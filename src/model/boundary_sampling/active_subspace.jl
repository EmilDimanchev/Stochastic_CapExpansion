using LinearAlgebra, Random, Statistics, Dates, DataFrames, CSV, YAML

#==========================================================================================

Part A: the active-subspace metric

Goal: a low-dimensional coordinate z for each capacity portfolio such that Euclidean
distance in z is a first-order estimate of how much EV / CVaR operating cost changes
between two portfolios, measured in "chord units" (fractions of the anchor gaps).

Steps (one function each below):
 1. Scale capacities to comparable units:  x̃ = x ./ cap_max   (cap_max = max over finals)
    so the gradient w.r.t. x̃ is  g = β .* cap_max.
 2. Per solve i, two gradients:
      g_EV,i   = Σ_s p_s g_{i,s}                    (expected-value gradient)
      g_CVaR,i = Σ_s w_s g_{i,s} / α                (CVaR subgradient; w = tail weights
                                                     from the same rule as compute_cvar)
    divided by the anchor gaps ΔC = C_R − C_E and ΔV = V_E − V_R so both are in chord units.
 3. M = (1/2N) Σ_i (g_EV,i g_EV,iᵀ + g_CVaR,i g_CVaR,iᵀ) = GᵀG, where G stacks all 2N
    gradients as rows, scaled by 1/sqrt(2N). For large d we never form M and take the thin
    SVD of G instead: same eigenvectors, eigenvalues = singular values squared.
 4. Keep the top r eigenpairs (smallest r capturing "Energy share" of the eigenvalue sum,
    capped at "r max").
 5. z = diag(sqrt(λ_1..λ_r)) W_rᵀ x̃.  Then ||z_a − z_b||² = Δx̃ᵀ (W_r Λ_r W_rᵀ) Δx̃ ≈ Δx̃ᵀ M Δx̃,
    the mean squared first-order change in (EV, CVaR) cost, in chord units.
 6. Fallback ("lengthscale_file"): z = x ./ (cap_max_MW .* length_scale) from a CSV.

==========================================================================================#

struct CapacityMetric
    mode::String                      # "active_subspace" or "lengthscale_file"
    var_names::Vector{String}
    cap_max::Vector{Float64}          # d; scaling used for x̃ (active subspace)
    scale::Vector{Float64}            # d; lengthscale mode: z = x ./ scale
    W::Matrix{Float64}                # d × r eigenvectors (active subspace only)
    λ_all::Vector{Float64}            # full eigenvalue spectrum, descending
    r::Int
    ΔC::Float64                       # anchor gaps (NaN in lengthscale mode)
    ΔV::Float64
    n_solves::Int
    α::Float64
    energy_share::Float64
end

#---- Step 1: capacity scale ---------------------------------------------------------------
# cap_max per variable over the given capacity matrix (rows = solves). A variable that is
# always ~0 MW gets a scale of 1, so dividing by it cannot blow up.
function capacity_scale(X::AbstractMatrix; tiny::Float64 = 1e-3)
    cap_max = vec(maximum(X; dims = 1))
    return [c < tiny ? 1.0 : c for c in cap_max]
end

#---- Step 2: per-solve EV / CVaR gradients in scaled units --------------------------------
# (fs is a FinalSolve; left untyped because the module loads files alphabetically and
#  gradient_bank.jl, which defines FinalSolve, comes after this file.)
function solve_gradients(fs, cap_max::Vector{Float64}, α::Real)
    G = fs.β .* cap_max                         # d × S, column s = g_s = β_s .* cap_max
    p = fs.p ./ sum(fs.p)
    g_ev = G * p
    g_cvar = G * cvar_tail_weights(fs.θ, fs.p, α) ./ α
    return g_ev, g_cvar
end

# Anchor gaps from the two extreme risk-weight solves:
#   E = Risk_Weight_1.0 (EV-optimal), R = Risk_Weight_0.0 (CVaR-optimal)
#   C = Inv + EV op cost,  V = Inv + CVaR op cost;  ΔC = C_R − C_E > 0, ΔV = V_E − V_R > 0
function anchor_gaps(bank::AbstractVector, α::Real)
    find(name) = (i = findfirst(fs -> fs.name == name, bank); i === nothing ? error("Anchor $name not found in the gradient bank; it is needed for the chord units.") : bank[i])
    E, R = find("Risk_Weight_1.0"), find("Risk_Weight_0.0")
    ΔC = (R.inv_cost + solve_ev(R)) - (E.inv_cost + solve_ev(E))
    ΔV = (E.inv_cost + solve_cvar(E, α)) - (R.inv_cost + solve_cvar(R, α))
    (ΔC > 0 && ΔV > 0) || error("Anchor gaps must be positive (got ΔC = $ΔC, ΔV = $ΔV).")
    return ΔC, ΔV
end

#---- Step 3: principal directions of M = GᵀG ----------------------------------------------
# Returns (λ, W) with λ descending and W's columns the matching unit eigenvectors.
# explicit = true forms the d × d matrix M; false uses the thin SVD of G (2N × d), which
# is what you want once d is in the thousands. Both give the same result (tested).
function principal_directions(G::AbstractMatrix; explicit::Bool)
    if explicit
        E = eigen(Symmetric(G' * G))
        order = sortperm(E.values; rev = true)
        return max.(E.values[order], 0.0), E.vectors[:, order]   # clip tiny negative round-off
    else
        F = svd(G)                                               # thin SVD: G = U Σ Vᵀ
        return F.S .^ 2, F.V                                     # eigen of GᵀG = V Σ² Vᵀ
    end
end

#---- Step 4: how many directions to keep --------------------------------------------------
function choose_rank(λ::AbstractVector, energy_share::Real, r_max::Int)
    total = sum(λ)
    total <= 0 && return 1
    share = cumsum(λ) ./ total
    r = something(findfirst(>=(energy_share - 1e-12), share), length(λ))
    return clamp(r, 1, r_max)
end

# Build the active-subspace metric from the solves in `bank` that have usable gradients.
function build_active_subspace(bank::AbstractVector, var_names::Vector{String}, cfg::AbstractDict, α::Real)
    usable = [fs for fs in bank if fs.ok && has_gradients(fs)]
    isempty(usable) && error("No solves with gradients available to build the active subspace.")
    X = reduce(vcat, [fs.x' for fs in usable])
    cap_max = capacity_scale(X)
    ΔC, ΔV = anchor_gaps(bank, α)

    # Stack the 2N chord-unit gradients as rows of G, scaled so that GᵀG = M
    N = length(usable)
    rows = Vector{Vector{Float64}}()
    for fs in usable
        g_ev, g_cvar = solve_gradients(fs, cap_max, α)
        push!(rows, g_ev ./ ΔC, g_cvar ./ ΔV)
    end
    G = reduce(vcat, [r' for r in rows]) ./ sqrt(2N)

    d = length(var_names)
    λ, W = principal_directions(G; explicit = d <= cfg["Explicit M max dim"])
    r = choose_rank(λ, cfg["Energy share"], cfg["r max"])

    @info("Active subspace built from $N solves (d = $d): keeping r = $r directions " *
          "(energy share target $(cfg["Energy share"]), r max $(cfg["r max"])).")
    @info("Eigenvalue spectrum (descending): " * join(string.(round.(λ; sigdigits = 4)), ", "))
    @info("Cumulative share: " * join(string.(round.(cumsum(λ) ./ sum(λ); digits = 4)), ", "))

    return CapacityMetric("active_subspace", var_names, cap_max, cap_max, W[:, 1:r], λ, r, ΔC, ΔV, N, Float64(α), Float64(cfg["Energy share"]))
end

#---- Step 6: fallback metric from a length-scale file -------------------------------------
# CSV columns: tech, cap_max_MW, length_scale (one row per capacity variable, any order).
function lengthscale_metric(csv_path::String, var_names::Vector{String})
    df = CSV.read(csv_path, DataFrame)
    lookup = Dict(string(t) => i for (i, t) in enumerate(df.tech))
    missing_vars = [v for v in var_names if !haskey(lookup, v)]
    isempty(missing_vars) || error("Metric file $csv_path has no row for: $(missing_vars)")
    rows = [lookup[v] for v in var_names]
    cap_max = Float64.(df.cap_max_MW[rows])
    scale = cap_max .* Float64.(df.length_scale[rows])
    d = length(var_names)
    return CapacityMetric("lengthscale_file", var_names, cap_max, scale, zeros(d, 0), Float64[], d, NaN, NaN, 0, NaN, NaN)
end

#---- Step 5: metric coordinates -----------------------------------------------------------
# X has one portfolio per row (N × d, original MW units); returns Z (N × r).
function metric_coords(m::CapacityMetric, X::AbstractMatrix)
    if m.mode == "active_subspace"
        Xt = X ./ m.cap_max'                                   # x̃
        return (Xt * m.W) .* sqrt.(m.λ_all[1:m.r])'            # z = diag(√λ) Wᵀ x̃, as rows
    else
        return X ./ m.scale'
    end
end
metric_coords(m::CapacityMetric, x::AbstractVector) = vec(metric_coords(m, reshape(x, 1, :)))

# Pick the metric named in the boundary settings ("Metric"). The active subspace needs no
# prior knowledge of the model. Only if no solve has usable gradients does it fall back to
# a length-scale file, and only when one was given explicitly; otherwise it stops with an error.
function build_metric(bank::AbstractVector, var_names::Vector{String}, cfg::AbstractDict, α::Real)
    file = cfg["Metric file"]
    if cfg["Metric"] == "active_subspace"
        if any(fs -> fs.ok && has_gradients(fs), bank)
            return build_active_subspace(bank, var_names, cfg, α)
        end
        file === nothing && error("Cannot build the active-subspace metric: no final solve has usable gradients " *
                                  "(every solve had a subproblem that did not solve to optimality), and no \"Metric file\" is set to fall back on.")
        @warn("No gradients available; falling back to the length-scale metric file $file.")
        return lengthscale_metric(file, var_names)
    elseif cfg["Metric"] == "lengthscale_file"
        file === nothing && error("\"Metric\" => \"lengthscale_file\" needs \"Metric file\" set to a CSV with columns tech, cap_max_MW, length_scale.")
        return lengthscale_metric(file, var_names)
    else
        error("Unknown metric \"$(cfg["Metric"])\"; use \"active_subspace\" or \"lengthscale_file\".")
    end
end

#---- Saving the metric --------------------------------------------------------------------
# <folder>/<tag>/eigenvalues.csv, eigenvectors.csv, scaling.csv, metadata.yaml
function write_metric(folder::String, m::CapacityMetric; tag::String = "wave_0")
    out = joinpath(folder, tag)
    mkpath(out)
    CSV.write(joinpath(out, "scaling.csv"), DataFrame(variable = m.var_names, cap_max = m.cap_max, scale = m.scale))
    if m.mode == "active_subspace"
        share = cumsum(m.λ_all) ./ sum(m.λ_all)
        CSV.write(joinpath(out, "eigenvalues.csv"), DataFrame(index = 1:length(m.λ_all), eigenvalue = m.λ_all, cumulative_share = share, kept = (1:length(m.λ_all)) .<= m.r))
        df_W = DataFrame(variable = m.var_names)
        for j in 1:m.r
            df_W[!, "w$j"] = m.W[:, j]
        end
        CSV.write(joinpath(out, "eigenvectors.csv"), df_W)
    end
    commit = try readchomp(`git -C $(@__DIR__) rev-parse HEAD`) catch; "unknown" end
    meta = Dict("mode" => m.mode, "r" => m.r, "n_solves" => m.n_solves, "d" => length(m.var_names),
                "anchor_gap_EV" => m.ΔC, "anchor_gap_CVaR" => m.ΔV, "alpha_tail_mass" => m.α,
                "energy_share" => m.energy_share, "git_commit" => commit, "written" => string(now()))
    YAML.write_file(joinpath(out, "metadata.yaml"), meta)
    return out
end
