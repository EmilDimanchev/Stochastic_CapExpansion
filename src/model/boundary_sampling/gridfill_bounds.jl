using LinearAlgebra, Random, Statistics, JuMP

#==========================================================================================

GridFill: the per-scenario cut bank and the op-cost bracket of a blend

Every subproblem solve at a capacity x_j gives one supporting hyperplane per scenario,
    θ_s(x) >= θ_s(x_j) + β_sᵀ (x − x_j) = a + gᵀx,
valid everywhere because each scenario's op cost is convex in capacity. The CutBank holds
these per scenario, from two sources:
  - the CutPool (live + archived optimality_cut_* rows of the master): every Benders
    iterate of every anchor / MGA / wing solve, not just the final ones;
  - FinalSolves: the bank's final solves and every evaluated GridFill sample.
All values are unscaled ($ and $/MW), in the [x; x_line] order of capacity_variable_names.

Bracket for a blend x = Σ w_j x_j of solved parents (θ_j their per-scenario costs):
  upper: θ_s(x) <= Σ w_j θ_{j,s} =: θ̄_s                  (convexity)
         EV <= p·θ̄,  CVaR <= CVaR(θ̄)  (tighter than Σ w_j CVaR_j)
  lower: θ_s(x) >= l_s := max over scenario-s cuts
         EV >= p·l,  CVaR >= CVaR(l)  (CVaR is monotone in the scenario costs)
The midpoint of each bracket is the aiming predictor; it is never used as a label.

==========================================================================================#

mutable struct CutBank
    a::Vector{Vector{Float64}}     # per scenario: cut intercepts, $
    G::Vector{Matrix{Float64}}     # per scenario: d × n_s cut gradients, $/MW
end

CutBank(S::Int, d::Int) = CutBank([Float64[] for _ in 1:S], [zeros(d, 0) for _ in 1:S])

n_cuts(cb::CutBank) = sum(length, cb.a)

# One cut per scenario from a solve at fs.x (only solves with gradients and all scenarios optimal)
function add_final_solve!(cb::CutBank, fs)
    (fs.ok && has_gradients(fs)) || return cb
    for s in eachindex(cb.a)
        g = fs.β[:, s]
        push!(cb.a[s], fs.θ[s] - dot(g, fs.x))
        cb.G[s] = hcat(cb.G[s], g)
    end
    return cb
end

# Every optimality cut in the pool, live or archived. The master row is
#     c·alpha[s,f,k] + Σ a_r x_r + Σ a_l x_line_l >= b      (scaled by "Scaling factor cost")
# with c = 1 for an optimal subproblem (c = 0 rows carry no information and are skipped).
# alpha is the scaled op cost of scenario (s,f,k), so θ_s(x) >= scale·(b − aᵀx)/c.
# cvar_tail_cuts_* repeat the same tangent with −ζ and are skipped.
function add_pool_cuts!(cb::CutBank, pool, MP::Model, settings::Dict)
    scale = settings["Scaling factor cost"]
    alpha = MP[:alpha]
    xvars = vcat(vec(collect(MP[:x])), vec(collect(MP[:x_line])))
    scen_of = Dict(alpha[I] => li for (li, I) in enumerate(CartesianIndices(size(alpha))))   # (s,f,k) column-major, same as θ
    pos_of = Dict(v => i for (i, v) in enumerate(xvars))
    d = length(xvars)
    a_new = [Float64[] for _ in eachindex(cb.a)]
    g_new = [Vector{Float64}[] for _ in eachindex(cb.a)]
    added = 0
    function add_row!(terms, b)
        s, c = 0, 0.0
        g = zeros(d)
        for (v, coef) in terms
            if haskey(scen_of, v)
                s, c = scen_of[v], coef
            elseif haskey(pos_of, v)
                g[pos_of[v]] = coef
            end
        end
        (s == 0 || c <= 0) && return
        push!(a_new[s], scale * b / c)
        push!(g_new[s], -scale .* g ./ c)
        added += 1
    end
    for (n, con) in pool.live
        startswith(n, "optimality_cut_") || continue
        obj = constraint_object(con)
        add_row!(obj.func.terms, MOI.constant(obj.set) - obj.func.constant)
    end
    for j in eachindex(pool.names)
        (pool.alive[j] && startswith(pool.names[j], "optimality_cut_")) || continue
        rng = pool.colptr[j]:(pool.colptr[j + 1] - 1)
        add_row!(((pool.vars[pool.rowval[k]], pool.nzval[k]) for k in rng), pool.b[j])
    end
    for s in eachindex(cb.a)
        isempty(a_new[s]) && continue
        append!(cb.a[s], a_new[s])
        cb.G[s] = hcat(cb.G[s], reduce(hcat, g_new[s]))
    end
    @info("GridFill cut bank: $added optimality cuts taken from the cut pool ($(length(pool.live)) live, $(n_archived(pool)) archived cuts in the pool)")
    return cb
end

# Per-scenario cut maximum at each column of X (d × m): S × m. −Inf where a scenario has no cut.
function cut_values(cb::CutBank, X::AbstractMatrix)
    l = fill(-Inf, length(cb.a), size(X, 2))
    for s in eachindex(cb.a)
        isempty(cb.a[s]) && continue
        l[s, :] = vec(maximum(cb.a[s] .+ cb.G[s]' * X; dims = 1))
    end
    return l
end

# Keep, per scenario, only the cuts that attain the maximum at one of the reference points
# (columns of Xref). Any subset of cuts still gives a valid lower bound; keeping those active
# around the parents keeps the bound almost as tight at a fraction of the cost.
function prune_cuts(cb::CutBank, Xref::AbstractMatrix)
    out = CutBank(length(cb.a), size(Xref, 1))
    for s in eachindex(cb.a)
        isempty(cb.a[s]) && continue
        keep = Int[]
        for cols in Iterators.partition(axes(Xref, 2), 500)        # bounds the n_s × chunk temporary
            append!(keep, vec(getindex.(argmax(cb.a[s] .+ cb.G[s]' * Xref[:, cols]; dims = 1), 1)))
        end
        unique!(sort!(keep))
        out.a[s] = cb.a[s][keep]
        out.G[s] = cb.G[s][:, keep]
    end
    return out
end

# CVaR of each column of L (S × m) under probabilities p
cvar_columns(L::AbstractMatrix, p::AbstractVector, α::Real) = [dot(cvar_tail_weights(c, p, α), c) / α for c in eachcol(L)]

# Bracket for m blends at once. Θ: S × n parent op costs; inv: n investment costs; parents:
# m × 3 parent indices; W: m × 3 weights; X: d × m blended capacities.
# Returns a NamedTuple of m-vectors (ub/lb/mid for EV and CVaR, and the exact inv).
function blend_bracket(cb::CutBank, Θ::AbstractMatrix, inv::AbstractVector, p::AbstractVector, α::Real,
                       parents::AbstractMatrix{Int}, W::AbstractMatrix, X::AbstractMatrix)
    m = size(parents, 1)
    pn = p ./ sum(p)
    θbar = zeros(size(Θ, 1), m)
    for i in 1:m, k in axes(parents, 2)
        θbar[:, i] .+= W[i, k] .* view(Θ, :, parents[i, k])
    end
    L = cut_values(cb, X)
    L = min.(L, θbar)                     # a cut can never exceed the true cost; guards numerical noise
    ub_ev, lb_ev = vec(pn' * θbar), vec(pn' * L)
    ub_cvar, lb_cvar = cvar_columns(θbar, pn, α), cvar_columns(L, pn, α)
    inv_b = [sum(W[i, k] * inv[parents[i, k]] for k in axes(parents, 2)) for i in 1:m]
    return (inv = inv_b, ub_ev = ub_ev, lb_ev = lb_ev, mid_ev = (ub_ev .+ lb_ev) ./ 2,
            ub_cvar = ub_cvar, lb_cvar = lb_cvar, mid_cvar = (ub_cvar .+ lb_cvar) ./ 2)
end
