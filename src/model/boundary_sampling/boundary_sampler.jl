using LinearAlgebra, Random, Statistics, Dates, DataFrames, CSV, YAML

#==========================================================================================

Part B: boundary sampler (Julia port of select_boundary_samples.py)

No model is trained. Everything is nearest-neighbour lookups and distance arithmetic in
the metric coordinates z from Part A (active_subspace.jl).

 1. Every final MGA solve is (almost always) a corner of the sampled hull, so the final
    solves are the boundary vertices. Each gets a region label (gradient_bank.jl).
 2. Find each vertex's K nearest neighbours in z.
 3. Candidates: for each region r with quota share q_r, n_r = round(q_r·N_new) picks are
    needed and n_r·(candidates per pick) candidates are generated. A candidate blends a
    base vertex from region r (weight w0 ~ U(lo, hi)) with m ∈ {1,2,3} of its neighbours
    (remaining 1 − w0 split by Dirichlet(1,…,1)). Blending NEARBY vertices keeps it just
    inside the boundary. Budgets are convex, so a blend of budget-feasible solves is
    budget-feasible.
 4. Greedy max-min ("farthest point first"), region by region in quota order: take the
    pool candidate farthest from everything solved or already picked, then shrink every
    candidate's distance to account for the new pick.
 5. Per pick, record lineage and bounds on its operating costs before it is solved:
      upper bound  blend = Σ w_j f_j        (op cost is convex in x: LP value function)
      lower bound  max_j [f_j + g_jᵀ(x − x_j)]  (each parent's Benders cut is a valid
                                                 under-estimator of a convex function)

==========================================================================================#

const BOUNDARY_REGION_ORDER = ["outer_shell", "ev_wing", "cvar_wing", "interior"]

# One candidate blend: parents[1] is the base vertex, weights sum to 1
struct BlendCandidate
    region::String
    parents::Vector{Int}
    weights::Vector{Float64}
end

# One selected interior portfolio, with everything needed to audit it after evaluation
struct BoundaryPick
    name::String
    wave::Int
    region::String                 # region of the base vertex
    parent_names::Vector{String}
    weights::Vector{Float64}
    dist_to_nearest_solve::Float64 # in z, to the reference set (finals + earlier waves)
    x::Vector{Float64}             # capacities to solve, original MW units
    z::Vector{Float64}
    blend_ev::Float64              # upper bounds (Σ w_j · parent value)
    blend_cvar::Float64
    blend_inv::Float64             # investment is linear, so this is exact
    lb_ev::Float64                 # first-order lower bounds (NaN without gradients)
    lb_cvar::Float64
end

#---- Distances ----------------------------------------------------------------------------
# Euclidean distances between the rows of A (n × r) and the rows of B (m × r)
function pairwise_distances(A::AbstractMatrix, B::AbstractMatrix)
    D2 = sum(abs2, A; dims = 2) .+ sum(abs2, B; dims = 2)' .- 2 .* (A * B')
    return sqrt.(max.(D2, 0.0))     # clip tiny negative round-off
end

# K nearest neighbours of every row of Z (excluding the row itself)
function knn_indices(Z::AbstractMatrix, K::Int)
    D = pairwise_distances(Z, Z)
    N = size(Z, 1)
    return [filter(!=(i), sortperm(D[i, :]))[1:min(K, N - 1)] for i in 1:N]
end

#---- Quotas -------------------------------------------------------------------------------
# Quotas are an ordered list of region => share pairs (order = selection order).
# round() ties to even, like Python's round in the reference script.
quota_counts(quotas, n_new::Int) = [(String(region), round(Int, share * n_new)) for (region, share) in quotas]

#---- Step 3: candidate generation ---------------------------------------------------------
function generate_candidates(rng::AbstractRNG, regions::Vector{String}, nbr::Vector{Vector{Int}}, counts, cfg::AbstractDict)
    lo, hi = cfg["Base weight range"]
    candidates = BlendCandidate[]
    for (region, n_r) in counts
        n_r == 0 && continue
        bases = findall(==(region), regions)
        if isempty(bases)
            @warn("No final solves in region $region; its $n_r picks are skipped.")
            continue
        end
        for _ in 1:(n_r * cfg["Candidates"])
            b = rand(rng, bases)                                  # base vertex, uniform in region
            m = min(rand(rng, 1:3), length(nbr[b]))               # 1 to 3 partners
            partners = nbr[b][randperm(rng, length(nbr[b]))[1:m]] # m distinct neighbours
            w0 = lo + (hi - lo) * rand(rng)                       # weight kept on the base
            e = randexp(rng, m)                                   # normalized Exp(1) draws
            w_rest = (1 - w0) .* e ./ sum(e)                      # = Dirichlet(1,…,1)·(1 − w0)
            push!(candidates, BlendCandidate(region, vcat(b, partners), vcat(w0, w_rest)))
        end
    end
    return candidates
end

# A blend's coordinates: the same weights applied to the parents' rows (z is linear in x)
blend_rows(M::AbstractMatrix, c::BlendCandidate) = vec(c.weights' * M[c.parents, :])

#---- Step 4: greedy max-min selection -----------------------------------------------------
# dmin[c] = distance from candidate c to the nearest reference point (solved portfolio) or
# already-picked candidate. Returns the picked candidate indices in pick order.
function greedy_maxmin(cand_Z::AbstractMatrix, cand_region::Vector{String}, ref_Z::AbstractMatrix, counts)
    dmin = vec(minimum(pairwise_distances(cand_Z, ref_Z); dims = 2))
    chosen = Int[]
    for (region, n_r) in counts
        pool = findall(==(region), cand_region)
        isempty(pool) && continue
        for _ in 1:n_r
            i = pool[argmax(dmin[pool])]        # farthest from everything so far
            push!(chosen, i)
            # the pick now counts as data: shrink EVERY candidate's distance to it
            dmin .= min.(dmin, vec(pairwise_distances(cand_Z, cand_Z[i:i, :])))
        end
    end
    return chosen
end

#---- Step 5: bounds -----------------------------------------------------------------------
# Lower bounds from the parents' Benders cuts. For parent j and scenario s:
#     θ_s(x) >= θ_{j,s} + β_{j,s}ᵀ (x − x_j)
# (written with unscaled β in MW; the same as g_{j,s}ᵀ(x̃ − x̃_j) because the cap_max
# scaling cancels). EV: the probability-weighted cut of each parent, max over parents.
# CVaR: the per-scenario max over parents l_s, then the CVaR tail mean of l. That is valid
# because CVaR is monotone in the scenario costs.
function cut_lower_bounds(parents::AbstractVector, x::Vector{Float64}, α::Real)
    all(has_gradients, parents) || return NaN, NaN
    p = parents[1].p ./ sum(parents[1].p)
    lin = [fs.θ .+ fs.β' * (x .- fs.x) for fs in parents]   # per parent: S-vector of cuts
    lb_ev = maximum(dot(p, l) for l in lin)
    l = reduce((a, b) -> max.(a, b), lin)
    lb_cvar = dot(cvar_tail_weights(l, p, α), l) / α
    return lb_ev, lb_cvar
end

#---- Top level ----------------------------------------------------------------------------
# finals:    vertices (Vector{FinalSolve}) with regions assigned; only `ok` ones are used.
# Z_finals:  their metric coordinates (rows aligned with finals).
# extra_ref_Z: coordinates of already-solved interior points (waves). They count for
#            distance but are never base vertices or parents.
function select_boundary_samples(finals::AbstractVector, Z_finals::AbstractMatrix, cfg::AbstractDict, rng::AbstractRNG, α::Real;
                                 n_new::Int, extra_ref_Z::AbstractMatrix = zeros(0, size(Z_finals, 2)),
                                 first_index::Int = 1, wave::Int = 1)
    keep = [fs.ok for fs in finals]
    finals, Z = finals[keep], Z_finals[keep, :]
    regions = [fs.region for fs in finals]
    X = reduce(vcat, [fs.x' for fs in finals])

    counts = quota_counts(cfg["Quotas"], n_new)
    nbr = knn_indices(Z, cfg["K neighbours"])
    candidates = generate_candidates(rng, regions, nbr, counts, cfg)
    isempty(candidates) && return BoundaryPick[]
    cand_Z = reduce(vcat, [blend_rows(Z, c)' for c in candidates])

    ref_Z = vcat(Z, extra_ref_Z)
    chosen = greedy_maxmin(cand_Z, [c.region for c in candidates], ref_Z, counts)

    ev = [solve_ev(fs) for fs in finals]
    cv = [solve_cvar(fs, α) for fs in finals]
    inv = [fs.inv_cost for fs in finals]

    picks = BoundaryPick[]
    for (k, i) in enumerate(chosen)
        c = candidates[i]
        w = c.weights
        x = blend_rows(X, c)                                  # Σ w_j x_j, MW
        lb_ev, lb_cvar = cut_lower_bounds(finals[c.parents], x, α)
        dist = minimum(pairwise_distances(cand_Z[i:i, :], ref_Z))
        push!(picks, BoundaryPick("Interior_$(first_index + k - 1)", wave, c.region, [fs.name for fs in finals[c.parents]], w,
                                  dist, x, cand_Z[i, :], dot(w, ev[c.parents]), dot(w, cv[c.parents]), dot(w, inv[c.parents]), lb_ev, lb_cvar))
    end
    return picks
end

#---- Diagnostics --------------------------------------------------------------------------
# Spread check from the spec: ratio of the picks' median distance-to-nearest-solve over the
# finals' median nearest-neighbour distance (should be ≈ 1), and the same ratio for the
# midpoints of random final pairs (should be larger: midpoints sit farther from the data).
function sampler_spread_diagnostics(picks::AbstractVector, Z_finals::AbstractMatrix, rng::AbstractRNG; n_midpoints::Int = 300)
    D = pairwise_distances(Z_finals, Z_finals)
    for i in 1:size(D, 1); D[i, i] = Inf; end
    nn_median = median(vec(minimum(D; dims = 2)))
    pick_median = median([p.dist_to_nearest_solve for p in picks])
    N = size(Z_finals, 1)
    mids = reduce(vcat, [((Z_finals[a, :] .+ Z_finals[b, :]) ./ 2)' for (a, b) in ((rand(rng, 1:N), rand(rng, 1:N)) for _ in 1:n_midpoints)])
    mid_median = median(vec(minimum(pairwise_distances(mids, Z_finals); dims = 2)))
    return (finals_nn_median = nn_median, picks_median = pick_median, ratio = pick_median / nn_median,
            midpoint_median = mid_median, midpoint_ratio = mid_median / nn_median)
end

# Spearman rank correlation (average ranks for ties), to compare two metrics' distances
function tied_ranks(v::AbstractVector)
    order = sortperm(v)
    ranks = zeros(length(v))
    i = 1
    while i <= length(v)
        j = i
        while j < length(v) && v[order[j + 1]] == v[order[i]]
            j += 1
        end
        ranks[order[i:j]] .= (i + j) / 2
        i = j + 1
    end
    return ranks
end
spearman_correlation(a::AbstractVector, b::AbstractVector) = cor(tied_ranks(a), tied_ranks(b))
