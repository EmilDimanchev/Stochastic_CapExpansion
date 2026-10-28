using LinearAlgebra, Random, Statistics, Distributed

#==========================================================================================

GridFill: cost-space targets, candidate blends and selection

Cost coordinates (anchors E = Risk_Weight_1.0, EV-optimal; R = Risk_Weight_0.0, CVaR-optimal):
    C = Inv + EV op cost,  V = Inv + CVaR op cost
    ĉ = (C − C_E)/(C_R − C_E),  v̂ = (V − V_R)/(V_E − V_R)
EV and CVaR always come from a solve's own per-scenario costs θ (solve_ev / solve_cvar),
never from the Summary's CVaR_SystemCost (a risk-weighted blend).

1. Targets: a grid of cell size g over (ĉ, v̂), restricted to cells above the risk-weight
   sweep (piecewise-linear through the Risk_Weight_* points) and with ĉ + v̂ <= Delta cap.
   A cell is a target if no solved point sat in it when the phase started; it stays open
   until "Per cell" solved points sit in it.
2. Candidates, one task per open cell (cell_candidates): aim points spread over the cell;
   random parent triples from ALL solved points; the barycentric coordinates of the aim in
   the triple's cost triangle are the blend weights. The weights are then re-aimed so the
   bracket MIDPOINT (gridfill_bounds.jl), not the raw blend, lands on the aim. Seeded per
   (Seed, cell, round), so results don't depend on how cells are spread over workers.
3. Selection (select_gridfill): greedy max-min distance in range-scaled capacity space with
   a per-cell quota, against every solved point and every earlier pick in any cell. The
   score is a function argument (novelty_score / bracket_novelty_score / your own).

==========================================================================================#

# Internal tuning (not in settings.yaml); override with run_gridfill_sampling!(...; tuning = ...).
const GRIDFILL_TUNING = (
    aims_per_cell = 6,        # aim points spread over each cell (uniform in its central 80%)
    triples_per_aim = 8000,   # random parent triples tested per aim point
    max_weight = 0.95,        # drop blends that are almost a copy of one parent
    reaim_pool = 300,         # most novel valid blends per cell that get re-aimed
    reaim_iters = 3,          # aim' += aim − predicted midpoint, this many times
    prune_blends = 2000,      # random parent blends used (with the solved points) to pick the cuts that matter
)

#---- Cost coordinates ---------------------------------------------------------------------
struct GridFillAnchors
    C_E::Float64
    C_R::Float64
    V_E::Float64
    V_R::Float64
end

function gridfill_anchors(solves::AbstractVector, α::Real)
    find(name) = (i = findfirst(fs -> fs.name == name, solves); i === nothing ? error("GridFill needs anchor $name among the solved points.") : solves[i])
    E, R = find("Risk_Weight_1.0"), find("Risk_Weight_0.0")
    a = GridFillAnchors(E.inv_cost + solve_ev(E), R.inv_cost + solve_ev(R), E.inv_cost + solve_cvar(E, α), R.inv_cost + solve_cvar(R, α))
    (a.C_R > a.C_E && a.V_E > a.V_R) || error("GridFill anchor gaps must be positive (C_R − C_E = $(a.C_R - a.C_E), V_E − V_R = $(a.V_E - a.V_R)).")
    return a
end

c_hat(a::GridFillAnchors, C) = (C .- a.C_E) ./ (a.C_R - a.C_E)
v_hat(a::GridFillAnchors, V) = (V .- a.V_R) ./ (a.V_E - a.V_R)

# n × 2 matrix of (ĉ, v̂) for solved points
cost_coords(solves::AbstractVector, a::GridFillAnchors, α::Real) =
    hcat(c_hat(a, [fs.inv_cost + solve_ev(fs) for fs in solves]), v_hat(a, [fs.inv_cost + solve_cvar(fs, α) for fs in solves]))

#---- Grid ---------------------------------------------------------------------------------
cell_of(k::AbstractVector, g::Real) = (floor(Int, k[1] / g), floor(Int, k[2] / g))

function cell_counts(K::AbstractMatrix, g::Real)
    counts = Dict{Tuple{Int,Int}, Int}()
    for k in eachrow(K)
        c = cell_of(k, g)
        counts[c] = get(counts, c, 0) + 1
    end
    return counts
end

# Lower frontier v̂(ĉ) of the risk-weight sweep: linear between its points, flat beyond them
function sweep_frontier(sweep::AbstractMatrix)
    o = sortperm(sweep[:, 1])
    cs, vs = sweep[o, 1], sweep[o, 2]
    return c -> c <= cs[1] ? vs[1] : c >= cs[end] ? vs[end] : begin
        i = searchsortedlast(cs, c)
        t = (c - cs[i]) / (cs[i + 1] - cs[i])
        (1 - t) * vs[i] + t * vs[i + 1]
    end
end

# Every cell of the near-optimal region: centre above the sweep frontier, (i+j+1)·g <= δ cap
function grid_cells(sweep::AbstractMatrix, g::Real, δcap::Real)
    frontier = sweep_frontier(sweep)
    n = floor(Int, δcap / g + 1e-9)
    return [(i, j) for i in 0:n-1 for j in 0:n-1 if (i + j + 1) * g <= δcap + 1e-9 && (j + 0.5) * g >= frontier((i + 0.5) * g)]
end

#---- 2-D hull (report unreachable cells) -------------------------------------------------
cross2(o, a, b) = (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1])

# Counter-clockwise convex hull of the rows of K (Andrew's monotone chain)
function convex_hull_2d(K::AbstractMatrix)
    P = sort([Tuple(r) for r in eachrow(K)])
    length(P) <= 2 && return P
    lower, upper = NTuple{2,Float64}[], NTuple{2,Float64}[]
    for p in P
        while length(lower) >= 2 && cross2(lower[end-1], lower[end], p) <= 0; pop!(lower); end
        push!(lower, p)
    end
    for p in reverse(P)
        while length(upper) >= 2 && cross2(upper[end-1], upper[end], p) <= 0; pop!(upper); end
        push!(upper, p)
    end
    return vcat(lower[1:end-1], upper[1:end-1])
end

in_hull(hull, t) = length(hull) >= 3 && all(cross2(hull[i], hull[mod1(i + 1, length(hull))], t) >= -1e-12 for i in eachindex(hull))

#---- Barycentric weights ------------------------------------------------------------------
# Weights of aim t[i, :] in the triangle K[tri[i, :], :], for every row i (m × 3).
# Degenerate triangles give non-finite weights (filtered by valid_weights).
function barycentric(K::AbstractMatrix, tri::AbstractMatrix{Int}, t::AbstractMatrix)
    m = size(tri, 1)
    W = zeros(m, 3)
    for i in 1:m
        a, b, c = view(K, tri[i, 1], :), view(K, tri[i, 2], :), view(K, tri[i, 3], :)
        v0, v1, v2 = b .- a, c .- a, view(t, i, :) .- a
        den = v0[1] * v1[2] - v1[1] * v0[2]
        w1 = (v2[1] * v1[2] - v1[1] * v2[2]) / den
        w2 = (v0[1] * v2[2] - v2[1] * v0[2]) / den
        W[i, 1], W[i, 2], W[i, 3] = 1 - w1 - w2, w1, w2
    end
    return W
end

valid_weights(W::AbstractMatrix, max_weight::Real) = [all(isfinite, w) && all(>=(0), w) && maximum(w) < max_weight for w in eachrow(W)]

#---- Context shipped to the candidate tasks ----------------------------------------------
struct GridFillContext
    X::Matrix{Float64}          # n × d solved capacities (parents and novelty references), MW
    Θ::Matrix{Float64}          # S × n per-scenario op costs
    inv::Vector{Float64}        # n investment costs
    K::Matrix{Float64}          # n × 2 solved (ĉ, v̂)
    scale::Vector{Float64}      # d capacity ranges for distances
    p::Vector{Float64}
    α::Float64
    anchors::GridFillAnchors
    cuts::CutBank
    g::Float64
    hull::Vector{NTuple{2,Float64}}
end

function GridFillContext(solves::AbstractVector, cuts::CutBank, anchors::GridFillAnchors, α::Real, g::Real)
    X = reduce(vcat, [fs.x' for fs in solves])
    K = cost_coords(solves, anchors, α)
    rng = vec(maximum(X; dims = 1) .- minimum(X; dims = 1))
    scale = [r > 0 ? r : 1.0 for r in rng]
    return GridFillContext(X, reduce(hcat, [fs.θ for fs in solves]), [fs.inv_cost for fs in solves], K, scale,
                           solves[1].p ./ sum(solves[1].p), Float64(α), anchors, cuts, Float64(g), convex_hull_2d(K))
end

# Predicted (ĉ, v̂) of the bracket midpoints (m × 2)
predicted_position(ctx::GridFillContext, br) = hcat(c_hat(ctx.anchors, br.inv .+ br.mid_ev), v_hat(ctx.anchors, br.inv .+ br.mid_cvar))

# m × d blended capacities Σ_k W[i,k] X[tri[i,k], :]
function blend_capacities(X::AbstractMatrix, tri::AbstractMatrix{Int}, W::AbstractMatrix)
    Y = zeros(size(tri, 1), size(X, 2))
    for i in axes(tri, 1), k in axes(tri, 2)
        Y[i, :] .+= W[i, k] .* view(X, tri[i, k], :)
    end
    return Y
end

# Distance from each row of Y (MW) to the nearest row of Ref (MW), range-scaled
nearest_distance(Y::AbstractMatrix, Ref::AbstractMatrix, scale::AbstractVector) =
    size(Ref, 1) == 0 ? fill(Inf, size(Y, 1)) : vec(minimum(pairwise_distances(Y ./ scale', Ref ./ scale'); dims = 2))

#---- Step 2: candidates for one cell ------------------------------------------------------
# status: :ok, :aim_missed (blends exist but no midpoint landed in the cell),
# :no_blend (no valid triangle around any aim, inside the parents' cost hull) or
# :outside_hull (cell centre outside the hull of the solved cost points: not fillable by blending).
struct GridFillCellResult
    cell::Tuple{Int,Int}
    status::Symbol
    tri::Matrix{Int}            # m × 3 parent indices (into the context's solves)
    W::Matrix{Float64}          # m × 3 weights
    aim::Matrix{Float64}        # m × 2 aim points (ĉ, v̂)
    X::Matrix{Float64}          # m × d blended capacities
    pred::Matrix{Float64}       # m × 2 predicted midpoint position
    bracket::NamedTuple         # m-vectors: inv, ub/lb/mid for EV and CVaR
    novelty::Vector{Float64}    # range-scaled distance to the nearest solved point
end

function empty_cell_result(cell, status, d)
    br = (inv = Float64[], ub_ev = Float64[], lb_ev = Float64[], mid_ev = Float64[], ub_cvar = Float64[], lb_cvar = Float64[], mid_cvar = Float64[])
    return GridFillCellResult(cell, status, zeros(Int, 0, 3), zeros(0, 3), zeros(0, 2), zeros(0, d), zeros(0, 2), br, Float64[])
end

subset_bracket(br, idx) = map(v -> v[idx], br)

function cell_candidates(ctx::GridFillContext, cell::Tuple{Int,Int}, tuning, n_keep::Int, seed::Integer, rnd::Integer)
    rng = Xoshiro(hash((seed, cell[1], cell[2], rnd)))
    g, n, d = ctx.g, size(ctx.X, 1), size(ctx.X, 2)
    tris, aims = Matrix{Int}[], Matrix{Float64}[]
    for _ in 1:tuning.aims_per_cell
        aim = (collect(cell) .+ 0.1 .+ 0.8 .* rand(rng, 2)) .* g
        tri = rand(rng, 1:n, tuning.triples_per_aim, 3)
        A = repeat(aim', size(tri, 1))
        ok = valid_weights(barycentric(ctx.K, tri, A), tuning.max_weight)
        push!(tris, tri[ok, :]); push!(aims, A[ok, :])
    end
    tri, aim = reduce(vcat, tris), reduce(vcat, aims)
    if isempty(tri)
        centre = ((cell[1] + 0.5) * g, (cell[2] + 0.5) * g)
        return empty_cell_result(cell, in_hull(ctx.hull, centre) ? :no_blend : :outside_hull, d)
    end

    # Pre-screen: only the most novel blends are worth re-aiming
    W = barycentric(ctx.K, tri, aim)
    nov = nearest_distance(blend_capacities(ctx.X, tri, W), ctx.X, ctx.scale)
    top = sortperm(nov; rev = true)[1:min(tuning.reaim_pool, length(nov))]
    tri, aim = tri[top, :], aim[top, :]

    # Re-aim: shift each blend's aim by the miss of its predicted midpoint, same triple
    aim_eff = copy(aim)
    local X, br, P
    for it in 0:tuning.reaim_iters
        if it > 0
            aim_eff .+= aim .- P
        end
        W = barycentric(ctx.K, tri, aim_eff)
        ok = valid_weights(W, tuning.max_weight)
        tri, aim, aim_eff, W = tri[ok, :], aim[ok, :], aim_eff[ok, :], W[ok, :]
        isempty(tri) && return empty_cell_result(cell, :aim_missed, d)
        X = blend_capacities(ctx.X, tri, W)
        br = blend_bracket(ctx.cuts, ctx.Θ, ctx.inv, ctx.p, ctx.α, tri, W, permutedims(X))
        P = predicted_position(ctx, br)
    end

    landed = [cell_of(P[i, :], g) == cell for i in axes(P, 1)]
    any(landed) || return empty_cell_result(cell, :aim_missed, d)
    nov = nearest_distance(X, ctx.X, ctx.scale)
    keep = filter(i -> landed[i], sortperm(nov; rev = true))
    keep = keep[1:min(n_keep, length(keep))]
    return GridFillCellResult(cell, :ok, tri[keep, :], W[keep, :], aim[keep, :], X[keep, :], P[keep, :], subset_bracket(br, keep), nov[keep])
end

#---- Step 2, distributed: one remotecall per worker over a static slice of the cells -------
# The context lives in a worker global for the duration of one round and is cleared right
# after, so nothing accumulates across rounds (no pmap / closures holding on to it).
const GRIDFILL_WORKER_CTX = Ref{Any}(nothing)

_set_gridfill_ctx!(ctx) = (GRIDFILL_WORKER_CTX[] = ctx; nothing)
_clear_gridfill_ctx!() = (GRIDFILL_WORKER_CTX[] = nothing; nothing)
_gridfill_cells_on_worker(cells, tuning, n_keep, seed, rnd) = [cell_candidates(GRIDFILL_WORKER_CTX[], c, tuning, n_keep, seed, rnd) for c in cells]

function run_candidate_stage(ctx::GridFillContext, cells::Vector{Tuple{Int,Int}}, tuning, n_keep::Int, seed::Integer, rnd::Integer; parallel::Bool = nworkers() > 0 && workers() != [myid()])
    parallel || return [cell_candidates(ctx, c, tuning, n_keep, seed, rnd) for c in cells]
    pids = workers()
    out = Vector{Vector{GridFillCellResult}}(undef, length(pids))
    try
        @sync for pid in pids
            @async remotecall_wait(_set_gridfill_ctx!, pid, ctx)
        end
        @sync for (w, pid) in enumerate(pids)
            @async out[w] = remotecall_fetch(_gridfill_cells_on_worker, pid, cells[w:length(pids):end], tuning, n_keep, seed, rnd)
        end
    finally
        @sync for pid in pids
            @async remotecall_wait(_clear_gridfill_ctx!, pid)
        end
        _gc_all_processes!()
    end
    results = reduce(vcat, out)
    return results[sortperm([r.cell for r in results])]
end

#---- Step 3: selection --------------------------------------------------------------------
# Scores take the current nearest distances and the pooled candidates; a larger score is
# picked first. novelty_score is pure max-min; bracket_novelty_score weights it by the
# bracket width in cost-space units (an acquisition-style score).
novelty_score(dmin, cands, ctx) = dmin
function bracket_novelty_score(dmin, cands, ctx)
    a = ctx.anchors
    width = (cands.ub_ev .- cands.lb_ev) ./ (a.C_R - a.C_E) .+ (cands.ub_cvar .- cands.lb_cvar) ./ (a.V_E - a.V_R)
    return dmin .* width
end

# Pool the cell results into flat arrays (cell index per row).
function pool_candidates(results::Vector{GridFillCellResult})
    ok = [r for r in results if r.status == :ok]
    isempty(ok) && return nothing
    stack_rows(f) = reduce(vcat, [f(r) for r in ok])
    return (cell = reduce(vcat, [fill(r.cell, size(r.tri, 1)) for r in ok]), tri = stack_rows(r -> r.tri), W = stack_rows(r -> r.W),
            aim = stack_rows(r -> r.aim), X = stack_rows(r -> r.X), pred = stack_rows(r -> r.pred), novelty = stack_rows(r -> r.novelty),
            inv = stack_rows(r -> r.bracket.inv), ub_ev = stack_rows(r -> r.bracket.ub_ev), lb_ev = stack_rows(r -> r.bracket.lb_ev),
            mid_ev = stack_rows(r -> r.bracket.mid_ev), ub_cvar = stack_rows(r -> r.bracket.ub_cvar), lb_cvar = stack_rows(r -> r.bracket.lb_cvar),
            mid_cvar = stack_rows(r -> r.bracket.mid_cvar))
end

# Greedy max-min with a per-cell quota (cell => picks still allowed). Distances count against
# every solved point (ctx.X) and every earlier pick, in any cell. Returns (row, distance at
# pick) in pick order: every prefix is itself a well-spread design.
function select_gridfill(cands, ctx::GridFillContext, quota::AbstractDict, n::Int; score = novelty_score)
    cands === nothing && return Tuple{Int,Float64}[]
    quota = copy(quota)
    Y = cands.X ./ ctx.scale'
    dmin = nearest_distance(cands.X, ctx.X, ctx.scale)
    picks = Tuple{Int,Float64}[]
    while length(picks) < n
        s = score(dmin, cands, ctx)
        s = [get(quota, cands.cell[i], 0) > 0 ? s[i] : -Inf for i in eachindex(s)]
        i = argmax(s)
        isfinite(s[i]) || break
        push!(picks, (i, dmin[i]))
        quota[cands.cell[i]] -= 1
        dmin .= min.(dmin, vec(pairwise_distances(Y, Y[i:i, :])))
    end
    return picks
end
