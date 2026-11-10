using LinearAlgebra, Random, Statistics, JuMP

#==========================================================================================

GridFill: the per-subproblem cut bank and the op-cost bracket of a blend

Every subproblem solve at a capacity x_j gives one supporting hyperplane per subproblem
(scenario s × time block b, see subproblem_blocks),
    θ_{s,b}(x) >= θ_{s,b}(x_j) + β_{s,b}ᵀ (x − x_j) = a + gᵀx,
valid everywhere because each block's op cost is convex in capacity. A scenario's cost is the
sum of its blocks', so its floor is the sum over blocks of each block's best cut. The CutBank
holds the cuts per subproblem, in groups g = s + (b−1)·S, from two sources:
  - the CutPool (live + archived optimality_cut_* rows of the master): every Benders
    iterate of every anchor / MGA / wing solve, not just the final ones;
  - FinalSolves: the bank's final solves and every evaluated GridFill sample.
All values are unscaled ($ and $/MW), in the [x; x_line] order of capacity_variable_names.

Bracket for a blend x = Σ w_j x_j of solved parents (θ_j their per-scenario costs):
  upper: θ_s(x) <= Σ w_j θ_{j,s} =: θ̄_s                  (convexity)
         EV <= p·θ̄,  CVaR <= CVaR(θ̄)  (tighter than Σ w_j CVaR_j)
  lower: θ_s(x) >= l_s := Σ_b max over the cuts of block b of scenario s
         EV >= p·l,  CVaR >= CVaR(l)  (CVaR is monotone in the scenario costs)
The aiming predictor is LB + f·(UB − LB) ("Aim fraction" f, see predicted_position); it is
never used as a label. The midpoint (f = 0.5) is kept as a lineage column.

==========================================================================================#

mutable struct CutBank
    a::Vector{Vector{Float64}}     # per group (scenario × block): cut intercepts, $
    G::Vector{Matrix{Float64}}     # per group: d × n_g cut gradients, $/MW
    n_blocks::Int                  # B; there are S·B groups, g = s + (b−1)·S
end

CutBank(a::Vector{Vector{Float64}}, G::Vector{<:AbstractMatrix}) = CutBank(a, Vector{Matrix{Float64}}(G), 1)
CutBank(S::Int, d::Int; n_blocks::Int = 1) = CutBank([Float64[] for _ in 1:S*n_blocks], [zeros(d, 0) for _ in 1:S*n_blocks], n_blocks)

n_cuts(cb::CutBank) = sum(length, cb.a)
n_scenarios(cb::CutBank) = length(cb.a) ÷ cb.n_blocks

# One cut per subproblem from a solve at fs.x (only solves with gradients and all scenarios optimal)
function add_final_solve!(cb::CutBank, fs)
    (fs.ok && has_gradients(fs)) || return cb
    length(fs.θb) == length(cb.a) || error("add_final_solve!: $(fs.name) has $(length(fs.θb)) subproblems, the cut bank $(length(cb.a)) groups.")
    for g in eachindex(cb.a)
        grad = fs.βb[:, g]
        push!(cb.a[g], fs.θb[g] - dot(grad, fs.x))
        cb.G[g] = hcat(cb.G[g], grad)
    end
    return cb
end

# Every optimality cut in the pool, live or archived, unscaled and in pool order. The master row is
#     c·alpha[s,f,k,b] + Σ a_r x_r + Σ a_l x_line_l >= b      (scaled by "Scaling factor cost")
# with c = 1 for an optimal subproblem (c = 0 rows carry no information and are skipped).
# alpha is the scaled op cost of block b of scenario (s,f,k), so θ_{s,b}(x) >= scale·(b − aᵀx)/c.
# Returns the cuts as (group, a, g, name), group = s + (b−1)·S over (s,f,k,b) column-major (the
# FinalSolve θb order), and the number of optimality rows skipped.
function pool_optimality_cuts(pool, MP::Model, settings::Dict)
    scale = settings["Scaling factor cost"]
    alpha = MP[:alpha]
    xvars = vcat(vec(collect(MP[:x])), vec(collect(MP[:x_line])))
    scen_of = Dict(alpha[I] => li for (li, I) in enumerate(CartesianIndices(size(alpha))))   # (s,f,k,b) column-major, same as θb
    pos_of = Dict(v => i for (i, v) in enumerate(xvars))
    d = length(xvars)
    cuts = NamedTuple{(:group, :a, :g, :name), Tuple{Int, Float64, Vector{Float64}, String}}[]
    skipped = 0
    function add_row!(terms, b, cut_name)
        s, c = 0, 0.0
        g = zeros(d)
        for (v, coef) in terms
            if haskey(scen_of, v)
                s, c = scen_of[v], coef
            elseif haskey(pos_of, v)
                g[pos_of[v]] = coef
            end
        end
        if s == 0 || c <= 0
            skipped += 1
            return
        end
        push!(cuts, (group = s, a = scale * b / c, g = -scale .* g ./ c, name = cut_name))
    end
    for (n, con) in pool.live
        startswith(n, "optimality_cut_") || continue
        obj = constraint_object(con)
        add_row!(obj.func.terms, MOI.constant(obj.set) - obj.func.constant, n)
    end
    for j in eachindex(pool.names)
        (pool.alive[j] && startswith(pool.names[j], "optimality_cut_")) || continue
        rng = pool.colptr[j]:(pool.colptr[j + 1] - 1)
        add_row!(((pool.vars[pool.rowval[k]], pool.nzval[k]) for k in rng), pool.b[j], pool.names[j])
    end
    return cuts, skipped
end

# Adds every optimality cut in the pool to the bank (see pool_optimality_cuts).
function add_pool_cuts!(cb::CutBank, pool, MP::Model, settings::Dict)
    cuts, _ = pool_optimality_cuts(pool, MP, settings)
    a_new = [Float64[] for _ in eachindex(cb.a)]
    g_new = [Vector{Float64}[] for _ in eachindex(cb.a)]
    length(cb.a) == length(MP[:alpha]) || error("add_pool_cuts!: the cut bank has $(length(cb.a)) groups, the master $(length(MP[:alpha])) subproblems.")
    for cut in cuts
        push!(a_new[cut.group], cut.a)
        push!(g_new[cut.group], cut.g)
    end
    for s in eachindex(cb.a)
        isempty(a_new[s]) && continue
        append!(cb.a[s], a_new[s])
        cb.G[s] = hcat(cb.G[s], reduce(hcat, g_new[s]))
    end
    @info("GridFill cut bank: $(length(cuts)) optimality cuts taken from the cut pool ($(length(pool.live)) live, $(n_archived(pool)) archived cuts in the pool)")
    return cb
end

# Per-group cut maximum at each column of X (d × m): (S·B) × m. −Inf where a group has no cut.
function group_cut_values(cb::CutBank, X::AbstractMatrix)
    l = fill(-Inf, length(cb.a), size(X, 2))
    for g in eachindex(cb.a)
        isempty(cb.a[g]) && continue
        l[g, :] = vec(maximum(cb.a[g] .+ cb.G[g]' * X; dims = 1))
    end
    return l
end

# Per-scenario floor at each column of X: S × m, the sum over a scenario's blocks of each block's
# best cut. −Inf where any block of the scenario has no cut.
function cut_values(cb::CutBank, X::AbstractMatrix)
    l = group_cut_values(cb, X)
    S = n_scenarios(cb)
    return dropdims(sum(reshape(l, S, cb.n_blocks, size(X, 2)); dims = 2); dims = 2)
end

# Keep, per group (scenario × block), only the cuts that attain the maximum at one of the reference
# points (columns of Xref). Any subset of cuts still gives a valid lower bound; keeping those active
# around the parents keeps the bound almost as tight at a fraction of the cost.
function prune_cuts(cb::CutBank, Xref::AbstractMatrix)
    out = CutBank(n_scenarios(cb), size(Xref, 1); n_blocks = cb.n_blocks)
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
