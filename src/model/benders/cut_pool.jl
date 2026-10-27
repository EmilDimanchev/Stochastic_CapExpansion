using SparseArrays

export CutPool, seed!, separate!, prune!, add_cuts!, trim_cuts!, count_live_cuts

#=====================================

Cut pool

One object owns every Benders cut (optimality_cut_* and cvar_tail_cuts_*) built against a
shared master problem, across every solve that reuses it (the risk-weight anchors and every
MGA vector's budget/cap solves). Each cut is either
  - live:     a constraint in MP right now (`live`, name => ConstraintRef), or
  - archived: removed from MP and stored as one column of a sparse matrix (below).

The Benders loops use three operations:
  seed!     - start of each solve: choose which cuts are live at the first MP solve
  separate! - every MP solve: solve, then bring back archived cuts the solution violates
  prune!    - after every MP solve: archive live cuts that haven't bound for `window` iterations

Archive layout. Every cut is a scalar `a'v >= b` over MP variables v (x, x_line, alpha, u, ζ).
Column j of the matrix Aᵀ holds the a of archived cut j over the pool's variable columns
(`vars`), and b[j] its RHS - an exact copy of the constraint, so restoring it is lossless.
The matrix is stored as raw CSC vectors (colptr/rowval/nzval) so archiving a cut is a cheap
append; restored cuts are only flagged dead in `alive` and are compacted away later. This
lets separate! check every archived cut at once with a single Aᵀ'v product instead of a
per-cut JuMP evaluation.

Groups. `group` is the MGA vector index (0 for the anchor sweep). `bound_cur` holds cuts
that have bound at some point during the current group, `bound_prev` the previous group's.
seed! starts a solve from bound_cur (later solves of the same group) or bound_prev (first
solve of a new group); separate! recovers anything else the solve needs.

Settings come from the "Cut settings" block of settings.yaml (see CutPool below).

=====================================#

is_cut_name(n::AbstractString) = startswith(n, "optimality_cut_") || startswith(n, "cvar_tail_cuts_")

mutable struct CutPool
    cfg::Dict{String, Any}                  # settings["Cut settings"]
    live::Dict{String, ConstraintRef}       # cuts currently in MP

    # Archive: column j is cut names[j], a'v >= b[j], with a in colptr/rowval/nzval
    vars::Vector{VariableRef}               # matrix row -> MP variable
    col::Dict{VariableRef, Int}             # MP variable -> matrix row
    colptr::Vector{Int}
    rowval::Vector{Int}
    nzval::Vector{Float64}
    b::Vector{Float64}
    names::Vector{String}
    alive::BitVector                        # false once a column has been restored or trimmed
    colof::Dict{String, Int}                # archived (alive) cut name -> column

    # Trim bookkeeping: name => (created, last_bound) on a monotone clock bumped once per
    # prune! call, so stamps are comparable across solves (see trim_cuts!)
    stats::Dict{String, Tuple{Int, Int}}
    clock::Int

    last_bound_iter::Dict{String, Int}      # per solve: iteration each live cut last bound (or went live)

    group::Int
    bound_cur::Set{String}
    bound_prev::Set{String}
end

# Hard-indexes every key so a missing setting fails loudly rather than silently defaulting.
function CutPool(settings::Dict)
    cfg = settings["Cut settings"]
    for key in ("Cuts retained", "Trim rule", "Prune window", "Binding tolerance", "Separation tolerance", "Cuts per round", "Max rounds")
        haskey(cfg, key) || error("Cut settings is missing \"$key\"")
    end
    return CutPool(cfg, Dict{String, ConstraintRef}(),
                   VariableRef[], Dict{VariableRef, Int}(), [1], Int[], Float64[], Float64[], String[], falses(0), Dict{String, Int}(),
                   Dict{String, Tuple{Int, Int}}(), 0, Dict{String, Int}(),
                   -1, Set{String}(), Set{String}())
end

prune_window(pool::CutPool, solve_type::String) = pool.cfg["Prune window"][solve_type]
n_archived(pool::CutPool) = length(pool.colof)

# Number of cuts currently live in MP (scans MP, so it also works without a pool).
function count_live_cuts(MP::Model)
    return count(con -> is_cut_name(name(con)), all_constraints(MP, include_variable_in_set_constraints = false))
end

# Register cuts just created by add_optimality_cuts! as live.
add_cuts!(pool::CutPool, refs::AbstractDict) = (merge!(pool.live, refs); pool)

#=====================================

Archive / restore (private)

=====================================#

function _row!(pool::CutPool, v::VariableRef)
    return get!(pool.col, v) do
        push!(pool.vars, v)
        length(pool.vars)
    end
end

# Removes live cuts from MP and appends each one's exact constraint (coefficients + RHS) as a
# column of the archive, so it can be rebuilt losslessly later by _restore!.
#
# The delete is batched (one delete(MP, ::Vector{ConstraintRef}) call after the loop) rather
# than called once per cut. Both HiGHS.jl and Gurobi.jl's MOI wrappers implement a scalar
# delete that reindexes every remaining constraint in MP (an O(total constraints) walk) on
# every call - and Gurobi's scalar path also forces a GRBupdatemodel sync per call - so
# deleting N cuts one at a time is O(N * total_constraints). The vector delete does that
# reindex/sync once for the whole batch.
#
# Cuts are removed outright rather than masked with a huge RHS: -1e6 rows next to ~1e5-scale
# rows is the bad row-bound ratio that drove level-set regularization to numerical infeasibility.
function _archive!(pool::CutPool, MP::Model, cut_names)
    cons = ConstraintRef[]
    for n in cut_names
        con = pop!(pool.live, n, nothing)
        con === nothing && continue
        obj = constraint_object(con)
        # Sort by row so each column is a valid CSC column (rows ascending)
        entries = sort!([(_row!(pool, v), c) for (v, c) in obj.func.terms]; by = first)
        for (r, c) in entries
            push!(pool.rowval, r)
            push!(pool.nzval, c)
        end
        push!(pool.colptr, length(pool.rowval) + 1)
        push!(pool.b, MOI.constant(obj.set) - obj.func.constant)
        push!(pool.names, n)
        push!(pool.alive, true)
        pool.colof[n] = length(pool.names)
        push!(cons, con)
    end
    isempty(cons) || delete(MP, cons)
    return length(cons)
end

# Rebuilds archived columns as constraints in MP under their original names and marks the
# columns dead. Compacts the archive once dead columns are the majority, so the separation
# product never does more than ~2x the necessary work.
function _restore!(pool::CutPool, MP::Model, cols)
    for j in cols
        pool.alive[j] || continue
        expr = AffExpr(0.0)
        for k in pool.colptr[j]:pool.colptr[j+1]-1
            add_to_expression!(expr, pool.nzval[k], pool.vars[pool.rowval[k]])
        end
        n = pool.names[j]
        pool.live[n] = add_constraint(MP, ScalarConstraint(expr, MOI.GreaterThan(pool.b[j])), n)
        pool.alive[j] = false
        delete!(pool.colof, n)
    end
    2 * n_archived(pool) < length(pool.names) && _compact!(pool)
    return pool
end

# Drops dead columns from the archive.
function _compact!(pool::CutPool)
    keep = findall(pool.alive)
    colptr, rowval, nzval = [1], Int[], Float64[]
    for j in keep
        rng = pool.colptr[j]:pool.colptr[j+1]-1
        append!(rowval, view(pool.rowval, rng))
        append!(nzval, view(pool.nzval, rng))
        push!(colptr, length(rowval) + 1)
    end
    pool.colptr, pool.rowval, pool.nzval = colptr, rowval, nzval
    pool.b, pool.names = pool.b[keep], pool.names[keep]
    pool.alive = trues(length(keep))
    pool.colof = Dict(n => j for (j, n) in enumerate(pool.names))
    return pool
end

# Relative violation of every archived cut at the variable values v (dead columns -> -Inf).
function archive_violations(pool::CutPool, v::AbstractVector{Float64})
    At = SparseMatrixCSC(length(pool.vars), length(pool.b), pool.colptr, pool.rowval, pool.nzval)
    viol = (pool.b .- At' * v) ./ max.(1.0, abs.(pool.b))
    viol[.!pool.alive] .= -Inf
    return viol
end

# Columns of the (at most n) most violated archived cuts with relative violation above tol,
# most violated first.
function most_violated(pool::CutPool, v::AbstractVector{Float64}, tol::Real, n::Int)
    viol = archive_violations(pool, v)
    violated = findall(>(tol), viol)
    return violated[partialsortperm(viol[violated], 1:min(n, length(violated)); rev = true)]
end

#=====================================

seed! / separate! / prune!

=====================================#

# Start of a solve. `group` is the MGA vector index (0 for the anchor sweep). On a new group
# the binding sets roll over; the solve then starts with exactly the cuts that bound during
# this group so far, or - for the first solve of a group - during the previous group.
# Everything else is archived and left for separate! to bring back if needed. Resets the
# per-solve grace bookkeeping, so every live cut gets a full prune window again.
function seed!(pool::CutPool, MP::Model, group::Int)
    new_group = group != pool.group
    if new_group
        pool.bound_prev, pool.bound_cur, pool.group = pool.bound_cur, Set{String}(), group
    end
    target = isempty(pool.bound_cur) ? pool.bound_prev : pool.bound_cur
    _archive!(pool, MP, [n for n in keys(pool.live) if !(n in target)])
    _restore!(pool, MP, [pool.colof[n] for n in target if haskey(pool.colof, n)])
    empty!(pool.last_bound_iter)
    @info("Cut pool seed ($(new_group ? "new group $group" : "same group $group")): $(length(pool.live)) live cuts at first MP solve, $(n_archived(pool)) archived")
    return length(pool.live)
end

# Solves MP, then checks every archived cut at the solution and restores those violated by
# more than "Separation tolerance" (relative to max(1,|b|)), most violated first, at most
# "Cuts per round" per round, re-solving after each round until nothing is violated or
# "Max rounds" rounds have run. The sub-problems should only see the point this returns.
# Returns the planning-model output of the last solve and a NamedTuple of separation stats.
function separate!(pool::CutPool, MP::Model, settings::Dict, risk_aversion_weight)
    tol, per_round, max_rounds = pool.cfg["Separation tolerance"], pool.cfg["Cuts per round"], pool.cfg["Max rounds"]
    output = run_planning_model(MP, settings, risk_aversion_weight)
    restored, rounds, t_check, t_extra_mp = 0, 0, 0.0, 0.0
    while rounds < max_rounds && n_archived(pool) > 0
        t = time()
        violated = most_violated(pool, value.(pool.vars), tol, per_round)
        t_check += time() - t
        isempty(violated) && break
        _restore!(pool, MP, violated)
        restored += length(violated)
        rounds += 1
        t_extra_mp += @elapsed output = run_planning_model(MP, settings, risk_aversion_weight)
    end
    if restored > 0
        @info("Cut separation: restored $restored cuts in $rounds round(s), check $(round(t_check; digits=3)) s, extra MP $(round(t_extra_mp; digits=3)) s")
    end
    return output, (restored = restored, rounds = rounds, t_check = t_check, t_extra_mp = t_extra_mp)
end

function _stamp!(pool::CutPool, n::String; bound::Bool = false)
    created, last_bound = get(pool.stats, n, (pool.clock, 0))
    pool.stats[n] = (created, bound ? pool.clock : last_bound)
end

# After an MP solve (MP must still hold that solve's duals). A live cut binds when its dual
# is at least "Binding tolerance"; binding cuts are stamped (for trim_cuts!) and recorded in
# this group's binding set. A cut first seen this solve (new, or brought back by seed!/
# separate!) starts its grace period here; one that hasn't bound for `window` iterations is
# archived.
function prune!(pool::CutPool, MP::Model, iteration::Int, window::Int)
    pool.clock += 1
    tol = pool.cfg["Binding tolerance"]
    to_archive = String[]
    for (n, con) in pool.live
        is_new = !haskey(pool.last_bound_iter, n)
        if is_new
            pool.last_bound_iter[n] = iteration
            _stamp!(pool, n)
        end
        if dual(con) >= tol
            pool.last_bound_iter[n] = iteration
            _stamp!(pool, n; bound = true)
            push!(pool.bound_cur, n)
        elseif !is_new && iteration - pool.last_bound_iter[n] >= window
            push!(to_archive, n)
        end
    end
    return _archive!(pool, MP, to_archive)
end

# Trims the combined pool - live cuts *and* archived ones - down to at most "Cuts retained",
# permanently discarding the rest (their stats go too). This only bounds memory: seed! no
# longer reactivates the whole archive, so the archive size doesn't set the MP size. Rule
# ("Trim rule"):
#   "recent_binding" - keep the cuts that bound most recently (ties: later-created first)
#   "earliest"       - keep the earliest-created cuts
# Cuts with no stats yet count as newest, so they survive. Always logs, including when it's a
# no-op, so a trim that never engages is visible. Returns (live cuts deleted, archived deleted).
function trim_cuts!(MP::Model, pool::CutPool)
    n_keep, rule = pool.cfg["Cuts retained"], pool.cfg["Trim rule"]
    n_live, n_arch = length(pool.live), n_archived(pool)
    if n_live + n_arch <= n_keep
        @info("Cut trim ($rule): $n_live live + $n_arch archived = $(n_live + n_arch) cuts <= Cuts retained ($n_keep); nothing trimmed")
        return (0, 0)
    end
    names = vcat(collect(keys(pool.live)), collect(keys(pool.colof)))
    stat(n) = get(pool.stats, n, (pool.clock, pool.clock))
    if rule == "recent_binding"
        sort!(names; by = n -> (stat(n)[2], stat(n)[1]), rev = true)
    elseif rule == "earliest"
        sort!(names; by = n -> stat(n)[1])
    else
        error("Unknown cut trim rule \"$rule\". Choose \"recent_binding\" or \"earliest\".")
    end
    live_to_delete = ConstraintRef[]
    n_arch_deleted = 0
    for n in names[n_keep+1:end]
        if haskey(pool.live, n)
            push!(live_to_delete, pop!(pool.live, n))
        else
            pool.alive[pop!(pool.colof, n)] = false
            n_arch_deleted += 1
        end
        delete!(pool.stats, n)
        delete!(pool.bound_cur, n)
        delete!(pool.bound_prev, n)
    end
    isempty(live_to_delete) || delete(MP, live_to_delete) # one batched delete, see _archive!
    _compact!(pool)
    @info("Cut trim ($rule): $n_live live + $n_arch archived -> kept $n_keep; deleted $(length(live_to_delete)) live and $n_arch_deleted archived cuts")
    return (length(live_to_delete), n_arch_deleted)
end
