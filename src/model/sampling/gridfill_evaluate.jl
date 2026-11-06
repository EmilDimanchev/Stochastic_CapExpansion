using LinearAlgebra, Random, Statistics, DataFrames, Distributed

#==========================================================================================

GridFill: batched evaluation with feedback (the whole phase)

run_gridfill_sampling! loops: candidates for every open cell → greedy max-min selection →
evaluate one batch of subproblems → feed the solves back (new parents, new cuts, new
occupancy) → re-select on what is still open. Every Summary label (op costs, prices,
emissions) comes from the subproblem solves; the predictions are lineage columns only.

==========================================================================================#

# Valid values of the "Method" key in Sampling settings
const SAMPLING_METHODS = ("None", "GridFill")

struct GridFillPick
    name::String
    batch::Int
    cell::Tuple{Int,Int}
    aim::Vector{Float64}           # (ĉ, v̂) the blend was aimed at
    parent_names::Vector{String}
    weights::Vector{Float64}
    x::Vector{Float64}             # capacities to solve, MW
    novelty::Float64               # range-scaled distance to the nearest solved point / earlier pick, at selection
    inv::Float64                   # exact (linear)
    ub_ev::Float64
    lb_ev::Float64
    mid_ev::Float64
    ub_cvar::Float64
    lb_cvar::Float64
    mid_cvar::Float64
    pred::Vector{Float64}          # (ĉ, v̂) predicted by LB + aim_frac·(UB − LB)
    layer::Int                     # layer the sample was generated in
    attempt::Int                   # attempt number for its target cell
end

function add_gridfill_columns!(df_syscost::DataFrame, pick::GridFillPick, g::Real)
    df_syscost[!, :Parents] = [join(pick.parent_names, ";")]
    df_syscost[!, :Weights] = [join(string.(round.(pick.weights; digits = 4)), ";")]
    df_syscost[!, :GF_Batch] = [pick.batch]
    df_syscost[!, :GF_Layer] = [pick.layer]
    df_syscost[!, :GF_Attempt] = [pick.attempt]
    df_syscost[!, :GF_Cell_c] = [pick.cell[1] * g]       # lower-left corner of the target cell
    df_syscost[!, :GF_Cell_v] = [pick.cell[2] * g]
    df_syscost[!, :GF_Aim_c] = [pick.aim[1]]
    df_syscost[!, :GF_Aim_v] = [pick.aim[2]]
    df_syscost[!, :GF_Novelty] = [pick.novelty]
    df_syscost[!, :GF_UB_EV] = [pick.ub_ev]
    df_syscost[!, :GF_LB_EV] = [pick.lb_ev]
    df_syscost[!, :GF_Mid_EV] = [pick.mid_ev]
    df_syscost[!, :GF_UB_CVaR] = [pick.ub_cvar]
    df_syscost[!, :GF_LB_CVaR] = [pick.lb_cvar]
    df_syscost[!, :GF_Mid_CVaR] = [pick.mid_cvar]
    df_syscost[!, :GF_Pred_c] = [pick.pred[1]]
    df_syscost[!, :GF_Pred_v] = [pick.pred[2]]
    return df_syscost
end

function make_gridfill_picks(cands, sel, solves::AbstractVector, batch::Int, first_index::Int, layer::Int, attempts::AbstractVector{Int})
    return [GridFillPick("GridFill_$(first_index + k - 1)", batch, cands.cell[i], cands.aim[i, :], [solves[j].name for j in cands.tri[i, :]],
                         cands.W[i, :], cands.X[i, :], dist, cands.inv[i], cands.ub_ev[i], cands.lb_ev[i], cands.mid_ev[i],
                         cands.ub_cvar[i], cands.lb_cvar[i], cands.mid_cvar[i], cands.pred[i, :], layer, attempts[k]) for (k, (i, dist)) in enumerate(sel)]
end

# One batch through run_subproblems_batch (each worker loops over the batch for its own
# scenarios). Rows go through make_results_mapping_dfs, so Summary columns match the MGA rows.
function evaluate_gridfill_batch(picks::Vector{GridFillPick}, SPs, inputs::Dict, settings::Dict, g::Real)
    R = length(inputs["Resources"])
    p = scenario_probabilities(inputs)
    α = settings["Value-at-Risk percent"]
    sps = run_subproblems_batch(SPs, inputs, settings, [pk.x[1:R] for pk in picks], [pk.x[R+1:end] for pk in picks]; minimal_payload = false)
    rows_cap, rows_syscost, rows_emissions = DataFrame[], DataFrame[], DataFrame[]
    solved = FinalSolve[]
    for (pk, sp) in zip(picks, sps)
        inv, inv_by_zone = investment_cost(pk.x, inputs)
        fs = final_solve_from_sp_outputs(pk.name, pk.x, sp, p, inv)
        ev, cvar = solve_ev(fs), solve_cvar(fs, α)
        df_cap, df_syscost, df_emissions = make_results_mapping_dfs(pk.x[1:R], pk.x[R+1:end], sp, inv, inv_by_zone, cvar, ev, inputs, settings)
        add_gridfill_columns!(df_syscost, pk, g)
        push!(rows_cap, df_cap); push!(rows_syscost, df_syscost); push!(rows_emissions, df_emissions)
        push!(solved, fs)
    end
    return rows_cap, rows_syscost, rows_emissions, solved
end

# Prediction quality of one batch: aimed prediction, midpoint and raw blend (UB) vs solved, relative
# op-cost errors and position errors in cell units, plus bracket violations.
function gridfill_batch_stats(picks::Vector{GridFillPick}, solved::AbstractVector, ctx, rel_tol::Real)   # ctx::GridFillContext (defined in gridfill_sampler.jl, loaded later)
    α, a, g = ctx.α, ctx.anchors, ctx.g
    ev = [solve_ev(fs) for fs in solved]
    cv = [solve_cvar(fs, α) for fs in solved]
    pos = cost_coords(solved, a, α)
    ub_pos = hcat(c_hat(a, [pk.inv + pk.ub_ev for pk in picks]), v_hat(a, [pk.inv + pk.ub_cvar for pk in picks]))
    mid_pos = reduce(vcat, [pk.pred' for pk in picks])
    tol(v) = rel_tol * max(abs(v), 1.0)
    viol = count(i -> ev[i] < picks[i].lb_ev - tol(ev[i]) || ev[i] > picks[i].ub_ev + tol(ev[i]) ||
                      cv[i] < picks[i].lb_cvar - tol(cv[i]) || cv[i] > picks[i].ub_cvar + tol(cv[i]), eachindex(picks))
    mrel(x, y) = median(abs.(x .- y) ./ abs.(y))
    aimed(lb, ub, f) = lb + f * (ub - lb)
    return (pred_ev = mrel([aimed(pk.lb_ev, pk.ub_ev, ctx.aim_frac.ev) for pk in picks], ev), pred_cvar = mrel([aimed(pk.lb_cvar, pk.ub_cvar, ctx.aim_frac.cvar) for pk in picks], cv),
            mid_ev = mrel([pk.mid_ev for pk in picks], ev), blend_ev = mrel([pk.ub_ev for pk in picks], ev),
            mid_cvar = mrel([pk.mid_cvar for pk in picks], cv), blend_cvar = mrel([pk.ub_cvar for pk in picks], cv),
            pred_cells = median(vec(maximum(abs.(mid_pos .- pos); dims = 2))) / g,
            blend_cells = median(vec(maximum(abs.(ub_pos .- pos); dims = 2))) / g,
            landed = count(i -> cell_of(pos[i, :], g) == picks[i].cell, eachindex(picks)),
            violations = viol, failed = count(fs -> !fs.ok, solved))
end

# Random parent blends (Dirichlet(1,1,1) over random triples), MW, d × m: reference points
# where the cut bank should stay tight (with the solved points themselves).
function random_parent_blends(X::AbstractMatrix, m::Int, rng::AbstractRNG)
    n = size(X, 1)
    tri = rand(rng, 1:n, m, 3)
    E = randexp(rng, m, 3)
    return permutedims(blend_capacities(X, tri, E ./ sum(E; dims = 2)))
end

#---- The whole phase -----------------------------------------------------------------------
# bank:      FinalSolves of the exploration phase (anchors, MGA, wings). GridFill samples are
#            appended, so they become parents for later batches.
# cut_pool:  the CutPool shared by every Benders solve on MP (nothing: bank cuts only).
# sampling:  the resolved "Sampling settings" block ("Samples", "Batch size", "Candidates",
#            "Write gradients", and the "GridFill" sub-block).
# checkpoint(rows_cap, rows_syscost, rows_emissions, labels) is called after every batch with
# all GridFill rows so far (the caller rewrites the Summary).
# pick_log: optional Vector that collects every GridFillPick, in order (used by the cut export
# to reproduce the logged bounds).
# Returns rows_cap, rows_syscost, rows_emissions, labels.
function run_gridfill_sampling!(bank::AbstractVector, cut_pool, MP, SPs, inputs::Dict, settings::Dict, sampling::AbstractDict;
                                checkpoint = (c, s, e, l) -> nothing, results_folder = nothing, pick_log = nothing,
                                tuning = GRIDFILL_TUNING, score = novelty_score,
                                evaluate = evaluate_gridfill_batch, parallel::Bool = nworkers() > 0 && workers() != [myid()])
    gcfg = sampling["GridFill"]
    g, δcap = Float64(gcfg["Cell size"]), Float64(gcfg["Delta cap"])
    # "Per cell" = "auto": fill in layers (every reachable target cell to 1 sample, then 2, ...)
    # until "Samples" is spent; the level rises once every open cell is full or exhausted.
    auto = gcfg["Per cell"] == "auto"
    per_cell = auto ? 1 : Int(gcfg["Per cell"])
    max_misses, max_dry = Int(get(gcfg, "Max misses", 3)), Int(get(gcfg, "Max dry rounds", 3))
    aim_frac = aim_fractions(get(gcfg, "Aim fraction", GRIDFILL_DEFAULT_AIM_FRAC))
    aim_margin = Float64(get(gcfg, "Aim margin", GRIDFILL_DEFAULT_AIM_MARGIN))
    frontier_depth, boundary_margin = Float64(get(gcfg, "Frontier depth", 1.0)), Float64(get(gcfg, "Boundary margin", 0.0))
    n_max, batch_size, n_keep = Int(sampling["Samples"]), Int(sampling["Batch size"]), Int(sampling["Candidates"])
    α, seed = settings["Value-at-Risk percent"], settings["Seed"]
    var_names = capacity_variable_names(inputs)

    solved() = [fs for fs in bank if fs.ok]
    solves = solved()
    anchors = gridfill_anchors(solves, α)
    S, d = length(solves[1].θ), length(solves[1].x)

    cuts = CutBank(S, d)
    cut_pool === nothing || add_pool_cuts!(cuts, cut_pool, MP, settings)
    foreach(fs -> add_final_solve!(cuts, fs), solves)

    sweep = cost_coords([fs for fs in solves if startswith(fs.name, "Risk_Weight")], anchors, α)
    cells = grid_cells(sweep, g, δcap)
    counts0 = cell_counts(cost_coords(solves, anchors, α), g)
    targets = [c for c in cells if get(counts0, c, 0) == 0]
    excluded = near_lower_frontier(targets, convex_hull_2d(cost_coords(solves, anchors, α)), g, boundary_margin)
    setdiff!(targets, excluded)
    @info("GridFill: $(length(cells)) cells of size $g in the near-optimal region (δ <= $δcap), $(length(cells) - length(targets) - length(excluded)) already occupied, $(length(excluded)) excluded by gf_boundary_margin = $boundary_margin, $(length(targets)) targets; $(auto ? "samples per cell raised in layers until" : "up to $per_cell samples per cell,") $n_max in total, batches of $batch_size; gf_max_misses = $max_misses, gf_max_dry = $max_dry, gf_aim_frac = EV $(aim_frac.ev), CVaR $(aim_frac.cvar), gf_aim_margin = $aim_margin.")

    rows_cap, rows_syscost, rows_emissions, labels = DataFrame[], DataFrame[], DataFrame[], String[]
    state = GridFillCellState()
    rnd, layer, picks_in_layer = 0, 1, 0
    while length(labels) < n_max && !isempty(targets)
        rnd += 1
        solves = solved()
        counts = cell_counts(cost_coords(solves, anchors, α), g)
        open_cells = [c for c in targets if get(counts, c, 0) < per_cell && !is_exhausted(state, c) && !(c in state.parked)]
        if isempty(open_cells)
            all(c -> is_exhausted(state, c), targets) && (@info("GridFill: every target cell is exhausted; stopping."); break)
            auto || (@info("GridFill: every target cell holds $per_cell samples, is exhausted or has no candidates."); break)
            # Parked cells get retried next layer, so a layer with no picks at all could repeat forever
            (picks_in_layer == 0 && !isempty(state.parked)) && (@info("GridFill: layer $layer produced no samples and its $(length(state.parked)) remaining cells have no candidates; stopping."); break)
            unpark_all!(state)
            per_cell += 1; layer += 1; picks_in_layer = 0
            @info("GridFill: every target cell holds $(per_cell - 1) samples, is exhausted or has no candidates; layer $layer, raising to $per_cell per cell.")
            continue
        end
        n_exhausted0 = length(state.exhausted)

        X = reduce(vcat, [fs.x' for fs in solves])
        Xref = hcat(permutedims(X), random_parent_blends(X, tuning.prune_blends, Xoshiro(hash((seed, rnd)))))
        ctx = GridFillContext(solves, prune_cuts(cuts, Xref), anchors, α, g; aim_frac = [aim_frac.ev, aim_frac.cvar], aim_margin = aim_margin)
        t_cand = @elapsed results = run_candidate_stage(ctx, open_cells, tuning, n_keep, seed, rnd; parallel = parallel)
        mark_outside_hull!(state, results)
        for r in results
            r.status in (:no_blend, :aim_missed) || continue   # no candidate this round: a miss if empty, else a dry round
            record_attempt!(state, r.cell)
            record_no_candidate!(state, r.cell, get(counts, r.cell, 0), max_misses, max_dry)
        end
        cands = pool_candidates(results)
        quota = Dict(c => per_cell - get(counts, c, 0) for c in open_cells)
        sel = select_gridfill(cands, ctx, quota, min(batch_size, n_max - length(labels)); score = score)
        if isempty(sel)
            # Every open cell was just charged a miss or a dry round, so this can't repeat forever
            @info("GridFill round $rnd: none of the $(length(open_cells)) open cells has a candidate; exhausted this round $(length(state.exhausted) - n_exhausted0) (total $(length(state.exhausted))).")
            continue
        end
        attempts = [record_attempt!(state, cands.cell[i]) for (i, _) in sel]
        foreach(((i, _),) -> push!(state.picked, cands.cell[i]), sel)
        picks_in_layer += length(sel)
        picks = make_gridfill_picks(cands, sel, solves, rnd, length(labels) + 1, layer, attempts)
        pick_log === nothing || append!(pick_log, picks)

        t_eval = @elapsed c, s, e, new_solves = evaluate(picks, SPs, inputs, settings, g)
        append!(rows_cap, c); append!(rows_syscost, s); append!(rows_emissions, e)
        append!(labels, [pk.name for pk in picks])
        for fs in new_solves
            push!(bank, fs)
            add_final_solve!(cuts, fs)
            if sampling["Write gradients"] && results_folder !== nothing
                write_gradient_file(joinpath(results_folder, "GridFill", fs.name), fs, var_names)
            end
        end

        st = gridfill_batch_stats(picks, new_solves, ctx, settings["Convergence tolerance"])
        counts = cell_counts(cost_coords(solved(), anchors, α), g)
        pos = cost_coords(new_solves, anchors, α)
        for (i, pk) in enumerate(picks)
            (new_solves[i].ok && cell_of(pos[i, :], g) == pk.cell) && continue
            record_miss!(state, pk.cell, max_misses; open = get(counts, pk.cell, 0) < per_cell)
        end
        reached = count(c -> get(counts, c, 0) > 0, targets)
        full = count(c -> get(counts, c, 0) >= per_cell, targets)
        n_status(sym) = count(r -> r.status == sym, results)
        @info("GridFill batch $rnd: $(length(picks)) samples (candidates $(round(t_cand; digits = 1)) s for $(length(open_cells)) open cells, evaluation $(round(t_eval; digits = 1)) s). " *
              "Median |error| predicted vs solved: EV $(round(100 * st.pred_ev; sigdigits = 3))%, CVaR $(round(100 * st.pred_cvar; sigdigits = 3))%, position $(round(st.pred_cells; sigdigits = 3)) cells; " *
              "midpoint: EV $(round(100 * st.mid_ev; sigdigits = 3))%, CVaR $(round(100 * st.mid_cvar; sigdigits = 3))%; " *
              "raw blend: EV $(round(100 * st.blend_ev; sigdigits = 3))%, CVaR $(round(100 * st.blend_cvar; sigdigits = 3))%, position $(round(st.blend_cells; sigdigits = 3)) cells. " *
              "Landed in target cell: $(st.landed)/$(length(picks)); bracket violations: $(st.violations); failed solves: $(st.failed). " *
              "Layer $layer (quota $per_cell). Target cells reached $reached/$(length(targets)), full ($per_cell per cell) $full, still empty $(length(targets) - reached), " *
              "exhausted this round $(length(state.exhausted) - n_exhausted0) (total $(length(state.exhausted))), parked $(length(state.parked)) " *
              "(this round: $(n_status(:outside_hull)) outside the parents' cost hull, $(n_status(:no_blend)) without a valid blend, $(n_status(:aim_missed)) missed).")
        try
            checkpoint(rows_cap, rows_syscost, rows_emissions, labels)
        catch err
            @warn("GridFill summary checkpoint failed", exception = (err, catch_backtrace()))
        end
    end

    fmt_cells(cs) = join(["($(round(c[1] * g; digits = 4)), $(round(c[2] * g; digits = 4)))" for c in cs], " ")
    out_of_hull, repeated_miss = exhausted_cells(state, :outside_hull), exhausted_cells(state, :missed)
    if !isempty(state.exhausted)
        @warn("GridFill: $(length(state.exhausted)) target cells exhausted (lower-left corners). " *
              "$(length(out_of_hull)) outside the parents' cost hull: $(fmt_cells(out_of_hull)). " *
              "$(length(repeated_miss)) after $max_misses misses: $(fmt_cells(repeated_miss)).")
    end
    final_counts = cell_counts(cost_coords(solved(), anchors, α), g)
    per_target = [get(final_counts, c, 0) for c in targets]
    final_hull = convex_hull_2d(cost_coords(solved(), anchors, α))
    labels_unreached = label_unreached_cells(state, final_counts, final_hull, g, frontier_depth)
    n_reached = count(>(0), per_target)
    n_reachable = length(targets) - length(out_of_hull) - length(labels_unreached.frontier) - length(labels_unreached.out_of_reach)
    @info("GridFill: coverage $n_reached/$(length(targets)) target cells raw, $n_reached/$n_reachable excluding $(length(out_of_hull)) outside the hull, " *
          "$(length(labels_unreached.frontier)) on the lower frontier (depth < $frontier_depth cells): $(fmt_cells(labels_unreached.frontier)), " *
          "and $(length(labels_unreached.out_of_reach)) out of reach near the upper edge: $(fmt_cells(labels_unreached.out_of_reach)).")
    @info("GridFill: final samples per target cell (samples: cells): " * join(["$k: $(count(==(k), per_target))" for k in sort(unique(per_target))], ", ") * ".")
    @info("GridFill: $(length(labels)) samples evaluated in $rnd rounds; cut bank holds $(n_cuts(cuts)) cuts.")
    return rows_cap, rows_syscost, rows_emissions, labels
end
