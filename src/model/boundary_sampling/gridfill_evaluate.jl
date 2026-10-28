using LinearAlgebra, Random, Statistics, DataFrames, Distributed

#==========================================================================================

GridFill: batched evaluation with feedback (the whole phase)

run_gridfill_sampling! loops: candidates for every open cell → greedy max-min selection →
evaluate one batch of subproblems → feed the solves back (new parents, new cuts, new
occupancy) → re-select on what is still open. Every Summary label (op costs, prices,
emissions) comes from the subproblem solves; the predictions are lineage columns only.

==========================================================================================#

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
    pred::Vector{Float64}          # (ĉ, v̂) of the midpoints
end

function add_gridfill_columns!(df_syscost::DataFrame, pick::GridFillPick, g::Real)
    df_syscost[!, :Parents] = [join(pick.parent_names, ";")]
    df_syscost[!, :Weights] = [join(string.(round.(pick.weights; digits = 4)), ";")]
    df_syscost[!, :GF_Batch] = [pick.batch]
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

function make_gridfill_picks(cands, sel, solves::AbstractVector, batch::Int, first_index::Int)
    return [GridFillPick("GridFill_$(first_index + k - 1)", batch, cands.cell[i], cands.aim[i, :], [solves[j].name for j in cands.tri[i, :]],
                         cands.W[i, :], cands.X[i, :], dist, cands.inv[i], cands.ub_ev[i], cands.lb_ev[i], cands.mid_ev[i],
                         cands.ub_cvar[i], cands.lb_cvar[i], cands.mid_cvar[i], cands.pred[i, :]) for (k, (i, dist)) in enumerate(sel)]
end

# One batch through run_subproblems_batch (each worker loops over the batch for its own
# scenarios). Rows are built exactly as evaluate_picks does, so Summary columns match.
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
        fs.region = "gridfill_sample"
        ev, cvar = solve_ev(fs), solve_cvar(fs, α)
        df_cap, df_syscost, df_emissions = make_results_mapping_dfs(pk.x[1:R], pk.x[R+1:end], sp, inv, inv_by_zone, cvar, ev, inputs, settings)
        add_gridfill_columns!(df_syscost, pk, g)
        push!(rows_cap, df_cap); push!(rows_syscost, df_syscost); push!(rows_emissions, df_emissions)
        push!(solved, fs)
    end
    return rows_cap, rows_syscost, rows_emissions, solved
end

# Prediction quality of one batch: midpoint vs solved and raw blend (UB) vs solved, relative
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
    return (mid_ev = mrel([pk.mid_ev for pk in picks], ev), blend_ev = mrel([pk.ub_ev for pk in picks], ev),
            mid_cvar = mrel([pk.mid_cvar for pk in picks], cv), blend_cvar = mrel([pk.ub_cvar for pk in picks], cv),
            mid_cells = median(vec(maximum(abs.(mid_pos .- pos); dims = 2))) / g,
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
# Returns rows_cap, rows_syscost, rows_emissions, labels.
function run_gridfill_sampling!(bank::AbstractVector, cut_pool, MP, SPs, inputs::Dict, settings::Dict, sampling::AbstractDict;
                                checkpoint = (c, s, e, l) -> nothing, results_folder = nothing,
                                tuning = GRIDFILL_TUNING, score = novelty_score,
                                evaluate = evaluate_gridfill_batch, parallel::Bool = nworkers() > 0 && workers() != [myid()])
    gcfg = sampling["GridFill"]
    g, δcap = Float64(gcfg["Cell size"]), Float64(gcfg["Delta cap"])
    # "Per cell" = "auto": fill in layers (every reachable target cell to 1 sample, then 2, ...)
    # until "Samples" is spent; the level only rises once a layer can't take more picks.
    auto = gcfg["Per cell"] == "auto"
    per_cell = auto ? 1 : Int(gcfg["Per cell"])
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
    @info("GridFill: $(length(cells)) cells of size $g in the near-optimal region (δ <= $δcap), $(length(cells) - length(targets)) already occupied, $(length(targets)) targets; $(auto ? "samples per cell raised in layers until" : "up to $per_cell samples per cell,") $n_max in total, batches of $batch_size.")

    rows_cap, rows_syscost, rows_emissions, labels = DataFrame[], DataFrame[], DataFrame[], String[]
    last_results = GridFillCellResult[]
    rnd = 0
    progress = false   # picks made since per_cell was last raised
    while length(labels) < n_max && !isempty(targets)
        rnd += 1
        solves = solved()
        counts = cell_counts(cost_coords(solves, anchors, α), g)
        open_cells = [c for c in targets if get(counts, c, 0) < per_cell]
        if isempty(open_cells)
            auto || (@info("GridFill: every target cell holds $per_cell samples."); break)
            per_cell += 1; progress = false
            @info("GridFill: every target cell holds $(per_cell - 1) samples; raising to $per_cell per cell.")
            continue
        end

        X = reduce(vcat, [fs.x' for fs in solves])
        Xref = hcat(permutedims(X), random_parent_blends(X, tuning.prune_blends, Xoshiro(hash((seed, rnd)))))
        ctx = GridFillContext(solves, prune_cuts(cuts, Xref), anchors, α, g)
        t_cand = @elapsed last_results = run_candidate_stage(ctx, open_cells, tuning, n_keep, seed, rnd; parallel = parallel)
        cands = pool_candidates(last_results)
        quota = Dict(c => per_cell - get(counts, c, 0) for c in open_cells)
        sel = select_gridfill(cands, ctx, quota, min(batch_size, n_max - length(labels)); score = score)
        if isempty(sel)
            if auto && progress
                per_cell += 1; progress = false
                @info("GridFill round $rnd: the $(length(open_cells)) open cells can't be reached by a blend; raising to $per_cell per cell.")
                continue
            end
            @info("GridFill round $rnd: no open cell can be reached by a blend; stopping.")
            break
        end
        progress = true
        picks = make_gridfill_picks(cands, sel, solves, rnd, length(labels) + 1)

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
        reached = count(c -> get(counts, c, 0) > 0, targets)
        full = count(c -> get(counts, c, 0) >= per_cell, targets)
        n_status(sym) = count(r -> r.status == sym, last_results)
        @info("GridFill batch $rnd: $(length(picks)) samples (candidates $(round(t_cand; digits = 1)) s for $(length(open_cells)) open cells, evaluation $(round(t_eval; digits = 1)) s). " *
              "Median |error| midpoint vs solved: EV $(round(100 * st.mid_ev; sigdigits = 3))%, CVaR $(round(100 * st.mid_cvar; sigdigits = 3))%, position $(round(st.mid_cells; sigdigits = 3)) cells; " *
              "raw blend: EV $(round(100 * st.blend_ev; sigdigits = 3))%, CVaR $(round(100 * st.blend_cvar; sigdigits = 3))%, position $(round(st.blend_cells; sigdigits = 3)) cells. " *
              "Landed in target cell: $(st.landed)/$(length(picks)); bracket violations: $(st.violations); failed solves: $(st.failed). " *
              "Target cells reached $reached/$(length(targets)), full ($per_cell per cell) $full, still empty $(length(targets) - reached) " *
              "(this round: $(n_status(:outside_hull)) outside the parents' cost hull, $(n_status(:no_blend)) without a valid blend, $(n_status(:aim_missed)) missed).")
        try
            checkpoint(rows_cap, rows_syscost, rows_emissions, labels)
        catch err
            @warn("GridFill summary checkpoint failed", exception = (err, catch_backtrace()))
        end
    end

    unreachable = [r.cell for r in last_results if r.status in (:outside_hull, :no_blend)]
    if !isempty(unreachable)
        cs, vs = first.(unreachable) .* g, last.(unreachable) .* g
        @warn("GridFill: $(length(unreachable)) target cells could not be filled by blending ($(count(r -> r.status == :outside_hull, last_results)) outside the parents' cost hull), " *
              "ĉ in [$(minimum(cs)), $(maximum(cs) + g)], v̂ in [$(minimum(vs)), $(maximum(vs) + g)]: " * join(["($(round(c; digits = 4)), $(round(v; digits = 4)))" for (c, v) in zip(cs, vs)], " "))
    end
    @info("GridFill: $(length(labels)) samples evaluated in $rnd rounds; cut bank holds $(n_cuts(cuts)) cuts.")
    return rows_cap, rows_syscost, rows_emissions, labels
end
