#=========

GridFill on the real model (src/model/sampling/gridfill_*.jl), without a Benders run:
  - add_pool_cuts! reads the master's optimality cuts (live and archived, as add_optimality_cuts!
    writes them) into exactly the cuts add_final_solve! builds from the same subproblem solves;
  - evaluate_gridfill_batch gives Summary rows with the per-scenario columns and the GF_ lineage
    columns, solved labels inside the bracket, and write_exploration_results! writes them.
Gas price uncertainty off (10 scenarios).

==========#

function test_gridfill_real_model(test_index = "gridfill_real")
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")
    settings = load_settings(inputs_folder)
    settings["Gas price uncertainty"] = false
    settings["Risk aversion flag"] = true          # also writes cvar_tail_cuts_*, which must be skipped
    inputs = load_input_data(inputs_folder, settings)
    configure_parallel_workers!(settings)
    SPs = build_all_subproblems(inputs, settings)
    R, L = length(inputs["Resources"]), inputs["Number of lines"]
    S, F, K = size(SPs)
    p, α = scenario_probabilities(inputs), settings["Value-at-Risk percent"]

    MP = build_planning_model(inputs, settings; risk_aversion_weight = 0.5)
    pool = CutPool(settings)
    rng = Xoshiro(4)
    ub = inputs["Capacity upper bounds"]
    solves = FinalSolve[]
    for j in 1:4
        cap, line = 0.3 .* ub .* rand(rng, R), 200 .* rand(rng, L)
        sp = run_all_subproblems(SPs, inputs, settings, cap, line; minimal_payload = true)
        obj = reshape([sp[s,f,k]["SP objective"] for s in 1:S, f in 1:F, k in 1:K], (S, F, K))
        cd = reshape([sp[s,f,k]["capacity dual"][r] for r in 1:R, s in 1:S, f in 1:F, k in 1:K], (R, S, F, K))
        ld = reshape([sp[s,f,k]["line dual"][l] for l in 1:L, s in 1:S, f in 1:F, k in 1:K], (L, S, F, K))
        co = reshape([sp[s,f,k]["coeff"] for s in 1:S, f in 1:F, k in 1:K], (S, F, K))
        add_cuts!(pool, add_optimality_cuts!(MP, obj, cd, ld, cap, line, co, inputs, settings, j, "gfcheck"))
        x = vcat(cap, line)
        fs = final_solve_from_sp_outputs("Point_$j", x, sp, p, investment_cost(x, inputs)[1])
        push!(solves, fs)
    end
    _archive!(pool, MP, [n for n in keys(pool.live) if occursin("_gfcheck_2_", n)])   # read both paths
    cb_pool = add_pool_cuts!(CutBank(S * F * K, R + L), pool, MP, settings)
    cb_ref = CutBank(S * F * K, R + L)
    foreach(fs -> add_final_solve!(cb_ref, fs), solves)
    @assert n_cuts(cb_pool) == n_cuts(cb_ref) == 4 * S * F * K "pool gave $(n_cuts(cb_pool)) cuts, expected $(4 * S * F * K)"
    Y = hcat([fs.x for fs in solves]..., [0.3 .* vcat(ub, 200) .* rand(rng, R + L) for _ in 1:5]...)
    @assert isapprox(cut_values(cb_pool, Y), cut_values(cb_ref, Y); rtol = 1e-8) "pool cuts differ from the (θ, β, x) cuts"

    # Two blends of the points through the GridFill evaluator
    Θ = reduce(hcat, [fs.θ for fs in solves]); inv = [fs.inv_cost for fs in solves]
    tri = [1 2 3; 2 3 4]; W = [0.5 0.3 0.2; 0.2 0.2 0.6]
    X = blend_capacities(reduce(vcat, [fs.x' for fs in solves]), tri, W)
    br = blend_bracket(cb_pool, Θ, inv, p, α, tri, W, permutedims(X))
    picks = [GridFillPick("GridFill_$i", 1, (i, i), [0.1, 0.2], [solves[j].name for j in tri[i, :]], W[i, :], X[i, :], 0.5,
                          br.inv[i], br.ub_ev[i], br.lb_ev[i], br.mid_ev[i], br.ub_cvar[i], br.lb_cvar[i], br.mid_cvar[i], [0.1, 0.2], 1, 1) for i in 1:2]
    c, s, e, solved = evaluate_gridfill_batch(picks, SPs, inputs, settings, 0.1)
    tol(v) = 1e-6 * abs(v)
    for (pk, fs) in zip(picks, solved)
        ev, cv = solve_ev(fs), solve_cvar(fs, α)
        @assert pk.lb_ev - tol(ev) <= ev <= pk.ub_ev + tol(ev) "$(pk.name): EV $ev outside [$(pk.lb_ev), $(pk.ub_ev)]"
        @assert pk.lb_cvar - tol(cv) <= cv <= pk.ub_cvar + tol(cv) "$(pk.name): CVaR $cv outside [$(pk.lb_cvar), $(pk.ub_cvar)]"
        @assert isapprox(fs.inv_cost, pk.inv; rtol = 1e-9) "investment cost should blend exactly"
        @info("[$(test_index)] $(pk.name): EV $(round(ev; sigdigits = 6)) in [$(round(pk.lb_ev; sigdigits = 6)), $(round(pk.ub_ev; sigdigits = 6))], midpoint error $(round(100 * abs(pk.mid_ev - ev) / ev; sigdigits = 3))% vs blend $(round(100 * abs(pk.ub_ev - ev) / ev; sigdigits = 3))%")
    end
    folder = joinpath("outputs", "Test_" * test_index, "Summary")
    write_exploration_results!(c, s, e, folder, [pk.name for pk in picks], "gridfill_check")
    df = CSV.read(joinpath(folder, "Summary_gridfill_check.csv"), DataFrame)
    @assert df.Index == ["GridFill_1", "GridFill_2"]
    for col in ("EV_OpCost", "CVaR_OpCost", "Parents", "Weights", "GF_Batch", "GF_Layer", "GF_Attempt", "GF_Mid_EV", "GF_LB_CVaR", "GF_Pred_v", "GF_Novelty")
        @assert col in names(df) "Summary is missing $col"
    end
    settings["Write all scenarios flag"] && @assert count(startswith("OperatingCost_D_"), names(df)) == S * F * K "per-scenario op cost columns missing"
    @assert isapprox(df.EV_OpCost, [solve_ev(fs) for fs in solved]; rtol = 1e-9) "Summary labels must be the solved values"
    @info("[$(test_index)] PASSED: pool cuts match, GridFill rows bracketed and written to $(folder).")
    return true
end
