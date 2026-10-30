#=========

Capacity bound tightening test: build_planning_model previously bounded every resource's
capacity variable by an arbitrary global x_ub = 1e6, unrelated to the input data (every
resource's own real Capacity_UB in Resources.csv/Resources_storage.csv is far smaller,
~1e5-2e5 for the laptop inputs) and needlessly widened the master problem's RHS/bound
coefficient range. build_planning_model now bounds each x[r] by its own
inputs["Capacity upper bounds"][r] instead. This test checks that binding is actually
wired up correctly (JuMP's upper_bound matches the real per-resource data, not 1e6) and
that the master problem still solves normally under a realistic synthetic cut with the
tighter bounds in place.

Note: an earlier version of this fix tried to additionally rescale the x/x_line
variables themselves (dividing their values while multiplying every cost coefficient
that touches them back up to compensate) to bring their magnitude closer to the
already-scaled cost terms. That was reverted: compensating a coefficient by exactly the
factor used to shrink the variable is a wash for any row that also contains other,
differently-scaled terms (e.g. probability-weighted alpha/cvar terms in budget
constraints and cuts) - it does not shrink the row's coefficient range, and in some rows
widens it. It also introduced a real bug (the small "avoid exact zero" lower bound on x
got divided down to 1e-8, which is exactly what HiGHS's "excessively small column
bounds" warning was flagging). Bound-tightening against real data has none of those
failure modes since it never touches a matrix coefficient, only a box bound.

==========#

function test_capacity_bound_tightening(test_index = "cap_bound")
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)

    P_s = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]
    S, F, K = length(P_s), length(P_f), length(P_k)
    G = inputs["Number of generation resources"]
    O = inputs["Number of storage resources"]
    R = G + O
    L = inputs["Number of lines"]
    cost_inv = inputs["Investment costs"]
    x_ub_real = inputs["Capacity upper bounds"]

    MP = build_planning_model(inputs, settings)

    bounds = [normalized_rhs(MP[:max_capacity][r]) for r in 1:R]
    @assert bounds == x_ub_real "max_capacity bound does not match inputs[\"Capacity upper bounds\"]: $(bounds) vs $(x_ub_real)"
    @assert all(bounds .< 1e6) "Expected every per-resource bound to be tighter than the old arbitrary 1e6 constant, got $(bounds)"
    @info("[$(test_index)] Confirmed max_capacity bounds match real per-resource Capacity_UB data (max $(maximum(bounds)), previously a flat 1e6 for every resource).")

    Random.seed!(1234)
    # Synthetic cut shaped like a real Benders cut: capacity duals centered on investment
    # cost (economically meaningful rather than degenerate) plus a baseline dispatch cost,
    # so the master problem actually trades off inv_cost against alpha instead of just
    # collapsing to a bound.
    x_prev = zeros(R)
    x_prev_line = zeros(L)
    cap_dual = (cost_inv .* (0.9 .+ 0.2 .* rand(R))) .* ones(R, S, F, K)
    line_dual = 100 .* rand(L, S, F, K)
    SP_obj = 1e6 .* (0.9 .+ 0.2 .* rand(S, F, K))
    coeffs = fill(1.0, S, F, K)

    add_optimality_cuts!(MP, SP_obj, cap_dual, line_dual, x_prev, x_prev_line, coeffs, inputs, settings, 1, "boundtest")
    output = run_planning_model(MP, settings, 0.5)

    @assert all(output["Capacity"] .<= x_ub_real .+ 1e-6) "Solved capacity exceeds the tightened per-resource bound"
    @info("[$(test_index)] PASSED: master problem solves normally with data-driven per-resource bounds.")

    return (bounds = bounds, capacity = output["Capacity"])
end

# Checks that wing_caps places every cap strictly inside the wing (never at the central
# edge or the lens end). Uses the CVaR-wing values from the test (14) cluster run: the
# central edge (13929.06) and the budget 1 / budget 3 lens ends recovered from that run's caps.
function test_wing_caps(test_index = "wing_caps")
    edge = 13929.05591446445
    span = 20000.0
    min_gap_frac = 0.03
    for (label, tight_end, expected) in (("Budget 1", 7686.983378547586, [11848.365069158828, 9767.674223853207]),
                                         ("Budget 3", 9186.767128728072, [12348.292985885657, 10767.530057306865]))
        caps = wing_caps(tight_end, edge, span, 2, min_gap_frac)
        @assert length(caps) == 2 "$label: expected 2 caps, got $(length(caps))"
        @assert all(isapprox.(caps, expected; rtol = 1e-9)) "$label: caps $caps do not match expected $expected"
        @assert all(tight_end .< caps .< edge) "$label: a cap lies on or outside the wing boundary"
        @assert issorted(caps; rev = true) "$label: caps not ordered loose -> tight"
        @info("[$(test_index)] $label: caps = $(round.(caps; digits=2)), all strictly inside ($(round(tight_end; digits=2)), $(round(edge; digits=2))).")
    end

    # Narrow wing: room for 3 gaps of min_gap_frac*span (=600) means only 2 caps fit even with n_max = 5
    caps = wing_caps(0.0, 1900.0, span, 5, min_gap_frac)
    @assert length(caps) == 2 "narrow wing: expected 2 caps, got $(length(caps))"
    @assert minimum(diff(vcat(1900.0, caps, 0.0)) .* -1) >= min_gap_frac * span - 1e-9 "narrow wing: gap smaller than minimum spacing"

    # Too narrow for any interior cap (fewer than 2 minimum gaps)
    @assert isempty(wing_caps(0.0, 1100.0, span, 2, min_gap_frac)) "too-narrow wing should return no caps"

    @info("[$(test_index)] PASSED: wing caps lie strictly inside the wing and respect the minimum spacing.")
    return true
end

#=========

Settings groups test (experiments/mga/run_settings.jl): every group block is present in both
input sets' settings.yaml, overrides replace only the keys they name (a nested block such as
Cut settings' "Prune window" merges), unknown keys and groups error out, and resolving never
touches the caller's settings.

==========#

function test_settings_groups(test_index = "settings_groups")
    for folder in ("Inputs_30d_1000scen_7tech_2z", "Inputs_30d_1000scen_7tech_2z_Della")
        settings = load_settings(joinpath("inputs", folder))
        @assert all(g -> haskey(settings, g), SETTINGS_GROUPS) "$folder is missing a settings group"
        CutPool(resolve_all_settings(settings)) # every Cut settings key is present
    end
    settings = load_settings(joinpath("inputs", "Inputs_30d_1000scen_7tech_2z"))
    before = deepcopy(settings)

    resolved = resolve_all_settings(settings, Dict("Budget settings" => Dict("Tighten budget" => true),
                                                   "Cut settings" => Dict("Prune window" => Dict("Budget" => 5), "Max rounds" => 9)))
    @assert resolved["Budget settings"]["Tighten budget"] == true
    @assert resolved["Budget settings"]["Budget type"] == settings["Budget settings"]["Budget type"] "non-overridden key changed"
    @assert resolved["Cut settings"]["Max rounds"] == 9
    window = resolved["Cut settings"]["Prune window"]
    @assert window["Budget"] == 5 && window["Anchor"] == settings["Cut settings"]["Prune window"]["Anchor"] "nested block should merge, not replace"
    @assert CutPool(resolved).cfg["Max rounds"] == 9 "the pool should see the resolved cut settings"
    @assert settings == before "resolving modified the caller's settings"

    errors(f) = try f(); false catch e; e isa ErrorException end
    @assert errors(() -> resolve_settings(settings, "Budget settings", Dict("Tighten budgte" => true))) "a misspelled key should error"
    @assert errors(() -> resolve_all_settings(settings, Dict("Budget setings" => Dict()))) "a misspelled group should error"
    @assert errors(() -> resolve_settings(settings, "Missing settings")) "a missing group should error"

    # Sampling settings: method sub-blocks merge, misspelled nested keys and the retired
    # Mapping / Boundary keys error out
    sampling = resolve_settings(settings, "Sampling settings", Dict("Method" => "GridFill", "GridFill" => Dict("Cell size" => 0.1)))
    @assert sampling["GridFill"]["Cell size"] == 0.1 && sampling["GridFill"]["Per cell"] == settings["Sampling settings"]["GridFill"]["Per cell"] "GridFill block should merge"
    @assert sampling["Method"] in SAMPLING_METHODS
    @assert boundary_cfg(sampling)["Candidates"] == sampling["Candidates"] && haskey(boundary_cfg(sampling), "K neighbours")
    @assert errors(() -> resolve_settings(settings, "Sampling settings", Dict("GridFill" => Dict("Cell sise" => 0.1)))) "a misspelled nested key should error"
    @assert errors(() -> resolve_settings(settings, "Mapping settings", Dict("Sampling" => true))) "the retired Mapping \"Sampling\" key should error"
    @assert errors(() -> resolve_settings(settings, "Mapping settings", Dict("Interior samples" => 5))) "the retired \"Interior samples\" key should error"
    @assert errors(() -> resolve_all_settings(settings, Dict("Boundary settings" => Dict()))) "the retired Boundary settings group should error"
    @info("[$(test_index)] settings groups resolve and validate as expected for both input sets.")
    return true
end

#=========

Batched subproblem test (src/model/benders/worker_functions.jl): run_subproblems_batch solves
several capacity points with one call per worker over its own scenarios. Its objectives and
duals must match run_all_subproblems point by point (gas price uncertainty off: 10 scenarios).

==========#

function test_subproblems_batch(test_index = "sp_batch")
    settings = load_settings(joinpath("inputs", "Inputs_30d_1000scen_7tech_2z"))
    settings["Gas price uncertainty"] = false
    inputs = load_input_data(joinpath("inputs", "Inputs_30d_1000scen_7tech_2z"), settings)
    configure_parallel_workers!(settings)
    SPs = build_all_subproblems(inputs, settings)
    R, L = length(inputs["Resources"]), inputs["Number of lines"]
    ub = inputs["Capacity upper bounds"]
    rng = Xoshiro(3)
    caps = [0.2 .* ub .* rand(rng, R) for _ in 1:3]
    lines = [100 .* rand(rng, L) for _ in 1:3]
    batch = run_subproblems_batch(SPs, inputs, settings, caps, lines; minimal_payload = false)
    for i in 1:3
        single = run_all_subproblems(SPs, inputs, settings, caps[i], lines[i]; minimal_payload = false)
        for I in CartesianIndices(single)
            a, b = batch[i][I], single[I]
            @assert isapprox(a["SP objective"], b["SP objective"]; rtol = 1e-6) "point $i scenario $(Tuple(I)): objective $(a["SP objective"]) vs $(b["SP objective"])"
            @assert isapprox(a["capacity dual"], b["capacity dual"]; rtol = 1e-5, atol = 1e-3) "point $i scenario $(Tuple(I)): capacity duals differ"
            @assert haskey(a, "Power price") "full payload expected"
        end
    end
    @info("[$(test_index)] PASSED: run_subproblems_batch matches run_all_subproblems on 3 points × $(length(SPs)) scenarios.")
    return true
end

#=========

GridFill on the real model (src/model/boundary_sampling/gridfill_*.jl), without a Benders run:
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
                          br.inv[i], br.ub_ev[i], br.lb_ev[i], br.mid_ev[i], br.ub_cvar[i], br.lb_cvar[i], br.mid_cvar[i], [0.1, 0.2]) for i in 1:2]
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
    for col in ("EV_OpCost", "CVaR_OpCost", "Parents", "Weights", "GF_Batch", "GF_Mid_EV", "GF_LB_CVaR", "GF_Pred_v", "GF_Novelty")
        @assert col in names(df) "Summary is missing $col"
    end
    settings["Write all scenarios flag"] && @assert count(startswith("OperatingCost_D_"), names(df)) == S * F * K "per-scenario op cost columns missing"
    @assert isapprox(df.EV_OpCost, [solve_ev(fs) for fs in solved]; rtol = 1e-9) "Summary labels must be the solved values"
    @info("[$(test_index)] PASSED: pool cuts match, GridFill rows bracketed and written to $(folder).")
    return true
end
