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
    inputs_folder = joinpath("inputs", "Inputs_2z_7tech")
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

Settings groups test (experiments/mga/run_settings.jl): every group block is present in the
base settings and the Della profile, overrides replace only the keys they name (a nested block such as
Cut settings' "Prune window" merges), unknown keys and groups error out, and resolving never
touches the caller's settings.

==========#

function test_settings_groups(test_index = "settings_groups")
    for profile in (nothing, "della")
        settings = load_settings(joinpath("inputs", "Inputs_2z_7tech"); profile = profile)
        @assert all(g -> haskey(settings, g), SETTINGS_GROUPS) "profile $profile is missing a settings group"
        CutPool(resolve_all_settings(settings)) # every Cut settings key is present
    end
    settings = load_settings(joinpath("inputs", "Inputs_2z_7tech"))
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
    # Mapping / Boundary / Delaunay keys error out
    sampling = resolve_settings(settings, "Sampling settings", Dict("Method" => "GridFill", "GridFill" => Dict("Cell size" => 0.1)))
    @assert sampling["GridFill"]["Cell size"] == 0.1 && sampling["GridFill"]["Per cell"] == settings["Sampling settings"]["GridFill"]["Per cell"] "GridFill block should merge"
    @assert sampling["Method"] in SAMPLING_METHODS
    @assert errors(() -> resolve_settings(settings, "Sampling settings", Dict("Boundary" => Dict()))) "the retired Boundary sampler block should error"
    @assert errors(() -> resolve_settings(settings, "Sampling settings", Dict("Delaunay" => Dict()))) "the retired Delaunay sampler block should error"
    @assert errors(() -> resolve_settings(settings, "Sampling settings", Dict("GridFill" => Dict("Cell sise" => 0.1)))) "a misspelled nested key should error"
    @assert errors(() -> resolve_settings(settings, "Mapping settings", Dict("Sampling" => true))) "the retired Mapping \"Sampling\" key should error"
    @assert errors(() -> resolve_settings(settings, "Mapping settings", Dict("Interior samples" => 5))) "the retired \"Interior samples\" key should error"
    @assert errors(() -> resolve_all_settings(settings, Dict("Boundary settings" => Dict()))) "the retired Boundary settings group should error"
    @info("[$(test_index)] settings groups resolve and validate as expected for both input sets.")
    return true
end

#=========

"Write case files" (Run settings): off skips each case's own folder, and leaves what feeds the
Summary unchanged. Solver-free: the solve results are a small synthetic dict, with eval SPs so
the *_eval files are covered too.

==========#

function test_case_file_flag(test_index = "case_file_flag")
    probs = [0.5, 0.5]
    sp(x) = Dict{String,Any}("Power price" => fill(x, 2, 1), "Cost by zone" => [x], "Emissions by zone" => [x], "Emissions" => x,
                             "SP objective" => x, "Generation" => fill(x, 1, 2), "Load shedding" => fill(0.0, 2, 1))
    SPs = reshape([sp(1.0), sp(2.0)], 2, 1, 1)
    inputs = Dict{String,Any}("Demand scenario probabilities" => probs, "Gas price scenario probabilities" => [1.0], "Weather scenario probabilities" => [1.0],
                              "Full Demand scenario probabilities" => probs, "Full Gas price scenario probabilities" => [1.0], "Full Weather scenario probabilities" => [1.0],
                              "Number of periods" => 2, "Number of lines" => 1, "Number of zones" => 1, "Resources" => ["Gas"], "Generation resources" => ["Gas"], "Time index" => 1:2)
    results = Dict{String,Any}("MP" => Dict{String,Any}("Capacity" => [10.0], "Line expansion" => [1.0], "Inv_cost" => 5.0, "Inv cost by zone" => [5.0]),
                               "SPs" => SPs, "CVaR" => 2.0, "Expected Value" => 1.5, "SPs_eval" => SPs, "CVaR_Eval" => 2.0, "Expected Value Eval" => 1.5,
                               "Gaps" => [0.1, 0.01], "Capacity per iteration" => [[9.0], [10.0]], "Time MP hist" => [0.1, 0.1], "Time SP hist" => [0.2, 0.2])
    settings(on) = Dict{String,Any}("Risk aversion weight" => 0.5, "Write all scenarios flag" => true, "Run settings" => Dict{String,Any}("Write case files" => on))

    root = mktempdir()
    folder_on, folder_off = joinpath(root, "on", "Case"), joinpath(root, "off", "Case")
    dfs_on = write_results_benders(results, inputs, settings(true), folder_on)
    dfs_off = write_results_benders(results, inputs, settings(false), folder_off)

    expected = ["Convergence.csv", "generation.csv", "generation_eval.csv", "nse.csv", "nse_eval.csv", "price.csv", "price_eval.csv"]
    @assert sort(readdir(folder_on)) == expected "on should write every case file, got $(readdir(folder_on))"
    @assert !isdir(folder_off) "off should not create the case folder"
    @assert all(a == b for (a, b) in zip(dfs_on, dfs_off)) "the flag should not change what feeds the Summary"
    @assert count(startswith("OperatingCost_D_"), names(dfs_off[2])) == 2 "per-scenario Summary columns follow \"Write all scenarios flag\", not this flag"

    # Gradients.csv lives in the case folders, so the flag also turns it off, without touching the caller's settings
    base = load_settings(joinpath("inputs", "Inputs_2z_7tech"))
    @assert base["Cut export settings"]["Export cuts"] == resolve_all_settings(base, Dict("Run settings" => Dict("Write case files" => false)))["Cut export settings"]["Export cuts"] "the flag should not change the cut export"
    gradients_on = Dict("Sampling settings" => Dict("Write gradients" => true))
    @assert resolve_all_settings(base, merge(gradients_on, Dict("Run settings" => Dict("Write case files" => true))))["Sampling settings"]["Write gradients"] "gradients stay on while case files are on"
    resolved = resolve_all_settings(base, merge(gradients_on, Dict("Run settings" => Dict("Write case files" => false))))
    @assert !resolved["Sampling settings"]["Write gradients"] "gradients should follow the case-file switch"
    @assert base["Sampling settings"]["Write gradients"] == false "resolving modified the caller's settings"

    @info("[$(test_index)] PASSED: case files follow the flag; Summary data and the cut export setting don't.")
    return true
end

#=========

Batched subproblem test (src/model/benders/worker_functions.jl): run_subproblems_batch solves
several capacity points with one call per worker over its own scenarios. Its objectives and
duals must match run_all_subproblems point by point (one gas and weather scenario: 5 scenarios).

==========#

function test_subproblems_batch(test_index = "sp_batch")
    settings = load_settings(joinpath("inputs", "Inputs_2z_7tech"))
    settings["Scenarios"]["Gas"]["Count"] = 1
    settings["Scenarios"]["Weather"]["Count"] = 1
    inputs = load_input_data(joinpath("inputs", "Inputs_2z_7tech"), settings)
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

Design CSV evaluation test:
  - read_design_csv maps the planner CSVs' <Tech>_Z<z>/Transmission columns (GW) onto
    capacity_variable_names (MW) and rejects a CSV with a missing or
    unknown capacity column.
  - evaluate_design_csv (batched, with a partial last batch) writes one Summary row per design
    with the Design_ metadata, and its EV/CVaR match a direct run_all_subproblems solve.
One gas and one weather scenario (5 demand scenarios). Investment_B in these CSVs reflects the
costs they were designed under, so it is not compared with current costs.

==========#

function test_design_csv_evaluation(test_index = "design_csv")
    inputs_folder = joinpath("inputs", "Inputs_2z_7tech")
    settings = load_settings(inputs_folder)
    settings["Scenarios"]["Gas"]["Count"] = 1
    settings["Scenarios"]["Weather"]["Count"] = 1
    inputs = load_input_data(inputs_folder, settings)
    R = length(inputs["Resources"])
    planner_file = joinpath("experiments", "designs", "run13_planner_designs.csv")
    ra_file = joinpath("experiments", "designs", "run13_RA_planner_designs.csv")

    # Parsing
    for (file, n) in ((planner_file, 500), (ra_file, 100))
        df = CSV.read(file, DataFrame)
        labels, X, meta = read_design_csv(df, inputs)
        @assert length(X) == n && all(length(x) == R + inputs["Number of lines"] for x in X) "$file: $(length(X)) designs"
        @assert labels == string.(df.point_id)
        @assert all("Design_" * c in names(meta) for c in ("Investment_B", "EV_Op_ceil_B", "CVaR_Op_ceil_B"))
    end
    labels, X, _ = read_design_csv(planner_file, inputs)
    names_x = capacity_variable_names(inputs)
    @assert labels[1] == "P0001"
    @assert X[1][findfirst(==("Gas-z1"), names_x)] ≈ 2840.5
    @assert X[1][findfirst(==("Line 1"), names_x)] ≈ 2041.6
    df = CSV.read(ra_file, DataFrame)
    for bad in (select(df, Not("Coal_Z2")), rename(df, "Wind_Z1" => "Wind_Z9"))
        threw = try read_design_csv(bad, inputs); false catch; true end
        @assert threw "a design CSV with a missing/unknown capacity column should error"
    end

    # Batched evaluation of 3 designs, batch size 2
    configure_parallel_workers!(settings)
    SPs = build_all_subproblems(inputs, settings)
    folder = joinpath("outputs", "Test_"*string(test_index), "Design_Evaluation", "Summary")
    done = evaluate_design_csv(SPs, inputs, settings, ra_file, folder; batch_size = 2, rows = 1:3)
    summary = CSV.read(joinpath(folder, "Summary_run13_RA_planner_designs.csv"), DataFrame)
    @assert done == ["RA0001", "RA0002", "RA0003"] && nrow(summary) == 3 && string.(summary.Index) == done
    @assert all(c in names(summary) for c in ("Design_Investment_B", "Design_ELCC_GasCC", "Design_Solve_OK"))
    _, X_ra, _ = read_design_csv(ra_file, inputs)
    for i in 1:3
        sp = run_all_subproblems(SPs, inputs, settings, X_ra[i][1:R], X_ra[i][R+1:end]; minimal_payload = false)
        ev, cvar = evaluate_subproblems(sp, inputs["Demand scenario probabilities"], inputs["Gas price scenario probabilities"],
                                        inputs["Weather scenario probabilities"], settings["Value-at-Risk percent"])
        @assert isapprox(summary[i, "EV_OpCost"], ev; rtol = 1e-6) "$(done[i]): EV $(summary[i, "EV_OpCost"]) vs direct $ev"
        @assert isapprox(summary[i, "CVaR_OpCost"], cvar; rtol = 1e-6) "$(done[i]): CVaR $(summary[i, "CVaR_OpCost"]) vs direct $cvar"
    end
    @info("[$(test_index)] PASSED: design CSVs parse (600 designs) and 3 batched evaluations match direct solves.")
    return true
end

#=========

Zonal cost test (src/model/benders/worker_functions.jl): "Cost by zone" must be weighted like the
objective, so the zones sum to "SP objective". Capacity is set low enough that generation,
non-served energy and CRM slack costs are all non-zero, and repeated with the CRM off.

==========#

function test_cost_by_zone(test_index = "cost_by_zone")
    inputs_folder = joinpath("inputs", "Inputs_2z_7tech")
    for crm_flag in (true, false)
        settings = load_settings(inputs_folder)
        settings["CRM flag"] = crm_flag
        inputs = load_input_data(inputs_folder, settings)
        SP = build_subproblem(inputs, settings, [5, 1, 1])
        set_parameter_value.(SP[:x], fill(1500.0, length(inputs["Resources"])))
        set_parameter_value.(SP[:x_line], fill(500.0, inputs["Number of lines"]))
        output = run_subproblem(SP, inputs, settings; minimal_payload = false)

        @assert output["coeff"] == 1 "subproblem did not solve to optimality"
        @assert sum(output["Load shedding"]) > 0 "test point should shed load"
        @assert sum(value.(SP[:g])) > 0 "test point should generate"
        crm_flag && @assert sum(output["CRM slack"]) > 0 "test point should use CRM slack"
        total = sum(output["Cost by zone"])
        @assert isapprox(total, output["SP objective"]; rtol = 1e-6) "CRM $crm_flag: zones sum to $total, objective is $(output["SP objective"])"
    end
    @info("[$(test_index)] PASSED: Cost by zone sums to the SP objective with the CRM on and off.")
    return true
end

#=========

Temporal decomposition test (src/model/benders/worker_functions.jl): "Subproblems per scenario"
splits each scenario evenly into blocks of whole representative periods. Periods don't interact, so
the merged per-scenario results must match one subproblem per scenario:
  - objectives match, cost by zone sums to the objective, and time series keep their full-T
    shapes. The zonal split and emissions aren't unique (equally cheap flows), so their
    differences are only reported
    (30 subproblems of 1 day, and 5 of 6 days)
  - the summed duals give valid cuts. The LP can be degenerate, so duals are checked as
    subgradients (theta(x') >= theta(x) + beta'(x' - x)) rather than for equality
  - a parallel build gives the same results as the serial one, through run_all_subproblems
    and run_subproblems_batch
  - a count that doesn't divide the periods, and a block that cuts through a period, are rejected

==========#

function test_temporal_blocks(test_index = "temporal_blocks")
    folder = joinpath("inputs", "Inputs_2z_7tech")
    settings = load_settings(folder)
    settings["Scenarios"]["Gas"]["Count"] = 1
    settings["Scenarios"]["Weather"]["Count"] = 2
    settings["Scenarios"]["Demand"]["Count"] = 3
    inputs = load_input_data(folder, settings)
    R, L = length(inputs["Resources"]), inputs["Number of lines"]
    ub = inputs["Capacity upper bounds"]
    rng = Xoshiro(7)
    caps = [0.2 .* ub .* rand(rng, R) for _ in 1:3]
    lines = [100 .* rand(rng, L) for _ in 1:3]

    with_blocks(n_blocks, parallel) = merge(settings, Dict("Subproblems per scenario" => n_blocks, "Parallel flag" => parallel))
    solve_points(SPs, s) = [run_all_subproblems(SPs, inputs, s, caps[i], lines[i]; minimal_payload = false) for i in 1:3]
    rel(a, b) = abs(a - b) / max(1.0, abs(b))

    # Every cut from point i, evaluated at point j, must stay below the solved objective at j
    function check_cuts(results, label)
        worst = 0.0
        for i in 1:3, j in 1:3, I in CartesianIndices(results[i])
            i == j && continue
            o = results[i][I]
            predicted = o["SP objective"] + o["capacity dual"]' * (caps[j] - caps[i]) + o["line dual"]' * (lines[j] - lines[i])
            worst = max(worst, (predicted - results[j][I]["SP objective"]) / max(1.0, abs(results[j][I]["SP objective"])))
        end
        @assert worst <= 1e-4 "$label: a cut overestimates the objective at another point by $(worst) (relative)"
        return worst
    end

    serial = with_blocks(1, false)
    reference = solve_points(build_all_subproblems(inputs, serial), serial)
    n_periods = inputs["Number of representative periods"]

    for n_blocks in (n_periods, 5)
        s = with_blocks(n_blocks, false)
        blocks = subproblem_blocks(inputs, s)
        @assert length(blocks) == n_blocks "expected $(n_blocks) blocks, got $(length(blocks))"
        @assert allequal(length.(blocks)) "blocks should all be the same length"
        @assert reduce(vcat, collect.(blocks)) == collect(1:inputs["Number of periods"]) "blocks don't cover every hour once, in order"
        SPs = build_all_subproblems(inputs, s)
        @assert size(SPs, 4) == length(blocks)
        results = solve_points(SPs, s)
        max_dual_gap, max_split_gap, max_emissions_gap = 0.0, 0.0, 0.0
        for i in 1:3, I in CartesianIndices(reference[i])
            a, b = results[i][I], reference[i][I]
            where_ = "$n_blocks subproblems per scenario, point $i, scenario $(Tuple(I))"
            @assert a["coeff"] == 1 "$where_: a block did not solve to optimality"
            @assert rel(a["SP objective"], b["SP objective"]) <= 1e-6 "$where_: objective $(a["SP objective"]) vs $(b["SP objective"])"
            # Zonal splits and emissions aren't unique (e.g. equally cheap flows), so only the total is exact
            @assert rel(sum(a["Cost by zone"]), a["SP objective"]) <= 1e-6 "$where_: cost by zone doesn't sum to the objective"
            max_split_gap = max(max_split_gap, maximum(abs.(a["Cost by zone"] - b["Cost by zone"])) / max(1.0, b["SP objective"]))
            max_emissions_gap = max(max_emissions_gap, rel(a["Emissions"], b["Emissions"]))
            for key in ("Generation", "Power price", "Load shedding", "Net Flow", "CRM dual", "CRM slack", "Storage charging", "Storage discharging")
                haskey(b, key) && @assert size(a[key]) == size(b[key]) "$where_: $key has shape $(size(a[key])), expected $(size(b[key]))"
            end
            max_dual_gap = max(max_dual_gap, maximum(abs.(a["capacity dual"] - b["capacity dual"])) / max(1.0, maximum(abs.(b["capacity dual"]))))
        end
        worst = check_cuts(results, "$n_blocks subproblems per scenario")
        @info("[$(test_index)] $n_blocks subproblems per scenario ($(length(blocks[1])) hours each): objectives match; worst cut overestimate $(worst). Differences from one block per scenario, which degeneracy allows: capacity duals $(max_dual_gap), zonal cost split $(max_split_gap), emissions $(max_emissions_gap) (relative)")
    end
    check_cuts(reference, "one block per scenario")

    # Parallel build: same blocks spread over the workers
    parallel = with_blocks(5, true)
    configure_parallel_workers!(parallel)
    serial_5 = solve_points(build_all_subproblems(inputs, with_blocks(5, false)), with_blocks(5, false))
    SPs_par = build_all_subproblems(inputs, parallel)
    on_workers = solve_points(SPs_par, parallel)
    batch = run_subproblems_batch(SPs_par, inputs, parallel, caps, lines; minimal_payload = false)
    for i in 1:3, I in CartesianIndices(serial_5[i]), (label, res) in (("run_all_subproblems", on_workers), ("run_subproblems_batch", batch))
        @assert rel(res[i][I]["SP objective"], serial_5[i][I]["SP objective"]) <= 1e-6 "parallel $label, point $i, scenario $(Tuple(I)): objective differs from serial"
        @assert size(res[i][I]["Generation"]) == size(serial_5[i][I]["Generation"]) "parallel $label: Generation shape differs"
    end
    check_cuts(on_workers, "parallel run_all_subproblems")
    check_cuts(batch, "parallel run_subproblems_batch")

    # Counts that don't split the periods evenly, or aren't positive integers, are rejected with the valid values
    for bad in (7, 0, n_periods + 1, 2.5, nothing, "12")
        message = try
            subproblem_blocks(inputs, with_blocks(bad, false))
            ""
        catch e
            sprint(showerror, e)
        end
        @assert occursin("Valid values: 1, 2, 3, 5, 6, 10, 15, 30", message) "Subproblems per scenario = $(repr(bad)) should be rejected with the valid values, got: $(message)"
    end
    @assert length(subproblem_blocks(inputs, delete!(copy(settings), "Subproblems per scenario"))) == 1 "a missing setting should mean one subproblem per scenario"

    # A block must start and end on period boundaries
    period_length = inputs["Set of last periods"][1]
    rejected = try
        build_subproblem(inputs, serial, [1, 1, 1]; hours = 2:(period_length + 1))
        false
    catch e
        occursin("period boundaries", sprint(showerror, e))
    end
    @assert rejected "a block cutting through a period should be rejected"

    @info("[$(test_index)] PASSED: temporal blocks reproduce one subproblem per scenario (30 and 5 subproblems per scenario, serial and parallel).")
    return true
end
