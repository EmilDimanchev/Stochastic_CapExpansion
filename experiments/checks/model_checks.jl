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
    @info("[$(test_index)] settings groups resolve and validate as expected for both input sets.")
    return true
end
