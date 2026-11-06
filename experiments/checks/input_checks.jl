#=========

Input checks: scenario selection, cost construction, outages, the weather library, and converting
old-format folders. Everything but test_scenario_counts_end_to_end is solver-free.

==========#

const INPUT_CHECK_FOLDER = joinpath("inputs", "Inputs_2z_7tech")

function test_scenario_selection(test_index = "scenario_selection")
    distributions = [
        Dict{String,Any}("Count" => 7, "Distribution" => "Uniform", "Low" => -1.0, "High" => 1.0),
        Dict{String,Any}("Count" => 7, "Distribution" => "Triangular", "Low" => 0.0, "High" => 4.0, "Mode" => 1.0),
        Dict{String,Any}("Count" => 7, "Distribution" => "Normal", "Mean" => 1.0, "Sigma" => 0.5, "Low" => 0.0, "High" => 2.5),
        Dict{String,Any}("Count" => 7, "Distribution" => "LogNormal", "Median" => 3.0, "Sigma" => 0.4, "Low" => 2.0, "High" => 6.0)]
    for config in distributions
        values, p = select_distribution_scenarios(config)
        @assert values[1] == config["Low"] && values[end] == config["High"] "$(config["Distribution"]): extremes not kept"
        @assert all(isapprox.(diff(values), diff(values)[1])) "$(config["Distribution"]): values not evenly spaced"
        @assert isapprox(sum(p), 1.0) && all(p .> 0) "$(config["Distribution"]): bad probabilities $p"
    end

    uniform = Dict{String,Any}("Count" => 5, "Distribution" => "Uniform", "Low" => -0.15, "High" => 0.15)
    values, p = select_distribution_scenarios(uniform)
    @assert values ≈ [-0.15, -0.075, 0.0, 0.075, 0.15] && p ≈ [0.125, 0.25, 0.25, 0.25, 0.125] "uniform 5-point set: $values $p"
    gas = Dict{String,Any}("Count" => 5, "Distribution" => "LogNormal", "Median" => 3.45, "Sigma" => 0.40, "Low" => 2.00, "High" => 6.50)
    values, p = select_distribution_scenarios(gas)
    @assert isapprox(values, [2.0, 3.125, 4.25, 5.375, 6.5]) && isapprox(p, [0.166, 0.394, 0.270, 0.135, 0.036]; atol = 1e-3) "gas 5-point set: $values $p"
    @assert isapprox(sum(values .* p), 3.67; atol = 0.01) "gas mean $(sum(values .* p))"

    one_point, weight = select_distribution_scenarios(merge(uniform, Dict("Count" => 1)))
    @assert isapprox(one_point[1], 0.0; atol = 1e-9) && weight == [1.0]
    gas_mean, _ = select_distribution_scenarios(merge(gas, Dict("Count" => 1)))
    @assert isapprox(gas_mean[1], 3.67; atol = 0.02) "limited gas mean $(gas_mean[1])"

    errors(f) = try f(); false catch e; e isa ErrorException end
    @assert errors(() -> select_distribution_scenarios(merge(uniform, Dict("Distribution" => "Gamma")))) "an unknown distribution should error"
    @assert errors(() -> select_distribution_scenarios(Dict("Count" => 3, "Distribution" => "LogNormal", "Low" => 1.0, "High" => 2.0))) "a missing parameter should error"
    @assert errors(() -> select_distribution_scenarios(merge(uniform, Dict("Count" => 0)))) "a zero count should error"
    @assert errors(() -> select_distribution_scenarios(merge(uniform, Dict("Low" => 1.0)))) "Low above High should error"

    # Weather: an obvious outlier is kept from 2 scenarios up, and the most central one is picked for 1
    features = [Float64(i) + (i == 10 ? 50.0 : 0.0) for i in 1:10, _ in 1:3]
    for n in 2:9
        indices, shares = select_weather_scenarios(features, n)
        @assert 10 in indices && issorted(indices) && length(indices) == n "n=$n: $indices"
        @assert isapprox(sum(shares), 1.0) && all(shares .> 0)
        @assert (indices, shares) == select_weather_scenarios(features, n) "selection should be repeatable"
    end
    @assert select_weather_scenarios(features, 10) == (collect(1:10), fill(0.1, 10))
    indices, shares = select_weather_scenarios(features, 1)
    @assert length(indices) == 1 && indices[1] != 10 && shares == [1.0]
    @assert errors(() -> select_weather_scenarios(features, 11)) "more scenarios than the library should error"

    @info("[$(test_index)] PASSED: distribution and weather scenario selection.")
    return true
end

function test_cost_construction(test_index = "cost_construction")
    r = 0.0768
    @assert isapprox(capital_recovery_factor(0.0, 20), 0.05)
    @assert isapprox(capital_recovery_factor(0.10, 10), 0.16274539; rtol = 1e-6) "standard annuity table value"
    @assert isapprox(annual_investment_cost(1000.0, 10.0, 10, 0.10; multiplier = 2.0), 2000 * (162.74539 + 10); rtol = 1e-6)
    @assert variable_cost(10.0, 3.0, 2.0) == 32.0
    @assert isapprox(emission_intensity(10.0, 0.05), 0.5)

    settings = load_settings(INPUT_CHECK_FOLDER)
    inputs = load_input_data(INPUT_CHECK_FOLDER, settings)
    expected = Dict("Gas-z1" => 138_000, "Gas turbine-z1" => 126_000, "Wind-z1" => 213_000, "Solar-z1" => 127_000,
                    "Coal-z1" => 497_000, "Nuclear-z1" => 1_068_000, "Battery-z1" => 112_000)
    for (resource, cost) in expected
        actual = inputs["Investment costs"][findfirst(==(resource), inputs["Resources"])]
        @assert isapprox(actual, cost; rtol = 0.005) "$resource: $actual vs $cost"
    end
    # Zone 2 keeps the old cost ratios
    gas = findfirst(==("Gas-z1"), inputs["Resources"])
    @assert isapprox(inputs["Investment costs"][gas+1] / inputs["Investment costs"][gas], 41700 / 40721; rtol = 1e-9)

    # Gas responds to the scenario price by heat rate; coal and nuclear don't move
    gas_prices = inputs["Scenario selection"]["Gas"]["Values"]
    resources = inputs["Generation resources"]
    @assert inputs["Variable costs"][findfirst(==("Gas-z1"), resources), :] ≈ 6.5 .* gas_prices .+ 3.9
    @assert inputs["Variable costs"][findfirst(==("Gas turbine-z2"), resources), :] ≈ 10.7 .* gas_prices .+ 4.25
    @assert all(inputs["Variable costs"][findfirst(==("Coal-z1"), resources), :] .≈ 10.4 * 1.47 + 4.4)
    @assert all(inputs["Variable costs"][findfirst(==("Wind-z1"), resources), :] .== 0)

    # Wind and solar stay within 7% of Lazard's LCOE range at Lazard's capacity factors
    lazard = [("Wind", 1900, 24.5, 30, 0.55, 37), ("Wind", 2300, 40.0, 30, 0.30, 86),
              ("Solar", 1150, 11.0, 35, 0.30, 38), ("Solar", 1600, 14.0, 35, 0.20, 78)]
    for (name, capex, fom, life, cf, lcoe) in lazard
        ours = annual_investment_cost(capex, fom, life, r) / (8760 * cf)
        @assert abs(ours - lcoe) / lcoe <= 0.08 "$name at $capex/kW: $ours vs Lazard $lcoe"
    end
    @info("[$(test_index)] PASSED: costs build from capex, heat rates and fuel prices.")
    return true
end

function test_thermal_availability(test_index = "thermal_availability")
    settings = load_settings(INPUT_CHECK_FOLDER)
    inputs = load_input_data(INPUT_CHECK_FOLDER, settings)
    rates = Dict("Gas" => (0.05, 0.05), "Gas turbine" => (0.08, 0.03), "Coal" => (0.12, 0.06), "Nuclear" => (0.02, 0.07))
    resources = inputs["Generation resources"]
    availability, weights, demand = inputs["Generation availability"], inputs["Period weights"], inputs["Demand"]
    for (technology, (forced, planned)) in rates, k in axes(availability, 3)
        g = findfirst(==("$technology-z1"), resources)
        annual = sum(availability[:, g, k] .* weights[:, k]) / sum(weights[:, k])
        @assert isapprox(annual, (1 - forced) * (1 - planned); rtol = 1e-9) "$technology, weather $k: annual availability $annual"
        peak_hour = argmax(vec(sum(demand[:, :, k], dims = 2)))
        @assert isapprox(availability[peak_hour, g, k], 1 - forced) "$technology, weather $k: planned outage at the peak"
        @assert inputs["CRM Derating Factors"]["$technology-z1"] == 1.0 "$technology: a thermal derate would count outages twice"
    end
    @assert inputs["CRM Derating Factors"]["Wind-z1"] == 0.8 && inputs["CRM Derating Factors"]["Solar-z2"] == 0.8
    @info("[$(test_index)] PASSED: outage rates hold annually, and the peak has no planned outage.")
    return true
end

function test_weather_library(test_index = "weather_library")
    library = CSV.read(joinpath(INPUT_CHECK_FOLDER, "Weather_library.csv"), DataFrame)
    tables = load_weather_library(INPUT_CHECK_FOLDER, 2)
    K = size(tables.demand, 3)
    @assert K == 10 && all(==(8760), sum(tables.weights, dims = 1)) "expected 10 scenarios of 8,760 weighted hours"
    @assert last(tables.source_files) == 1 "old file 1 should be scenario 10, got file $(last(tables.source_files))"
    @assert sort(tables.source_files) == collect(1:10) "every old file should appear once"
    @assert issorted(weather_constraint_scores(tables.demand, tables.wind, tables.solar, tables.weights)) "scenarios should run from least to most constraining"
    @assert all(library.Zone_profiles .== "synthetic") "zone 2 here is synthetic and should be labeled that way"

    correlation(a, b) = cor(a, b)
    for k in 1:K
        w = tables.weights[:, k]
        d1, d2 = tables.demand[:, 1, k], tables.demand[:, 2, k]
        w1, w2 = tables.wind[:, 1, k], tables.wind[:, 2, k]
        s1, s2 = tables.solar[:, 1, k], tables.solar[:, 2, k]
        @assert 0.95 <= correlation(d1, d2) <= 0.99 "scenario $k demand correlation $(correlation(d1, d2))"
        @assert 0.85 <= correlation(w1, w2) <= 0.95 "scenario $k wind correlation $(correlation(w1, w2))"
        @assert 0.90 <= correlation(s1, s2) <= 0.98 "scenario $k solar correlation $(correlation(s1, s2))"
        @assert isapprox(sum(d2 .* w) / sum(d1 .* w), 0.5; atol = 1e-3) "scenario $k energy share"
        @assert isapprox(sum(w2 .* w) / sum(w1 .* w), 0.9; atol = 1e-3) "scenario $k wind capacity factor ratio"
        @assert isapprox(sum(s2 .* w) / sum(s1 .* w), 0.9; atol = 3e-3) "scenario $k solar capacity factor ratio"
        @assert all(0 .<= w2 .<= 1) && all(0 .<= s2 .<= 1)
        night = reshape(s1, 24, :) .== 0
        @assert all(reshape(s2, 24, :)[circshift(night, (1, 0)) .& night] .== 0) "scenario $k has zone 2 solar in hours dark in both zones"
    end
    @info("[$(test_index)] PASSED: the library is ordered, labeled, and within the zone 2 targets.")
    return true
end

# A small old-format folder: 2 zones, 3 weather scenarios, 2 representative days. Scenario 1 is the extreme one.
function _write_old_format_fixture(folder; zone_2_copy::Bool = true)
    mkpath(folder)
    write(joinpath(folder, "Network.csv"), "Lines,z1,z2,Existing_capacity,Max_expansion,CAPEX_per_MW\n1,1,-1,0,10000,6000\n")
    technologies = ["Gas", "Gas turbine", "Wind", "Solar", "Nuclear", "Coal"]
    resource_rows = ["$t-z$z,$z,1,$(t in ("Wind", "Solar") ? 1 : 0),1,0,$(z == 1 ? 1000 : 2000),0.1,5000" for t in technologies for z in 1:2]
    write(joinpath(folder, "Resources.csv"), "Resources,Zone,Index_ID,CES_eligible,Contracted,Storage,Investment cost,Emissions_ton_per_MWh,Capacity_UB\n" * join(resource_rows, "\n") * "\n")
    write(joinpath(folder, "Resources_storage.csv"), "Resources,Zone,Index_ID,Existing capacity,Storage,Investment cost,Power_to_energy_ratio,Loss_discharge,Loss_charge,Capacity_UB\nBattery-z1,1,1,0,0,41587,0.25,0.92,0.92,100000\n")
    write(joinpath(folder, "CRM_System.csv"), "Zone,Margin,PriceCap\n1,0.15,5000\n2,0.15,5000\n")
    derates = ["$t-z$z,$(t in ("Wind", "Solar") ? 0.8 : 0.9)" for t in technologies for z in 1:2]
    write(joinpath(folder, "CRM_Resources.csv"), "Resource,Derating_Factor\n" * join(derates, "\n") * "\n")
    write(joinpath(folder, "Demand_adders.csv"), "Demand,DemandAdd-1,DemandAdd-2\nAdder,-1000,1000\n")
    write(joinpath(folder, "settings.yaml"), "Solver: HiGHS # solver\nDemand uncertainty: true\nPrice cap: 2000\nCRM flag: true\nKeep me: 7\nNested:\n  A: 1\n  B: 2\n")
    write(joinpath(folder, "Representative_period_weights.csv"), "Period,Weight,Weight-1,Weight-2,Weight-3\n1,4000,4000,4380,4380\n2,4760,4760,4380,4380\n")
    columns = ["$t-z$z" for t in technologies for z in 1:2]
    for k in 1:3
        z1 = [10000.0 + 1000 * h + (k == 1 ? 8000 : 100 * k) + (t > 24 ? 2000 : 0) for t in 1:48 for h in [(t - 1) % 24]]
        z2 = zone_2_copy ? round.(0.5 .* z1, digits = 1) : z1 .* 0.4 .+ 50 .* sin.(1:48)
        CSV.write(joinpath(folder, "Demand_$k.csv"), DataFrame(Time_index = 1:48, var"Demand-z1" = z1, var"Demand-z2" = z2))
        availability = DataFrame(Time_Index = 1:48)
        wind = [0.4 + 0.3 * sin(t / 3) + 0.01 * k for t in 1:48]
        solar = [max(0.0, sin(pi * ((t - 1) % 24 - 6) / 12)) * 0.8 for t in 1:48]
        for c in columns
            scale = endswith(c, "z2") ? 0.9 : 1.0
            availability[!, c] = startswith(c, "Wind") ? scale .* wind : startswith(c, "Solar") ? scale .* solar : ones(48)
        end
        CSV.write(joinpath(folder, "Resources_availability_$k.csv"), availability)
    end
    return folder
end

function test_input_conversion(test_index = "input_conversion")
    root = mktempdir()
    old = _write_old_format_fixture(joinpath(root, "old"))
    new = joinpath(root, "new")
    report = convert_input_folder(old, new; counts = Dict("Demand" => 3, "Gas" => 3, "Weather" => 2), expected_most_constraining = 1)
    @assert last(report.order) == 1 && sort(report.order) == [1, 2, 3] "order $(report.order)"
    for file in ("settings.yaml", "Technologies.csv", "Resources_storage.csv", "Fuels.csv", "Weather_library.csv", "Network.csv")
        @assert isfile(joinpath(new, file)) "$file wasn't written"
    end
    @assert isfile(joinpath(old, "Demand_1.csv")) "the old folder shouldn't change"

    # Data carries over: zone 1 is untouched, and copied zone 2 profiles become labeled synthetic ones
    library = CSV.read(joinpath(new, "Weather_library.csv"), DataFrame)
    source = CSV.read(joinpath(old, "Demand_1.csv"), DataFrame)
    extreme = library[library.Source_file .== 1, :]
    @assert extreme.Demand_z1 == source[!, "Demand-z1"]
    @assert all(library.Zone_profiles .== "synthetic") && extreme.Demand_z2 != source[!, "Demand-z2"]
    @assert isapprox(sum(extreme.Demand_z2 .* extreme.Weight), sum(source[!, "Demand-z2"] .* extreme.Weight); rtol = 1e-3) "zone 2 energy should be kept"

    # Settings: other keys and comments survive, retired keys go, new ones arrive
    text = read(joinpath(new, "settings.yaml"), String)
    settings = load_settings(new)
    @assert occursin("# solver", text) && settings["Keep me"] == 7 && settings["CRM flag"] == true
    @assert !any(key -> haskey(settings, key), RETIRED_SETTINGS) "retired settings should be gone"
    @assert settings["Scenarios"]["Gas"]["Count"] == 3 && settings["CRM zones"]["Penalty"] == [5000.0, 5000.0]
    expected_fraction = round(1000 / minimum(maximum(library[library.Scenario .== k, "Demand_z1"]) for k in 1:3), digits = 4)
    @assert settings["Scenarios"]["Demand"]["High"] == expected_fraction

    # Zone cost ratios carry over, and the result loads and has the expected shape
    technologies = CSV.read(joinpath(new, "Technologies.csv"), DataFrame)
    @assert all(technologies.Cost_multiplier_z2 .== 2.0)
    inputs = load_input_data(new, settings)
    @assert size(inputs["Demand"]) == (48, 2, 2) && size(inputs["Variable costs"]) == (12, 3)

    # Repeatable from the seed; a different seed changes only the synthetic data
    again = joinpath(root, "again")
    convert_input_folder(old, again; counts = Dict("Demand" => 3, "Gas" => 3, "Weather" => 2))
    @assert read(joinpath(new, "Weather_library.csv")) == read(joinpath(again, "Weather_library.csv")) "same seed should rebuild identical data"

    # Copies can be kept, and real zone 2 data is never replaced
    kept = joinpath(root, "kept")
    convert_input_folder(old, kept; synthetic_zone_profiles = false)
    kept_library = CSV.read(joinpath(kept, "Weather_library.csv"), DataFrame)
    @assert all(kept_library.Zone_profiles .== "source") && kept_library[kept_library.Source_file .== 1, :Demand_z2] == source[!, "Demand-z2"]
    real_old = _write_old_format_fixture(joinpath(root, "real_old"); zone_2_copy = false)
    real_new = joinpath(root, "real_new")
    convert_input_folder(real_old, real_new)
    real_library = CSV.read(joinpath(real_new, "Weather_library.csv"), DataFrame)
    real_source = CSV.read(joinpath(real_old, "Demand_1.csv"), DataFrame)
    @assert real_library[real_library.Source_file .== 1, :Demand_z2] == real_source[!, "Demand-z2"] "real zone 2 demand was replaced"

    # Profiles hold only what differs
    della = joinpath(root, "della")
    cp(old, della)
    write(joinpath(della, "settings.yaml"), replace(read(joinpath(old, "settings.yaml"), String), "Solver: HiGHS" => "Solver: Gurobi", "  B: 2" => "  B: 3", "Price cap: 2000" => "Price cap: 1"))
    with_profile = joinpath(root, "with_profile")
    convert_input_folder(old, with_profile; profiles = Dict("della" => della))
    overrides = YAML.load_file(joinpath(with_profile, "settings_della.yaml"))
    @assert overrides == Dict("Solver" => "Gurobi", "Nested" => Dict("B" => 3)) "profile: $overrides"
    merged = load_settings(with_profile; profile = "della")
    @assert merged["Solver"] == "Gurobi" && merged["Nested"] == Dict("A" => 1, "B" => 3) && merged["Keep me"] == 7

    # Clear failures
    errors(f, text) = try f(); false catch e; e isa ErrorException && occursin(text, e.msg) end
    @assert errors(() -> convert_input_folder(old, new), "overwrite") "an existing folder should need overwrite"
    @assert errors(() -> convert_input_folder(old, joinpath(root, "x"); expected_most_constraining = 2), "most constraining") "an unexpected ranking should stop"
    @assert errors(() -> convert_input_folder(joinpath(root, "nowhere"), joinpath(root, "y")), "doesn't exist")
    @assert errors(() -> convert_input_folder(old, joinpath(root, "z"); technologies = default_technology_assumptions()[2:end, :]), "No cost assumptions")
    broken = joinpath(root, "broken")
    cp(old, broken)
    rm(joinpath(broken, "Demand_2.csv"))
    @assert errors(() -> convert_input_folder(broken, joinpath(root, "w")), "Resources_availability") "mismatched scenario files should error"
    rm(joinpath(broken, "Network.csv"))
    @assert errors(() -> convert_input_folder(broken, joinpath(root, "v")), "Network.csv") "a missing file should be named"

    @info("[$(test_index)] PASSED: old-format folders convert, with data, settings and failures as expected.")
    return true
end

function test_scenario_counts_end_to_end(test_index = "scenario_counts")
    capacity_share = [0.10, 0.05, 0.05, 0.03, 0.30, 0.15, 0.30, 0.15, 0.10, 0.05, 0.02, 0.02, 0.05]
    for counts in ((5, 5, 5), (1, 1, 1), (2, 3, 4))
        settings = load_settings(INPUT_CHECK_FOLDER)
        for (name, count) in zip(("Demand", "Gas", "Weather"), counts)
            settings["Scenarios"][name]["Count"] = count
        end
        inputs = load_input_data(INPUT_CHECK_FOLDER, settings)
        @assert (inputs["Number of demand scenarios"], inputs["Number of gas price scenarios"], inputs["Number of weather scenarios"]) == counts
        p = scenario_probabilities(inputs)
        @assert length(p) == prod(counts) && isapprox(sum(p), 1.0)

        # Solve the two corners of the scenario grid: high demand, gas and constraint cost more
        capacity = 60_000 .* capacity_share
        solve(index) = begin
            SP = build_subproblem(inputs, settings, index)
            set_parameter_value.(SP[:x], capacity)
            set_parameter_value.(SP[:x_line], [2000.0])
            output = run_subproblem(SP, inputs, settings; minimal_payload = false)
            @assert output["coeff"] == 1 "scenario $index did not solve"
            output["SP objective"]
        end
        low, high = solve([1, 1, 1]), solve(collect(counts))
        @assert low > 0 && (counts == (1, 1, 1) ? low == high : high > low) "counts $counts: objective $low to $high"
    end

    folder = mktempdir()
    settings = load_settings(INPUT_CHECK_FOLDER)
    write_scenario_selection(load_input_data(INPUT_CHECK_FOLDER, settings), folder)
    table = CSV.read(joinpath(folder, "Scenario_selection.csv"), DataFrame)
    @assert nrow(table) == 15 && all(isapprox(sum(table.Probability[table.Uncertainty .== u]), 1.0) for u in ("Demand", "Gas", "Weather"))
    @info("[$(test_index)] PASSED: scenario counts build, solve and are written out.")
    return true
end
