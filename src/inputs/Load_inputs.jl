#==========

Builds the `inputs` dictionary from one input folder (see Convert_inputs.jl for the file layout):
costs come from Technologies.csv and Fuels.csv, gas and demand scenarios from the distributions in
settings.yaml, and weather scenarios from Weather_library.csv.

==========#

# Reads Weather_library.csv into arrays over the whole library: demand, wind and solar are
# T × Z × K, and weights are T × K in hours per model hour.
function load_weather_library(inputs_path, Z::Int)
    library = CSV.read(joinpath(inputs_path, "Weather_library.csv"), DataFrame)
    sort!(library, [:Scenario, :Period, :Hour])
    scenarios = sort(unique(library.Scenario))
    K = length(scenarios)
    K > 0 || error("Weather_library.csv in $inputs_path has no scenarios.")
    scenarios == collect(1:K) || error("Weather_library.csv scenarios must be numbered 1 to $K.")
    period_length = maximum(library.Hour)
    n_periods = maximum(library.Period)
    T = period_length * n_periods
    nrow(library) == T * K || error("Weather_library.csv should have $(T * K) rows ($K scenarios × $T hours), found $(nrow(library)).")

    column(name) = reshape(Float64.(library[!, name]), T, K)
    zonal(prefix) = cat([reshape(column("$(prefix)_z$z"), T, 1, K) for z in 1:Z]...; dims = 2)
    # Any <Name>_z<N> column other than demand is an availability profile a technology can name
    profile_names = unique(m[1] for m in match.(r"^(.+)_z\d+$", names(library)) if m !== nothing && m[1] != "Demand")
    return (demand = zonal("Demand"), wind = zonal("Wind"), solar = zonal("Solar"),
            profiles = Dict(name => zonal(name) for name in profile_names),
            weights = column("Weight") ./ period_length, period_length = period_length, n_periods = n_periods,
            source_files = [first(library.Source_file[library.Scenario .== k]) for k in 1:K])
end

# Planned outages go in low-demand periods: each period's factor is proportional to its gap below the
# year's peak and averages 1 across the year's weighted hours, so a technology's annual rate is kept.
function planned_outage_factors(system_demand::AbstractVector, weights::AbstractVector, period_length::Int)
    period_peak = vec(maximum(reshape(system_demand, period_length, :), dims = 1))
    period_weight = vec(sum(reshape(weights, period_length, :), dims = 1))
    gap = maximum(period_peak) .- period_peak
    weighted_gap = sum(gap .* period_weight) / sum(period_weight)
    weighted_gap > 0 || error("Every period has the same peak, so planned outages can't be placed.")
    return gap ./ weighted_gap
end

# Hourly availability (T × G × K) of each generator: wind and solar from the library, and thermal
# technologies from their forced and planned outage rates.
function generation_availability(technologies, resource_technology, resource_zone, library, selected)
    T, G, K = size(library.demand, 1), length(resource_technology), length(selected)
    availability = ones(T, G, K)
    for (index, k) in enumerate(selected), g in 1:G
        row = technologies[resource_technology[g], :]
        if !ismissing(row.Profile)
            haskey(library.profiles, row.Profile) || error("Unknown Profile \"$(row.Profile)\" for $(row.Technology). Valid options: $(join(sort!(collect(keys(library.profiles))), ", ")), or blank.")
            profiles = library.profiles[row.Profile]
            availability[:, g, index] = profiles[:, resource_zone[g], k]
        else
            factors = planned_outage_factors(vec(sum(library.demand[:, :, k], dims = 2)), library.weights[:, k], library.period_length)
            planned = row.Planned_outage_rate .* factors
            maximum(planned) <= 1 || error("$(row.Technology): planned outages exceed 100% in a period; lower Planned_outage_rate.")
            availability[:, g, index] = (1 - row.Forced_outage_rate) .* (1 .- repeat(planned, inner = library.period_length))
        end
    end
    return availability
end

# Capacity bounds (MW) of technology row `t` in zone `z`. Capacity_UB_z<N> columns override the
# technology-wide Capacity_UB, a blank meaning the technology isn't offered there, and
# Capacity_LB_z<N> columns give capacity that must stay (0 if absent).
function technology_capacity_bounds(technologies, t::Int, z::Int)
    upper = "Capacity_UB_z$z" in names(technologies) ? coalesce(technologies[t, "Capacity_UB_z$z"], 0.0) : technologies.Capacity_UB[t]
    lower = "Capacity_LB_z$z" in names(technologies) ? coalesce(technologies[t, "Capacity_LB_z$z"], 0.0) : 0.0
    lower <= upper || error("$(technologies.Technology[t]) in zone $z: Capacity_LB ($lower) is above Capacity_UB ($upper).")
    return Float64(lower), Float64(upper)
end

function load_input_data(inputs_path, settings)

    # Collect data here
    inputs = Dict{String, Any}()

    lines_input = CSV.read(joinpath(inputs_path, "Network.csv"), DataFrame, header=true)
    inputs["Number of lines"] = size(lines_input)[1]
    # Collect a list of column names matching z1 z2 etc
    zone_columns = filter(x -> occursin(r"^z\d+$", x), names(lines_input))
    Z = length(zone_columns)
    inputs["Number of zones"] = Z
    inputs["Line capacities"] = Array(lines_input[:, "Existing_capacity"])
    inputs["Max expansion"] = Array(lines_input[:, "Max_expansion"])
    inputs["CAPEX per MW"] = Array(lines_input[:, "CAPEX_per_MW"])
    inputs["Zone map"] = Array(lines_input[:, zone_columns])

    # ~~~
    # Technologies: one row per technology, expanded here to one resource per zone
    # ~~~
    technologies = CSV.read(joinpath(inputs_path, "Technologies.csv"), DataFrame, header=true, stringtype=String)
    fuels = CSV.read(joinpath(inputs_path, "Fuels.csv"), DataFrame, header=true, stringtype=String)
    storage_input = CSV.read(joinpath(inputs_path, "Resources_storage.csv"), DataFrame, header=true, stringtype=String)
    discount_rate = settings["Discount rate"]

    # One resource per zone the technology is offered in
    offered = [(t, z) for t in 1:nrow(technologies) for z in 1:Z if technology_capacity_bounds(technologies, t, z)[2] > 0]
    resource_technology = first.(offered)
    resource_zone = last.(offered)
    gen_resources = ["$(technologies.Technology[t])-z$z" for (t, z) in zip(resource_technology, resource_zone)]
    G = length(gen_resources)
    O = nrow(storage_input)
    storage_flag = O > 0

    inputs["Generation resources"] = gen_resources
    inputs["Number of generation resources"] = G
    inputs["Number of storage resources"] = O
    inputs["Number of resources"] = G + O

    gen_zone_map = zeros(Int, G, Z)
    for g in 1:G
        gen_zone_map[g, resource_zone[g]] = 1
    end
    stor_zone_map = zeros(Int, O, Z)
    for o in 1:O
        stor_zone_map[o, storage_input[o, "Zone"]] = 1
    end
    # Full resource zone map: R × Z  (generation rows first, then storage rows)
    inputs["Resource zone map"] = vcat(gen_zone_map, stor_zone_map)

    # ~~~
    # Scenarios
    # ~~~
    scenario_settings = settings["Scenarios"]
    library = load_weather_library(inputs_path, Z)

    demand_fractions, demand_probabilities = select_distribution_scenarios(scenario_settings["Demand"])
    gas_prices, gas_probabilities = select_distribution_scenarios(scenario_settings["Gas"])
    features = weather_features(library.demand, library.wind, library.solar, library.weights; period_length = library.period_length)
    selected, weather_probabilities = select_weather_scenarios(features, scenario_settings["Weather"]["Count"])

    inputs["Scenario selection"] = Dict(
        "Demand" => Dict("Values" => demand_fractions, "Probabilities" => demand_probabilities),
        "Gas" => Dict("Values" => gas_prices, "Probabilities" => gas_probabilities),
        "Weather" => Dict("Values" => float.(selected), "Probabilities" => weather_probabilities,
                          "Source files" => library.source_files[selected]))
    inputs["Demand scenario probabilities"] = demand_probabilities
    inputs["Gas price scenario probabilities"] = gas_probabilities
    inputs["Weather scenario probabilities"] = weather_probabilities
    inputs["Number of demand scenarios"] = length(demand_probabilities)
    inputs["Number of gas price scenarios"] = length(gas_probabilities)
    inputs["Number of weather scenarios"] = length(weather_probabilities)

    # Demand: T × Z × K from the library, plus a flat MW adder sized to each zone. A zone's reference
    # peak is its lowest annual peak across the whole library, so it doesn't change with the selection.
    T = size(library.demand, 1)
    inputs["Time index"] = collect(1:T)
    inputs["Number of periods"] = T
    inputs["Demand"] = library.demand[:, :, selected]
    reference_peak = [minimum(maximum(library.demand[:, z, k]) for k in axes(library.demand, 3)) for z in 1:Z]
    inputs["Demand reference peak"] = reference_peak
    inputs["Demand adders"] = demand_fractions * reference_peak' # S × Z

    # Representative periods
    inputs["Number of representative periods"] = library.n_periods
    inputs["Set of first periods"] = collect(1:library.period_length:T)
    inputs["Set of last periods"] = collect(library.period_length:library.period_length:T)
    inputs["Period weights"] = library.weights[:, selected]
    inputs["Generation availability"] = generation_availability(technologies, resource_technology, resource_zone, library, selected)

    # ~~~
    # Costs
    # ~~~
    fuel_row(name) = begin
        match_row = findfirst(==(name), fuels.Fuel)
        match_row === nothing && error("Fuel \"$name\" is not in Fuels.csv.")
        fuels[match_row, :]
    end
    fuel_price(row, gas_price) = ismissing(row.Fuel) ? 0.0 : row.Fuel == "Gas" ? gas_price : fuel_row(row.Fuel).Price_per_MMBtu
    # Capture_rate is optional: the share of a technology's fuel carbon that isn't emitted
    captured(row) = "Capture_rate" in names(technologies) ? coalesce(row.Capture_rate, 0.0) : 0.0
    carbon(row) = ismissing(row.Fuel) ? 0.0 : fuel_row(row.Fuel).Carbon_t_per_MMBtu * (1 - captured(row))
    heat_rate(row) = ismissing(row.Heat_rate_MMBtu_per_MWh) ? 0.0 : row.Heat_rate_MMBtu_per_MWh

    # $/MWh, G × F
    inputs["Variable costs"] = [variable_cost(heat_rate(technologies[resource_technology[g], :]),
                                              fuel_price(technologies[resource_technology[g], :], price),
                                              technologies.Variable_OM_per_MWh[resource_technology[g]])
                                for g in 1:G, price in gas_prices]
    # Energy_limit_CF is optional: the share of full-capacity energy a technology can produce in each
    # representative period (see generation_energy_limit in the subproblem). Blank or 1 means no limit.
    inputs["Energy limit fractions"] = "Energy_limit_CF" in names(technologies) ?
        Float64.(coalesce.(technologies.Energy_limit_CF[resource_technology], 1.0)) : ones(G)
    all(0 .< inputs["Energy limit fractions"] .<= 1) || error("Energy_limit_CF must be above 0 and at most 1.")
    inputs["CO2 emission intensities"] = [emission_intensity(heat_rate(technologies[t, :]), carbon(technologies[t, :])) for t in resource_technology]

    # $/MW-yr, generation rows first
    gen_investment = [annual_investment_cost(technologies.Capex_per_kW[t], technologies.Fixed_OM_per_kW_yr[t], technologies.Life_years[t], discount_rate;
                                             multiplier = technologies[t, "Cost_multiplier_z$z"])
                      for (t, z) in zip(resource_technology, resource_zone)]
    stor_investment = [annual_investment_cost(storage_input.Capex_per_kW[o], storage_input.Fixed_OM_per_kW_yr[o], storage_input.Life_years[o], discount_rate) for o in 1:O]
    inputs["Investment costs"] = vcat(gen_investment, stor_investment)

    # ~~~
    # CRM: outages are already in availability, so thermal derates are 1.0
    # ~~~
    if settings["CRM flag"]
        inputs["CRM Margin"] = Float64.(settings["CRM zones"]["Margin"])
        inputs["CRM Price Cap"] = Float64.(settings["CRM zones"]["Penalty"])
        inputs["CRM Derating Factors"] = Dict(gen_resources[g] => technologies.CRM_derate[resource_technology[g]] for g in 1:G)
    end

    # ~~~
    # Storage and capacity limits
    # ~~~
    gen_bounds = [technology_capacity_bounds(technologies, t, z) for (t, z) in offered]
    x_min, x_max = first.(gen_bounds), last.(gen_bounds)
    if storage_flag
        inputs["Storage power to energy ratio"] = Float64.(storage_input[!,"Power_to_energy_ratio"])
        inputs["Storage discharge loss"] = Float64.(storage_input[!,"Loss_discharge"])
        inputs["Storage charge loss"] = Float64.(storage_input[!,"Loss_charge"])
        inputs["Resources"] = vcat(gen_resources, storage_input[1:end,1])
        inputs["Capacity upper bounds"] = vcat(x_max, storage_input[1:end,"Capacity_UB"])
        storage_min = "Capacity_LB" in names(storage_input) ? Float64.(coalesce.(storage_input[!,"Capacity_LB"], 0.0)) : zeros(O)
        inputs["Capacity lower bounds"] = vcat(x_min, storage_min)
    else
        inputs["Resources"] = gen_resources
        inputs["Capacity upper bounds"] = x_max
        inputs["Capacity lower bounds"] = x_min
    end
    inputs["Storage flag"] = storage_flag
    inputs["Storage"] = storage_input

    return inputs
end
