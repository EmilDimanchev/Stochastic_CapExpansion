#==========

Converts an old-format input folder (Demand_N.csv, Resources_availability_N.csv, Variable_cost.csv,
Representative_period_weights.csv, ...) into the layout load_input_data reads:

  settings.yaml           base settings, with Scenarios, Discount rate and CRM zones
  settings_<profile>.yaml only the keys that differ for a profile, e.g. settings_della.yaml
  Technologies.csv        one row per technology: costs, outages, CRM derate, zone cost multipliers
  Resources_storage.csv   storage technologies, with capex-based costs
  Fuels.csv               carbon intensity and fixed fuel prices
  Weather_library.csv     long format, one row per scenario and hour, ordered least to most constraining
  Network.csv             unchanged

Usage, from the repo root after including experiments/MGA_Tests.jl:

  convert_input_folder("inputs/Old_folder", "inputs/New_folder")
  convert_input_folder("inputs/Old_folder", "inputs/New_folder"; profiles = Dict("della" => "inputs/Old_folder_Della"))

What carries over: demand, wind and solar profiles, period weights, network, capacity limits, zone
cost ratios, eligibility flags, CRM margins and penalties, and every other setting.

What doesn't: old $/MW-yr and $/MWh costs can't be turned back into capex and heat rates, so costs
come from `technologies` and `fuels` (defaults below). The function reports old against new costs.
Old fuel scenario tables are replaced by the `gas` distribution, and old thermal availability
columns by outage rates.

Zones whose demand, wind or solar is an exact scaled copy of zone 1 get a synthetic profile by
default, labeled in the Zone_profiles column. Set synthetic_zone_profiles = false to keep the copies.

==========#

# Lazard LCOE+ v18 (June 2025) midpoints. Outage rates: forced rates for gas and coal follow recent
# NERC figures and nuclear an ISO-NE class table; the gas turbine rate and planned rates are estimates.
function default_technology_assumptions()
    return DataFrame(
        Technology              = ["Gas", "Gas turbine", "Wind", "Solar", "Nuclear", "Coal"],
        Capex_per_kW            = [1400.0, 1300.0, 2100.0, 1375.0, 11920.0, 5307.5],
        Fixed_OM_per_kW_yr      = [17.75, 13.5, 32.25, 12.5, 147.0, 67.6],
        Life_years              = [30, 30, 30, 35, 70, 40],
        Fuel                    = ["Gas", "Gas", missing, missing, "Uranium", "Coal"],
        Heat_rate_MMBtu_per_MWh = [6.5, 10.7, missing, missing, 10.45, 10.4],
        Variable_OM_per_MWh     = [3.9, 4.25, 0.0, 0.0, 4.8, 4.4],
        Profile                 = [missing, missing, "Wind", "Solar", missing, missing],
        Forced_outage_rate      = [0.05, 0.08, 0.0, 0.0, 0.02, 0.12],
        Planned_outage_rate     = [0.05, 0.03, 0.0, 0.0, 0.07, 0.06])
end

default_fuel_assumptions() = DataFrame(Fuel = ["Gas", "Coal", "Uranium"],
                                       Carbon_t_per_MMBtu = [0.053, 0.095, 0.0],
                                       Price_per_MMBtu = [missing, 1.47, 0.85]) # gas comes from Scenarios in settings

default_storage_assumptions() = (Capex_per_kW = 904.5, Fixed_OM_per_kW_yr = 22.04, Life_years = 20) # 4-hour battery, 20-year life assumed

default_gas_distribution() = Dict{String,Any}("Distribution" => "LogNormal", "Median" => 3.45, "Sigma" => 0.40, "Low" => 2.00, "High" => 6.50)

# Settings the new format replaces
const RETIRED_SETTINGS = ("Demand uncertainty", "Gas price uncertainty", "Weather uncertainty", "Price cap")

_read_old(folder, file) = CSV.read(joinpath(folder, file), DataFrame, header = true, stringtype = String)

# Numbers of the Demand_N.csv / Resources_availability_N.csv files in `folder`, after checking the
# folder is in the old format
function old_format_scenarios(folder::AbstractString)
    isdir(folder) || error("Input folder $folder doesn't exist.")
    files = readdir(folder)
    numbers(pattern) = sort([parse(Int, m[1]) for m in match.(pattern, files) if m !== nothing])
    scenarios = numbers(r"^Demand_(\d+)\.csv$")
    required = ["Network.csv", "Resources.csv", "Resources_storage.csv", "CRM_System.csv", "CRM_Resources.csv",
                "Demand_adders.csv", "Representative_period_weights.csv", "settings.yaml"]
    missing_files = setdiff(required, files)
    isempty(scenarios) && push!(missing_files, "Demand_N.csv")
    isempty(missing_files) || error("$folder isn't an old-format folder this can convert; missing: $(join(missing_files, ", ")).")
    scenarios == numbers(r"^Resources_availability_(\d+)\.csv$") || error("$folder needs one Resources_availability_N.csv per Demand_N.csv.")
    return scenarios
end

# True when `profile` is a constant multiple of `reference`, within rounding
function is_scaled_copy(profile::AbstractVector, reference::AbstractVector; tolerance = 1e-3)
    ratio = sum(profile) / sum(reference)
    return maximum(abs.(profile .- ratio .* reference)) <= tolerance * sum(profile) / length(profile)
end

function _ar1_noise(rng, n::Int; persistence = 0.8)
    noise = zeros(n)
    noise[1] = randn(rng)
    for t in 2:n
        noise[t] = persistence * noise[t-1] + sqrt(1 - persistence^2) * randn(rng)
    end
    return noise
end

# Scales `profile` so its weighted total matches `target`, keeping values at or below `limit`
function _rescale_to_total(profile::AbstractVector, weights::AbstractVector, target::Real; limit = Inf)
    scaled = copy(profile)
    for _ in 1:50
        total = sum(scaled .* weights)
        total > 0 || return scaled
        isapprox(total, target; rtol = 1e-9) && break
        scaled = min.(scaled .* (target / total), limit)
    end
    return scaled
end

# Applies `f(day, rng)` to each representative day. Days aren't consecutive, so shifts wrap within a day.
function _by_day(f, profile::AbstractVector, period_length::Int, rng)
    return reduce(vcat, [f(day, rng) for day in eachcol(reshape(profile, period_length, :))])
end

# Synthetic profiles for a zone that's a scaled copy of zone 1. Each keeps the zone's weighted total,
# and noise levels are set to hit the hourly correlation targets checked in test_weather_library.
function synthetic_demand(copied, weights, period_length, rng; lag = 1, noise = 0.010)
    shaped = _by_day(copied, period_length, rng) do day, r
        circshift(day, lag) .* (1 .+ noise .* _ar1_noise(r, length(day)))
    end
    return round.(_rescale_to_total(shaped, weights, sum(copied .* weights)), digits = 1)
end

function synthetic_wind(copied, weights, period_length, rng; noise = 0.075)
    shaped = _by_day(copied, period_length, rng) do day, r
        clamp.(circshift(day, rand(r, 1:2)) .+ noise .* _ar1_noise(r, length(day)), 0.0, 1.0)
    end
    return _rescale_to_total(shaped, weights, sum(copied .* weights); limit = 1.0)
end

function synthetic_solar(copied, weights, period_length, rng; lag = 1, noise = 0.10)
    shaped = _by_day(copied, period_length, rng) do day, r
        circshift(day, lag) .* clamp.(1 .- noise .* abs.(_ar1_noise(r, length(day))), 0.0, 1.0)
    end
    return _rescale_to_total(shaped, weights, sum(copied .* weights); limit = 1.0)
end

"""
    weather_constraint_scores(demand, wind, solar, weights; wind_mw, solar_mw, hours) -> scores

How constraining each library scenario is: average residual load (demand minus wind and solar) over
its `hours` highest residual-load weighted hours. The reference fleet is split across zones by demand
share and is only used for ranking. Arguments are T × Z × K (demand, wind, solar) and T × K (weights).
"""
function weather_constraint_scores(demand, wind, solar, weights; wind_mw = 20_000.0, solar_mw = 17_000.0, hours = 100.0)
    K, Z = size(demand, 3), size(demand, 2)
    scores = zeros(K)
    for k in 1:K
        share = vec(sum(demand[:, :, k], dims = 1)) ./ sum(demand[:, :, k])
        residual = vec(sum(demand[:, :, k], dims = 2)) .-
                   sum(share[z] .* (wind_mw .* wind[:, z, k] .+ solar_mw .* solar[:, z, k]) for z in 1:Z)
        total, covered = 0.0, 0.0
        for t in sortperm(residual; rev = true)
            used = min(weights[t, k], hours - covered)
            total += used * residual[t]
            covered += used
            covered >= hours && break
        end
        scores[k] = total / covered
    end
    return scores
end

# Reads the old scenario files into a library, ordered least to most constraining
function build_weather_library(old_folder, scenarios, Z; synthetic_zone_profiles = true, seed = 1234, expected_most_constraining = nothing)
    weight_table = _read_old(old_folder, "Representative_period_weights.csv")
    n_periods = nrow(weight_table)
    rng = Xoshiro(seed)
    T = nrow(_read_old(old_folder, "Demand_$(first(scenarios)).csv"))
    period_length = T ÷ n_periods
    T == period_length * n_periods || error("$T hours don't divide evenly into $n_periods representative periods.")
    K = length(scenarios)
    demand, wind, solar, weights = zeros(T, Z, K), zeros(T, Z, K), zeros(T, Z, K), zeros(T, K)
    profile_source = fill("source", K)

    for (k, scenario) in enumerate(scenarios)
        demand_file = _read_old(old_folder, "Demand_$scenario.csv")
        availability_file = _read_old(old_folder, "Resources_availability_$scenario.csv")
        weights[:, k] = repeat(Float64.(weight_table[!, "Weight-$scenario"]) ./ period_length, inner = period_length)
        for z in 1:Z
            demand[:, z, k] = Float64.(demand_file[!, "Demand-z$z"])
            wind[:, z, k] = Float64.(availability_file[!, "Wind-z$z"])
            solar[:, z, k] = Float64.(availability_file[!, "Solar-z$z"])
        end
        synthetic_zone_profiles || continue
        for z in 2:Z
            if is_scaled_copy(demand[:, z, k], demand[:, 1, k])
                demand[:, z, k] = synthetic_demand(demand[:, z, k], weights[:, k], period_length, rng)
                profile_source[k] = "synthetic"
            end
            if is_scaled_copy(wind[:, z, k], wind[:, 1, k])
                wind[:, z, k] = synthetic_wind(wind[:, z, k], weights[:, k], period_length, rng)
                profile_source[k] = "synthetic"
            end
            if is_scaled_copy(solar[:, z, k], solar[:, 1, k])
                solar[:, z, k] = synthetic_solar(solar[:, z, k], weights[:, k], period_length, rng)
                profile_source[k] = "synthetic"
            end
        end
    end

    scores = weather_constraint_scores(demand, wind, solar, weights)
    order = sortperm(scores)
    if expected_most_constraining !== nothing && scenarios[last(order)] != expected_most_constraining
        error("Old scenario $(scenarios[last(order)]) ranks most constraining, not $expected_most_constraining as expected. Scores: $(round.(scores, digits = 0))")
    end

    rows = DataFrame(Scenario = Int[], Source_file = Int[], Zone_profiles = String[], Period = Int[], Hour = Int[], Weight = Float64[])
    for z in 1:Z, prefix in ("Demand", "Wind", "Solar")
        rows[!, "$(prefix)_z$z"] = Float64[]
    end
    for (new_index, k) in enumerate(order), t in 1:T
        row = Dict{String,Any}("Scenario" => new_index, "Source_file" => scenarios[k], "Zone_profiles" => profile_source[k],
                               "Period" => (t - 1) ÷ period_length + 1, "Hour" => (t - 1) % period_length + 1,
                               "Weight" => weights[t, k] * period_length)
        for z in 1:Z
            row["Demand_z$z"], row["Wind_z$z"], row["Solar_z$z"] = demand[t, z, k], wind[t, z, k], solar[t, z, k]
        end
        push!(rows, row; cols = :subset)
    end
    return rows, scenarios[order], scores[order]
end

# Technologies.csv: costs from the assumption table, everything site-specific from the old files
function build_technologies(old_folder, technologies, Z)
    resources = _read_old(old_folder, "Resources.csv")
    "Zone" in names(resources) || error("Resources.csv has no Zone column; only multi-zone folders can be converted.")
    derates = _read_old(old_folder, "CRM_Resources.csv")
    base_name(resource) = replace(resource, r"-z\d+$" => "")
    names_old = unique(base_name.(resources.Resources))
    unknown = setdiff(names_old, technologies.Technology)
    isempty(unknown) || error("No cost assumptions for $(join(unknown, ", ")). Add them to the `technologies` table.")

    table = technologies[indexin(names_old, technologies.Technology), :]
    table.CES_eligible = [first(resources.CES_eligible[base_name.(resources.Resources) .== t]) for t in names_old]
    table.Contracted = [first(resources.Contracted[base_name.(resources.Resources) .== t]) for t in names_old]
    table.Capacity_UB = [maximum(resources.Capacity_UB[base_name.(resources.Resources) .== t]) for t in names_old]
    # Outages already reduce thermal availability, so a derate would count them twice
    table.CRM_derate = [ismissing(profile) ? 1.0 : derates.Derating_Factor[findfirst(==("$(t)-z1"), derates.Resource)]
                        for (t, profile) in zip(names_old, table.Profile)]
    for z in 1:Z
        table[!, "Cost_multiplier_z$z"] = [begin
            cost(zone) = resources[resources.Resources .== "$(t)-z$zone", "Investment cost"][1]
            cost(z) / cost(1)
        end for t in names_old]
    end
    return table
end

function _new_settings_block(old_folder, scenarios, Z; gas, counts, discount_rate, demand_low, demand_high)
    crm = _read_old(old_folder, "CRM_System.csv")
    block = Dict{String,Any}(
        "Scenarios" => Dict{String,Any}(
            "Demand" => Dict{String,Any}("Count" => counts["Demand"], "Distribution" => "Uniform", "Low" => demand_low, "High" => demand_high),
            "Gas" => merge(Dict{String,Any}("Count" => counts["Gas"]), gas),
            "Weather" => Dict{String,Any}("Count" => min(counts["Weather"], length(scenarios)))),
        "Discount rate" => discount_rate,
        "CRM zones" => Dict{String,Any}("Margin" => Float64.(crm.Margin[1:Z]), "Penalty" => Float64.(crm.PriceCap[1:Z])))
    return block
end

# Old settings text minus retired keys, keeping comments and layout, plus the new block
function _write_settings(old_folder, new_folder, block)
    kept = filter(readlines(joinpath(old_folder, "settings.yaml"))) do line
        !any(key -> startswith(line, key * ":"), RETIRED_SETTINGS)
    end
    settings_path = joinpath(new_folder, "settings.yaml")
    open(settings_path, "w") do io
        println(io, join(kept, "\n"))
        println(io, "\n# Scenarios: Count points evenly spaced from Low to High, weighted by the Distribution")
        println(io, "# (Uniform, Triangular, Normal or LogNormal). Demand values are fractions of each zone's reference peak; Gas values are \$/MMBtu.")
        YAML.write(io, block)
    end
end

function _settings_difference(base::AbstractDict, other::AbstractDict)
    difference = Dict{String,Any}()
    for (key, value) in other
        key in RETIRED_SETTINGS && continue
        if value isa AbstractDict && get(base, key, nothing) isa AbstractDict
            nested = _settings_difference(base[key], value)
            isempty(nested) || (difference[key] = nested)
        elseif !isequal(get(base, key, nothing), value)
            difference[key] = value
        end
    end
    return difference
end

"""
    convert_input_folder(old_folder, new_folder; kwargs...) -> report

Writes `old_folder`'s data to `new_folder` in the current input layout (see the top of this file) and
returns the library order and an old-against-new cost table. The old folder isn't changed.

Keywords: `technologies`, `fuels`, `storage`, `gas` (assumption tables), `counts` (default scenario
counts written to settings), `discount_rate`, `synthetic_zone_profiles`, `seed`, `profiles`
(Dict of profile name => old folder whose settings differ, e.g. "della"), `expected_most_constraining`
(old scenario that must rank last, or the conversion stops), and `overwrite`.
"""
function convert_input_folder(old_folder::AbstractString, new_folder::AbstractString;
                              technologies = default_technology_assumptions(), fuels = default_fuel_assumptions(),
                              storage = default_storage_assumptions(), gas = default_gas_distribution(),
                              counts = Dict("Demand" => 5, "Gas" => 5, "Weather" => 5), discount_rate = 0.0768,
                              synthetic_zone_profiles::Bool = true, seed::Integer = 1234,
                              profiles::AbstractDict = Dict(), expected_most_constraining = nothing, overwrite::Bool = false)
    scenarios = old_format_scenarios(old_folder)
    if isdir(new_folder) && !isempty(readdir(new_folder))
        overwrite || error("$new_folder already has files. Pass overwrite = true to replace them.")
    end
    mkpath(new_folder)

    network = _read_old(old_folder, "Network.csv")
    Z = count(c -> occursin(r"^z\d+$", c), names(network))
    Z >= 1 || error("Network.csv in $old_folder has no zone columns.")

    library, order, scores = build_weather_library(old_folder, scenarios, Z; synthetic_zone_profiles = synthetic_zone_profiles,
                                                   seed = seed, expected_most_constraining = expected_most_constraining)
    technology_table = build_technologies(old_folder, technologies, Z)

    old_storage = _read_old(old_folder, "Resources_storage.csv")
    storage_table = select(old_storage, Not("Investment cost"))
    if "Index_ID" in names(storage_table)
        select!(storage_table, Not("Index_ID"))
    end
    storage_table.Capex_per_kW .= storage.Capex_per_kW
    storage_table.Fixed_OM_per_kW_yr .= storage.Fixed_OM_per_kW_yr
    storage_table.Life_years .= storage.Life_years

    # The adders are stored in MW; settings hold them as fractions of zone 1's reference peak
    adders = Float64.(collect(_read_old(old_folder, "Demand_adders.csv")[1, 2:end]))
    reference = minimum(maximum(library[library.Scenario .== k, "Demand_z1"]) for k in unique(library.Scenario))
    block = _new_settings_block(old_folder, scenarios, Z; gas = gas, counts = counts, discount_rate = discount_rate,
                                demand_low = round(minimum(adders) / reference, digits = 4), demand_high = round(maximum(adders) / reference, digits = 4))

    CSV.write(joinpath(new_folder, "Weather_library.csv"), library)
    CSV.write(joinpath(new_folder, "Technologies.csv"), technology_table)
    CSV.write(joinpath(new_folder, "Resources_storage.csv"), storage_table)
    CSV.write(joinpath(new_folder, "Fuels.csv"), fuels)
    cp(joinpath(old_folder, "Network.csv"), joinpath(new_folder, "Network.csv"); force = true)
    _write_settings(old_folder, new_folder, block)

    base_settings = YAML.load_file(joinpath(old_folder, "settings.yaml"))
    for (profile, profile_folder) in profiles
        difference = _settings_difference(base_settings, YAML.load_file(joinpath(profile_folder, "settings.yaml")))
        YAML.write_file(joinpath(new_folder, "settings_$(profile).yaml"), difference)
    end

    # Report costs that changed, since old values can't carry over
    resources = _read_old(old_folder, "Resources.csv")
    new_cost = [begin
        t = technology_table[findfirst(==(replace(r, r"-z\d+$" => "")), technology_table.Technology), :]
        annual_investment_cost(t.Capex_per_kW, t.Fixed_OM_per_kW_yr, t.Life_years, discount_rate; multiplier = t[Symbol("Cost_multiplier_z", z)])
    end for (r, z) in zip(resources.Resources, resources.Zone)]
    costs = DataFrame(Resource = resources.Resources, Old_per_MW_yr = Float64.(resources[!, "Investment cost"]), New_per_MW_yr = new_cost)
    @info("Converted $old_folder to $new_folder: $(length(scenarios)) weather scenarios, least to most constraining (old file numbers): $order.")
    @info("Costs now come from the assumption tables, not the old folder:\n$(costs)")
    return (order = order, scores = scores, costs = costs)
end
