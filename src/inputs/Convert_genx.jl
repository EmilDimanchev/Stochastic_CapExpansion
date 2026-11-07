#==========

Converts a GenX 0.3.7 case folder (Generators_data.csv, Generators_variability.csv, Load_data.csv,
Network.csv, Capacity_reserve_margin.csv) into the layout load_input_data reads. See
Conversion_notes.md in the output folder for what a converted case keeps and drops.

Usage, from the repo root after including experiments/MGA_Tests.jl:

  convert_genx_case("path/to/genx_case", "inputs/New_folder")

Resources group into technology rows (one resource per zone), which keeps existing fleets, new-build
options and zone costs apart:

  fixed fleet (New_Build = -1)      Capacity_LB = Capacity_UB = existing MW, paying fixed O&M
  retirable fleet (New_Build = 0)   Capacity_UB = existing MW, paying fixed O&M only on what stays
  new build (New_Build = 1)         Capacity_UB = build limit, paying the case's annual cost
  wind and solar clusters           tiers by levelized cost, so the supply curve keeps its shape

Annual costs reproduce the case's $/MW-yr exactly, which embed its 2.6-3.6% real WACC. Capex_per_kW is
therefore the capex that gives that annual cost at the settings' Discount rate, not a quoted capex.

The source has one weather year, so other library years are synthetic: each week is drawn from nearby
weeks of the source year, keeping every zone's demand and resource together, then wind and solar are
scaled by an annual anomaly. Demand keeps each zone's annual energy, since the Demand scenarios cover
that. Each year is then clustered down to representative weeks.

==========#

const GENX_UNLIMITED_MW = 500_000.0 # the model needs a finite bound; used where the case has none

# GenX technology => model technology. `nothing` has no equivalent in the model and is dropped.
const GENX_TECHNOLOGIES = Dict{String,Union{Nothing,String}}(
    "Conventional Hydroelectric" => "Hydro", "Small Hydroelectric" => "Hydro", "Biomass" => "Biomass",
    "Geothermal" => "Geothermal", "distributed_generation" => "Distributed",
    "Conventional Steam Coal" => "Coal existing", "Nuclear" => "Nuclear existing", "Petroleum Liquids" => "Oil existing",
    "Natural Gas Fired Combined Cycle" => "Gas CC existing", "Natural Gas Fired Combustion Turbine" => "Gas CT existing",
    "Natural Gas Internal Combustion Engine" => "Gas CT existing", "Natural Gas Steam Turbine" => "Gas ST existing",
    "Onshore Wind Turbine" => "Wind existing", "Solar Photovoltaic" => "Solar existing", "Offshore Wind Turbine" => "Offshore existing",
    "NaturalGas_CCAvgCF_Moderate" => "Gas CC new", "NaturalGas_CTAvgCF_Moderate" => "Gas CT new",
    "NaturalGas_CCCCSAvgCF_Conservative" => "Gas CCS new",
    "LandbasedWind_Class3_Moderate_" => "Wind", "UtilityPV_Class1_Moderate_" => "Solar",
    "OffShoreWind_Class3_Moderate_fixed_1" => "Offshore fixed", "OffShoreWind_Class12_Moderate_floating_1" => "Offshore floating",
    "Batteries" => "Battery existing", "Hydroelectric Pumped Storage" => "Pumped hydro", "Battery_*_Moderate" => "Battery",
    "Coal_CCS90AvgCF_Moderate" => nothing, "Nuclear_Nuclear_Moderate" => nothing, # listed, but no capacity and not buildable
    "trans_light_duty" => nothing, "comm_water_heat" => nothing, "res_water_heat" => nothing,
    "comm_space_heat_cool" => nothing, "res_space_heat_cool" => nothing) # flexible demand

# GenX has no outages here (thermal availability is 1), and the model represents them through rates.
# Same rates as default_technology_assumptions; steam turbine, oil and biomass are estimates. (forced, planned)
# Hydro's 0.20 is the fleet's unavailable share from all causes (about 80% availability, ORNL/DOE U.S.
# Hydropower Market Report, 2023 edition), so it isn't split into forced and planned.
const GENX_OUTAGE_RATES = Dict("Hydro" => (0.20, 0.0), "Gas CC" => (0.05, 0.05), "Gas CCS" => (0.05, 0.05), "Gas CT" => (0.08, 0.03), "Gas ST" => (0.08, 0.05),
                               "Coal" => (0.12, 0.06), "Nuclear" => (0.02, 0.07), "Oil" => (0.08, 0.03), "Biomass" => (0.05, 0.05))

const GENX_STORAGE_GROUPS = ("Battery", "Battery existing", "Pumped hydro")

# Share of full-capacity energy available per representative period. Hydro's is the 2016-2025 mean of EIA's
# fleet capacity factor (Electric Power Monthly, Table 6.07.B), which stands in for water and operating limits.
const GENX_ENERGY_LIMITS = Dict("Hydro" => 0.38)

_genx_fuel(fuel::AbstractString) = begin
    token = last(split(fuel, "_"))
    get(Dict("naturalgas" => "Gas", "coal" => "Coal", "uranium" => "Uranium", "distillate" => "Oil", "biomass" => "Biomass"), token, missing)
end

_weighted_mean(values, weights) = sum(values .* weights) / sum(weights)

function _genx_profile_matrix(variability, resources, weights)
    profile = zeros(nrow(variability))
    for (resource, weight) in zip(resources, weights)
        profile .+= weight .* variability[!, resource]
    end
    return profile ./ sum(weights)
end

# Where a resource with no usable profile takes one from, in order of preference. Distributed generation
# is assumed to be solar: the case gives it no fuel, cost or real profile to say otherwise.
const GENX_PROFILE_DONORS = Dict("Wind" => ["Wind"], "Wind existing" => ["Wind"], "Solar" => ["Solar"], "Solar existing" => ["Solar"],
                                 "Offshore fixed" => ["Offshore fixed", "Offshore floating"], "Offshore floating" => ["Offshore floating", "Offshore fixed"],
                                 "Offshore existing" => ["Offshore fixed", "Offshore floating"], "Distributed" => ["Solar"])

"""
    repair_flat_availability!(generators, variability) -> repairs

Wind, solar and distributed availability of exactly 1 in every hour isn't a real profile. Replaces each such column
with the build-limit-weighted mean of real profiles: the same zone first, then any zone, trying each
donor group in GENX_PROFILE_DONORS order. Returns one row per repair. Other resources are left alone,
since a flat profile can be genuine for them (biomass, geothermal) and there's no donor for the rest.
"""
function repair_flat_availability!(generators, variability)
    weight(r) = r.Max_Cap_MW > 0 ? r.Max_Cap_MW : GENX_UNLIMITED_MW
    candidates = filter(r -> !ismissing(r.Group) && haskey(GENX_PROFILE_DONORS, r.Group), generators)
    is_flat = Dict(r => all(==(1), variability[!, r]) for r in candidates.Resource)
    repairs = DataFrame(Resource = String[], Group = String[], Zone = Int[], Donor_group = String[], Donor_zone = String[], Mean_availability = Float64[])
    for flat in eachrow(candidates)
        is_flat[flat.Resource] || continue
        donors, donor_group, donor_zone = nothing, "", ""
        for same_zone in (true, false), group in GENX_PROFILE_DONORS[flat.Group]
            found = filter(r -> r.Group == group && !is_flat[r.Resource] && (!same_zone || r.Zone == flat.Zone), candidates)
            if nrow(found) > 0
                donors, donor_group, donor_zone = found, group, same_zone ? string(flat.Zone) : "all"
                break
            end
        end
        donors === nothing && error("$(flat.Resource) has a flat profile and no resource to take a real one from.")
        profile = _genx_profile_matrix(variability, donors.Resource, weight.(eachrow(donors)))
        variability[!, flat.Resource] = profile
        push!(repairs, (flat.Resource, flat.Group, flat.Zone, donor_group, donor_zone, mean(profile)))
    end
    return repairs
end

# Assigns a zone's wind or solar clusters to tiers by levelized cost, in steps of the zone's peak
# demand: a cluster goes where the middle of its capacity falls.
function genx_cost_tiers(levelized_cost, capacity, zone_peak; tier_limits = (1.0, 3.0))
    tiers = zeros(Int, length(levelized_cost))
    cumulative = 0.0
    for i in sortperm(levelized_cost)
        midpoint = (cumulative + capacity[i] / 2) / zone_peak
        tiers[i] = something(findfirst(>(midpoint), tier_limits), length(tier_limits) + 1)
        cumulative += capacity[i]
    end
    return tiers
end

"""
    select_representative_weeks(features, n; forced = Int[]) -> (weeks, counts)

Chooses `n` rows of `features` (weeks, already scaled) by k-medoids and returns them in order, with how
many weeks each stands for. `forced` weeks are always kept, so extremes survive. Deterministic.
"""
function select_representative_weeks(features::AbstractMatrix, n::Integer; forced = Int[])
    N = size(features, 1)
    (1 <= n <= N) || error("Representative weeks must be between 1 and $N, got $n.")
    length(unique(forced)) <= n || error("$(length(unique(forced))) forced weeks don't fit in $n representative weeks.")
    distance = [sqrt(sum(abs2, features[i, :] .- features[j, :])) for i in 1:N, j in 1:N]
    medoids = isempty(forced) ? [argmin(vec(sum(distance, dims = 2)))] : unique(forced)
    n_fixed = isempty(forced) ? 0 : length(medoids)
    while length(medoids) < n # farthest-first start
        push!(medoids, argmax(vec(minimum(distance[:, medoids], dims = 2))))
    end
    for _ in 1:100
        assignment = [argmin(distance[i, medoids]) for i in 1:N]
        updated = copy(medoids)
        for c in n_fixed+1:n
            members = findall(==(c), assignment)
            updated[c] = members[argmin(vec(sum(distance[members, members], dims = 2)))]
        end
        updated == medoids && break
        medoids = updated
    end
    assignment = [argmin(distance[i, medoids]) for i in 1:N]
    order = sortperm(medoids)
    return medoids[order], [count(==(c), assignment) for c in order]
end

# Weekly features for clustering: demand level and peak, mean wind and solar per zone, each z-scored
function _week_features(demand, wind, solar, period_length)
    weekly(f, data) = [f(view(data, (w-1)*period_length+1:w*period_length, z)) for w in 1:size(data, 1)÷period_length, z in 1:size(data, 2)]
    raw = hcat(weekly(mean, demand), weekly(maximum, demand), weekly(mean, wind), weekly(mean, solar))
    keep = [std(raw[:, j]) > 1e-9 for j in axes(raw, 2)]
    return mapslices(c -> (c .- mean(c)) ./ std(c), raw[:, keep], dims = 1)
end

# Weeks that must survive clustering: highest system demand, highest residual load, lowest renewables
function _extreme_weeks(demand, wind, solar, period_length)
    n_weeks = size(demand, 1) ÷ period_length
    weeks(t) = (t - 1) ÷ period_length + 1
    system = vec(sum(demand, dims = 2))
    share = vec(sum(demand, dims = 1)) ./ sum(demand)
    renewable = (wind * share .+ solar * share) ./ 2
    residual = system .* (1 .- renewable)
    weekly_renewable = [mean(renewable[(w-1)*period_length+1:w*period_length]) for w in 1:n_weeks]
    return unique([weeks(argmax(system)), weeks(argmax(residual)), argmin(weekly_renewable)])
end

# One synthetic year: each week comes from within `window` weeks of itself in the source year
function _synthetic_year(demand, profiles, period_length, rng; window = 3, wind_sigma = 0.06, solar_sigma = 0.03)
    n_weeks = size(demand, 1) ÷ period_length
    source_week = [mod1(w + rand(rng, -window:window), n_weeks) for w in 1:n_weeks]
    hours = reduce(vcat, [collect((w-1)*period_length+1:w*period_length) for w in source_week])
    wind_factor, solar_factor = 1 + wind_sigma * randn(rng), 1 + solar_sigma * randn(rng)
    new_demand = demand[hours, :]
    new_demand .*= (sum(demand, dims = 1) ./ sum(new_demand, dims = 1)) # keep each zone's annual energy
    factor(name) = occursin("Wind", name) || occursin("Offshore", name) ? wind_factor : occursin("Solar", name) || name == "Distributed" ? solar_factor : 1.0
    new_profiles = Dict(name => clamp.(p[hours, :] .* factor(name), 0.0, 1.0) for (name, p) in profiles)
    return new_demand, new_profiles
end

# Builds the Weather_library table, least to most constraining
function build_genx_library(demand, profiles, period_length; rep_weeks, weather_years, seed)
    rng = Xoshiro(seed)
    years = Any[(demand, profiles)]
    for _ in 2:weather_years
        push!(years, _synthetic_year(demand, profiles, period_length, rng))
    end
    Z = size(demand, 2)
    T, K = size(demand, 1), length(years)
    peak = maximum(sum(demand, dims = 2))
    scores = weather_constraint_scores(cat([y[1] for y in years]...; dims = 3), cat([y[2]["Wind"] for y in years]...; dims = 3),
                                       cat([y[2]["Solar"] for y in years]...; dims = 3), ones(T, K);
                                       wind_mw = 0.3 * peak, solar_mw = 0.3 * peak)
    names_sorted = sort!(collect(keys(profiles)))
    tables = DataFrame[]
    selected_weeks = Vector{Int}[]
    for (new_index, k) in enumerate(sortperm(scores))
        year_demand, year_profiles = years[k]
        wind, solar = year_profiles["Wind"], year_profiles["Solar"]
        features = _week_features(year_demand, wind, solar, period_length)
        weeks, counts = select_representative_weeks(features, rep_weeks; forced = _extreme_weeks(year_demand, wind, solar, period_length))
        push!(selected_weeks, weeks)
        hours = reduce(vcat, [collect((w-1)*period_length+1:w*period_length) for w in weeks])
        weights = repeat(Float64.(counts), inner = period_length)
        table = DataFrame(Scenario = new_index, Source_file = 1, Zone_profiles = k == 1 ? "source" : "synthetic",
                          Period = repeat(1:length(weeks), inner = period_length), Hour = repeat(1:period_length, outer = length(weeks)),
                          Weight = repeat(Float64.(counts) .* period_length, inner = period_length))
        # Clustering shifts annual totals a little, so each zone's demand and resource energy is put back
        for z in 1:Z
            table[!, "Demand_z$z"] = round.(_rescale_to_total(year_demand[hours, z], weights, sum(year_demand[:, z])), digits = 1)
            for name in names_sorted
                profile = year_profiles[name]
                table[!, "$(name)_z$z"] = round.(_rescale_to_total(profile[hours, z], weights, sum(profile[:, z]); limit = 1.0), digits = 4)
            end
        end
        push!(tables, table)
    end
    return reduce(vcat, tables), scores[sortperm(scores)], selected_weeks
end

# Technologies.csv and the hourly profiles it names, from the generator table
function build_genx_technologies(generators, variability, zone_peaks, discount_rate; tier_limits = (1.0, 3.0))
    Z = length(zone_peaks)
    rows = filter(r -> !ismissing(r.Group) && !(r.Group in GENX_STORAGE_GROUPS), generators)
    # Split wind and solar build options into cost tiers within each zone
    rows.Tier = fill(0, nrow(rows))
    build_limit(r) = r.Max_Cap_MW > 0 ? r.Max_Cap_MW : GENX_UNLIMITED_MW
    for group in ("Wind", "Solar"), z in 1:Z
        members = findall(r -> r.Group == group && r.Zone == z, eachrow(rows))
        isempty(members) && continue
        limits = build_limit.(eachrow(rows[members, :]))
        hourly_mean = [mean(variability[!, r]) for r in rows.Resource[members]]
        annual = rows.Inv_Cost_per_MWyr[members] .+ rows.Fixed_OM_Cost_per_MWyr[members]
        rows.Tier[members] = genx_cost_tiers(annual ./ max.(hourly_mean, 1e-3), limits, zone_peaks[z]; tier_limits = tier_limits)
    end
    rows.Technology = [r.Tier > 0 ? "$(r.Group) T$(r.Tier)" : r.Group for r in eachrow(rows)]

    # Wind and Solar are always written, for weather features: the zone's build options, else what it has
    profiles = Dict{String,Matrix{Float64}}()
    technology_rows = DataFrame()
    for technology_rows_group in groupby(rows, :Technology, sort = true)
        name = first(technology_rows_group.Technology)
        zone_rows = [technology_rows_group[technology_rows_group.Zone .== z, :] for z in 1:Z]
        is_new = first(technology_rows_group.New_Build) == 1
        weight(frame) = is_new ? build_limit.(eachrow(frame)) : frame.Existing_Cap_MW
        offered = [nrow(f) > 0 && sum(weight(f)) > 0 for f in zone_rows]
        any(offered) || continue
        zones = findall(offered)

        annual_inv = [_weighted_mean(zone_rows[z].Inv_Cost_per_MWyr, weight(zone_rows[z])) for z in zones]
        annual_fom = [_weighted_mean(zone_rows[z].Fixed_OM_Cost_per_MWyr, weight(zone_rows[z])) for z in zones]
        life = is_new ? first(technology_rows_group.cap_recovery_years) : 30.0
        life = life > 0 ? life : 30.0
        capex = is_new ? mean(annual_inv) / capital_recovery_factor(discount_rate, life) / 1000 : 0.0
        fom = mean(annual_fom) / 1000
        reference = 1000 * (capex * capital_recovery_factor(discount_rate, life) + fom)
        multiplier = Vector{Union{Missing,Float64}}(missing, Z)
        multiplier[zones] = reference > 0 ? (annual_inv .+ annual_fom) ./ reference : fill(1.0, length(zones))

        # Profile: blank when availability is always 1, in which case outage rates apply
        profile_name = replace(name, " " => "_")
        zone_profiles = zeros(nrow(variability), Z)
        for z in zones
            zone_profiles[:, z] = _genx_profile_matrix(variability, zone_rows[z].Resource, weight(zone_rows[z]))
        end
        has_profile = any(zone_profiles[:, zones] .< 1 - 1e-9)
        has_profile && (profiles[profile_name] = zone_profiles)

        all_weights = vcat([weight(zone_rows[z]) for z in zones]...)
        all_rows = reduce(vcat, [zone_rows[z] for z in zones])
        fuel = _genx_fuel(first(technology_rows_group.Fuel))
        capture = _weighted_mean(all_rows.CO2_Capture_Rate, all_weights)
        heat_rate = ismissing(fuel) ? missing : _weighted_mean(all_rows.Heat_Rate_MMBTU_per_MWh, all_weights)
        carbon = ismissing(fuel) ? 0.0 : Dict("Gas" => 0.053, "Coal" => 0.095, "Uranium" => 0.0, "Oil" => 0.073, "Biomass" => 0.0)[fuel]
        capture_cost = ismissing(fuel) ? 0.0 : heat_rate * carbon * capture * _weighted_mean(all_rows.CO2_Capture_Cost_per_Metric_Ton, all_weights)
        forced, planned = has_profile ? (0.0, 0.0) : get(GENX_OUTAGE_RATES, replace(name, r" (existing|new)$" => ""), (0.0, 0.0))

        row = Dict{String,Any}("Technology" => name, "Capex_per_kW" => capex, "Fixed_OM_per_kW_yr" => fom, "Life_years" => life,
                               "Fuel" => fuel, "Heat_rate_MMBtu_per_MWh" => heat_rate,
                               "Variable_OM_per_MWh" => _weighted_mean(all_rows.Var_OM_Cost_per_MWh, all_weights) + capture_cost,
                               "Capture_rate" => capture, "Energy_limit_CF" => get(GENX_ENERGY_LIMITS, name, missing), "Profile" => has_profile ? profile_name : missing,
                               "Forced_outage_rate" => forced, "Planned_outage_rate" => planned,
                               "CES_eligible" => ismissing(fuel) || fuel in ("Uranium", "Biomass") || capture > 0 ? 1 : 0, "Contracted" => 1,
                               "CRM_derate" => has_profile ? _weighted_mean(all_rows.CapRes_1, all_weights) : 1.0)
        for z in 1:Z
            row["Cost_multiplier_z$z"] = multiplier[z]
            is_offered = offered[z]
            must_stay = is_offered && first(zone_rows[z].New_Build) == -1
            total = is_offered ? sum(weight(zone_rows[z])) : 0.0
            row["Capacity_UB_z$z"] = is_offered ? min(total, GENX_UNLIMITED_MW) : missing
            row["Capacity_LB_z$z"] = must_stay ? total : (is_offered ? 0.0 : missing)
        end
        order = vcat(["Technology", "Capex_per_kW", "Fixed_OM_per_kW_yr", "Life_years", "Fuel", "Heat_rate_MMBtu_per_MWh", "Variable_OM_per_MWh",
                      "Capture_rate", "Energy_limit_CF", "Profile", "Forced_outage_rate", "Planned_outage_rate", "CES_eligible", "Contracted", "CRM_derate"],
                     ["Cost_multiplier_z$z" for z in 1:Z], ["Capacity_LB_z$z" for z in 1:Z], ["Capacity_UB_z$z" for z in 1:Z])
        technology_rows = vcat(technology_rows, DataFrame([key => [row[key]] for key in order]); cols = :union)
    end

    # Weather features need Wind and Solar for every zone
    for (base, group, existing) in (("Wind", "Wind", "Wind existing"), ("Solar", "Solar", "Solar existing"))
        base_profile = zeros(nrow(variability), Z)
        for z in 1:Z
            for candidate in (group, existing)
                members = filter(r -> r.Group == candidate && r.Zone == z, rows)
                if nrow(members) > 0 && sum(candidate == group ? build_limit.(eachrow(members)) : members.Existing_Cap_MW) > 0
                    base_profile[:, z] = _genx_profile_matrix(variability, members.Resource, candidate == group ? build_limit.(eachrow(members)) : members.Existing_Cap_MW)
                    break
                end
            end
        end
        profiles[base] = base_profile
    end
    return technology_rows, profiles
end

# Resources_storage.csv: one row per storage kind per zone
function build_genx_storage(generators, discount_rate; new_duration_hours = 4.0)
    table = DataFrame()
    storage = filter(r -> !ismissing(r.Group) && r.Group in GENX_STORAGE_GROUPS, generators)
    for group in groupby(storage, [:Group, :Zone], sort = true)
        name, zone = first(group.Group), first(group.Zone)
        is_new = name == "Battery"
        weight = is_new ? ones(nrow(group)) : group.Existing_Cap_MW
        sum(weight) > 0 || continue
        duration = is_new ? new_duration_hours : sum(group.Existing_Cap_MWh) / sum(group.Existing_Cap_MW)
        life = is_new ? first(group.cap_recovery_years) : 30.0
        annual_inv = is_new ? _weighted_mean(group.Inv_Cost_per_MWyr .+ duration .* group.Inv_Cost_per_MWhyr, weight) : 0.0
        fom = _weighted_mean(group.Fixed_OM_Cost_per_MWyr .+ duration .* group.Fixed_OM_Cost_per_MWhyr, weight)
        capacity = is_new ? GENX_UNLIMITED_MW : sum(group.Existing_Cap_MW)
        append!(table, DataFrame(Resources = "$(name)-z$(zone)", Zone = zone, Storage = 0, Power_to_energy_ratio = 1 / duration,
                                 Loss_discharge = _weighted_mean(group.Eff_Down, weight), Loss_charge = _weighted_mean(group.Eff_Up, weight),
                                 Capacity_LB = is_new ? 0.0 : capacity, Capacity_UB = capacity,
                                 Capex_per_kW = annual_inv / capital_recovery_factor(discount_rate, life) / 1000,
                                 Fixed_OM_per_kW_yr = fom / 1000, Life_years = life))
    end
    return table
end

function _genx_network(old_network, Z)
    table = DataFrame(Lines = old_network.Network_Lines)
    for z in 1:Z
        table[!, "z$z"] = old_network[!, "z$z"]
    end
    table.Existing_capacity = Float64.(old_network.Line_Max_Flow_MW)
    table.Max_expansion = Float64.(old_network.Line_Max_Reinforcement_MW)
    table.CAPEX_per_MW = Float64.(old_network.Line_Reinforcement_Cost_per_MWyr) # already annual, as the master problem expects
    return table
end

# Settings from a current folder's settings.yaml, with the case's values for what the case sets
function _write_genx_settings(template_folder, new_folder, block, nse_costs, nse_shares)
    lines = readlines(joinpath(template_folder, "settings.yaml"))
    cut = findfirst(line -> startswith(line, "# Scenarios:"), lines)
    kept = cut === nothing ? lines : lines[1:cut-1]
    replaced = map(kept) do line
        startswith(line, "CRM flag:") ? "CRM flag: false # off in the GenX case; margins below are its values" :
        startswith(line, "Cost of nonserved energy:") ? "Cost of nonserved energy: [$(join(nse_costs, ","))]" :
        startswith(line, "Max nse:") ? "Max nse: [$(join(nse_shares, ","))]" : line
    end
    open(joinpath(new_folder, "settings.yaml"), "w") do io
        println(io, join(replaced, "\n"))
        println(io, "# Scenarios: Count points evenly spaced from Low to High, weighted by the Distribution")
        println(io, "# (Uniform, Triangular, Normal or LogNormal). Demand values are fractions of each zone's reference peak; Gas values are \$/MMBtu.")
        YAML.write(io, block)
    end
    della = joinpath(template_folder, "settings_della.yaml")
    isfile(della) && cp(della, joinpath(new_folder, "settings_della.yaml"); force = true)
end

function _write_genx_notes(new_folder, source, rep_weeks, weather_years, scores, dropped, repairs)
    repair_counts = join(["$(nrow(g)) $(first(g.Group))" for g in groupby(repairs, :Group, sort = true)], ", ")
    borrowed = count(==("all"), repairs.Donor_zone)
    open(joinpath(new_folder, "Conversion_notes.md"), "w") do io
        println(io, """
        # Converted from GenX 0.3.7

        Source: `$(source)`. Built by `convert_genx_case` (`src/inputs/Convert_genx.jl`).

        ## Representation
        - $(rep_weeks) representative weeks per weather year, including peak demand, peak residual load and lowest renewable weeks. Each zone's annual energy is restored after clustering.
        - $(weather_years) weather years: scenario labels in `Zone_profiles` show which is the source and which are synthetic. Constraint scores, least to most: $(round.(scores, digits = 0)).
        - Synthetic years are block-resampled from the one source year with annual wind and solar anomalies. They aren't observed weather.

        ## Assumptions added
        - Fuel: repo gas distribution, coal 1.47 and uranium 0.85 replace the case's prices (about \$14.50 gas, \$21 coal). Oil (\$38.82) and biomass (\$5) keep the case's values.
        - Thermal outages: the case's availability is 1, so rates come from `GENX_OUTAGE_RATES`.
        - New batteries are fixed at 4 hours (the case allowed 1 to 10), with energy costs included. CRM penalty is the repo's 5000.

        ## Kept from the case
        Annual investment and fixed O&M by zone, heat rates and variable O&M (capacity-weighted across zones per technology), existing fleets, build limits, transmission, demand, availability profiles, NSE costs, CRM margins (flag off, as in the case).

        ## Availability repaired
        The source has wind, solar and distributed resources at exactly 1.0 in every hour, most of them absent from the top-level Generators_variability.csv. These $(nrow(repairs)) resources ($(repair_counts)) now use the mean of real profiles in their zone, and $(borrowed) of them, with none in the zone, use a mean across zones. Fixed offshore mostly borrows floating offshore profiles, since few fixed clusters have real ones. `Availability_repairs.csv` lists each one. They're estimates, and wind and solar tiers are ranked on the repaired values.

        Distributed generation is assumed to be solar and uses its zone's new-build solar profile.

        ## Hydro limits
        The source has hydro flat at 1.0. It's now capped at 80% of capacity in any hour (fleet availability, ORNL/DOE U.S. Hydropower Market Report, 2023 edition) and at 38% of full-capacity energy in each representative week (2016-2025 mean capacity factor, EIA Electric Power Monthly Table 6.07.B). The same values apply to every zone and season, and to conventional and small hydro alike, so this has no seasonal or regional shape. The energy limit uses a dispatch result as a stand-in for water and operating limits.

        ## Not represented
        - Dropped resources: $(join(dropped, ", ")).
        - Unit commitment, ramping and minimum power, line losses, zone-specific fuel prices, storage variable O&M, and moving hydro energy between weeks.
        - Build-limit policy: the case's inactive minimum capacity and CO2 inputs aren't carried over.
        """)
    end
end

"""
    convert_genx_case(genx_folder, new_folder; kwargs...) -> report

Writes a GenX 0.3.7 case to `new_folder` in the current input layout (see the top of this file). The
source folder isn't changed. Scenario assumptions are keywords: `gas`, `fuels`, `demand_low`,
`demand_high`, `counts`. Others: `rep_weeks`, `weather_years`, `period_length`, `seed`,
`discount_rate`, `tier_limits` (multiples of zone peak), `template_folder` (settings to build on), and
`overwrite`. Returns the technology and storage tables, week selection, constraint scores and the
availability repairs (see repair_flat_availability!).
"""
function convert_genx_case(genx_folder::AbstractString, new_folder::AbstractString;
                           rep_weeks::Integer = 12, weather_years::Integer = 5, period_length::Integer = 168, seed::Integer = 1234,
                           discount_rate = 0.0768, gas = default_gas_distribution(), demand_low = -0.15, demand_high = 0.15,
                           counts = Dict("Demand" => 5, "Gas" => 5, "Weather" => 5), tier_limits = (1.0, 3.0),
                           fuels = DataFrame(Fuel = ["Gas", "Coal", "Uranium", "Oil", "Biomass"], Carbon_t_per_MMBtu = [0.053, 0.095, 0.0, 0.073, 0.0],
                                             Price_per_MMBtu = [missing, 1.47, 0.85, 38.82, 5.0]),
                           template_folder::AbstractString = joinpath("inputs", "Inputs_2z_7tech"), overwrite::Bool = false)
    for file in ("Generators_data.csv", "Generators_variability.csv", "Load_data.csv", "Network.csv", "Capacity_reserve_margin.csv")
        isfile(joinpath(genx_folder, file)) || error("$genx_folder isn't a GenX case this can convert; missing $file.")
    end
    if isdir(new_folder) && !isempty(readdir(new_folder))
        overwrite || error("$new_folder already has files. Pass overwrite = true to replace them.")
    end
    mkpath(new_folder)

    read_genx(file) = CSV.read(joinpath(genx_folder, file), DataFrame; stringtype = String)
    generators = read_genx("Generators_data.csv")
    unknown = setdiff(unique(generators.technology), keys(GENX_TECHNOLOGIES))
    isempty(unknown) || error("No model technology for GenX technologies: $(join(unknown, ", ")). Add them to GENX_TECHNOLOGIES.")
    generators.Group = [something(GENX_TECHNOLOGIES[t], missing) for t in generators.technology]
    dropped = sort(unique(generators.technology[ismissing.(generators.Group)]))

    # The time-domain-reduced files are what GenX ran, and can hold resources the top-level ones lack.
    # They're only usable here if the periods are still in calendar order, since years resample nearby weeks.
    tdr = joinpath(genx_folder, "TDR_Results")
    time_folder = isdir(tdr) ? joinpath("TDR_Results") : ""
    if isdir(tdr)
        period_map = CSV.read(joinpath(tdr, "Period_map.csv"), DataFrame)
        period_map.Rep_Period_Index == period_map.Period_Index || error("$tdr has clustered periods. convert_genx_case needs every period, in order.")
    end
    variability = read_genx(joinpath(time_folder, "Generators_variability.csv"))
    load = read_genx(joinpath(time_folder, "Load_data.csv"))
    missing_profiles = setdiff(generators.Resource, names(variability))
    isempty(missing_profiles) || error("No availability profile for: $(join(first(missing_profiles, 5), ", ")) ($(length(missing_profiles)) in all).")
    network = read_genx("Network.csv")
    Z = count(c -> occursin(r"^z\d+$", c), names(network))
    demand = reduce(hcat, [Float64.(load[!, "Load_MW_z$z"]) for z in 1:Z])
    size(demand, 1) % period_length == 0 || error("$(size(demand, 1)) hours don't divide into periods of $period_length.")
    zone_peaks = vec(maximum(demand, dims = 1))
    repairs = repair_flat_availability!(generators, variability)

    technologies, profiles = build_genx_technologies(generators, variability, zone_peaks, discount_rate; tier_limits = tier_limits)
    storage = build_genx_storage(generators, discount_rate)
    library, scores, weeks = build_genx_library(demand, profiles, period_length; rep_weeks = rep_weeks, weather_years = weather_years, seed = seed)

    segments = load[.!ismissing.(load.Demand_Segment), :]
    nse_costs = Float64.(segments[!, "\$/MWh"])
    shares = Float64.(segments.Max_Demand_Curtailment)
    nse_shares = vcat(1 - sum(shares[2:end]), shares[2:end]) # the model's segments are shares summing to all demand

    margins = read_genx("Capacity_reserve_margin.csv").CapRes_1
    block = Dict{String,Any}(
        "Discount rate" => discount_rate,
        "Scenarios" => Dict{String,Any}(
            "Demand" => Dict{String,Any}("Count" => counts["Demand"], "Distribution" => "Uniform", "Low" => demand_low, "High" => demand_high),
            "Gas" => merge(Dict{String,Any}("Count" => counts["Gas"]), gas),
            "Weather" => Dict{String,Any}("Count" => min(counts["Weather"], weather_years))),
        "CRM zones" => Dict{String,Any}("Margin" => Float64.(margins[1:Z]), "Penalty" => fill(5000.0, Z)))

    CSV.write(joinpath(new_folder, "Weather_library.csv"), library)
    CSV.write(joinpath(new_folder, "Technologies.csv"), technologies)
    CSV.write(joinpath(new_folder, "Resources_storage.csv"), storage)
    CSV.write(joinpath(new_folder, "Fuels.csv"), fuels)
    CSV.write(joinpath(new_folder, "Network.csv"), _genx_network(network, Z))
    _write_genx_settings(template_folder, new_folder, block, nse_costs, nse_shares)
    CSV.write(joinpath(new_folder, "Availability_repairs.csv"), repairs)
    _write_genx_notes(new_folder, genx_folder, rep_weeks, weather_years, scores, dropped, repairs)
    @info("Converted $genx_folder to $new_folder: $(nrow(technologies)) technologies, $(nrow(storage)) storage resources, $weather_years weather years, $rep_weeks weeks each.")
    return (technologies = technologies, storage = storage, weeks = weeks, scores = scores, dropped = dropped, repairs = repairs)
end
