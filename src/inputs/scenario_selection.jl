using Distributions, Statistics, Combinatorics

#==========

Chooses the scenarios a run uses. Selection is deterministic, with no random seed.

Gas and demand are value progressions: a distribution limited to Low-High in settings.yaml.
Weather is a library of whole years, so scenarios are picked from it rather than interpolated.

==========#

function _scenario_parameter(config, name::String, key::String)
    haskey(config, key) || error("Scenario distribution $name needs \"$key\" in settings.yaml.")
    return config[key]
end

# The distribution named in a settings block (before limiting it to Low-High)
function scenario_distribution(config)
    name = _scenario_parameter(config, "(unnamed)", "Distribution")
    low, high = config["Low"], config["High"]
    name == "Uniform" && return Uniform(low, high)
    name == "Triangular" && return TriangularDist(low, high, _scenario_parameter(config, name, "Mode"))
    name == "Normal" && return Normal(_scenario_parameter(config, name, "Mean"), _scenario_parameter(config, name, "Sigma"))
    name == "LogNormal" && return LogNormal(log(_scenario_parameter(config, name, "Median")), _scenario_parameter(config, name, "Sigma"))
    error("Unknown scenario distribution \"$name\". Valid options: Uniform, Triangular, Normal, LogNormal.")
end

"""
    select_distribution_scenarios(config) -> (values, probabilities)

Spaces `Count` points evenly in value from `Low` to `High`, so both extremes are always kept. Each
point gets the probability the distribution assigns to values nearest it, normalized over Low-High
so the set sums to 1. `Count` of 1 returns the mean of the limited distribution.
"""
function select_distribution_scenarios(config)
    for key in ("Count", "Low", "High")
        _scenario_parameter(config, "block", key)
    end
    n, low, high = config["Count"], config["Low"], config["High"]
    (n isa Integer && n >= 1) || error("Scenario Count must be a positive integer, got $n.")
    low < high || error("Scenario Low ($low) must be below High ($high).")
    distribution = scenario_distribution(config)
    if n == 1
        return [_limited_mean(distribution, low, high)], [1.0]
    end
    values = collect(range(low, high; length = n))
    edges = vcat(low, (values[1:end-1] .+ values[2:end]) ./ 2, high)
    mass = diff(cdf.(distribution, edges))
    sum(mass) > 0 || error("The distribution puts no probability between Low ($low) and High ($high).")
    return values, mass ./ sum(mass)
end

# Numerical mean of `distribution` limited to low-high. Not every limited distribution has a closed form.
function _limited_mean(distribution, low, high; cells::Int = 20_000)
    edges = range(low, high; length = cells + 1)
    mass = diff(cdf.(Ref(distribution), edges))
    return sum(((edges[1:end-1] .+ edges[2:end]) ./ 2) .* mass) / sum(mass)
end

# Scales each feature by its interquartile range around the median, so one extreme scenario doesn't
# compress the distances between all the others. Rows are scenarios.
function _robust_scale(features::AbstractMatrix)
    scaled = float.(features)
    for j in axes(scaled, 2)
        low, high = quantile(scaled[:, j], (0.25, 0.75))
        # The floor stops a feature that's nearly identical across scenarios from turning rounding into distance
        spread = max(high - low, 0.01 * abs(median(scaled[:, j])), eps())
        scaled[:, j] = (scaled[:, j] .- median(scaled[:, j])) ./ spread
    end
    return scaled
end

"""
    select_weather_scenarios(features, n) -> (indices, probabilities)

Chooses `n` rows of `features` (library scenarios, equally likely). Tries every subset and keeps the
one with the smallest worst-case distance from any scenario to its nearest selected one, breaking
ties on total distance. Each selected scenario gets the share of scenarios assigned to it. `n` of 1
picks the most central scenario instead, since the worst-case rule would favor the outlier.
Indices stay sorted, so the library's least-to-most-constraining order carries through.
"""
function select_weather_scenarios(features::AbstractMatrix, n::Integer)
    K = size(features, 1)
    (1 <= n <= K) || error("Weather Count must be between 1 and $K (the library size), got $n.")
    scaled = _robust_scale(features)
    distance = [sqrt(sum(abs2, scaled[i, :] .- scaled[j, :])) for i in 1:K, j in 1:K]
    if n == 1
        return [argmin(vec(sum(distance, dims = 2)))], [1.0]
    end
    best, best_score = Int[], (Inf, Inf)
    for subset in combinations(1:K, n)
        nearest = minimum(distance[:, subset], dims = 2)
        score = (maximum(nearest), sum(nearest))
        if score < best_score
            best, best_score = subset, score
        end
    end
    assigned = [argmin(distance[i, best]) for i in 1:K]
    return best, [count(==(j), assigned) / K for j in eachindex(best)]
end

# Features that describe how a library scenario stresses the system, including the score that orders
# the library, so a selection spans least to most constraining. Arguments are T × Z × K (demand, wind,
# solar) and T × K (hourly weights).
function weather_features(demand, wind, solar, weights; period_length::Int = 24)
    K = size(demand, 3)
    features = zeros(K, 6)
    constraint_scores = weather_constraint_scores(demand, wind, solar, weights)
    for k in 1:K
        system_demand = vec(sum(demand[:, :, k], dims = 2))
        wind_mean, solar_mean = vec(mean(wind[:, :, k], dims = 2)), vec(mean(solar[:, :, k], dims = 2))
        total_hours = sum(weights[:, k])
        daily_renewable = mean(reshape((wind_mean .+ solar_mean) ./ 2, period_length, :), dims = 1)
        features[k, :] = [maximum(system_demand), sum(system_demand .* weights[:, k]),
                          sum(wind_mean .* weights[:, k]) / total_hours, sum(solar_mean .* weights[:, k]) / total_hours,
                          minimum(daily_renewable), constraint_scores[k]]
    end
    return features
end

# Flat format of everything selected, for Scenario_selection.csv
function scenario_selection_table(selection)
    rows = NamedTuple{(:Uncertainty, :Scenario, :Value, :Unit, :Probability),Tuple{String,Int,Float64,String,Float64}}[]
    for (name, unit) in (("Demand", "fraction of zone reference peak"), ("Gas", "\$/MMBtu"), ("Weather", "library scenario"))
        for (i, (value, probability)) in enumerate(zip(selection[name]["Values"], selection[name]["Probabilities"]))
            push!(rows, (Uncertainty = name, Scenario = i, Value = value, Unit = unit, Probability = probability))
        end
    end
    return DataFrame(rows)
end

function write_scenario_selection(inputs::Dict, folder::String)
    mkpath(folder)
    CSV.write(joinpath(folder, "Scenario_selection.csv"), scenario_selection_table(inputs["Scenario selection"]))
end
