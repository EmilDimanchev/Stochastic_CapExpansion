#-------------------------------------------------------------------------------------------
# Settings groups for run_stochastic_exploration_risk_pareto.
#
# Each group is a block in settings.yaml - the base values for that machine/input set, with the
# meaning of every key documented next to it there:
#   "Budget settings", "Cap settings", "Mapping settings", "Sampling settings",
#   "Run settings", "Cut settings"
# A run overrides individual keys through the matching keyword argument of the driver
# (budget_settings = Dict("Tighten budget" => true), cut_settings = Dict("Max rounds" => 10), ...).
# Unknown keys error out, so a typo can't silently fall back to the base value. A nested
# block (e.g. Cut settings' "Prune window", Sampling settings' "GridFill") can be overridden
# partially; the rest is kept.
#-------------------------------------------------------------------------------------------

const SETTINGS_GROUPS = ("Budget settings", "Cap settings", "Mapping settings", "Sampling settings", "Run settings", "Cut settings")

# Base block settings[group] with `overrides` applied.
function resolve_settings(settings::AbstractDict, group::String, overrides::AbstractDict = Dict())
    haskey(settings, group) || error("settings.yaml has no \"$group\" block.")
    base = Dict{String,Any}(settings[group])
    unknown = setdiff(keys(overrides), keys(base))
    isempty(unknown) || error("Unknown key(s) $(collect(unknown)) in $group. Valid keys: $(sort(collect(keys(base)))).")
    for (k, v) in overrides
        if base[k] isa AbstractDict && v isa AbstractDict
            bad = setdiff(keys(v), keys(base[k]))
            isempty(bad) || error("Unknown key(s) $(collect(bad)) in $group → $k. Valid keys: $(sort(collect(keys(base[k])))).")
            base[k] = merge(Dict{String,Any}(base[k]), v)
        else
            base[k] = v
        end
    end
    return base
end

# A copy of `settings` with every group resolved against its overrides (group name => Dict).
# Everything downstream - the Benders loops, CutPool, the samplers - reads its group from the
# returned dict, so a run's overrides reach them without touching the caller's settings.
function resolve_all_settings(settings::AbstractDict, overrides::AbstractDict = Dict())
    unknown = setdiff(keys(overrides), SETTINGS_GROUPS)
    isempty(unknown) || error("Unknown settings group(s) $(collect(unknown)). Valid groups: $(collect(SETTINGS_GROUPS)).")
    resolved = copy(settings)
    for group in SETTINGS_GROUPS
        resolved[group] = resolve_settings(settings, group, get(overrides, group, Dict()))
    end
    return resolved
end
