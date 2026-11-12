# Replaces the keys named in `overrides`, merging nested blocks so sibling keys are kept
function merge_settings_profile(base::AbstractDict, overrides::AbstractDict)
    merged = copy(base)
    for (key, value) in overrides
        merged[key] = (value isa AbstractDict && get(base, key, nothing) isa AbstractDict) ? merge_settings_profile(base[key], value) : value
    end
    return merged
end

# `profile` applies settings_<profile>.yaml, which holds only the keys that differ from settings.yaml
function load_settings(inputpath::String; profile::Union{Nothing,String} = nothing)
    settings = YAML.load_file(joinpath(inputpath, "settings.yaml"))
    if profile !== nothing
        profile_path = joinpath(inputpath, "settings_$(profile).yaml")
        isfile(profile_path) || error("No settings profile \"$profile\": $profile_path doesn't exist.")
        overrides = YAML.load_file(profile_path)
        overrides === nothing || (settings = merge_settings_profile(settings, overrides))
    end
    check_master_settings(settings)
    return settings
end

# Replaces the keys named in `overrides` like a profile does (nested blocks merged, siblings kept),
# but errors on any key, at any depth, that `settings` doesn't already have, so a typo or a key from
# before a block was introduced can't silently add an unused entry.
function apply_settings_overrides(settings::AbstractDict, overrides::AbstractDict; path::String = "settings")
    for (key, value) in overrides
        if !haskey(settings, key)
            hint = path == "settings" && haskey(MOVED_TO_MASTER_SETTINGS, key) ?
                " It is now \"$(MOVED_TO_MASTER_SETTINGS[key])\" in the \"$MASTER_SETTINGS\" block: pass Dict(\"$MASTER_SETTINGS\" => Dict(\"$(MOVED_TO_MASTER_SETTINGS[key])\" => ...))." : ""
            error("Unknown key $(repr(key)) in $path.$hint")
        end
        value isa AbstractDict && settings[key] isa AbstractDict && apply_settings_overrides(settings[key], value; path = "$path → $key")
    end
    return merge_settings_profile(settings, overrides)
end

#-------------------------------------------------------------------------------------------
# "Master problem settings": how the Benders master (investment) problem is solved and
# stabilized. Read through master_settings(settings); load_settings checks the block.
#-------------------------------------------------------------------------------------------

const MASTER_SETTINGS = "Master problem settings"
const MASTER_SETTINGS_KEYS = ("Threads", "Regularization flag", "Regularization strategy", "LB stall limit", "LB stall tolerance")
const REGULARIZATION_STRATEGIES = ("Level Set", "QTR")

# Top-level keys that moved into the block (old name => name in the block)
const MOVED_TO_MASTER_SETTINGS = Dict(
    "MP threads" => "Threads",
    "Regularization flag" => "Regularization flag",
    "Regularization strategy" => "Regularization strategy",
    "LB stall limit" => "LB stall limit",
    "LB stall tolerance" => "LB stall tolerance")

# The values the code fell back on before the block existed; only convert_input_folder uses these,
# to complete the block for an old folder that never set some of the keys.
const MASTER_SETTINGS_DEFAULTS = Dict{String,Any}("Threads" => 0, "Regularization flag" => true,
    "Regularization strategy" => "Level Set", "LB stall limit" => 10, "LB stall tolerance" => 1e-6)

master_settings(settings::AbstractDict) = settings[MASTER_SETTINGS]

function check_master_settings(settings::AbstractDict)
    moved = [k for k in keys(MOVED_TO_MASTER_SETTINGS) if haskey(settings, k)]
    isempty(moved) || error("Settings key(s) $(moved) moved into the \"$MASTER_SETTINGS\" block (as " *
        join(("\"$(MOVED_TO_MASTER_SETTINGS[k])\"" for k in moved), ", ") * "). Move them there, or run nest_master_settings on an old settings dict.")
    haskey(settings, MASTER_SETTINGS) || error("settings.yaml has no \"$MASTER_SETTINGS\" block (keys: $(collect(MASTER_SETTINGS_KEYS))).")
    block = settings[MASTER_SETTINGS]
    missing_keys = String[k for k in MASTER_SETTINGS_KEYS if !haskey(block, k)]
    unknown = String[k for k in keys(block) if !(k in MASTER_SETTINGS_KEYS)]
    isempty(missing_keys) && isempty(unknown) || error("\"$MASTER_SETTINGS\": missing $(missing_keys), unknown $(unknown). Valid keys: $(collect(MASTER_SETTINGS_KEYS)).")
    block["Regularization strategy"] in REGULARIZATION_STRATEGIES || error("\"$MASTER_SETTINGS\" → \"Regularization strategy\" must be one of $(REGULARIZATION_STRATEGIES), got $(repr(block["Regularization strategy"])).")
    return nothing
end

# Moves old top-level keys into the block, for settings written before it existed
# (convert_input_folder uses this on old folders). Keys the old dict doesn't have are left out.
function nest_master_settings(settings::AbstractDict)
    nested = Dict{String,Any}(k => v for (k, v) in settings if !haskey(MOVED_TO_MASTER_SETTINGS, k))
    block = Dict{String,Any}(MOVED_TO_MASTER_SETTINGS[k] => v for (k, v) in settings if haskey(MOVED_TO_MASTER_SETTINGS, k))
    isempty(block) || (nested[MASTER_SETTINGS] = merge(get(settings, MASTER_SETTINGS, Dict{String,Any}()), block))
    return nested
end
