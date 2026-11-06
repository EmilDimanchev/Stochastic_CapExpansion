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
    profile === nothing && return settings
    profile_path = joinpath(inputpath, "settings_$(profile).yaml")
    isfile(profile_path) || error("No settings profile \"$profile\": $profile_path doesn't exist.")
    overrides = YAML.load_file(profile_path)
    return overrides === nothing ? settings : merge_settings_profile(settings, overrides)
end
