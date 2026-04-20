function load_settings(inputpath::String)
    settings_path = joinpath(inputpath, "settings.yaml")

    settings = YAML.load_file(settings_path)
    return settings
end