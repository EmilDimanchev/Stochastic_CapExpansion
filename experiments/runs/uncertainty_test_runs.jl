function weather_gridfill_test_della(test_index; n_samples::Int = 2000, sampling_overrides::AbstractDict = Dict())
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della")
    mkpath(summary_folder)
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)
    configure_parallel_workers!(settings)
    SPs = build_all_subproblems(inputs, settings)
    results_folder = joinpath(results_folder, "GridFill")
    settings["Gas price uncertainty"] = false
    settings["Weather uncertainty"] =  true
    run_stochastic_exploration_risk_pareto(SPs, inputs, settings, joinpath(results_folder, "Pareto5"), summary_folder;
        budget_settings = Dict("Budget type" => "Transformed", "Tighten budget" => true),
        cap_settings = Dict("Cap wings" => true, "Cap levels" => [1, 3]),
        mapping_settings = Dict("Mapping" => false),
        run_settings = Dict("Summary name" => "GridFill"),
        sampling_settings = merge(Dict{String,Any}("Method" => "GridFill", "Samples" => n_samples, "Write gradients" => true), sampling_overrides))
end