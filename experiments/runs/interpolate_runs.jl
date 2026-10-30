#=========

Run interpolate evaluation

==========#

function run_interpolate_evaluation_laptop(test_index)
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")#joinpath("inputs", "Inputs_30repdays_ext_1000scen_7techs")
    results_folder = joinpath("outputs", "Test_"*string(test_index), "Interpolate_Evaluation")
    summary_folder = joinpath(results_folder, "Summary")
    interp_file = joinpath("experiments", "Della_experiments", "Mapping","mapping_budget_2p_MGCA.csv")

    if !isdir(results_folder)
        mkpath(results_folder)
    end
    if !isdir(summary_folder)
        mkpath(summary_folder)
    end
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)
    interpolate_df = CSV.read(interp_file, DataFrame, header=true)

    configure_parallel_workers!(settings)

    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)

    evaluate_interpolates(SPs, inputs, settings, interpolate_df, summary_folder)

end



function run_interpolate_evaluation_della(test_index, interp_file)
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della")#joinpath("inputs", "Inputs_30repdays_ext_1000scen_7techs")
    results_folder = joinpath("outputs", "Test_"*string(test_index), "Interpolate_Evaluation")
    summary_folder = joinpath(results_folder, "Summary")

    if !isdir(results_folder)
        mkpath(results_folder)
    end
    if !isdir(summary_folder)
        mkpath(summary_folder)
    end
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)
    interpolate_df = CSV.read(interp_file, DataFrame, header=true)

    configure_parallel_workers!(settings)

    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)

    evaluate_interpolates(SPs, inputs, settings, interpolate_df, summary_folder)

end



#=========

Evaluate planner design CSVs (see read_design_csv), e.g.
  run_design_csv_evaluation_laptop("design_csv", "experiments/designs/run13_RA_planner_designs.csv"; rows = 1:10)

==========#

function run_design_csv_evaluation(inputs_folder, test_index, design_file; batch_size = 25, rows = nothing)
    summary_folder = joinpath("outputs", "Test_"*string(test_index), "Design_Evaluation", "Summary")
    mkpath(summary_folder)
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)

    configure_parallel_workers!(settings)
    SPs = build_all_subproblems(inputs, settings)

    evaluate_design_csv(SPs, inputs, settings, design_file, summary_folder; batch_size = batch_size, rows = rows)
end

run_design_csv_evaluation_laptop(test_index, design_file; kwargs...) =
    run_design_csv_evaluation(joinpath("inputs", "Inputs_30d_1000scen_7tech_2z"), test_index, design_file; kwargs...)

run_design_csv_evaluation_della(test_index, design_file; kwargs...) =
    run_design_csv_evaluation(joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della"), test_index, design_file; kwargs...)
