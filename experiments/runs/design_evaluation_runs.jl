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
