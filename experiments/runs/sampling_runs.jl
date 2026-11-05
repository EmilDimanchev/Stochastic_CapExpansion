#======================

GridFill sampler: entry points (src/model/sampling)

  test_gridfill_units()                 - solver-free unit tests (tests/gridfill_tests.jl)
  gridfill_test_laptop("gf_1")          - end-to-end pareto run with Sampling "Method" => "GridFill"
  gridfill_test_della("gf_1")           - the same on the Della inputs
  weather_gridfill_test_della("gf_1")   - Della run with weather instead of gas price uncertainty

=======================#

test_gridfill_units(test_index = "gridfill_units") = include(joinpath(@__DIR__, "..", "..", "tests", "gridfill_tests.jl"))

# Laptop inputs, workers and subproblems. `fast` turns off the gas price uncertainty set
# (10 instead of 100 scenarios), keeping demand uncertainty on.
function _sampling_laptop_setup(; fast::Bool = false)
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")
    settings = load_settings(inputs_folder)
    fast && (settings["Gas price uncertainty"] = false)
    inputs = load_input_data(inputs_folder, settings)
    configure_parallel_workers!(settings)
    SPs = build_all_subproblems(inputs, settings)
    return inputs, settings, SPs
end

# End-to-end: full risk-pareto + MGA + wing caps, then GridFill on the finals. Laptop defaults
# are small (coarse 0.1 cells, batches of 8) so the whole loop, feedback included, runs in a
# few minutes; `fast` turns off the gas price uncertainty set (demand stays stochastic).
function gridfill_test_laptop(test_index; n_samples::Int = 24, fast::Bool = true, sampling_overrides::AbstractDict = Dict())
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    mkpath(summary_folder)
    inputs, settings, SPs = _sampling_laptop_setup(; fast = fast)
    run_stochastic_exploration_risk_pareto(SPs, inputs, settings, joinpath(results_folder, "Pareto5"), summary_folder;
        budget_settings = Dict("Budget type" => "Transformed", "Tighten budget" => true),
        cap_settings = Dict("Cap wings" => true),
        mapping_settings = Dict("Mapping" => false),
        run_settings = Dict("Summary name" => "pareto"),
        sampling_settings = merge(Dict{String,Any}("Method" => "GridFill", "Samples" => n_samples, "Batch size" => 8,
                                                   "GridFill" => Dict("Cell size" => 0.1)), sampling_overrides))
end

function gridfill_test_della(test_index; n_samples::Int = 2000, sampling_overrides::AbstractDict = Dict())
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della")
    mkpath(summary_folder)
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)
    configure_parallel_workers!(settings)
    SPs = build_all_subproblems(inputs, settings)
    results_folder = joinpath(results_folder, "GridFill")
    run_stochastic_exploration_risk_pareto(SPs, inputs, settings, joinpath(results_folder, "Pareto5"), summary_folder;
        budget_settings = Dict("Budget type" => "Transformed", "Tighten budget" => true),
        cap_settings = Dict("Cap wings" => true, "Cap levels" => [1, 3]),
        mapping_settings = Dict("Mapping" => false),
        run_settings = Dict("Summary name" => "GridFill"),
        sampling_settings = merge(Dict{String,Any}("Method" => "GridFill", "Samples" => n_samples, "Write gradients" => true), sampling_overrides))
end

# The uncertainty flags are read by load_input_data, so they must be set before it runs.
function weather_gridfill_test_della(test_index; n_samples::Int = 2000, sampling_overrides::AbstractDict = Dict())
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della")
    mkpath(summary_folder)
    settings = load_settings(inputs_folder)
    settings["Gas price uncertainty"] = false
    settings["Weather uncertainty"] = true
    inputs = load_input_data(inputs_folder, settings)
    configure_parallel_workers!(settings)
    SPs = build_all_subproblems(inputs, settings)
    results_folder = joinpath(results_folder, "GridFill")
    run_stochastic_exploration_risk_pareto(SPs, inputs, settings, joinpath(results_folder, "Pareto5"), summary_folder;
        budget_settings = Dict("Budget type" => "Transformed", "Tighten budget" => true),
        cap_settings = Dict("Cap wings" => true, "Cap levels" => [1, 3]),
        mapping_settings = Dict("Mapping" => false),
        run_settings = Dict("Summary name" => "GridFill"),
        sampling_settings = merge(Dict{String,Any}("Method" => "GridFill", "Samples" => n_samples, "Write gradients" => true), sampling_overrides))
end
