#======================

GridFill sampler: entry points (src/model/sampling)

  test_gridfill_units()                 - solver-free unit tests (tests/gridfill_tests.jl)
  gridfill_test_laptop("gf_1")          - end-to-end pareto run with Sampling "Method" => "GridFill"
  gridfill_test_della("gf_1")           - the same on the Della inputs
  weather_gridfill_test_della("gf_1")   - Della run with weather instead of gas price uncertainty
  test_cut_export_units()               - solver-free unit tests of the cut export (tests/cut_export_tests.jl)
  test_memory_log_units()               - log_memory_status and worker BLAS threads (tests/memory_log_tests.jl)
  test_run_log_units()                  - "Log style" loggers, ITER|/SOLVE| lines and parse_run_log, plus a tiny della/laptop run (tests/run_log_tests.jl)
  test_nearest_distance_memory()        - memory/time of nearest_distance, full vs chunked, at Della sizes (needs ~3 GB)
  cut_export_test_laptop("cx_1")        - gridfill_test_laptop with the cut export on

=======================#

test_gridfill_units(test_index = "gridfill_units") = include(joinpath(@__DIR__, "..", "..", "tests", "gridfill_tests.jl"))
test_cut_export_units(test_index = "cut_export_units") = include(joinpath(@__DIR__, "..", "..", "tests", "cut_export_tests.jl"))
test_memory_log_units(test_index = "memory_log_units") = include(joinpath(@__DIR__, "..", "..", "tests", "memory_log_tests.jl"))
test_run_log_units(test_index = "run_log_units") = include(joinpath(@__DIR__, "..", "..", "tests", "run_log_tests.jl"))
test_nearest_distance_memory(test_index = "nearest_distance_memory") = include(joinpath(@__DIR__, "..", "..", "tests", "nearest_distance_memory.jl"))

# Laptop inputs, workers and subproblems. `fast` uses one gas and one weather scenario, keeping
# demand uncertainty on.
function _sampling_laptop_setup(; fast::Bool = false)
    inputs_folder = joinpath("inputs", "Inputs_2z_7tech")
    settings = load_settings(inputs_folder)
    if fast
        settings["Scenarios"]["Gas"]["Count"] = 1
        settings["Scenarios"]["Weather"]["Count"] = 1
    end
    inputs = load_input_data(inputs_folder, settings)
    configure_parallel_workers!(settings)
    SPs = build_all_subproblems(inputs, settings)
    return inputs, settings, SPs
end

# End-to-end: full risk-pareto + MGA + wing caps, then GridFill on the finals. Laptop defaults
# are small (coarse 0.1 cells, batches of 8) so the whole loop, feedback included, runs in a
# few minutes; `fast` uses one gas and one weather scenario (demand stays stochastic).
function gridfill_test_laptop(test_index; n_samples::Int = 24, fast::Bool = true, sampling_overrides::AbstractDict = Dict(), cut_export_overrides::AbstractDict = Dict(), run_overrides::AbstractDict = Dict())
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    mkpath(summary_folder)
    inputs, settings, SPs = _sampling_laptop_setup(; fast = fast)
    run_stochastic_exploration_risk_pareto(SPs, inputs, settings, joinpath(results_folder, "Pareto5"), summary_folder;
        budget_settings = Dict("Budget type" => "Transformed", "Tighten budget" => true),
        cap_settings = Dict("Cap wings" => true),
        mapping_settings = Dict("Mapping" => false),
        run_settings = merge(Dict{String,Any}("Summary name" => "pareto"), run_overrides),
        sampling_settings = merge(Dict{String,Any}("Method" => "GridFill", "Samples" => n_samples, "Batch size" => 8,
                                                   "GridFill" => Dict("Cell size" => 0.1)), sampling_overrides),
        cut_export_settings = cut_export_overrides)
end

# The same end-to-end run with the cut export on. Writes outputs/Test_<index>/Summary/cut_archive_pareto.h5.
cut_export_test_laptop(test_index; cut_export_overrides::AbstractDict = Dict(), kwargs...) =
    gridfill_test_laptop(test_index; cut_export_overrides = merge(Dict{String,Any}("Export cuts" => true), cut_export_overrides), kwargs...)

# `inputs_folder` picks the case and `settings_overrides` replaces top-level settings keys (e.g.
# "Subproblems per scenario") after the Della profile is loaded, before the inputs are read.
function gridfill_test_della(test_index; n_samples::Int = 2000, inputs_folder = joinpath("inputs", "Inputs_2z_7tech"), settings_overrides::AbstractDict = Dict(), sampling_overrides::AbstractDict = Dict(), cut_export_overrides::AbstractDict = Dict(), run_overrides::AbstractDict = Dict())
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    mkpath(summary_folder)
    settings = load_settings(inputs_folder; profile = "della")
    for (key, value) in settings_overrides
        haskey(settings, key) || error("Unknown settings key in settings_overrides: $(repr(key))")
        settings[key] = value
    end
    inputs = load_input_data(inputs_folder, settings)
    configure_parallel_workers!(settings)
    SPs = build_all_subproblems(inputs, settings)
    results_folder = joinpath(results_folder, "GridFill")
    run_stochastic_exploration_risk_pareto(SPs, inputs, settings, results_folder, summary_folder;
        budget_settings = Dict("Budget type" => "Transformed", "Tighten budget" => true),
        cap_settings = Dict("Cap wings" => true, "Cap levels" => [1, 3]),
        mapping_settings = Dict("Mapping" => false),
        run_settings = merge(Dict{String,Any}("Summary name" => "GridFill"), run_overrides),
        sampling_settings = merge(Dict{String,Any}("Method" => "GridFill", "Samples" => n_samples), sampling_overrides),
        cut_export_settings = cut_export_overrides)
end

# Scenario counts are read by load_input_data, so they must be set before it runs.
function weather_gridfill_test_della(test_index; n_samples::Int = 2000, sampling_overrides::AbstractDict = Dict())
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    inputs_folder = joinpath("inputs", "Inputs_2z_7tech")
    mkpath(summary_folder)
    settings = load_settings(inputs_folder; profile = "della")
    settings["Scenarios"]["Gas"]["Count"] = 1
    settings["Scenarios"]["Weather"]["Count"] = 10
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
