using Distributed

@everywhere include("../src/Stochastic_CapExpansion.jl")

@everywhere using .Stochastic_CapExpansion
@everywhere using Revise, JuMP, Gurobi, HiGHS, Ipopt, DataFrames, CSV, YAML, Random, LinearAlgebra, Combinatorics, Dates, Distributions
using Clarabel

#=========================

Entry point for the MGA / risk-pareto experiments. Include this file, then call any runner
or check with a test_index. The code lives in the files below:

  mga/     reusable pieces the runners build on
  runs/    laptop and Della entry points (call with a test_index)
  checks/  self-contained diagnostic tests (test_* functions)

=========================#

# Core helpers and the main pareto/MGA driver
include("mga/worker_utils.jl")            # memory logging, parallel worker setup
include("mga/budget_schedules.jl")        # budget floor/tightening schedule, wing caps
include("mga/run_settings.jl")            # default_*_settings + resolve_settings
include("mga/risk_pareto.jl")             # run_stochastic_exploration_risk_pareto, run_base_mga
include("mga/design_evaluation.jl")       # evaluate_subproblems, read_design_csv, evaluate_design_csv

# Experiment runners
include("runs/mapping_runs.jl")           # mapping_test_*
include("runs/risk_pareto_runs.jl")       # risk_pareto_*, wing_caps_test_*, warmstart_verification_laptop, simple_comp
include("runs/design_evaluation_runs.jl") # run_design_csv_evaluation_*
include("runs/sampling_runs.jl")          # test_gridfill_units, gridfill_test_*, weather_gridfill_test_della, test_cut_export_units, cut_export_test_laptop

# Diagnostic checks
include("checks/cut_checks.jl")           # test_cut_management, test_cut_trim_rules
include("checks/model_checks.jl")         # test_capacity_bound_tightening, test_wing_caps, test_settings_groups, test_case_file_flag, test_subproblems_batch, test_design_csv_evaluation, test_cost_by_zone, test_temporal_blocks, test_worker_only_subproblems
include("checks/input_checks.jl")         # test_scenario_selection, test_cost_construction, test_thermal_availability, test_weather_library, test_input_conversion, test_scenario_counts_end_to_end, test_genx_conversion
include("checks/sampling_checks.jl")      # test_gridfill_real_model
