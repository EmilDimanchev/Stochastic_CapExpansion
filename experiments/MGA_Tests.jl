using Distributed

@everywhere include("../src/Stochastic_CapExpansion.jl")

@everywhere using .Stochastic_CapExpansion
@everywhere using Revise, JuMP, Gurobi, HiGHS, Ipopt, DataFrames, CSV, YAML, Random, LinearAlgebra, Combinatorics, Dates, Distributions, Surrogates, DelaunayTriangulation
using Clarabel
using Serialization   # caches the run-13 gradient bank (run13_boundary_checks)

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
include("mga/interior_sampling.jl")       # MGCA problem, Delaunay interior sampling
include("mga/interpolate_evaluation.jl")  # evaluate_interpolates, evaluate_subproblems

# Experiment runners
include("runs/mapping_runs.jl")           # mapping_test_*
include("runs/risk_pareto_runs.jl")       # risk_pareto_*, wing_caps_test_*, warmstart_verification_laptop, simple_comp
include("runs/interpolate_runs.jl")       # run_interpolate_evaluation_*
include("runs/boundary_sampling_runs.jl") # boundary_sampling_test_*, gradient_check, run13_boundary_checks

# Diagnostic checks
include("checks/cut_checks.jl")           # test_cut_management, test_cut_trim_rules
include("checks/model_checks.jl")         # test_capacity_bound_tightening, test_wing_caps
include("checks/sampler_checks.jl")       # test_hit_and_run_sampler, test_sample_interior_delaunay
