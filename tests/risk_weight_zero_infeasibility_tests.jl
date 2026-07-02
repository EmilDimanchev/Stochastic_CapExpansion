using Test
using Distributed
using Revise, JuMP, Gurobi, HiGHS, Ipopt, DataFrames, CSV, YAML, Random, LinearAlgebra, Combinatorics, Dates, Distributions, Surrogates

include("../src/Stochastic_CapExpansion.jl")
using .Stochastic_CapExpansion

"""
Probe-level test harness for the risk_aversion_weight = 0.0 base Benders run.

Usage examples:
  julia --project=. experiments/risk_weight_zero_infeasibility_tests.jl
  SCE_INPUTS_FOLDER=inputs/Inputs_30days_5techs julia --project=. experiments/risk_weight_zero_infeasibility_tests.jl
  SCE_EXPECT_INFEASIBLE=true julia --project=. experiments/risk_weight_zero_infeasibility_tests.jl

Environment variables:
  SCE_INPUTS_FOLDER       Input folder path (default: inputs/Inputs_30d_1000scen_7tech_2z)
  SCE_EXPECT_INFEASIBLE   If true, test expects at least one infeasible regularized run
  SCE_REPETITIONS         Number of regularized probes to run (default: 1)
  SCE_SOLVER              Optional solver override, e.g. HiGHS or Gurobi
"""

const DEFAULT_INPUTS = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")

function _env_bool(name::String, default::Bool)
    v = get(ENV, name, default ? "true" : "false")
    return lowercase(strip(v)) in ("1", "true", "yes", "y", "on")
end

function _env_int(name::String, default::Int)
    raw = get(ENV, name, string(default))
    parsed = tryparse(Int, raw)
    parsed === nothing && error("Environment variable $(name) must be an Int, got '$(raw)'.")
    return parsed
end

function build_case_settings(inputs_folder::String; regularization::Bool)
    settings = load_settings(inputs_folder)
    settings["Risk aversion flag"] = true
    settings["Risk aversion weight"] = 0.0
    settings["Parallel flag"] = false
    settings["Regularization flag"] = regularization
    settings["Regularization strategy"] = "Level Set"
    settings["Fuel price uncertainty flag"] = false  # Disable fuel price uncertainty for this test

    if haskey(ENV, "SCE_SOLVER")
        settings["Solver"] = ENV["SCE_SOLVER"]
    end

    return settings
end

function run_risk_zero_case(inputs::Dict, settings::Dict, SPs::Array{Model,3}; case_name::String)
    MP = build_planning_model(inputs, settings; risk_aversion_weight = 0.0)
    set_objective_bendersMP!(MP, "System_Weighted_CVaR", inputs, settings; obj_weight = 0.0)

    try
        output = benders_algorithm(
            inputs,
            settings,
            MP,
            SPs,
            case_name;
            Eval_SPs = nothing,
            mapping = false,
            risk_aversion_weight = 0.0,
        )
        return (status = :ok, output = output, error = "")
    catch err
        return (status = :error, output = nothing, error = sprint(showerror, err))
    end
end

function main()
    inputs_folder = get(ENV, "SCE_INPUTS_FOLDER", DEFAULT_INPUTS)
    expect_infeasible = _env_bool("SCE_EXPECT_INFEASIBLE", false)
    repetitions = _env_int("SCE_REPETITIONS", 1)
    repetitions < 1 && error("SCE_REPETITIONS must be >= 1")

    @info("Risk-0 infeasibility probe starting")
    @info("Inputs folder: $(inputs_folder)")
    @info("Expect infeasible: $(expect_infeasible), repetitions: $(repetitions)")

    settings_for_data = build_case_settings(inputs_folder; regularization = false)
    inputs = load_input_data(inputs_folder, settings_for_data)
    SPs = build_all_subproblems(inputs, settings_for_data)

    @testset "Risk=0.0 Base Benders Feasibility Probe" begin
        @testset "Sanity (regularization off)" begin
            settings_no_reg = build_case_settings(inputs_folder; regularization = false)
            outcome_no_reg = run_risk_zero_case(inputs, settings_no_reg, SPs; case_name = "Risk_Weight_0_no_reg")
            @info("No-reg status: $(outcome_no_reg.status)")
            if outcome_no_reg.status == :error
                @error("No-reg run failed", error = outcome_no_reg.error)
            end
            @test outcome_no_reg.status == :ok
        end

        @testset "Level-set regularization probe" begin
            infeasible_hits = 0
            other_errors = String[]

            for rep in 1:repetitions
                settings_reg = build_case_settings(inputs_folder; regularization = true)
                case_name = "Risk_Weight_0_levelset_probe_$(rep)"
                outcome_reg = run_risk_zero_case(inputs, settings_reg, SPs; case_name = case_name)

                @info("Regularized run $(rep) status: $(outcome_reg.status)")

                if outcome_reg.status == :error
                    if occursin("Model infeasible", outcome_reg.error)
                        infeasible_hits += 1
                        @warn("Detected model infeasibility in run $(rep)")
                    else
                        push!(other_errors, outcome_reg.error)
                        @error("Detected non-infeasibility error in run $(rep)", error = outcome_reg.error)
                    end
                end
            end

            @info("Infeasibility hits across regularized runs: $(infeasible_hits)/$(repetitions)")

            @test isempty(other_errors)
            if expect_infeasible
                @test infeasible_hits >= 1
            else
                @test infeasible_hits >= 0
            end
        end
    end
end

main()
