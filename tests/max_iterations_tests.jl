using Test
using Distributed

"""
Tests for the "Max iterations" setting (max_benders_iterations in src/model/benders/algorithm.jl).

Usage, from a REPL that already ran include("experiments/MGA_Tests.jl"):
  include("tests/max_iterations_tests.jl")
"""

isdefined(Main, :run_stochastic_exploration_risk_pareto) || include(joinpath(@__DIR__, "..", "experiments", "MGA_Tests.jl"))

@testset "Max iterations" begin
    @testset "Setting parsing" begin
        @test max_benders_iterations(Dict{String,Any}()) == 150 # settings files without the key
        @test max_benders_iterations(Dict{String,Any}("Max iterations" => 7)) == 7
        @test_throws ErrorException max_benders_iterations(Dict{String,Any}("Max iterations" => 0))
        @test_throws ErrorException max_benders_iterations(Dict{String,Any}("Max iterations" => 10.5))
        @test_throws ErrorException max_benders_iterations(Dict{String,Any}("Max iterations" => nothing))
    end

    @testset "Base settings files set it" begin
        for folder in ("Inputs_2z_7tech", "Inputs_16z_GenX")
            @test load_settings(joinpath(@__DIR__, "..", "inputs", folder))["Max iterations"] == 150
        end
    end

    @testset "Benders stops at the cap" begin
        folder = joinpath(@__DIR__, "..", "inputs", "Inputs_2z_7tech")
        settings = load_settings(folder)
        settings["Parallel flag"] = false
        settings["Scenarios"]["Gas"]["Count"] = 2
        settings["Scenarios"]["Weather"]["Count"] = 1
        settings["Scenarios"]["Demand"]["Count"] = 1
        settings["Risk aversion flag"] = true
        settings["Max iterations"] = 2
        inputs = load_input_data(folder, settings)
        SPs = build_all_subproblems(inputs, settings)
        MP = build_planning_model(inputs, settings; risk_aversion_weight = 0.5)
        set_objective_bendersMP!(MP, "System_Weighted_CVaR", inputs, settings; obj_weight = 0.5)

        err = try
            benders_algorithm(inputs, settings, MP, SPs, "max_iter_test"; risk_aversion_weight = 0.5)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test err !== nothing && occursin("within 2 iterations", sprint(showerror, err))
    end
end
