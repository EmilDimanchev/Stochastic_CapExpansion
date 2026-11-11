using Test
using Distributed

"""
Tests that benders_algorithm's LB and UB describe the objective the master actually minimizes
(objective_risk_weight in src/model/benders/master_planning.jl, bound_risk_weight in algorithm.jl).

The bug: with "Risk aversion flag" false (the laptop settings.yaml default), set_objective_bendersMP!
still adds the CVaR terms and minimizes inv + w·E + (1 - w)·CVaR, but the bounds were built from the
flag, so they compared inv + E at the MP solution with the best inv + E seen so far. Those aren't a
lower and upper bound on the same thing, and on the 2/2/2 laptop case at w = 0.5 the run stopped with
"Negative gap detected" at iteration 15 (LB 124574.6, UB 123345.4).

Usage, from a REPL that already ran include("experiments/MGA_Tests.jl"):
  include("tests/bound_objective_tests.jl")
"""

isdefined(Main, :run_stochastic_exploration_risk_pareto) || include(joinpath(@__DIR__, "..", "experiments", "MGA_Tests.jl"))

function _bound_test_case(; gas = 2, weather = 2, demand = 2)
    folder = joinpath(@__DIR__, "..", "inputs", "Inputs_2z_7tech")
    settings = load_settings(folder)
    settings["Parallel flag"] = false
    settings["Risk aversion flag"] = false # the settings.yaml default; the case that failed
    settings["Scenarios"]["Gas"]["Count"] = gas
    settings["Scenarios"]["Weather"]["Count"] = weather
    settings["Scenarios"]["Demand"]["Count"] = demand
    inputs = load_input_data(folder, settings)
    return settings, inputs
end

@testset "Benders bounds follow the master's objective" begin
    settings, inputs = _bound_test_case()

    @testset "The master records its objective's risk weight" begin
        MP = build_planning_model(inputs, settings; risk_aversion_weight = 0.5)
        @test objective_risk_weight(MP) == 1.0 # flag off: expected cost only
        set_objective_bendersMP!(MP, "System_Weighted_CVaR", inputs, settings; obj_weight = 0.3)
        @test objective_risk_weight(MP) == 0.3
        @test settings["Risk aversion flag"] == false # not mutated behind the caller's back
        set_objective_bendersMP!(MP, "System_Expected", inputs, settings)
        @test objective_risk_weight(MP) == 1.0
        set_objective_bendersMP!(MP, "Capacity", inputs, settings; set_coeffs = ones(length(MP[:x]) + length(MP[:x_line])))
        @test objective_risk_weight(MP) === nothing
        @test_throws ErrorException bound_risk_weight(MP, nothing) # not a cost objective

        set_objective_bendersMP!(MP, "System_Weighted_CVaR", inputs, settings; obj_weight = 0.3)
        @test bound_risk_weight(MP, nothing) == 0.3
        @test bound_risk_weight(MP, 0.3) == 0.3
        @test_throws ErrorException bound_risk_weight(MP, 0.5) # caller and master disagree
    end

    @testset "CVaR outputs are reported whenever the master has CVaR terms" begin
        MP = build_planning_model(inputs, settings)
        set_objective_bendersMP!(MP, "System_Weighted_CVaR", inputs, settings; obj_weight = 0.5)
        optimize!(MP)
        out = write_outputs(MP, settings; risk_aversion_weight = 0.5)
        @test out["CVaR term"] == value(MP[:cvar_term])*settings["Scaling factor cost"]
        @test out["Risk adjusted system cost"] ≈ out["Planning objective"] rtol = 1e-9
    end

    @testset "Weighted anchor converges with the flag off" begin
        # Failed with "Negative gap detected" at iteration 15 before the fix
        settings["Convergence tolerance"] = 1e-3
        SPs = build_all_subproblems(inputs, settings)
        MP = build_planning_model(inputs, settings)
        set_objective_bendersMP!(MP, "System_Weighted_CVaR", inputs, settings; obj_weight = 0.5)
        out = benders_algorithm(inputs, settings, MP, SPs, "bound_test"; cut_pool = CutPool(settings))
        @test all(out["Gaps"] .>= -1e-6)
        @test out["Gaps"][end] <= 1e-3
        @test out["LB_hist"][end] <= out["UB_hist"][end]
    end
end
