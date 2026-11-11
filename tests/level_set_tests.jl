using Test
using Distributed

"""
Tests for the level-set regularization row (set_level_set_row! and level_set_regularization in
src/model/benders/master_planning.jl).

Usage, from a REPL that already ran include("experiments/MGA_Tests.jl"):
  include("tests/level_set_tests.jl")
"""

isdefined(Main, :run_stochastic_exploration_risk_pareto) || include(joinpath(@__DIR__, "..", "experiments", "MGA_Tests.jl"))

# Coefficients of an affine function, keyed by variable, without zeros
_coefs(f) = Dict(v => c for (v, c) in f.terms if c != 0)
_row_matches_eObj(MP) = _coefs(constraint_object(MP[:cLevel_set]).func) == _coefs(MP[:eObj])
_n_level_rows(MP) = count(c -> name(c) == "cLevel_set", all_constraints(MP, AffExpr, MOI.LessThan{Float64}))

@testset "Level-set row" begin
    folder = joinpath(@__DIR__, "..", "inputs", "Inputs_2z_7tech")
    settings = load_settings(folder)
    settings["Parallel flag"] = false
    settings["Risk aversion flag"] = true
    settings["Scenarios"]["Gas"]["Count"] = 2
    settings["Scenarios"]["Weather"]["Count"] = 1
    settings["Scenarios"]["Demand"]["Count"] = 1
    inputs = load_input_data(folder, settings)
    MP = build_planning_model(inputs, settings; risk_aversion_weight = 0.5)

    @testset "Follows the objective across anchors" begin
        set_objective_bendersMP!(MP, "System_Weighted_CVaR", inputs, settings; obj_weight = 0.5)
        row = set_level_set_row!(MP, 10.0)
        @test _row_matches_eObj(MP)
        @test set_level_set_row!(MP, 20.0) == row # same objective: RHS updated in place
        @test normalized_rhs(row) == 20.0

        set_objective_bendersMP!(MP, "System_Weighted_CVaR", inputs, settings; obj_weight = 1.0)
        @test !_row_matches_eObj(MP) # the old row no longer describes the objective...
        set_level_set_row!(MP, 30.0)
        @test _row_matches_eObj(MP) # ...so it is rebuilt
        @test _n_level_rows(MP) == 1
        @test !is_valid(MP, row)
        @test normalized_rhs(MP[:cLevel_set]) == 30.0
    end

    @testset "Regularized solve leaves the master as it found it" begin
        set_objective_bendersMP!(MP, "System_Weighted_CVaR", inputs, settings; obj_weight = 0.5)
        optimize!(MP)
        LB = objective_value(MP)
        out = level_set_regularization(MP, 2LB + 1.0, LB, 0.5, settings)
        @test out !== nothing
        @test _row_matches_eObj(MP)
        @test normalized_rhs(MP[:cLevel_set]) == Inf # relaxed, not a large finite RHS
        @test objective_function(MP) == MP[:eObj]
        if settings["Solver"] == "HiGHS"
            @test get_optimizer_attribute(MP, "primal_feasibility_tolerance") == 1e-3
            @test get_optimizer_attribute(MP, "optimality_tolerance") == 1e-3
            @test get_optimizer_attribute(MP, "user_bound_scale") == -6
        end
        optimize!(MP) # an ordinary solve with the relaxed row
        @test termination_status(MP) == MOI.OPTIMAL
        @test objective_value(MP) ≈ LB rtol = 1e-6
    end
end
