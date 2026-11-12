using Test
using Distributed

"""
Tests for the level-set regularization row (add_level_set_row!, level_set_regularization), the master
basis bookkeeping kept across a level-set solve (row_basis_statuses, save_/restore_master_basis!) and
run_planning_model's retry, all in src/model/benders/master_planning.jl.
HiGHS only: the Gurobi basis calls themselves are checked on Della (see the plan's verification notes).

Usage, from a REPL that already ran include("experiments/MGA_Tests.jl"):
  include("tests/level_set_tests.jl")
"""

isdefined(Main, :run_stochastic_exploration_risk_pareto) || include(joinpath(@__DIR__, "..", "experiments", "MGA_Tests.jl"))

# Coefficients of an affine function, keyed by variable, without zeros
_coefs(f) = Dict(v => c for (v, c) in f.terms if c != 0)
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

    @testset "Built from the current objective, deleted after the solve" begin
        for w in (0.5, 1.0) # one MP across anchors: each pass must cap the current eObj
            set_objective_bendersMP!(MP, "System_Weighted_CVaR", inputs, settings; obj_weight = w)
            row = add_level_set_row!(MP, 10.0)
            @test _coefs(constraint_object(row).func) == _coefs(MP[:eObj])
            @test normalized_rhs(row) == 10.0
            @test_throws ErrorException add_level_set_row!(MP, 20.0) # never two rows
            delete_level_set_row!(MP)
            @test !is_valid(MP, row)
            @test !haskey(MP.obj_dict, :cLevel_set)
            @test _n_level_rows(MP) == 0
        end
    end

    @testset "Regularized solve leaves the master as it found it" begin
        set_objective_bendersMP!(MP, "System_Weighted_CVaR", inputs, settings; obj_weight = 0.5)
        optimize!(MP)
        LB = objective_value(MP)
        n_rows = length(master_rows(MP))
        out = level_set_regularization(MP, 2LB + 1.0, LB, 0.5, settings)
        @test out !== nothing
        @test !haskey(MP.obj_dict, :cLevel_set) # deleted, not relaxed to an Inf RHS
        @test length(master_rows(MP)) == n_rows
        @test objective_function(MP) == MP[:eObj]
        if settings["Solver"] == "HiGHS"
            @test get_optimizer_attribute(MP, "primal_feasibility_tolerance") == 1e-3
            @test get_optimizer_attribute(MP, "optimality_tolerance") == 1e-3
            @test get_optimizer_attribute(MP, "user_bound_scale") == -6
        end
        optimize!(MP)
        @test termination_status(MP) == MOI.OPTIMAL
        @test objective_value(MP) ≈ LB rtol = 1e-6
        # A second pass works too (the first one's row is gone)
        @test level_set_regularization(MP, 2LB + 1.0, LB, 0.5, settings) !== nothing
    end

    @testset "Saved basis maps onto the rows present at restore" begin
        x = first(all_variables(MP))
        kept = @constraint(MP, x >= -1e6)
        pruned = @constraint(MP, x <= 1e9)
        saved = Dict(kept => -1, pruned => -1)
        delete(MP, pruned)
        added = @constraint(MP, 2x >= -1e6)
        statuses = Dict(row_basis_statuses(saved, master_rows(MP)))
        @test statuses[kept] == -1    # saved status kept
        @test statuses[added] == 0    # new row: basic slack
        @test !haskey(statuses, pruned) # deleted row dropped
        @test length(statuses) == length(master_rows(MP))
        delete(MP, [kept, added])

        # HiGHS rebuilds the solver on its level-set swap, so nothing is saved or restored
        if settings["Solver"] != "Gurobi"
            optimize!(MP)
            save_master_basis!(MP, settings)
            @test !haskey(MP.ext, :saved_basis)
            @test restore_master_basis!(MP) === nothing
        end
    end

    @testset "run_planning_model returns an output or fails clearly" begin
        out = run_planning_model(MP, settings, 0.5)
        @test haskey(out, "Expected alpha")
        if settings["Solver"] == "HiGHS"
            # A row that cuts off the last optimum, so HiGHS has to iterate, and no time to do it: the
            # first solve and the retry both stop early
            x = MP[:x][1]
            dummy = @constraint(MP, x >= value(x) + 1.0)
            set_optimizer_attribute(MP, "time_limit", 0.0)
            @test_throws ErrorException run_planning_model(MP, settings, 0.5)
            set_optimizer_attribute(MP, "time_limit", Inf)
            delete(MP, dummy)
            @test haskey(run_planning_model(MP, settings, 0.5), "Expected alpha")
        end
    end
end
