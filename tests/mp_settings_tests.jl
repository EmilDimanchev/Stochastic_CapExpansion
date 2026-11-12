using Test
using Distributed

"""
Tests for the master/subproblem solver attributes (master_gurobi_attributes, level_set_gurobi_attributes,
subproblem_gurobi_attributes), the cut pool's "Seed scope", and PHASE| markers.

There is no Gurobi license on the laptop, so the Gurobi attribute lists are checked as data. The
restore after a level-set solve is checked on HiGHS in tests/level_set_tests.jl.

Usage, from a REPL that already ran include("experiments/MGA_Tests.jl"):
  include("tests/mp_settings_tests.jl")
"""

isdefined(Main, :run_stochastic_exploration_risk_pareto) || include(joinpath(@__DIR__, "..", "experiments", "MGA_Tests.jl"))

@testset "Solver attributes and seed scope" begin
    settings = load_settings(joinpath(@__DIR__, "..", "inputs", "Inputs_2z_7tech"))

    @testset "Gurobi attribute lists" begin
        master = Dict(master_gurobi_attributes(settings))
        @test master["Method"] == 3       # concurrent: uses more than one core
        @test master["LPWarmStart"] == 1  # simplex warm-starts; barrier isn't forced through the basis
        @test master["Crossover"] == -1   # a barrier win still leaves a basis
        @test master["Threads"] == master_settings(settings)["Threads"] == 0

        level = Dict(level_set_gurobi_attributes(settings))
        @test level["Method"] == 2 && level["LPWarmStart"] == 0
        @test level["Crossover"] == 0     # interior point on purpose
        # Every attribute the level-set solve changes is one the restore sets back
        @test issubset(keys(level), keys(master))

        sp = Dict(subproblem_gurobi_attributes(settings))
        @test sp["Threads"] == settings["Subproblem threads"] == 1
        @test sp["Method"] == -1 && sp["Crossover"] == 1

        s = apply_settings_overrides(settings, Dict("Master problem settings" => Dict("Threads" => 64), "Subproblem threads" => 2))
        @test Dict(master_gurobi_attributes(s))["Threads"] == 64
        @test Dict(level_set_gurobi_attributes(s))["Threads"] == 64
        @test Dict(subproblem_gurobi_attributes(s))["Threads"] == 2
    end

    @testset "Seed scope" begin
        MP, pool, inputs, s, n_sub = _cut_test_setup(4)
        all_names = collect(keys(pool.live))
        @test pool.cfg["Seed scope"] == "group"

        # A solve: MP solve, then prune! records what binds
        run_planning_model(MP, s, s["Risk aversion weight"])
        prune!(pool, MP, 1, 2)
        binding = Set(n for (n, con) in pool.live if dual(con) >= pool.cfg["Binding tolerance"])
        @test !isempty(binding)
        @test pool.bound_last == binding
        @test binding ⊆ pool.bound_cur

        # A later prune! replaces bound_last rather than adding to it
        pool.bound_cur = union(pool.bound_cur, Set(all_names[1:5]))
        prune!(pool, MP, 2, 2)
        @test pool.bound_last == binding # same solution, same binding set

        # "group": the group's binding set (unchanged behavior)
        pool.group = 0
        seed!(pool, MP, 0)
        @test Set(keys(pool.live)) == pool.bound_cur

        # "last solve": only what bound at the most recent prune!
        pool.cfg["Seed scope"] = "last solve"
        seed!(pool, MP, 0)
        @test Set(keys(pool.live)) == binding
        @test pool.bound_cur ⊇ binding # group bookkeeping still kept

        # "none": nothing live; every cut is archived and none is lost
        pool.cfg["Seed scope"] = "none"
        seed!(pool, MP, 1)
        @test isempty(pool.live) && count_live_cuts(MP) == 0
        @test n_archived(pool) == length(all_names)

        # Separation alone brings back what the solution violates
        output, sep = separate!(pool, MP, s, s["Risk aversion weight"])
        @test sep.restored > 0

        bad = deepcopy(s)
        bad["Cut settings"]["Seed scope"] = "everything"
        @test_throws ErrorException CutPool(bad)
    end

    @testset "PHASE| markers" begin
        buf = IOBuffer()
        with_logger(CompactLogger(buf)) do
            log_phase("MP solve (10 live cuts)")
        end
        @test occursin("PHASE| MP solve (10 live cuts)", String(take!(buf)))
        with_logger(run_logger("laptop", buf)) do
            log_phase("hidden in the laptop style")
        end
        @test isempty(take!(buf))
    end
end
