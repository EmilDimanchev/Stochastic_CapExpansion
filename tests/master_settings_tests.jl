using Test
using Distributed

"""
Tests for the "Master problem settings" block (master_settings, check_master_settings,
nest_master_settings, apply_settings_overrides in src/settings/load_settings.jl).

Usage, from a REPL that already ran include("experiments/MGA_Tests.jl"):
  include("tests/master_settings_tests.jl")
"""

isdefined(Main, :run_stochastic_exploration_risk_pareto) || include(joinpath(@__DIR__, "..", "experiments", "MGA_Tests.jl"))

_errors_with(f, text) = try
    f()
    false
catch e
    occursin(text, sprint(showerror, e))
end

@testset "Master problem settings" begin
    inputs_root = joinpath(@__DIR__, "..", "inputs")

    @testset "Every input folder and profile loads with a complete block" begin
        for folder in ("Inputs_2z_7tech", "Inputs_16z_GenX"), profile in (nothing, "della")
            settings = load_settings(joinpath(inputs_root, folder); profile = profile)
            @test Set(keys(master_settings(settings))) == Set(MASTER_SETTINGS_KEYS)
            @test !any(k -> haskey(settings, k), keys(MOVED_TO_MASTER_SETTINGS))
        end
        @test master_settings(load_settings(joinpath(inputs_root, "Inputs_16z_GenX"); profile = "della"))["Threads"] == 32
    end

    base = load_settings(joinpath(inputs_root, "Inputs_2z_7tech"))

    @testset "Old top-level keys and bad blocks fail clearly" begin
        old = merge(base, Dict("MP threads" => 8))
        @test _errors_with(() -> check_master_settings(old), "\"Threads\"")
        @test _errors_with(() -> check_master_settings(delete!(copy(base), MASTER_SETTINGS)), "no \"Master problem settings\" block")
        typo = merge_settings_profile(base, Dict(MASTER_SETTINGS => Dict("Thread" => 8)))
        @test _errors_with(() -> check_master_settings(typo), "unknown [\"Thread\"]")
        bad = merge_settings_profile(base, Dict(MASTER_SETTINGS => Dict("Regularization strategy" => "level set")))
        @test _errors_with(() -> check_master_settings(bad), "Regularization strategy")
    end

    @testset "nest_master_settings moves old keys into the block" begin
        nested = nest_master_settings(Dict{String,Any}("MP threads" => 4, "Regularization flag" => false, "Solver" => "HiGHS"))
        @test nested == Dict("Solver" => "HiGHS", MASTER_SETTINGS => Dict("Threads" => 4, "Regularization flag" => false))
        @test nest_master_settings(Dict{String,Any}("Solver" => "HiGHS")) == Dict("Solver" => "HiGHS") # nothing to move
    end

    @testset "apply_settings_overrides merges blocks and rejects unknown keys" begin
        s = apply_settings_overrides(base, Dict(MASTER_SETTINGS => Dict("Threads" => 32), "Workers" => 3))
        @test master_settings(s)["Threads"] == 32 && s["Workers"] == 3
        @test master_settings(s)["Regularization flag"] == master_settings(base)["Regularization flag"] # siblings kept
        @test master_settings(base)["Threads"] == 0 # caller's dict untouched
        @test _errors_with(() -> apply_settings_overrides(base, Dict("MP threads" => 32)), "\"Threads\" in the \"Master problem settings\" block")
        @test _errors_with(() -> apply_settings_overrides(base, Dict(MASTER_SETTINGS => Dict("Thread" => 32))), "Unknown key \"Thread\"")
    end
end
