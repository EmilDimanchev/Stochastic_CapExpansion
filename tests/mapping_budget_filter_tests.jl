using Test
using Distributed
using Logging
using DataFrames

"""
Standalone probe for the budget-violation filter added to sample_interior_delauney
(experiments/MGA_Tests.jl). Confirms that mapped solutions whose transformed cost
(transform_sys*(inv+ev) + transform_cvar*(inv+cvar)) exceeds budget + "Mapping Gap
Threshold" are dropped before Delaunay interior sampling, while solutions at or under
the cutoff are kept.

Usage:
  julia --project=. tests/mapping_budget_filter_tests.jl
"""

include(joinpath(@__DIR__, "..", "experiments", "MGA_Tests.jl"))

function main()
    budget = 1.0
    transform_sys = 1.0
    transform_cvar = 1.0
    settings = Dict("Mapping Gap Threshold" => 0.2, "Seed" => 42)
    cutoff = budget + settings["Mapping Gap Threshold"] # 1.2

    @testset "Mapping budget filter (sample_interior_delauney)" begin

        @testset "Gross budget violation is dropped and excluded from sampling" begin
            # Four solutions inside/at the unit square, all with transformed cost 1.0 <= cutoff,
            # plus one gross outlier (transformed cost 62.0) mirroring the reported bug
            # ("budgets of 31 or something when the normal threshold is 1").
            points = [0.0 0.0; 1.0 0.0; 0.0 1.0; 1.0 1.0; 100.0 100.0]
            cost_df = DataFrame(
                Investment_Cost = [0.0, 0.0, 0.0, 0.0, 0.0],
                EV_OpCost       = [0.5, 0.4, 0.6, 0.6, 31.0],
                CVaR_OpCost     = [0.5, 0.6, 0.4, 0.6, 31.0],
            )

            logs, samples = Test.collect_test_logs() do
                sample_interior_delauney(points, cost_df, 30, settings, budget, transform_sys, transform_cvar)
            end

            drop_logs = filter(l -> occursin("Dropping", l.message), logs)
            @test length(drop_logs) == 1
            @test occursin("Dropping 1 of 5", drop_logs[1].message)

            @test length(samples) == 30
            # The outlier row's capacities are [100, 100]; if it had leaked into the
            # triangulation, some barycentric sample would land far outside the unit
            # square. Every sample must stay a convex combination of the four kept
            # (unit-square) points.
            for s in samples
                @test all(0.0 - 1e-9 .<= s .<= 1.0 + 1e-9)
            end
        end

        @testset "Solution exactly at the cutoff is kept (<=, not <)" begin
            # Same four unit-square points; the fourth is exactly at the cutoff
            # (transformed cost == 1.2) and must not be dropped.
            points = [0.0 0.0; 1.0 0.0; 0.0 1.0; 1.0 1.0]
            cost_df = DataFrame(
                Investment_Cost = [0.0, 0.0, 0.0, 0.0],
                EV_OpCost       = [0.5, 0.4, 0.6, 0.6],
                CVaR_OpCost     = [0.5, 0.6, 0.4, 0.6],
            )
            @test cost_df[4, "Investment_Cost"] + cost_df[4, "EV_OpCost"] == 0.6
            @test transform_sys * (0.6) + transform_cvar * (0.6) == cutoff

            logs, samples = Test.collect_test_logs() do
                sample_interior_delauney(points, cost_df, 10, settings, budget, transform_sys, transform_cvar)
            end

            drop_logs = filter(l -> occursin("Dropping", l.message), logs)
            @test isempty(drop_logs)
            @test length(samples) == 10
        end
    end
end

main()
