using Test
using Distributed

"""
Unit tests for log_memory_status (src/results/Memory_log.jl).

Usage, from a REPL that already ran include("experiments/MGA_Tests.jl"):
  test_memory_log_units()
"""

isdefined(Main, :run_stochastic_exploration_risk_pareto) || include(joinpath(@__DIR__, "..", "experiments", "MGA_Tests.jl"))

@testset "Memory log" begin
    @testset "Main process only" begin
        r = log_memory_status("unit test, main only"; include_workers = false)
        @test r.main.rss > 0 && r.main.peak_rss > 0 && r.main.live > 0
        @test isempty(r.workers) && r.worker_rss_total == 0
        r2 = log_memory_status("unit test, second call"; include_workers = false)
        @test 0 <= r2.gc_share <= 1           # share of wall time since the first call
    end

    @testset "With workers" begin
        added = nworkers() < 2 || workers() == [myid()] ? addprocs(2) : Int[]
        try
            configure_parallel_workers!(Dict{String,Any}("Parallel flag" => true, "Workers" => nworkers()))
            @test all(remotecall_fetch(() -> LinearAlgebra.BLAS.get_num_threads(), pid) == 1 for pid in workers())
            r = log_memory_status("unit test, with workers")
            @test length(r.workers) == nworkers() == length(r.pids)
            @test all(m.rss > 0 && m.live > 0 for m in r.workers)
            @test r.worker_rss_total >= r.worker_rss_max > 0
            @test r.worker_rss_max_pid in workers()
        finally
            isempty(added) || rmprocs(added)
        end
    end

    @testset "Allocation-gated forced GC" begin
        added = nworkers() < 2 || workers() == [myid()] ? addprocs(2) : Int[]
        try
            configure_parallel_workers!(Dict{String,Any}("Parallel flag" => true, "Workers" => nworkers()))
            @test gc_all_processes!(Dict{String,Any}()) == nprocs()                                   # unset: always
            @test gc_all_processes!(Dict{String,Any}("Worker GC threshold MB" => nothing)) == nprocs()
            gated = Dict{String,Any}("Worker GC threshold MB" => 64)
            gc_all_processes!(Dict{String,Any}())                                                     # reset the counters
            @test gc_all_processes!(gated) == 0                                                       # nothing allocated since
            junk = [zeros(10^6) for _ in 1:10]                                                        # 80 MB on the main process
            junk = nothing
            @test gc_all_processes!(gated) == 1                                                       # only main collects
            remotecall_wait(() -> (global _gc_test_junk = [zeros(10^6) for _ in 1:10]; _gc_test_junk = nothing; nothing), workers()[1])
            @test gc_all_processes!(gated) == 1                                                       # only that worker
            @test gc_all_processes!(gated) == 0
        finally
            isempty(added) || rmprocs(added)
        end
    end

    @testset "Line format" begin
        buf = IOBuffer()
        Base.CoreLogging.with_logger(Base.CoreLogging.SimpleLogger(buf)) do
            log_memory_status("format check"; include_workers = false)
        end
        line = String(take!(buf))
        @test occursin("MEM|format check|main rss", line) && occursin("live", line)
        @test occursin("|free ", line) && occursin("main GC", line)
    end
end
