using Test
using Distributed
using Logging

"""
Tests for the "Log style" setting and the run log (src/results/Run_log.jl).

Usage, from a REPL that already ran include("experiments/MGA_Tests.jl"):
  test_run_log_units()                # unit tests + a small della/laptop end-to-end run
  include("tests/run_log_tests.jl")   # the same

The end-to-end part runs the risk-pareto driver on 2 scenarios with 1 MGA vector (a minute or two).
Set RUN_LOG_E2E = false before including to skip it.
"""

isdefined(Main, :run_stochastic_exploration_risk_pareto) || include(joinpath(@__DIR__, "..", "experiments", "MGA_Tests.jl"))

# Messages logged by `f` under the `style` logger, as text
function captured_log(f, style::AbstractString; start = time())
    buf = IOBuffer()
    logger = style == "della" ? CompactLogger(buf; start = start) : run_logger(style, buf)
    with_logger(f, logger)
    return String(take!(buf))
end

fake_sep(restored, rounds) = (restored = restored, rounds = rounds, t_check = 0.0, t_extra_mp = 0.0)

@testset "Run log" begin
    @testset "Log style setting" begin
        @test run_log_style(Dict{String,Any}()) == "laptop"
        @test run_log_style(Dict{String,Any}("Log style" => "della")) == "della"
        @test_throws ErrorException run_log_style(Dict{String,Any}("Log style" => "cluster"))
        folder = joinpath(@__DIR__, "..", "inputs", "Inputs_2z_7tech")
        @test run_log_style(load_settings(folder)) == "laptop"
        @test run_log_style(load_settings(folder; profile = "della")) == "della"
    end

    @testset "Laptop logger shows detail as Info" begin
        text = captured_log("laptop") do
            log_detail("detail message")
            @info("info message")
        end
        @test occursin("Info: detail message", text)
        @test occursin("Info: info message", text)
        @test !occursin("LogLevel", text)
        @test !occursin("Debug", captured_log(() -> @debug("hidden"), "laptop"))
    end

    @testset "Della logger: stamped single lines, detail hidden" begin
        text = captured_log("della"; start = time() - 3725) do
            log_detail("detail message")
            @info("info message")
        end
        @test !occursin("detail message", text)
        @test text == "[+01:02:05] info message\n"

        text = captured_log("della") do
            try
                error("boom")
            catch e
                @warn("it failed", exception = (e, catch_backtrace()))
            end
        end
        lines = split(rstrip(text), '\n')
        @test startswith(lines[1], "[+00:00:00] ┌ Warning: it failed")
        @test occursin("boom", text)
        @test all(startswith(l, " "^12) for l in lines[2:end]) # continuation lines indented under the stamp
    end

    @testset "ITER| and SOLVE| lines" begin
        anchor = iter_line("Risk_Weight_0.5", 3; gap = 0.0531, LB = 51.0021, UB = 53.7104, time_mp = 0.114, time_sp = 2.981,
                           sep = fake_sep(405, 1), live = 1234, time_reg = 0.4)
        @test anchor == "ITER| Risk_Weight_0.5 | it 3 | gap 5.31% | LB 51.0021 UB 53.7104 | MP 0.114 s (sep +405 in 1 rnd) | reg 0.400 s | SP 2.981 s | live 1234"

        mga = iter_line("Random_1_Budget_1_CVaR_Cap_1", 2; gap = ["Transformed" => 0.0005, "CVaR" => 0.0016],
                        UB = ["Transformed" => 51.1242, "CVaR" => 100971.53], time_mp = 0.15, time_sp = 2.9, sep = fake_sep(4, 1), live = 150)
        @test mga == "ITER| Random_1_Budget_1_CVaR_Cap_1 | it 2 | gap T 0.05% CVaR 0.16% | UB T 51.1242 CVaR 100971.5 | MP 0.150 s (sep +4 in 1 rnd) | SP 2.900 s | live 150"

        out = Dict{String,Any}("Gaps" => [[0.084, 0.01], [0.0011, 0.0004]], "Time SP hist" => [3.1, 3.0], "Time MP hist" => [0.2, 0.1],
                               "Time final" => 1.5, "Elapsed" => 8.0, "Sep hist" => [fake_sep(2400, 5), fake_sep(12, 1)],
                               "Live cuts at start" => 935, "MP" => Dict("Inv_cost" => 6.125e9), "Expected Value" => 5.616e9, "CVaR" => 1.042e10)
        budgets = Dict("Transformed" => Dict("budget" => 52.097, "transform_cvar" => 0.1, "transform_sys" => 0.1), "System_Expected" => 11.59)
        line = solve_line("MGA 1/50 L1/6", "Random_1_Budget_1", out; pool_size = 9800, budgets = budgets)
        @test line == "SOLVE| MGA 1/50 L1/6 | Random_1_Budget_1 | budget SysEV≤11.59 T≤52.097 | 2 it, gap 0.11% | " *
                      "wall 8.0 s (MP 0.3, SP 6.1, final 1.5) | cuts 935 live at start, +2.4k restored, pool 9.8k | inv 6.125 EV 5.616 CVaR 10.42 (\$B/yr)"

        out_anchor = merge(out, Dict{String,Any}("Gaps" => [0.4, 0.0008], "Time Reg hist" => [0.4, 0.0]))
        line = solve_line("A 1/11 w=0.5", "Risk_Weight_0.5", out_anchor; pool_size = 10, weight = 0.5)
        @test occursin("(MP 0.3, reg 0.4, SP 6.1, final 1.5)", line)
        @test occursin("2 it, gap 0.08%", line)
        @test occursin("risk-adj 14.14", line) # 6.125 + 0.5 * 5.616 + 0.5 * 10.42
    end

    @testset "parse_run_log round-trips both line types" begin
        path = tempname()
        open(path, "w") do io
            with_logger(CompactLogger(io; start = time() - 61)) do
                @info(iter_line("Random_1_Budget_1", 1; gap = ["Transformed" => 0.084], UB = ["Transformed" => 56.47], time_mp = 0.19,
                                time_sp = 3.14, sep = fake_sep(2400, 5), live = 3300))
                @info(iter_line("Risk_Weight_0.5", 2; gap = 0.12, LB = 47.1, UB = 52.7, time_mp = 0.1, time_sp = 1.66,
                                sep = fake_sep(0, 0), live = 125, time_reg = 0.38))
                @info("=== something else")
                @info("SOLVE| A 1/11 w=0.5 | Risk_Weight_0.5 | 14 it, gap 0.08% | wall 18.7 s (MP 1.2, reg 0.4, SP 15.3, final 1.8) | cuts 0 live at start, +0 restored, pool 1.8k | inv 6.1 EV 5.6 CVaR 10 risk-adj 14 (\$B/yr)")
            end
        end
        laptop_line = "[ Info: ITER| Risk_Weight_0.5 | it 3 | gap 1.1% | LB 47.1 UB 52.7 | MP 0.10 s (sep +3 in 1 rnd) | SP 1.60 s | live 130"
        open(io -> println(io, laptop_line), path, "a")

        log = parse_run_log(path)
        @test nrow(log.iters) == 3 && nrow(log.solves) == 1
        r = log.iters[1, :]
        @test (r.elapsed_s, r.label, r.iteration, r.gap_pct, r.mp_s, r.sep_restored, r.sep_rounds, r.sp_s, r.live) ==
              (61, "Random_1_Budget_1", 1, 8.4, 0.19, 2400, 5, 3.14, 3300)
        @test ismissing(r.reg_s)
        @test log.iters[2, :reg_s] == 0.38
        @test ismissing(log.iters[3, :elapsed_s]) && log.iters[3, :iteration] == 3
        s = log.solves[1, :]
        @test (s.tag, s.label, s.iterations, s.gap_pct, s.wall_s, s.mp_s, s.reg_s, s.sp_s, s.final_s) ==
              ("A 1/11 w=0.5", "Risk_Weight_0.5", 14, 0.08, 18.7, 1.2, 0.4, 15.3, 1.8)
    end

    @testset "show_conflicts: names only in della" begin
        m = Model()
        @variable(m, x >= 0)
        @constraint(m, c_test, x <= 1)
        text = captured_log(() -> show_conflicts([c_test, LowerBoundRef(x)]), "della")
        @test occursin("Conflicting constraints (2): c_test, x ≥ 0", text)
        @test !occursin("x ≤ 1", text) # the expression of a named constraint isn't printed
    end
end

#==========================================================================================
End to end: the risk-pareto driver on a tiny problem, in each style
==========================================================================================#

# Runs the driver with `style`, its log going to a file; returns the log text
function run_log_e2e(style::AbstractString; capacity_exploration::Bool)
    folder = joinpath(@__DIR__, "..", "inputs", "Inputs_2z_7tech")
    settings = load_settings(folder)
    settings["Parallel flag"] = false
    settings["Scenarios"]["Gas"]["Count"] = 2
    settings["Scenarios"]["Weather"]["Count"] = 1
    settings["Scenarios"]["Demand"]["Count"] = 1
    settings["Log style"] = style
    inputs = load_input_data(folder, settings)
    SPs = build_all_subproblems(inputs, settings)
    out_dir = mktempdir()
    log_path = joinpath(out_dir, "run.log")
    previous = global_logger()
    try
        open(log_path, "w") do io
            configure_run_output!(settings; stream = io)
            run_stochastic_exploration_risk_pareto(SPs, inputs, settings, joinpath(out_dir, "results"), joinpath(out_dir, "summary");
                run_settings = Dict("MGA vectors" => 1, "Capacity exploration" => capacity_exploration, "Summary name" => "run_log_e2e"),
                budget_settings = Dict("Tighten budget" => true, "Tighten step size" => 0.5),
                cap_settings = Dict("Number of caps" => 1, "Cap levels" => [1]))
        end
    finally
        global_logger(previous)
    end
    return log_path
end

if !isdefined(Main, :RUN_LOG_E2E) || Main.RUN_LOG_E2E
    @testset "Run log end to end" begin
        @testset "della" begin
            path = run_log_e2e("della"; capacity_exploration = true)
            text = read(path, String)
            lines = filter(!isempty, split(text, '\n'))
            @test all(l -> startswith(l, "[+") || startswith(l, " "), lines) # every line stamped or a continuation
            for hidden in ("Running economic dispatch", "*** Iteration", "Added budget constraint", "Objective set", "New capacity mix", "Cut separation")
                @test !occursin(hidden, text)
            end
            for shown in ("RUN| started", "RUN| Summary", "=== Anchors: 11 risk weights", "=== MGA vector 1/1", "=== Final summary write", "=== Done in")
                @test occursin(shown, text)
            end

            log = parse_run_log(path)
            @test nrow(log.solves) >= 12 # 11 anchors + at least the first budget level
            @test allunique(log.solves.label)
            for s in eachrow(log.solves)
                its = log.iters[log.iters.label .== s.label, :]
                @test nrow(its) == s.iterations
                @test its.iteration == 1:s.iterations
                # Per-iteration times add up to the solve's totals (rounded to 0.001 s / 0.1 s)
                @test isapprox(sum(its.mp_s), s.mp_s; atol = 0.05 + 0.0005 * nrow(its))
                @test isapprox(sum(its.sp_s), s.sp_s; atol = 0.05 + 0.0005 * nrow(its))
                @test s.mp_s + s.sp_s + coalesce(s.reg_s, 0.0) + s.final_s <= s.wall_s + 0.2
            end
            @test any(startswith.(log.solves.tag, "MGA 1/1 L1/"))
            @test any(occursin.("cap 1/1", log.solves.tag))
            println("\n---- della log (first 40 lines) ----")
            foreach(println, first(lines, 40))
        end

        @testset "laptop keeps the detail lines" begin
            text = read(run_log_e2e("laptop"; capacity_exploration = false), String)
            for shown in ("Running economic dispatch", "*** Iteration", "Objective set to minimize weighted CVaR", "ITER| Risk_Weight_0.5", "SOLVE| A 1/11")
                @test occursin(shown, text)
            end
            @test !occursin("LogLevel", text)
        end
    end
end
