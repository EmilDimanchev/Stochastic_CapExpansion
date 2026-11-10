using Test
using Distributed
using LinearAlgebra, Random, Statistics, DataFrames, JuMP
import HDF5

"""
Solver-free unit tests for the cut export (src/results/Export_cuts.jl).

Usage:
  julia --project=. tests/cut_export_tests.jl
or, from a REPL that already ran include("experiments/MGA_Tests.jl"):
  test_cut_export_units()

Synthetic data (tests/sampling_toy_helpers.jl): each scenario's op cost is convex piecewise-
linear in capacity with known gradients, so every bound can be checked against the true value.
"""

isdefined(Main, :run_stochastic_exploration_risk_pareto) || include(joinpath(@__DIR__, "..", "experiments", "MGA_Tests.jl"))
include(joinpath(@__DIR__, "sampling_toy_helpers.jl"))   # ToyCost, toy_cost, toy_eval

const CX_D = 4        # first-stage variables
const CX_ALPHA = 0.1
const CX_FOLDER = mktempdir()

const CX_SETTINGS = Dict{String,Any}("Export cuts" => true, "3-parent blends" => 1500, "2-parent blends" => 400, "Sparse blends" => 200,
                                     "Max sparse parents" => 6, "Held-out blends" => 500, "Max repair rounds" => 10, "Clean rounds required" => 2, "Gap tolerance" => 1e-4,
                                     "Near-active tolerance" => 1e-5, "Duplicate tolerance" => 1e-8, "Float32 tolerance" => 1e-5,
                                     "Validity tolerance" => 1e-6, "Sparse density threshold" => 0.3, "Target MB" => 25, "Hard cap MB" => 50)
cx_settings(; kwargs...) = merge(CX_SETTINGS, Dict{String,Any}(string(k) => v for (k, v) in kwargs))

# True cost of scenario s at each column of X
toy_eval_all(tc, X, s) = [toy_eval(tc, X[:, j])[1][s] for j in axes(X, 2)]

cx_solve(name, x, tc, p) = (tv = toy_eval(tc, x); FinalSolve(name, copy(x), tv[1], tv[2], p, 0.0, true))

function cx_bank(; seed = 3, n = 40, S = 6, pieces = 25, spread = 30.0)
    rng = Xoshiro(seed)
    tc = toy_cost(rng, S, pieces, CX_D)
    p = fill(1 / S, S)
    return [cx_solve("Plan_$i", spread .* rand(rng, CX_D), tc, p) for i in 1:n], tc, p
end

cx_input(bank; pool_cuts = [], skipped = 0, picks = GridFillPick[], seed = 7) = build_cut_export_input(pool_cuts, skipped, bank, picks;
    run_id = "toy", var_name = ["v$i" for i in 1:length(bank[1].x)], var_resource = ["r$i" for i in 1:length(bank[1].x)], var_zone = fill("z1", length(bank[1].x)),
    var_units = fill("MW", length(bank[1].x)), scenario_label = ["D_$(s)_F_1_W_1" for s in 1:length(bank[1].θ)], cvar_alpha = CX_ALPHA, seed = seed, scaling_factor_cost = 1e5)

# Extra cuts for the pool, on top of the ones the plans give: tangents at random capacities
function cx_pool_cuts(tc, p; n = 60, seed = 9)
    rng = Xoshiro(seed)
    S = length(p)
    cuts = NamedTuple{(:group, :a, :g, :name), Tuple{Int, Float64, Vector{Float64}, String}}[]
    for j in 1:n
        x = 30 .* rand(rng, CX_D)
        θ, β = toy_eval(tc, x)
        for s in 1:S
            push!(cuts, (group = s, a = θ[s] - dot(β[:, s], x), g = β[:, s], name = "optimality_cut_Risk_Weight_0.5_$(j)_[$s,1,1,1]"))
        end
    end
    return cuts
end

@testset "Cut export" begin

    @testset "Deduplication drops only duplicates and dominated intercepts" begin
        beta = [1.0 1.0 2.0 2.0 2.0 3.0;
                0.0 0.0 1.0 1.0 1.0 1.0]
        alpha = [5.0, 7.0, 3.0, 3.0, 3.0 + 1e-12, 1.0]
        kept = dedupe_cuts(alpha, beta, ones(2), 1e-8)
        @test kept == [2, 5, 6]     # 7 beats 5 at the same slope; cuts 3-5 share a slope; cut 6 stays

        # A real difference in slope or intercept must survive
        beta2 = [1.0 1.0 1.0001; 0.0 0.0 0.0]
        @test dedupe_cuts([100.0, 200.0, 100.0], beta2, ones(2), 1e-8) == [2, 3]

        # Near-duplicates in (α, β·xscale), differing by much less than the tolerance
        beta3 = [1.0 1.0 + 1e-13; 2.0 2.0]
        @test length(dedupe_cuts([50.0, 50.0], beta3, ones(2), 1e-8)) == 1
    end

    @testset "Probe blends stay inside the plans' convex hull" begin
        X = [0.0 10.0 0.0 4.0; 0.0 0.0 10.0 4.0]    # n = 2, N = 4
        cfg = cut_export_config(cx_settings())
        rng = Xoshiro(1)
        probes = build_probe_blends(X, cfg, 1.0, rng)
        @test size(probes, 2) == 4 + 1500 + 400 + 200
        @test probes[:, 1:4] == X                    # every plan is a probe
        @test all(probes .>= -1e-12) && all(sum(probes; dims = 1) .<= 10 + 1e-9)
        @test size(build_probe_blends(X, cfg, 0.5, Xoshiro(1)), 2) == 4 + 750 + 200 + 100
        sparse_points = random_sparse_blends(X, 50, 6, rng)
        @test size(sparse_points) == (2, 50)
    end

    @testset "Chunked evaluation matches the direct computation" begin
        rng = Xoshiro(2)
        alpha, beta, X = rand(rng, 200), randn(rng, CX_D, 200), rand(rng, CX_D, 1300)
        direct = alpha .+ beta' * X
        values, cuts = best_cut_at(alpha, beta, X)
        @test values ≈ vec(maximum(direct; dims = 1))
        @test [direct[cuts[j], j] for j in axes(X, 2)] ≈ values
        @test findall(active_cut_mask(alpha, beta, X, 0.0)) == sort(unique(cuts))
        near = 0.01
        expected = [any(direct[k, j] >= values[j] - near * max(abs(values[j]), 1.0) for j in axes(X, 2)) for k in eachindex(alpha)]
        @test active_cut_mask(alpha, beta, X, near) == expected
    end

    @testset "Pruning: fewer cuts, a valid bound, gap below tolerance at fresh blends" begin
        bank, tc, p = cx_bank()
        input = cx_input(bank; pool_cuts = cx_pool_cuts(tc, p))
        cfg = cut_export_config(cx_settings())
        X = input.plan_x
        probes = build_probe_blends(X, cfg, 1.0, Xoshiro(1))
        held_rng = Xoshiro(2)
        sets = [build_held_out_blends(X, cfg, held_rng) for _ in 1:cfg.max_repair_rounds]
        test = build_held_out_blends(X, cfg, held_rng)
        kept, coverage = prune_cut_table(input.cuts, probes, sets, test, ones(CX_D), cfg.near_active_tol, cfg.gap_tol, cfg.duplicate_tol, cfg.clean_rounds)
        @test sum(length, kept) < n_table_cuts(input.cuts)
        @test coverage.max_gap <= cfg.gap_tol
        fresh = random_sparse_blends(X, 300, 6, Xoshiro(5))
        for s in eachindex(kept)
            pruned, _ = best_cut_at(input.cuts.alpha[s][kept[s]], input.cuts.beta[s][:, kept[s]], fresh)
            full, _ = best_cut_at(input.cuts.alpha[s], input.cuts.beta[s], fresh)
            @test all(pruned .<= full .+ 1e-9)                     # a subset can't exceed the full max
            @test all(pruned .<= toy_eval_all(tc, fresh, s) .+ 1e-9)   # and stays below the true cost
        end
    end

    @testset "Repair recovers a cut the probes missed" begin
        # Two cuts in one variable. The second is only the max for x > 5. Probes sit at x <= 1.
        alpha, beta = [10.0, -40.0], [0.0 10.0]
        probes = reshape(collect(0.0:0.1:1.0), 1, :)
        missed = falses(2)
        missed[active_cut_mask(alpha, beta, probes, 0.0)] .= true
        @test missed == [true, false]
        held_out = [reshape([8.0, 9.0], 1, :)]
        kept, rounds = prune_scenario(alpha, beta, probes, held_out, ones(1), 0.0, 1e-4, 1e-8, 1)
        @test kept == [1, 2] && rounds == 1
        # One clean set stops the repair only if one is all that's required; a later set can
        # still expose a gap when more clean sets are required
        clean_then_gap = [reshape([0.5], 1, :), reshape([9.0], 1, :), reshape([0.5], 1, :)]
        kept2, rounds2 = prune_scenario(alpha, beta, probes, clean_then_gap, ones(1), 0.0, 1e-4, 1e-8, 1)
        @test rounds2 == 1 && kept2 == [1]
        kept3, rounds3 = prune_scenario(alpha, beta, probes, clean_then_gap, ones(1), 0.0, 1e-4, 1e-8, 2)
        @test rounds3 == 3 && kept3 == [1, 2]
    end

    @testset "Collection: units, order, cut names and batches" begin
        bank, tc, p = cx_bank(n = 3)
        @test parse_cut_name("optimality_cut_Random_1_Budget_2_CVaR_Cap_1_14_[3,1,1]") == ("Random_1_Budget_2_CVaR_Cap_1", 14)
        @test parse_cut_name("optimality_cut_Risk_Weight_0.5_2_[1,1,1]") == ("Risk_Weight_0.5", 2)
        pool = cx_pool_cuts(tc, p; n = 2)
        input = cx_input(bank; pool_cuts = pool)
        S = length(p)
        @test all(length.(input.cuts.alpha) .== 2 + 3)
        @test input.sources == ["Risk_Weight_0.5", "Plan_1", "Plan_2", "Plan_3"]
        @test input.cuts.iteration[1] == [1, 2, 0, 0, 0] && input.cuts.batch[1] == zeros(Int, 5)
        # The tangent built from a plan is exact at that plan
        for s in 1:S, (j, fs) in enumerate(bank)
            k = 2 + j
            @test input.cuts.alpha[s][k] + dot(input.cuts.beta[s][:, k], fs.x) ≈ fs.θ[s]
        end
        # A failed plan contributes no cuts but stays in the plan list
        bank[2].ok = false
        input2 = cx_input(bank)
        @test all(length.(input2.cuts.alpha) .== 2)
        @test input2.plan_ok == [true, false, true] && length(input2.plan_ids) == 3
    end

    @testset "Pool reader matches the master's rows" begin
        settings = load_settings(joinpath(@__DIR__, "..", "inputs", "Inputs_2z_7tech"))
        scale = settings["Scaling factor cost"]
        S, F = 3, 2
        bank, _, _ = cx_bank(S = S * F, n = 4)
        MP = Model()
        @variable(MP, x[1:3] >= 0); @variable(MP, x_line[1:1] >= 0)
        @variable(MP, alpha[1:S, 1:F, 1:1, 1:1] >= 0)     # one block per scenario
        pool = CutPool(settings)
        for (it, fs) in enumerate(bank)
            θ = reshape(fs.θ ./ scale, S, F, 1, 1); β = reshape(fs.β ./ scale, 4, S, F, 1, 1)
            xp = fs.x
            coeff(s) = (it == 4 && s == 1) ? 0 : 1      # the last solve's first scenario "failed"
            opt = @constraint(MP, [s in 1:S, f in 1:F, k in 1:1, b in 1:1], coeff(s) * alpha[s,f,k,b] >= θ[s,f,k,b] + sum(β[r,s,f,k,b]*(x[r] - xp[r]) for r in 1:3) + β[4,s,f,k,b]*(x_line[1] - xp[4]), base_name = "optimality_cut_t_$(it)_")
            add_cuts!(pool, Dict(name(c) => c for c in vec(opt)))
        end
        _archive!(pool, MP, [n for n in keys(pool.live) if occursin("_t_2_", n)])
        cuts, skipped = pool_optimality_cuts(pool, MP, settings)
        @test skipped == F                      # the scenarios with s == 1 in solve 4
        @test length(cuts) == 4 * S * F - F
        for cut in cuts
            _, it = parse_cut_name(cut.name)
            fs = bank[it]
            @test cut.a + dot(cut.g, fs.x) ≈ fs.θ[cut.group] rtol = 1e-9     # unscaled, right scenario, right column order
            @test cut.g ≈ fs.β[:, cut.group] rtol = 1e-9
        end
    end

    @testset "File: layout, CSR and dense give the same bound, float32 falls back" begin
        bank, tc, p = cx_bank()
        input = cx_input(bank; pool_cuts = cx_pool_cuts(tc, p))
        out = export_cut_archive(input, joinpath(CX_FOLDER, "layout"), cx_settings(); file_stem = "a")
        @test isfile(out.path) && !isfile(joinpath(CX_FOLDER, "layout", "a_raw.jls"))
        HDF5.h5open(out.path, "r") do f
            K = length(read(f["cuts/alpha"]))
            @test eltype(f["cuts/scenario"]) == Int16 && eltype(f["cuts/alpha"]) == Float64
            @test size(read(f["cuts/beta"])) == (CX_D, K)       # K×n in h5py, column-major here
            @test read(f["meta/beta_format"]) == "dense"
            @test read(f["plans/id"]) == input.plan_ids
            @test read(f["vars/name"]) == input.var_name
            @test read(f["meta/raw_cuts_total"]) == n_table_cuts(input.cuts)
            @test read(f["meta/pruned_cuts_total"]) == K
            @test read(f["meta/scenario_probability"]) == input.p
            @test !isempty(read(f["meta/validation_report"])) && !isempty(read(f["meta/budget_tradeoff"]))
            @test length(HDF5.get_create_properties(f["cuts/beta"]).filters) == 2     # shuffle + gzip
        end

        # Same export with CSR forced on must give the same bounds
        csr = export_cut_archive(input, joinpath(CX_FOLDER, "csr"), cx_settings(; var"Sparse density threshold" = 2.0); file_stem = "a")
        HDF5.h5open(f -> @test(read(f["meta/beta_format"]) == "csr"), csr.path, "r")
        Y = input.plan_x
        @test cut_values(CutBank(csr.report.written.cuts.alpha, csr.report.written.cuts.beta), Y) ≈ cut_values(CutBank(out.report.written.cuts.alpha, out.report.written.cuts.beta), Y)

        # Large capacities close together under an affine cost: α and β·x nearly cancel and LB is
        # small, so Float32 rounding of either β or x is a large relative change
        big_beta = [1234567.3, -1000000.7, 500000.3]
        x0 = [1.1e6, 1.3e6, 0.7e6]
        big_p = [0.5, 0.5]
        big = map(1:5) do i
            x = x0 .+ 1e-3 .* rand(Xoshiro(i), 3)
            FinalSolve("Big_$i", x, [1e4, 2e4] .+ dot(big_beta, x .- x0), [big_beta big_beta], big_p, 0.0, true)
        end
        big_input = cx_input(big)
        kept = [collect(eachindex(a)) for a in big_input.cuts.alpha]
        @test float32_lb_change(big_input.cuts, kept, big_input.plan_x; round_beta = true, round_x = true) > 1e-5
        choice = choose_float_types(big_input.cuts, kept, big_input.plan_x, 1e-5)
        @test choice.beta == Float64 && choice.x == Float64
        @test float32_lb_change(big_input.cuts, kept, big_input.plan_x; round_beta = false, round_x = false) == 0
        # Ordinary numbers stay in Float32
        ordinary = choose_float_types(input.cuts, [collect(eachindex(a)) for a in input.cuts.alpha], input.plan_x, 1e-5)
        @test ordinary.beta == Float32 && ordinary.x == Float32
    end

    @testset "Blocks: the file stores scenario and block, and floors sum the block maxima" begin
        S, B = 4, 3
        rng = Xoshiro(21)
        tc = toy_cost(rng, S * B, 15, CX_D)
        p = fill(1 / S, S)
        bank = [toy_block_solve("Plan_$i", 30 .* rand(rng, CX_D), tc, p, B) for i in 1:25]
        input = cx_input(bank)
        @test input.cuts.n_blocks == B && length(input.cuts.alpha) == S * B
        out = export_cut_archive(input, joinpath(CX_FOLDER, "blocks"), cx_settings(); file_stem = "b")
        scenario, block, alpha, beta = HDF5.h5open(out.path, "r") do f
            @test read(f["meta/n_blocks"]) == B
            read(f["cuts/scenario"]), read(f["cuts/block"]), read(f["cuts/alpha"]), Float64.(read(f["cuts/beta"]))
        end
        @test sort(unique(block)) == 1:B && sort(unique(scenario)) == 1:S
        # The floor a reader computes from the raw datasets: sum over blocks of each block's max
        file_floor(X, s) = sum(vec(maximum(alpha[r] .+ beta[:, r]' * X; dims = 1)) for r in (findall((scenario .== s) .& (block .== b)) for b in 1:B))
        X = input.plan_x
        Y = 30 .* rand(Xoshiro(5), CX_D, 200)
        for s in 1:S
            @test file_floor(X, s) ≈ input.plan_θ[s, :] rtol = 1e-3       # pruned to the 1e-4 gap tolerance per block
            truth = [sum(toy_eval(tc, Y[:, j])[1][s .+ (0:B-1) .* S]) for j in axes(Y, 2)]
            @test all(file_floor(Y, s) .<= truth .+ 1e-6 .* abs.(truth))      # valid elsewhere, up to Float32 storage of beta
        end
        # Unpruned, the plans' own block tangents make the floor exact there; the package's reader agrees with the raw one
        @test cut_values(table_cut_bank(input.cuts), X) ≈ input.plan_θ rtol = 1e-9
        @test cut_values(table_cut_bank(out.report.written.cuts), X) ≈ reduce(vcat, [file_floor(X, s)' for s in 1:S]) rtol = 1e-9
        summed = CutBank(S, CX_D); foreach(fs -> add_final_solve!(summed, FinalSolve(fs.name, fs.x, fs.θ, fs.β, fs.p, 0.0, true)), bank)
        @test all(cut_values(table_cut_bank(input.cuts), Y) .>= cut_values(summed, Y) .- 1e-9)
    end

    @testset "Budget: shrink order, tradeoff reported, and a cap that can't be met writes nothing" begin
        bank, tc, p = cx_bank()
        input = cx_input(bank; pool_cuts = cx_pool_cuts(tc, p))
        steps = budget_steps(cut_export_config(cx_settings()))
        @test first(steps).near_active_tol == 1e-5 && steps[2].near_active_tol == 0.0   # shrinking lowers it
        @test issorted([s.density for s in steps]; rev = true) && last(steps).gap_tol == 1e-3
        @test all(s.gap_tol == 1e-4 for s in steps[1:end-1])    # the gap tolerance goes last

        over_target = export_cut_archive(input, joinpath(CX_FOLDER, "target"), cx_settings(; var"Target MB" = 1e-6, var"Hard cap MB" = 50); file_stem = "a")
        @test length(over_target.tradeoff) == length(steps)     # tried every step, kept the smallest
        @test isfile(over_target.path) && !isfile(over_target.path * ".tmp")

        folder = joinpath(CX_FOLDER, "cap")
        @test_throws ErrorException export_cut_archive(input, folder, cx_settings(; var"Target MB" = 1e-6, var"Hard cap MB" = 1e-6); file_stem = "a")
        @test !isfile(joinpath(folder, "a.h5")) && !isfile(joinpath(folder, "a.h5.tmp"))
        @test isfile(joinpath(folder, "a_raw.jls"))              # the raw copy survives for a retry

        # ...and that retry works without the model
        retry_out = export_cut_archive_from_raw(joinpath(folder, "a_raw.jls"), folder, cx_settings(); file_stem = "a")
        @test isfile(retry_out.path) && !isfile(joinpath(folder, "a_raw.jls"))
    end

    @testset "Validation: replay matches the sampler exactly, bounds are valid and tight" begin
        # A real GridFill loop on the toy model, so the logged bounds come from the sampler itself
        rng = Xoshiro(11)
        S, n0 = 20, 80
        tc = toy_cost(rng, S, 4, CX_D)
        p = fill(1 / S, S)
        cinv = [0.9, 0.6, 1.1, 0.4]
        solve(name, x) = (tv = toy_eval(tc, x); FinalSolve(name, copy(x), tv[1], tv[2], p, dot(cinv, x), true))
        bank = [solve("Random_$i", 30 .* rand(rng, CX_D)) for i in 1:n0]
        C(fs) = fs.inv_cost + solve_ev(fs)
        V(fs) = fs.inv_cost + solve_cvar(fs, CX_ALPHA)
        iE, iR = argmin(C.(bank)), argmin(V.(bank))
        iM = argmin(C.(bank) .+ V.(bank))
        bank[iE].name, bank[iR].name = "Risk_Weight_1.0", "Risk_Weight_0.0"
        iM in (iE, iR) || (bank[iM].name = "Risk_Weight_0.5")
        sampling = Dict{String,Any}("Samples" => 24, "Batch size" => 8, "Candidates" => 30, "Write gradients" => false,
                                    "GridFill" => Dict{String,Any}("Cell size" => 0.1, "Delta cap" => 1.1, "Per cell" => 2))
        evaluator = (picks, SPs, inputs, settings, g) -> begin
            solved = [solve(pk.name, pk.x) for pk in picks]
            rows = [DataFrame(Name = [pk.name]) for pk in picks]
            rows, copy(rows), copy(rows), solved
        end
        tuning = merge(GRIDFILL_TUNING, (triples_per_aim = 3000, prune_blends = 300))
        picks = GridFillPick[]
        run_gridfill_sampling!(bank, nothing, nothing, nothing, Dict{String,Any}("Resources" => ["r1", "r2", "r3"], "Number of lines" => 1),
                               Dict{String,Any}("Value-at-Risk percent" => CX_ALPHA, "Seed" => 7, "Convergence tolerance" => 1e-9), sampling;
                               tuning = tuning, evaluate = evaluator, parallel = false, pick_log = picks)
        @test length(picks) > 8 && length(unique(pk.batch for pk in picks)) >= 2

        input = cx_input(bank; picks = picks, seed = 7)
        @test all(input.plan_batch[n0+1:end] .>= 1) && all(input.plan_batch[1:n0] .== 0)
        replay = replay_logged_bounds(input; prune_blends = tuning.prune_blends)
        @test replay.lb_ev == [pk.lb_ev for pk in picks]         # exact, not approximate
        @test replay.lb_cvar == [pk.lb_cvar for pk in picks]

        out = export_cut_archive(input, joinpath(CX_FOLDER, "validate"), cx_settings(); file_stem = "a", prune_blends = tuning.prune_blends)
        @test out.report.replay_diff.ev == (0.0, 0.0) && out.report.replay_diff.cvar == (0.0, 0.0)
        @test out.report.violations_raw == 0 && out.report.violations_file == 0
        @test out.report.max_slack < 1e-5                        # each plan's own cut is exact at that plan

        # The full-archive floor at a pick's batch is at least the sampler's pruned floor...
        full = pick_floors(input.cuts, input)
        @test all(full.ev .>= [pk.lb_ev for pk in picks] .- 1e-9)
        # ...and the file's batch filter removes the cuts generated at later batches
        written = out.report.written.cuts
        first_batch = minimum(pk.batch for pk in picks)
        @test all(all(cut_bank_before_batch(written, first_batch).a[s] .== written.alpha[s][written.batch[s] .< first_batch]) for s in eachindex(written.alpha))
    end

    @testset "An inexact cut above a plan's cost is reported" begin
        bank, tc, p = cx_bank(n = 8)
        bank[3].θ[2] -= 50.0      # the plan now looks cheaper than the tangents say
        input = cx_input(bank)
        out = export_cut_archive(input, joinpath(CX_FOLDER, "violate"), cx_settings(); file_stem = "a")
        @test out.report.violations_raw >= 1
    end

    @testset "Settings: a missing key fails loudly" begin
        @test_throws KeyError cut_export_config(delete!(copy(CX_SETTINGS), "Gap tolerance"))
    end
end
