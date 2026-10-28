using Test
using Distributed
using LinearAlgebra, Random, Statistics, DataFrames, JuMP

"""
Solver-free unit tests for the GridFill sampler (src/model/boundary_sampling/gridfill_*.jl).

Usage:
  julia --project=. tests/gridfill_tests.jl
or, from a REPL that already ran include("experiments/MGA_Tests.jl"):
  test_gridfill_units()

Synthetic data (tests/sampling_toy_helpers.jl): each scenario's op cost is convex piecewise-
linear in capacity with known gradients; investment cost is linear. The anchors are the
cheapest sampled points in C = Inv + EV (E) and V = Inv + CVaR (R).
"""

isdefined(Main, :run_stochastic_exploration_risk_pareto) || include(joinpath(@__DIR__, "..", "experiments", "MGA_Tests.jl"))
include(joinpath(@__DIR__, "sampling_toy_helpers.jl"))   # ToyCost, toy_cost, toy_eval

const GF_ALPHA = 0.1
const GF_CINV = [0.9, 0.6, 1.1, 0.4]       # $/MW investment cost of the 4 toy capacities

gf_solve(name, x, tc, p) = (tv = toy_eval(tc, x); FinalSolve(name, "", copy(x), tv[1], tv[2], p, dot(GF_CINV, x), true))
gf_C(fs) = fs.inv_cost + solve_ev(fs)
gf_V(fs) = fs.inv_cost + solve_cvar(fs, GF_ALPHA)

# Random solves; the cheapest in C / V / their mean become the Risk_Weight_1.0 / 0.0 / 0.5 sweep.
function gf_bank(seed = 11; n = 80, S = 20, d = 4)
    rng = Xoshiro(seed)
    tc = toy_cost(rng, S, 4, d)
    p = fill(1 / S, S)
    bank = [gf_solve("Random_$i", 30 .* rand(rng, d), tc, p) for i in 1:n]
    iE, iR = argmin(gf_C.(bank)), argmin(gf_V.(bank))
    iM = argmin(gf_C.(bank) .+ gf_V.(bank))
    bank[iE].name, bank[iR].name = "Risk_Weight_1.0", "Risk_Weight_0.0"
    iM in (iE, iR) || (bank[iM].name = "Risk_Weight_0.5")
    return bank, tc, p
end

# Stand-in for evaluate_gridfill_batch: solves the toy model at each pick
function gf_toy_evaluator(tc, p)
    return (picks, SPs, inputs, settings, g) -> begin
        solved = [(fs = gf_solve(pk.name, pk.x, tc, p); fs.region = "gridfill_sample"; fs) for pk in picks]
        rows = [DataFrame(Name = [pk.name]) for pk in picks]
        rows, copy(rows), copy(rows), solved
    end
end

gf_sampling(; samples = 30, batch = 8, cell = 0.1, per_cell = 2, candidates = 30) =
    Dict{String,Any}("Samples" => samples, "Batch size" => batch, "Candidates" => candidates, "Write gradients" => false,
                     "GridFill" => Dict{String,Any}("Cell size" => cell, "Delta cap" => 1.1, "Per cell" => per_cell))
const GF_INPUTS = Dict{String,Any}("Resources" => ["r1", "r2", "r3"], "Number of lines" => 1)
const GF_SETTINGS = Dict{String,Any}("Value-at-Risk percent" => GF_ALPHA, "Seed" => 7, "Convergence tolerance" => 1e-9)
const GF_FAST = merge(GRIDFILL_TUNING, (triples_per_aim = 3000, prune_blends = 300))

function gf_context(bank; g = 0.1)
    a = gridfill_anchors(bank, GF_ALPHA)
    cb = CutBank(length(bank[1].θ), length(bank[1].x))
    foreach(fs -> add_final_solve!(cb, fs), bank)
    return GridFillContext(bank, cb, a, GF_ALPHA, g)
end

@testset "GridFill sampler" begin
    bank, tc, p = gf_bank()
    a = gridfill_anchors(bank, GF_ALPHA)

    @testset "Cost coordinates and anchors" begin
        K = cost_coords(bank, a, GF_ALPHA)
        E = findfirst(fs -> fs.name == "Risk_Weight_1.0", bank)
        R = findfirst(fs -> fs.name == "Risk_Weight_0.0", bank)
        @test K[E, :] ≈ [0.0, 1.0] atol = 1e-12
        @test K[R, :] ≈ [1.0, 0.0] atol = 1e-12
        @test all(K[:, 1] .>= -1e-12)            # E is the cheapest in C
        @test all(K[:, 2] .>= -1e-12)            # R is the cheapest in V
    end

    @testset "Barycentric weights" begin
        K = [0.0 0.0; 1.0 0.0; 0.0 1.0; 2.0 2.0]
        tri = [1 2 3; 1 2 3; 1 1 2]
        t = [0.2 0.3; 0.8 0.8; 0.5 0.0]
        W = barycentric(K, tri, t)
        @test W[1, :] ≈ [0.5, 0.2, 0.3]
        @test (W[1, :]' * K[tri[1, :], :])' ≈ t[1, :]
        @test valid_weights(W, 0.95) == [true, false, false]   # outside the triangle; degenerate triple
        @test valid_weights(reshape([0.96, 0.02, 0.02], 1, 3), 0.95) == [false]
    end

    @testset "Bracket: LB <= true <= UB, CVaR(θ̄) tighter, LB monotone in cuts" begin
        rng = Xoshiro(5)
        n = length(bank)
        Θ = reduce(hcat, [fs.θ for fs in bank])
        X = reduce(vcat, [fs.x' for fs in bank])
        inv = [fs.inv_cost for fs in bank]
        tri = rand(rng, 1:n, 200, 3)
        E = randexp(rng, 200, 3); W = E ./ sum(E; dims = 2)
        Xb = blend_capacities(X, tri, W)
        cb_all = gf_context(bank).cuts
        cb_half = CutBank(length(p), 4); foreach(fs -> add_final_solve!(cb_half, fs), bank[1:2:end])
        br = blend_bracket(cb_all, Θ, inv, p, GF_ALPHA, tri, W, permutedims(Xb))
        br_half = blend_bracket(cb_half, Θ, inv, p, GF_ALPHA, tri, W, permutedims(Xb))
        for i in 1:200
            θ, _ = toy_eval(tc, Xb[i, :])
            ev, cv = dot(p, θ), dot(cvar_tail_weights(θ, p, GF_ALPHA), θ) / GF_ALPHA
            @test br.lb_ev[i] <= ev + 1e-9 && ev <= br.ub_ev[i] + 1e-9
            @test br.lb_cvar[i] <= cv + 1e-9 && cv <= br.ub_cvar[i] + 1e-9
            @test br.ub_cvar[i] <= sum(W[i, k] * solve_cvar(bank[tri[i, k]], GF_ALPHA) for k in 1:3) + 1e-9
            @test br.lb_ev[i] >= br_half.lb_ev[i] - 1e-9 && br.lb_cvar[i] >= br_half.lb_cvar[i] - 1e-9
            @test br.inv[i] ≈ dot(GF_CINV, Xb[i, :])
        end
        # Pruning to the cuts active at the solves and random blends keeps a valid (slightly looser) bound
        pr = prune_cuts(cb_all, hcat(permutedims(X), random_parent_blends(X, 300, Xoshiro(1))))
        @test n_cuts(pr) <= n_cuts(cb_all)
        br_pr = blend_bracket(pr, Θ, inv, p, GF_ALPHA, tri, W, permutedims(Xb))
        @test all(br_pr.lb_ev .<= br.lb_ev .+ 1e-9)
        @test all(cut_values(pr, permutedims(X)) .≈ cut_values(cb_all, permutedims(X)))   # exact at the reference points
    end

    @testset "Pool cuts match the same cuts built from (θ, β, x)" begin
        settings = load_settings(joinpath(@__DIR__, "..", "inputs", "Inputs_30d_1000scen_7tech_2z"))
        scale = settings["Scaling factor cost"]
        S, F = 5, 4                                 # 20 scenarios, flattened (s,f) column-major
        MP = Model()
        @variable(MP, x[1:3] >= 0); @variable(MP, x_line[1:1] >= 0)
        @variable(MP, alpha[1:S, 1:F, 1:1] >= 0); @variable(MP, u[1:S, 1:F, 1:1] >= 0); @variable(MP, ζ)
        pool = CutPool(settings)
        cb_ref = CutBank(S * F, 4)
        for (it, fs) in enumerate(bank[1:6])
            add_final_solve!(cb_ref, fs)
            θ = reshape(fs.θ ./ scale, S, F, 1); β = reshape(fs.β ./ scale, 4, S, F, 1)
            xp = fs.x
            opt = @constraint(MP, [s in 1:S, f in 1:F, k in 1:1], alpha[s,f,k] >= θ[s,f,k] + sum(β[r,s,f,k]*(x[r] - xp[r]) for r in 1:3) + β[4,s,f,k]*(x_line[1] - xp[4]), base_name = "optimality_cut_t_$(it)_")
            tail = @constraint(MP, [s in 1:S, f in 1:F, k in 1:1], u[s,f,k] >= θ[s,f,k] + sum(β[r,s,f,k]*(x[r] - xp[r]) for r in 1:3) + β[4,s,f,k]*(x_line[1] - xp[4]) - ζ, base_name = "cvar_tail_cuts_t_$(it)_")
            add_cuts!(pool, Dict(name(c) => c for c in vcat(vec(opt), vec(tail))))
        end
        # Archive half of them: both the live and the archived path are read
        to_archive = [n for n in keys(pool.live) if occursin("_t_2_", n) || occursin("_t_5_", n)]
        _archive!(pool, MP, to_archive)
        cb_pool = CutBank(S * F, 4)
        add_pool_cuts!(cb_pool, pool, MP, settings)
        @test n_cuts(cb_pool) == n_cuts(cb_ref) == 6 * S * F       # tail cuts skipped
        Y = permutedims(reduce(vcat, [fs.x' for fs in bank]))
        @test cut_values(cb_pool, Y) ≈ cut_values(cb_ref, Y) rtol = 1e-9
    end

    @testset "Cells, frontier and unreachable cells" begin
        sweep = [0.0 1.0; 0.5 0.3; 1.0 0.0]
        g, δcap = 0.1, 1.1
        cells = grid_cells(sweep, g, δcap)
        fr = sweep_frontier(sweep)
        @test fr(0.25) ≈ 0.65 && fr(2.0) == 0.0 && fr(-1.0) == 1.0
        all_cells = [(i, j) for i in 0:10 for j in 0:10]
        inside(c) = (c[1] + c[2] + 1) * g <= δcap + 1e-9 && (c[2] + 0.5) * g >= fr((c[1] + 0.5) * g)
        @test Set(cells) == Set(filter(inside, all_cells))
        @test cell_counts([0.05 0.95; 0.07 0.91; 0.55 0.31], g) == Dict((0, 9) => 2, (5, 3) => 1)

        hull = convex_hull_2d([0.0 0.0; 1.0 0.0; 0.0 1.0; 0.2 0.2])
        @test length(hull) == 3 && in_hull(hull, (0.2, 0.2)) && !in_hull(hull, (0.8, 0.8))

        ctx = gf_context(bank)
        far = (40, 40)                                             # nowhere near any solve
        r = cell_candidates(ctx, far, GF_FAST, 10, 7, 1)
        @test r.status == :outside_hull && isempty(r.tri)
    end

    @testset "Candidates: re-aimed midpoints land in the cell" begin
        ctx = gf_context(bank)
        counts = cell_counts(ctx.K, ctx.g)
        targets = [c for c in grid_cells(cost_coords([fs for fs in bank if startswith(fs.name, "Risk_Weight")], a, GF_ALPHA), 0.1, 1.1) if !haskey(counts, c)]
        @test !isempty(targets)
        with = [cell_candidates(ctx, c, GF_FAST, 1000, 7, 1) for c in targets]
        without = [cell_candidates(ctx, c, merge(GF_FAST, (reaim_iters = 0,)), 1000, 7, 1) for c in targets]
        n_with, n_without = sum(r -> size(r.tri, 1), with), sum(r -> size(r.tri, 1), without)
        @info("GridFill toy: $(length(targets)) target cells; candidates landing in their cell with re-aim $n_with, without $n_without")
        @test n_with > 0 && n_with >= n_without
        for r in with
            r.status == :ok || continue
            @test all(cell_of(r.pred[i, :], ctx.g) == r.cell for i in axes(r.pred, 1))
            @test all(r.W .>= 0) && all(abs.(sum(r.W; dims = 2) .- 1) .< 1e-9)
            br = blend_bracket(ctx.cuts, ctx.Θ, ctx.inv, ctx.p, ctx.α, r.tri, r.W, permutedims(r.X))
            @test predicted_position(ctx, br) ≈ r.pred
            @test r.X ≈ blend_capacities(ctx.X, r.tri, r.W)
        end
        # Seeded per (seed, cell, round): order of the cells doesn't matter
        sh = shuffle(Xoshiro(2), targets)
        res_a = run_candidate_stage(ctx, targets, GF_FAST, 20, 7, 1; parallel = false)
        res_b = run_candidate_stage(ctx, sh, GF_FAST, 20, 7, 1; parallel = false)
        @test [r.cell for r in res_b[sortperm([r.cell for r in res_b])]] == [r.cell for r in res_a]
        @test all(r1.tri == r2.tri && r1.W == r2.W for (r1, r2) in zip(res_a, res_b[sortperm([r.cell for r in res_b])]))
        @test cell_candidates(ctx, targets[1], GF_FAST, 20, 7, 2).tri != res_a[1].tri || res_a[1].status != :ok
    end

    @testset "Selection: quota, spread, order" begin
        ctx = gf_context(bank)
        counts = cell_counts(ctx.K, ctx.g)
        targets = [c for c in grid_cells(cost_coords([fs for fs in bank if startswith(fs.name, "Risk_Weight")], a, GF_ALPHA), 0.1, 1.1) if !haskey(counts, c)]
        cands = pool_candidates(run_candidate_stage(ctx, targets, GF_FAST, 30, 7, 1; parallel = false))
        @test cands !== nothing
        quota = Dict(c => 2 for c in targets)
        sel = select_gridfill(cands, ctx, quota, 1000)
        picked_cells = [cands.cell[i] for (i, _) in sel]
        @test all(count(==(c), picked_cells) <= 2 for c in unique(picked_cells))
        @test length(sel) == sum(min(2, count(==(c), cands.cell)) for c in unique(cands.cell))
        @test issorted([d for (_, d) in sel]; rev = true)     # max-min: each pick's distance is <= the previous one's
        @test select_gridfill(cands, ctx, quota, 5) == sel[1:5]  # the pick order is a priority list
        # No collapse to the centroid: picks are more spread out than random candidates of the same size
        Y = cands.X ./ ctx.scale'
        mind(idx) = (D = pairwise_distances(Y[idx, :], Y[idx, :]); minimum(D[i, j] for i in eachindex(idx) for j in eachindex(idx) if i < j))
        k = min(10, length(sel))
        picks = [i for (i, _) in sel[1:k]]
        rnd = mean(mind(randperm(Xoshiro(s), size(Y, 1))[1:k]) for s in 1:20)
        @test mind(picks) >= rnd
        @test length(select_gridfill(cands, ctx, quota, 1000; score = bracket_novelty_score)) == length(sel)
    end

    @testset "Distributed candidate stage matches serial, context cleared" begin
        ctx = gf_context(bank)
        cells = [(i, j) for i in 0:9 for j in 0:9 if i + j <= 10][1:12]
        serial = run_candidate_stage(ctx, cells, GF_FAST, 10, 7, 1; parallel = false)
        added = nworkers() < 2 || workers() == [myid()] ? addprocs(2) : Int[]
        try
            configure_parallel_workers!(Dict{String,Any}("Parallel flag" => true, "Workers" => nworkers()))
            par = run_candidate_stage(ctx, cells, GF_FAST, 10, 7, 1; parallel = true)
            @test [r.cell for r in par] == [r.cell for r in serial]
            @test all(r1.status == r2.status && r1.tri == r2.tri && r1.W ≈ r2.W for (r1, r2) in zip(par, serial))
            @test all(remotecall_fetch(() -> GRIDFILL_WORKER_CTX[] === nothing, pid) for pid in workers())
        finally
            isempty(added) || rmprocs(added)
        end
    end

    @testset "Full loop with feedback (toy evaluator)" begin
        bank1, tc1, p1 = gf_bank()
        n0 = length(bank1)
        sampling = gf_sampling(samples = 30, batch = 8)
        ckpts = Int[]
        c, s, e, labels = run_gridfill_sampling!(bank1, nothing, nothing, nothing, GF_INPUTS, GF_SETTINGS, sampling;
                                                 checkpoint = (c, s, e, l) -> push!(ckpts, length(l)), tuning = GF_FAST,
                                                 evaluate = gf_toy_evaluator(tc1, p1), parallel = false)
        @test 0 < length(labels) <= 30 && length(bank1) == n0 + length(labels)
        @test labels == ["GridFill_$k" for k in eachindex(labels)]
        @test issorted(ckpts; lt = <=) && ckpts[end] == length(labels) && all(diff([0; ckpts]) .<= 8)   # one checkpoint per batch
        @test all(fs.region == "gridfill_sample" for fs in bank1[n0+1:end])
        # Every new sample sits within the per-cell quota of the cells it was aimed at, and labels
        # come from solves (bracket holds for each sample against its own solve)
        gf = bank1[n0+1:end]
        a1 = gridfill_anchors(bank1, GF_ALPHA)
        K0 = cost_coords(bank1[1:n0], a1, GF_ALPHA)
        Kn = cost_coords(gf, a1, GF_ALPHA)
        occupied0 = keys(cell_counts(K0, 0.1))
        landed = [cell_of(Kn[i, :], 0.1) for i in axes(Kn, 1)]
        @info("GridFill toy loop: $(length(labels)) samples, $(count(c -> !(c in occupied0), landed)) landed in initially empty cells")
        @test count(c -> !(c in occupied0), landed) >= length(labels) ÷ 2
        # Reproducible
        bank2, _, _ = gf_bank()
        _, _, _, labels2 = run_gridfill_sampling!(bank2, nothing, nothing, nothing, GF_INPUTS, GF_SETTINGS, sampling;
                                                  tuning = GF_FAST, evaluate = gf_toy_evaluator(tc1, p1), parallel = false)
        @test labels2 == labels && all(f1.x == f2.x for (f1, f2) in zip(bank1[n0+1:end], bank2[n0+1:end]))
    end

    @testset "Per cell \"auto\": layers until Samples is spent" begin
        # Fixed quota of 1 stops once every reachable target cell holds a sample
        bank1, tc, p = gf_bank()
        n0 = length(bank1)
        _, _, _, l1 = run_gridfill_sampling!(bank1, nothing, nothing, nothing, GF_INPUTS, GF_SETTINGS, gf_sampling(samples = 500, per_cell = 1);
                                             tuning = GF_FAST, evaluate = gf_toy_evaluator(tc, p), parallel = false)
        @test 0 < length(l1) < 500
        # "auto" runs the same first layer, then keeps going to exactly Samples
        n_total = length(l1) + 12
        bank2, _, _ = gf_bank()
        _, _, _, l2 = run_gridfill_sampling!(bank2, nothing, nothing, nothing, GF_INPUTS, GF_SETTINGS, gf_sampling(samples = n_total, per_cell = "auto");
                                             tuning = GF_FAST, evaluate = gf_toy_evaluator(tc, p), parallel = false)
        @test length(l2) == n_total
        @test all(f1.x == f2.x for (f1, f2) in zip(bank1[n0+1:end], bank2[n0+1:n0+length(l1)]))
    end
end
