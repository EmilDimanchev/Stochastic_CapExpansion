using Test
using Distributed
using LinearAlgebra, Random, Statistics, DataFrames, JuMP

"""
Solver-free unit tests for the GridFill sampler (src/model/sampling/gridfill_*.jl).

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

gf_solve(name, x, tc, p) = (tv = toy_eval(tc, x); FinalSolve(name, copy(x), tv[1], tv[2], p, dot(GF_CINV, x), true))
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
        solved = [gf_solve(pk.name, pk.x, tc, p) for pk in picks]
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

# Tuning whose triple count can be switched off mid-run, to force rounds with no candidates
mutable struct GFSwitchableTuning
    aims_per_cell::Int
    triples_per_aim::Int
    max_weight::Float64
    reaim_pool::Int
    reaim_iters::Int
    prune_blends::Int
end
GFSwitchableTuning(t::NamedTuple) = GFSwitchableTuning(t.aims_per_cell, t.triples_per_aim, t.max_weight, t.reaim_pool, t.reaim_iters, t.prune_blends)

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
        settings = load_settings(joinpath(@__DIR__, "..", "inputs", "Inputs_2z_7tech"))
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
        @test all(startswith(fs.name, "GridFill_") for fs in bank1[n0+1:end])
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

    @testset "Aim: LB + f·(UB − LB), bracket columns unchanged" begin
        ctx = gf_context(bank)
        @test ctx.aim_frac == (ev = 0.10, cvar = 0.13)
        ctx_mid = GridFillContext(bank, ctx.cuts, a, GF_ALPHA, 0.1; aim_frac = 0.5)
        tri = [1 2 3; 4 5 6]; W = [0.5 0.3 0.2; 0.2 0.2 0.6]
        X = blend_capacities(ctx.X, tri, W)
        br = blend_bracket(ctx.cuts, ctx.Θ, ctx.inv, ctx.p, ctx.α, tri, W, permutedims(X))
        @test br.mid_ev ≈ (br.ub_ev .+ br.lb_ev) ./ 2
        @test predicted_position(ctx_mid, br) ≈ hcat(c_hat(a, br.inv .+ br.mid_ev), v_hat(a, br.inv .+ br.mid_cvar))
        @test predicted_position(ctx, br) ≈ hcat(c_hat(a, br.inv .+ br.lb_ev .+ 0.10 .* (br.ub_ev .- br.lb_ev)),
                                                 v_hat(a, br.inv .+ br.lb_cvar .+ 0.13 .* (br.ub_cvar .- br.lb_cvar)))
        # Each metric uses its own fraction; one number applies to both
        ctx_split = GridFillContext(bank, ctx.cuts, a, GF_ALPHA, 0.1; aim_frac = [0.5, 0.0])
        @test predicted_position(ctx_split, br) ≈ hcat(c_hat(a, br.inv .+ br.mid_ev), v_hat(a, br.inv .+ br.lb_cvar))
        @test aim_fractions(0.3) == (ev = 0.3, cvar = 0.3) && aim_fractions([0.1, 0.2]) == (ev = 0.1, cvar = 0.2)
        @test_throws ErrorException aim_fractions([0.1, 0.2, 0.3])
        @test_throws ErrorException aim_fractions(1.5)
    end

    @testset "Aim margin bounds the aim box" begin
        ctx = gf_context(bank)
        counts = cell_counts(ctx.K, ctx.g)
        targets = [c for c in grid_cells(cost_coords([fs for fs in bank if startswith(fs.name, "Risk_Weight")], a, GF_ALPHA), 0.1, 1.1) if !haskey(counts, c)]
        for margin in (0.1, 0.2, 0.5)
            m_ctx = GridFillContext(bank, ctx.cuts, a, GF_ALPHA, 0.1; aim_margin = margin)
            results = [r for r in (cell_candidates(m_ctx, c, GF_FAST, 1000, 7, 1) for c in targets) if r.status == :ok]
            @test !isempty(results)
            fracs = reduce(vcat, [r.aim ./ 0.1 .- [r.cell[1] r.cell[2]] for r in results])
            @test minimum(fracs) >= 0.1 - 1e-9 && maximum(fracs) <= 1 - margin + 1e-9
        end
        @test_throws ErrorException GridFillContext(bank, ctx.cuts, a, GF_ALPHA, 0.1; aim_margin = 0.9)
    end

    @testset "Hull depth, nearest edge and unreached labels" begin
        hull = convex_hull_2d([0.0 0.0; 4.0 0.0; 4.0 4.0; 0.0 4.0])   # square, counter-clockwise
        e = hull_nearest_edge(hull, (2.0, 0.5))
        @test e.depth ≈ 0.5 && e.lower
        e = hull_nearest_edge(hull, (2.0, 3.5))
        @test e.depth ≈ 0.5 && !e.lower
        @test hull_nearest_edge(hull, (2.0, -1.0)).depth ≈ -1.0
        @test hull_nearest_edge(hull[1:2], (0.0, 0.0)).depth == -Inf
        # g = 1: centres at 0.5 cells above the lower edge, 1.5 above it, and 0.5 below the upper edge
        @test near_lower_frontier([(2, 0), (2, 1), (2, 3)], hull, 1.0, 1.0) == [(2, 0)]
        @test isempty(near_lower_frontier([(2, 0)], hull, 1.0, 0.0))
        st = GridFillCellState()
        for c in [(2, 0), (1, 1), (2, 3), (1, 0), (3, 3)]
            st.exhausted[c] = :missed
        end
        st.exhausted[(0, 3)] = :outside_hull
        push!(st.picked, (1, 0))                  # a sample was aimed here: a real gap, not unreachable
        counts = Dict((3, 3) => 1)                # holds a sample: not unreached
        lab = label_unreached_cells(st, counts, hull, 1.0, 1.0)
        @test lab.frontier == [(2, 0)] && lab.out_of_reach == [(2, 3)]     # (1, 1) is 1.5 deep on the lower side
    end

    @testset "Cell state: attempts, misses, exhaustion" begin
        st = GridFillCellState()
        @test record_attempt!(st, (1, 1)) == 1 && record_attempt!(st, (1, 1)) == 2 && record_attempt!(st, (2, 2)) == 1
        @test !record_miss!(st, (1, 1), 3) && !record_miss!(st, (1, 1), 3) && !is_exhausted(st, (1, 1))
        @test !record_miss!(st, (1, 1), 3; open = false) && !is_exhausted(st, (1, 1))   # at quota: not exhausted yet
        @test record_miss!(st, (1, 1), 3) && is_exhausted(st, (1, 1))
        @test !record_miss!(st, (1, 1), 3)                                              # only newly exhausted counts
        # Outside the hull: exhausted on the first occurrence, with its own reason
        far = cell_candidates(gf_context(bank), (40, 40), GF_FAST, 10, 7, 1)
        mark_outside_hull!(st, [far])
        @test exhausted_cells(st, :outside_hull) == [(40, 40)] && exhausted_cells(st, :missed) == [(1, 1)]
        # Solved misses stay cumulative: hits in between don't reset the count
        st2 = GridFillCellState()
        record_miss!(st2, (1, 1), 3); record_attempt!(st2, (1, 1)); record_miss!(st2, (1, 1), 3)
        @test record_miss!(st2, (1, 1), 3) && is_exhausted(st2, (1, 1))
    end

    @testset "No candidate: empty cells take a miss, cells with samples are parked" begin
        st = GridFillCellState()
        # Empty cell: three rounds exhaust it, nothing parked
        @test !record_no_candidate!(st, (1, 1), 0, 3, 3) && !record_no_candidate!(st, (1, 1), 0, 3, 3)
        @test record_no_candidate!(st, (1, 1), 0, 3, 3) && is_exhausted(st, (1, 1)) && isempty(st.parked)
        # Holding samples: no misses charged, parked at the cap, never exhausted
        for _ in 1:3
            @test !record_no_candidate!(st, (2, 2), 4, 3, 3)
        end
        @test (2, 2) in st.parked && !is_exhausted(st, (2, 2)) && get(st.misses, (2, 2), 0) == 0
        unpark_all!(st)
        @test isempty(st.parked) && isempty(st.dry) && is_exhausted(st, (1, 1))
    end

    # Rows carry the GF_* columns, solves optionally reported as failed (every pick then misses)
    gf_column_evaluator(tc, p; ok = true) = (picks, SPs, inputs, settings, g) -> begin
        solved = [(fs = gf_solve(pk.name, pk.x, tc, p); fs.ok = ok; fs) for pk in picks]
        rows = [add_gridfill_columns!(DataFrame(Name = [pk.name]), pk, g) for pk in picks]
        rows, rows, rows, solved
    end

    @testset "Retry cap: cells that always miss are dropped and the loop stops" begin
        bank1, tc, p = gf_bank()
        _, s, _, labels = run_gridfill_sampling!(bank1, nothing, nothing, nothing, GF_INPUTS, GF_SETTINGS, gf_sampling(samples = 5000, per_cell = "auto");
                                                 tuning = GF_FAST, evaluate = gf_column_evaluator(tc, p; ok = false), parallel = false)
        @test 0 < length(labels) < 5000          # auto only stops early once every target is exhausted
        df = reduce(vcat, s)
        @test maximum(df.GF_Attempt) <= 3 && all(==(1), df.GF_Layer)
        # A lower cap gives up sooner
        bank2, _, _ = gf_bank()
        sampling = gf_sampling(samples = 5000, per_cell = "auto"); sampling["GridFill"]["Max misses"] = 1
        _, s2, _, l2 = run_gridfill_sampling!(bank2, nothing, nothing, nothing, GF_INPUTS, GF_SETTINGS, sampling;
                                              tuning = GF_FAST, evaluate = gf_column_evaluator(tc, p; ok = false), parallel = false)
        @test 0 < length(l2) < length(labels) && maximum(reduce(vcat, s2).GF_Attempt) == 1
    end

    @testset "Termination: stops when a whole layer has no candidates" begin
        # Candidates exist for the first batch only; afterwards every open cell stays dry
        bank1, tc, p = gf_bank()
        tuning = GFSwitchableTuning(GF_FAST)
        base = gf_toy_evaluator(tc, p)
        evaluate = (picks, SPs, inputs, settings, g) -> (tuning.triples_per_aim = 0; base(picks, SPs, inputs, settings, g))
        _, _, _, labels = run_gridfill_sampling!(bank1, nothing, nothing, nothing, GF_INPUTS, GF_SETTINGS, gf_sampling(samples = 5000, per_cell = "auto");
                                                 tuning = tuning, evaluate = evaluate, parallel = false)
        @test length(labels) == 8               # one batch, then no samples, and the loop ended
    end

    @testset "Layers advance past exhausted cells; GF_Layer and GF_Attempt" begin
        bank1, tc, p = gf_bank()
        _, _, _, l1 = run_gridfill_sampling!(bank1, nothing, nothing, nothing, GF_INPUTS, GF_SETTINGS, gf_sampling(samples = 500, per_cell = 1);
                                             tuning = GF_FAST, evaluate = gf_toy_evaluator(tc, p), parallel = false)
        n_total = length(l1) + 12
        bank2, _, _ = gf_bank()
        _, s, _, l2 = run_gridfill_sampling!(bank2, nothing, nothing, nothing, GF_INPUTS, GF_SETTINGS, gf_sampling(samples = n_total, per_cell = "auto");
                                             tuning = GF_FAST, evaluate = gf_column_evaluator(tc, p), parallel = false)
        df = reduce(vcat, s)
        @test length(l2) == n_total
        @test issorted(df.GF_Layer) && maximum(df.GF_Layer) >= 2
        @test count(==(1), df.GF_Layer) == length(l1)
        for sub in groupby(df, [:GF_Cell_c, :GF_Cell_v])
            @test all(diff(sub.GF_Attempt) .> 0)   # attempts count up per target cell
        end
    end
    @testset "nearest_distance in row chunks matches the full distance matrix" begin
        full_nearest(Y, R, sc) = size(R, 1) == 0 ? fill(Inf, size(Y, 1)) : vec(minimum(pairwise_distances(Y ./ sc', R ./ sc'); dims = 2))
        rng = Xoshiro(5)
        d = 14
        sc = 1 .+ 100 .* rand(rng, d)
        R = 1000 .* rand(rng, 60, d)
        for (m, chunk) in [(50, 1), (50, 7), (5, 1024), (2048, 1024), (2500, 1024), (0, 16)]
            Y = 1000 .* rand(rng, m, d)
            got, want = nearest_distance(Y, R, sc; chunk = chunk), full_nearest(Y, R, sc)
            @test length(got) == m
            if got != want
                @info("nearest_distance chunk = $chunk, m = $m: not bit-identical, max difference $(maximum(abs.(got .- want)))")
            end
            @test got ≈ want rtol = 1e-12
        end
        @test nearest_distance(rand(rng, 10, d), zeros(0, d), sc) == fill(Inf, 10)
        @test_throws ErrorException nearest_distance(rand(rng, 3, d), R, sc; chunk = 0)

        # Peak memory is measured in fresh processes by test_nearest_distance_memory(): @allocated
        # counts every block, so it is the same for both versions and can't show the peak.
    end

    @testset "Valid-blend count per cell (sizes the candidate-stage temporaries)" begin
        ctx = gf_context(bank)
        counts = cell_counts(ctx.K, ctx.g)
        targets = [c for c in grid_cells(cost_coords([fs for fs in bank if startswith(fs.name, "Risk_Weight")], a, GF_ALPHA), 0.1, 1.1) if !haskey(counts, c)]
        res = [cell_candidates(ctx, c, GRIDFILL_TUNING, 30, 7, 1) for c in targets]
        drawn = GRIDFILL_TUNING.aims_per_cell * GRIDFILL_TUNING.triples_per_aim
        @test all(0 <= r.n_valid <= drawn for r in res)
        @test all(r.n_valid > 0 for r in res if r.status == :ok)
        @test all(r.n_valid == 0 for r in res if r.status in (:no_blend, :outside_hull))
        fr = [r.n_valid / drawn for r in res if r.n_valid > 0]
        @info("GridFill toy ($(length(bank)) plans): valid blends per cell, as a share of the $drawn drawn: median $(round(median(fr); digits = 3)), max $(round(maximum(fr); digits = 3))")
    end
end
