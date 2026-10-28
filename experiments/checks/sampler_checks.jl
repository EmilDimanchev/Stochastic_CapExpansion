#=========

Hit-and-run sampler test: sample_interior_distributed() previously only supported "CVT"
(k-means) and "Random" (uniform-box + convex-hull projection via an LP/QP solve) sample
methods. This test exercises the new "HitAndRun" method, which walks the same
convex-hull-of-points polytope (expressed in weight-space as l>=0, sum(l)=1, l<=0.95,
exactly as the existing MGCA projection does) via a hit-and-run MCMC chain instead,
needing no solver at all.

It also exercises a real bug fix bundled into the same change: the previous `if
nworkers() > 0` branch condition in sample_interior_distributed was always true in
practice (Distributed.nworkers() returns 1, not 0, with no extra processes), so the
serial `else` branch was dead code containing two further bugs (an undefined `pids`
reference, and a sampling loop that never wrote into `samples[idx]`). The condition is
now `nprocs() > 1`, so this test deliberately runs the sampler once before adding any
workers (serial path) and once after (distributed path).

=========#

function test_hit_and_run_sampler(test_index = "hit_and_run")
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")
    settings = load_settings(inputs_folder)
    settings["Sampling settings"]["Delaunay"]["Interior sample method"] = "HitAndRun"
    settings["Workers"] = 2

    points = [1.0 0.0 0.0 0.0;
              0.0 1.0 0.0 0.0;
              0.0 0.0 1.0 0.0;
              0.0 0.0 0.0 1.0;
              0.5 0.5 0.0 0.0;
              0.0 0.0 0.5 0.5]
    rows, cols = size(points)
    lo = vec(minimum(points, dims=1))
    hi = vec(maximum(points, dims=1))
    n_samples = 40

    # --- l-space checks on the chain in isolation ---
    rng = MersenneTwister(42)
    ls = _run_hit_and_run_chain(rows, n_samples; rng = rng, burn_in = 200, thin = 20)
    for l in ls
        @assert all(l .>= -1e-9) && all(l .<= 0.95 + 1e-9) "l out of [0, 0.95] bounds"
        @assert abs(sum(l) - 1.0) < 1e-6 "l does not sum to 1"
    end
    @assert norm(ls[end] - ls[1]) > 1e-3 "chain did not move in l-space - suspect a self-loop or burn-in/thin no-op bug"
    @info("[$(test_index)] l-space chain: all $(n_samples) samples satisfy l>=0, l<=0.95, sum(l)=1; chain moved (||l_end - l_1|| = $(round(norm(ls[end]-ls[1]); digits=4))).")

    # --- serial path (nprocs()==1 -> the fixed `else` branch) ---
    @assert nprocs() == 1 "expected a clean single-process session for the serial-path assertions"
    samples_serial = sample_interior_distributed(points, n_samples, settings)
    @assert length(samples_serial) == n_samples
    for x in samples_serial
        @assert all(x .>= lo .- 1e-6) && all(x .<= hi .+ 1e-6) "capacity-space sample outside bounding box of points"
    end
    spread_serial = maximum(norm(samples_serial[i] - samples_serial[1]) for i in 2:n_samples)
    @assert spread_serial > 1e-3 "serial samples show no spread - suspect samples collapsed to the start point"
    @info("[$(test_index)] Serial path: $(n_samples) samples all within bounding box of points; spread=$(round(spread_serial; digits=4)).")

    # --- distributed path (nprocs()>1 -> the distributed branch, independent chain per worker) ---
    configure_parallel_workers!(settings)
    @assert nprocs() > 1
    samples_dist = sample_interior_distributed(points, n_samples, settings)
    @assert length(samples_dist) == n_samples
    for x in samples_dist
        @assert all(x .>= lo .- 1e-6) && all(x .<= hi .+ 1e-6) "capacity-space sample outside bounding box of points"
    end
    spread_dist = maximum(norm(samples_dist[i] - samples_dist[1]) for i in 2:n_samples)
    @assert spread_dist > 1e-3 "distributed samples show no spread"
    # n_samples=40 over Workers=2 divides evenly (per_worker=20), so index 1 and 21 are
    # each the first sample drawn by a different worker's chain - distinct seeds (Seed +
    # myid()) should make these different chains, not identical draws.
    @assert samples_dist[1] != samples_dist[21] "first samples from worker 1 and worker 2 are identical - suspect identical RNG seeds across workers"
    @info("[$(test_index)] Distributed path ($(nworkers()) workers): $(n_samples) samples all within bounding box of points; spread=$(round(spread_dist; digits=4)); per-worker chains are distinct.")

    rmprocs(workers())
    @info("[$(test_index)] PASSED: hit-and-run sampler correct in l-space and capacity-space, on both serial and distributed paths.")

    return (spread_serial = spread_serial, spread_dist = spread_dist)
end

function test_sample_interior_delaunay(test_index = "sample_interior_delaunay")
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")
    settings = load_settings(inputs_folder)

    # capacity-space points (same shape/style as test_hit_and_run_sampler)
    points = [1.0 0.0 0.0 0.0;
              0.0 1.0 0.0 0.0;
              0.0 0.0 1.0 0.0;
              0.0 0.0 0.0 1.0;
              0.5 0.5 0.0 0.0;
              0.0 0.0 0.5 0.5]
    # matching cost-space points (Investment_Cost+EV_OpCost, Investment_Cost+CVaR_OpCost), row i
    # aligned with row i of `points` above - a convex pentagon plus one interior point so the
    # triangulation has several simplices of clearly different areas to rank.
    cost_points = [0.0 0.0;
                   4.0 0.0;
                   4.0 3.0;
                   2.0 5.0;
                   0.0 3.0;
                   2.0 1.5]
    lo = vec(minimum(points, dims=1))
    hi = vec(maximum(points, dims=1))
    n_samples = 20

    @assert nprocs() == 1 "expected a clean single-process session - this sampler must never dispatch to workers"

    # --- unit check on the triangulation+ranking helper in isolation ---
    unit_square = [0.0 0.0; 1.0 0.0; 0.0 1.0; 1.0 1.0]
    ranked, elapsed = _build_cost_hull_simplices(unit_square)
    @assert length(ranked) == 2 "unit square should triangulate into exactly 2 simplices"
    @assert all(abs(area - 0.5) < 1e-9 for (_, area) in ranked) "unit square simplices should each have area 0.5, got areas=$([a for (_,a) in ranked])"
    @assert elapsed >= 0.0
    @info("[$(test_index)] _build_cost_hull_simplices unit check: 2 simplices, areas=$([round(a; digits=4) for (_,a) in ranked]).")

    # --- CycleSampling mode: single triangulation, cycles through top-fraction pool with fresh
    # random Dirichlet draws each pass, so re-running with the same Seed must reproduce exactly ---
    settings["Sampling settings"]["Delaunay"]["Repeat mode"] = "CycleSampling"
    settings["Seed"] = 4321
    samples_cycle_a = sample_interior_distributed(points, cost_points, n_samples, settings)
    samples_cycle_b = sample_interior_distributed(points, cost_points, n_samples, settings)
    @assert length(samples_cycle_a) == n_samples
    @assert all(samples_cycle_a[i] == samples_cycle_b[i] for i in 1:n_samples) "CycleSampling is not reproducible under a fixed Seed"
    for x in samples_cycle_a
        @assert all(x .>= lo .- 1e-6) && all(x .<= hi .+ 1e-6) "capacity-space sample outside bounding box of points"
    end
    spread_cycle = maximum(norm(samples_cycle_a[i] - samples_cycle_a[1]) for i in 2:n_samples)
    @assert spread_cycle > 1e-6 "CycleSampling samples show no spread - suspect Dirichlet draws collapsed to the same point"
    @info("[$(test_index)] CycleSampling: $(n_samples) samples within bounding box of points, reproducible under fixed Seed, spread=$(round(spread_cycle; digits=4)).")

    # --- ProgressiveRefinement mode: fixed centroids only, so it is fully deterministic
    # regardless of Seed/rng ---
    settings["Sampling settings"]["Delaunay"]["Repeat mode"] = "ProgressiveRefinement"
    samples_prog_a = sample_interior_distributed(points, cost_points, n_samples, settings)
    settings["Seed"] = 9999
    samples_prog_b = sample_interior_distributed(points, cost_points, n_samples, settings)
    @assert length(samples_prog_a) == n_samples
    @assert all(samples_prog_a[i] == samples_prog_b[i] for i in 1:n_samples) "ProgressiveRefinement should be deterministic (centroids only) regardless of Seed"
    for x in samples_prog_a
        @assert all(x .>= lo .- 1e-6) && all(x .<= hi .+ 1e-6) "capacity-space sample outside bounding box of points"
    end
    spread_prog = maximum(norm(samples_prog_a[i] - samples_prog_a[1]) for i in 2:n_samples)
    @assert spread_prog > 1e-6 "ProgressiveRefinement samples show no spread"
    @info("[$(test_index)] ProgressiveRefinement: $(n_samples) samples within bounding box of points, deterministic across Seeds, spread=$(round(spread_prog; digits=4)).")

    @assert nprocs() == 1 "sampler must not have spawned workers"
    @info("[$(test_index)] PASSED: Delaunay interior sampling correct for both CycleSampling and ProgressiveRefinement modes.")

    return (spread_cycle = spread_cycle, spread_prog = spread_prog)
end
