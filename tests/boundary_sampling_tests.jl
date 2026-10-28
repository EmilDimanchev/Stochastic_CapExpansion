using Test
using Distributed
using LinearAlgebra, Random, Statistics

"""
Solver-free unit tests for the boundary interior sampler (src/model/boundary_sampling).

Usage:
  julia --project=. tests/boundary_sampling_tests.jl
or, from a REPL that already ran include("experiments/MGA_Tests.jl"):
  test_boundary_sampling_units()

Synthetic data: each scenario's op cost is a convex piecewise-linear function of capacity,
θ_s(x) = max_l (a_sl + c_slᵀ x), so its exact gradient (the active piece's c) is known and
the bounds can be checked against the true value of any blend.
"""

isdefined(Main, :run_stochastic_exploration_risk_pareto) || include(joinpath(@__DIR__, "..", "experiments", "MGA_Tests.jl"))

include(joinpath(@__DIR__, "sampling_toy_helpers.jl"))   # ToyCost, toy_cost, toy_eval, toy_bank

toy_cfg() = boundary_cfg(resolve_settings(load_settings(joinpath(@__DIR__, "..", "inputs", "Inputs_30d_1000scen_7tech_2z")), "Sampling settings", Dict("Candidates" => 10)))
unit_metric(d) = CapacityMetric("lengthscale_file", ["v$i" for i in 1:d], ones(d), fill(5.0, d), zeros(d, 0), Float64[], d, NaN, NaN, 0, NaN, NaN)
X_of(bank) = reduce(vcat, [fs.x' for fs in bank])

# compute_cvar exactly as it was before the tail-weight refactor
function old_compute_cvar(SP_obj, P, P_f, P_k, alpha)
    sp = vec(SP_obj)
    pr = vec([P[s]*P_f[f]*P_k[k] for s in 1:size(SP_obj,1), f in 1:size(SP_obj,2), k in 1:size(SP_obj,3)])
    pr = pr ./ sum(pr)
    idx = sortperm(sp; rev = true)
    sps, prs = sp[idx], pr[idx]
    rem, acc = alpha, 0.0
    for i in eachindex(sps)
        rem <= 1e-14 && break
        w = min(prs[i], rem)
        acc += sps[i] * w
        rem -= w
    end
    return acc / alpha
end

function main()
    rng = Xoshiro(42)

    @testset "Boundary sampler" begin

        @testset "CVaR tail weights match compute_cvar (and compute_cvar is unchanged)" begin
            for α in (0.05, 0.033, 0.5)
                θ = 1e6 .* rand(rng, 100)
                P_s, P_f, P_k = fill(0.1, 10), fill(0.1, 10), [1.0]
                θ3 = reshape(θ, 10, 10, 1)
                w = cvar_tail_weights(θ, vec([P_s[s]*P_f[f] for s in 1:10, f in 1:10]), α)
                @test all(w .>= 0) && sum(w) ≈ α
                @test dot(w, θ) / α ≈ compute_cvar(θ3, P_s, P_f, P_k, α)
                @test compute_cvar(θ3, P_s, P_f, P_k, α) == old_compute_cvar(θ3, P_s, P_f, P_k, α)
            end
            # unequal probabilities, fractional boundary scenario
            θ, p = [5.0, 1.0, 9.0, 3.0], [0.1, 0.4, 0.05, 0.45]
            w = cvar_tail_weights(θ, p, 0.1)
            @test w ≈ [0.05, 0.0, 0.05, 0.0]
        end

        @testset "Region labels" begin
            bank, _ = toy_bank(rng)
            region(n) = bank[findfirst(fs -> fs.name == n, bank)].region
            @test region("Risk_Weight_0.5") == "anchor"
            @test region("Random_3_Budget_1") == "outer_shell"
            @test region("Random_3_Budget_2") == "interior"
            @test region("Random_3_Budget_1_System_Expected_Cap_2") == "ev_wing"
            @test region("Random_3_Budget_1_CVaR_Cap_1") == "cvar_wing"
            # no Budget tag (no tightening) = loosest budget
            b2 = [FinalSolve("Random_$k", "", zeros(2), zeros(1), zeros(2, 1), [1.0], 0.0, true) for k in 1:3]
            @test all(fs.region == "outer_shell" for fs in assign_regions!(b2))
        end

        @testset "Thin SVD matches explicit M" begin
            for (n, d) in ((20, 8), (60, 30))
                G = randn(rng, n, d)
                λe, We = principal_directions(G; explicit = true)
                λs, Ws = principal_directions(G; explicit = false)
                k = min(n, d)
                @test λe[1:k] ≈ λs[1:k] rtol = 1e-8
                @test abs.(diag(We[:, 1:k]' * Ws[:, 1:k])) ≈ ones(k) rtol = 1e-6   # same directions up to sign
            end
            @test choose_rank([5.0, 3.0, 1.0, 1.0], 0.8, 50) == 2
            @test choose_rank([5.0, 3.0, 1.0, 1.0], 0.95, 50) == 4
            @test choose_rank([5.0, 3.0, 1.0, 1.0], 0.95, 3) == 3
        end

        @testset "Active subspace build" begin
            bank, _ = toy_bank(rng)
            S = length(bank[1].θ)
            # anchors with positive gaps: E (EV-optimal) has low EV / high CVaR, R the reverse
            E = bank[findfirst(fs -> fs.name == "Risk_Weight_1.0", bank)]
            R = bank[findfirst(fs -> fs.name == "Risk_Weight_0.0", bank)]
            E.θ = vcat(zeros(S - 1), 200.0); E.inv_cost = 0.0
            R.θ = fill(30.0, S); R.inv_cost = 0.0
            cfg = toy_cfg()
            m = build_active_subspace(bank, ["v$i" for i in 1:5], cfg, 0.05)
            @test m.ΔC ≈ 20.0 && m.ΔV ≈ 170.0
            @test 1 <= m.r <= 5 && issorted(m.λ_all; rev = true)
            # explicit vs SVD path inside the build give the same metric distances
            m_svd = build_active_subspace(bank, ["v$i" for i in 1:5], merge(cfg, Dict("Explicit M max dim" => 0)), 0.05)
            Z1, Z2 = metric_coords(m, X_of(bank)), metric_coords(m_svd, X_of(bank))
            @test pairwise_distances(Z1, Z1) ≈ pairwise_distances(Z2, Z2) rtol = 1e-6
            # with r = d, ||Δz||² = Δx̃ᵀ M Δx̃
            m_full = build_active_subspace(bank, ["v$i" for i in 1:5], merge(cfg, Dict("Energy share" => 1.0)), 0.05)
            M = m_full.W * Diagonal(m_full.λ_all[1:m_full.r]) * m_full.W'
            Δ = (bank[4].x .- bank[5].x) ./ m_full.cap_max
            @test norm(metric_coords(m_full, bank[4].x) .- metric_coords(m_full, bank[5].x))^2 ≈ dot(Δ, M * Δ) rtol = 1e-6
        end

        @testset "Candidates: weights, blends, linear z" begin
            bank, _ = toy_bank(rng)
            cfg = toy_cfg()
            Z = metric_coords(unit_metric(5), X_of(bank))
            nbr = knn_indices(Z, cfg["K neighbours"])
            counts = quota_counts(cfg["Quotas"], 40)
            cands = generate_candidates(Xoshiro(1), [fs.region for fs in bank], nbr, counts, cfg)
            lo, hi = cfg["Base weight range"]
            for c in cands
                @test all(c.weights .>= 0) && sum(c.weights) ≈ 1
                @test lo <= c.weights[1] <= hi
                @test 2 <= length(c.parents) <= 4 && allunique(c.parents)
                @test all(p in nbr[c.parents[1]] for p in c.parents[2:end])
                @test bank[c.parents[1]].region == c.region
                x = blend_rows(X_of(bank), c)
                @test x ≈ sum(c.weights[j] .* bank[c.parents[j]].x for j in eachindex(c.parents))
                @test blend_rows(Z, c) ≈ metric_coords(unit_metric(5), x)
            end
        end

        @testset "Quotas, reproducibility, z linearity and bounds of picks" begin
            bank, tc = toy_bank(rng)
            cfg = toy_cfg()
            @test quota_counts(cfg["Quotas"], 300) == [("outer_shell", 150), ("ev_wing", 45), ("cvar_wing", 45), ("interior", 60)]
            n = 40                          # 20 / 6 / 6 / 8 picks
            m = unit_metric(5)
            Z = metric_coords(m, X_of(bank))
            picks = select_boundary_samples(bank, Z, cfg, Xoshiro(0), 0.05; n_new = n)
            @test [count(p -> p.region == r, picks) for r in BOUNDARY_REGION_ORDER] == [20, 6, 6, 8]
            @test all(p.dist_to_nearest_solve >= 0 for p in picks)
            again = select_boundary_samples(bank, Z, cfg, Xoshiro(0), 0.05; n_new = n)
            @test [p.x for p in picks] == [p.x for p in again]
            other = select_boundary_samples(bank, Z, cfg, Xoshiro(1), 0.05; n_new = n)
            @test [p.x for p in picks] != [p.x for p in other]
            byname = Dict(fs.name => fs for fs in bank)
            for pk in picks
                parents = [byname[n] for n in pk.parent_names]
                @test pk.x ≈ sum(pk.weights[j] .* parents[j].x for j in eachindex(parents))
                @test pk.z ≈ metric_coords(m, pk.x)
                @test pk.blend_inv ≈ sum(pk.x)                            # toy investment = sum(x), linear
                θ, _ = toy_eval(tc, pk.x)                                 # true op cost of the blend
                p = parents[1].p
                ev, cv = dot(p, θ), dot(cvar_tail_weights(θ, p, 0.05), θ) / 0.05
                @test pk.lb_ev <= ev + 1e-9 && ev <= pk.blend_ev + 1e-9
                @test pk.lb_cvar <= cv + 1e-9 && cv <= pk.blend_cvar + 1e-9
            end
            # waves: solved interior points enter the reference set, so new picks keep away from them
            Zs = reduce(vcat, [p.z' for p in picks])
            wave2 = select_boundary_samples(bank, Z, cfg, Xoshiro(7), 0.05; n_new = n, extra_ref_Z = Zs, first_index = 41, wave = 2)
            @test wave2[1].name == "Interior_41" && all(p.wave == 2 for p in wave2)
            @test all(minimum(pairwise_distances(reshape(p.z, 1, :), Zs)) >= p.dist_to_nearest_solve - 1e-9 for p in wave2)
        end

        @testset "Spearman helper" begin
            @test tied_ranks([1.0, 2.0, 2.0, 3.0]) == [1.0, 2.5, 2.5, 4.0]
            a = rand(rng, 50)
            @test spearman_correlation(a, exp.(a)) ≈ 1.0
            @test spearman_correlation(a, -a) ≈ -1.0
        end
    end
end

main()
