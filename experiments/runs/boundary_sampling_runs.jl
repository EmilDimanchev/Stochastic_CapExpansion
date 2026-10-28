#======================

Interior samplers: entry points and checks (src/model/boundary_sampling)

  test_boundary_sampling_units()        - solver-free unit tests (tests/boundary_sampling_tests.jl)
  test_gridfill_units()                 - solver-free unit tests (tests/gridfill_tests.jl)
  gradient_check("grad_chk")            - finite-difference check of the Benders gradients β_s
  run13_boundary_checks("run13")        - metric + sampler + bracket checks on run 13's finals
  boundary_sampling_test_laptop("bs_1") - end-to-end pareto run with Sampling "Method" => "Boundary"
  gridfill_test_laptop("gf_1")          - end-to-end pareto run with Sampling "Method" => "GridFill"

=======================#

##const RUN13_SUMMARY = joinpath(@__DIR__, "..", "..", ".claude", "files for claude", "Summary_pareto (11).csv")
#const RUN13_METRIC = joinpath(@__DIR__, "..", "..", ".claude", "files for claude", "run13_metric.csv")   # run-13 length scales, for comparison only

test_boundary_sampling_units(test_index = "boundary_units") = include(joinpath(@__DIR__, "..", "..", "tests", "boundary_sampling_tests.jl"))
test_gridfill_units(test_index = "gridfill_units") = include(joinpath(@__DIR__, "..", "..", "tests", "gridfill_tests.jl"))

# Laptop inputs, workers and subproblems, shared by the checks below. `fast` turns off the
# gas price uncertainty set (10 instead of 100 scenarios), keeping demand uncertainty on.
function _boundary_laptop_setup(; fast::Bool = false)
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")
    settings = load_settings(inputs_folder)
    fast && (settings["Gas price uncertainty"] = false)
    inputs = load_input_data(inputs_folder, settings)
    configure_parallel_workers!(settings)
    SPs = build_all_subproblems(inputs, settings)
    return inputs, settings, SPs
end

# End-to-end: full risk-pareto + MGA + wing caps, then the boundary sampler on the finals.
function boundary_sampling_test_laptop(test_index; n_samples::Int = 20, fast::Bool = false, sampling_overrides::AbstractDict = Dict())
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    mkpath(summary_folder)
    inputs, settings, SPs = _boundary_laptop_setup(; fast = fast)
    run_stochastic_exploration_risk_pareto(SPs, inputs, settings, joinpath(results_folder, "Pareto5"), summary_folder;
        budget_settings = Dict("Budget type" => "Transformed", "Tighten budget" => true),
        cap_settings = Dict("Cap wings" => true),
        mapping_settings = Dict("Mapping" => false),
        run_settings = Dict("Summary name" => "pareto"),
        sampling_settings = merge(Dict{String,Any}("Method" => "Boundary", "Samples" => n_samples, "Batch size" => n_samples, "Write gradients" => true), sampling_overrides))
end

function boundary_sampling_test_della(test_index; n_samples::Int = 2000, sampling_overrides::AbstractDict = Dict())
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della")
    if !isdir(results_folder)
        mkpath(results_folder)
    end
    if !isdir(summary_folder)
        mkpath(summary_folder)
    end
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)

    configure_parallel_workers!(settings)

    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)
    results_folder = joinpath(results_folder, "Boundary")
    run_stochastic_exploration_risk_pareto(SPs, inputs, settings, joinpath(results_folder, "Pareto5"), summary_folder;
        budget_settings = Dict("Budget type" => "Transformed", "Tighten budget" => true),
        cap_settings = Dict("Cap wings" => true, "Cap levels" => [1, 3]),
        mapping_settings = Dict("Mapping" => false),
        run_settings = Dict("Summary name" => "Boundary"),
        sampling_settings = merge(Dict{String,Any}("Method" => "Boundary", "Samples" => n_samples, "Batch size" => n_samples, "Write gradients" => true), sampling_overrides))
end

# End-to-end: full risk-pareto + MGA + wing caps, then GridFill on the finals. Laptop defaults
# are small (coarse 0.1 cells, batches of 8) so the whole loop, feedback included, runs in a
# few minutes; `fast` turns off the gas price uncertainty set (demand stays stochastic).
function gridfill_test_laptop(test_index; n_samples::Int = 24, fast::Bool = true, sampling_overrides::AbstractDict = Dict())
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    mkpath(summary_folder)
    inputs, settings, SPs = _boundary_laptop_setup(; fast = fast)
    run_stochastic_exploration_risk_pareto(SPs, inputs, settings, joinpath(results_folder, "Pareto5"), summary_folder;
        budget_settings = Dict("Budget type" => "Transformed", "Tighten budget" => true),
        cap_settings = Dict("Cap wings" => true),
        mapping_settings = Dict("Mapping" => false),
        run_settings = Dict("Summary name" => "pareto"),
        sampling_settings = merge(Dict{String,Any}("Method" => "GridFill", "Samples" => n_samples, "Batch size" => 8,
                                                   "GridFill" => Dict("Cell size" => 0.1)), sampling_overrides))
end

function gridfill_test_della(test_index; n_samples::Int = 2000, sampling_overrides::AbstractDict = Dict())
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della")
    mkpath(summary_folder)
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)
    configure_parallel_workers!(settings)
    SPs = build_all_subproblems(inputs, settings)
    results_folder = joinpath(results_folder, "GridFill")
    run_stochastic_exploration_risk_pareto(SPs, inputs, settings, joinpath(results_folder, "Pareto5"), summary_folder;
        budget_settings = Dict("Budget type" => "Transformed", "Tighten budget" => true),
        cap_settings = Dict("Cap wings" => true, "Cap levels" => [1, 3]),
        mapping_settings = Dict("Mapping" => false),
        run_settings = Dict("Summary name" => "GridFill"),
        sampling_settings = merge(Dict{String,Any}("Method" => "GridFill", "Samples" => n_samples, "Write gradients" => true), sampling_overrides))
end

# Gradient check: for a few final solves of run 13 and a couple of scenarios, perturb each
# capacity variable by ±h (h = h_frac · cap_max), re-solve that one scenario subproblem, and
# compare the one-sided differences with β_s. The op cost is convex in x, so a correct β is a
# subgradient: D⁻ <= β <= D⁺. Away from kinks D⁻ ≈ D⁺ ≈ β. Solves the main-process copies of
# the subproblems directly (one scenario at a time), so the worker caches are untouched.
function gradient_check(test_index = "grad_chk"; summary_csv::String = RUN13_SUMMARY, h_frac::Float64 = 1e-3,
                        solve_names = ["Risk_Weight_1.0", "Risk_Weight_0.0", "Random_1_Budget_1", "Random_1_Budget_1_CVaR_Cap_1"])
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    mkpath(results_folder)
    inputs, settings, SPs = _boundary_laptop_setup()
    R = length(inputs["Resources"])
    var_names = capacity_variable_names(inputs)
    α = settings["Value-at-Risk percent"]
    bank = load_finals_from_summary(summary_csv, inputs)
    cap_max = capacity_scale(reduce(vcat, [fs.x' for fs in bank]))
    keys_flat = vec([(s, f, k) for s in 1:size(SPs, 1), f in 1:size(SPs, 2), k in 1:size(SPs, 3)])

    function solve_one(key, x)
        SP = SPs[key...]
        set_parameter_value.(SP[:x], x[1:R])
        set_parameter_value.(SP[:x_line], x[R+1:end])
        out = run_subproblem(SP, inputs, settings; minimal_payload = true)
        return out["SP objective"], vcat(out["capacity dual"], out["line dual"]), out["coeff"] == 1
    end

    df = DataFrame(Solve = String[], Scenario = String[], Variable = String[], Beta = Float64[], D_minus = Float64[], D_plus = Float64[],
                   Rel_mismatch_central = Float64[], Kink = Bool[], Subgradient_ok = Bool[])
    for name in solve_names
        i = findfirst(fs -> fs.name == name, bank)
        i === nothing && (@warn("Solve $name not in $summary_csv, skipped"); continue)
        fs = bank[i]
        # two scenarios: the worst one (in the CVaR tail) and the median one
        order = sortperm(fs.θ; rev = true)
        for si in unique([order[1], order[div(length(order), 2)]])
            key = keys_flat[si]
            θ0, β, ok = solve_one(key, fs.x)
            ok || (@warn("Base solve not optimal for $name $key"); continue)
            for v in eachindex(var_names)
                h = h_frac * cap_max[v]
                xp = copy(fs.x); xp[v] += h
                θp, _, okp = solve_one(key, xp)
                D_plus = (θp - θ0) / h
                if fs.x[v] >= h
                    xm = copy(fs.x); xm[v] -= h
                    θm, _, _ = solve_one(key, xm)
                    D_minus = (θ0 - θm) / h
                else
                    D_minus = -Inf     # at the lower bound: no backward step, only D⁺ is checked
                end
                # tolerance: LP objective accuracy (~1e-7 relative) turned into a slope error, plus 0.1% of |β|
                tol = 1e-7 * max(abs(θ0), 1.0) / h + 1e-3 * abs(β[v])
                central = isfinite(D_minus) ? (D_plus + D_minus) / 2 : D_plus
                kink = isfinite(D_minus) && (D_plus - D_minus) > 1e-2 * max(abs(D_plus), abs(D_minus), 1.0)
                sub_ok = okp && (D_minus - tol <= β[v] <= D_plus + tol)
                push!(df, (name, string(key), var_names[v], β[v], D_minus, D_plus, abs(central - β[v]) / max(abs(central), 1.0), kink, sub_ok))
            end
        end
    end
    CSV.write(joinpath(results_folder, "gradient_check.csv"), df)
    bad = df[.!df.Subgradient_ok, :]
    @info("Gradient check: $(nrow(df)) comparisons, $(count(df.Kink)) at kinks (D⁻ ≠ D⁺), $(nrow(bad)) violate D⁻ <= β <= D⁺. Median relative mismatch to the central difference (non-kink): $(isempty(df[.!df.Kink, :]) ? NaN : median(df[.!df.Kink, :Rel_mismatch_central])).")
    nrow(bad) > 0 && show(stdout, bad; allrows = true)
    return df
end

# Metric, sampler and bracket checks on run 13 (d = 14), against the length-scale metric.
# Run 13 saved no gradients, so they are recovered once by re-solving the laptop subproblems
# at each final capacity and cached in outputs/Test_<index>/run13_bank.jls.
function run13_boundary_checks(test_index = "run13"; summary_csv::String = RUN13_SUMMARY, n_samples::Int = 300, seed::Int = 0, evaluate::Bool = true)
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    mkpath(summary_folder)
    inputs, settings, SPs = _boundary_laptop_setup()
    var_names = capacity_variable_names(inputs)
    α = settings["Value-at-Risk percent"]
    cfg = boundary_cfg(resolve_settings(settings, "Sampling settings", Dict("Boundary" => Dict("Metric file" => RUN13_METRIC))))

    cache = joinpath(results_folder, "run13_bank.jls")
    if isfile(cache)
        bank = deserialize(cache)
        @info("Loaded cached gradient bank from $cache")
    else
        bank = load_finals_from_summary(summary_csv, inputs)
        @time recompute_gradients!(bank, SPs, inputs, settings)
        serialize(cache, bank)
    end
    finals = [fs for fs in bank if fs.ok]
    @info("$(length(finals)) / $(length(bank)) final solves usable (all scenarios optimal)")
    X = reduce(vcat, [fs.x' for fs in finals])

    # ---- Metric check -------------------------------------------------------------------
    m_as = build_active_subspace(finals, var_names, cfg, α)
    write_metric(joinpath(results_folder, "Metric"), m_as; tag = "active_subspace")
    m_ls = lengthscale_metric(cfg["Metric file"], var_names)
    println("\nActive subspace: r = $(m_as.r) of d = $(length(var_names)); anchor gaps ΔC = $(m_as.ΔC), ΔV = $(m_as.ΔV)")
    println("Spectrum: ", round.(m_as.λ_all; sigdigits = 3))
    println("Cumulative share: ", round.(cumsum(m_as.λ_all) ./ sum(m_as.λ_all); digits = 3))
    for j in 1:min(m_as.r, 5)
        top = sortperm(abs.(m_as.W[:, j]); rev = true)[1:4]
        println("  direction $j (λ = $(round(m_as.λ_all[j]; sigdigits = 3))): ", join(["$(var_names[t]) $(round(m_as.W[t, j]; digits = 2))" for t in top], ", "))
    end
    # Per-variable sensitivity: how far z moves per unit of x̃ along each variable
    sens_as = vec(sqrt.(sum((m_as.W .^ 2) .* m_as.λ_all[1:m_as.r]'; dims = 2)))
    sens_ls = m_as.cap_max ./ m_ls.scale          # = (bank cap_max / file cap_max) / length scale
    by_as = sortperm(sens_as; rev = true)
    by_ls = sortperm(sens_ls; rev = true)
    println("\nMost sensitive variables, active subspace: ", join(var_names[by_as[1:6]], ", "))
    println("Most sensitive variables, length scales:   ", join(var_names[by_ls[1:6]], ", "))
    println("Spearman correlation of per-variable sensitivities: ", round(spearman_correlation(sens_as, sens_ls); digits = 3))

    Z_as, Z_ls = metric_coords(m_as, X), metric_coords(m_ls, X)
    upper(D) = [D[i, j] for j in 1:size(D, 1) for i in 1:j-1]
    ρ = spearman_correlation(upper(pairwise_distances(Z_as, Z_as)), upper(pairwise_distances(Z_ls, Z_ls)))
    println("Spearman correlation of pairwise distances (active subspace vs length scale): ", round(ρ; digits = 3))

    # ---- Sampler check under both metrics -----------------------------------------------
    picks_by_metric = Dict{String, Any}()
    for (label, m, Z) in (("active_subspace", m_as, Z_as), ("lengthscale_file", m_ls, Z_ls))
        rng = Xoshiro(seed)
        picks = select_boundary_samples(finals, Z, cfg, rng, α; n_new = n_samples)
        spread = sampler_spread_diagnostics(picks, Z, Xoshiro(seed + 1))
        counts = Dict(r => count(p -> p.region == r, picks) for r in BOUNDARY_REGION_ORDER)
        println("\n[$label] picks per base region: ", join(["$r $(counts[r])" for r in BOUNDARY_REGION_ORDER], ", "))
        println("[$label] median pick distance / median final NN distance = $(round(spread.ratio; digits = 3)); random-pair midpoints = $(round(spread.midpoint_ratio; digits = 3))")
        picks_by_metric[label] = picks
    end

    # ---- Evaluate the active-subspace picks and check their brackets ----------------------
    brackets = nothing
    if evaluate
        picks = picks_by_metric["active_subspace"]
        rows_cap, rows_syscost, rows_emissions, solved = evaluate_picks(picks, SPs, inputs, settings)
        write_exploration_results!(rows_cap, rows_syscost, rows_emissions, summary_folder, [p.name for p in picks], "run13_interior")
        brackets = check_brackets(picks, solved, α, settings["Convergence tolerance"])
        CSV.write(joinpath(summary_folder, "Interior_brackets_run13.csv"), brackets)
    end
    return (metric = m_as, metric_ls = m_ls, picks = picks_by_metric, brackets = brackets, spearman = ρ)
end
