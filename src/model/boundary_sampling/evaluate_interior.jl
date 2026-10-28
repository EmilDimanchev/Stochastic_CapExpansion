using LinearAlgebra, Random, Statistics, DataFrames, CSV

#==========================================================================================

Evaluating the boundary picks and checking their bounds

evaluate_picks: each pick goes through the same fixed-capacity path as the existing
interior sampler: run_all_subproblems (parallel over scenarios on the workers that own them,
minimal_payload = false so prices/costs by zone are available), then
make_results_mapping_dfs for the Summary row. Investment cost is computed directly; it is
linear in capacity, so no master solve is needed. Each solved pick is also stored as a
FinalSolve (with its own gradients) so it can feed later waves or a metric refit.

run_boundary_sampling!: the whole interior phase (metric → select → evaluate, per wave).

==========================================================================================#

# Lineage + bound columns appended to a pick's Summary row (missing for every MGA row)
function add_lineage_columns!(df_syscost::DataFrame, pick)
    df_syscost[!, :Base_Region] = [pick.region]
    df_syscost[!, :Parents] = [join(pick.parent_names, ";")]
    df_syscost[!, :Weights] = [join(string.(round.(pick.weights; digits = 4)), ";")]
    df_syscost[!, :Wave] = [pick.wave]
    df_syscost[!, :Dist_To_Nearest_Solve] = [pick.dist_to_nearest_solve]
    df_syscost[!, :Blend_EV_OpCost] = [pick.blend_ev]
    df_syscost[!, :Blend_CVaR_OpCost] = [pick.blend_cvar]
    df_syscost[!, :Blend_Investment_Cost] = [pick.blend_inv]
    df_syscost[!, :LB_EV_OpCost] = [pick.lb_ev]
    df_syscost[!, :LB_CVaR_OpCost] = [pick.lb_cvar]
    return df_syscost
end

function evaluate_picks(picks::AbstractVector, SPs, inputs::Dict, settings::Dict)
    R = length(inputs["Resources"])
    p = scenario_probabilities(inputs)
    α = settings["Value-at-Risk percent"]
    rows_cap, rows_syscost, rows_emissions = DataFrame[], DataFrame[], DataFrame[]
    solved = FinalSolve[]
    for pick in picks
        cap, line = pick.x[1:R], pick.x[R+1:end]
        sp = run_all_subproblems(SPs, inputs, settings, cap, line; minimal_payload = false)
        inv, inv_by_zone = investment_cost(pick.x, inputs)
        fs = final_solve_from_sp_outputs(pick.name, pick.x, sp, p, inv)
        fs.region = "interior_sample"
        ev, cvar = solve_ev(fs), solve_cvar(fs, α)      # same rule as evaluate_subproblems
        df_cap, df_syscost, df_emissions = make_results_mapping_dfs(cap, line, sp, inv, inv_by_zone, cvar, ev, inputs, settings)
        add_lineage_columns!(df_syscost, pick)
        push!(rows_cap, df_cap); push!(rows_syscost, df_syscost); push!(rows_emissions, df_emissions)
        push!(solved, fs)
        @info("$(pick.name) ($(pick.region)): EV op = $(ev) in [$(pick.lb_ev), $(pick.blend_ev)], CVaR op = $(cvar) in [$(pick.lb_cvar), $(pick.blend_cvar)]")
    end
    return rows_cap, rows_syscost, rows_emissions, solved
end

#---- Bracket check: LB − tol <= solved <= UB + tol ---------------------------------------
# A violation means the convexity assumption fails, a gradient is wrong, or a subproblem did
# not converge. tol is relative (rel_tol · |value|).
function check_brackets(picks::AbstractVector, solved::AbstractVector, α::Real, rel_tol::Real)
    df = DataFrame(Index = String[], Base_Region = String[],
                   LB_EV = Float64[], EV = Float64[], UB_EV = Float64[],
                   LB_CVaR = Float64[], CVaR = Float64[], UB_CVaR = Float64[],
                   Width_EV_rel = Float64[], Width_CVaR_rel = Float64[],
                   Inv_solved = Float64[], Inv_blend = Float64[], Violation = String[])
    for (pick, fs) in zip(picks, solved)
        ev, cvar = solve_ev(fs), solve_cvar(fs, α)
        tol(v) = rel_tol * max(abs(v), 1.0)
        issues = String[]
        !isnan(pick.lb_ev) && ev < pick.lb_ev - tol(ev) && push!(issues, "EV<LB")
        ev > pick.blend_ev + tol(ev) && push!(issues, "EV>UB")
        !isnan(pick.lb_cvar) && cvar < pick.lb_cvar - tol(cvar) && push!(issues, "CVaR<LB")
        cvar > pick.blend_cvar + tol(cvar) && push!(issues, "CVaR>UB")
        abs(fs.inv_cost - pick.blend_inv) > tol(fs.inv_cost) && push!(issues, "Inv!=blend")
        fs.ok || push!(issues, "SP not optimal")
        push!(df, (pick.name, pick.region, pick.lb_ev, ev, pick.blend_ev, pick.lb_cvar, cvar, pick.blend_cvar,
                   (pick.blend_ev - pick.lb_ev) / abs(ev), (pick.blend_cvar - pick.lb_cvar) / abs(cvar),
                   fs.inv_cost, pick.blend_inv, join(issues, ";")))
        isempty(issues) || @warn("Bracket violation for $(pick.name): $(join(issues, ", ")) (EV $(pick.lb_ev) <= $ev <= $(pick.blend_ev); CVaR $(pick.lb_cvar) <= $cvar <= $(pick.blend_cvar))")
    end
    n_bad = count(!isempty, df.Violation)
    if nrow(df) > 0
        @info("Bracket check: $n_bad / $(nrow(df)) picks violate their bounds. Relative bracket width EV: median $(round(median(df.Width_EV_rel); sigdigits = 3)), max $(round(maximum(df.Width_EV_rel); sigdigits = 3)); CVaR: median $(round(median(df.Width_CVaR_rel); sigdigits = 3)), max $(round(maximum(df.Width_CVaR_rel); sigdigits = 3)).")
    end
    return df
end

#---- The whole interior phase -------------------------------------------------------------
# The boundary sampler's own config: the "Boundary" sub-block of Sampling settings plus the
# shared candidate pool size.
boundary_cfg(sampling::AbstractDict) = merge(Dict{String,Any}(sampling["Boundary"]), Dict{String,Any}("Candidates" => sampling["Candidates"]))

# bank: FinalSolves recorded during the MGA phase. sampling: the resolved "Sampling settings"
# block; "Samples" picks in total, split into waves of at most "Batch size". The seed is
# settings["Seed"], as for the Delaunay sampler. Returns the Summary rows and labels for the
# evaluated picks (to be appended after the MGA rows), plus the picks and solved points.
function run_boundary_sampling!(bank::AbstractVector, SPs, inputs::Dict, settings::Dict, sampling::AbstractDict, results_folder::String, summary_folder::String, summary_name::String)
    cfg = boundary_cfg(sampling)
    n_samples = sampling["Samples"]
    α = settings["Value-at-Risk percent"]
    var_names = capacity_variable_names(inputs)
    assign_regions!(bank)
    finals = [fs for fs in bank if fs.region != "interior_sample"]
    X_of(solves) = reduce(vcat, [fs.x' for fs in solves])
    for region in vcat(BOUNDARY_REGION_ORDER, "anchor")
        @info("Final solves in region $region: $(count(fs -> fs.region == region, finals))")
    end

    rng = Xoshiro(settings["Seed"])      # one seeded stream for the whole phase: reproducible
    metric = build_metric(finals, var_names, cfg, α)
    write_metric(joinpath(results_folder, "Metric"), metric; tag = "wave_0")

    W = cld(n_samples, sampling["Batch size"])
    n_total = n_samples
    rows_cap, rows_syscost, rows_emissions = DataFrame[], DataFrame[], DataFrame[]
    picks_all, solved_all = BoundaryPick[], FinalSolve[]
    for wave in 1:W
        n_wave = div(n_total, W) + (wave <= n_total % W ? 1 : 0)
        Z_finals = metric_coords(metric, X_of(finals))
        Z_solved = isempty(solved_all) ? zeros(0, size(Z_finals, 2)) : metric_coords(metric, X_of(solved_all))
        picks = select_boundary_samples(finals, Z_finals, cfg, rng, α; n_new = n_wave, extra_ref_Z = Z_solved,
                                        first_index = length(picks_all) + 1, wave = wave)
        for region in BOUNDARY_REGION_ORDER
            @info("Wave $wave picks with base region $region: $(count(p -> p.region == region, picks))")
        end
        spread = sampler_spread_diagnostics(picks, Z_finals[[fs.ok for fs in finals], :], rng)
        @info("Wave $wave spread: median pick distance / median final NN distance = $(round(spread.ratio; digits = 3)); random-pair midpoints: $(round(spread.midpoint_ratio; digits = 3))")

        c, s, e, solved = evaluate_picks(picks, SPs, inputs, settings)
        append!(rows_cap, c); append!(rows_syscost, s); append!(rows_emissions, e)
        append!(picks_all, picks); append!(solved_all, solved)
        if sampling["Write gradients"]
            for fs in solved
                write_gradient_file(joinpath(results_folder, "Interior", fs.name), fs, var_names)
            end
        end

        # Optional refit: rebuild the subspace with the interior gradients included
        if cfg["Refit metric each wave"] && wave < W && metric.mode == "active_subspace"
            metric = build_metric(vcat(finals, solved_all), var_names, cfg, α)
            write_metric(joinpath(results_folder, "Metric"), metric; tag = "wave_$wave")
        end
    end
    append!(bank, solved_all)

    brackets = check_brackets(picks_all, solved_all, α, settings["Convergence tolerance"])
    mkpath(summary_folder)
    CSV.write(joinpath(summary_folder, "Interior_brackets_" * summary_name * ".csv"), brackets)
    labels = [p.name for p in picks_all]
    return rows_cap, rows_syscost, rows_emissions, labels, picks_all, solved_all
end
