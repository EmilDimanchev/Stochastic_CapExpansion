using Distributed

@everywhere include("../src/Stochastic_CapExpansion.jl")

@everywhere using .Stochastic_CapExpansion
@everywhere using Revise, JuMP, Gurobi, HiGHS, Ipopt, DataFrames, CSV, YAML, Random, LinearAlgebra, Combinatorics, Dates, Distributions, Surrogates, DelaunayTriangulation
using Clarabel

#=========================

Utility functions for logging and managing memory, configuring parallel workers, and running stochastic exploration with Benders and MGA

=========================#

function log_result_memory!(label::String, output::Dict)
    size_mb = round(Base.summarysize(output) / 1024^2; digits = 2)
    @info("Approx memory for $(label): $(size_mb) MiB")
end

function release_heavy_payload!(output::Dict)
    for key in ("SPs", "SPs_eval", "All MP outputs per iteration", "All SP outputs per iteration")
        if haskey(output, key)
            delete!(output, key)
        end
    end
    GC.gc(false)
    return nothing
end

function configure_parallel_workers!(settings::Dict)
    if settings["Parallel flag"]
        desired_workers = haskey(settings, "Workers") ? settings["Workers"] : max(Sys.CPU_THREADS - 1, 1)
        wkers = nworkers() == 1 ? 0 : nworkers()
        if wkers < desired_workers
            println("Current workers: ", wkers, ". Adding ", desired_workers - wkers, " workers for parallel processing...")
            addprocs(desired_workers - wkers)
        end

        project_root = abspath(joinpath(@__DIR__, ".."))
        init_expr = quote
            if !isdefined(Main, :Stochastic_CapExpansion)
                include(joinpath($project_root, "src", "Stochastic_CapExpansion.jl"))
            end
            using .Stochastic_CapExpansion
            using Revise, JuMP, Ipopt, Gurobi, HiGHS, DataFrames, CSV, YAML, Random, LinearAlgebra, Combinatorics, Dates, Distributions
            nothing
        end
        @sync for pid in workers()
            @async remotecall_wait(Core.eval, pid, Main, init_expr)
        end

        settings["Workers"] = nworkers()
        @info("Running with parallelization using $(nworkers()) distributed workers")
    else
        @info("Running without parallelization")
    end
end



#===============================

Helpers for main pareto + exploration run



================================#


function compute_budget_floor(extreme_values::Dict, transform_cvar::Float64, transform_sys::Float64, floor_offset::Float64)
    transformed_vals = [transform_sys*ev["System Expected"] + transform_cvar*ev["CVaR"] for ev in values(extreme_values)]
    isempty(transformed_vals) && error("Need at least one risk-pareto anchor point to compute a budget floor.")
    return minimum(transformed_vals) + floor_offset
end

function compute_budget_tighten_schedule(budget_start::Float64, budget_floor::Float64, step_size::Float64)
    step_size <= 0 && error("Budget tighten step size must be positive.")
    budget_floor > budget_start && error("Budget floor ($budget_floor) must be <= starting budget ($budget_start) — check budget_offset vs. the anchor points' own transformed costs.")
    n_steps = ceil(Int, (budget_start - budget_floor) / step_size)
    return [max(budget_floor, budget_start - step_size*k) for k in 0:n_steps]
end


#-------------------------------------------------------------------------------------------
# Helper 1: piecewise-linear interpolation (Base Julia has no built-in interp function).
# Given sorted breakpoints xs with values ys, return the linearly interpolated y at x.
# Values of x outside [xs[1], xs[end]] are clamped to the end values.
#-------------------------------------------------------------------------------------------
function linear_interp(xs::Vector{Float64}, ys::Vector{Float64}, x::Real)
    # Clamp to the ends of the curve
    if x <= xs[1]
        return ys[1]
    elseif x >= xs[end]
        return ys[end]
    end

    # Index i of the last breakpoint with xs[i] ≤ x, so x lies in [xs[i], xs[i+1])
    i = searchsortedlast(xs, x)

    # Guard against a zero-length segment (two breakpoints at the same x)
    if xs[i+1] == xs[i]
        return ys[i]
    end

    # Fraction of the way from xs[i] to xs[i+1], then blend the two y values
    t = (x - xs[i]) / (xs[i+1] - xs[i])
    return ys[i] + t * (ys[i+1] - ys[i])
end


#-------------------------------------------------------------------------------------------
# Helper 2: lens end. Walk the polyline vertices in the given order (from one anchor inward)
# and return the arc position of the first point where T ≤ B. Within the crossing segment
# T is linear, so the exact crossing point is found by linear interpolation.
# Returns `nothing` if no point on the polyline satisfies T ≤ B.
#-------------------------------------------------------------------------------------------
function lens_end(B::Real, order::Vector{Int}, Tv::Vector{Float64}, arc::Vector{Float64})
    for k in 1:length(order)-1
        a = order[k]      # current vertex
        b = order[k+1]    # next vertex inward

        # The current vertex already satisfies the budget: the lens starts here
        if Tv[a] <= B
            return arc[a]
        end

        # The budget line is crossed inside segment a -> b: solve for the crossing point
        if Tv[b] <= B
            t = (Tv[a] - B) / (Tv[a] - Tv[b])
            return arc[a] + t * (arc[b] - arc[a])
        end
    end

    # Only the final vertex is left to check
    if Tv[order[end]] <= B
        return arc[order[end]]
    end
    return nothing
end


#-------------------------------------------------------------------------------------------
# Helper 3: cap values for one wing. Places n evenly spaced caps between the lens end
# (tight_end) and the central edge (edge), excluding the lens end itself. n is limited by
# n_max and by the minimum spacing min_gap_frac * span. Caps are returned loose -> tight,
# which is the order to solve them in.
#-------------------------------------------------------------------------------------------
function wing_caps(tight_end::Real, edge::Real, span::Real, n_max::Int, min_gap_frac::Real)
    width = edge - tight_end

    # How many caps fit while keeping adjacent caps at least min_gap_frac * span apart
    n = min(n_max, Int(floor(width / (min_gap_frac * span))))
    if n <= 0
        return Float64[]
    end

    # i = n gives the loosest cap (at the central edge); i = 1 the tightest
    caps = Float64[]
    for i in n:-1:1
        push!(caps, tight_end + width * i / n)
    end
    return caps
end


#-------------------------------------------------------------------------------------------
# Main function.
#
# Inputs
# - extreme_values: Dict(risk_weight => Dict("System Expected" => …, "CVaR" => …)) for ALL
#   risk-weight solves. Intermediate weights are required to locate the lens ends and the
#   tangency point. The key convention does not matter: the EV-optimal solve is the one with
#   the lowest System Expected.
# - budgets: transformed budget levels, in the same units as k_ev*SystemExpected + k_cvar*CVaR.
# - k_ev, k_cvar: weights of the transformed budget.
# - n_max: maximum number of caps per wing per level.
# - min_gap_frac: minimum spacing between adjacent caps, as a fraction of that metric's span
#   (0.03 ≈ 2.5× the Benders overshoot in the reference run).
# - edge_pullback: moves each central edge from the tangency toward its anchor, as a fraction
#   of that metric's span. 0.0 = edges exactly at the tangency.
#
# Returns Dict(B => Dict("System_Expected" => caps, "CVaR" => caps)). Iterate over `budgets`
# (not over the Dict) to keep level order.
#
# Feasibility: every polyline point is a convex combination of solved portfolios, so by
# convexity of the operational problem some feasible portfolio is at least as good on both
# metrics. Each capped region {T ≤ B, cap} therefore contains a feasible point.
#-------------------------------------------------------------------------------------------
function compute_budget_schedule_wings(extreme_values::Dict, budgets::AbstractVector;
                                       k_ev::Real, k_cvar::Real,
                                       n_max::Int = 3, min_gap_frac::Real = 0.03,
                                       edge_pullback::Real = 0.05)

    # --- Block 1: collect the risk-weight solves as (System Expected, CVaR) pairs and sort
    #     them by System Expected, so the polyline runs EV-optimal -> CVaR-optimal.
    pts = Tuple{Float64,Float64}[]
    for v in values(extreme_values)
        push!(pts, (Float64(v["System Expected"]), Float64(v["CVaR"])))
    end
    pts = sort(pts; by = first)
    if length(pts) < 3
        error("need intermediate risk-weight solves to locate lens ends and the tangency point")
    end
    ev = [first(p) for p in pts]   # System Expected at each vertex (increasing)
    cv = [last(p) for p in pts]    # CVaR at each vertex (decreasing along a proper frontier)
    n_pts = length(pts)

    # --- Block 2: spans of each metric between the two extreme solves. These normalize the
    #     axes for arc length and set the minimum cap spacing.
    ev_span = ev[end] - ev[1]
    cv_span = cv[1] - cv[end]
    if !(ev_span > 0 && cv_span > 0)
        error("extreme solves must trade off EV against CVaR")
    end

    # --- Block 3: normalized arc length along the polyline (0 at EV-optimal, 1 at CVaR-optimal).
    #     Each segment length is measured with both axes scaled by their spans.
    seg = Float64[]
    for i in 1:n_pts-1
        push!(seg, hypot((ev[i+1] - ev[i]) / ev_span, (cv[i+1] - cv[i]) / cv_span))
    end
    arc = vcat(0.0, cumsum(seg)) ./ sum(seg)

    # --- Block 4: transformed budget value T at each polyline vertex.
    Tv = [k_ev * ev[i] + k_cvar * cv[i] for i in 1:n_pts]

    # --- Block 5: central edges at the tangency point (vertex with the lowest T), optionally
    #     pulled back toward each anchor. These are the same for every budget level.
    i_tan = argmin(Tv)
    ev_central_edge = ev[i_tan] - edge_pullback * ev_span   # EV-wing caps stop here
    cv_central_edge = cv[i_tan] - edge_pullback * cv_span   # CVaR-wing caps stop here

    # --- Block 6: for each budget level, find both lens ends and place the wing caps.
    order_from_ev   = collect(1:n_pts)       # walk inward from the EV-optimal anchor
    order_from_cvar = collect(n_pts:-1:1)    # walk inward from the CVaR-optimal anchor
    schedule = Dict{Float64,Dict{String,Vector{Float64}}}()

    for B in budgets
        s_ev = lens_end(B, order_from_ev, Tv, arc)
        s_cv = lens_end(B, order_from_cvar, Tv, arc)

        # EV wing: caps on System Expected, only if the lens end lies outside the central edge
        caps_ev = Float64[]
        if !isnothing(s_ev)
            ev_end = linear_interp(arc, ev, s_ev)
            if ev_end < ev_central_edge
                caps_ev = wing_caps(ev_end, ev_central_edge, ev_span, n_max, min_gap_frac)
            end
        end

        # CVaR wing: caps on CVaR, only if the lens end lies outside the central edge
        caps_cv = Float64[]
        if !isnothing(s_cv)
            cv_end = linear_interp(arc, cv, s_cv)
            if cv_end < cv_central_edge
                caps_cv = wing_caps(cv_end, cv_central_edge, cv_span, n_max, min_gap_frac)
            end
        end

        schedule[Float64(B)] = Dict("System_Expected" => caps_ev, "CVaR" => caps_cv)
    end
    return schedule
end

#============================


Main MGA functions for running stochastic exploration with Benders and MGA, including base runs and iterative exploration with budget constraints, cut management, and result logging


==============================#


function run_stochastic_exploration_risk_pareto(SPs::Array{Model, 3}, inputs::Dict, settings::Dict, results_folder::String, summary_folder::String; budget_multiplier::Float64 = 1.10, budget_offset::Float64 = 1.0, vector_set::Union{AbstractVector, Nothing} = nothing, summary_name::String = "new_setup", Eval_SPs = nothing, mapping = false, n_samples = 100, budget_type = "Transformed", tighten_budget::Bool = false, budget_tighten_step_size::Float64 = 0.1, floor_offset::Float64 = 0.01, cap_wings::Bool = true)

    #configure_parallel_workers!(settings)

    # Result containers
    results_cap = []
    results_syscost = []
    results_emissions = []
    run_labels = []
    gaps = []
    cuts_to_keep = []

    # Model Settings
    iterations = settings["Iterations"]
    risk_aversion_weights = [[i for i in 0.5:-0.1:0.1]; [i for i in 0.6:0.1:0.9]; 0.0; 1.0]
    risk_aversion_weight = settings["Risk aversion weight"] ### set as anchor point
    VaR_percent = settings["Value-at-Risk percent"] ### Currently set consistently across runs and ahead of time
    settings["Risk aversion flag"] = true
    scaling = settings["Scaling factor cost"]
    # Create and set expected value model
    extreme_values = Dict()
    # initialize quantities in function memory (not subsequent if/for statements)
    budget_val_transform = 0.0
    transform_sys = 0.0
    transform_cvar = 0.0

    P_s = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]
    R = length(inputs["Resources"])

    MP = build_planning_model(inputs, settings; risk_aversion_weight = risk_aversion_weight)

    # Shared across every risk weight / MGA direction solved against this MP so that cuts
    # pruned while solving one problem get a fresh look (and are rebuilt losslessly, since
    # deactivate_cuts now archives the exact constraint rather than discarding it) when a
    # later problem starts.
    cut_archive = Dict{String, Any}()

    for (i, risk) in enumerate(risk_aversion_weights)
        # Base Runs
        case_name = "Risk_Weight_"*string(risk)
        set_objective_bendersMP!(MP, "System_Weighted_CVaR", inputs, settings; obj_weight = risk)
        output = benders_algorithm(inputs, settings, MP, SPs, case_name; Eval_SPs = Eval_SPs, mapping = mapping, risk_aversion_weight = risk, cut_archive = cut_archive)
        log_result_memory!(case_name*" output", output)
        gap_cvar = output["Gaps"]
        push!(run_labels, case_name)
        push!(gaps, gap_cvar)
        # Write results
        results_destination = joinpath(results_folder, case_name)
        df_cap, df_syscost, df_emissions = write_results_benders(output, inputs, settings, results_destination)
        if mapping && risk >= 0.00001 # > 0 with tolerance for machine precision
        # do not map intermediate risk = 0 solutions because they are unbounded in EV, which can distort the outputs
            temp_df_cap, temp_df_syscost, temp_df_emissions = write_mapping_results(output, inputs, settings)
            push!(results_cap, temp_df_cap)
            push!(results_syscost, temp_df_syscost)
            push!(results_emissions, temp_df_emissions)
        else
            push!(results_cap, df_cap)
            push!(results_syscost, df_syscost)
            push!(results_emissions, df_emissions)
        end
        extreme_values[risk] = Dict("CVaR" => output["CVaR"]/scaling, "System Expected" => (output["Expected Value"] + output["MP"]["Inv_cost"])/scaling)

        @info("Risk weight $risk solution has investment cost of $(output["MP"]["Inv_cost"])")
        @info("Expected value system cost of " * "Risk weight $risk" * " solution: $(output["Expected Value"] + output["MP"]["Inv_cost"])")
        @info("Risk adjusted system cost of " * "Risk weight $risk" * " solution: $((1-risk_aversion_weight)*output["CVaR"] + risk_aversion_weight*output["Expected Value"]+ output["MP"]["Inv_cost"])")
        
        if settings["Cut deactivation strategy"] == "in mga"
            
            cuts_to_keep = [name(con) for con in all_constraints(MP, include_variable_in_set_constraints=false) if ((startswith(string(con), "optimality_cut_") && (split(string(con), "_")[3] == string(risk))) || (startswith(string(con), "cvar_tail_cuts_") && (split(string(con), "_")[3] == string(risk))))]
        
            if mapping
                cuts_to_keep = filter(cut -> parse(Int, split(string(cut), "_")[end-1]) <= output["first_write"] + 10, cuts_to_keep)
            end
            cuts_to_keep = manage_cuts(MP, cuts_to_keep)
        end

        release_heavy_payload!(output)
    end

    if settings["Capacity Exploration"]
        budgets = Dict()
        constraint_dict = Dict()
        # Wing caps only make sense with a tightening schedule on the Transformed budget
        do_cap_wings = cap_wings && tighten_budget && budget_type == "Transformed"
        cap_schedule = Dict{Float64,Dict{String,Vector{Float64}}}()

        if tighten_budget && budget_type != "Transformed"
            error("tighten_budget=true is only supported for budget_type == \"Transformed\" (got \"$budget_type\").")
        end

        if mapping && budget_type != "Transformed"
            error("mapping=true is only supported for budget_type == \"Transformed\" (got \"$budget_type\") - sample_interior_delauney filters against the Transformed budget.")
        end

        if budget_type == "Transformed"
            transform_cvar = 1/(extreme_values[1.0]["CVaR"] - extreme_values[0.0]["CVaR"])
            transform_sys = 1/(extreme_values[0.0]["System Expected"] - extreme_values[1.0]["System Expected"])
            budget_val_transform = transform_cvar*extreme_values[0.0]["CVaR"] + transform_sys*extreme_values[1.0]["System Expected"] + budget_offset
            @info("Budget value for transformed budget constraint: ", budget_val_transform)
            budgets, constraint_dict = add_budget_constraint_bendersMP(MP, budget_val_transform, "Transformed", budgets; extreme_values = extreme_values, constraint_dict = constraint_dict)
            if tighten_budget
                budget_floor = compute_budget_floor(extreme_values, transform_cvar, transform_sys, floor_offset)
                budget_schedule = compute_budget_tighten_schedule(budget_val_transform, budget_floor, budget_tighten_step_size)
                @info("Budget tightening enabled: $(length(budget_schedule)) levels per MGA vector, from $(budget_val_transform) down to $(budget_floor) in steps of $(budget_tighten_step_size).")
                if do_cap_wings
                    cap_schedule = compute_budget_schedule_wings(extreme_values, budget_schedule; k_ev = transform_sys, k_cvar = transform_cvar)
                end
            end
        elseif budget_type == "Box"
            budgets, constraint_dict = add_budget_constraint_bendersMP(MP, extreme_values[1.0]["CVaR"], "CVaR", budgets; constraint_dict = constraint_dict)
            budgets, constraint_dict = add_budget_constraint_bendersMP(MP, extreme_values[0.0]["System Expected"]*budget_multiplier, "System_Expected", budgets; constraint_dict = constraint_dict)
        else
            error("Unknown budget type: $budget_type")
        end

        if settings["Cut deactivation strategy"] == "in mga"
            cuts = [name(con) for con in all_constraints(MP, include_variable_in_set_constraints=false) if (startswith(string(con), "optimality_cut_") || startswith(string(con), "cvar_tail_cuts_"))]
            cuts_to_keep = copy(cuts)

            if mapping
                cuts_to_keep = filter(cut -> parse(Int, split(string(cut), "_")[end-1]) <= output["first_write"] + 10, cuts_to_keep)
            end

        end
        vectors = vector_set !== nothing ? vector_set : generate_weights(iterations, length(MP[:x])+length(MP[:x_line]), settings["Vector Type"], settings)
        #@info("Keeping $(length(cuts_to_keep)) cuts for MGA iterations")
        for iteration in 1:iterations
            set_objective_bendersMP!(MP, "Capacity", inputs, settings; set_coeffs = vectors[iteration])
            if settings["Cut deactivation strategy"] == "in mga"
                cuts_to_keep = manage_cuts(MP, cuts_to_keep)
            end
            #cuts_to_keep = manage_cuts(MP, cuts_to_keep)
            # introduce new variable to track if this is the not the first tightening. Used to indicate mapping
            map_bool = false

            levels = tighten_budget ? budget_schedule : [budget_val_transform]
            local output_random, avg_time_mp
            for (level, budget_level_val) in enumerate(levels)
                if tighten_budget
                    budgets = update_budget_constraint_bendersMP!(MP, budget_level_val, "Transformed", budgets)
                end
                run_name = tighten_budget ? "Random_$(iteration)_Budget_$(level)" : "Random_"*string(iteration)
                map_bool = true
                try
                    output_random = mga_benders(inputs, settings, MP, SPs, budgets, run_name; Eval_SPs = Eval_SPs, mapping = map_bool, cut_archive = cut_archive)
                    log_result_memory!(run_name*" output", output_random)
                    avg_time_mp = mean(output_random["Time MP hist"])
                    gap = output_random["Gaps"]
                    push!(run_labels, run_name)
                    push!(gaps, gap)
                    results_destination = joinpath(results_folder, run_name)
                    df_cap, df_syscost, df_emissions = write_results_benders(output_random, inputs, settings, results_destination; budgets = budgets)
                    if mapping
                        temp_df_cap, temp_df_syscost, temp_df_emissions = write_mapping_results(output_random, inputs, settings)
                        push!(results_cap, temp_df_cap)
                        push!(results_syscost, temp_df_syscost)
                        push!(results_emissions, temp_df_emissions)
                    else
                        push!(results_cap, df_cap)
                        push!(results_syscost, df_syscost)
                        push!(results_emissions, df_emissions)
                    end
                    release_heavy_payload!(output_random)

                    if do_cap_wings
                        schedule = cap_schedule[budget_level_val]
                        for key in keys(schedule) # iterates through cvar and ev wings
                            caps = schedule[key]
                            isempty(caps) && continue
                            for (cap_level, cap) in enumerate(caps) # loose -> tight
                                if !haskey(budgets, key) # first cap of this wing: add the constraint
                                    budgets, constraint_dict = add_budget_constraint_bendersMP(MP, cap, key, budgets; constraint_dict = constraint_dict)
                                else
                                    budgets = update_budget_constraint_bendersMP!(MP, cap, key, budgets)
                                end

                                cap_run_name = "Random_$(iteration)_Budget_$(level)_$(key)_Cap_$(cap_level)"

                                try
                                    cap_output = mga_benders(inputs, settings, MP, SPs, budgets, cap_run_name; Eval_SPs = Eval_SPs, mapping = true, cut_archive = cut_archive)
                                    log_result_memory!(cap_run_name*" output", cap_output)
                                    avg_time_mp = mean(cap_output["Time MP hist"])
                                    push!(run_labels, cap_run_name)
                                    push!(gaps, cap_output["Gaps"])
                                    results_destination = joinpath(results_folder, cap_run_name)
                                    df_cap, df_syscost, df_emissions = write_results_benders(cap_output, inputs, settings, results_destination; budgets = budgets)
                                    if mapping
                                        temp_df_cap, temp_df_syscost, temp_df_emissions = write_mapping_results(cap_output, inputs, settings)
                                        push!(results_cap, temp_df_cap)
                                        push!(results_syscost, temp_df_syscost)
                                        push!(results_emissions, temp_df_emissions)
                                    else
                                        push!(results_cap, df_cap)
                                        push!(results_syscost, df_syscost)
                                        push!(results_emissions, df_emissions)
                                    end
                                    release_heavy_payload!(cap_output)
                                catch e
                                    @warn("Cap run $cap_run_name failed, skipping remaining caps for this wing", exception = (e, catch_backtrace()))
                                    break
                                end
                            end
                            # Remove this wing's constraint so it does not affect the other wing or later runs
                            con_ref = constraint_dict[key]
                            con_sym = Symbol(name(con_ref))
                            delete(MP, con_ref)
                            unregister(MP, con_sym)
                            delete!(budgets, key)
                            delete!(constraint_dict, key)
                        end
                    end
                catch e
                    @warn("Budget too tight (or error), terminating tightening", exception = (e, catch_backtrace()))
                    break
                end
            end
            if settings["Cut deactivation strategy"] == "in mga"
                if length(cuts_to_keep) < settings["Cuts retained"] && !mapping && avg_time_mp < 3*balanced_avg_time_mp
                    push!(cuts_to_keep, [name(con) for con in all_constraints(MP, include_variable_in_set_constraints=false) if startswith(string(con), "optimality_cut_") || startswith(string(con), "cvar_tail_cuts_")]...)
                end
            end

        end
        #write_gaps!(gaps, run_labels, joinpath(results_path, "Gaps"))

        # Map interior after exterior mapping
        if mapping
            all_caps = Matrix(vcat(results_cap...))
            all_costs = vcat(results_syscost...)
            @time samples = sample_interior_delauney(all_caps, all_costs, n_samples, settings, budget_val_transform, transform_sys, transform_cvar, scaling)
            outputs_mp = run_distributed_sampling(samples)
            @time for (i, sample) in enumerate(samples)
                outputs_sp = run_all_subproblems(SPs, inputs, settings, sample[1:R], sample[R+1:end]; minimal_payload=false)
                ev, cvar = evaluate_subproblems(outputs_sp, P_s, P_f, P_k, VaR_percent)
                @info("Sample $i: Investment cost = $(outputs_mp[i]["Inv_cost"]), Expected value = $ev, CVaR = $(cvar)")
                temp_df_cap, temp_df_syscost, temp_df_emissions = make_results_mapping_dfs(sample[1:R], sample[R+1:end], outputs_sp, outputs_mp[i]["Inv_cost"], outputs_mp[i]["Inv cost by zone"], cvar, ev, inputs, settings)
                push!(results_cap, temp_df_cap)
                push!(results_syscost, temp_df_syscost)
                push!(results_emissions, temp_df_emissions)
            end
        end
        write_exploration_results!(results_cap, results_syscost, results_emissions, summary_folder, run_labels, summary_name; mapping = mapping)
        return vectors
    end
end

function run_base_mga(SPs, new_inputs::Dict, settings::Dict, results_path::String, summary_folder::String; budget_multiplier::Float64 = 1.1, vector_set::Union{AbstractVector, Nothing} = nothing, scenario::Int = -1)
    configure_parallel_workers!(settings)

    # Result containers
    outputs = []
    labels = [] 
    results_cap = []
    results_syscost = []
    results_emissions = []
    gaps = []
    budget = 0.0

    # Model Settings
    iterations = settings["Iterations"]
    

    

    MP = build_planning_model(new_inputs, settings) #### When loaded with one scenario weighted, this is equivalent to a deterministic model with that scenario selected
    SP_one_scen = build_all_subproblems(new_inputs, settings)
    output = benders_algorithm(new_inputs, settings, MP, SP_one_scen, "OneScenarioLC"; Eval_SPs = SPs)
    log_result_memory!("OneScenarioLC output", output)
    gap = output["Gaps"]
    push!(labels, "OneScenarioLC")
    push!(gaps, gap)
    results_destination = joinpath(results_path,"CostOptimal")
    df_cap, df_syscost, df_emissions = write_results_benders(output, new_inputs, settings, results_destination)
    release_heavy_payload!(output)
    push!(results_cap, df_cap)
    push!(results_syscost, df_syscost)
    push!(results_emissions, df_emissions)
    lc_value = df_syscost[1, :EV_SystemCost]
    if budget_multiplier <= 10
        budget = (lc_value) * (budget_multiplier)
    end
    
    vectors = vector_set !== nothing ? vector_set : generate_weights(iterations, length(MP[:x])+length(MP[:x_line]), settings["Vector Type"], settings)
    @info("Using budget of ", budget, " for Base MGA test")
    percent_over_lc = round((budget - lc_value)/lc_value * 100, digits=2)
    @info("This budget is ", percent_over_lc, "% over the cost optimal solution")
    budgets = Dict()
    budgets, _ = add_budget_constraint_bendersMP(MP, budget/settings["Scaling factor cost"], "System_Expected", budgets)
    cuts_to_keep = [name(con) for (F, S) in list_of_constraint_types(MP) for con in all_constraints(MP, F, S) if startswith(string(con), "optimality_cut_") || startswith(string(con), "cvar_tail_cuts_")]

    for iteration in 1:iterations
        set_objective_bendersMP!(MP, "Capacity", new_inputs, settings; set_coeffs = vectors[iteration])
        output_random = mga_benders(new_inputs, settings, MP, SP_one_scen, budgets, "Base_MGA_"*string(iteration); Eval_SPs = SPs)
        log_result_memory!("Base_MGA_"*string(iteration)*" output", output_random)
        #cuts_to_keep = manage_cuts(MP, cuts_to_keep)
        gap = output_random["Gaps"]
        results_destination = joinpath(results_path,"Base_MGA_"*string(iteration))
        df_cap, df_syscost, df_emissions = write_results_benders(output_random, new_inputs, settings, results_destination; budgets = budgets)
        release_heavy_payload!(output_random)
        push!(results_cap, df_cap)
        push!(results_syscost, df_syscost)
        push!(results_emissions, df_emissions)
        push!(labels, "Random_"*string(iteration))
        push!(gaps, gap)
    end
    #write_gaps!(gaps, labels, joinpath(results_path, "Gaps"))
    write_exploration_results!(results_cap, results_syscost, results_emissions, summary_folder, labels, string(scenario))
end



#======================

Interior sampling functions


=======================#

function make_mgca_problem(points)
    rows, cols = size(points)
    points = max.(points, 0.0) # ensure all points are non-negative for convex combination

    model = Model(Ipopt.Optimizer)
    @variable(model, l[1:rows] >= 0)
    @constraint(model, c_sum, sum(l[i] for i in 1:rows) == 1)
    @variable(model, x[1:cols])
    @constraint(model, x == points' * l)
    @constraint(model, c_max_l, l .<= 0.95) # do not replicate individual points
    @objective(model, Min, 0)
    set_silent(model)
    return model
end

function sample_interior(points, num_samples, settings)
    return sample_interior_distributed(points, num_samples, settings)
end

# takes in points, cost_df, num samples and settings
# transforms cost_df into a 2d array of inv+ev and inv+cvar columns, with each row being a point in order of solve
# sends those to sample_interior_distributed which will return the evaluated mps

function sample_interior_delauney(points, cost_df::DataFrame, num_samples, settings, budget::Float64, transform_sys::Float64, transform_cvar::Float64, scaling::Float64)
    inv = cost_df[!,"Investment_Cost"] ./ scaling
    cvar = cost_df[!,"CVaR_OpCost"] ./ scaling
    ev = cost_df[!,"EV_OpCost"] ./ scaling

    # Drop solutions the Benders mapping recording pass let through despite violating the
    # Transformed budget constraint (it only checks the Benders convergence gap, not the
    # budget itself) - see algorithm.jl:432. Filtered against the loosest (untightened)
    # budget since results from every tightening level are combined here.
    transformed_cost = transform_sys .* (inv .+ ev) .+ transform_cvar .* (cvar)
    cutoff = budget + settings["Mapping Gap Threshold"]
    @info("Budget in sampling: $(budget); cutoff: $(cutoff)")
    keep = transformed_cost .<= cutoff

    n_dropped = length(keep) - sum(keep)
    if n_dropped > 0
        @info("Dropping $(n_dropped) of $(length(keep)) mapped solutions with transformed cost above budget + Mapping Gap Threshold ($(cutoff)) before Delaunay sampling.")
    end



    points = points[keep, :]
    inv, cvar, ev = inv[keep], cvar[keep], ev[keep]
    cost_points = Matrix([inv+ev inv+cvar])

    return sample_interior_simplex(points, cost_points, num_samples, settings)
end

#=================


Mapping test functions - run with mapping flag on and evaluation SPs to see how well the mapping performs in approximating the true performance of solutions across iterations


===================#



function mapping_test_laptop(test_index)

    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")#joinpath("inputs", "Inputs_30repdays_ext_1000scen_7techs")
    results_folder = joinpath("outputs", "Test_"*string(test_index), "Mapping_Test")
    summary_folder = joinpath(results_folder, "Summary")

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

    vector_set = nothing
    summary_name = "mapping"
    Eval_SPs = nothing
    mapping = true
    budget_multiplier = 1.001
    n_samples = 10

    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, results_folder, summary_folder; budget_multiplier = 1.001, vector_set = nothing, summary_name = "mapping", Eval_SPs = nothing, mapping = true, n_samples = 50)
    rmprocs(workers())

end

function mapping_test_della(test_index)

    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della")#joinpath("inputs", "Inputs_30repdays_ext_1000scen_7techs")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
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
    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, joinpath(results_folder, "Mapping_Test"), summary_folder; budget_multiplier = 1.001, vector_set = nothing, summary_name = "mapping", Eval_SPs = nothing, mapping=true, n_samples = settings["Interior Samples"])

end

function mapping_test_della_001(test_index)

    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della")#joinpath("inputs", "Inputs_30repdays_ext_1000scen_7techs")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
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
    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, joinpath(results_folder, "Mapping_Test"), summary_folder; budget_multiplier = 1.001, vector_set = nothing, summary_name = "mapping", Eval_SPs = nothing, mapping=true, n_samples = settings["Interior Samples"])

end

function risk_pareto_test_della(test_index, budget_type; budget_multiplier = 1.0)

    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della")#joinpath("inputs", "Inputs_30repdays_ext_1000scen_7techs")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
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
    results_folder = joinpath(results_folder, "Pareto5")
    #budget_multiplier = 1.00
    vector_set = nothing
    summary_name = "pareto"
    Eval_SPs = nothing
    mapping=false
    n_samples = settings["Interior Samples"]
    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, results_folder , summary_folder; budget_multiplier = budget_multiplier, vector_set = nothing, summary_name = "pareto", Eval_SPs = nothing, mapping=false, n_samples = settings["Interior Samples"], budget_type = budget_type)

end

function risk_pareto_test_della_deterministic(test_index, budget_type; budget_multiplier=1.0)

    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della")#joinpath("inputs", "Inputs_30repdays_ext_1000scen_7techs")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
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
    results_folder = joinpath(results_folder, "Pareto5")
    #budget_multiplier = 1.00
    vector_set = nothing
    summary_name = "pareto"
    Eval_SPs = nothing
    mapping=false
    n_samples = settings["Interior Samples"]
    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, results_folder , summary_folder; budget_multiplier = budget_multiplier, vector_set = nothing, summary_name = "pareto", Eval_SPs = nothing, mapping=false, n_samples = settings["Interior Samples"], budget_type = budget_type)

end

function risk_pareto_test_laptop_5(test_index)

    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")#joinpath("inputs", "Inputs_30repdays_ext_1000scen_7techs")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
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
    results_folder = joinpath(results_folder, "Pareto5")
    budget_multiplier = 0.00
    vector_set = nothing
    summary_name = "pareto"
    Eval_SPs = nothing
    mapping=true
    n_samples = settings["Interior Samples"]
    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, results_folder , summary_folder; budget_multiplier = budget_multiplier, vector_set = nothing, summary_name = "pareto", Eval_SPs = nothing, mapping=mapping, n_samples = settings["Interior Samples"], budget_type = "Transformed")

end

function risk_pareto_tighten_test_laptop(test_index)

    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
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
    results_folder = joinpath(results_folder, "Pareto5")
    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, results_folder, summary_folder; vector_set = nothing, summary_name = "pareto", Eval_SPs = nothing, mapping = true, n_samples = settings["Interior Samples"], budget_type = "Transformed", tighten_budget = true)

end

function risk_pareto_tighten_test_della(test_index)

    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
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
    results_folder = joinpath(results_folder, "Pareto5")
    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, results_folder, summary_folder; vector_set = nothing, summary_name = "pareto", Eval_SPs = nothing, mapping = true, n_samples = settings["Interior Samples"], budget_type = "Transformed", tighten_budget = true)

end

function risk_pareto_test_no_tighten_della(test_index)

    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
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
    results_folder = joinpath(results_folder, "Pareto5")
    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, results_folder, summary_folder; vector_set = nothing, summary_name = "pareto", Eval_SPs = nothing, mapping = false, n_samples = settings["Interior Samples"], budget_type = "Transformed", tighten_budget = false)

end

function wing_caps_test_della(test_index)

    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
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
    results_folder = joinpath(results_folder, "Pareto5")
    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, results_folder, summary_folder; vector_set = nothing, summary_name = "pareto", Eval_SPs = nothing, mapping = false, n_samples = settings["Interior Samples"], budget_type = "Transformed", tighten_budget = true, cap_wings = true)

end

# Minimal diagnostic run for comparing solver warm-start behavior before/after changes to
# Method/Crossover settings (see src/model/benders/master_planning.jl,
# subproblems_economic_dispatch.jl). Skips the full MGA/budget-tightening wrapper - one
# fixed objective, one straight-through benders_algorithm call - so it's fast to re-run
# after each settings change, unlike the full pareto/tighten harness above. Turns on
# "Warmstart debug flag" so each master/subproblem solve logs its simplex/barrier
# iteration count via _log_solver_work! (algorithm.jl) - a more direct signal that warm
# starting is reducing solver work than wall-clock time alone, which is noisy on a laptop.
function warmstart_verification_laptop(test_index)
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    if !isdir(results_folder)
        mkpath(results_folder)
    end
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)
    settings["Warmstart debug flag"] = true

    configure_parallel_workers!(settings)

    MP = build_planning_model(inputs, settings)
    set_objective_bendersMP!(MP, "System_Expected", inputs, settings)
    SPs = build_all_subproblems(inputs, settings)

    elapsed = @elapsed output = benders_algorithm(inputs, settings, MP, SPs, "Warmstart_verification")

    @info("Warmstart verification run complete in $(round(elapsed; digits=2)) seconds")
    @info("Final planning objective: $(output["MP"]["Planning objective"])")
    @info("Final gap: $(output["Gaps"][end])")
    @info("Time MP hist: $(output["Time MP hist"])")
    @info("Time SP hist: $(output["Time SP hist"])")
    @info("Time Reg hist: $(output["Time Reg hist"])")

    return output
end


function simple_comp(test_index)
    inputs_folder = joinpath("inputs","Inputs_30d_1000scen_7tech_2z")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    if !isdir(results_folder)
        mkpath(results_folder)
    end
    if !isdir(summary_folder)
        mkpath(summary_folder)
    end
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)
    
    inputs["Output Demand scenario probabilities"] = inputs["Demand scenario probabilities"] #Establish base weights ahead of time
    inputs["Output Gas price scenario probabilities"] = inputs["Gas price scenario probabilities"]
    inputs["Output Weather scenario probabilities"] = inputs["Weather scenario probabilities"]

    configure_parallel_workers!(settings)

    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)

    #outputs_mixed, vectors = run_stochastic_exploration(SPs, inputs, settings, joinpath(results_folder, "Both_Flipped"), summary_folder; budget_multiplier=1.10, vector_set=nothing)#vectors
    #outputs_exp, vectors = run_stochastic_exploration_single_type(SPs, inputs, settings, joinpath(results_folder, "Expected"), summary_folder; type = "System_Expected", vector_set = nothing, budget_multiplier=1.10)
    #outputs_cvar, vectors = run_stochastic_exploration_single_type(SPs, inputs, settings, joinpath(results_folder, "CVaR"), summary_folder; type ="System_Weighted_CVaR",vector_set = nothing, budget_multiplier=1.10)
    
    @info("Running base MGA test for mean scenario")
    
    # set all uncertainties to false and risk aversion to false for this test, since we are just running with one scenario selected
    settings["Risk aversion flag"] = false
    settings["Demand uncertainty"] = false
    settings["Gas price uncertainty"] = false
    settings["Weather uncertainty"] = false
    new_inputs = load_input_data(inputs_folder, settings)
    new_inputs["Full Demand scenario probabilities"] = inputs["Demand scenario probabilities"]
    new_inputs["Full Gas price scenario probabilities"] = inputs["Gas price scenario probabilities"]
    new_inputs["Full Weather scenario probabilities"] = inputs["Weather scenario probabilities"]

    i = 0
    results_path = joinpath(results_folder, "Results_Base_MGA", "Scenario_"*string(i))
    _ = run_base_mga(SPs, new_inputs, settings, results_path, summary_folder; budget_multiplier=1.20, vector_set=nothing, scenario=i)

end


###########

# Interpolation evaluation functions

#######

function evaluate_interpolates(SPs::Array{Model, 3}, inputs::Dict, settings::Dict, interpolate_df::DataFrame, results_folder::String)
    # result containers
    results_cap = []
    results_syscost = []
    results_emissions = []

    P_s = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]
    VaR_percent = settings["Value-at-Risk percent"]
    Z = inputs["Number of zones"]
    costs_by_resource = inputs["Investment costs"]
    cost_of_lines = inputs["CAPEX per MW"]
    n_lines = length(cost_of_lines)

    # isolate capacity vectors
    index_indx = findfirst(==("Index"), names(interpolate_df))
    labels = Vector(interpolate_df[:, index_indx])
    select!(interpolate_df, Not(:Index))
    investment_indx = findfirst(==("Investment_Cost"), names(interpolate_df))
    cap_vectors = interpolate_df[:, 1:investment_indx-n_lines-1]
    line_vectors = interpolate_df[:, investment_indx-n_lines:investment_indx-1]
    # get costs
    cap_inv_costs = cap_vectors .* reshape(costs_by_resource, 1, :)
    line_inv_costs = line_vectors .* reshape(cost_of_lines, 1, :)
    inv_costs = hcat(cap_inv_costs, line_inv_costs)
    tot_inv_costs = sum.(eachrow(Matrix(inv_costs)))
    
    col_by_zone = collect(findall(x -> occursin("z"*string(z), x), names(inv_costs)) for z in 1:Z)
    cost_by_zone = [sum.(eachrow(inv_costs[:, col_by_zone[z]])) for z in 1:Z]

    for i in 1:nrow(interpolate_df)
        @info("Evaluating interpolate ", labels[i])
        caps = Vector(cap_vectors[i, :])
        lines = Vector(line_vectors[i, :])
        outputs_sp = run_all_subproblems(SPs, inputs, settings, caps, lines; minimal_payload=false)
        ev, cvar = evaluate_subproblems(outputs_sp, P_s, P_f, P_k, VaR_percent)

        costs_by_zone_it = [cost_by_zone[z][i] for z in 1:Z]

        df_cap, df_syscost, df_emissions = make_results_mapping_dfs(caps, lines, outputs_sp, tot_inv_costs[i], costs_by_zone_it, cvar, ev, inputs, settings)
        push!(results_cap, df_cap)
        push!(results_syscost, df_syscost)
        push!(results_emissions, df_emissions)
    end
    write_exploration_results!(results_cap, results_syscost, results_emissions, results_folder, labels, "Interpolates"; mapping=true)
end

function evaluate_subproblems(outputs_sp::Array{Dict{String, Any}, 3}, P_s::Vector, P_f::Vector, P_k::Vector, VaR_percent::Float64)
    S = length(P_s)
    F = length(P_f)
    K = length(P_k)

    sp_obj_per_iter = reshape([outputs_sp[s,f,k]["SP objective"] for s in 1:S, f in 1:F, k in 1:K], (S, F, K))
    ev = sum(P_s[s]*P_f[f]*P_k[k]*sp_obj_per_iter[s,f,k] for s in 1:S, f in 1:F, k in 1:K)
    cvar = compute_cvar(sp_obj_per_iter, P_s, P_f, P_k, VaR_percent)
    return ev, cvar
end

#=========

Run interpolate evaluation

==========#

function run_interpolate_evaluation_laptop(test_index)
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")#joinpath("inputs", "Inputs_30repdays_ext_1000scen_7techs")
    results_folder = joinpath("outputs", "Test_"*string(test_index), "Interpolate_Evaluation")
    summary_folder = joinpath(results_folder, "Summary")
    interp_file = joinpath("experiments", "Della_experiments", "Mapping","mapping_budget_2p_MGCA.csv")

    if !isdir(results_folder)
        mkpath(results_folder)
    end
    if !isdir(summary_folder)
        mkpath(summary_folder)
    end
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)
    interpolate_df = CSV.read(interp_file, DataFrame, header=true)

    configure_parallel_workers!(settings)

    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)

    evaluate_interpolates(SPs, inputs, settings, interpolate_df, summary_folder)

end



function run_interpolate_evaluation_della(test_index, interp_file)
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della")#joinpath("inputs", "Inputs_30repdays_ext_1000scen_7techs")
    results_folder = joinpath("outputs", "Test_"*string(test_index), "Interpolate_Evaluation")
    summary_folder = joinpath(results_folder, "Summary")

    if !isdir(results_folder)
        mkpath(results_folder)
    end
    if !isdir(summary_folder)
        mkpath(summary_folder)
    end
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)
    interpolate_df = CSV.read(interp_file, DataFrame, header=true)

    configure_parallel_workers!(settings)

    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)

    evaluate_interpolates(SPs, inputs, settings, interpolate_df, summary_folder)

end

#=========

Cut management tests: exercises add_optimality_cuts!/deactivate_cuts/reactivate_cuts
directly against a real (laptop-sized) MP, without running a full Benders loop. Checks
that a deactivate -> reactivate round trip reproduces cuts exactly (same RHS), and times
three deletion patterns at production-scale batch sizes so a regression back to any of the
slower paths would show up as a large gap: (a) constraint_by_name-per-call + scalar delete
(the original pattern), (b) cut_refs-cache-per-call + scalar delete (isolates the cache win
alone), (c) cut_refs-cache + batched delete via deactivate_cuts itself (the production
function - isolates the batch-delete win on top of the cache).

==========#

function test_cut_management(test_index = "cut_mgmt")
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)

    MP = build_planning_model(inputs, settings)

    P_s = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]
    S, F, K = length(P_s), length(P_f), length(P_k)
    G = inputs["Number of generation resources"]
    O = inputs["Number of storage resources"]
    R = G + O
    L = inputs["Number of lines"]

    x_prev = zeros(R)
    x_prev_line = zeros(L)
    coeffs = fill(1.0 / (S*F*K), S, F, K)
    SP_obj = zeros(S, F, K)
    cap_dual = zeros(R, S, F, K)
    line_dual = zeros(L, S, F, K)

    cut_refs = Dict{String, ConstraintRef}()
    n_iterations = 6
    for it in 1:n_iterations
        new_refs = add_optimality_cuts!(MP, SP_obj, cap_dual, line_dual, x_prev, x_prev_line, coeffs, inputs, settings, it, "cuttest")
        merge!(cut_refs, new_refs)
    end

    total_cuts = length(cut_refs)
    @assert total_cuts == n_iterations * S * F * K "expected $(n_iterations*S*F*K) cuts, got $(total_cuts)"
    @info("Added $(total_cuts) cuts across $(n_iterations) iterations ($(S)x$(F)x$(K) scenarios/iteration).")

    # --- Correctness: deactivate + reactivate round trip ---
    cut_archive = Dict{String, Any}()
    all_names = collect(keys(cut_refs))
    to_deactivate = all_names[1:min(50, length(all_names))]
    rhs_before = Dict(n => normalized_rhs(cut_refs[n]) for n in to_deactivate)

    deactivate_cuts(MP, to_deactivate, cut_archive, cut_refs)
    @assert length(cut_archive) == length(to_deactivate)
    @assert all(n -> !haskey(cut_refs, n), to_deactivate)
    @assert all(n -> constraint_by_name(MP, n) === nothing, to_deactivate)
    @info("Deactivated $(length(to_deactivate)) cuts; archived and removed from MP as expected.")

    reactivated = reactivate_cuts(MP, cut_archive, to_deactivate)
    @assert isempty(cut_archive)
    @assert length(reactivated) == length(to_deactivate)
    merge!(cut_refs, reactivated)

    rhs_after = Dict(n => normalized_rhs(cut_refs[n]) for n in to_deactivate)
    @assert rhs_before == rhs_after "reactivated cuts' RHS did not match the originals"
    @info("Reactivated all $(length(to_deactivate)) cuts; RHS values match the originals exactly.")

    # --- Self-healing fallback: a cut live in MP but absent from cut_refs (simulates a
    # cut carried over from a previous problem on a reused MP) should still resolve. ---
    probe_name = all_names[end]
    delete!(cut_refs, probe_name)
    con = constraint_by_name(MP, probe_name)
    @assert con !== nothing
    cut_refs[probe_name] = con # mirrors the fallback populate step in algorithm.jl
    @info("Fallback lookup for an untracked-but-live cut resolved correctly.")

    # --- Performance: constraint_by_name/scalar (a) vs cut_refs-cache/scalar (b) vs
    # cut_refs-cache/batched deactivate_cuts (c), at a production-scale batch size. n_perf
    # is intentionally much larger than the old N=300: the batch-delete win comes from
    # avoiding an O(total constraints) reindex/resync per deleted cut in HiGHS.jl/Gurobi.jl,
    # which only shows up clearly once the batch size is in the thousands (matching real
    # cuts_to_remove sizes given `Cut deactivation threshold: 8`/`Cuts retained: 50000`).
    n_perf = min(3000, S*F*K*n_iterations)
    refill_iter = Ref(n_iterations)
    refill!(needed) = begin
        iters = ceil(Int, needed / (S*F*K)) + 1
        for _ in 1:iters
            refill_iter[] += 1
            new_refs = add_optimality_cuts!(MP, SP_obj, cap_dual, line_dual, x_prev, x_prev_line, coeffs, inputs, settings, refill_iter[], "cuttest_refill")
            merge!(cut_refs, new_refs)
        end
    end

    # (a) constraint_by_name-per-call + scalar delete (the original pattern)
    perf_batch_a = collect(keys(cut_refs))[1:n_perf]
    t_a = @elapsed begin
        for n in perf_batch_a
            con = constraint_by_name(MP, n)
            delete(MP, con)
        end
    end
    @info("(a) constraint_by_name + scalar delete: removed $(n_perf) cuts in $(round(t_a; digits=4))s")
    for n in perf_batch_a
        delete!(cut_refs, n)
    end

    # (b) cut_refs cache + scalar delete (isolates the cache win alone, no batching)
    refill!(n_perf)
    perf_batch_b = collect(keys(cut_refs))[1:n_perf]
    t_b = @elapsed begin
        for n in perf_batch_b
            delete(MP, cut_refs[n])
        end
    end
    @info("(b) cut_refs cache + scalar delete: removed $(n_perf) cuts in $(round(t_b; digits=4))s")
    for n in perf_batch_b
        delete!(cut_refs, n)
    end

    # (c) cut_refs cache + batched delete, via the production deactivate_cuts function
    refill!(n_perf)
    perf_batch_c = collect(keys(cut_refs))[1:n_perf]
    cut_archive_c = Dict{String, Any}()
    t_c = @elapsed deactivate_cuts(MP, perf_batch_c, cut_archive_c, cut_refs)
    @info("(c) cut_refs cache + batched delete (deactivate_cuts): removed $(n_perf) cuts in $(round(t_c; digits=4))s")

    speedup_cache = t_a / max(t_b, 1e-9)
    speedup_batch = t_b / max(t_c, 1e-9)
    speedup_total = t_a / max(t_c, 1e-9)
    @info("[$(test_index)] Speedup from cut_refs cache alone: $(round(speedup_cache; digits=1))x; from batching on top of the cache: $(round(speedup_batch; digits=1))x; combined: $(round(speedup_total; digits=1))x")

    return (n_perf = n_perf, scalar_by_name_time = t_a, scalar_cached_time = t_b, batched_time = t_c,
            speedup_cache = speedup_cache, speedup_batch = speedup_batch, speedup_total = speedup_total,
            total_cuts = total_cuts)

end

#=========

Capacity bound tightening test: build_planning_model previously bounded every resource's
capacity variable by an arbitrary global x_ub = 1e6, unrelated to the input data (every
resource's own real Capacity_UB in Resources.csv/Resources_storage.csv is far smaller,
~1e5-2e5 for the laptop inputs) and needlessly widened the master problem's RHS/bound
coefficient range. build_planning_model now bounds each x[r] by its own
inputs["Capacity upper bounds"][r] instead. This test checks that binding is actually
wired up correctly (JuMP's upper_bound matches the real per-resource data, not 1e6) and
that the master problem still solves normally under a realistic synthetic cut with the
tighter bounds in place.

Note: an earlier version of this fix tried to additionally rescale the x/x_line
variables themselves (dividing their values while multiplying every cost coefficient
that touches them back up to compensate) to bring their magnitude closer to the
already-scaled cost terms. That was reverted: compensating a coefficient by exactly the
factor used to shrink the variable is a wash for any row that also contains other,
differently-scaled terms (e.g. probability-weighted alpha/cvar terms in budget
constraints and cuts) - it does not shrink the row's coefficient range, and in some rows
widens it. It also introduced a real bug (the small "avoid exact zero" lower bound on x
got divided down to 1e-8, which is exactly what HiGHS's "excessively small column
bounds" warning was flagging). Bound-tightening against real data has none of those
failure modes since it never touches a matrix coefficient, only a box bound.

==========#

function test_capacity_bound_tightening(test_index = "cap_bound")
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)

    P_s = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]
    S, F, K = length(P_s), length(P_f), length(P_k)
    G = inputs["Number of generation resources"]
    O = inputs["Number of storage resources"]
    R = G + O
    L = inputs["Number of lines"]
    cost_inv = inputs["Investment costs"]
    x_ub_real = inputs["Capacity upper bounds"]

    MP = build_planning_model(inputs, settings)

    bounds = [normalized_rhs(MP[:max_capacity][r]) for r in 1:R]
    @assert bounds == x_ub_real "max_capacity bound does not match inputs[\"Capacity upper bounds\"]: $(bounds) vs $(x_ub_real)"
    @assert all(bounds .< 1e6) "Expected every per-resource bound to be tighter than the old arbitrary 1e6 constant, got $(bounds)"
    @info("[$(test_index)] Confirmed max_capacity bounds match real per-resource Capacity_UB data (max $(maximum(bounds)), previously a flat 1e6 for every resource).")

    Random.seed!(1234)
    # Synthetic cut shaped like a real Benders cut: capacity duals centered on investment
    # cost (economically meaningful rather than degenerate) plus a baseline dispatch cost,
    # so the master problem actually trades off inv_cost against alpha instead of just
    # collapsing to a bound.
    x_prev = zeros(R)
    x_prev_line = zeros(L)
    cap_dual = (cost_inv .* (0.9 .+ 0.2 .* rand(R))) .* ones(R, S, F, K)
    line_dual = 100 .* rand(L, S, F, K)
    SP_obj = 1e6 .* (0.9 .+ 0.2 .* rand(S, F, K))
    coeffs = fill(1.0, S, F, K)

    add_optimality_cuts!(MP, SP_obj, cap_dual, line_dual, x_prev, x_prev_line, coeffs, inputs, settings, 1, "boundtest")
    output = run_planning_model(MP, settings, 0.5)

    @assert all(output["Capacity"] .<= x_ub_real .+ 1e-6) "Solved capacity exceeds the tightened per-resource bound"
    @info("[$(test_index)] PASSED: master problem solves normally with data-driven per-resource bounds.")

    return (bounds = bounds, capacity = output["Capacity"])
end

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
    settings["Sample Method"] = "HitAndRun"
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
    settings["Delaunay Repeat Mode"] = "CycleSampling"
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
    settings["Delaunay Repeat Mode"] = "ProgressiveRefinement"
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