#============================


Main MGA functions for running stochastic exploration with Benders and MGA, including base runs and iterative exploration with budget constraints, cut management, and result logging


==============================#

const SAMPLING_METHODS = ("None", "Delaunay", "Boundary", "GridFill")

# Each *_settings keyword overrides keys of the matching settings.yaml block (see run_settings.jl).
function run_stochastic_exploration_risk_pareto(SPs::Array{Model, 3}, inputs::Dict, settings::Dict, results_folder::String, summary_folder::String; budget_settings::AbstractDict = Dict(), cap_settings::AbstractDict = Dict(), mapping_settings::AbstractDict = Dict(), run_settings::AbstractDict = Dict(), sampling_settings::AbstractDict = Dict(), cut_settings::AbstractDict = Dict())

    settings = resolve_all_settings(settings, Dict("Budget settings" => budget_settings, "Cap settings" => cap_settings, "Mapping settings" => mapping_settings,
                                                   "Run settings" => run_settings, "Sampling settings" => sampling_settings, "Cut settings" => cut_settings))
    budget_settings = settings["Budget settings"]
    sampling_settings = settings["Sampling settings"]
    cap_settings = settings["Cap settings"]
    mapping_settings = settings["Mapping settings"]
    run_settings = settings["Run settings"]

    budget_type = budget_settings["Budget type"]
    budget_offset = budget_settings["Budget offset"]
    budget_multiplier = budget_settings["Budget multiplier"]
    tighten_budget = budget_settings["Tighten budget"]
    budget_tighten_step_size = budget_settings["Tighten step size"]
    floor_offset = budget_settings["Floor offset"]

    cap_wings = cap_settings["Cap wings"]
    n_caps = cap_settings["Number of caps"]
    cap_levels = cap_settings["Cap levels"]

    mapping = mapping_settings["Mapping"]
    if mapping
        @info("Mapping on")
    end
    sampling_method = sampling_settings["Method"]
    sampling_method in SAMPLING_METHODS || error("Unknown sampling \"Method\" $(sampling_method); use one of $(SAMPLING_METHODS).")
    n_samples = sampling_settings["Samples"]

    vector_set = run_settings["Vector set"]
    summary_name = run_settings["Summary name"]
    Eval_SPs = run_settings["Eval SPs"]
    run_log = run_settings["Run log"]

    #configure_parallel_workers!(settings)

    # Optional per-run diagnostics (used by test_cut_trim_rules / test_cut_pool_ab): SP rounds
    # (= Benders iterations), MP solves including separation re-solves, live cuts at the first
    # MP solve, MP and separation time, and the pool size (live + archived) before/after.
    log_run!(label, out, archive_before) = run_log === nothing || begin
        sep = out["Sep hist"]
        entry = (label = label, iterations = length(out["Gaps"]), sp_rounds = length(out["Time SP hist"]),
                 mp_solves = length(sep) + sum(s.rounds for s in sep), live_at_start = out["Live cuts at start"],
                 mp_time = sum(out["Time MP hist"]), sep_time = sum(s.t_check for s in sep),
                 archive_before = archive_before, archive_after = length(cut_pool.live) + n_archived(cut_pool))
        push!(run_log, entry)
        @info("TRIMRUN|$(entry.label)|iters=$(entry.iterations)|sp_rounds=$(entry.sp_rounds)|mp_solves=$(entry.mp_solves)|live_at_start=$(entry.live_at_start)|mp_time=$(round(entry.mp_time; digits=3))|sep_time=$(round(entry.sep_time; digits=3))|cuts_before=$(entry.archive_before)|archive_after=$(entry.archive_after)")
    end
    pool_size() = length(cut_pool.live) + n_archived(cut_pool)

    # Result containers
    results_cap = []
    results_syscost = []
    results_emissions = []
    run_labels = []
    gaps = []

    # Model Settings
    iterations = run_settings["MGA vectors"]
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

    # Shared across every risk weight / MGA direction solved against this MP, so cuts built
    # for one problem can be brought back (seed!/separate!) when a later problem needs them.
    cut_pool = CutPool(settings)

    # Gradient bank for the Boundary / GridFill samplers: one FinalSolve (capacity, per-scenario
    # op cost and capacity gradients) per final solve. The converged subproblem outputs in
    # out["SPs"] already hold these duals, so this copies them before release_heavy_payload!
    # drops them.
    bank = FinalSolve[]
    use_bank = sampling_method in ("Boundary", "GridFill")
    var_names = capacity_variable_names(inputs)
    bank_record!(label, out, folder) = use_bank && begin
        fs = record_final_solve!(bank, label, out, inputs)
        sampling_settings["Write gradients"] && write_gradient_file(folder, fs, var_names)
    end

    for (i, risk) in enumerate(risk_aversion_weights)
        # Base Runs
        case_name = "Risk_Weight_"*string(risk)
        set_objective_bendersMP!(MP, "System_Weighted_CVaR", inputs, settings; obj_weight = risk)
        trim_cuts!(MP, cut_pool)
        archive_before = pool_size()
        output = benders_algorithm(inputs, settings, MP, SPs, case_name; Eval_SPs = Eval_SPs, mapping = mapping, risk_aversion_weight = risk, cut_pool = cut_pool, solve_type = "Anchor", mga_vector = 0)
        log_run!(case_name, output, archive_before)
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


        bank_record!(case_name, output, results_destination)
        release_heavy_payload!(output)
    end

    if run_settings["Capacity exploration"]
        checkpoint_time = 0.0 # cumulative seconds spent on incremental summary writes
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
                    cap_schedule = compute_budget_schedule_wings(extreme_values, budget_schedule; k_ev = transform_sys, k_cvar = transform_cvar, n_max = n_caps, )
                end
            end
        elseif budget_type == "Box"
            budgets, constraint_dict = add_budget_constraint_bendersMP(MP, extreme_values[1.0]["CVaR"], "CVaR", budgets; constraint_dict = constraint_dict)
            budgets, constraint_dict = add_budget_constraint_bendersMP(MP, extreme_values[0.0]["System Expected"]*budget_multiplier, "System_Expected", budgets; constraint_dict = constraint_dict)
        else
            error("Unknown budget type: $budget_type")
        end

        vectors = vector_set !== nothing ? vector_set : generate_weights(iterations, length(MP[:x])+length(MP[:x_line]), run_settings["Vector type"], settings)
        for iteration in 1:iterations
            set_objective_bendersMP!(MP, "Capacity", inputs, settings; set_coeffs = vectors[iteration])
            # Cuts (live + archived) accumulate freely across this vector's runs; trim the pool once per new vector.
            trim_cuts!(MP, cut_pool)

            levels = tighten_budget ? budget_schedule : [budget_val_transform]
            local output_random
            for (level, budget_level_val) in enumerate(levels)
                if tighten_budget
                    budgets = update_budget_constraint_bendersMP!(MP, budget_level_val, "Transformed", budgets)
                end
                run_name = tighten_budget ? "Random_$(iteration)_Budget_$(level)" : "Random_"*string(iteration)
                try
                    archive_before = pool_size()
                    output_random = mga_benders(inputs, settings, MP, SPs, budgets, run_name; Eval_SPs = Eval_SPs, mapping = mapping, cut_pool = cut_pool, solve_type = "Budget", mga_vector = iteration)
                    log_run!(run_name, output_random, archive_before)
                    log_result_memory!(run_name*" output", output_random)
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
                    bank_record!(run_name, output_random, results_destination)
                    release_heavy_payload!(output_random)

                    if do_cap_wings && level in cap_levels # only carry out on budgets asked for
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
                                    archive_before = pool_size()
                                    cap_output = mga_benders(inputs, settings, MP, SPs, budgets, cap_run_name; Eval_SPs = Eval_SPs, mapping = mapping, cut_pool = cut_pool, solve_type = key == "CVaR" ? "CVaR cap" : "EV cap", mga_vector = iteration)
                                    log_run!(cap_run_name, cap_output, archive_before)
                                    log_result_memory!(cap_run_name*" output", cap_output)
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
                                    bank_record!(cap_run_name, cap_output, results_destination)
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
                # Checkpoint the summary after every budget level so a late failure does not lose the run
                try
                    if !isempty(results_cap)
                        t_ckpt = @elapsed write_exploration_results!(results_cap, results_syscost, results_emissions, summary_folder, run_labels, summary_name; mapping = mapping)
                        checkpoint_time += t_ckpt
                        @info("Summary checkpoint after $run_name: $(round(t_ckpt; digits=3)) s ($(sum(nrow, results_cap)) rows), cumulative $(round(checkpoint_time; digits=3)) s")
                    end
                catch e
                    @warn("Summary checkpoint failed", exception = (e, catch_backtrace()))
                end
            end
        end
        #write_gaps!(gaps, run_labels, joinpath(results_path, "Gaps"))

        # Map interior after exterior mapping
        if sampling_method == "Boundary"
            # Interior_k rows carry their own labels plus lineage/bound columns; their
            # bracket check is written to Summary/Interior_brackets_<summary name>.csv
            @time rows_cap, rows_syscost, rows_emissions, interior_labels, _, _ = run_boundary_sampling!(bank, SPs, inputs, settings, sampling_settings, results_folder, summary_folder, summary_name)
            append!(results_cap, rows_cap)
            append!(results_syscost, rows_syscost)
            append!(results_emissions, rows_emissions)
            append!(run_labels, interior_labels)
        elseif sampling_method == "GridFill"
            # GridFill_k rows carry lineage + prediction columns (GF_*). The Summary is rewritten
            # after every batch with the MGA rows plus the samples solved so far.
            n_mga = length(results_cap)
            checkpoint(c, s, e, l) = write_exploration_results!(vcat(results_cap[1:n_mga], c), vcat(results_syscost[1:n_mga], s), vcat(results_emissions[1:n_mga], e),
                                                                summary_folder, vcat(run_labels, l), summary_name; mapping = mapping)
            @time rows_cap, rows_syscost, rows_emissions, gf_labels = run_gridfill_sampling!(bank, cut_pool, MP, SPs, inputs, settings, sampling_settings; checkpoint = checkpoint, results_folder = results_folder)
            append!(results_cap, rows_cap)
            append!(results_syscost, rows_syscost)
            append!(results_emissions, rows_emissions)
            append!(run_labels, gf_labels)
        elseif sampling_method == "Delaunay"
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
        t_final = @elapsed write_exploration_results!(results_cap, results_syscost, results_emissions, summary_folder, run_labels, summary_name; mapping = mapping)
        @info("Final summary write: $(round(t_final; digits=3)) s; total checkpoint time $(round(checkpoint_time; digits=3)) s")
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
    iterations = settings["Run settings"]["MGA vectors"]
    

    

    MP = build_planning_model(new_inputs, settings) #### When loaded with one scenario weighted, this is equivalent to a deterministic model with that scenario selected
    SP_one_scen = build_all_subproblems(new_inputs, settings)
    cut_pool = CutPool(settings)
    output = benders_algorithm(new_inputs, settings, MP, SP_one_scen, "OneScenarioLC"; Eval_SPs = SPs, cut_pool = cut_pool)
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
    
    vectors = vector_set !== nothing ? vector_set : generate_weights(iterations, length(MP[:x])+length(MP[:x_line]), settings["Run settings"]["Vector type"], settings)
    @info("Using budget of ", budget, " for Base MGA test")
    percent_over_lc = round((budget - lc_value)/lc_value * 100, digits=2)
    @info("This budget is ", percent_over_lc, "% over the cost optimal solution")
    budgets = Dict()
    budgets, _ = add_budget_constraint_bendersMP(MP, budget/settings["Scaling factor cost"], "System_Expected", budgets)

    for iteration in 1:iterations
        set_objective_bendersMP!(MP, "Capacity", new_inputs, settings; set_coeffs = vectors[iteration])
        output_random = mga_benders(new_inputs, settings, MP, SP_one_scen, budgets, "Base_MGA_"*string(iteration); Eval_SPs = SPs, cut_pool = cut_pool, mga_vector = iteration)
        log_result_memory!("Base_MGA_"*string(iteration)*" output", output_random)
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
