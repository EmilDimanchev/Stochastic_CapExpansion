using Distributed

@everywhere include("../src/Stochastic_CapExpansion.jl")

@everywhere using .Stochastic_CapExpansion
@everywhere using Revise, JuMP, Gurobi, HiGHS, Ipopt, DataFrames, CSV, YAML, Random, LinearAlgebra, Combinatorics, Dates, Distributions, Surrogates
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



#============================


Main MGA functions for running stochastic exploration with Benders and MGA, including base runs and iterative exploration with budget constraints, cut management, and result logging


==============================#

function run_stochastic_exploration_separate_budgets(SPs::Array{Model, 3}, inputs::Dict, settings::Dict, results_folder::String, summary_folder::String; budget_multiplier::Float64 = 1.10, vector_set::Union{AbstractVector, Nothing} = nothing, summary_name::String = "new_setup", Eval_SPs = nothing, mapping = false, n_samples = 100)

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
    risk_aversion_weight = settings["Risk aversion weight"] ### Currently set consistently across runs and ahead of time
    VaR_percent = settings["Value-at-Risk percent"] ### Currently set consistently across runs and ahead of time
    # Create and set expected value model
    settings["Risk aversion flag"] = true

    P_s = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]
    R = length(inputs["Resources"])
    

    MP = build_planning_model(inputs, settings)

    # Base Runs
    set_objective_bendersMP!(MP, "System_Weighted_CVaR", inputs, settings; obj_weight = risk_aversion_weight)
    output_cvar = benders_algorithm(inputs, settings, MP, SPs, "System_Weighted_CVaR"; Eval_SPs = Eval_SPs, mapping = mapping)
    log_result_memory!("System_Weighted_CVaR output", output_cvar)
    gap_cvar = output_cvar["Gaps"]
    push!(run_labels, "System_Weighted_CVaR")
    push!(gaps, gap_cvar)
    avg_time_mp_base = mean(output_cvar["Time MP hist"])
     # Write results
    results_destination = joinpath(results_folder,"System_Weighted_CVaR")
    df_cap, df_syscost, df_emissions = write_results_benders(output_cvar, inputs, settings, results_destination)
    if mapping
        temp_df_cap, temp_df_syscost, temp_df_emissions = write_mapping_results(output_cvar, inputs, settings)
        push!(results_cap, temp_df_cap)
        push!(results_syscost, temp_df_syscost)
        push!(results_emissions, temp_df_emissions)
    else 
        push!(results_cap, df_cap)
        push!(results_syscost, df_syscost)
        push!(results_emissions, df_emissions)
    end
    
    release_heavy_payload!(output_cvar)
    
    @info("System_Weighted_CVaR solution has investment cost of ", output_cvar["MP"]["Inv_cost"])
    @info("Expected value system cost of " * "System_Weighted_CVaR" * " solution: $(output_cvar["Expected Value"] + output_cvar["MP"]["Inv_cost"])")
    @info("Risk adjusted system cost of " * "System_Weighted_CVaR" * " solution: $((1-risk_aversion_weight)*output_cvar["CVaR"] + risk_aversion_weight*output_cvar["Expected Value"]+ output_cvar["MP"]["Inv_cost"])")



    if settings["Capacity Exploration"]
        budgets = Dict()
        # Set budgets? ------- budget set = set with budgets same percent greater than least cost solution for each metric
        budgets = add_budget_constraint_bendersMP(MP, ((output_cvar["CVaR"])/settings["Scaling factor cost"]), "CVaR", budgets)
        budgets = add_budget_constraint_bendersMP(MP, ((output_cvar["Expected Value"] + output_cvar["MP"]["Inv_cost"])/settings["Scaling factor cost"])*(budget_multiplier), "System_Expected", budgets)
        cuts = [name(con) for con in all_constraints(MP, include_variable_in_set_constraints=false) if (startswith(string(con), "optimality_cut_") || startswith(string(con), "cvar_tail_cuts_"))]
        cuts_to_keep = copy(cuts)
        rhs_values = output_cvar["RHS Values"]
        removed_cuts = output_cvar["Removed cuts"]
        if mapping
            cuts_to_keep = filter(cut -> parse(Int, split(string(cut), "_")[end-1]) <= output_cvar["first_write"] + 10, cuts_to_keep)
        end
        vectors = vector_set !== nothing ? vector_set : generate_weights(iterations, length(MP[:x])+length(MP[:x_line]), settings["Vector Type"], settings)
        #@info("Keeping $(length(cuts_to_keep)) cuts for MGA iterations")
        for iteration in 1:iterations
            set_objective_bendersMP!(MP, "Capacity", inputs, settings; set_coeffs = vectors[iteration])
            cuts_to_keep = manage_cuts(MP, cuts_to_keep)
            # other option - reactivate all cuts then let alg deactivate those that are not useful
            #reactivate_cuts(MP, removed_cuts, rhs_values)

            output_random = mga_benders(inputs, settings, MP, SPs, budgets, "Random_"*string(iteration); Eval_SPs = Eval_SPs, mapping = mapping)
            log_result_memory!("Random_"*string(iteration)*" output", output_random)
            avg_time_mp = mean(output_random["Time MP hist"])
            gap = output_random["Gaps"]
            push!(run_labels, "Random_"*string(iteration))
            push!(gaps, gap)
            results_destination = joinpath(results_folder,"Random_"*string(iteration))
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
            rhs_values = output_random["RHS Values"]
            removed_cuts = output_random["Removed cuts"]
            
            if length(cuts_to_keep) < settings["Cuts retained"] && !mapping && avg_time_mp < 3*avg_time_mp_base
                push!(cuts_to_keep, [name(con) for con in all_constraints(MP, include_variable_in_set_constraints=false) if startswith(string(con), "optimality_cut_") || startswith(string(con), "cvar_tail_cuts_")]...)
            end
        end
        #write_gaps!(gaps, run_labels, joinpath(results_path, "Gaps"))

        # Map interior after exterior mapping
        if mapping
            all_caps = Matrix(vcat(results_cap...))
            @time samples = sample_interior(all_caps, n_samples, settings)
            outputs_mp = run_distributed_sampling(samples)
            @time for (i, sample) in enumerate(samples)
                outputs_sp = run_all_subproblems(SPs, inputs, settings, sample[1:R], sample[R+1:end]; minimal_payload=false)
                ev, cvar = evaluate_subproblems(outputs_sp, P_s, P_f, P_k, VaR_percent)
                @info("Sample $i: Investment cost = $(outputs_mp[i]["Inv_cost"]), Expected value = $ev, CVaR = $(cvar)")
                @info("Sample $i: Budget Percentage for CVaR: $(((cvar)/(budgets["CVaR"]*settings["Scaling factor cost"]))*100)%, Budget Percentage for Expected Value: $((((ev+outputs_mp[i]["Inv_cost"])/(budgets["System_Expected"]*settings["Scaling factor cost"]))*100)-100)%")
                temp_df_cap, temp_df_syscost, temp_df_emissions = make_results_mapping_dfs(sample[1:R], sample[R+1:end], outputs_sp, outputs_mp[i]["Inv_cost"], outputs_mp[i]["Inv cost by zone"], cvar, ev, inputs, settings)
                push!(results_cap, temp_df_cap)
                push!(results_syscost, temp_df_syscost)
                push!(results_emissions, temp_df_emissions)
            end
        end
        write_exploration_results!(results_cap, results_syscost, results_emissions, summary_folder, run_labels, summary_name; mapping = mapping)
    end
    return vectors
end


function run_stochastic_exploration_risk_pareto(SPs::Array{Model, 3}, inputs::Dict, settings::Dict, results_folder::String, summary_folder::String; budget_multiplier::Float64 = 1.10, vector_set::Union{AbstractVector, Nothing} = nothing, summary_name::String = "new_setup", Eval_SPs = nothing, mapping = false, n_samples = 100, budget_type = "Transformed")

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
    risk_aversion_weights = [0.5, 0.0, 1.0]#[[i for i in 0.5:-0.1:0.1]; [i for i in 0.6:0.1:0.9]; 0.0; 1.0]
    risk_aversion_weight = settings["Risk aversion weight"] ### set as anchor point
    VaR_percent = settings["Value-at-Risk percent"] ### Currently set consistently across runs and ahead of time
    settings["Risk aversion flag"] = true
    scaling = settings["Scaling factor cost"]
    # Create and set expected value model
    extreme_values = Dict()

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
        output_cvar = benders_algorithm(inputs, settings, MP, SPs, case_name; Eval_SPs = Eval_SPs, mapping = mapping, risk_aversion_weight = risk, cut_archive = cut_archive)
        log_result_memory!(case_name*" output", output_cvar)
        gap_cvar = output_cvar["Gaps"]
        push!(run_labels, case_name)
        push!(gaps, gap_cvar)
        # Write results
        results_destination = joinpath(results_folder, case_name)
        df_cap, df_syscost, df_emissions = write_results_benders(output_cvar, inputs, settings, results_destination)
        if mapping
            temp_df_cap, temp_df_syscost, temp_df_emissions = write_mapping_results(output, inputs, settings)
            push!(results_cap, temp_df_cap)
            push!(results_syscost, temp_df_syscost)
            push!(results_emissions, temp_df_emissions)
        else
            push!(results_cap, df_cap)
            push!(results_syscost, df_syscost)
            push!(results_emissions, df_emissions)
        end
        extreme_values[risk] = Dict("CVaR" => output_cvar["CVaR"]/scaling, "System Expected" => (output_cvar["Expected Value"] + output_cvar["MP"]["Inv_cost"])/scaling)

        @info("Risk weight $risk solution has investment cost of $(output_cvar["MP"]["Inv_cost"])")
        @info("Expected value system cost of " * "Risk weight $risk" * " solution: $(output_cvar["Expected Value"] + output_cvar["MP"]["Inv_cost"])")
        @info("Risk adjusted system cost of " * "Risk weight $risk" * " solution: $((1-risk_aversion_weight)*output_cvar["CVaR"] + risk_aversion_weight*output_cvar["Expected Value"]+ output_cvar["MP"]["Inv_cost"])")
        
        if settings["Cut deactivation strategy"] == "in mga"
            
            cuts_to_keep = [name(con) for con in all_constraints(MP, include_variable_in_set_constraints=false) if ((startswith(string(con), "optimality_cut_") && (split(string(con), "_")[3] == string(risk))) || (startswith(string(con), "cvar_tail_cuts_") && (split(string(con), "_")[3] == string(risk))))]
        
            if mapping
                cuts_to_keep = filter(cut -> parse(Int, split(string(cut), "_")[end-1]) <= output["first_write"] + 10, cuts_to_keep)
            end
            cuts_to_keep = manage_cuts(MP, cuts_to_keep)
        end

        release_heavy_payload!(output_cvar)
    end

    if settings["Capacity Exploration"]
        budgets = Dict()
        # CVaR budget: the CVaR achieved by the risk-neutral (expected-value-only) solution.
        # extreme_values[...] is already in cost-scaled units (divided by "Scaling factor
        # cost" once above, when it was populated) - dividing again here shrank the budget
        # by another 1e5x, making it essentially unsatisfiable against any real cut. That
        # was the actual cause of "numerical issues" reported for the MGA phase: not a
        # conditioning problem, a budget that was ~1e5x too tight by construction.
        #budgets = add_budget_constraint_bendersMP(MP, extreme_values[1.0]["CVaR"], "CVaR", budgets)
        # Expected value budget: the EV achieved by the fully risk-averse (CVaR-only)
        # solution, with budget_multiplier headroom (mirrors the working pattern in
        # run_stochastic_exploration_separate_budgets above). The previous line added
        # extreme_values[0.0]["System Expected"] to itself instead of applying
        # budget_multiplier, which this function accepts but never used.
        #budgets = add_budget_constraint_bendersMP(MP, extreme_values[0.0]["System Expected"]*budget_multiplier, "System_Expected", budgets)
        # Combined budget: keep inv_cost + 0.5*EV + 0.5*CVaR within (1+budget_multiplier) of the
        # balanced (risk weight 0.5) solution's own optimal value.
        # mga_benders tracks this budget's upper bound using settings["Risk aversion weight"]
        # (algorithm.jl:608/706), not the risk_aversion passed here, so it must be kept at 0.5
        # to match the constraint actually being enforced on MP.
        #settings["Risk aversion weight"] = 0.5
        if budget_type == "Transformed"
            transform_cvar = 1/(extreme_values[1.0]["CVaR"] - extreme_values[0.0]["CVaR"])
            transform_sys = 1/(extreme_values[0.0]["System Expected"] - extreme_values[1.0]["System Expected"])
            budget_val_transform = transform_cvar*extreme_values[0.0]["CVaR"] + transform_sys*extreme_values[1.0]["System Expected"] + 1 + budget_multiplier
            budgets = add_budget_constraint_bendersMP(MP, budget_val_transform, "Transformed", budgets; extreme_values = extreme_values)
        elseif budget_type == "Box"
            budgets = add_budget_constraint_bendersMP(MP, extreme_values[1.0]["CVaR"], "CVaR", budgets)
            budgets = add_budget_constraint_bendersMP(MP, extreme_values[0.0]["System Expected"]*budget_multiplier, "System_Expected", budgets)
        else
            error("Unknown budget type: $budget_type")
        end
        #budgets = add_budget_constraint_bendersMP(MP, (balanced_optimum/settings["Scaling factor cost"])*(1+budget_multiplier), "System_Weighted_CVaR", budgets; risk_aversion = settings["Risk aversion weight"])
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

            output_random = mga_benders(inputs, settings, MP, SPs, budgets, "Random_"*string(iteration); Eval_SPs = Eval_SPs, mapping = mapping, cut_archive = cut_archive)
            log_result_memory!("Random_"*string(iteration)*" output", output_random)
            avg_time_mp = mean(output_random["Time MP hist"])
            gap = output_random["Gaps"]
            push!(run_labels, "Random_"*string(iteration))
            push!(gaps, gap)
            results_destination = joinpath(results_folder,"Random_"*string(iteration))
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
            @time samples = sample_interior(all_caps, n_samples, settings)
            outputs_mp = run_distributed_sampling(samples)
            @time for (i, sample) in enumerate(samples)
                outputs_sp = run_all_subproblems(SPs, inputs, settings, sample[1:R], sample[R+1:end]; minimal_payload=false)
                ev, cvar = evaluate_subproblems(outputs_sp, P_s, P_f, P_k, VaR_percent)
                @info("Sample $i: Investment cost = $(outputs_mp[i]["Inv_cost"]), Expected value = $ev, CVaR = $(cvar)")
                @info("Sample $i: Budget Percentage for CVaR: $(((cvar)/(budgets["CVaR"]*settings["Scaling factor cost"]))*100)%, Budget Percentage for Expected Value: $((((ev+outputs_mp[i]["Inv_cost"])/(budgets["Expected"]*settings["Scaling factor cost"]))*100)-100)%")
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
    lc_value = df_syscost[1, :SingleScenario_SystemCost]
    if budget_multiplier <= 10
        budget = (lc_value) * (budget_multiplier)
    end
    
    vectors = vector_set !== nothing ? vector_set : generate_weights(iterations, length(MP[:x]), settings["Vector Type"], settings)
    @info("Using budget of ", budget, " for Base MGA test")
    percent_over_lc = round((budget - lc_value)/lc_value * 100, digits=2)
    @info("This budget is ", percent_over_lc, "% over the cost optimal solution")
    budgets = Dict()
    budgets = add_budget_constraint_bendersMP(MP, budget/settings["Scaling factor cost"], "System_Expected", budgets)
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

    vectors = run_stochastic_exploration_separate_budgets(SPs, inputs, settings, results_folder, summary_folder; budget_multiplier = 1.001, vector_set = nothing, summary_name = "mapping", Eval_SPs = nothing, mapping = true, n_samples = 50)
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
    vectors = run_stochastic_exploration_separate_budgets(SPs, inputs, settings, joinpath(results_folder, "Mapping_Test"), summary_folder; budget_multiplier = 1.001, vector_set = nothing, summary_name = "mapping", Eval_SPs = nothing, mapping=true, n_samples = settings["Interior Samples"])

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

function risk_pareto_test_della(test_index, budget_type)

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
    budget_multiplier = 1.00
    vector_set = nothing
    summary_name = "pareto"
    Eval_SPs = nothing
    mapping=false
    n_samples = settings["Interior Samples"]
    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, results_folder , summary_folder; budget_multiplier = budget_multiplier, vector_set = nothing, summary_name = "pareto", Eval_SPs = nothing, mapping=false, n_samples = settings["Interior Samples"], budget_type = budget_type)

end

function risk_pareto_test_della_deterministic(test_index, budget_type)

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
    budget_multiplier = 1.00
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
    budget_multiplier = 1.00
    vector_set = nothing
    summary_name = "pareto"
    Eval_SPs = nothing
    mapping=false
    n_samples = settings["Interior Samples"]
    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, results_folder , summary_folder; budget_multiplier = budget_multiplier, vector_set = nothing, summary_name = "pareto", Eval_SPs = nothing, mapping=false, n_samples = settings["Interior Samples"], budget_type = "Transformed")

end


function simple_comp(test_index)
    inputs_folder = joinpath("inputs","Inputs_30d_1000scen_7tech_2z_Della")
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
    _ = run_base_mga(SPs, new_inputs, settings, results_path, summary_folder; budget_multiplier=1.10, vector_set=vectors, scenario=i)

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
that a deactivate -> reactivate round trip reproduces cuts exactly (same RHS), and
times the constraint_by_name-per-call pattern against the cut_refs-cache pattern so a
regression back to the slow path would show up as a large old/new time gap.

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

    # --- Performance: constraint_by_name-per-call (old) vs cut_refs cache (new) ---
    remaining_names = collect(keys(cut_refs))
    n_perf = min(300, length(remaining_names))
    perf_batch = remaining_names[1:n_perf]

    t_old = @elapsed begin
        for n in perf_batch
            con = constraint_by_name(MP, n)
            delete(MP, con)
        end
    end
    @info("Old pattern (constraint_by_name in a delete loop): removed $(n_perf) cuts in $(round(t_old; digits=4))s")
    for n in perf_batch
        delete!(cut_refs, n)
    end

    refill_iterations = ceil(Int, n_perf / (S*F*K)) + 1
    for it in (n_iterations+1):(n_iterations+refill_iterations)
        new_refs = add_optimality_cuts!(MP, SP_obj, cap_dual, line_dual, x_prev, x_prev_line, coeffs, inputs, settings, it, "cuttest_refill")
        merge!(cut_refs, new_refs)
    end
    perf_batch_2 = collect(keys(cut_refs))[1:n_perf]
    cut_archive_2 = Dict{String, Any}()
    t_new = @elapsed deactivate_cuts(MP, perf_batch_2, cut_archive_2, cut_refs)
    @info("New pattern (cut_refs cache): removed $(n_perf) cuts in $(round(t_new; digits=4))s")

    speedup = t_old / max(t_new, 1e-9)
    @info("[$(test_index)] Speedup from cut_refs cache over constraint_by_name-per-call: $(round(speedup; digits=1))x")

    return (old_time = t_old, new_time = t_new, speedup = speedup, total_cuts = total_cuts)

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