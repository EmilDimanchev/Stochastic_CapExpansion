using Distributed

@everywhere include("../src/Stochastic_CapExpansion.jl")

@everywhere using .Stochastic_CapExpansion
@everywhere using Revise, JuMP, Gurobi, HiGHS, Ipopt, DataFrames, CSV, YAML, Random, LinearAlgebra, Combinatorics, Dates, Distributions, Surrogates


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
        #cuts_to_keep = copy(cuts)
        rhs_values = output_cvar["RHS Values"]
        removed_cuts = output_cvar["Removed cuts"]
        if mapping
            #cuts_to_keep = filter(cut -> parse(Int, split(string(cut), "_")[end-1]) <= output_cvar["first_write"] + 10, cuts_to_keep)
        end
        vectors = vector_set !== nothing ? vector_set : generate_weights(iterations, length(MP[:x])+length(MP[:x_line]), settings["Vector Type"], settings)
        #@info("Keeping $(length(cuts_to_keep)) cuts for MGA iterations")
        for iteration in 1:iterations
            set_objective_bendersMP!(MP, "Capacity", inputs, settings; set_coeffs = vectors[iteration])
            #cuts_to_keep = manage_cuts(MP, cuts_to_keep)
            # other option - reactivate all cuts then let alg deactivate those that are not useful
            reactivate_cuts(MP, removed_cuts, rhs_values)

            output_random = mga_benders(inputs, settings, MP, SPs, budgets, "Random_"*string(iteration); Eval_SPs = Eval_SPs, mapping = mapping)
            log_result_memory!("Random_"*string(iteration)*" output", output_random)
            
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
            
            #if length(cuts_to_keep) < settings["Cuts retained"] && !mapping
             #   push!(cuts_to_keep, [name(con) for con in all_constraints(MP, include_variable_in_set_constraints=false) if startswith(string(con), "optimality_cut_") || startswith(string(con), "cvar_tail_cuts_")]...)
           # end
        end
        #write_gaps!(gaps, run_labels, joinpath(results_path, "Gaps"))

        # Map interior after exterior mapping
        if mapping
            all_caps = Matrix(vcat(results_cap...))
            @time samples = sample_interior(all_caps, n_samples, settings)
            outputs_mp = run_distributed_sampling(samples)
            for (i, sample) in enumerate(samples)
                outputs_sp = run_all_subproblems(SPs, inputs, settings, sample[1:R], sample[R+1:end]; minimal_payload=false)
                ev, cvar = evaluate_subproblems(outputs_sp, P_s, P_f, P_k, VaR_percent)
                @info("Sample $i: Investment cost = $(outputs_mp[i]["Inv_cost"]), Expected value = $ev, CVaR = $(cvar)")
                @info("Sample $i: Budget Percentage for CVaR: $(((cvar)/(budgets["CVaR"]*settings["Scaling factor cost"])*100)*100)%, Budget Percentage for Expected Value: $((((ev+outputs_mp[i]["Inv_cost"])/(budgets["System_Expected"]*settings["Scaling factor cost"]))*100)-100)%")
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

    # isolate capacity vectors
    index_indx = findfirst(==("Index"), names(interpolate_df))
    labels = Vector(interpolate_df[:, index_indx])
    select!(interpolate_df, Not(:Index))
    investment_indx = findfirst(==("Investment_Cost"), names(interpolate_df))
    cap_vectors = interpolate_df[:, 1:investment_indx-1]
    # get costs
    inv_costs = cap_vectors .* reshape(costs_by_resource, 1, :)
    tot_inv_costs = sum.(eachrow(Matrix(inv_costs)))
    
    col_by_zone = collect(findall(x -> occursin("z"*string(z), x), names(inv_costs)) for z in 1:Z)
    cost_by_zone = [sum.(eachrow(inv_costs[:, col_by_zone[z]])) for z in 1:Z]

    for i in 1:nrow(interpolate_df)
        @info("Evaluating interpolate ", labels[i])
        caps = Vector(cap_vectors[i, :])
        outputs_sp = run_all_subproblems(SPs, inputs, settings, caps, minimal_payload=false)
        ev, cvar = evaluate_subproblems(outputs_sp, P_s, P_f, P_k, VaR_percent)

        costs_by_zone_it = [cost_by_zone[z][i] for z in 1:Z]

        df_cap, df_syscost, df_emissions = make_results_mapping_dfs(caps, outputs_sp, tot_inv_costs[i], costs_by_zone_it, cvar, ev, inputs, settings)
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