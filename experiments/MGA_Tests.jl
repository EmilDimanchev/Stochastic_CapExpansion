using Distributed

@everywhere include("../src/Stochastic_CapExpansion.jl")

@everywhere using .Stochastic_CapExpansion
@everywhere using Revise, JuMP, Gurobi, HiGHS, DataFrames, CSV, YAML, Random, LinearAlgebra, Combinatorics, Dates, Distributions

function configure_parallel_workers!(settings::Dict)
    if settings["Parallel flag"]
        desired_workers = haskey(settings, "Workers") ? settings["Workers"] : max(Sys.CPU_THREADS - 1, 1)
        if nworkers() < desired_workers
            addprocs(desired_workers - nworkers())
        end

        project_root = abspath(joinpath(@__DIR__, ".."))
        init_expr = quote
            if !isdefined(Main, :Stochastic_CapExpansion)
                include(joinpath($project_root, "src", "Stochastic_CapExpansion.jl"))
            end
            using .Stochastic_CapExpansion
            using Revise, JuMP, Gurobi, HiGHS, DataFrames, CSV, YAML, Random, LinearAlgebra, Combinatorics, Dates, Distributions
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

function run_stochastic_exploration(SPs::Array{Model, 3}, inputs::Dict, settings::Dict, results_folder::String, summary_folder::String; budget_multiplier::Float64 = 1.10, vector_set::Union{AbstractVector, Nothing} = nothing, summary_name::String = "dual_constraint", Eval_SPs = nothing)

    configure_parallel_workers!(settings)

      # Result containers
    results_cap = []
    results_syscost = []
    results_emissions = []
    outputs = []
    run_labels = []
    gaps = []

    # Model Settings
    iterations = settings["Iterations"]
    risk_aversion_weight = settings["Risk aversion weight"] ### Currently set consistently across runs and ahead of time
    value_at_risk_percent = settings["Value-at-Risk percent"] ### Currently set consistently across runs and ahead of time
    # Create and set expected value model
    settings["Risk aversion flag"] = false

    MP = build_planning_model(inputs, settings)

    
    
    output_exp = benders_algorithm(inputs, settings, MP, SPs, "System_Expected"; Eval_SPs = Eval_SPs)
    gap_exp = output_exp["Gaps"]
    push!(gaps, gap_exp)
    push!(run_labels, "System_Expected")
    # Write results
    results_destination = joinpath(results_folder,"System_Expected")
    df_cap, df_syscost, df_emissions = write_results_benders(output_exp, inputs, settings, results_destination)
    push!(results_cap, df_cap)
    push!(results_syscost, df_syscost)
    push!(results_emissions, df_emissions)
    @info("System_Expected solution has investment cost of ", output_exp["MP"]["Inv_cost"])
    @info("Expected value system cost of " * "System_Expected" * " solution: $(output_exp["Expected Value"] + output_exp["MP"]["Inv_cost"])")
    @info("Risk adjusted system cost of " * "System_Expected" * " solution: $((1-risk_aversion_weight)*output_exp["CVaR"] + risk_aversion_weight*output_exp["Expected Value"]+ output_exp["MP"]["Inv_cost"])")

  
    # Base Runs
    settings["Risk aversion flag"] = true
    set_objective_bendersMP!(MP, "System_Weighted_CVaR", inputs, settings; obj_weight = risk_aversion_weight)
    output_cvar = benders_algorithm(inputs, settings, MP, SPs, "System_Weighted_CVaR"; Eval_SPs = Eval_SPs)
    gap_cvar = output_cvar["Gaps"]
    push!(run_labels, "System_Weighted_CVaR")
    push!(gaps, gap_cvar)
     # Write results
    results_destination = joinpath(results_folder,"System_Weighted_CVaR")
    df_cap, df_syscost, df_emissions = write_results_benders(output_cvar, inputs, settings, results_destination)
    push!(results_cap, df_cap)
    push!(results_syscost, df_syscost)
    push!(results_emissions, df_emissions)
    @info("System_Weighted_CVaR solution has investment cost of ", output_cvar["MP"]["Inv_cost"])
    @info("Expected value system cost of " * "System_Weighted_CVaR" * " solution: $(output_cvar["Expected Value"] + output_cvar["MP"]["Inv_cost"])")
    @info("Risk adjusted system cost of " * "System_Weighted_CVaR" * " solution: $((1-risk_aversion_weight)*output_cvar["CVaR"] + risk_aversion_weight*output_cvar["Expected Value"]+ output_cvar["MP"]["Inv_cost"])")



    if settings["Capacity Exploration"]
        budgets = Dict()
        # Set budgets? ------- budget set = set with budgets same percent greater than least cost solution for each metric
        #budgets = add_budget_constraint_bendersMP(MP, ((output_exp["MP"]["Inv_cost"] + output_exp["Expected Value"])/settings["Scaling factor cost"])*budget_multiplier, "System_Expected", budgets)
        #budgets = add_budget_constraint_bendersMP(MP, (((1-risk_aversion_weight)*output_cvar["CVaR"] + (risk_aversion_weight)*output_cvar["Expected Value"] + output_cvar["MP"]["Inv_cost"])/settings["Scaling factor cost"])*budget_multiplier, "System_Weighted_CVaR", budgets; risk_aversion=risk_aversion_weight)
        budgets = add_budget_constraint_bendersMP(MP, ((output_cvar["MP"]["Inv_cost"] + output_cvar["Expected Value"])/settings["Scaling factor cost"]), "System_Expected", budgets)
        budgets = add_budget_constraint_bendersMP(MP, (((1-risk_aversion_weight)*output_exp["CVaR"] + (risk_aversion_weight)*output_exp["Expected Value"] + output_exp["MP"]["Inv_cost"])/settings["Scaling factor cost"]), "System_Weighted_CVaR", budgets; risk_aversion=risk_aversion_weight)
        cuts_to_keep = [name(con) for con in all_constraints(MP, include_variable_in_set_constraints=false) if startswith(string(con), "optimality_cut_") || startswith(string(con), "cvar_tail_cuts_")]
        vectors = vector_set !== nothing ? vector_set : generate_weights(iterations, length(MP[:x]), settings["Vector Type"], settings)
        @info("Keeping $(length(cuts_to_keep)) cuts for MGA iterations")
        for iteration in 1:iterations
            set_objective_bendersMP!(MP, "Capacity", inputs, settings; set_coeffs = vectors[iteration])
            output_random = mga_benders(inputs, settings, MP, SPs, budgets, "Random_"*string(iteration); Eval_SPs = Eval_SPs)
            cuts_to_keep = manage_cuts(MP, cuts_to_keep)
            gap = output_random["Gaps"]
            push!(run_labels, "Random_"*string(iteration))
            push!(gaps, gap)
            results_destination = joinpath(results_folder,"Random_"*string(iteration))
            df_cap, df_syscost, df_emissions = write_results_benders(output_random, inputs, settings, results_destination; budgets = budgets)
            push!(results_cap, df_cap)
            push!(results_syscost, df_syscost)
            push!(results_emissions, df_emissions)
        end
        #write_gaps!(gaps, run_labels, joinpath(results_path, "Gaps"))
        write_exploration_results!(results_cap, results_syscost, results_emissions, joinpath(summary_folder), run_labels, summary_name)
    end
    return outputs, vectors
end


function run_stochastic_exploration_single_type(SPs::Array{Model, 3}, inputs::Dict, settings::Dict, results_path::String, summary_folder::String; type::String = "System_Expected", vector_set::Union{AbstractVector, Nothing} = nothing, budget_multiplier::Float64 = 1.10, Eval_SPs=nothing)

    configure_parallel_workers!(settings)

    # ~~~
    # Define paths
    # ~~~

    # Load everything
    # ~~~
    # Set up experiment
    # ~~~

    # Result containers
    results_cap = []
    results_syscost = []
    results_emissions = []
    outputs = []
    run_labels = []
    gaps = []

    

    # Model Settings
    iterations = settings["Iterations"]
    risk_aversion_weight = settings["Risk aversion weight"] ### Currently set consistently across runs and ahead of time
    value_at_risk_percent = settings["Value-at-Risk percent"] ### Currently set consistently across runs and ahead of time

    # Base Runs

    # Set type of model
    if type == "System_Expected"
        settings["Risk aversion flag"] = false
    elseif type == "System_Weighted_CVaR"
        settings["Risk aversion flag"] = true
    end
    MP = build_planning_model(inputs, settings)
    output = benders_algorithm(inputs, settings, MP, SPs, type;  Eval_SPs = Eval_SPs)
    gap = output["Gaps"]
    push!(outputs, output)
    push!(run_labels, type)
    push!(gaps, gap)
    # Write results
    results_destination = joinpath(results_path,type)
    df_cap, df_syscost, df_emissions = write_results_benders(output, inputs, settings, results_destination)
    push!(results_cap, df_cap)
    push!(results_syscost, df_syscost)
    push!(results_emissions, df_emissions)
    @info(type, " solution has investment cost of ", output["MP"]["Inv_cost"])
    @info("Expected value system cost of " * type * " solution: $(output["Expected Value"] + output["MP"]["Inv_cost"])")
    @info("Risk adjusted system cost of " * type * " solution: $(output["CVaR"] + output["MP"]["Inv_cost"])")

    if settings["Capacity Exploration"]
        budgets = Dict()
        # Set budgets?
        if type == "System_Expected"
            budgets = add_budget_constraint_bendersMP(MP, ((output["MP"]["Inv_cost"] + output["Expected Value"])/settings["Scaling factor cost"])*budget_multiplier, "System_Expected", budgets)
        elseif type == "System_Weighted_CVaR"
            budgets = add_budget_constraint_bendersMP(MP, (((1-risk_aversion_weight)*output["CVaR"] + (risk_aversion_weight)*output["Expected Value"] + output["MP"]["Inv_cost"]))/settings["Scaling factor cost"]*budget_multiplier, "System_Weighted_CVaR", budgets; risk_aversion=risk_aversion_weight)
        end
        cuts_to_keep = [name(con) for (F, S) in list_of_constraint_types(MP) for con in all_constraints(MP, F, S) if startswith(string(con), "optimality_cut_") || startswith(string(con), "cvar_tail_cuts_")]
        
        vectors = vector_set !== nothing ? vector_set : generate_weights(iterations, length(MP[:x]), settings["Vector Type"], settings)
        for iteration in 1:iterations
            set_objective_bendersMP!(MP, "Capacity", inputs, settings; set_coeffs = vectors[iteration])
            output_random = mga_benders(inputs, settings, MP, SPs, budgets, "Random_"*string(iteration); Eval_SPs = Eval_SPs)
            #cuts_to_keep = manage_cuts(MP, cuts_to_keep)
            gap = output_random["Gaps"]
            push!(outputs, output_random)
            push!(run_labels, "Random_"*string(iteration))
            push!(gaps, gap)
            results_destination = joinpath(results_path,"Random_"*string(iteration))
            df_cap, df_syscost, df_emissions = write_results_benders(output_random, inputs, settings, results_destination; budgets = budgets)
            push!(results_cap, df_cap)
            push!(results_syscost, df_syscost)
            push!(results_emissions, df_emissions)
        end
        #write_gaps!(gaps, run_labels, joinpath(results_path, "Gaps"))
        write_exploration_results!(results_cap, results_syscost, results_emissions, joinpath(summary_folder), run_labels, type)
    end
    return outputs, vectors
    
end



#### Base MGA test - run base model with one scenario selected
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
    gap = output["Gaps"]
    push!(outputs, output)
    push!(labels, "OneScenarioLC")
    push!(gaps, gap)
    results_destination = joinpath(results_path,"CostOptimal")
    df_cap, df_syscost, df_emissions = write_results_benders(output, new_inputs, settings, results_destination)
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
        #cuts_to_keep = manage_cuts(MP, cuts_to_keep)
        gap = output_random["Gaps"]
        results_destination = joinpath(results_path,"Base_MGA_"*string(iteration))
        df_cap, df_syscost, df_emissions = write_results_benders(output_random, new_inputs, settings, results_destination; budgets = budgets)
        push!(results_cap, df_cap)
        push!(results_syscost, df_syscost)
        push!(results_emissions, df_emissions)
        push!(outputs, output_random)
        push!(labels, "Random_"*string(iteration))
        push!(gaps, gap)
    end
    #write_gaps!(gaps, labels, joinpath(results_path, "Gaps"))
    write_exploration_results!(results_cap, results_syscost, results_emissions, summary_folder, labels, string(scenario))
end

function simple_comp()
    inputs_folder = joinpath("inputs","Inputs_30repdays_ext_1000scen_7techs")
    test_index = 10
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
    outputs_cvar, vectors = run_stochastic_exploration_single_type(SPs, inputs, settings, joinpath(results_folder, "CVaR"), summary_folder; type ="System_Weighted_CVaR",vector_set = nothing, budget_multiplier=1.10)
    
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
    #_ = run_base_mga(SPs, new_inputs, settings, results_path, summary_folder; budget_multiplier=1.10, vector_set=nothing, scenario=i)

end


function test_stability_laptop(test_index)
    inputs_folder = joinpath("inputs", "Inputs_30repdays_ext_1000scen_7techs")
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
    probabilities = [inputs["Demand scenario probabilities"], inputs["Gas price scenario probabilities"], inputs["Weather scenario probabilities"]]
    
    distribution_types = [["Gaussian", "Gaussian", "Gaussian"], ["LogNormal", "LogNormal", "LogNormal"], ["Gaussian", "LogNormal", "Gaussian"], ["LogNormal", "Gaussian", "LogNormal"]]
    names = ["Gaussian", "LogNormal", "Gaussian-LogNormal", "LogNormal-Gaussian"]

    # Set evaluation probabilities
    inputs["Full Demand scenario probabilities"] = inputs["Demand scenario probabilities"]
    inputs["Full Gas price scenario probabilities"] = inputs["Gas price scenario probabilities"]
    inputs["Full Weather scenario probabilities"] = inputs["Weather scenario probabilities"]

    # Base Run

    configure_parallel_workers!(settings)

    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)

    _, vectors = run_stochastic_exploration(SPs, inputs, settings, joinpath(results_folder, "Flat"), summary_folder; budget_multiplier=1.10, vector_set=nothing, summary_name = "Flat", Eval_SPs = SPs)

    # Set new probabilities in inputs and run again
    for (i, probs) in enumerate(distribution_types)
        println("************* Beginning run with new probabilities - ", names[i], " **************")
        new_probabilities = generate_probabilities(probabilities, probs)
        inputs["Demand scenario probabilities"] = new_probabilities[1]
        inputs["Gas price scenario probabilities"] = new_probabilities[2]
        #inputs["Weather scenario probabilities"] = new_probabilities[3]

        _, vectors = run_stochastic_exploration(SPs, inputs, settings, joinpath(results_folder, names[i]), summary_folder; budget_multiplier=1.10, vector_set=vectors, summary_name = names[i], Eval_SPs = SPs)
    end

end

function test_stability_della(test_index)
    inputs_folder = "/home/ml6802/Stochastic_CapExpansion/inputs/Inputs_30repdays_ext_1000scen_7techs_Della"#joinpath("inputs", "Inputs_30repdays_ext_1000scen_7techs")
    results_folder = "/home/ml6802/Stochastic_CapExpansion/outputs/Test_"*string(test_index) #joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    if !isdir(results_folder)
        mkpath(results_folder)
    end
    if !isdir(summary_folder)
        mkpath(summary_folder)
    end
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)
    probabilities = [inputs["Demand scenario probabilities"], inputs["Gas price scenario probabilities"], inputs["Weather scenario probabilities"]]
    
    distribution_types = [["Gaussian", "Gaussian", "Gaussian"], ["LogNormal", "LogNormal", "LogNormal"], ["Gaussian", "LogNormal", "Gaussian"], ["LogNormal", "Gaussian", "LogNormal"]]
    names = ["Gaussian", "LogNormal", "Gaussian-LogNormal", "LogNormal-Gaussian"]

    # Set evaluation probabilities
    inputs["Full Demand scenario probabilities"] = inputs["Demand scenario probabilities"]
    inputs["Full Gas price scenario probabilities"] = inputs["Gas price scenario probabilities"]
    inputs["Full Weather scenario probabilities"] = inputs["Weather scenario probabilities"]

    # Base Run

    configure_parallel_workers!(settings)

    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)

    _, vectors = run_stochastic_exploration(SPs, inputs, settings, joinpath(results_folder, "Flat"), summary_folder; budget_multiplier=1.10, vector_set=nothing, summary_name = "Flat")

    # Set new probabilities in inputs and run again
    for (i, probs) in enumerate(distribution_types)
        println("************* Beginning run with new probabilities - ", names[i], " **************")
        new_probabilities = generate_probabilities(probabilities, probs)
        inputs["Demand scenario probabilities"] = new_probabilities[1]
        inputs["Gas price scenario probabilities"] = new_probabilities[2]
        #inputs["Weather scenario probabilities"] = new_probabilities[3]

        _, vectors = run_stochastic_exploration(SPs, inputs, settings, joinpath(results_folder, names[i]), summary_folder; budget_multiplier=1.10, vector_set=vectors, summary_name = names[i], Eval_SPs = SPs)
    end

end

function generate_probabilities(Probablities::AbstractVector, Distribution_Type::AbstractVector)
    result = Vector{Any}(undef, length(Distribution_Type))
    parameters_normal = Dict("mean" => 5.5, "std" => 2.0)
    parameters_LogNormal = Dict("mean" => 2, "std" => 1)
    for (i, dist) in enumerate(Distribution_Type)
        if dist == "Gaussian"
            result[i] = generate_normalized_distribution(length(Probablities[i]), dist, parameters_normal)
        elseif dist == "LogNormal"
            result[i] = generate_normalized_distribution(length(Probablities[i]), dist, parameters_LogNormal)
        else
            error("Distribution type not supported")
        end
    end

    return result
end

function generate_normalized_distribution(num_scenarios::Int, distribution_type::String, parameters::Dict)
    if distribution_type == "Gaussian"
        dist = Normal(parameters["mean"], parameters["std"]) # Example parameters, can be adjusted
        samples = pdf.(dist, collect(1:num_scenarios))
        normalized_samples = samples ./ sum(samples)
        return normalized_samples
    elseif distribution_type == "LogNormal"
        dist = LogNormal(parameters["mean"], parameters["std"]) # Example parameters, can be adjusted
        samples = pdf.(dist, collect(1:num_scenarios))
        normalized_samples = samples ./ sum(samples)
        return normalized_samples
    else
        error("Distribution type not supported")
    end
end