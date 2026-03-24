
include("../src/Stochastic_CapExpansion.jl")

using .Stochastic_CapExpansion
using Revise, JuMP, Gurobi, DataFrames, CSV, YAML, Random, LinearAlgebra, Combinatorics, Dates

function run_stochastic_exploration(SPs::Array{Model, 3}, inputs::Dict, settings::Dict, results_path::String, summary_folder::String; budget_multiplier::Float64 = 1.10, vector_set::Union{AbstractVector, Nothing} = nothing)

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

    # Create and set expected value model
    settings["Risk aversion flag"] = false
    MP = build_planning_model(inputs, settings)
    output_exp = benders_algorithm(inputs, settings, MP, SPs)
    gap_exp = output_exp["Gaps"]
    push!(outputs, output_exp)
    push!(gaps, gap_exp)
    push!(run_labels, type)
    # Write results
    results_destination = joinpath(results_path,type)
    df_cap, df_syscost, df_emissions = write_results(output_exp, inputs, settings, results_destination)
    push!(results_cap, df_cap)
    push!(results_syscost, df_syscost)
    push!(results_emissions, df_emissions)
    @info(type, " solution has investment cost of ", output_exp["MP"]["Inv_cost"])
    @info("Expected value system cost of " * type * " solution: $(output_exp["Expected Value"] + output_exp["MP"]["Inv_cost"])")
    @info("Risk adjusted system cost of " * type * " solution: $((1-risk_aversion_weight)*output_exp["CVaR"] + risk_aversion_weight*output_exp["Expected Value"]+ output_exp["MP"]["Inv_cost"])")

    settings["Risk aversion flag"] = true
    set_objective_bendersMP!(MP, "System_Weighted_CVaR", inputs, settings; obj_weight = risk_aversion_weight)
    output_cvar = benders_algorithm(inputs, settings, MP, SPs)
    gap_cvar = output_cvar["Gaps"]
    push!(outputs, output_cvar)
    push!(run_labels, "System_Weighted_CVaR")
    push!(gaps, gap_cvar)
     # Write results
    results_destination = joinpath(results_path,type)
    df_cap, df_syscost, df_emissions = write_results_benders(output_cvar, inputs, settings, results_destination)
    push!(results_cap, df_cap)
    push!(results_syscost, df_syscost)
    push!(results_emissions, df_emissions)
    @info(type, " solution has investment cost of ", output_cvar["MP"]["Inv_cost"])
    @info("Expected value system cost of " * type * " solution: $(output_cvar["Expected Value"] + output_cvar["MP"]["Inv_cost"])")
    @info("Risk adjusted system cost of " * type * " solution: $((1-risk_aversion_weight)*output_cvar["CVaR"] + risk_aversion_weight*output_cvar["Expected Value"]+ output_cvar["MP"]["Inv_cost"])")

    if settings["Capacity Exploration"]
        budgets = Dict()
        # Set budgets? ------- budget set = set with budgets same percent greater than least cost solution for each metric
        budgets = add_budget_constraint_bendersMP(MP, (output_exp["MP"]["Inv_cost"] + output_exp["Expected Value"])*budget_multiplier, "System_Expected", budgets)
        budgets = add_budget_constraint_bendersMP(MP, ((1-risk_aversion_weight)*output_cvar["CVaR"] + (risk_aversion_weight)*output_cvar["Expected Value"] + output_cvar["MP"]["Inv_cost"])*budget_multiplier, "System_Weighted_CVaR", budgets; risk_aversion=risk_aversion_weight)

        vectors = vector_set !== nothing ? vector_set : generate_weights(iterations, length(MP[:x]), settings["Vector Type"], settings)
        for iteration in 1:iterations
            set_objective_bendersMP!(MP, "Capacity", inputs, settings; set_coeffs = vectors[iteration])
            output_random = mga_benders(inputs, settings, MP, SPs, budgets)
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
        write_gaps!(gaps, run_labels, joinpath(results_path, "Gaps"))
        write_exploration_results!(results_cap, results_syscost, results_emissions, joinpath(summary_folder), run_labels, type)
    end
    return outputs, vectors
end


function run_stochastic_exploration_single_type(SPs::Array{Model, 3}, inputs::Dict, settings::Dict, results_path::String, summary_folder::String; type::String = "System_Expected", vector_set::Union{AbstractVector, Nothing} = nothing, budget_multiplier::Float64 = 1.10)

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
    output = benders_algorithm(inputs, settings, MP, SPs)
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
            budgets = add_budget_constraint_bendersMP(MP, (output["MP"]["Inv_cost"] + output["Expected Value"])*budget_multiplier, "System_Expected", budgets)
        elseif type == "System_Weighted_CVaR"
            budgets = add_budget_constraint_bendersMP(MP, ((1-risk_aversion_weight)*output["CVaR"] + (risk_aversion_weight)*output["Expected Value"] + output["MP"]["Inv_cost"])*budget_multiplier, "System_Weighted_CVaR", budgets; risk_aversion=risk_aversion_weight)
        end

        vectors = vector_set !== nothing ? vector_set : generate_weights(iterations, length(MP[:x]), settings["Vector Type"], settings)
        for iteration in 1:iterations
            set_objective_bendersMP!(MP, "Capacity", inputs, settings; set_coeffs = vectors[iteration])
            output_random = mga_benders(inputs, settings, MP, SPs, budgets)
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
        write_gaps!(gaps, run_labels, joinpath(results_path, "Gaps"))
        write_exploration_results!(results_cap, results_syscost, results_emissions, joinpath(summary_folder), run_labels, type)
    end
    return outputs, vectors
    
end



#### Base MGA test - run base model with one scenario selected
function run_base_mga(base_SPs, inputs::Dict, settings::Dict, results_path::String, summary_folder::String; budget_multiplier::Float64 = 1.1, vector_set::Union{AbstractVector, Nothing} = nothing, scenario::Int = -1)
    # Result containers
    outputs = []
    labels = [] 
    results_cap = []
    results_syscost = []
    results_emissions = []
    budget = 0.0

    # Model Settings
    iterations = settings["Iterations"]
    

    

    MP = build_planning_model(inputs, settings) #### When loaded with one scenario weighted, this is equivalent to a deterministic model with that scenario selected
    SP_one_scen = build_all_subproblems(inputs, settings)
    output = benders_algorithm(inputs, settings, MP, SP_one_scen; eval_SPs = base_SPs)
    gap = output["Gaps"]
    push!(outputs, output)
    push!(labels, "OneScenarioLC")
    push!(gaps, gap)
    results_destination = joinpath(results_path,"CostOptimal")
    df_cap, df_syscost, df_emissions = write_results_benders(output, inputs, settings, results_destination)
    push!(results_cap, df_cap)
    push!(results_syscost, df_syscost)
    push!(results_emissions, df_emissions)
    lc_value = output["Investment cost"] + output["Expected cost"]
    @info("System cost of cost optimal solution: ", lc_value)
    if budget_multiplier <= 10
        budget = lc_value * (budget_multiplier)
    end
    
    vectors = vector_set !== nothing ? vector_set : generate_weights(iterations, length(MP[:x]), settings["Vector Type"], settings)
    @info("Using budget of ", budget, " for Base MGA test")
    percent_over_lc = round((budget - lc_value)/lc_value * 100, digits=2)
    @info("This budget is ", percent_over_lc, "% over the cost optimal solution")
    budgets = Dict()
    budgets = add_budget_constraint_bendersMP(MP, budget, "System_Expected", budgets)

    for iteration in 1:iterations
        set_objective_bendersMP!(MP, "Capacity", inputs, settings; set_coeffs = vectors[iteration])
        output_random = mga_benders(inputs, settings, MP, SPs, budgets; Eval_SPs = base_SPs)
        gap = output_random["Gaps"]
        results_destination = joinpath(results_path,"Base_MGA_"*string(iteration))
        df_cap, df_syscost, df_emissions = write_results_benders(output_random, inputs, settings, results_destination; budgets = budgets)
        push!(results_cap, df_cap)
        push!(results_syscost, df_syscost)
        push!(results_emissions, df_emissions)
        push!(outputs, output_random)
        push!(labels, "Random_"*string(iteration))
        push!(gaps, gap)
    end
    write_gaps!(gaps, labels, joinpath(results_path, "Gaps"))
    write_exploration_results!(results_cap, results_syscost, results_syscost_risk, results_emissions, investment_costs, results_oper_all, summary_folder, labels, string(scenario))
end

function simple_comp()
    inputs_folder = joinpath("inputs","Inputs_30repdays_ext_1000scen_7techs")
    test_index = 4
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

    if settings["Parallel flag"]
        settings["Threads"] = Threads.nthreads()
        @info("Running with parallelization using $(Threads.nthreads()) threads")
    else
        @info("Running without parallelization")
    end

    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)


    outputs_exp, vectors = run_stochastic_exploration_single_type(SPs, inputs, settings, results_folder, summary_folder; type = "System_Expected", vector_set = nothing, budget_multiplier=1.10)
    outputs_cvar, vectors = run_stochastic_exploration_single_type(SPs, inputs, settings, results_folder, summary_folder; type ="System_Weighted_CVaR",vector_set = vectors, budget_multiplier=1.10)
    outputs_mixed, vectors = run_stochastic_exploration(SPs, inputs, settings, results_folder, summary_folder; budget_multiplier=1.10, vector_set=vectors)
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