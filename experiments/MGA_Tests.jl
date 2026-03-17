
include("../src/Stochastic_CapExpansion.jl")

using .Stochastic_CapExpansion

using JuMP, Gurobi, DataFrames, CSV, Random, LinearAlgebra
using Revise, YAML

function run_stochastic_exploration(inputs::Dict, settings::Dict, results_path::String, summary_folder::String; budget_multiplier::Float64 = 1.10)

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
    results_syscost_risk = []
    results_emissions = []
    results_oper_all = []
    outputs = []
    investment_costs = []
    run_labels = []
    scenario_names = ["S1", "S2", "S3", "S4"]
    

    # Model Settings
    iterations = settings["Iterations"]
    risk_aversion_weight = settings["Risk aversion weight"] ### Currently set consistently across runs and ahead of time
    value_at_risk_percent = settings["Value-at-Risk percent"] ### Currently set consistently across runs and ahead of time

    # Base Runs

    # Build base stochastic model and run for lowest first stage cost solution
    model = build_optimization_model(inputs, settings, "System_Expected")
    output_exp, model = run_optimization_model(model, inputs, settings)
    push!(outputs, output_exp)
    push!(run_labels, "System_Expected")
    # Write results
    results_destination = joinpath(results_path,"EV")
    df_cap, df_syscost, df_syscost_risk, df_emissions, df_oper_all = write_results(output_exp, inputs, settings, results_destination)
    push!(results_cap, df_cap)
    push!(results_syscost, df_syscost)
    push!(results_syscost_risk, df_syscost_risk)
    push!(results_emissions, df_emissions)
    push!(results_oper_all, df_oper_all)
    push!(investment_costs, output_exp["Investment cost"])
    @info("Expected value system cost of EV solution: ", output_exp["Expected cost"])

    #model = build_optimization_model(inputs, settings, "CVaR", 0.75)
    set_objective!(model, "System_Weighted_CVaR"; obj_weight=risk_aversion_weight)

    output_cvar, model = run_optimization_model(model, inputs, settings)
    results_destination = joinpath(results_path,"CVaR")
    push!(outputs, output_cvar)
    push!(run_labels, "CVaR_"*string(risk_aversion_weight)*"_"*string(value_at_risk_percent))
    df_cap, df_syscost, df_syscost_risk, df_emissions, df_oper_all = write_results(output_cvar, inputs, settings, results_destination)
    push!(results_cap, df_cap)
    push!(results_syscost, df_syscost)
    push!(results_syscost_risk, df_syscost_risk)
    push!(results_emissions, df_emissions)
    push!(results_oper_all, df_oper_all)
    push!(investment_costs, output_cvar["Investment cost"])
    @info("Expected value system cost of Weighted CVaR solution: ", output_cvar["Expected cost"])


    
    @info("Risk adjusted system cost of EV solution: ", output_exp["Risk adjusted system cost"])
    @info("Risk adjusted system cost of CVaR solution: ", output_cvar["Risk adjusted system cost"])

    if settings["Capacity Exploration"]
        # Set budgets?
        add_budget_constraint!(model, (output_exp["Investment cost"] + output_exp["Expected cost"])*budget_multiplier, "System_Expected")
        add_budget_constraint!(model, (output_cvar["Risk adjusted system cost"])*budget_multiplier, "System_Weighted_CVaR"; risk_aversion=risk_aversion_weight)

        vectors = generate_weights(iterations, length(model[:x]), settings["Vector Type"])
        for iteration in 1:iterations
            set_objective!(model, "Capacity"; set_coeffs = vectors[iteration])
            output_random, model = run_optimization_model(model, inputs, settings)
            push!(outputs, output_random)
            push!(run_labels, "Random_"*string(iteration))
            results_destination = joinpath(results_path,"Random_"*string(iteration))
            df_cap, df_syscost, df_syscost_risk, df_emissions, df_oper_all = write_results(output_random, inputs, settings, results_destination)
            push!(results_cap, df_cap)
            push!(results_syscost, df_syscost)
            push!(results_syscost_risk, df_syscost_risk)
            push!(results_emissions, df_emissions)
            push!(results_oper_all, df_oper_all)
            push!(investment_costs, output_random["Investment cost"])
        end

        write_exploration_results!(results_cap, results_syscost, results_syscost_risk, results_emissions, investment_costs, results_oper_all, joinpath(summary_folder), run_labels, "stoch")
    end
    return outputs, vectors
    
end


function run_stochastic_exploration_single_type(inputs::Dict, settings::Dict, results_path::String, summary_folder::String; type::String = "System_Expected", vector_set::Union{AbstractVector, Nothing} = nothing, budget_multiplier::Float64 = 1.10)

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
    results_syscost_risk = []
    results_emissions = []
    results_oper_all = []
    outputs = []
    investment_costs = []
    run_labels = []
    scenario_names = ["S1", "S2", "S3", "S4"]
    

    # Model Settings
    iterations = settings["Iterations"]
    risk_aversion_weight = settings["Risk aversion weight"] ### Currently set consistently across runs and ahead of time
    value_at_risk_percent = settings["Value-at-Risk percent"] ### Currently set consistently across runs and ahead of time

    # Base Runs

    # Build base stochastic model and run for lowest first stage cost solution
    model = build_optimization_model(inputs, settings, type)
    @info(objective_function(model))
    output, model = run_optimization_model(model, inputs, settings)
    push!(outputs, output)
    push!(run_labels, type)
    # Write results
    results_destination = joinpath(results_path,type)
    df_cap, df_syscost, df_syscost_risk, df_emissions, df_oper_all = write_results(output, inputs, settings, results_destination)
    push!(results_cap, df_cap)
    push!(results_syscost, df_syscost)
    push!(results_syscost_risk, df_syscost_risk)
    push!(results_emissions, df_emissions)
    push!(results_oper_all, df_oper_all)
    push!(investment_costs, output["Investment cost"])
    @info(type, " solution has investment cost of ", output["Investment cost"])
    @info("Expected value system cost of " * type * " solution: ", output["Investment cost"] + output["Expected cost"])
    @info("Risk adjusted system cost of " * type * " solution: ", output["Risk adjusted system cost"])

    if settings["Capacity Exploration"]
        # Set budgets?
        if type == "System_Expected"
            add_budget_constraint!(model, (output["Investment cost"] + output["Expected cost"])*budget_multiplier, "System_Expected")
        elseif type == "System_Weighted_CVaR"
            add_budget_constraint!(model, (output["Risk adjusted system cost"])*budget_multiplier, "System_Weighted_CVaR"; risk_aversion=risk_aversion_weight)
        end

        vectors = vector_set !== nothing ? vector_set : generate_weights(iterations, length(model[:x]), settings["Vector Type"])
        for iteration in 1:iterations
            set_objective!(model, "Capacity"; set_coeffs = vectors[iteration])
            output_random, model = run_optimization_model(model, inputs, settings)
            push!(outputs, output_random)
            push!(run_labels, "Random_"*string(iteration))
            results_destination = joinpath(results_path,"Random_"*string(iteration))
            df_cap, df_syscost, df_syscost_risk, df_emissions, df_oper_all = write_results(output_random, inputs, settings, results_destination)
            push!(results_cap, df_cap)
            push!(results_syscost, df_syscost)
            push!(results_syscost_risk, df_syscost_risk)
            push!(results_emissions, df_emissions)
            push!(results_oper_all, df_oper_all)
            push!(investment_costs, output_random["Investment cost"])
        end

        write_exploration_results!(results_cap, results_syscost, results_syscost_risk, results_emissions, investment_costs, results_oper_all, joinpath(summary_folder), run_labels, type)
    end
    return outputs, vectors
    
end



#### Base MGA test - run base model with one scenario selected
function run_base_mga(inputs::Dict, settings::Dict, results_path::String, summary_folder::String; budget_multiplier::Float64 = 1.1, vector_set::Union{AbstractVector, Nothing} = nothing, scenario::Int = -1)
    # Result containers
    outputs = []
    labels = [] 
    results_cap = []
    results_syscost = []
    results_syscost_risk = []
    results_emissions = []
    results_oper_all = []
    investment_costs = []
    budget = 0.0

    # Model Settings
    iterations = settings["Iterations"]

    model = build_optimization_model(inputs, settings, "System_Expected") #### When loaded with one scenario weighted, this is equivalent to a deterministic model with that scenario selected
    output, model = run_optimization_model(model, inputs, settings)
    push!(outputs, output)
    push!(labels, "OneScenarioLC")
    results_destination = joinpath(results_path,"CostOptimal")
    df_cap, df_syscost, df_syscost_risk, df_emissions, df_oper_all = write_results(output, inputs, settings, results_destination)
    push!(results_cap, df_cap)
    push!(results_syscost, df_syscost)
    push!(results_syscost_risk, df_syscost_risk)
    push!(results_emissions, df_emissions)
    push!(results_oper_all, df_oper_all)
    push!(investment_costs, output["Investment cost"])
    lc_value = output["Investment cost"] + output["Expected cost"]
    @info("System cost of cost optimal solution: ", lc_value)
    if budget_multiplier <= 10
        budget = lc_value * (budget_multiplier)
    end
    
    vectors = vector_set !== nothing ? vector_set : generate_weights(iterations, length(model[:x]), settings["Vector Type"])
    @info("Using budget of ", budget, " for Base MGA test")
    percent_over_lc = round((budget - lc_value)/lc_value * 100, digits=2)
    @info("This budget is ", percent_over_lc, "% over the cost optimal solution")
    add_budget_constraint!(model, budget, "System_Expected")

    for iteration in 1:iterations
        set_objective!(model, "Capacity"; set_coeffs = vectors[iteration])
        output_random, model = run_optimization_model(model, inputs, settings)
        results_destination = joinpath(results_path,"Base_MGA_"*string(iteration))
        df_cap, df_syscost, df_syscost_risk, df_emissions, df_oper_all = write_results(output_random, inputs, settings, results_destination)
        push!(results_cap, df_cap)
        push!(results_syscost, df_syscost)
        push!(results_syscost_risk, df_syscost_risk)
        push!(results_emissions, df_emissions)
        push!(outputs, output_random)
        push!(labels, "Random_"*string(iteration))
        push!(results_oper_all, df_oper_all)
        push!(investment_costs, output_random["Investment cost"])
    end

    write_exploration_results!(results_cap, results_syscost, results_syscost_risk, results_emissions, investment_costs, results_oper_all, summary_folder, labels, string(scenario))
end

function complex_comp()
    budget_exp = []
    budget_cvar = []
    
    vectors = []
    EXPERIMENTS = ["Stochastic Exploration", "Base MGA"] ### Options: "Stochastic Exploration", "Base MGA"
    inputs_folder = joinpath("inputs","Inputs_30days_5techs")
    test_index = 1
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
    if settings["Seed"] !== nothing
        Random.seed!(settings["Seed"])
    end
    a = [1.0 0.0]
    b = [0.0 1.0]
    probabilities_scenarios = [[a; a], [a; b], [b; a], [b; b]]
    println("Demand and fuel price scenario probabilities: ", probabilities_scenarios)
    
    if "Stochastic Exploration" in EXPERIMENTS
        results_path = joinpath(results_folder,"Results_Stochastic_Exploration/")
        outputs, vectors = run_stochastic_exploration(inputs, settings, results_path, summary_folder)
        budget_exp = collect(outputs[1]["Operating cost"][i,j] for i in 1:size(outputs[1]["Operating cost"], 1) for j in 1:size(outputs[1]["Operating cost"], 2)) .+ outputs[1]["Investment cost"]
        budget_cvar = collect(outputs[2]["Operating cost"][i,j] for i in 1:size(outputs[2]["Operating cost"], 1) for j in 1:size(outputs[2]["Operating cost"], 2)) .+ outputs[2]["Investment cost"]
    end
    @info("Expected value solution budgets: ", budget_exp)
    @info("CVaR solution budgets: ", budget_cvar)
    if "Base MGA" in EXPERIMENTS
        for (i, scenario) in enumerate(probabilities_scenarios)
            budget = budget_exp[i]
            println("Using budget of ", budget, " for Base MGA test, which is the expected cost of the Expected Value solution from the stochastic exploration experiment")
            println("Testing with demand scenario probabilities: ", scenario[1,:], " and fuel price scenario probabilities: ", scenario[2,:])
            inputs["Demand scenario probabilities"] = scenario[1,:]
            inputs["Fuel price scenario probabilities"] = scenario[2, :]
            results_path = joinpath(results_folder, "Results_Base_MGA", "Scenario_"*string(i))
            outputs = run_base_mga(inputs, settings, results_path, summary_folder; budget_multiplier=budget_exp[i], vector_set=vectors, scenario=i)
        end
    end
end

function simple_comp()
    inputs_folder = joinpath("inputs","Inputs_30days_5techs")
    test_index = 2
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
    if settings["Seed"] !== nothing
        Random.seed!(settings["Seed"])
    end
    settings["Demand risk flag"] = false
    settings["Fuel risk flag"] = false
    """
    a = [1.0 0.0]
    b = [0.0 1.0]
    probabilities_scenarios = [[a; a], [a; b], [b; a], [b; b]]
    println("Demand and fuel price scenario probabilities: ", probabilities_scenarios)
    #outputs_exp, vectors = run_stochastic_exploration_single_type(inputs, settings, results_folder, summary_folder; type = "System_Expected", vector_set = nothing, budget_multiplier=1.10)
    #outputs_cvar, vectors = run_stochastic_exploration_single_type(inputs, settings, results_folder, summary_folder; type ="System_Weighted_CVaR",vector_set = nothing, budget_multiplier=1.10)
    outputs_mixed, vectors = run_stochastic_exploration(inputs, settings, results_folder, summary_folder; budget_multiplier=1.10)
"""
    #for (i, scenario) in enumerate(probabilities_scenarios)
     #   println("Testing with demand scenario probabilities: ", scenario[1,:], " and fuel price scenario probabilities: ", scenario[2,:])
     #   inputs["Demand scenario probabilities"] = scenario[1,:]
    #    inputs["Fuel price scenario probabilities"] = scenario[2, :]
        i = 0
        results_path = joinpath(results_folder, "Results_Base_MGA", "Scenario_"*string(i))
        _ = run_base_mga(inputs, settings, results_path, summary_folder; budget_multiplier=1.10, vector_set=nothing, scenario=i)
    #end

end

simple_comp()