
include("../src/Stochastic_CapExpansion.jl")

using .Stochastic_CapExpansion

using JuMP, Gurobi, DataFrames, CSV, Random, LinearAlgebra
using Revise, YAML

function run_stochastic_exploration(inputs_path::String, results_path::String)

    # ~~~
    # Define paths
    # ~~~

    # Load everything
    # ~~~
    # Set up experiment
    # ~~~

    # Model settings
    settings = load_settings(inputs_path)
    
    # Load data
    inputs = load_input_data(inputs_path, settings)

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
    model = build_optimization_model(inputs, settings, "Expectation")
    output_exp, model = run_optimization_model(model, inputs, settings)
    push!(outputs, output_exp)
    push!(run_labels, "Expectation")
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
    set_objective!(model, "Weighted CVaR"; obj_weight=risk_aversion_weight)

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
        add_budget_constraint!(model, (output_cvar["Investment cost"] + output_cvar["Expected cost"]), "System_Expected")
        add_budget_constraint!(model, (output_cvar["Risk adjusted system cost"]), "System_CVaR"; risk_aversion=risk_aversion_weight)

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

        write_exploration_results!(results_cap, results_syscost, results_syscost_risk, results_emissions, investment_costs, results_oper_all, results_path, run_labels)
    end
    return outputs
    
end

#### Base MGA test - run base model with one scenario selected
function run_base_mga(inputs_path::String, results_path::String; budget::Float64 = .1)
    settings = load_settings(inputs_path)
    inputs = load_input_data(inputs_path, settings)

    # Result containers
    outputs = []
    labels = [] 
    results_cap = []
    results_syscost = []
    results_syscost_risk = []
    results_emissions = []
    results_oper_all = []
    investment_costs = []
    scenario_names = ["S1", "S2", "S3", "S4"]

    # Model Settings
    iterations = settings["Iterations"]

    model = build_optimization_model(inputs, settings, "Expectation") #### When loaded with one scenario weighted, this is equivalent to a deterministic model with that scenario selected
    output, model = run_optimization_model(model, inputs, settings)
    push!(outputs, output)
    push!(labels, "OneScenarioLC")
    results_destination = joinpath(results_path,"CostOptimal_S1")
    df_cap, df_syscost, df_syscost_risk, df_emissions, df_oper_all = write_results(output, inputs, settings, results_destination)
    push!(results_cap, df_cap)
    push!(results_syscost, df_syscost)
    push!(results_syscost_risk, df_syscost_risk)
    push!(results_emissions, df_emissions)
    push!(results_oper_all, df_oper_all)
    push!(investment_costs, output["Investment cost"])
    if budget <= 2
        budget = (output["Investment cost"] + output["Expected cost"])*(1 + budget)
    end

    vectors = generate_weights(iterations, length(model[:x]), settings["Vector Type"])
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

    write_exploration_results!(results_cap, results_syscost, results_syscost_risk, results_emissions, investment_costs, results_oper_all, results_path, labels)

    return outputs
end

function main()
    budget_exp_s1 = 0.0
    budget_cvar_s1 = 0.0
    Random.seed!(1234)
    EXPERIMENTS = ["Stochastic Exploration","Base MGA"]
    if "Stochastic Exploration" in EXPERIMENTS
        inputs_folder = "Inputs_30days_5techs"
        inputs_path = joinpath("inputs",inputs_folder)
        results_path = joinpath("outputs","Results_Stochastic_Exploration/")
        outputs = run_stochastic_exploration(inputs_path, results_path)
        budget_exp_s1 = outputs[1]["Investment cost"] + outputs[1]["Operating cost"][1,1]
        budget_cvar_s1 = outputs[2]["Investment cost"] + outputs[2]["Operating cost"][1,1]
    end
    if "Base MGA" in EXPERIMENTS
        budget = budget_exp_s1
        @info("Using budget of ", budget, " for Base MGA test, which is the expected cost of the Expected Value solution from the stochastic exploration experiment")
        inputs_folder = "Inputs_30days_5techs_OneScenario"
        inputs_path = joinpath("inputs",inputs_folder)
        results_path = joinpath("outputs","Results_Base_MGA/")
        outputs = run_base_mga(inputs_path, results_path)
    end
end

main()