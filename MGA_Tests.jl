using Revise
includet("Stochastic_CapExpansion.jl")

using Pkg
using .Stochastic_CapExpansion
using CSV, DataFrames
using JuMP





function run_stochastic_exploration()

    # ~~~
    # Define paths
    # ~~~

    inputs_folder = "Inputs_30days_5techs"

    results_folder = "Results_Stochastic_Exploration"

    inputs_path = string("./Inputs/",inputs_folder)

    results_path = "./"

    # Load everything
    # ~~~
    # Set up experiment
    # ~~~

    # Model settings
    settings = Dict{String, Any}()
    settings["Results folder path"] = results_path
    # Risk
    settings["Demand risk flag"] = true
    settings["Fuel risk flag"] = true
    settings["Risk aversion flag"] = true
    settings["Risk aversion weight"] = 0.5
    settings["Value-at-Risk percent"] = 0.01
    # Policy instruments
    settings["Price cap"] = 2000
    settings["Investment tax credits flag"] = false
    settings["Production tax credits flag"] = false
    settings["Standard flag"] = false
    settings["Carbon cap flag"] = false
    settings["Carbon tax flag"] = false
    # Policy parameters
    # Investment tax credits
    settings["Endogenous tax credits flag"] = false
    settings["Tax credits"] = [0, 0, 0, 0, 0]
    # Production tax credit
    settings["Production credits"] = [0,0,0,0]
    # CES
    settings["Standard mandate"] = 0.80
    # Carbon tax 
    settings["Carbon tax"] = 0
    # Cap and trade
    settings["Carbon cap"] = 25e3 # Value used if activated
    settings["Contracts"] = false
    settings["Budget Type"] = "System_Expectation"
    settings["Vector Type"] = "random"
    settings["Iterations"] = 6

    # Load data
    inputs = load_input_data(inputs_path, settings)
    results_cap = []
    results_syscost = []
    results_syscost_risk = []
    results_emissions = []
    run_labels = []
    iterations = settings["Iterations"]

    # Run model

    # Build base stochastic model and run for lowest first stage cost solution (risk aversion weight = 0)
    model = build_optimization_model(inputs, settings, "Expectation", 0.75)
    output_exp, model = run_optimization_model(model, inputs, settings)
    # Write results
    results_destination = string(results_path,results_folder,"/EV")
    df_cap, df_syscost, df_syscost_risk, df_emissions = write_results(output_exp, inputs, settings, results_destination)
    push!(results_cap, df_cap)
    push!(results_syscost, df_syscost)
    push!(results_syscost_risk, df_syscost_risk)
    push!(results_emissions, df_emissions)

    @info("Expected value system cost of EV solution: ", output_exp["Expected cost"])

    #model = build_optimization_model(inputs, settings, "CVaR", 0.75)
    set_objective!(model, "Weighted CVaR"; obj_weight=0.75)

    output_cvar, model = run_optimization_model(model, inputs, settings)
    results_destination = string(results_path,results_folder,"/CVaR")
    df_cap, df_syscost, df_syscost_risk, df_emissions = write_results(output_cvar, inputs, settings, results_destination)
    push!(results_cap, df_cap)
    push!(results_syscost, df_syscost)
    push!(results_syscost_risk, df_syscost_risk)
    push!(results_emissions, df_emissions)


    # Set budgets?
    @info("Risk adjusted system cost of EV solution: ", output_exp["Risk adjusted system cost"])
    @info("Risk adjusted system cost of CVaR solution: ", output_cvar["Risk adjusted system cost"])
    
    #add_budget_constraint!(model, 0.5*(output_cvar["Risk adjusted system cost"] + output_exp["Risk adjusted system cost"]) , "System_CVaR"; risk_aversion=0.75)
    add_budget_constraint!(model, (output_cvar["Investment cost"] + output_cvar["Expected cost"]), "System_Expected")
    output_cvar_w_expbudget, model = run_optimization_model(model, inputs, settings)
    results_destination = string(results_path,results_folder,"/CVaR_w_exp_budget")
    df_cap, df_syscost, df_syscost_risk, df_emissions = write_results(output_cvar_w_expbudget, inputs, settings, results_destination)
    push!(results_cap, df_cap)
    push!(results_syscost, df_syscost)
    push!(results_syscost_risk, df_syscost_risk)
    push!(results_emissions, df_emissions)
    add_budget_constraint!(model, (output_cvar["Risk adjusted system cost"]), "System_CVaR"; risk_aversion=0.75)


    vectors = generate_weights(iterations, length(model[:x]), settings["Vector Type"])

    for iteration in 1:iterations
        set_objective!(model, "Capacity"; set_coeffs = vectors[iteration])
        output_random, model = run_optimization_model(model, inputs, settings)
        results_destination = string(results_path,results_folder,"/Random_",iteration)
        df_cap, df_syscost, df_syscost_risk, df_emissions = write_results(output_random, inputs, settings, results_destination)
        push!(results_cap, df_cap)
        push!(results_syscost, df_syscost)
        push!(results_syscost_risk, df_syscost_risk)
        push!(results_emissions, df_emissions)
    end

    write_exploration_results!(results_cap, results_syscost, results_syscost_risk, results_emissions, results_folder)
    
end

run_stochastic_exploration()

