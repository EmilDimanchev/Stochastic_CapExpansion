using CSV, DataFrames

include("./Write_results.jl")
include("./Load_inputs.jl")
include("./Stochastic_expansion.jl")

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

    # Load data
    inputs = load_input_data(inputs_path, settings)
    # Run model

    # Build base stochastic model and run for lowest first stage cost solution (risk aversion weight = 0)
    model = build_optimization_model(inputs, settings, "Expectation", 0.5)
    output_exp, model = run_optimization_model(model, inputs, settings)
    # Write results
    results_destination = string(results_path,results_folder,"/EV")
    write_results(output_exp, inputs, settings, results_destination)

    # Max Risk Aversion
    set_objective!(model, "Weighted CVaR"; obj_weight=1.0)
    output_cvar, model = run_optimization_model(model, inputs, settings)
    results_destination = string(results_path,results_folder,"/CVaR_max")
    write_results(output_cvar, inputs, settings, results_destination)


    # Set budgets?
    println("Risk adjusted system cost of EV solution: ", output_exp["Risk adjusted system cost"])
    println("Risk adjusted system cost of CVaR solution: ", output_cvar["Risk adjusted system cost"])
    println("Setting budget between the two solutions and running random exploration. Budget: ", 0.5*(output_cvar["Risk adjusted system cost"] + output_exp["Risk adjusted system cost"]))
    add_budget_constraint!(model, 0.5*(output_cvar["Risk adjusted system cost"] + output_exp["Risk adjusted system cost"]) , "System_CVaR"; risk_aversion=0.75)
    add_budget_constraint!(model, output_cvar["Investment cost"], "Investment"; risk_aversion=0.0)

    # add nse cap constraint
    println(sum(output_exp["Load shedding"][:,:,:]))
    println("Load shedding in CVaR solution: ", sum(output_cvar["Load shedding"][:,:,:]))
    min_nse = min(sum(output_exp["Load shedding"][:,:,:]), sum(output_cvar["Load shedding"][:,:,:]))
    println("Setting NSE cap at ", min_nse, " and running random exploration.")
    add_nse_cap!(model, min_nse)

    for iteration in 1:5
        set_objective!(model, "Capacity_Random")
        output_random, model = run_optimization_model(model, inputs, settings)
        results_destination = string(results_path,results_folder,"/Random_",iteration)
        write_results(output_random, inputs, settings, results_destination)
    end


end

run_stochastic_exploration()

