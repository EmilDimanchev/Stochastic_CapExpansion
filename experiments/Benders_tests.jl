
include("../src/Stochastic_CapExpansion.jl")

using .Stochastic_CapExpansion
using Revise, JuMP, Gurobi, DataFrames, CSV, YAML, Random, LinearAlgebra, Combinatorics, Dates

INPUTS_PATH = "/Users/mike_/Documents/Code/Stochastic_CapExpansion/inputs/Inputs_30repdays_ext_1000scen_7techs"

function benders_test_compare()
    # Load inputs and settings
    
    settings = load_settings(INPUTS_PATH)
    inputs = load_input_data(INPUTS_PATH, settings)
    settings["Results folder path"] = "/Users/mike_/Documents/Code/Stochastic_CapExpansion/outputs/Benders_test_base"
    if !isdir(settings["Results folder path"])
        mkdir(settings["Results folder path"])
    end
    names = inputs["Resources"]
    
    # Run Old Benders algorithm
    results = run_benders_algorithm(inputs, settings)
    CSV.write(joinpath(settings["Results folder path"], "capacity_mix_per_iteration_old.csv"), DataFrame(stack(results["Capacity per iteration"], dims=1), names))

    # Build models
    MP = build_planning_model(inputs, settings)
    SPs = build_all_subproblems(inputs, settings)
    # Run New Benders algorithm
    results_new = benders_algorithm(inputs, settings, MP, SPs)
    CSV.write(joinpath(settings["Results folder path"], "capacity_mix_per_iteration_new.csv"), DataFrame(stack(results_new["Capacity per iteration"], dims=1), names))

end

function run_old_benders()
    # Load inputs and settings
    
    settings = load_settings(INPUTS_PATH)
    inputs = load_input_data(INPUTS_PATH, settings)
    settings["Results folder path"] = "/Users/mike_/Documents/Code/Stochastic_CapExpansion/outputs/Benders_test_base"
    if !isdir(settings["Results folder path"])
        mkdir(settings["Results folder path"])
    end
    names = inputs["Resources"]
    
    # Run Old Benders algorithm
    results = run_benders_algorithm(inputs, settings)
    CSV.write(joinpath(settings["Results folder path"], "capacity_mix_per_iteration_old.csv"), DataFrame(stack(results["Capacity per iteration"], dims=1), names))

end

function run_new_benders()
    # Load inputs and settings
    
    settings = load_settings(INPUTS_PATH)
    inputs = load_input_data(INPUTS_PATH, settings)
    inputs["Output Demand scenario probabilities"] = inputs["Demand scenario probabilities"]
    inputs["Output Fuel price scenario probabilities"] = inputs["Gas price scenario probabilities"]
    inputs["Output Weather scenario probabilities"] = inputs["Weather scenario probabilities"]
    if settings["Parallel flag"]
        settings["Threads"] = Threads.nthreads()
        @info("Running with parallelization using $(Threads.nthreads()) threads")
    else
        @info("Running without parallelization")
    end
    
    settings["Results folder path"] = "/Users/mike_/Documents/Code/Stochastic_CapExpansion/outputs/Benders_test_base"
    if !isdir(settings["Results folder path"])
        mkdir(settings["Results folder path"])
    end
    
    # Build models
    MP = build_planning_model(inputs, settings)
    SPs = build_all_subproblems(inputs, settings)
    # Run New Benders algorithm
    results_new = benders_algorithm(inputs, settings, MP, SPs)
    
    df_cap, df_syscost, df_emissions = write_results_benders(results_new, inputs, settings, settings["Results folder path"])
    write_exploration_results!([df_cap], [df_syscost], [df_emissions], settings["Results folder path"], ["Benders test new"], scenario = "Benders_test_new")
    
end

function run_single_scenario_benders()
    # Load inputs and settings
    
    settings = load_settings(INPUTS_PATH)
    inputs = load_input_data(INPUTS_PATH, settings)
    settings["Results folder path"] = "/Users/mike_/Documents/Code/Stochastic_CapExpansion/outputs/Benders_test_base"
    if !isdir(settings["Results folder path"])
        mkdir(settings["Results folder path"])
    end
    names = inputs["Resources"]
    
    # Build models
    
    eval_SPs = build_all_subproblems(inputs, settings)
    
    # set all uncertainties to false and risk aversion to false for this test, since we are just running with one scenario selected
    settings["Risk aversion flag"] = false
    settings["Demand uncertainty"] = false
    settings["Gas price uncertainty"] = false
    settings["Weather uncertainty"] = false
    new_inputs = load_input_data(INPUTS_PATH, settings) # need to reload inputs after changing settings to turn off uncertainties

    MP = build_planning_model(new_inputs, settings)
    SPs = build_all_subproblems(new_inputs, settings) # need to rebuild SPs after changing settings to turn off uncertainties
    new_inputs["Full Demand scenario probabilities"] = inputs["Demand scenario probabilities"]
    new_inputs["Full Gas price scenario probabilities"] = inputs["Gas price scenario probabilities"]
    new_inputs["Full Weather scenario probabilities"] = inputs["Weather scenario probabilities"]

    # Run New Benders algorithm
    results_new, gaps = benders_algorithm(new_inputs, settings, MP, SPs; Eval_SPs=eval_SPs)
    CSV.write(joinpath(settings["Results folder path"], "capacity_mix_per_iteration_single_scenario.csv"), DataFrame(stack(results_new["Capacity per iteration"], dims=1), names))
    return results_new
end

#benders_test_compare()
#run_new_benders()
#run_old_benders()
