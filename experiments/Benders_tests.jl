include("../src/Stochastic_CapExpansion.jl")
include("MGA_Tests.jl")

using .Stochastic_CapExpansion
using Revise, JuMP, Gurobi, DataFrames, CSV, YAML, Random, LinearAlgebra, Combinatorics, Dates

INPUTS_PATH = "../Inputs/Inputs_30repdays_ext_1000scen_7techs"

function benders_test_compare()
    # Load inputs and settings
    
    settings = load_settings(INPUTS_PATH)
    inputs = load_input_data(INPUTS_PATH, settings)
    configure_parallel_workers!(settings)
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
    configure_parallel_workers!(settings)
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
    configure_parallel_workers!(settings)
    inputs["Output Demand scenario probabilities"] = inputs["Demand scenario probabilities"]
    inputs["Output Fuel price scenario probabilities"] = inputs["Gas price scenario probabilities"]
    inputs["Output Weather scenario probabilities"] = inputs["Weather scenario probabilities"]
    
    settings["Results folder path"] = "/Users/mike_/Documents/Code/Stochastic_CapExpansion/outputs/Benders_test_base"
    if !isdir(settings["Results folder path"])
        mkdir(settings["Results folder path"])
    end
    
    # Build models
    MP = build_planning_model(inputs, settings)
    SPs = build_all_subproblems(inputs, settings)
    # Run New Benders algorithm
    results_new = benders_algorithm(inputs, settings, MP, SPs, "Benders_test_new")
    
    df_cap, df_syscost, df_emissions = write_results_benders(results_new, inputs, settings, settings["Results folder path"])
    write_exploration_results!([df_cap], [df_syscost], [df_emissions], settings["Results folder path"], ["Benders test new"], scenario = "Benders_test_new")
    
end

function run_single_scenario_benders()
    # Load inputs and settings
    
    settings = load_settings(INPUTS_PATH)
    inputs = load_input_data(INPUTS_PATH, settings)
    configure_parallel_workers!(settings)
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
    results_new, gaps = benders_algorithm(new_inputs, settings, MP, SPs, "Benders_test_single_scenario"; Eval_SPs=eval_SPs)
    CSV.write(joinpath(settings["Results folder path"], "capacity_mix_per_iteration_single_scenario.csv"), DataFrame(stack(results_new["Capacity per iteration"], dims=1), names))
    return results_new
end

function run_eval_SP_validation()
    inputs_folder = joinpath("inputs","Inputs_30repdays_ext_1000scen_7techs")
    test_index = 5
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
    
    inputs["Output Demand scenario probabilities"] = inputs["Demand scenario probabilities"] #Establish base weights ahead of time
    inputs["Output Gas price scenario probabilities"] = inputs["Gas price scenario probabilities"]
    inputs["Output Weather scenario probabilities"] = inputs["Weather scenario probabilities"]
    inputs["Full Demand scenario probabilities"] = inputs["Demand scenario probabilities"]
    inputs["Full Gas price scenario probabilities"] = inputs["Gas price scenario probabilities"]
    inputs["Full Weather scenario probabilities"] = inputs["Weather scenario probabilities"]

    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)

    outputs_exp, vectors = run_stochastic_exploration_single_type(SPs, inputs, settings, joinpath(results_folder, "Expected"), summary_folder; type = "System_Expected", vector_set = nothing, budget_multiplier=1.10, Eval_SPs = SPs)

end

function run_risk_pareto_frontier(test_index)
    inputs_folder = joinpath("inputs","Inputs_30repdays_ext_1000scen_7techs_Della")

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

    inputs["Output Demand scenario probabilities"] = inputs["Demand scenario probabilities"] #Establish base weights ahead of time
    inputs["Output Gas price scenario probabilities"] = inputs["Gas price scenario probabilities"]
    inputs["Output Weather scenario probabilities"] = inputs["Weather scenario probabilities"]

    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)
   

    risk_aversion_weights = [0.0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 1.0]

    results_cap = []
    results_syscost = []
    results_emissions = []
    outputs = []
    run_labels = []
    gaps = []
    settings["Risk aversion flag"] = true
    MP = build_planning_model(inputs, settings)
    for risk_aversion_weight in risk_aversion_weights
        settings["Risk aversion weight"] = risk_aversion_weight
        set_objective_bendersMP!(MP, "System_Weighted_CVaR", inputs, settings; obj_weight=risk_aversion_weight)
        output = benders_algorithm(inputs, settings, MP, SPs, "Benders_test_risk_aversion_weight_"*string(risk_aversion_weight))
        gap = output["Gaps"]
        df_cap, df_syscost, df_emissions = write_results_benders(output, inputs, settings, joinpath(results_folder, "Benders_risk_aversion_weight_"*string(risk_aversion_weight)))
        push!(results_cap, df_cap)
        push!(results_syscost, df_syscost)
        push!(results_emissions, df_emissions)
        push!(outputs, output)
        push!(run_labels, "Risk aversion weight "*string(risk_aversion_weight))
        push!(gaps, gap)
    end

    write_exploration_results!(results_cap, results_syscost, results_emissions, joinpath(results_folder, "Summary"), run_labels,"Risk aversion weights")


end

#benders_test_compare()
#run_new_benders()
#run_old_benders()
