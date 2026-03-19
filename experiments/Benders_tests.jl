
include("../src/Stochastic_CapExpansion.jl")

using .Stochastic_CapExpansion
using Revise, JuMP, Gurobi, DataFrames, CSV, YAML, Random, LinearAlgebra, Combinatorics, Dates, Distributed, DistributedArrays

INPUTS_PATH = "/Users/mike_/Documents/Code/Stochastic_CapExpansion/inputs/Inputs_30repdays_ext_1000scen_7techs"

function benders_test_compare()
    # Load inputs and settings
    
    settings = load_settings(INPUTS_PATH)
    inputs = load_input_data(INPUTS_PATH, settings)
    settings["Results folder path"] = "/Users/mike_/Documents/Code/Stochastic_CapExpansion/outputs/Benders_test_base"
    if !isdir(settings["Results folder path"])
        mkdir(settings["Results folder path"])
    end
    names = inputs["Generation resources"]
    
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
    names = inputs["Generation resources"]
    
    # Run Old Benders algorithm
    results = run_benders_algorithm(inputs, settings)
    CSV.write(joinpath(settings["Results folder path"], "capacity_mix_per_iteration_old.csv"), DataFrame(stack(results["Capacity per iteration"], dims=1), names))

end

function run_new_benders()
    # Load inputs and settings
    
    settings = load_settings(INPUTS_PATH)
    inputs = load_input_data(INPUTS_PATH, settings)
    settings["Results folder path"] = "/Users/mike_/Documents/Code/Stochastic_CapExpansion/outputs/Benders_test_base"
    if !isdir(settings["Results folder path"])
        mkdir(settings["Results folder path"])
    end
    names = inputs["Generation resources"]
    
    # Build models
    MP = build_planning_model(inputs, settings)
    SPs = build_all_subproblems(inputs, settings)
    # Run New Benders algorithm
    results_new = benders_algorithm(inputs, settings, MP, SPs)
    CSV.write(joinpath(settings["Results folder path"], "capacity_mix_per_iteration_new.csv"), DataFrame(stack(results_new["Capacity per iteration"], dims=1), names))
end

function run_new_benders_distributed()
    # Load inputs and settings
    
    settings = load_settings(INPUTS_PATH)
    inputs = load_input_data(INPUTS_PATH, settings)
    addprocs(settings["Workers"])
    @everywhere begin
        include("../src/Stochastic_CapExpansion.jl")

        using .Stochastic_CapExpansion
        using Revise, JuMP, Gurobi, DataFrames, CSV, YAML, Random, LinearAlgebra, Combinatorics, Dates

        INPUTS_PATH = "/Users/mike_/Documents/Code/Stochastic_CapExpansion/inputs/Inputs_30repdays_ext_1000scen_7techs"
        settings = load_settings(INPUTS_PATH)
        inputs = load_input_data(INPUTS_PATH, settings)
    end


    
    settings["Results folder path"] = "/Users/mike_/Documents/Code/Stochastic_CapExpansion/outputs/Benders_test_base"
    if !isdir(settings["Results folder path"])
        mkdir(settings["Results folder path"])
    end
    names = inputs["Generation resources"]
    
    # Build models
    MP = build_planning_model(inputs, settings)
    SPs = build_all_subproblems(inputs, settings)
    # Run New Benders algorithm
    results_new = benders_algorithm(inputs, settings, MP, SPs)
    CSV.write(joinpath(settings["Results folder path"], "capacity_mix_per_iteration_new.csv"), DataFrame(stack(results_new["Capacity per iteration"], dims=1), names))
end

#benders_test_compare()
run_new_benders()
#run_old_benders()