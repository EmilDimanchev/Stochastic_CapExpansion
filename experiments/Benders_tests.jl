
include("../src/Stochastic_CapExpansion.jl")

using .Stochastic_CapExpansion
using Revise, JuMP, Gurobi, DataFrames, CSV, YAML, Random, LinearAlgebra, Combinatorics, Dates

INPUTS_PATH = "/Users/mike_/Documents/Code/Stochastic_CapExpansion/inputs/Inputs_30repdays_ext_1000scen_7techs"

function benders_test()
    # Load inputs and settings
    
    settings = load_settings(INPUTS_PATH)
    inputs = load_input_data(INPUTS_PATH, settings)
    settings["Results folder path"] = "/Users/mike_/Documents/Code/Stochastic_CapExpansion/outputs/Benders_test_base"

    # Run Benders algorithm
    results = run_benders_algorithm(inputs, settings)

    return results
end

benders_test()
