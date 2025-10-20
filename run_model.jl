# Run equilibrium model
using CSV
using DataFrames

include("./Write_results.jl")
include("./Load_inputs.jl")
include("./Stochastic_expansion.jl")

# ~~~
# Define paths
# ~~~

inputs_folder = "Inputs_30days_5techs"

results_folder = "Results"

inputs_path = string("./Inputs/",inputs_folder)

results_path = "/Users/ed0400/Modeling/Stochastic_optimization/"

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

output = optimization_model(inputs, settings)

# Write results
results_destination = string(results_path,results_folder)
write_results(output, inputs, settings, results_destination)
