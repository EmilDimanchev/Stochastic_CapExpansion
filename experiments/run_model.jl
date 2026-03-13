# Run equilibrium model
include("../src/Stochastic_CapExpansion.jl")
using .Stochastic_CapExpansion

# ~~~
# Define paths
# ~~~

inputs_folder = "Inputs_30days_5techs"

results_folder = "Results"

inputs_path = string("./Inputs/",inputs_folder)

results_path = "./"

# ~~~
# Set up experiment
# ~~~

# Model settings
settings = load_settings(inputs_path)

# Load data
inputs = load_input_data(inputs_path, settings)
# Run model

output = optimization_model(inputs, settings)


# Write results
results_destination = string(results_path,results_folder)
write_results(output, inputs, settings, results_destination)
