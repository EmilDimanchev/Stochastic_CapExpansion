# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define paths
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


inputs_folder = "Inputs_30repdays_ext_1000scen_7techs"


results_folder = "results"

inputs_path = string("./Inputs/",inputs_folder)

results_path = string("./",results_folder)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set up experiment
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

settings = Dict{String, Any}()
settings["Results folder path"] = results_path
settings["Scaling factor cost"] = 1e5
settings["Scaling factor demand"] = 1
settings["Convergence tolerance"] = 1e-5

settings["Cost of nonserved energy"] = [9e3,1500,1100,400]
settings["Max nse"] = [0.917,0.05,0.03,0.003]

# Uncertainty representation
settings["Demand uncertainty"] = true
settings["Gas price uncertainty"] = true
settings["Weather uncertainty"] = false

# Risk preferences
settings["Risk aversion flag"] = true
settings["Risk aversion weight"] = 0.5
settings["Value-at-Risk percent"] = 0.05


export settings, inputs_path