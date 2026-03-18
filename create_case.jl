# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define paths
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# inputs_folder = "Inputs_4repdays_4scen"
# inputs_folder = "Inputs_30repdays_ext_4scen"
# inputs_folder = "Inputs_30repdays_ext_100scen"
inputs_folder = "Inputs_30repdays_ext_1000scen"
# inputs_folder = "Inputs_30repdays_ext_1000scen_7techs"
# inputs_folder = "Inputs_mays"

results_folder = "112624_raopt_demandfuel_test"

inputs_path = string("./Inputs/",inputs_folder)

results_path = "/Users/ed0400/MIT Dropbox/Emil Dimanchev/1_Projects/1. Postdoc/1. Paper on decomposed risk model/Results/"


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set up experiment
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

settings = Dict{String, Any}()
settings["Results folder path"] = results_path
settings["Scaling factor cost"] = 1e5
settings["Scaling factor demand"] = 1
settings["INV SD tolerance"] = 1e-5 # Can help with numerical issues
settings["Convergence tolerance"] = 1e-5
# settings["Cost of nonserved energy"] = [2000,] # $/MWh, the highest price is the price cap
# settings["Max nse"] = [1,]
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

# Contracts
settings["Contracts"] = false
settings["Consumers contracting"] = 1
settings["Contract volume limit flag"] = false
settings["Volume limit"] = 50e3
# If modeling common contracts
settings["Call options flag"] = true
settings["Futures flag"] = true
settings["Wind futures flag"] = false
settings["Solar futures flag"] = false
settings["Futures strike"] = [10.0,15.0, 30, 50]
settings["Wind futures strike"] = [35]
settings["Solar futures strike"] = [35]


export settings, inputs_path