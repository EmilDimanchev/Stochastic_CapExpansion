## Running Instructions:
After cloning follow these running instructions.
Always cd into this directory. <cd Stochastic_CapExpansion>

# Instantiation (First Time Only):
1. open julia <julia --project=.>
2. move to Pkg environment for project <]>
3. instantiate - <instantiate>

# Run case:
1. open julia <julia --project=.>
2. include test-file <include("experiments/test-file.jl")>

Please note that all test files must have their first lines as:
include("Stochastic_CapExpansion.jl")
using .Stochastic_CapExpansion
