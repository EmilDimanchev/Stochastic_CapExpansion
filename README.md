## Running Instructions:
After cloning follow these running instructions.
Always cd into this directory. <cd Stochastic_CapExpansion>

# Instantiation (First Time Only):
1. open julia - `julia --project=.`
2. move to Pkg environment for project - `]`
3. instantiate - `instantiate`

# Run case:
1. open julia with set number of threads for shared memory parallelism- `julia --project=. -t8`
2. include test-file - `include("experiments/MGA_Tests.jl")`
3. run desired function in terminal - `simple_comp()`


Please note that all test files must have their first lines as:

`include("Stochastic_CapExpansion.jl")`

`using .Stochastic_CapExpansion`
