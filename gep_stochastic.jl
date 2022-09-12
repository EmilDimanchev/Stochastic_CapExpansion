# A Generation Expansion model
# Stochastic version
# Contact: Emil Dimanchev
# ~~~

using DataFrames
using JuMP
using Complementarity
using CSV
using Gurobi

using Statistics

# ~~~
# Settings
# ~~~

# Uncertainty handling: "deterministic" or "stochastic"
uncertainty = "stochastic"
# Uncertainty handling: "random" or "distributed"
scenario_selection = "random"

# Policy
# Carbon constraint
co2_cap_flag = true


# ~~~
# Folder paths
# ~~~

if Sys.isunix()
    sep = "/"
elseif Sys.iswindows()
    sep = "\U005c"
end

working_path = pwd()

# Define input and output paths
inputs_path = string(working_path, sep, "Inputs", sep, "Inputs_course_3techs")
results_path = string(working_path, sep, "Results", sep, "Results_stochastic_031422")

if !(isdir(results_path))
    mkdir(results_path)
end

# ~~~
# Load data
# ~~~

if uncertainty == "stochastic"
    demand_input = CSV.read(string(inputs_path,sep,"Demand.csv"), DataFrame, header=true)
elseif uncertainty == "deterministic"
    demand_input = CSV.read(string(inputs_path,sep,"Demand.csv"), DataFrame, header=true, select=["Time_index", "Load"])
end

resources_input = CSV.read(string(inputs_path,sep,"Resources.csv"), DataFrame, header=true)
resource_avail_input = CSV.read(string(inputs_path,sep,"Resources_availability.csv"), DataFrame, header=true)
time_index = demand_input[:,1]

# ~~~
# Model
# ~~~

# Sets
T = size(demand_input)[1] # number of time steps
if uncertainty == "stochastic"
    W_all = size(demand_input[:,2:end])[2] # number of scenarios
end

R = size(resources_input)[1] # number of supply technologies
resources = resources_input[:,1]

# Parameters
inv_cost = resources_input[:,"Investment cost"] # $/MW-336days
var_cost = resources_input[:, "Operating cost"] # $/MWh
availability = Matrix(resource_avail_input[:, 2:end])
co2_factors = resources_input[:,"Emissions_ton_per_MWh"] # ton per MWh
price_cap = 1000 # $/MWh
carbon_cap = 300e3 # tCO2, arbitrary

demand = Array(demand_input[:,2:end])

# Select scenarios
W = 3 # in sample size
demand_sample = zeros(T,W)

if scenario_selection == "random"
    sample = randperm(W_all)[1:W]   
end
if scenario_selection == "distributed"
    sample = collect(1:W_all)[begin:5:end]
end

for (i,s) in zip(1:W, sample)
    demand_sample[:,i] = demand[:,s]
end
P = ones(W)*1/W # parameters scenario for probabilities

gep = Model(Gurobi.Optimizer)

# Supply
# Power supply, aka generation but could also be from storage in the future
@variable(gep, g[r in 1:R, t in 1:T, s in 1:W] >= 0) # amount of power supply, MW

# Capacity
@variable(gep, x[r in 1:R] >= 0) # Capacity, MW

# Non-served energy
@variable(gep, nse[t in 1:T, s in 1:W] >= 0)


@objective(gep, Min, sum(x[r]*inv_cost[r] for r in 1:R) + sum(P[s]*g[r,t,s]*var_cost[r] for r in 1:R, t in 1:T, s in 1:W) + sum(P[s]*price_cap*nse[t,s] for t in 1:T, s in 1:W))

@constraint(gep, PowerBalance[t in 1:T, s in 1:W], sum(g[r,t,s] for r in 1:R) + nse[t,s] == demand_sample[t,s])
@constraint(gep, CapacityLimit[r in 1:R, t in 1:T, s in 1:W], g[r,t,s] <= x[r]*availability[t,r])

if co2_cap_flag
    @constraint(gep, emissions_cap[s in 1:W], sum(g[r,t,s]*co2_factors[r] for r in 1:R, t in 1:T) <= carbon_cap)
end

JuMP.optimize!(gep)

# ~~~
# Process results
# ~~~

# Report Results
obj = objective_value(gep)
gen = value.(g)
nse_all = value.(nse)
cap = value.(x)

# Collect results
scenario_for_results = 1
df_cap = DataFrame(Resource = resources, Capacity = cap)
df_gen = DataFrame(transpose(gen[:,:,scenario_for_results]), resources)
insertcols!(df_gen, 1, :Time => time_index)
df_price = DataFrame(Time = time_index, Price = dual.(PowerBalance)[:,scenario_for_results])
df_nse = DataFrame(Time = time_index, Non_served_energy = nse_all[:,scenario_for_results])

co2_tot = zeros(3)
for s in 1:W
    co2_tot[s] = sum(gen[i,t,s]*co2_factors[i] for i in 1:R, t in 1:T)
end
co2_tot
dual.(emissions_cap)

# ~~~
# Write output files
# ~~~

# CSV.write(string(results_path,sep,"capacity.csv"), df_cap)
# CSV.write(string(results_path,sep,"generation.csv"), df_gen)
# CSV.write(string(results_path,sep,"price.csv"), df_price)
# CSV.write(string(results_path,sep,"nse.csv"), df_nse)
