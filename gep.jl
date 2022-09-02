# A Generation Expansion Planning model

using DataFrames
using JuMP
using CSV
using GLPK
using Ipopt

# ~~~
# Settings
# ~~~

# Policy
# Carbon constraint
co2_cap_flag = false

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
inputs_path = string(working_path, sep, "Inputs", sep, "Inputs_2weeks_storage")
results_path = string(working_path, sep, "Results", sep, "Results_090122")

if !(isdir(results_path))
    mkdir(results_path)
end

# ~~~
# Load data
# ~~~

demand_input = CSV.read(string(inputs_path,sep,"Demand.csv"), DataFrame, header=true, select=["Time_index", "Load"])
resources_input = CSV.read(string(inputs_path,sep,"Resources.csv"), DataFrame, header=true)
resource_avail_input = CSV.read(string(inputs_path,sep,"Resources_availability.csv"), DataFrame, header=true)
time_index = demand_input[:,1]

# ~~~
# Model
# ~~~

# Sets
# Time steps (number of hours)
T = size(demand_input)[1] # number of time steps
# Number of supply technologies
R = size(resources_input)[1] 
resources = resources_input[:,1]
# Number of generators
G = length(resources_input[(resources_input[!,:Generation].==1),:][!,:Index_ID])
# Number of storage technologies
S = length(resources_input[(resources_input[!,:Storage].==1),:][!,:Index_ID])

# Parameters
cost_inv = resources_input[:,"Investment cost"] # $/MW-336days
cost_var = resources_input[:, "Operating cost"] # $/MWh
co2_factors = resources_input[:,"Emissions_ton_per_MWh"] # ton per MWh
availability = Matrix(resource_avail_input[:, 2:end])
price_cap = 15000 # $/MWh
carbon_cap = 300e3 # tCO2, arbitrary
# Storage
# Current parameters assume battery storage
r = 1/8 # Power to energy ratio
# Single-trip efficiecy
eff_down = 0.9
eff_up = 0.9

demand = Array(demand_input[:,2])

gep = Model(Ipopt.Optimizer)
set_optimizer_attribute(gep, "print_level", 2)

# Variables

# Supply
# Power supply, aka generation but could also be from storage in the future
@variable(gep, g[r in 1:G, t in 1:T] >= 0) # amount of power supply, MW

# Capacity
@variable(gep, x[r in 1:R] >= 0) # Capacity, MW

# Non-served energy
@variable(gep, nse[t in 1:T] >= 0)

# Storage
@variable(gep, discharge[t in 1:T] >= 0)
@variable(gep, charge[t in 1:T] >= 0)
# State of charge
@variable(gep, e[t in 1:T] >= 0)

# Objective 

@objective(gep, Min, sum(x[r]*cost_inv[r] for r in 1:R) + sum(g[r,t]*cost_var[r] for r in (resources_input[(resources_input[!,:Generation].==1),:][!,:Index_ID]), t in 1:T) + sum(price_cap*nse[t] for t in 1:T))

# ~~~
# Constraints
# ~~~

# Energy balance
@constraint(gep, power_balance[t in 1:T], sum(g[r,t] for r in 1:G) + discharge[t] - charge[t] + nse[t] == demand[t])

# Generation capacity limit
@constraint(gep, capacity_limit[r=resources_input[(resources_input[!,:Generation].==1),:][!,:Index_ID], t in 1:T], x[r]*availability[t,r] >= g[r,t])

# Storage constraints
# Loop through storage technologies
for r in (resources_input[(resources_input[!,:Storage].==1),:][!,:Index_ID])
    # Wrap first and last periods
    @constraint(gep, state_of_charge_start[t in 1:1], e[t] == e[T] - (1/eff_down)*discharge[t] + eff_up*charge[t])
    # Energy balance for the remaining periods
    @constraint(gep, state_of_charge[t in 2:T], e[t] == e[t-1] - (1/eff_down)*discharge[t] + eff_up*charge[t])
    @constraint(gep, energy_limit[t in 1:T], e[t] <= (1/r)*x[r])
    @constraint(gep, charge_limit_total[t in 1:T], charge[t] <= (1/eff_up)*x[r])
    @constraint(gep, charge_limit[t in 1:T], charge[t] <= (1/r)*x[r] - e[t])
    @constraint(gep, discharge_limit_total[t in 1:T], discharge[t] <= eff_down*x[r])
    @constraint(gep, discharge_limit[t in 1:T], discharge[t] <= e[t])
    @constraint(gep, charge_discharge_balance[t in 1:T], (1/eff_down)*discharge[t] + eff_up*charge[t] <= x[r])
end

if co2_cap_flag
    @constraint(gep, emissions_cap, sum(g[r,t]*co2_factors[r] for r in (resources_input[(resources_input[!,:Generation].==1),:][!,:Index_ID]), t in 1:T) <= carbon_cap)
end

JuMP.optimize!(gep)
demand
# Report Results
gen = value.(g)
c = value.(charge)
d = value.(discharge)
e = value.(e)
maximum(d)
maximum(c)
maximum(e)
cap = value.(x)
dual.(gep[:state_of_charge])

# nse_all = value.(nse)
# m = dual.(capacity_limit)
# # Collect results
# df_cap = DataFrame(Resource = resources, Capacity = cap)
# df_gen = DataFrame(transpose(gen), resources)
# insertcols!(df_gen, 1, :Time => time_index)
# df_m = DataFrame(transpose(m), resources)
# insertcols!(df_m, 1, :Time => time_index)
# df_price = DataFrame(Time = time_index, Price = dual.(power_balance))
# df_nse = DataFrame(Time = time_index, Non_served_energy = nse_all)
# dual.(emissions_cap)
# # Write output files
# CSV.write(string(results_path,sep,"capacity.csv"), df_cap)
# CSV.write(string(results_path,sep,"generation.csv"), df_gen)
# CSV.write(string(results_path,sep,"price.csv"), df_price)
# CSV.write(string(results_path,sep,"nse.csv"), df_nse)
# CSV.write(string(results_path,sep,"capacity_rent.csv"), df_m)

# co2_tot = sum(gen[i,t]*co2_factors[i] for i in 1:R, t in 1:T)
# println("Emissions equal ", round(co2_tot), " tons")
# objective_value(gep)