# A Generation Expansion Planning model
# Stochastic version
# Contact: Emil Dimanchev

using DataFrames
using JuMP
using CSV
using Gurobi

# ~~~
# Settings
# ~~~

# Choose stochastic or deterministic
stochastic = true

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

if stochastic
    demand_input = CSV.read(string(inputs_path,sep,"Demand.csv"), DataFrame, header=true)
else
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
if stochastic
    S = size(demand_input[:,2:end])[2] # number of scenarios
end

R = size(resources_input)[1] # number of supply technologies
resources = resources_input[:,1]

# Parameters
cost_inv = resources_input[:,"Investment cost"] # $/MW-336days
cost_var = resources_input[:, "Operating cost"] # $/MWh
availability = Matrix(resource_avail_input[:, 2:end])
co2_factors = resources_input[:,"Emissions_ton_per_MWh"] # ton per MWh
price_cap = 1000 # $/MWh
carbon_cap = 300e3 # tCO2, arbitrary

demand = Array(demand_input[:,2:end])

# Define scenario probabilities
P = ones(S)*1/S # uniform distribution

gep = Model(Gurobi.Optimizer)


# Supply
# Power supply, aka generation but could also be from storage in the future
@variable(gep, g[r in 1:R, t in 1:T, s in 1:S] >= 0) # amount of power supply, MW

# Capacity
@variable(gep, x[r in 1:R] >= 0) # Capacity, MW

# Non-served energy
@variable(gep, nse[t in 1:T, s in 1:S] >= 0)


@objective(gep, Min, sum(x[r]*cost_inv[r] for r in 1:R) + sum(P[s]*g[r,t,s]*cost_var[r] for r in 1:R, t in 1:T, s in 1:S) + sum(P[s]*price_cap*nse[t,s] for t in 1:T, s in 1:S))

@constraint(gep, PowerBalance[t in 1:T, s in 1:S], sum(g[r,t,s] for r in 1:R) + nse[t,s] == demand[t,s])
@constraint(gep, CapacityLimit[r in 1:R, t in 1:T, s in 1:S], g[r,t,s] <= x[r]*availability[t,r])

if co2_cap_flag
    @constraint(gep, emissions_cap[s in 1:S], sum(g[r,t,s]*co2_factors[r] for r in 1:R, t in 1:T) <= carbon_cap)
end

JuMP.optimize!(gep)

# Report Results
obj = objective_value(gep)
gen = value.(g)
nse_all = value.(nse)
cap = value.(x)
co2_price = dual.(emissions_cap) 

# Collect results (the code needs work, currently only compiles results for one scenario)
scenario_for_results = 1

df_cap = DataFrame(Resource = resources, Capacity = cap)
df_gen = DataFrame(transpose(gen[:,:,scenario_for_results]), resources)
insertcols!(df_gen, 1, :Time => time_index)
df_price = DataFrame(Time = time_index, Price = dual.(PowerBalance)[:,scenario_for_results])
df_nse = DataFrame(Time = time_index, Non_served_energy = nse_all[:,scenario_for_results])

co2_tot = zeros(3)
for s in 1:S
    co2_tot[s] = sum(gen[i,t,s]*co2_factors[i] for i in 1:R, t in 1:T)
end


# Write output files (the code below might need work)
# CSV.write(string(results_path,sep,"capacity.csv"), df_cap)
# CSV.write(string(results_path,sep,"generation.csv"), df_gen)
# CSV.write(string(results_path,sep,"price.csv"), df_price)
# CSV.write(string(results_path,sep,"nse.csv"), df_nse)

