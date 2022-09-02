# A Generation Expansion model
# Stochastic version
# Contact: Emil Dimanchev
# ~~~

using DataFrames
using JuMP
using Complementarity
using CSV
using Ipopt
using Random
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
    S_all = size(demand_input[:,2:end])[2] # number of scenarios
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
S = 3 # in sample size
demand_sample = zeros(T,S)

if scenario_selection == "random"
    sample = randperm(S_all)[1:S]   
end
if scenario_selection == "distributed"
    sample = collect(1:S_all)[begin:5:end]
end

for (i,s) in zip(1:S, sample)
    demand_sample[:,i] = demand[:,s]
end
P = ones(S)*1/S # parameters scenario for probabilities

gep = Model(Ipopt.Optimizer)
set_optimizer_attribute(gep, "print_level", 2)
# Supply
# Power supply, aka generation but could also be from storage in the future
@variable(gep, g[r in 1:R, t in 1:T, s in 1:S] >= 0) # amount of power supply, MW

# Capacity
@variable(gep, c[r in 1:R] >= 0) # Capacity, MW

# Non-served energy
@variable(gep, nse[t in 1:T, s in 1:S] >= 0)


@objective(gep, Min, sum(c[r]*inv_cost[r] for r in 1:R) + sum(P[s]*g[r,t,s]*var_cost[r] for r in 1:R, t in 1:T, s in 1:S) + sum(P[s]*price_cap*nse[t,s] for t in 1:T, s in 1:S))

@constraint(gep, PowerBalance[t in 1:T, s in 1:S], sum(g[r,t,s] for r in 1:R) + nse[t,s] == demand_sample[t,s])
@constraint(gep, CapacityLimit[r in 1:R, t in 1:T, s in 1:S], g[r,t,s] <= c[r]*availability[t,r])

if co2_cap_flag
    @constraint(gep, emissions_cap[s in 1:S], sum(g[r,t,s]*co2_factors[r] for r in 1:R, t in 1:T) <= carbon_cap)
end

JuMP.optimize!(gep)

# Report Results
obj = objective_value(gep)
gen = getvalue.(g)
nse_all = getvalue.(nse)
cap = getvalue.(c)

# Collect results
scenario_for_results = 1
df_cap = DataFrame(Resource = resources, Capacity = cap)
df_gen = DataFrame(transpose(gen[:,:,scenario_for_results]), resources)
insertcols!(df_gen, 1, :Time => time_index)
df_price = DataFrame(Time = time_index, Price = getdual.(PowerBalance)[:,scenario_for_results])
df_nse = DataFrame(Time = time_index, Non_served_energy = nse_all[:,scenario_for_results])

co2_tot = zeros(3)
for s in 1:S
    co2_tot[s] = sum(gen[i,t,s]*co2_factors[i] for i in 1:R, t in 1:T)
end
co2_tot
getdual.(emissions_cap)
# df_co2_price = DataFrame(Time = time_index, Price = getdual.(power_balance))

# Write output files
# CSV.write(string(results_path,sep,"capacity.csv"), df_cap)
# CSV.write(string(results_path,sep,"generation.csv"), df_gen)
# CSV.write(string(results_path,sep,"price.csv"), df_price)
# CSV.write(string(results_path,sep,"nse.csv"), df_nse)


# ~~~
# Out of sample analysis
# ~~~

# function generation_model(capacity, demand)
#     # Model with exogenous capacity and deterministic demand
#     gep_test = Model(Ipopt.Optimizer)
#     # Supply
#     # Power supply, aka generation but could also be from storage in the future
#     @variable(gep_test, g[r in 1:R, t in 1:T] >= 0) # amount of power supply, MW

#     # Non-served energy
#     @variable(gep_test, nse[t in 1:T] >= 0)


#     @objective(gep_test, Min, sum(capacity[r]*inv_cost[r] for r in 1:R) + sum(g[r,t]*var_cost[r] for r in 1:R, t in 1:T) + sum(price_cap*nse[t] for t in 1:T))

#     @constraint(gep_test, PowerBalance[t in 1:T], sum(g[r,t] for r in 1:R) + nse[t] == demand[t])
#     @constraint(gep_test, CapacityLimit[r in 1:R, t in 1:T], g[r,t] <= capacity[r]*availability[t,r])

#     JuMP.optimize!(gep_test)

#     return objective_value(gep_test)
# end

# Out of sample testing
# tests = S_all-S
# scenarios_all = collect(1:S_all)
# sample_test = scenarios_all[(!in).(scenarios_all,Ref(sample))]

# obj_test = zeros(tests)

# for (i, s) in zip(1:tests, sample_test)
#     println("Scenario test # ", i)
#     obj_test[i] = generation_model(cap, demand[:,s])
# end
# CSV.write(string(results_path,sep,"out_of_sample_obj.csv"), DataFrame(Objective_values = obj_test))

# mean_error = mean(obj.-obj_test)
# std_error = std(obj.-obj_test)
# println("Average error is ", mean_error)
# println("Standard deviation of error is ", std_error)