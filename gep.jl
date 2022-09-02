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
co2_cap_flag = true

# Uncertainty handling: "deterministic" or "stochastic"
# NOTE: stochastic is implemented in another script
uncertainty = "deterministic"

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
results_path = string(working_path, sep, "Results", sep, "Results_031422")

if !(isdir(results_path))
    mkdir(results_path)
end

# ~~~
# Load data
# ~~~

if uncertainty == "deterministic"
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

R = size(resources_input)[1] # number of supply technologies
resources = resources_input[:,1]

# Parameters
cost_inv = resources_input[:,"Investment cost"] # $/MW-336days
cost_var = resources_input[:, "Operating cost"] # $/MWh
co2_factors = resources_input[:,"Emissions_ton_per_MWh"] # ton per MWh
availability = Matrix(resource_avail_input[:, 2:end])
price_cap = 1000 # $/MWh
carbon_cap = 300e3 # tCO2, arbitrary

demand = Array(demand_input[:,2])

gep = Model(Ipopt.Optimizer)
set_optimizer_attribute(gep, "print_level", 2)

# Supply
# Power supply, aka generation but could also be from storage in the future
@variable(gep, g[r in 1:R, t in 1:T] >= 0) # amount of power supply, MW

# Capacity
@variable(gep, x[r in 1:R] >= 0) # Capacity, MW

# Non-served energy
@variable(gep, nse[t in 1:T] >= 0)

@objective(gep, Min, sum(x[r]*cost_inv[r] for r in 1:R) + sum(g[r,t]*cost_var[r] for r in 1:R for t in 1:T) + sum(price_cap*nse[t] for t in 1:T))

@constraint(gep, power_balance[t in 1:T], sum(g[r,t] for r in 1:R) + nse[t] == demand[t])
@constraint(gep, capacity_limit[r in 1:R, t in 1:T], x[r]*availability[t,r] >= g[r,t])

if co2_cap_flag
    @constraint(gep, emissions_cap, sum(g[r,t]*co2_factors[r] for r in 1:R, t in 1:T) <= carbon_cap)
end

JuMP.optimize!(gep)

# Report Results
gen = getvalue.(g)
cap = getvalue.(x)
nse_all = getvalue.(nse)
m = getdual.(capacity_limit)
# Collect results
df_cap = DataFrame(Resource = resources, Capacity = cap)
df_gen = DataFrame(transpose(gen), resources)
insertcols!(df_gen, 1, :Time => time_index)
df_m = DataFrame(transpose(m), resources)
insertcols!(df_m, 1, :Time => time_index)
df_price = DataFrame(Time = time_index, Price = getdual.(power_balance))
df_nse = DataFrame(Time = time_index, Non_served_energy = nse_all)
getdual.(emissions_cap)
# Write output files
CSV.write(string(results_path,sep,"capacity.csv"), df_cap)
CSV.write(string(results_path,sep,"generation.csv"), df_gen)
CSV.write(string(results_path,sep,"price.csv"), df_price)
CSV.write(string(results_path,sep,"nse.csv"), df_nse)
CSV.write(string(results_path,sep,"capacity_rent.csv"), df_m)

co2_tot = sum(gen[i,t]*co2_factors[i] for i in 1:R, t in 1:T)
println("Emissions equal ", round(co2_tot), " tons")
objective_value(gep)