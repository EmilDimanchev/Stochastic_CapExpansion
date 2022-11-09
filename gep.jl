# A Generation Expansion Planning model

using DataFrames
using JuMP
using CSV
using Gurobi
import Plots
using Plots
using GR
using PyPlot
#using PlotlyJS
using Pkg


# ~~~
# Settings
# ~~~

# Policy]
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
results_path = string(working_path, sep, "Results", sep, "Results_091222")

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
# Model parameters
# ~~~

# Sets
# Time steps (number of hours)
T = size(demand_input)[1] # number of time steps
# Number of supply technologies
R = size(resources_input)[1] 
resources = resources_input[:,1]
generators = resources[(resources_input[(resources_input[!,:Generation].==1),:][!,:Index_ID])]
# Number of generators
G = length(resources_input[(resources_input[!,:Generation].==1),:][!,:Index_ID])
# Number of storage technologies
S = length(resources_input[(resources_input[!,:Storage].==1),:][!,:Index_ID])

# Parameters
cost_inv = resources_input[:,"Investment cost"] # $/MW-336days
cost_var = resources_input[:, "Operating cost"] # $/MWh
co2_factors = resources_input[:,"Emissions_ton_per_MWh"] # ton per MWh
availability = Matrix(resource_avail_input[:, 2:end])
price_cap = 70 # $/MWh
carbon_cap = 3e3 # tCO2, arbitrary
# Storage
# Current parameters assume battery storage
p_e_ratio = 1/8 # Power to energy ratio
# Single-trip efficiecy
eff_down = 0.9
eff_up = 0.9

demand = Array(demand_input[:,2])

# ~~~
# Build model
# ~~~

gep = Model(Gurobi.Optimizer)

# Variables

# Supply
# Power supply, aka generation but could also be from storage in the future
@variable(gep, g[r in 1:G, t in 1:T] >= 0) # amount of power supply, MW

# Capacity
@variable(gep, x[r in 1:R] >= 0) # Capacity, MW

# Non-served energy
@variable(gep, nse[t in 1:T] >= 0)

# Storage
@variable(gep, discharge[s in 1:S, t in 1:T] >= 0)
@variable(gep, charge[s in 1:S, t in 1:T] >= 0)
# State of charge
@variable(gep, e[s in 1:S, t in 1:T] >= 0)

# Objective 

@objective(gep, Min, sum(x[r]*cost_inv[r] for r in 1:R) + sum(g[r,t]*cost_var[r] for r in (resources_input[(resources_input[!,:Generation].==1),:][!,:Index_ID]), t in 1:T) + sum(price_cap*nse[t] for t in 1:T))

# ~~~
# Constraints
# ~~~

# Energy balance
@constraint(gep, power_balance[t in 1:T], sum(g[r,t] for r in 1:G) + sum(discharge[s,t] - charge[s,t] for s in 1:S) + nse[t] == demand[t])

# Generation capacity limit
@constraint(gep, capacity_limit[r=resources_input[(resources_input[!,:Generation].==1),:][!,:Index_ID], t in 1:T], x[r]*availability[t,r] >= g[r,t])

# Storage constraints
# Loop through storage technologies
for r in (resources_input[(resources_input[!,:Storage].==1),:][!,:Index_ID])
    # Wrap first and last periods
    @constraint(gep, state_of_charge_start[s in 1:S, t in 1:1], e[s,t] == e[s,T] - (1/eff_down)*discharge[s,t] + eff_up*charge[s,t])
    # Energy balance for the remaining periods
    @constraint(gep, state_of_charge[s in 1:S, t in 2:T], e[s,t] == e[t-1] - (1/eff_down)*discharge[s,t] + eff_up*charge[s,t])
    @constraint(gep, energy_limit[s in 1:S, t in 1:T], e[s,t] <= (1/p_e_ratio)*x[r])
    @constraint(gep, charge_limit_total[s in 1:S, t in 1:T], charge[s,t] <= (1/eff_up)*x[r])
    @constraint(gep, charge_limit[s in 1:S, t in 1:T], charge[s,t] <= (1/p_e_ratio)*x[r] - e[s,t])
    @constraint(gep, discharge_limit_total[s in 1:S, t in 1:T], discharge[s,t] <= eff_down*x[r])
    @constraint(gep, discharge_limit[s in 1:S, t in 1:T], discharge[s,t] <= e[s,t])
    @constraint(gep, charge_discharge_balance[s in 1:S, t in 1:T], (1/eff_down)*discharge[s,t] + eff_up*charge[s,t] <= x[r])
end

if co2_cap_flag
    @constraint(gep, emissions_cap, sum(g[r,t]*co2_factors[r] for r in (resources_input[(resources_input[!,:Generation].==1),:][!,:Index_ID]), t in 1:T) <= carbon_cap)
end

JuMP.optimize!(gep)

# ~~~
# Process results
# ~~~

# Report results
gen = value.(g)
c = value.(charge)
d = value.(discharge)
e = value.(e)
maximum(d)
maximum(c)
maximum(e)
cap = value.(x)

nse_all = value.(nse)
m = dual.(capacity_limit)

# Collect results
df_cap = DataFrame(Resource = resources, Capacity = cap)
df_gen = DataFrame(transpose(gen), generators)
insertcols!(df_gen, 1, :Time => time_index)
df_m = DataFrame(transpose(Array(m)), generators)
insertcols!(df_m, 1, :Time => time_index)
df_price = DataFrame(Time = time_index, Price = dual.(power_balance))
df_nse = DataFrame(Time = time_index, Non_served_energy = nse_all)
if co2_cap_flag
    co2_price = dual.(emissions_cap)
    println("CO2 price is ", round(co2_price))
end
# Write output files
CSV.write(string(results_path,sep,"capacity.csv"), df_cap)
CSV.write(string(results_path,sep,"generation.csv"), df_gen)
CSV.write(string(results_path,sep,"price.csv"), df_price)
CSV.write(string(results_path,sep,"nse.csv"), df_nse)
CSV.write(string(results_path,sep,"capacity_rent.csv"), df_m)

co2_tot = sum(gen[i,t]*co2_factors[i] for i in 1:G, t in 1:T)
println("Emissions equal ", round(co2_tot), " tons")
# objective_value(gep)


#Ploting the ouput capacities

display(Plots.plot(df_cap.Resource, df_cap.Capacity,title ="Production capacity",ylabel="Capacity [MWh]" ,seriestype =[:bar], palette = cgrad(:greens), fill=0, alpha=0.6))

#Estimate revenues per tecnology
revenues = zeros(5)#[nuclear, gas, wind, solar, batteries]

for r in 1:R-1
    revenues[r] = sum((df_price[i,2]-cost_var[r])*df_gen[i,r+1] for i in 1:T)
end


#Estimation of revenue

df_rev = DataFrame(Resource = resources, Revenue = revenues)

CSV.write(string(results_path,sep,"revenue.csv"), df_rev)
display(Plots.plot(df_rev.Resource, df_rev.Revenue, title = "Revenue per resource", ylabel = "Revenue [USD]" ,seriestype =[:bar], palette = cgrad(:blues), fill=0, alpha=0.6))


#Estimate and compile the emissions
CO2_price = 70 #[$/tCO2]
CO2_emission_factor = [0,0.4,0,0,0] #[tCO2/MWh] [Nuclear, gas, wind, solar, batteries]
emissions = zeros(T,R)

for j in 1:R-1
    for i in 1:T
    emissions[i,j] = df_gen[i,j+1] * CO2_emission_factor[j]
    end
end

df_emission = DataFrame(emissions,resources)
insertcols!(df_emission, 1, :Time => time_index)

CSV.write(string(results_path,sep,"emissions_per_tec.csv"), df_emission)
display(Plots.plot(collect(1:T), emissions ,title = "Emission per Technology",ylabel = "[Ton CO2]", seriestype =[:bar], palette = cgrad(:reds), fill=0, alpha=0.2, layout = (3,2))) #, label =("Nuclear", "Gas","Wind","Solar","Batteries")

Sum_Emissions = sum(df_emission.Gas + df_emission.Nuclear + df_emission.Wind +df_emission.Solar + df_emission.Batteries)
print("Total Emissions [ton CO2]: ")
display(Sum_Emissions)

#Non served energy, just for fun

display(Plots.plot(df_nse.Time, df_nse.Non_served_energy, xlabel = "Time [h]", ylabel="[MWh]", title ="Non Served Energy", label = "NSE"))


SumNSE = sum(df_nse.Non_served_energy) 
print("Non served energy [MWh]: ")
display(SumNSE)


#Revenue on batteries
#display(transpose(e))
c_trans = transpose(c)
d_trans =transpose(d)

df_storage = DataFrame(Charge = vec(c_trans), Discharge =vec(d_trans))#,Discharge=transpose(d), StateOfCharge = transpose(e))
insertcols!(df_storage,1,:Time => time_index)
#display(df_price)

revenue_storage = zeros(T)
for i in 1:T
    revenue_storage[i] = df_storage.Charge[i]*df_price.Price[i] + df_storage.Discharge[i]*df_price.Price[i]
end
print("Investment cost storage: ",cost_inv[5])
tot_revenue_storage = sum(revenue_storage) - cost_inv[5]
print("Storage Revenue: ", tot_revenue_storage)




#QUESTIONS
# Struggel how i can get values out of the dataframes df_price for a for loop etc. How should I do this
#Is this the right way to define variables
#Is this the type of plot you was thinking of? 
#Shoud we implement this at the end of the script as it is now, or should it be located elsewhere?
#How can I implement the new taxes in the most realistic way(maybe some logic so we can choose to include them or not)? And what are the different values telling us? If NSE grows --> risk increas? economic benificial decreas?




