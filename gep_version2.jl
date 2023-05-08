

using DataFrames
using JuMP
using CSV
using Gurobi
using Random
using Statistics
import Plots
using Plots; theme(:bright)
using RDatasets
using GR
using PyPlot
using Pkg

# ~~~
# Settings
# ~~~∏

# Choose stochastic or deterministic
stochastic = true
risk_aversion = false
# Policy
# Carbon constraint
co2_cap_flag = true
#CO2-tax policy
CO2_tax_flag = false

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
inputs_path = string(working_path, sep, "Inputs", sep, "Inputs_course_all_techs_annual_new_version") #"Inputs_course_3techs", "Inputs_course_all_techs_annual_2041"
results_path = string(working_path, sep, "Results", sep, "Results_stochastic_new_version")

if !(isdir(results_path))
    mkdir(results_path)
end

# ~~~
# Load data
# ~~~

if stochastic
    demand_input_z1 = CSV.read(string(inputs_path,sep,"Demand_z1.csv"), DataFrame, header=true)
    demand_input_z2 = CSV.read(string(inputs_path,sep,"Demand_z2.csv"), DataFrame, header=true)
    #demand_input_z3 = CSV.read(string(inputs_path,sep,"Demand_z3.csv"), DataFrame, header=true)
else
    demand_input_z1 = CSV.read(string(inputs_path,sep,"Demand_z1.csv"), DataFrame, header=true, select=["Time_index", "Load"])
    demand_input_z2 = CSV.read(string(inputs_path,sep,"Demand_z2.csv"), DataFrame, header=true, select=["Time_index", "Load"])
    #demand_input_z3 = CSV.read(string(inputs_path,sep,"Demand_z3.csv"), DataFrame, header=true, select=["Time_index", "Load"])
end

resources_input = CSV.read(string(inputs_path,sep,"Resources_updated_new_version.csv"), DataFrame, header=true)
resource_avail_input_z1 = CSV.read(string(inputs_path,sep,"Resources_availability_z1_new_version.csv"), DataFrame, header=true)
resource_avail_input_z2 = CSV.read(string(inputs_path,sep,"Resources_availability_z2_new_version.csv"), DataFrame, header=true)
#resource_avail_input_z3 = CSV.read(string(inputs_path,sep,"Resources_availability_z3.csv"), DataFrame, header=true)

time_index = demand_input_z1[:,1]

#Transmission and zone input values
transmission_input = CSV.read(string(inputs_path,sep,"trans_cap_and_zones_new_version.csv"), DataFrame, header=true)

# ~~~
# Model
# ~~~

## Sets
T = size(demand_input_z1)[1] # number of time steps
if stochastic
    S = size(demand_input_z1[:,2:end])[2] # number of scenarios
else
    S = 1
end
R_all = size(resources_input)[1]
R =  length(resources_input[(resources_input[!,:Generation].==1),:][!,:Index_ID]) #size(resources_input)[1]# number of supply technologies
#resources = resources_input[:,1:5]

#Number of storage technologies
R_S = length(resources_input[(resources_input[!,:Storage].==1),:][!,:Index_ID])
#transmission and zones
L = size(transmission_input[:,3])[1]
Z = size(transmission_input[:,2])[1]
#subset_Z =[1,3] # to exclude nse from zone 2 which has no demand
## Parameters 
cost_inv = resources_input[1:R_all,"Investment cost"] # €/MW-336days
cost_var = resources_input[1:R_all, "Operating cost"] # €/MWh
CO2_Tax = resources_input[1:R_all, "CO2_tax_per_ton_CO2"]


#cost_inv_transmission = zeros(R_all,Z) #Cost variables for transmission investemnet in the north sea grid Z2
#cost_inv_transmission[:,1] = resources_input[1:R_all,"Investment cost transmission z1"]
#cost_inv_transmission[:,2] = resources_input[1:R_all,"Investment cost transmission z2"] 
#cost_inv_transmission[:,3] = resources_input[1:R_all,"Investment cost transmission z3"]

availability = zeros(T,R_all,Z)
availability[:,:,1] = Matrix(resource_avail_input_z1[:, 2:end])
availability[:,:,2] = Matrix(resource_avail_input_z2[:, 2:end])
#availability[:,:,3] = Matrix(resource_avail_input_z3[:, 2:9]) 

co2_factors = resources_input[:,"Emissions_ton_per_MWh"] # ton per MWh
price_cap = 3000 # €/MWh
carbon_cap = 72775794.4  # tCO2, -90% from 2005 emission levels

# Storage
# Current parameters assume battery storage
p_e_ratio = 1/4 # Power to energy ratio
# Single-trip efficiecy
eff_down = 0.9
eff_up = 0.9

if risk_aversion
    # CVaR
    # Define parameters
    α = 1/4 # parameter for VaR set to 1/4 to look at worst case scenario which is the high demand scenario
    γ = 0.5 # parameter for degree of risk aversion, 1 means max/perfect risk aversion
else
    γ = 0
end

demand = zeros(T,S,Z)
demand[:,:,1] = Array(demand_input_z1[:,2:end])
demand[:,:,2] = Array(demand_input_z2[:,2:end])
#demand[:,:,3] = Array(demand_input_z3[:,2:end])

# Define scenario probabilities
P = ones(S)*1/S # uniform distribution

#Transmission parameters
MaxTransCapacity = transmission_input[:,7]
MinTransCapacity = transmission_input[:,8]
LineLoss = transmission_input[:,9]
ZoneMap = transmission_input[:,5:6]


gep = Model(Gurobi.Optimizer)

# Supply
# Power supply, aka generation but could also be from storage in the future
@variable(gep, g[r in 1:R, t in 1:T, s in 1:S, z in 1:Z] >= 0) # amount of power supply, MW

# Capacity
@variable(gep, x[r in 1:R_all, z in 1:Z] >= 0) # Capacity, MW

# Non-served energy
@variable(gep, nse[t in 1:T, s in 1:S, z in 1:Z] >= 0)

# Storage
@variable(gep, discharge[r in 1:R_all, t in 1:T, s in 1:S, z in 1:Z] >= 0)
@variable(gep, charge[r in 1:R_all, t in 1:T, s in 1:S, z in 1:Z] >= 0)
# State of charge
@variable(gep, e[r in R:R_all, t in 1:T, s in 1:S, z in 1:Z] >= 0)

##Variables, Zones and Transmission
@variable(gep, MaxTransCap[l in 1:L] >= 0)
#PF for every Transmission line for every time step t
@variable(gep,Flow[t in 1:T, l in 1:L])

# Risk aversion
if risk_aversion
    # Auxiliary varliables for CVaR
    @variable(gep, u[s in 1:S] >= 0) # loss relative to VaR, €/MW
    @variable(gep, VaR) # VaR variable, €/MW
    @constraint(gep, u_expression[s in 1:S], u[s] >= sum(g[r,t,s,z]*cost_var[r] for r in 1:R, t in 1:T, z in 1:Z) + sum(price_cap*nse[t,s,z] for r in 1:R, t in 1:T, z in 1:Z) - VaR)
end 

# Objective function
if CO2_tax_flag
    if risk_aversion
        @objective(gep, Min, sum(x[r,z]*cost_inv[r] for r in 1:R_all, z in 1:Z) +(1-γ)*(sum(P[s]*g[r,t,s,z]*cost_var[r] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z) + sum(P[s]*price_cap*nse[t,s,z] for t in 1:T, s in 1:S, z in 1:Z)+ sum(P[s]*CO2_Tax[r]*co2_factors[r]*g[r,t,s,z] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z)) + γ*(VaR + 1/α*sum(P[s]*u[s] for s in 1:S)))
    else
        # Risk neutral
        @objective(gep, Min, sum(x[r,z]*cost_inv[r] for r in 1:R_all, z in 1:Z) + sum(P[s]*g[r,t,s,z]*cost_var[r] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z) + sum(P[s]*price_cap*nse[t,s,z] for t in 1:T, s in 1:S, z in 1:Z)+ sum(P[s]*CO2_Tax[r]*co2_factors[r]*g[r,t,s,z] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z))
    end
else
    if risk_aversion
        @objective(gep, Min, sum(x[r,z]*cost_inv[r] for r in 1:R_all, z in 1:Z) +(1-γ)*(sum(P[s]*g[r,t,s,z]*cost_var[r] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z) + sum(P[s]*price_cap*nse[t,s,z] for t in 1:T, s in 1:S, z in 1:Z)) + γ*(VaR + 1/α*sum(P[s]*u[s] for s in 1:S)))
    else
        # Risk neutral
        @objective(gep, Min, sum(x[r,z]*cost_inv[r] for r in 1:R_all, z in 1:Z) + sum(P[s]*g[r,t,s,z]*cost_var[r] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z) + sum(P[s]*price_cap*nse[t,s,z] for t in 1:T, s in 1:S, z in 1:Z))
    end
end


@constraint(gep, PowerBalance[t in 1:T, s in 1:S, z in 1:Z], sum(g[r,t,s,z] for r in 1:R) + nse[t,s,z]+ sum(discharge[r,t,s,z] - charge[r,t,s,z] for r in (resources_input[(resources_input[!,:Storage].==1),:][!,:Index_ID]))  + sum(Flow[t,l]*ZoneMap[l,z] for l in 1:L) == demand[t,s,z])

#@constraint(gep, PowerBalance_z3[t in 1:T,s in 1:S,z in Z], sum(g[r,t,s,z] for r in 1:R) + sum(Flow[t,l]*ZoneMap[l,z] for l in 1:L) == demand[t,s,z])

@constraint(gep, CapacityLimit[r in 1:R, t in 1:T, s in 1:S, z in 1:Z], g[r,t,s,z] <= x[r,z]*availability[t,r,z])


#CAPACITY CONSTRAINTS ZONE 1
@constraint(gep, x[1,1] <= 8500)
@constraint(gep, x[4,1] <= 108045.33)
@constraint(gep, x[5,1] <= 43906.825)
@constraint(gep, x[6,1] == 0)
@constraint(gep, x[7,1] <= 423697.5)

#CAPACITY CONSTRAINTS ZONE 2 NOW REPRESENTIN UK
@constraint(gep, x[1,2] <= 6000)
@constraint(gep, x[4,2] <= 80170.994)
@constraint(gep, x[5,2] <= 88526)
@constraint(gep, x[7,2] <= 371750)
@constraint(gep, x[8,2] == 0)
@constraint(gep, x[9,2] == 0)

#CAPACITY CONSTRAINTS ZONE 3 NOW REPRESENTING THE OFFSHORE GRID
#@constraint(gep, x[1,3] == 0)
#@constraint(gep, x[2,3] == 0)
#@constraint(gep, x[3,3] == 0)
#@constraint(gep, x[4,3] == 0)
@constraint(gep, x[8,1] <= 85779.12)
#@constraint(gep, x[7,3] == 0)
#@constraint(gep, x[8,3] == 0)



#Storage constraints
# Loop through storage technologies
for r in (resources_input[(resources_input[!,:Storage].==1),:][!,:Index_ID])
    # Wrap first and last periods
    @constraint(gep, state_of_charge_start[t in 1:1, s in 1:S, z in 1:Z], e[r,t,s,z] == e[r,T,s,z] - (1/eff_down)*discharge[r,T,s,z] + eff_up*charge[r,T,s,z])
    # Energy balance for the remaining periods
    @constraint(gep, state_of_charge[t in 2:T, s in 1:S,z in 1:Z], e[r,t,s,z] == e[r,t-1,s,z] - (1/eff_down)*discharge[r,t-1,s,z] + eff_up*charge[r,t-1,s,z]) 
    @constraint(gep, energy_limit[t in 1:T, s in 1:S, z in 1:Z], e[r,t,s,z] <= (1/p_e_ratio)*x[r,z])
    @constraint(gep, charge_limit_total[t in 1:T, s in 1:S, z in 1:Z], charge[r,t,s,z] <= (1/eff_up)*x[r,z])
    @constraint(gep, charge_limit[t in 1:T, s in 1:S, z in 1:Z], charge[r,t,s,z] <= (1/p_e_ratio)*x[r,z] - e[r,t,s,z])
    @constraint(gep, discharge_limit_total[t in 1:T, s in 1:S, z in 1:Z], discharge[r,t,s,z] <= eff_down*x[r,z])
    @constraint(gep, discharge_limit[t in 1:T, s in 1:S, z in 1:Z], discharge[r,t,s,z] <= e[r,t,s,z])
    @constraint(gep, charge_discharge_balance[t in 1:T, s in 1:S, z in 1:Z], (1/eff_down)*discharge[r,t,s,z] + eff_up*charge[r,t,s,z] <= x[r,z])
end

# Climate policy
if co2_cap_flag
    @constraint(gep, emissions_cap[s in 1:S], sum(g[r,t,s,z]*co2_factors[r] for r in 1:R, t in 1:T, z in 1:Z) <= carbon_cap)
end

## Constraints, Transmission and zones
@constraint(gep, MaxFlowPos[l in 1:L,t in 1:T], Flow[t,l] <= MaxTransCapacity[l])
@constraint(gep, MaxFlowNeg[l in 1:L,t in 1:T], Flow[t,l] >= MinTransCapacity[l])


JuMP.optimize!(gep)

# Report Results
obj = objective_value(gep)
gen = value.(g)
nse_all_new = value.(nse)
cap = value.(x)
if co2_cap_flag
    co2_price = dual.(emissions_cap)
end 


# Collect results (the code needs work, currently only compiles results for one scenario)
scenario_for_results = 1
df_cap = DataFrames.DataFrameStyle()
df_cap = DataFrame(Resource = resources_input[:,1], Capacity_z1_CN = cap[1:R+1,1], Capacity_z2_GB =cap[1:R+1,2])#, Capacity_z3_NS = cap[1:R+1,3])
df_gen_z1_CN = DataFrame(transpose(gen[:,:,scenario_for_results,1]), resources_input[1:R,1])
insertcols!(df_gen_z1_CN, 1, :Time => time_index)
df_gen_z2_SK = DataFrame(transpose(gen[:,:,scenario_for_results,2]), resources_input[1:R,1])
insertcols!(df_gen_z2_SK, 1, :Time => time_index)




#Flow on transmissions lines
df_Flow_new = DataFrame(CN_to_GB = value.(Flow[:,1,1]),CN_to_NS = value.(Flow[:,2,1]))#, GB_to_NS = value.(Flow[:,3,1]))
#insertcols!(df_Flow, 1, :Time => time_index)



df_emission_sum =zeros(R+1,Z)
df_gen_sum = zeros(R,Z)
for r in 1:R
    for z in 1:Z
        df_gen_sum[r,z] = sum(gen[r,:,scenario_for_results,z])
        df_emission_sum[r,z] = sum(gen[r,:,scenario_for_results,z]*co2_factors[r])
    end
end
Scenario_stack = ["Base Demand","Low Demand","High Demand", "Very High Demand"]
Max_price_all_scen =  zeros(S,Z)
Mean_price_all_scen = zeros(S,Z)
nse_sumScen =zeros(S,2)
for s in 1:S
    for z in 1:Z
        if z < Z
            Max_price_all_scen[s,z]=maximum(dual.(PowerBalance)[:,s,z])
            Mean_price_all_scen[s,z]=mean(dual.(PowerBalance)[:,s,z])
        else
            Max_price_all_scen[s,z]=maximum(dual.(PowerBalance)[:,s,1])
            Mean_price_all_scen[s,z]=mean(dual.(PowerBalance)[:,s,1])
        end
    end
end
for s in 1:S
    for z in 1:2
        nse_sumScen[s,z] = sum(nse_all_new[:,s,z])
    end
end


df_nse_sum_all_scen = DataFrame(Scenario=Scenario_stack[:], Z1_NSE = nse_sumScen[:,1], Z2_NSE = nse_sumScen[:,2])

emission_sum_per_scen = zeros(S,3)
for s in 1:S
    emission_sum_per_scen[s,1] = sum(gen[2,:,s,:] * co2_factors[2])
    emission_sum_per_scen[s,2] = sum(gen[3,:,s,:] * co2_factors[3])
    emission_sum_per_scen[s,3] = emission_sum_per_scen[s,1]+emission_sum_per_scen[s,2]
end
df_sum_em_per_scen = DataFrame(Scenarios = Scenario_stack[:], Coal = emission_sum_per_scen[:,1], Gas = emission_sum_per_scen[:,2], Total = emission_sum_per_scen[:,3])

df_Total_system_cost = DataFrame(Total_System_Cost = obj)
