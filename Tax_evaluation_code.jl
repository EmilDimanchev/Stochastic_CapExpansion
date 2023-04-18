


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
co2_cap_flag = false
#CO2-tax policy
CO2_tax_flag = true
#Choose tax
Ground_rent_tax = true
High_price_tax =false
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
inputs_path = string(working_path, sep, "Inputs", sep, "Inputs_course_all_techs_annual_2041")
inputs_path_tax_cost = string(working_path, sep, "Results", sep, "Results_stochastic_031422")
results_path = string(working_path, sep, "Results", sep, "Results_With_Tax")

if !(isdir(results_path))
    mkdir(results_path)
end

# ~~~
# Load data
# ~~~

if stochastic
    demand_input_z1 = CSV.read(string(inputs_path,sep,"Demand_z1.csv"), DataFrame, header=true)
    demand_input_z2 = CSV.read(string(inputs_path,sep,"Demand_z2.csv"), DataFrame, header=true)
    demand_input_z3 = CSV.read(string(inputs_path,sep,"Demand_z3.csv"), DataFrame, header=true)
else
    demand_input_z1 = CSV.read(string(inputs_path,sep,"Demand_z1.csv"), DataFrame, header=true, select=["Time_index", "Load"])
    demand_input_z2 = CSV.read(string(inputs_path,sep,"Demand_z2.csv"), DataFrame, header=true, select=["Time_index", "Load"])
    demand_input_z3 = CSV.read(string(inputs_path,sep,"Demand_z3.csv"), DataFrame, header=true, select=["Time_index", "Load"])
end

resources_input = CSV.read(string(inputs_path,sep,"Resources_updated.csv"), DataFrame, header=true)
resource_avail_input_z1 = CSV.read(string(inputs_path,sep,"Resources_availability_z1.csv"), DataFrame, header=true)
resource_avail_input_z2 = CSV.read(string(inputs_path,sep,"Resources_availability_z2.csv"), DataFrame, header=true)
resource_avail_input_z3 = CSV.read(string(inputs_path,sep,"Resources_availability_z3.csv"), DataFrame, header=true)

time_index = demand_input_z1[:,1]

#Transmission and zone input values
transmission_input = CSV.read(string(inputs_path,sep,"trans_cap_and_zones.csv"), DataFrame, header=true)

#Tax costs input values
if Ground_rent_tax
    GRtax_cost_input_z1 = CSV.read(string(inputs_path_tax_cost,sep,"Ground_Rent_Tax_Cost_Z1.csv"), DataFrame, header=true)
    GRtax_cost_input_z2 = CSV.read(string(inputs_path_tax_cost,sep,"Ground_Rent_Tax_Cost_Z2.csv"), DataFrame, header=true)
    GRtax_cost_input_z3 = CSV.read(string(inputs_path_tax_cost,sep,"Ground_Rent_Tax_Cost_Z3.csv"), DataFrame, header=true)
else
    GRtax_cost_input_z1 = zeros(T,R_all)
    GRtax_cost_input_z2 = zeros(T,R_all)
    GRtax_cost_input_z3 = zeros(T,R_all)
end
if High_price_tax
    HPtax_cost_input_z1 = CSV.read(string(inputs_path_tax_cost,sep,"High_price_Tax_Cost_Z1.csv"), DataFrame, header=true)
    HPtax_cost_input_z2 = CSV.read(string(inputs_path_tax_cost,sep,"High_price_Tax_Cost_Z2.csv"), DataFrame, header=true)
    HPtax_cost_input_z3 = CSV.read(string(inputs_path_tax_cost,sep,"High_price_Tax_Cost_Z3.csv"), DataFrame, header=true)
else
    HPtax_cost_input_z1 = zeros(T,R_all)
    HPtax_cost_input_z2 = zeros(T,R_all)
    HPtax_cost_input_z3 = zeros(T,R_all)
end

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
resources = resources_input[:,1:5]

#Number of storage technologies
R_S = length(resources_input[(resources_input[!,:Storage].==1),:][!,:Index_ID])
#transmission and zones
L = size(transmission_input[:,3])[1]
Z = size(transmission_input[:,2])[1]

## Parameters 
cost_inv = resources_input[1:R_all,"Investment cost"] # €/MW-336days
cost_var = resources_input[1:R_all, "Operating cost"] # €/MWh
CO2_Tax = resources_input[1:R_all, "CO2_tax_per_ton_CO2"]

availability = zeros(T,R_all,Z)
availability[:,:,1] = Matrix(resource_avail_input_z1[:, 2:9])
availability[:,:,2] = Matrix(resource_avail_input_z2[:, 2:9])
availability[:,:,3] = Matrix(resource_avail_input_z3[:, 2:9]) 

co2_factors = resources_input[:,"Emissions_ton_per_MWh"] # ton per MWh
price_cap = 3000 # €/MWh
carbon_cap = 1.7642609241406734e9*0.2  # tCO2, -80% from Deterministic CO2 emisiions without CO2 cap

# Storage
# Current parameters assume battery storage
p_e_ratio = 1/4 # Power to energy ratio
# Single-trip efficiecy
eff_down = 0.9
eff_up = 0.9

if risk_aversion
    # CVaR
    # Define parameters
    α = 1/3 # parameter for VaR
    γ = 0.5 # parameter for degree of risk aversion, 1 means max/perfect risk aversion
else
    γ = 0
end

demand = zeros(T,S,Z)
demand[:,:,1] = Array(demand_input_z1[:,2:end])
demand[:,:,2] = Array(demand_input_z2[:,2:end])
demand[:,:,3] = Array(demand_input_z3[:,2:end])

# Define scenario probabilities
P = ones(S)*1/S # uniform distribution

#Transmission parameters
MaxTransCapacity = transmission_input[:,8]
LineLoss = transmission_input[:,9]
ZoneMap = transmission_input[:,5:7]

#define tax costs
GR_tax_cost = zeros(T,R_all,Z)
GR_tax_cost[:,:,1] = Array(GRtax_cost_input_z1[1:T,1:end])
GR_tax_cost[:,:,2] = Array(GRtax_cost_input_z2[:,1:end])
GR_tax_cost[:,:,3] = Array(GRtax_cost_input_z3[:,1:end])

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
@variable(gep, e[r in 7:R_all, t in 1:T, s in 1:S, z in 1:Z] >= 0)

##Variables, Zones and Transmission
@variable(gep, MaxTransCap[l in 1:L] >= 0)
#PF for every Transmission line for every time step t
@variable(gep,Flow[t in 1:T, l in 1:L])

# Risk aversion
if risk_aversion
    # Auxiliary varliables for CVaR
    @variable(gep, u[s in 1:S] >= 0) # loss relative to VaR, €/MW
    @variable(gep, VaR) # VaR variable, €/MW
    @constraint(gep, u_expression[s in 1:S, z in 1:Z], u[s] >= sum(g[r,t,s,z]*cost_var[r] + price_cap*nse[t,s,z] for r in 1:R, t in 1:T) - VaR)
end 

# Objective function
if CO2_tax_flag
    if risk_aversion
        @objective(gep, Min, sum(x[r,z]*cost_inv[r] for r in 1:R_all, z in 1:Z) +(1-γ)*(sum(P[s]*g[r,t,s,z]*cost_var[r] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z) + sum(P[s]*price_cap*nse[t,s,z] for t in 1:T, s in 1:S, z in 1:Z)) + γ*(VaR + 1/α*sum(P[s]*u[s] for s in 1:S)) + sum(CO2_Tax[r]*co2_factors[r]*g[r,t,s,z] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z) + sum(GR_tax_cost[t,r,z] for t in 1:T, r in 1:R_all, z in 1:Z))
    else
        # Risk neutral
        @objective(gep, Min, sum(x[r,z]*cost_inv[r] for r in 1:R_all, z in 1:Z) + sum(P[s]*g[r,t,s,z]*cost_var[r] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z) + sum(P[s]*price_cap*nse[t,s,z] for t in 1:T, s in 1:S, z in 1:Z)+ sum(CO2_Tax[r]*co2_factors[r]*g[r,t,s,z] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z) + sum(GR_tax_cost[t,r,z] for t in 1:T, r in 1:R_all, z in 1:Z))
    end
else
    if risk_aversion
        @objective(gep, Min, sum(x[r,z]*cost_inv[r] for r in 1:R_all, z in 1:Z) +(1-γ)*(sum(P[s]*g[r,t,s,z]*cost_var[r] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z) + sum(P[s]*price_cap*nse[t,s,z] for t in 1:T, s in 1:S, z in 1:Z)) + γ*(VaR + 1/α*sum(P[s]*u[s] for s in 1:S))+ sum(GR_tax_cost[t,r,z] for t in 1:T, r in 1:R_all, z in 1:Z))
    else
        # Risk neutral
        @objective(gep, Min, sum(x[r,z]*cost_inv[r] for r in 1:R_all, z in 1:Z) + sum(P[s]*g[r,t,s,z]*cost_var[r] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z) + sum(P[s]*price_cap*nse[t,s,z] for t in 1:T, s in 1:S, z in 1:Z)+ sum(GR_tax_cost[t,r,z] for t in 1:T, r in 1:R_all, z in 1:Z))
    end
end


@constraint(gep, PowerBalance[t in 1:T, s in 1:S, z in 1:Z], sum(g[r,t,s,z] for r in 1:R) + nse[t,s,z] + sum(discharge[r,t,s,z] - charge[r,t,s,z] for r in (resources_input[(resources_input[!,:Storage].==1),:][!,:Index_ID]))  + sum(Flow[t,l]*ZoneMap[l,z] for l in 1:L) == demand[t,s,z])
@constraint(gep, CapacityLimit[r in 1:R, t in 1:T, s in 1:S, z in 1:Z], g[r,t,s,z] <= x[r,z]*availability[t,r,z])

#Nuclear constraints
@constraint(gep, x[1,1]<= 8500)
@constraint(gep, x[1,2] <= 7000)
@constraint(gep, x[1,3] <= 6000)


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
@constraint(gep, MaxFlowNeg[l in 1:L,t in 1:T], Flow[t,l] >= -MaxTransCapacity[l])

JuMP.optimize!(gep)

# Report Results
obj = objective_value(gep)
gen = value.(g)
nse_all = value.(nse)
cap = value.(x)
if co2_cap_flag
    co2_price = dual.(emissions_cap)
end 