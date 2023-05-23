
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

#Choose Ground Rent tax
Ground_rent_tax = true
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
inputs_path = string(working_path, sep, "Inputs", sep, "Inputs_course_all_techs_annual_2041") #"Inputs_course_3techs", "Inputs_course_all_techs_annual_2041"
inputs_path_tax_cost = string(working_path, sep, "Results", sep, "Results_stochastic_031422")
results_path = string(working_path, sep, "Results", sep, "Results_stochastic_031422")

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

#Tax costs input values
if Ground_rent_tax
    GRtax_cost_input_z1 = CSV.read(string(inputs_path_tax_cost,sep,"Ground_Rent_Tax_Cost_Z1_RN_co2_cap.csv"), DataFrame, header=true)
    GRtax_cost_input_z2 = CSV.read(string(inputs_path_tax_cost,sep,"Ground_Rent_Tax_Cost_Z2_RN_co2_cap.csv"), DataFrame, header=true)
    GRtax_cost_input_z3 = CSV.read(string(inputs_path_tax_cost,sep,"Ground_Rent_Tax_Cost_Z3_RN_co2_cap.csv"), DataFrame, header=true)
else
    GRtax_cost_input_z1 = zeros(S,R_all)
    GRtax_cost_input_z2 = zeros(S,R_all)
    GRtax_cost_input_z3 = zeros(S,R_all)
end

#Transmission and zone input values
transmission_input = CSV.read(string(inputs_path,sep,"trans_cap_and_zones.csv"), DataFrame, header=true)

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
#subset_Z =[1,3] # to exclude nse from zone 2 which has no demand
## Parameters 
cost_inv = resources_input[1:R_all,"Investment cost"] # €/MW-336days
cost_var = resources_input[1:R_all, "Operating cost"] # €/MWh
CO2_Tax = resources_input[1:R_all, "CO2_tax_per_ton_CO2"]


cost_inv_transmission = zeros(R_all,Z) #Cost variables for transmission investemnet in the north sea grid Z2
cost_inv_transmission[:,1] = resources_input[1:R_all,"Investment cost transmission z1"]
cost_inv_transmission[:,2] = resources_input[1:R_all,"Investment cost transmission z2"] 
cost_inv_transmission[:,3] = resources_input[1:R_all,"Investment cost transmission z3"]

availability = zeros(T,R_all,Z)
availability[:,:,1] = Matrix(resource_avail_input_z1[:, 2:9])
availability[:,:,2] = Matrix(resource_avail_input_z2[:, 2:9])
availability[:,:,3] = Matrix(resource_avail_input_z3[:, 2:9]) 

#define tax costs
GR_tax_cost = zeros(S,R_all,Z)
GR_tax_cost[:,:,1] = Array(GRtax_cost_input_z1[:,1:end])
GR_tax_cost[:,:,2] = Array(GRtax_cost_input_z2[:,1:end])
GR_tax_cost[:,:,3] = Array(GRtax_cost_input_z3[:,1:end])

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
demand[:,:,3] = Array(demand_input_z3[:,2:end])

# Define scenario probabilities
P = ones(S)*1/S # uniform distribution

#Transmission parameters
MaxTransCapacity = transmission_input[:,8]
MinTransCapacity = transmission_input[:,9]
LineLoss = transmission_input[:,10]
ZoneMap = transmission_input[:,5:7]


gep = Model(Gurobi.Optimizer)

# Supply
# Power supply, aka generation but could also be from storage in the future
@variable(gep, g[r in 1:R, t in 1:T, s in 1:S, z in 1:Z] >= 0) # amount of power supply, MW

# Capacity
@variable(gep, x[r in 1:R_all, z in 1:Z] >= 0) # Capacity, MW

# Non-served energy
@variable(gep, nse[t in 1:T, s in 1:S, z in 1:Z-1] >= 0)

# Storage
@variable(gep, discharge[r in 1:R_all, t in 1:T, s in 1:S, z in 1:Z-1] >= 0)
@variable(gep, charge[r in 1:R_all, t in 1:T, s in 1:S, z in 1:Z-1] >= 0)
# State of charge
@variable(gep, e[r in R:R_all, t in 1:T, s in 1:S, z in 1:Z-1] >= 0)

##Variables, Zones and Transmission
@variable(gep, MaxTransCap[l in 1:L] >= 0)
#PF for every Transmission line for every time step t
@variable(gep,Flow[t in 1:T, l in 1:L])

# Risk aversion
if risk_aversion
    # Auxiliary varliables for CVaR
    @variable(gep, u[s in 1:S] >= 0) # loss relative to VaR, €/MW
    @variable(gep, VaR) # VaR variable, €/MW
    @constraint(gep, u_expression[s in 1:S], u[s] >= sum(g[r,t,s,z]*cost_var[r] for r in 1:R, t in 1:T, z in 1:Z) + sum(price_cap*nse[t,s,z] for r in 1:R, t in 1:T, z in 1:Z-1) - VaR)
end 

# Objective function
if CO2_tax_flag
    if risk_aversion
        @objective(gep, Min, sum(x[r,z]*cost_inv[r]+x[r,z]*cost_inv_transmission[r,z] for r in 1:R_all, z in 1:Z) +(1-γ)*(sum(P[s]*g[r,t,s,z]*cost_var[r] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z) + sum(P[s]*price_cap*nse[t,s,z] for t in 1:T, s in 1:S, z in 1:Z-1)+ sum(P[s]*CO2_Tax[r]*co2_factors[r]*g[r,t,s,z] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z)+sum(P[s]*GR_tax_cost[s,r,z]*x[r,z]*availability[t,r,z] for s in 1:S, r in 1:R, z in 1:Z)) + γ*(VaR + 1/α*sum(P[s]*u[s] for s in 1:S)))
    else
        # Risk neutral
        @objective(gep, Min, sum(x[r,z]*cost_inv[r]+x[r,z]*cost_inv_transmission[r,z] for r in 1:R_all, z in 1:Z) + sum(P[s]*g[r,t,s,z]*cost_var[r] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z) + sum(P[s]*price_cap*nse[t,s,z] for t in 1:T, s in 1:S, z in 1:Z-1)+ sum(P[s]*CO2_Tax[r]*co2_factors[r]*g[r,t,s,z] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z)+sum(P[s]*GR_tax_cost[s,r,z]*x[r,z]*availability[t,r,z] for s in 1:S, r in 1:R, z in 1:Z))
    end
else
    if risk_aversion
        @objective(gep, Min, sum(x[r,z]*cost_inv[r]+x[r,z]*cost_inv_transmission[r,z] for r in 1:R_all, z in 1:Z) +(1-γ)*(sum(P[s]*g[r,t,s,z]*cost_var[r] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z) + sum(P[s]*price_cap*nse[t,s,z] for t in 1:T, s in 1:S, z in 1:Z-1)+sum(P[s]*GR_tax_cost[s,r,z]*x[r,z]*availability[t,r,z] for s in 1:S, r in 1:R, z in 1:Z)) + γ*(VaR + 1/α*sum(P[s]*u[s] for s in 1:S)))
    else
        # Risk neutral
        @objective(gep, Min, sum(x[r,z]*cost_inv[r]+x[r,z]*cost_inv_transmission[r,z] for r in 1:R_all, z in 1:Z) + sum(P[s]*g[r,t,s,z]*cost_var[r] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z) + sum(P[s]*price_cap*nse[t,s,z] for t in 1:T, s in 1:S, z in 1:Z-1)+sum(P[s]*GR_tax_cost[s,r,z]*x[r,z]*availability[t,r,z] for t in 1:T, s in 1:S, r in 1:R, z in 1:Z))
    end
end


@constraint(gep, PowerBalance[t in 1:T, s in 1:S, z in 1:Z-1], sum(g[r,t,s,z] for r in 1:R) + nse[t,s,z]+ sum(discharge[r,t,s,z] - charge[r,t,s,z] for r in (resources_input[(resources_input[!,:Storage].==1),:][!,:Index_ID]))  + sum(Flow[t,l]*ZoneMap[l,z] for l in 1:L) == demand[t,s,z])

@constraint(gep, PowerBalance_z3[t in 1:T,s in 1:S,z in Z], sum(g[r,t,s,z] for r in 1:R) + sum(Flow[t,l]*ZoneMap[l,z] for l in 1:L) == demand[t,s,z])

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

#CAPACITY CONSTRAINTS ZONE 3 NOW REPRESENTING THE OFFSHORE GRID
@constraint(gep, x[1,3] == 0)
@constraint(gep, x[2,3] == 0)
@constraint(gep, x[3,3] == 0)
@constraint(gep, x[4,3] == 0)
@constraint(gep, x[5,3] <= 85779.12)
@constraint(gep, x[7,3] == 0)
@constraint(gep, x[8,3] == 0)


#CAPACITY CONSTRAINTS to represent RA revenue
#@constraint(gep, x[1,1] == 0)
#@constraint(gep, x[2,1] == 139360.15)
#@constraint(gep, x[3,1] == 55058.992)
#@constraint(gep, x[4,1] == 0)
#@constraint(gep, x[5,1] == 11634.448)
#@constraint(gep, x[6,1] == 0)
#@constraint(gep, x[7,1] == 114141.51)
#@constraint(gep, x[8,1] == 0)

#CAPACITY CONSTRAINTS ZONE 2 NOW REPRESENTIN UK

#@constraint(gep, x[1,2] == 0)
#@constraint(gep, x[2,2] == 48359.113)
#@constraint(gep, x[3,2] == 55648.062)
#@constraint(gep, x[4,2] == 80170.994)
#@constraint(gep, x[5,2] == 0)
#@constraint(gep, x[6,2] == 0)
#@constraint(gep, x[7,2] == 26117.319)
#@constraint(gep, x[8,2] == 0)

#CAPACITY CONSTRAINTS ZONE 3 NOW REPRESENTING THE OFFSHORE GRID
#@constraint(gep, x[1,3] == 0)
#@constraint(gep, x[2,3] == 0)
#@constraint(gep, x[3,3] == 0)
#@constraint(gep, x[4,3] == 0)
#@constraint(gep, x[5,3] == 0)
#@constraint(gep, x[7,3] == 0)
#@constraint(gep, x[8,3] == 0)


#Storage constraints
# Loop through storage technologies
for r in (resources_input[(resources_input[!,:Storage].==1),:][!,:Index_ID])
    # Wrap first and last periods
    @constraint(gep, state_of_charge_start[t in 1:1, s in 1:S, z in 1:Z-1], e[r,t,s,z] == e[r,T,s,z] - (1/eff_down)*discharge[r,T,s,z] + eff_up*charge[r,T,s,z])
    # Energy balance for the remaining periods
    @constraint(gep, state_of_charge[t in 2:T, s in 1:S,z in 1:Z-1], e[r,t,s,z] == e[r,t-1,s,z] - (1/eff_down)*discharge[r,t-1,s,z] + eff_up*charge[r,t-1,s,z]) 
    @constraint(gep, energy_limit[t in 1:T, s in 1:S, z in 1:Z-1], e[r,t,s,z] <= (1/p_e_ratio)*x[r,z])
    @constraint(gep, charge_limit_total[t in 1:T, s in 1:S, z in 1:Z-1], charge[r,t,s,z] <= (1/eff_up)*x[r,z])
    @constraint(gep, charge_limit[t in 1:T, s in 1:S, z in 1:Z-1], charge[r,t,s,z] <= (1/p_e_ratio)*x[r,z] - e[r,t,s,z])
    @constraint(gep, discharge_limit_total[t in 1:T, s in 1:S, z in 1:Z-1], discharge[r,t,s,z] <= eff_down*x[r,z])
    @constraint(gep, discharge_limit[t in 1:T, s in 1:S, z in 1:Z-1], discharge[r,t,s,z] <= e[r,t,s,z])
    @constraint(gep, charge_discharge_balance[t in 1:T, s in 1:S, z in 1:Z-1], (1/eff_down)*discharge[r,t,s,z] + eff_up*charge[r,t,s,z] <= x[r,z])
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
obj_tax = objective_value(gep)
gen_tax = value.(g)
nse_all_tax = value.(nse)
cap_tax = value.(x)
if co2_cap_flag
    co2_price_tax = dual.(emissions_cap)
end 
df_cap_GR_Tax = DataFrames.DataFrameStyle()
df_cap_GR_Tax = DataFrame(Resource = resources[1:R+1,1], Capacity_z1_CN = cap_tax[1:R+1,1], Capacity_z2_GB =cap_tax[1:R+1,2], Capacity_z3_NS = cap_tax[1:R+1,3])

CSV.write(string(results_path,sep,"GR_Tax_Capacity_RN_co2_cap.csv"), df_cap_GR_Tax)