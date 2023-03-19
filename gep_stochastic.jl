# A Generation Expansion Planning model
# Stochastic version
# Contact: Emil Dimanchev
# Updated and further developed model by Lars Skjelbred Nygaard

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
inputs_path = string(working_path, sep, "Inputs", sep, "Inputs_course_all_techs_annual_2041") #"Inputs_course_3techs", "Inputs_course_all_techs_annual_2041"
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

## Parameters 
cost_inv = resources_input[1:R_all,"Investment cost"] # €/MW-336days
cost_var = resources_input[1:R_all, "Operating cost"] # €/MWh
CO2_Tax = resources_input[1:R_all, "CO2_tax_per_ton_CO2"]

availability = zeros(T,R_all,Z)
availability[:,:,1] = Matrix(resource_avail_input_z1[:, 2:9])
availability[:,:,2] = Matrix(resource_avail_input_z2[:, 2:9])
availability[:,:,3] = Matrix(resource_avail_input_z3[:, 2:9]) 

co2_factors = resources_input[:,"Emissions_ton_per_MWh"] # ton per MWh
price_cap = 1000 # €/MWh
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
        @objective(gep, Min, sum(x[r,z]*cost_inv[r] for r in 1:R_all, z in 1:Z) +(1-γ)*(sum(P[s]*g[r,t,s,z]*cost_var[r] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z) + sum(P[s]*price_cap*nse[t,s,z] for t in 1:T, s in 1:S, z in 1:Z)) + γ*(VaR + 1/α*sum(P[s]*u[s] for s in 1:S)) + sum(CO2_Tax[r]*co2_factors[r]*g[r,t,s,z] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z))
    else
        # Risk neutral
        @objective(gep, Min, sum(x[r,z]*cost_inv[r] for r in 1:R_all, z in 1:Z) + sum(P[s]*g[r,t,s,z]*cost_var[r] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z) + sum(P[s]*price_cap*nse[t,s,z] for t in 1:T, s in 1:S, z in 1:Z)+ sum(CO2_Tax[r]*co2_factors[r]*g[r,t,s,z] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z))
    end
else
    if risk_aversion
        @objective(gep, Min, sum(x[r,z]*cost_inv[r] for r in 1:R_all, z in 1:Z) +(1-γ)*(sum(P[s]*g[r,t,s,z]*cost_var[r] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z) + sum(P[s]*price_cap*nse[t,s,z] for t in 1:T, s in 1:S, z in 1:Z)) + γ*(VaR + 1/α*sum(P[s]*u[s] for s in 1:S)))
    else
        # Risk neutral
        @objective(gep, Min, sum(x[r,z]*cost_inv[r] for r in 1:R_all, z in 1:Z) + sum(P[s]*g[r,t,s,z]*cost_var[r] for r in 1:R, t in 1:T, s in 1:S, z in 1:Z) + sum(P[s]*price_cap*nse[t,s,z] for t in 1:T, s in 1:S, z in 1:Z))
    end
end


@constraint(gep, PowerBalance[t in 1:T, s in 1:S, z in 1:Z], sum(g[r,t,s,z] for r in 1:R) + nse[t,s,z] + sum(discharge[r,t,s,z] - charge[r,t,s,z] for r in (resources_input[(resources_input[!,:Storage].==1),:][!,:Index_ID]))  + sum(Flow[t,l]*ZoneMap[l,z] for l in 1:L) == demand[t,s,z])
@constraint(gep, CapacityLimit[r in 1:R, t in 1:T, s in 1:S, z in 1:Z], g[r,t,s,z] <= x[r,z]*availability[t,r,z])
@constraint(gep, x[1,1]<= 8500)
@constraint(gep, x[1,2] <= 7000)
@constraint(gep, x[1,3] <= 6000)
#@constraint(gep, x[3,2] == 10000)
#@constraint(gep, x[4,2] == 10000)
#@constraint(gep, x[5,2] == 0)
#@constraint(gep, x[6,2] == 10000)
#@constraint(gep, x[7,2] == 10000)
#@constraint(gep, x[2,2] == 0)

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

#Transmission losses



JuMP.optimize!(gep)

# Report Results
obj = objective_value(gep)
gen = value.(g)
nse_all = value.(nse)
cap = value.(x)
if co2_cap_flag
    co2_price = dual.(emissions_cap)
end 

# Collect results (the code needs work, currently only compiles results for one scenario)
scenario_for_results = 1
df_cap = DataFrames.DataFrameStyle()
df_cap = DataFrame(Resource = resources[1:R+1,1], Capacity_z1_CN = cap[1:R+1,1], Capacity_z2_SK =cap[1:R+1,2], Capacity_z3_GB = cap[1:R+1,3])
df_gen_z1_CN = DataFrame(transpose(gen[:,:,scenario_for_results,1]), resources[1:R,1])
insertcols!(df_gen_z1_CN, 1, :Time => time_index)
df_gen_z2_SK = DataFrame(transpose(gen[:,:,scenario_for_results,2]), resources[1:R,1])
insertcols!(df_gen_z2_SK, 1, :Time => time_index)
df_gen_z3_GB = DataFrame(transpose(gen[:,:,scenario_for_results,3]), resources[1:R,1])
insertcols!(df_gen_z3_GB, 1, :Time => time_index)
df_gen_all = DataFrame(z1 = df_gen_z1_CN, z2 = df_gen_z2_SK, z3 = df_gen_z3_GB)
df_price = DataFrame(Time = time_index, Price_z1 = dual.(PowerBalance)[:,scenario_for_results,1], Price_z2 = dual.(PowerBalance)[:,scenario_for_results,2], Price_z3 = dual.(PowerBalance)[:,scenario_for_results,3])
df_nse = DataFrame(Time = time_index, Non_served_energy_z1 = nse_all[:,scenario_for_results,1],Non_served_energy_z2 = nse_all[:,scenario_for_results,2],Non_served_energy_z3 = nse_all[:,scenario_for_results,3])

co2_tot = zeros(S,Z)
for s in 1:S
    for z in 1:Z
        co2_tot[s,z] = sum(gen[i,t,s,z]*co2_factors[i] for i in 1:R, t in 1:T)
    end
end
co2_expected = sum(P[s]*co2_tot[s] for s in 1:S)


#Estimate Revenues per resource
# and Revenues on storage
revenues = zeros(R_all,S,Z)
Price = dual.(PowerBalance)
m = dual.(CapacityLimit)    
for s in 1:S
    for r in 1:R_all
        for z in 1:Z
            if r == R_all-R_S+1
                revenues[r,s,z] = sum(value.(discharge[r,t,s,z])*eff_down*Price[t,s,z]-value.(charge[r,t,s,z])*Price[t,s,z]*eff_up for t in 1:T)
            else
                revenues[r,s,z] = sum(m[r,t,s,z]*gen[r,t,s,z]*-1 for t in 1:T)
            end
        end 
    end           
end


#Check profitt
profitt = zeros(R_all,S,Z)
for s in 1:S
    for r in 1:R_all
        for z in 1:Z
            if revenues[r,s,z] > abs(1)
                profitt[r,s,z] = revenues[r,s,z] - cost_inv[r]*cap[r,z]
            end
        end
    end
end


#Apply to dataframe
df_rev = DataFrame(Resource = resources[1:R_all,1], Revenue_z1_CN = revenues[1:R_all,scenario_for_results,1],Revenue_z2_SK = revenues[1:R_all,scenario_for_results,2],Revenue_z3_GB = revenues[1:R_all,scenario_for_results,3])


#Estimation of emissions
emissions = zeros(T,R,Z)
for t in 1:T
    for r in 1:R
        for z in 1:Z
            emissions[t,r,z] = gen[r,t,scenario_for_results,z] * co2_factors[r]
        end
    end
end


df_emission = DataFrame(Emissions_z1_CN = sum(emissions[:,1:R,1]),Emissions_z2_SK = sum(emissions[:,1:R,2]),Emissions_z3_GB = sum(emissions[:,1:R,3]))
#insertcols!(df_emission, 1, :Time => time_index)
Sum_Emissions = sum(df_emission[1,:])
print("Total Emissions [ton CO2]: ")
display(Sum_Emissions)

#None served Energy
SumNSE = sum(df_nse.Non_served_energy_z1 + df_nse.Non_served_energy_z2 + df_nse.Non_served_energy_z3) 
print("Non served energy [MWh]: ")
display(SumNSE)


#Flow on transmissions lines
df_Flow = DataFrame(CN_to_SK = value.(Flow[:,1]),CN_to_GB = value.(Flow[:,2]), SK_to_GB = value.(Flow[:,3]))
insertcols!(df_Flow, 1, :Time => time_index)


#Charging and discharging of battery
df_charge = DataFrame(z1 = value.(charge[R_all,:,scenario_for_results,1]),z2 = value.(charge[R_all,:,scenario_for_results,2]),z3 = value.(charge[R_all,:,scenario_for_results,3]))
df_discharge = DataFrame(z1 = value.(discharge[R_all,:,scenario_for_results,1]),z2 = value.(discharge[R_all,:,scenario_for_results,2]),z3 = value.(discharge[R_all,:,scenario_for_results,3]))

#df_soc = DataFrame(z1 = value.(e[R_all,:,scenario_for_results,1]),z2 = value.(e[R_all,:,scenario_for_results,2]),z3 = value.(e[R_all,:,scenario_for_results,3]))
#df_stateOfCharge = DataFrame(z1 = value.(e[R_all,:,scenario_for_results,1]))

max_prices = zeros(Z)
average_price =zeros(Z)
max_prices[1] = maximum(dual.(PowerBalance)[:,scenario_for_results,1])
max_prices[2] = maximum(dual.(PowerBalance)[:,scenario_for_results,2])
max_prices[3] = maximum(dual.(PowerBalance)[:,scenario_for_results,3])
println("Max price z1 [€/MWh]: ", max_prices[1], ", Max price z2 [€/MWh]: ", max_prices[2],", Max price z3 [€/MWh]: ", max_prices[3])
average_price[1] = mean(dual.(PowerBalance)[:,scenario_for_results,1])
average_price[2] = mean(dual.(PowerBalance)[:,scenario_for_results,2])
average_price[3] = mean(dual.(PowerBalance)[:,scenario_for_results,3])
println("Mean price z1 [€/MWh]: ", average_price[1], ", Mean price z2 [€/MWh]: ", average_price[2],", Mean price z3 [€/MWh]: ", average_price[3])


#Write to CSV files

CSV.write(string(results_path,sep,"capacity_RN_co2_cap.csv"), df_cap)
CSV.write(string(results_path,sep,"generation_z1_CN_RN_co2_cap.csv"), df_gen_z1_CN)
CSV.write(string(results_path,sep,"generation_z2_SK_RN_co2_cap.csv"), df_gen_z2_SK)
CSV.write(string(results_path,sep,"generation_z3_GB_RN_co2_cap.csv"), df_gen_z3_GB)
CSV.write(string(results_path,sep,"generation_all_RN_co2_cap.csv"), df_gen_all)
CSV.write(string(results_path,sep,"price_RN_co2_cap.csv"), df_price[:,:])
CSV.write(string(results_path,sep,"nse_RN_co2_cap.csv"), df_nse)
CSV.write(string(results_path,sep,"revenue_RN_co2_cap.csv"), df_rev)
CSV.write(string(results_path,sep,"emissions_per_zone_RN_co2_cap.csv"), df_emission)
CSV.write(string(results_path,sep,"Transmission_Flows_RN_co2_cap.csv"), df_Flow)
CSV.write(string(results_path,sep,"charge_RN_co2_cap.csv"), df_charge)
CSV.write(string(results_path,sep,"discharge_RN_co2_cap.csv"), df_discharge)

#Filenames for scenarios and policies
#_D_No_co2_policies, _RN_No_co2_policies, _RA_No_co2_policies
#_D_co2_cap, _RN_co2_cap, _RA_co2_cap
#_D_co2_tax, _RN_co2_tax, _RA_co2_tax
#_D_both_co2_cap_and_tax, _RN_both_co2_cap_and_tax, _RA_both_co2_cap_and_tax

#Filenames differnt levels of risk aversion
#_RA_0_no_co2_policy,_RA_0.25_no_co2_policy,_RA_0.5_no_co2_policy,_RA_0.75_no_co2_policy, _RA_1_no_co2_policy
