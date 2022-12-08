# A Generation Expansion Planning model
# Stochastic version
# Contact: Emil Dimanchev

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
# ~~~

# Choose stochastic or deterministic
stochastic = false
risk_aversion = false
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
inputs_path = string(working_path, sep, "Inputs", sep, "Inputs_course_all_techs_annual_2041") #"Inputs_course_3techs", "Inputs_course_all_techs_annual_2041"
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
else
    S = 1
end
R_all = size(resources_input)[1]
R =  length(resources_input[(resources_input[!,:Generation].==1),:][!,:Index_ID]) #size(resources_input)[1]# number of supply technologies
resources = resources_input[:,1:5]

#Number of storage technologies
R_S = length(resources_input[(resources_input[!,:Storage].==1),:][!,:Index_ID])

# Parameters
cost_inv = resources_input[1:R+R_S,"Investment cost"] # $/MW-336days
cost_var = resources_input[1:R+R_S, "Operating cost"] # $/MWh
availability = Matrix(resource_avail_input[:, 2:7])
co2_factors = resources_input[:,"Emissions_ton_per_MWh"] # ton per MWh
price_cap = 1000 # $/MWh
carbon_cap =  114956929.4 # tCO2, arbitrary

# Storage
# Current parameters assume battery storage
p_e_ratio = 1/8 # Power to energy ratio
# Single-trip efficiecy
eff_down = 0.9
eff_up = 0.9

if risk_aversion
    # CVaR
    # Define parameters
    α = 1/3 # parameter for VaR
    γ = 0.5 # parameter for degree of risk aversion, 1 means no risk aversion
else
    γ = 1
end
demand = Array(demand_input[:,2:end])

# Define scenario probabilities
P = ones(S)*1/S # uniform distribution

gep = Model(Gurobi.Optimizer)


# Supply
# Power supply, aka generation but could also be from storage in the future
@variable(gep, g[r in 1:R, t in 1:T, s in 1:S] >= 0) # amount of power supply, MW

# Capacity
@variable(gep, x[r in 1:R_all] >= 0) # Capacity, MW

# Non-served energy
@variable(gep, nse[t in 1:T, s in 1:S] >= 0)

# Storage
@variable(gep, discharge[r_s in 1:R_S, t in 1:T, s in 1:S] >= 0)
@variable(gep, charge[r_s in 1:R_S, t in 1:T, s in 1:S] >= 0)
# State of charge
@variable(gep, e[r_s in 1:R_S, t in 1:T, s in 1:S] >= 0)

# Risk aversion
if risk_aversion
    # Auxiliary varliables for CVaR
    @variable(gep, u[s in 1:S] >= 0) # loss relative to VaR, $/MW
    @variable(gep, VaR) # VaR variable, $/MW
    @constraint(gep, u_expression[s in 1:S], u[s] >= sum(g[r,t,s]*cost_var[r] + price_cap*nse[t,s] for r in 1:R, t in 1:T) - VaR)
end 

# Objective function
if risk_aversion
    @objective(gep, Min, sum(x[r]*cost_inv[r] for r in 1:R) + γ*(sum(P[s]*g[r,t,s]*cost_var[r] for r in 1:R, t in 1:T, s in 1:S) + sum(P[s]*price_cap*nse[t,s] for t in 1:T, s in 1:S)) + (1-γ)*(VaR + 1/α*sum(P[s]*u[s] for s in 1:S)))
else
    # Risk neutral
    @objective(gep, Min, sum(x[r]*cost_inv[r] for r in 1:R) + sum(P[s]*g[r,t,s]*cost_var[r] for r in 1:R, t in 1:T, s in 1:S) + sum(P[s]*price_cap*nse[t,s] for t in 1:T, s in 1:S))
end

@constraint(gep, PowerBalance[t in 1:T, s in 1:S], sum(g[r,t,s] for r in 1:R) + nse[t,s] + sum(discharge[r_s,t,s] for r_s in 1:R_S) -sum(charge[r_s,t,s] for r_s in 1:R_S) == demand[t,s])
@constraint(gep, CapacityLimit[r in 1:R, t in 1:T, s in 1:S], g[r,t,s] <= x[r]*availability[t,r])
@constraint(gep, x[1]<= 4000)

#Storage constraints
# Loop through storage technologies
for r in (resources_input[(resources_input[!,:Storage].==1),:][!,:Index_ID])
    # Wrap first and last periods
    @constraint(gep, state_of_charge_start[r_s in 1:R_S, t in 1:1, s in 1:S], e[r_s,t,s] == e[r_s,T,s] - (1/eff_down)*discharge[r_s,t,s] + eff_up*charge[r_s,t,s])
    # Energy balance for the remaining periods
    @constraint(gep, state_of_charge[r_s in 1:R_S, t in 2:T, s in 1:S], e[r_s,t,s] == e[r_s,t-1,s] - (1/eff_down)*discharge[r_s,t,s] + eff_up*charge[r_s,t,s]) 
    @constraint(gep, energy_limit[r_s in 1:R_S, t in 1:T, s in 1:S], e[r_s,t,s] <= (1/p_e_ratio)*x[r_s])
    @constraint(gep, charge_limit_total[r_s in 1:R_S, t in 1:T, s in 1:S], charge[r_s,t,s] <= (1/eff_up)*x[r_s])
    @constraint(gep, charge_limit[r_s in 1:R_S, t in 1:T, s in 1:S], charge[r_s,t,s] <= (1/p_e_ratio)*x[r_s] - e[r_s,t,s])
    @constraint(gep, discharge_limit_total[r_s in 1:R_S, t in 1:T, s in 1:S], discharge[r_s,t,s] <= eff_down*x[r_s])
    @constraint(gep, discharge_limit[r_s in 1:R_S, t in 1:T, s in 1:S], discharge[r_s,t,s] <= e[r_s,t,s])
    @constraint(gep, charge_discharge_balance[r_s in 1:R_S, t in 1:T, s in 1:S], (1/eff_down)*discharge[r_s,t,s] + eff_up*charge[r_s,t,s] <= x[r_s])
end



# Climate policy

if co2_cap_flag
    @constraint(gep, emissions_cap[s in 1:S], sum(g[r,t,s]*co2_factors[r] for r in 1:R, t in 1:T) <= carbon_cap)
end

JuMP.optimize!(gep)

# Report Results
obj = objective_value(gep)
gen = value.(g)
nse_all = value.(nse)
cap = value.(x)
if co2_cap_flag
    co2_price = dual.(emissions_cap)
end 
#c = value(charge)
#d = value(discharge)

# Collect results (the code needs work, currently only compiles results for one scenario)
scenario_for_results = 1

df_cap = DataFrame(Resource = resources[1:R+1,1], Capacity = cap[1:R+1])
df_gen = DataFrame(transpose(gen[:,:,scenario_for_results]), resources[1:R,1])
insertcols!(df_gen, 1, :Time => time_index)
df_price = DataFrame(Time = time_index, Price = dual.(PowerBalance)[:,scenario_for_results])
df_nse = DataFrame(Time = time_index, Non_served_energy = nse_all[:,scenario_for_results])

co2_tot = zeros(S)
for s in 1:S
    co2_tot[s] = sum(gen[i,t,s]*co2_factors[i] for i in 1:R, t in 1:T)
end
co2_expected = sum(P[s]*co2_tot[s] for s in 1:S)


CSV.write(string(results_path,sep,"capacity.csv"), df_cap)
CSV.write(string(results_path,sep,"generation.csv"), df_gen)
CSV.write(string(results_path,sep,"price.csv"), df_price[:,:])
CSV.write(string(results_path,sep,"nse.csv"), df_nse)

#Estimate Net income per tecnology

revenues = zeros(R,S)
co2_cost_gas =zeros(S)
m = dual.(CapacityLimit)     #[Nuclear, Coal, gas, wind, solar]
for s in 1:S
    for r in 1:R
        if co2_cap_flag
            revenues[r,s] = sum((m[r,t,s]-co2_price[1]*co2_factors[r])*gen[r,t,s]*-1  for t in 1:T) 
        else
            revenues[r,s] = sum(m[r,t,s]*gen[r,t,s]*-1  for t in 1:T) 
        end
        if revenues[r,s] > abs(1)
           revenues[r,s] = revenues[r,s] - cost_inv[r]
        end
    end           
end

#Apply to dataframe
df_rev = DataFrame(Resource = resources[1:R,1], Revenue = revenues[1:R,1])
CSV.write(string(results_path,sep,"revenue.csv"), df_rev)

#Estimation of emissions
emissions = zeros(T,R)
for j in 1:R
    for i in 1:T
        emissions[i,j] = df_gen[i,j+1] * co2_factors[j]
    end
end

df_emission = DataFrame(emissions,resources[1:R,1])
insertcols!(df_emission, 1, :Time => time_index)
CSV.write(string(results_path,sep,"emissions_per_tec.csv"), df_emission)
Sum_Emissions = sum(df_emission.Gas + df_emission.Wind +df_emission.Solar +df_emission.Nuclear + df_emission.Coal)
print("Total Emissions [ton CO2]: ")
display(Sum_Emissions)

#None served Energy
SumNSE = sum(df_nse.Non_served_energy) 
print("Non served energy [MWh]: ")
display(SumNSE)
print("CO2 cap policy level" )
display(Sum_Emissions*0.2)

#Earnings on batteries
Price = dual.(PowerBalance)
revenue_storage = zeros(S,R_S)
if df_cap[6,2] > 0
    for s in 1:S
        for r_s in 1:R_S
            revenue_storage[s,r_s] = sum(value.(discharge[r_s,t,s])*eff_up*Price[t,s]-value.(charge[r_s,t,s])*Price[t,s] for t in 1:T)
        end
    end
end







#c_trans = transpose(c)
#d_trans =transpose(d)

#df_storage = DataFrame(Charge = vec(c_trans), Discharge =vec(d_trans))#,Discharge=transpose(d), StateOfCharge = transpose(e))
#insertcols!(df_storage,1,:Time => time_index)
#display(df_storage)



#
#for i in 1:T
    #revenue_storage[i] = -df_storage.Charge[i]*df_price.Price[i] + df_storage.Discharge[i]*df_price.Price[i]
#end
#display(cost_inv)
#print("Investment cost storage: ",cost_inv[2])
#tot_revenue_storage = sum(revenue_storage) - cost_inv[2]

#print("Storage Revenue: ", tot_revenue_storage)""""





