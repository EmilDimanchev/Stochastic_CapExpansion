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
using Plots
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
else
    S = 1
end

R = size(resources_input)[1] # number of supply technologies
resources = resources_input[:,1]

# Parameters
cost_inv = resources_input[:,"Investment cost"] # $/MW-336days
cost_var = resources_input[:, "Operating cost"] # $/MWh
availability = Matrix(resource_avail_input[:, 2:end])
co2_factors = resources_input[:,"Emissions_ton_per_MWh"] # ton per MWh
price_cap = 70 # $/MWh
carbon_cap = 300e3 # tCO2, arbitrary

if risk_aversion
    # CVaR
    # Define parameters
    α = 1/3 # parameter for VaR
    γ = 0.5 # parameter for degree of risk aversion, 1 means no risk aversion
end

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

@constraint(gep, PowerBalance[t in 1:T, s in 1:S], sum(g[r,t,s] for r in 1:R) + nse[t,s] == demand[t,s])
@constraint(gep, CapacityLimit[r in 1:R, t in 1:T, s in 1:S], g[r,t,s] <= x[r]*availability[t,r])

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
co2_expected = sum(P[s]*co2_tot[s] for s in 1:S)


CSV.write(string(results_path,sep,"capacity.csv"), df_cap)
CSV.write(string(results_path,sep,"generation.csv"), df_gen)
CSV.write(string(results_path,sep,"price.csv"), df_price)
CSV.write(string(results_path,sep,"nse.csv"), df_nse)


#ploting the recomended output
display(Plots.plot(df_cap.Resource, df_cap.Capacity,title ="Production capacity",ylabel="Capacity [MWh]" ,seriestype =[:bar], palette = cgrad(:greens), fill=0, alpha=0.6))

#Estimate revenues per tecnology
revenues = zeros(3)     #[gas, wind, solar]

for r in 1:R-1
    revenues[r] = sum((df_price[i,2]-var_cost[r])*df_gen[i,r+1] for i in 1:T)
end

#Estimation of revenue per technology

df_rev = DataFrame(Resource = resources, Revenue = revenues)
CSV.write(string(results_path,sep,"revenue.csv"), df_rev)
display(Plots.plot(df_rev.Resource, df_rev.Revenue, title = "Revenue per resource", ylabel = "Revenue [USD]" ,seriestype =[:bar], palette = cgrad(:blues), fill=0, alpha=0.6))


#Estimation of emissions

CO2_price = 70 #[$/tCO2]
CO2_emission_factor = [0.4,0,0] #[tCO2/MWh] [gas, wind, solar]
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

Sum_Emissions = sum(df_emission.Gas + df_emission.Wind +df_emission.Solar)
print("Total Emissions [ton CO2]: ")
display(Sum_Emissions)

#None served Energy

display(Plots.plot(df_nse.Time, df_nse.Non_served_energy, xlabel = "Time [h]", ylabel="[MWh]", title ="Non Served Energy", label = "NSE"))


SumNSE = sum(df_nse.Non_served_energy) 
print("Non served energy [MWh]: ")
display(SumNSE)


#Revenue on batteries
c_trans = transpose(c)
d_trans =transpose(d)

df_storage = DataFrame(Charge = vec(c_trans), Discharge =vec(d_trans))#,Discharge=transpose(d), StateOfCharge = transpose(e))
insertcols!(df_storage,1,:Time => time_index)
#display(df_storage)

revenue_storage = zeros(T)
for i in 1:T
    revenue_storage[i] = df_storage.Charge[i]*df_price.Price[i] + df_storage.Discharge[i]*df_price.Price[i]
end
#display(cost_inv)
print("Investment cost storage: ",cost_inv[5])
tot_revenue_storage = sum(revenue_storage) - cost_inv[5]

print("Storage Revenue: ", tot_revenue_storage)
