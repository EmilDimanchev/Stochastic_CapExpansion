using DataFrames
using JuMP
using CSV
using Gurobi
using Plots; theme(:dark)
using gep.jl
using RDatasets

#seriestype =[:bar], palette = cgrad(:greens), fill=0, alpha=0.6
#fillcolor=[:red,:green,:blue,:yellow, :black]

#revenue per tecnology
"""
r_vec() = Vector{Float64}(undef,5) #[nuclear, gas, wind, solar, batteries]
i=1

while i <= T
    r_vec(1) = r_vec(1,1) + (df_price[i] - cost_var[1])*df_gen[2,i];
    r_vec(2) = r_vec(1,2) + (df_price[i] - cost_var[2])*df_gen[3,i];
    r_vec(3) = r_vec(1,3) + (df_price[i] - cost_var[3])*df_gen[4,i];
    r_vec(4) = r_vec(1,4) + (df_price[i] - cost_var[4])*df_gen[5,i];
    #possible to add for batteries here aswell
    i += 1
    
end

x = ["Nuclear","Gas","Wind","Solar","Bateries"]
Plots.plot(x,r_array,seriestype =[:bar], palette = cgrad(:grays), fill=0, alpha=0.6)



#Estimate and compile emissions for every hour


@variable(gep, revenues[t in 1:R, r in 1:R] >= 0)
 
@objective(gep, max, sum((df_price[i]-cost_inv[i]-cost_var[i])*df_gen[i] for i in 1:T)) # Compute the revenues for each tecnology


#r_vec = Vector{Float64}(undef,5) #[nuclear, gas, wind, solar, batteries]
#for i in range(length(resources_input[(resources_input[!,:Generation].==1),:][!,:Index_ID]))
    r_vec = r_vec[1] + (df_price[i] - cost_var[1])*df_gen[2,i];
    r_vec = r_vec[2] + (df_price[i] - cost_var[2])*df_gen[3,i];
    r_vec = r_vec[3] + (df_price[i] - cost_var[3])*df_gen[4,i];
    r_vec = r_vec[4] + (df_price[i] - cost_var[4])*df_gen[5,i];
    #possible to add for batteries here aswell  
#end

df_rev = DataFrame(Resource = resources, revenue = revenues)

CSV.write(string(results_path,sep,"revenue.csv"), df_rev)
Plots.plot(df_rev.Resource, df_rev.revenue ,seriestype =[:bar], palette = cgrad(:grays), fill=0, alpha=0.6)


 #Estimate and compile the emissions
 CO2_price = 70 #[$/tCO2]
 CO2_emission_factor = [0,0.4,0,0,0] #[tCO2/MWh] [Nuclear, gas, wind, solar, batteries]
@variable(gep, emissions[])
@objective(gep, Min, )

#for i in 1:range(5)
    emissions[i] = value(gen[i]) * CO2_emission_factor[i]
#end

 df_emission = DataFrame(Resource =resources, Emission = emissions)

 CSV.write(string(results_path,sep,"emissions_per_tec.csv"), df_emission)
 Plots.plot(df_emission.Resource, df_emissons.emissions ,seriestype =[:bar], palette = cgrad(:grays), fill=0, alpha=0.6)


#Non served energy, just for fun

Plots.plot(df_nse.Time, df_nse.Non_served_energy)

SumNSE = sum(df_nse.Non_served_energy) 
display(SumNSE)



#NEW
revenues = [0,0,0,0,0]
#@variable(Vector{Float64}(undef,5), revenues[t in 1:R, r in 1:R] >= 0) #should I define this as a matrix and not the gep model
 

#r_vec = Vector{Float64}(undef,5) #[nuclear, gas, wind, solar, batteries]
for i in 1:T
    revenues[1] = revenues[1] + (df_price[i] - cost_var[1]- cost_inv[1])*df_gen[i,1];
    revenues[2] = revenues[2] + (df_price[i] - cost_var[2] - cost_inv[2])*df_gen[i,2];
    revenues[3] = revenues[3] + (df_price[i] - cost_var[3]-cost_inv[3])*df_gen[i,3];
    revenues[4] = revenues[4] + (df_price[i] - cost_var[4]-cost_inv[4])*df_gen[i,4];
    #possible to add for batteries here aswell  
end

df_rev = DataFrame(Resource = resources, Revenue = revenues)

CSV.write(string(results_path,sep,"revenue.csv"), df_rev)
Plots.plot(df_rev.Resource, df_rev.Revenue ,seriestype =[:bar], palette = cgrad(:grays), fill=0, alpha=0.6)


#Estimate and compile the emissions
CO2_price = 70 #[$/tCO2]
CO2_emission_factor = [0,0.4,0,0,0] #[tCO2/MWh] [Nuclear, gas, wind, solar, batteries]
@variable(Matrix, emissions[]) 

for j in 1:R
    for i in 1:T
    emissions[i,j] = value(gen[i,j]) * CO2_emission_factor[j]
    end
end
 df_emission = DataFrame(Resource =resources, Emission = emissions)

 CSV.write(string(results_path,sep,"emissions_per_tec.csv"), df_emission)
 Plots.plot(df_emission.Resource, df_emissons.emissions ,seriestype =[:bar], palette = cgrad(:grays), fill=0, alpha=0.6)



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

print("Storage Revenue: ", tot_revenue_storage) """

iris = dataset("datasets","iris")

