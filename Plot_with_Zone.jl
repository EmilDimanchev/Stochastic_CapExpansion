include("gep_stochastic.jl")
using Plots; theme(:bright)
using RDatasets
using StatsPlots
using Distributions
using StatsBase
using DataFrames
using JuMP
using CSV
using Gurobi
using Random
using Statistics



if Sys.isunix()
    sep = "/"
elseif Sys.iswindows()
    sep = "\U005c"
end

working_path = pwd()

# Define input and output paths
inputs_path_plots = string(working_path, sep, "Results", sep, "Results_stochastic_031422") #"Inputs_course_3techs", "Inputs_course_all_techs_annual_2041"
results_path_plots = string(working_path, sep, "Results", sep, "Plots_with_zones")

if !(isdir(results_path_plots))
    mkdir(results_path_plots)
end

Palette = cgrad(:Paired_12)
Palette_dark = cgrad(:RdGy)
Palette_bright = cgrad(:seaborn_pastel)


#Plotting capacities

Zone_stack =["","Z1","Z2","Z3"]
Zones = ["Z1 CN" , "Z2 SK" , "Z3 GB"]


Capacities = CSV.read(string(inputs_path_plots,sep,"capacity_RA_co2_cap.csv"), DataFrame, header=true)
plt_cap_all = zeros(R_all,Z)
for i in 1:Z
    plt_cap_all[:,i] = Capacities[:,i+1]
end
plt_cap = transpose(plt_cap_all)
Grouped_barPlot_Cap = StatsPlots.groupedbar([3;plt_cap[1:Z,:]],bar_position =:stack, bar_width=0.9, xlim = (1.47,6), xticks =(1:8,Zone_stack),title ="Capacities for each Zone",label = ["Nuclear" "Coal" "Gas" "Onshore Wind" "Offshore Wind Fixed" "Offshore Wind Float" "Solar" "Storage"], ylabel = "MW",xlabel ="Zones", color =[Palette[3] Palette_dark[11] Palette_dark[8] Palette[1] Palette[2] Palette[9] Palette[11] Palette[6]])
Plots.savefig(Grouped_barPlot_Cap,"Results/Plots_with_zones/grouped_Cap_with_zones_RA_CO2_Cap.pdf")


#Plotting Flows
Flows_plt = CSV.read(string(inputs_path_plots,sep,"Transmission_Flows_RA_co2_tax.csv"), DataFrame, header=true)
hours = value.(Flows_plt[:,1])

plt_flow_CN_to_SK = Plots.plot(hours[1:300],Flows_plt[1:300,2],tickfontsize =10,titel = "Flow on transmission line CN-SK",label = "Flow CN-SK",xlabel = "Time [h]",ylabel = "MW")
plt_flow_CN_to_GB = Plots.plot(hours[1:300],Flows_plt[1:300,3],tickfontsize =10,titel = "Flow on transmission line CN-GB",label = "Flow CN-GB",xlabel = "Time [h]",ylabel = "MW")
plt_flow_SK_to_GB = Plots.plot(hours[1:300],Flows_plt[1:300,4],tickfontsize =10,titel = "Flow on transmission line SK-GB",label = "Flow SK-GB",xlabel = "Time [h]",ylabel = "MW")
Plots.savefig(plt_flow_CN_to_SK,"Results/Plots_with_zones/Transmission_flow_CN_to_SK.pdf")
Plots.savefig(plt_flow_CN_to_GB,"Results/Plots_with_zones/Transmission_flow_CN_to_GB.pdf")
Plots.savefig(plt_flow_SK_to_GB,"Results/Plots_with_zones/Transmission_flow_SK_to_GB.pdf")
Combined_Plots_flow=Plots.plot(plt_flow_CN_to_SK,plt_flow_CN_to_GB,plt_flow_SK_to_GB, layout =(3,1),plot_title ="Flow on all Transmission lines")
Plots.savefig(Combined_Plots_flow,"Results/Plots_with_zones/Transmission_flow_combined_RA_co2_cap.pdf")


#Plotting emissions
Emission_per_zone = CSV.read(string(inputs_path_plots,sep,"emissions_per_zone_RA_co2_cap.csv"), DataFrame, header=true)
em_zone = zeros(Z)
for i in 1:Z
 em_zone[i] = Emission_per_zone[1,i]
end
plt_Emission_zone = Plots.plot(Zones,em_zone,label = "Ton CO2",title = "Emissions per zone",xlabel="Zone",ylabel="Ton CO2",seriestype =[:bar])
Plots.savefig(plt_Emission_zone,"Results/Plots_with_zones/Emissions_per_zone_RA_co2_cap.pdf")


#Plotting generation
Generation_z1 = CSV.read(string(inputs_path_plots,sep,"generation_z1_CN_RA_co2_cap.csv"), DataFrame, header=true)
Generation_z2 = CSV.read(string(inputs_path_plots,sep,"generation_z2_SK_RA_co2_cap.csv"), DataFrame, header=true)
Generation_z3 = CSV.read(string(inputs_path_plots,sep,"generation_z3_GB_RA_co2_cap.csv"), DataFrame, header=true)
Gen_all = zeros(R,Z)
for i in 1:R
    Gen_all[i,1] = sum(Generation_z1[:,i+1]) 
    Gen_all[i,2] = sum(Generation_z2[:,i+1]) 
    Gen_all[i,3] = sum(Generation_z3[:,i+1])
end 
plt_gen_all = transpose(Gen_all)
Grouped_barplt_gen = StatsPlots.groupedbar([3;plt_gen_all[1:Z,:]],bar_position =:stack, bar_width=0.9, xlim = (1.47,6), xticks =(1:8,Zone_stack),title ="Generation for each Zone",label = ["Nuclear" "Coal" "Gas" "Onshore Wind" "Offshore Wind Fixed" "Offshore Wind Float" "Solar"], ylabel = "MWh",xlabel ="Zones", color =[Palette[3] Palette_dark[11] Palette_dark[8] Palette[1] Palette[2] Palette[9] Palette[11] Palette[6]])
Plots.savefig(Grouped_barplt_gen,"Results/Plots_with_zones/Generation_per_zone_RA_co2_cap.pdf")

#Plotting Revenues
Revenues_all_zones = CSV.read(string(inputs_path_plots,sep,"revenue_RA_co2_cap.csv"), DataFrame, header=true)
plt_rev_all = zeros(R_all,Z)
for i in 1:Z
    plt_rev_all[:,i] = Revenues_all_zones[:,i+1]/1000000
end

plt_rev = transpose(plt_rev_all)
Grouped_barplt_rev = StatsPlots.groupedbar([3;plt_rev[1:Z,:]],bar_position =:stack, bar_width=0.9, xlim = (1.47,6), xticks =(1:8,Zone_stack),title ="Revenues for each Zone",label = ["Nuclear" "Coal" "Gas" "Onshore Wind" "Offshore Wind Fixed" "Offshore Wind Float" "Solar" "Storage"], ylabel = " Million EUR",xlabel ="Zones", color =[Palette[3] Palette_dark[11] Palette_dark[8] Palette[1] Palette[2] Palette[9] Palette[11] Palette[6]])
Plots.savefig(Grouped_barplt_rev,"Results/Plots_with_zones/Revenue_per_zone_RA_co2_cap.pdf")


#Plotting Power Price
Price_all_zones = CSV.read(string(inputs_path_plots,sep,"price_RA_co2_tax.csv"), DataFrame, header=true)
plt_price_z1 = Plots.plot(hours[1:300],Price_all_zones[1:300,2],tickfontsize =10,titel = "Flow on transmission line CN-SK",label = "Power Price Z1",xlabel = "Time [h]",ylabel = "EUR/MWh")
plt_price_z2 = Plots.plot(hours[1:300],Price_all_zones[1:300,3],tickfontsize =10,titel = "Flow on transmission line CN-SK",label = "Power Price Z2",xlabel = "Time [h]",ylabel = "EUR/MWh")
plt_price_z3 = Plots.plot(hours[1:300],Price_all_zones[1:300,4],tickfontsize =10,titel = "Flow on transmission line CN-SK",label = "Power Price Z3",xlabel = "Time [h]",ylabel = "EUR/MWh")
plt_comb_price = Combined_Plots_flow=Plots.plot(plt_price_z1,plt_price_z2,plt_price_z3, layout =(3,1),plot_title ="Power Price for each Zone")
Plots.savefig(plt_comb_price,"Results/Plots_with_zones/Power_Price_per_zone_RA_co2_cap.pdf")

#Plot chanrging and discharging of batteries
charge_for_plt = CSV.read(string(inputs_path_plots,sep,"charge_RA_co2_cap.csv"), DataFrame, header=true)
discharge_for_plt = CSV.read(string(inputs_path_plots,sep,"discharge_RA_co2_cap.csv"), DataFrame, header=true)
charge_all = zeros(T,Z)
discharge_all =zeros(T,Z)
plt_charge = zeros(Z)
for i in 1:Z
    charge_all[:,i] = charge_for_plt[:,i]
    discharge_all[:,i] = discharge_for_plt[:,i]
end
plt_charge_z1 = Plots.plot(hours[400:500],charge_all[400:500,1],tickfontsize =10,titel = "Charging of Battery",label = "Charge",xlabel = "Time [h]",ylabel = "MW")
plt_charge_z2 = Plots.plot(hours[1:5000],charge_all[1:5000,2],tickfontsize =10,titel = "Charging of Battery",label = "Charge",xlabel = "Time [h]",ylabel = "MW")
plt_charge_z3 = Plots.plot(hours[1:500],charge_all[1:500,3],tickfontsize =10,titel = "Charging of Battery",label = "Charge",xlabel = "Time [h]",ylabel = "MW")

plt_discharge_z1 = Plots.plot(hours[400:500],discharge_all[400:500,1],tickfontsize =10,titel = "Discharging of Battery",label = "Discharge",xlabel = "Time [h]",ylabel = "MW")
plt_discharge_z2 = Plots.plot(hours[1:5000],discharge_all[1:5000,2],tickfontsize =10,titel = "Discharging of Battery",label = "Discharge",xlabel = "Time [h]",ylabel = "MW")
plt_discharge_z3 = Plots.plot(hours[1:500],discharge_all[1:500,3],tickfontsize =10,titel = "Discharging of Battery",label = "Discharge",xlabel = "Time [h]",ylabel = "MW")

plt_ch_dch_comb_z1 = Combined_Plots_flow=Plots.plot(plt_charge_z1,plt_discharge_z1, layout =(2,1),plot_title ="Charging and Discharging of battery in Z1")
plt_ch_dch_comb_z2 = Combined_Plots_flow=Plots.plot(plt_charge_z2,plt_discharge_z2, layout =(2,1),plot_title ="Charging and Discharging of battery in Z2")
plt_ch_dch_comb_z3 = Combined_Plots_flow=Plots.plot(plt_charge_z3,plt_discharge_z3, layout =(2,1),plot_title ="Charging and Discharging of battery in Z3")





