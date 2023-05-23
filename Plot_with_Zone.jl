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
results_path_plots = string(working_path, sep, "Results", sep, "Results_plots_csv")
results_dc_plots = string(working_path, sep, "Results", sep, "Results_DC_plots")
inputs_path = string(working_path, sep, "Inputs", sep, "Inputs_course_all_techs_annual_2041")
if !(isdir(results_path_plots))
    mkdir(results_path_plots)
end

Palette = cgrad(:Paired_12)
Palette_dark = cgrad(:RdGy)
Palette_bright = cgrad(:seaborn_pastel)
Palette_green =cgrad(:YlGn_8)
palette_clear = cgrad(:Set2_6)


#Plotting capacities

Zone_stack =["","Z1","Z2","Z3"]
Zones = ["Z1 CN" , "Z2 UK" , "Z3 NS"]


Capacities = CSV.read(string(inputs_path_plots,sep,"capacity_RA_co2_cap_sensitivity_analysis.csv"), DataFrame, header=true)
plt_cap_all = zeros(8,3)
for i in 1:3
    plt_cap_all[:,i] = Capacities[:,i+1]/1000
end
plt_cap = transpose(plt_cap_all)
Grouped_barPlot_Cap = StatsPlots.groupedbar([3;plt_cap[1:3,:]],bar_position =:stack, bar_width=0.9, xlim = (1.47,6), xticks =(1:8,Zone_stack),title ="Capacities for each Zone",label = ["Nuclear" "Coal" "Gas" "Onshore Wind" "Offshore Wind Fixed" "Offshore Wind Float" "Solar" "Storage"],legendfont =font(8),legendfontsize=6, ylabel = "GW",xlabel ="Zones", color =[Palette[3] Palette_dark[11] Palette_dark[8] Palette[1] Palette[2] Palette[9] Palette[11] Palette[6]])
Plots.savefig(Grouped_barPlot_Cap,"Results/RA_Plots_with_zones/grouped_Cap_RA_CO2_Cap_Sensitivity_analysis.pdf")



#Plotting Flows
Flows_plt = CSV.read(string(inputs_path_plots,sep,"Transmission_Flows_RA_co2_tax_and_Nuclear_const.csv"), DataFrame, header=true)
hours = value.(Flows_plt[:,1])

plt_flow_CN_to_UK = Plots.plot(hours[1:8760],Flows_plt[1:8760,2],tickfontsize =10,titel = "Flow on transmission line CN-SK",label = "Flow Z1-Z2",xlabel = "Time [h]",ylabel = "MW",ylim=(-5000,5000))
plt_flow_CN_to_NS = Plots.plot(hours[1:8760],Flows_plt[1:8760,3],tickfontsize =10,titel = "Flow on transmission line CN-GB",label = "Flow Z1-Z3",xlabel = "Time [h]",ylabel = "MW")#,ylim=(-10000,90000))
#plt_flow_SK_to_GB = Plots.plot(hours[4000:4500],Flows_plt[4000:4500,4],tickfontsize =10,titel = "Flow on transmission line SK-GB",label = "Flow Z2-Z3",xlabel = "Time [h]",ylabel = "MW",ylim=(-3000,3000))
Plots.savefig(plt_flow_CN_to_UK,"Results/RA_Plots_with_zones/Transmission_flow_CN_to_UK.pdf")
Plots.savefig(plt_flow_CN_to_GB,"Results/RA_Plots_with_zones/Transmission_flow_CN_to_NS.pdf")
#Plots.savefig(plt_flow_SK_to_GB,"Results/RA_Plots_with_zones/Transmission_flow_SK_to_GB.pdf")
Combined_Plots_flow=Plots.plot(plt_flow_CN_to_UK,plt_flow_CN_to_NS, layout =(2,1),plot_title ="Flow on all Transmission lines")
Plots.savefig(Combined_Plots_flow,"Results/RA_Plots_with_zones/Transmission_flow_combined_RA_co2_tax_and_Nuclear_const,.pdf")


#Plotting emissions
Emission_per_zone = CSV.read(string(inputs_path_plots,sep,"emissions_per_zone_RA_co2_tax_and_Nuclear_const.csv"), DataFrame, header=true)
em_zone = zeros(1,Z)
for i in 1:Z
    em_zone[i] = Emission_per_zone[1,i]/(1000000)
end
trans_em_zone =transpose(em_zone)
plt_Emission_zone = Plots.plot(Zones,trans_em_zone,label = "MtCO2",title = "Emissions per zone",xlabel="Zone",ylabel="MtCO2",seriestype =[:bar],color =[Palette_dark[8] Palette_dark[8] Palette_dark[8]])
Plots.savefig(plt_Emission_zone,"Results/RA_Plots_with_zones/Emissions_per_zone_RA_co2_tax_and_Nuclear_const.pdf")


#Plotting generation
Generation_z1 = CSV.read(string(inputs_path_plots,sep,"generation_z1_CN_RA_co2_tax_and_Nuclear_const.csv"), DataFrame, header=true)
Generation_z2 = CSV.read(string(inputs_path_plots,sep,"generation_z2_SK_RA_co2_tax_and_Nuclear_const.csv"), DataFrame, header=true)
Generation_z3 = CSV.read(string(inputs_path_plots,sep,"generation_z3_GB_RA_co2_tax_and_Nuclear_const.csv"), DataFrame, header=true)
Gen_all = zeros(R,Z)
for i in 1:R
    Gen_all[i,1] = sum(Generation_z1[:,i+1])/1000000 
    Gen_all[i,2] = sum(Generation_z2[:,i+1])/1000000 
    Gen_all[i,3] = sum(Generation_z3[:,i+1])/1000000
end 
plt_gen_all = transpose(Gen_all)
Grouped_barplt_gen = StatsPlots.groupedbar([3;plt_gen_all[1:Z,:]],bar_position =:stack, bar_width=0.9, xlim = (1.47,6), xticks =(1:8,Zone_stack),title ="Generation for each Zone",label = ["Nuclear" "Coal" "Gas" "Onshore Wind" "Offshore Wind Fixed" "Offshore Wind Float" "Solar"], ylabel = "TWh",xlabel ="Zones", color =[Palette[3] Palette_dark[11] Palette_dark[8] Palette[1] Palette[2] Palette[9] Palette[11] Palette[6]])
Plots.savefig(Grouped_barplt_gen,"Results/RA_Plots_with_zones/Generation_per_zone_D_No_co2_policies.pdf")

#Plotting Revenues
Revenues_all_zones = CSV.read(string(inputs_path_plots,sep,"revenue_RA_1_Capacity_In_D_for_REV_no_policy.csv"), DataFrame, header=true)
plt_rev_all = zeros(R_all,Z)
for i in 1:Z
    plt_rev_all[:,i] = Revenues_all_zones[:,i+1]/1e9
end

plt_rev = transpose(plt_rev_all)
Grouped_barplt_rev = StatsPlots.groupedbar([3;plt_rev[1:Z,:]],bar_position =:stack, bar_width=0.9, ylim=(0,2.5),xlim = (1.47,6), xticks =(1:8,Zone_stack),title ="Revenues for each Zone",label = ["Nuclear" "Coal" "Gas" "Onshore Wind" "Offshore Wind Fixed" "Offshore Wind Float" "Solar" "Storage"], ylabel = " Billion €",xlabel ="Zones", color =[Palette[3] Palette_dark[11] Palette_dark[8] Palette[1] Palette[2] Palette[9] Palette[11] Palette[6]])
Plots.savefig(Grouped_barplt_rev,"Results/RA_Plots_with_zones/Revenue_per_zone_RA_co2_tax_and_Nuclear_const.pdf")


#Plotting Power Price
Price_all_zones = CSV.read(string(inputs_path_plots,sep,"price_RA_co2_tax_and_Nuclear_const.csv"), DataFrame, header=true)
plt_price_z1 = Plots.plot(hours[1:8760],Price_all_zones[1:8760,2],tickfontsize =10,titel = "Power Price Z1",label = "Power Price Z1",xlabel = "Time [h]",ylabel = "EUR/MWh")
plt_price_z2 = Plots.plot(hours[1:8760],Price_all_zones[1:8760,3],tickfontsize =10,titel = "Power Price Z2",label = "Power Price Z2",xlabel = "Time [h]",ylabel = "EUR/MWh")
#plt_price_z3 = Plots.plot(hours[1:8760],Price_all_zones[1:8760,4],tickfontsize =10,titel = "Power price Z3",label = "Power Price Z3",xlabel = "Time [h]",ylabel = "EUR/MWh")
plt_comb_price = Plots.plot(plt_price_z1,plt_price_z2, layout =(2,1),plot_title ="Power Price for Zone 1 and 2")
Plots.savefig(plt_comb_price,"Results/RA_Plots_with_zones/Power_Price_per_zone_RA_co2_tax_and_Nuclear_const.pdf")

#Plot chanrging and discharging of batteries
charge_for_plt = CSV.read(string(inputs_path_plots,sep,"charge_RA_co2_tax_and_Nuclear_const.csv"), DataFrame, header=true)
discharge_for_plt = CSV.read(string(inputs_path_plots,sep,"discharge_RA_co2_tax_and_Nuclear_const.csv"), DataFrame, header=true)
charge_all = zeros(T,Z-1)
discharge_all =zeros(T,Z-1)
plt_charge = zeros(Z-1)
for i in 1:Z-1
    charge_all[:,i] = charge_for_plt[:,i]
    discharge_all[:,i] = discharge_for_plt[:,i]
end
plt_charge_z1 = Plots.plot(hours[400:500],charge_all[400:500,1],tickfontsize =10,titel = "Charging of Battery",label = "Charge",xlabel = "Time [h]",ylabel = "MW")
plt_charge_z2 = Plots.plot(hours[1:5000],charge_all[1:5000,2],tickfontsize =10,titel = "Charging of Battery",label = "Charge",xlabel = "Time [h]",ylabel = "MW")
#plt_charge_z3 = Plots.plot(hours[1:500],charge_all[1:500,3],tickfontsize =10,titel = "Charging of Battery",label = "Charge",xlabel = "Time [h]",ylabel = "MW")

plt_discharge_z1 = Plots.plot(hours[400:500],discharge_all[400:500,1],tickfontsize =10,titel = "Discharging of Battery",label = "Discharge",xlabel = "Time [h]",ylabel = "MW")
plt_discharge_z2 = Plots.plot(hours[1:5000],discharge_all[1:5000,2],tickfontsize =10,titel = "Discharging of Battery",label = "Discharge",xlabel = "Time [h]",ylabel = "MW")
#plt_discharge_z3 = Plots.plot(hours[1:500],discharge_all[1:500,3],tickfontsize =10,titel = "Discharging of Battery",label = "Discharge",xlabel = "Time [h]",ylabel = "MW")

plt_ch_dch_comb_z1 = Combined_Plots_flow=Plots.plot(plt_charge_z1,plt_discharge_z1, layout =(2,1),plot_title ="Charging and Discharging of battery in Z1")
plt_ch_dch_comb_z2 = Combined_Plots_flow=Plots.plot(plt_charge_z2,plt_discharge_z2, layout =(2,1),plot_title ="Charging and Discharging of battery in Z2")
#plt_ch_dch_comb_z3 = Combined_Plots_flow=Plots.plot(plt_charge_z3,plt_discharge_z3, layout =(2,1),plot_title ="Charging and Discharging of battery in Z3")


#Plotting NSE
nse = CSV.read(string(inputs_path_plots,sep,"nse_RA_no policy_VH_demand.csv"), DataFrame, header=true)
plt_nse_all = zeros(T,Z-1)
for i in 1:Z-1
    plt_nse_all[:,i] = nse[:,i]/1000
end

plt_nse = transpose(plt_nse_all)
plt_nse_z1 = Plots.plot(hours[200:300],plt_nse[1,200:300],xlabel = "Time [h]",ylabel = "GWh",label = "NSE Z1")
plt_nse_z2 = Plots.plot(hours[200:300],plt_nse[2,200:300],xlabel = "Time [h]",ylabel = "GWh",label = "NSE Z2")
#plt_nse_z3 = Plots.plot(hours[1:8760],plt_nse[3,:],xlabel = "Time [h]",ylabel = "MWh",label = "NSE Z3")
plt_comb_nse = Plots.plot(plt_nse_z1,plt_nse_z2, layout =(2,1),plot_title ="Non Served Energy for each Zone")
Plots.savefig(plt_comb_nse,"Results/RA_Plots_with_zones/NSE_per_zone_VH_demand_RA_no_policy_h200_200.pdf")


#PLOTS FOR DIFFERENT DEGREES OF RISK AVERSION
RA_stack =["γ=0","γ=0.25","γ=0.5","γ=0.75","γ=1"]
demand_stack = ["Base", "Low","High","Very High"]
# Capacities NO CLIMATE POLICIES
RA_cap_0 =  CSV.read(string(inputs_path_plots,sep,"capacity_RN_No_co2_policies.csv"), DataFrame, header=true)
RA_cap_025= CSV.read(string(inputs_path_plots,sep,"capacity_RA_0.25_no_co2_policy.csv"), DataFrame, header=true)
RA_cap_05 = CSV.read(string(inputs_path_plots,sep,"capacity_RA_No_co2_policies.csv"), DataFrame, header=true)
RA_cap_075 =CSV.read(string(inputs_path_plots,sep,"capacity_RA_0.75_no_co2_policy.csv"), DataFrame, header=true)
RA_cap_1 = CSV.read(string(inputs_path_plots,sep,"capacity_RA_1_no_co2_policy.csv"), DataFrame, header=true)
RA_cap_sum = zeros(8, 5)
for r in 1:8
    RA_cap_sum[r,1] = sum(RA_cap_0[r,2:end])/1000
    RA_cap_sum[r,2] = sum(RA_cap_025[r,2:end])/1000
    RA_cap_sum[r,3] = sum(RA_cap_05[r,2:end])/1000
    RA_cap_sum[r,4] = sum(RA_cap_075[r,2:end])/1000
    RA_cap_sum[r,5] = sum(RA_cap_1[r,2:end])/1000
end
RA_cap_sum_trans = transpose(RA_cap_sum)
df_RA_cap_dif = DataFrame(Resource = resources_input[:,1],Gamma_0=RA_cap_sum[:,1],Gamma_025=RA_cap_sum[:,2],Gamma_05=RA_cap_sum[:,3],Gamma_075=RA_cap_sum[:,4],Gamma_1=RA_cap_sum[:,5])
CSV.write(string(results_path_plots,sep,"Capacity_dif_RA_no_co2_policy.csv"), df_RA_cap_dif)
RA_cap_dif=StatsPlots.groupedbar(RA_cap_sum_trans[1:end,1:end],bar_position =:stack, bar_width=0.9, xlim = (0.2,7), xticks =(1:8,RA_stack),title ="Capacities for Different levels of Risk Aversion",label = ["Nuclear" "Coal" "Gas" "Onshore Wind" "Offshore Wind Fixed" "Offshore Wind Float" "Solar" "Storage"],legendfont =font(8),legendfontsize=6, ylabel = "GW",xlabel ="Level of Risk Aversion", color =[Palette[3] Palette_dark[11] Palette_dark[8] Palette[1] Palette[2] Palette[9] Palette[11] Palette[6]])
Plots.savefig(RA_cap_dif,"Results/New_Differnt_RA_Plots_with_zones/GroupedCapacity_dif_RA_no_co2_policy.pdf")


# Capacities CO2 CAP POLICY
RA_cap_CO2_cap_0 =  CSV.read(string(inputs_path_plots,sep,"capacity_RN_co2_cap.csv"), DataFrame, header=true)
RA_cap_CO2_cap_025= CSV.read(string(inputs_path_plots,sep,"capacity_RA_0.25_co2_cap.csv"), DataFrame, header=true)
RA_cap_CO2_cap_05 = CSV.read(string(inputs_path_plots,sep,"capacity_RA_co2_cap.csv"), DataFrame, header=true)
RA_cap_CO2_cap_075 =CSV.read(string(inputs_path_plots,sep,"capacity_RA_0.75_co2_cap.csv"), DataFrame, header=true)
RA_cap_CO2_cap_1 = CSV.read(string(inputs_path_plots,sep,"capacity_RA_1_co2_cap.csv"), DataFrame, header=true)
RA_cap_CO2_cap_sum = zeros(8, 5)
for r in 1:8
    RA_cap_CO2_cap_sum[r,1] = sum(RA_cap_CO2_cap_0[r,2:end])/1000
    RA_cap_CO2_cap_sum[r,2] = sum(RA_cap_CO2_cap_025[r,2:end])/1000
    RA_cap_CO2_cap_sum[r,3] = sum(RA_cap_CO2_cap_05[r,2:end])/1000
    RA_cap_CO2_cap_sum[r,4] = sum(RA_cap_CO2_cap_075[r,2:end])/1000
    RA_cap_CO2_cap_sum[r,5] = sum(RA_cap_CO2_cap_1[r,2:end])/1000
end
RA_cap_CO2_cap_sum_trans = transpose(RA_cap_CO2_cap_sum)
df_RA_cap_CO2_cap_dif = DataFrame(Resource = resources_input[:,1],Gamma_0=RA_cap_CO2_cap_sum[:,1],Gamma_025=RA_cap_CO2_cap_sum[:,2],Gamma_05=RA_cap_CO2_cap_sum[:,3],Gamma_075=RA_cap_CO2_cap_sum[:,4],Gamma_1=RA_cap_CO2_cap_sum[:,5])
CSV.write(string(results_path_plots,sep,"Capacity_dif_RA_co2_Cap_policy.csv"), df_RA_cap_CO2_cap_dif)
RA_cap_CO2_cap_dif=StatsPlots.groupedbar(RA_cap_CO2_cap_sum_trans[1:end,1:end],bar_position =:stack, bar_width=0.9, xlim = (0.2,7.7), xticks =(1:8,RA_stack),title ="Capacities for Different levels of Risk Aversion",label = ["Nuclear" "Coal" "Gas" "Onshore Wind" "Offshore Wind Fixed" "Offshore Wind Float" "Solar" "Storage"],legendfont =font(8),legendfontsize=5, ylabel = "GW",xlabel ="Level of Risk Aversion", color =[Palette[3] Palette_dark[11] Palette_dark[8] Palette[1] Palette[2] Palette[9] Palette[11] Palette[6]])
Plots.savefig(RA_cap_CO2_cap_dif,"Results/New_Differnt_RA_Plots_with_zones/GroupedCapacity_dif_RA_co2_cap_policy.pdf")


# Capacities CO2 TAX POLICY
RA_cap_tax_0 =  CSV.read(string(inputs_path_plots,sep,"capacity_RN_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)
RA_cap_tax_025= CSV.read(string(inputs_path_plots,sep,"capacity_RA_0.25_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)
RA_cap_tax_05 = CSV.read(string(inputs_path_plots,sep,"capacity_RA_co2_tax_and_Nuclear_const.csv"), DataFrame, header=true)
RA_cap_tax_075 =CSV.read(string(inputs_path_plots,sep,"capacity_RA_0.75_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)
RA_cap_tax_1 = CSV.read(string(inputs_path_plots,sep,"capacity_RA_1_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)
RA_cap_tax_sum = zeros(R_all, 5)
for r in 1:8
    RA_cap_tax_sum[r,1] = sum(RA_cap_tax_0[r,2:end])/1000
    RA_cap_tax_sum[r,2] = sum(RA_cap_tax_025[r,2:end])/1000
    RA_cap_tax_sum[r,3] = sum(RA_cap_tax_05[r,2:end])/1000
    RA_cap_tax_sum[r,4] = sum(RA_cap_tax_075[r,2:end])/1000
    RA_cap_tax_sum[r,5] = sum(RA_cap_tax_1[r,2:end])/1000
end
RA_cap_tax_sum_trans = transpose(RA_cap_tax_sum)
df_RA_cap_tax_dif = DataFrame(Resource = resources_input[:,1],Gamma_0=RA_cap_tax_sum[:,1],Gamma_025=RA_cap_tax_sum[:,2],Gamma_05=RA_cap_tax_sum[:,3],Gamma_075=RA_cap_tax_sum[:,4],Gamma_1=RA_cap_tax_sum[:,5])
CSV.write(string(results_path_plots,sep,"Capacity_dif_RA_co2_tax.csv"), df_RA_cap_tax_dif)
RA_cap_tax_dif=StatsPlots.groupedbar(RA_cap_tax_sum_trans[1:end,1:end],bar_position =:stack, bar_width=0.9, xlim = (0.2,7), xticks =(1:8,RA_stack),title ="Capacities for Different levels of Risk Aversion",label = ["Nuclear" "Coal" "Gas" "Onshore Wind" "Offshore Wind Fixed" "Offshore Wind Float" "Solar" "Storage"],legendfont =font(8),legendfontsize=7, ylabel = "GW",xlabel ="Level of Risk Aversion", color =[Palette[3] Palette_dark[11] Palette_dark[8] Palette[1] Palette[2] Palette[9] Palette[11] Palette[6]])
Plots.savefig(RA_cap_tax_dif,"Results/New_Differnt_RA_Plots_with_zones/GroupedCapacity_dif_RA_co2_Tax_policy.pdf")



# NON SERVED ENERGY FOR DIFFERENT RA

# NSE NO CLIMATE POLICIES
RA_nse_0 =  CSV.read(string(inputs_path_plots,sep,"NSE_sum_all_scen_RN_No_co2_policies.csv"), DataFrame, header=true)
RA_nse_025= CSV.read(string(inputs_path_plots,sep,"NSE_sum_all_scen_RA_0.25_no_co2_policy.csv"), DataFrame, header=true)
RA_nse_05 = CSV.read(string(inputs_path_plots,sep,"NSE_sum_all_scen_RA_No_co2_policies.csv"), DataFrame, header=true)
RA_nse_075 =CSV.read(string(inputs_path_plots,sep,"NSE_sum_all_scen_RA_0.75_no_co2_policy.csv"), DataFrame, header=true)
RA_nse_1 = CSV.read(string(inputs_path_plots,sep,"NSE_sum_all_scen_RA_1_no_co2_policy.csv"), DataFrame, header=true)

RA_NSE_sum = zeros(5,4)
for s in 1:4
    RA_NSE_sum[1,s] = sum(RA_nse_0[s,2:end])/1000
    RA_NSE_sum[2,s] = sum(RA_nse_025[s,2:end])/1000
    RA_NSE_sum[3,s] = sum(RA_nse_05[s,2:end])/1000
    RA_NSE_sum[4,s] = sum(RA_nse_075[s,2:end])/1000
    RA_NSE_sum[5,s] = sum(RA_nse_1[s,2:end])/1000
end

df_RA_NSE_sum = DataFrame(Scenario = RA_nse_0[:,1], Gamma_0 =RA_NSE_sum[1,:], Gamma_025 =RA_NSE_sum[2,:], Gamma_05 =RA_NSE_sum[3,:], Gamma_075 =RA_NSE_sum[4,:], Gamma_1 =RA_NSE_sum[5,:])
CSV.write(string(results_path_plots,sep,"NSE_dif_RA_no_climate_policy.csv"), df_RA_NSE_sum)
RA_NSE_PLOT=StatsPlots.groupedbar(RA_NSE_sum[1:end,1:end],bar_position =:stack, bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,RA_stack),title ="Non-Served Energy",label = ["Base Demand" "Low Demand" "High Demand" "Very High Demand" ],legendfont =font(8),legendfontsize=6, ylabel = "GWh",xlabel ="Level of Risk Aversion", color =[Palette_bright[1] Palette_bright[9] Palette_bright[3] Palette_bright[4]])
Plots.savefig(RA_NSE_PLOT,"Results/New_Differnt_RA_Plots_with_zones/Grouped_NSE_dif_RA_No_climate_policy.pdf")


# NSE CO2 Cap Policy
RA_nse_co2_cap_0 =  CSV.read(string(inputs_path_plots,sep,"NSE_sum_all_scen_RN_co2_cap.csv"), DataFrame, header=true)
RA_nse_co2_cap_025= CSV.read(string(inputs_path_plots,sep,"NSE_sum_all_scen_RA_0.25_co2_cap.csv"), DataFrame, header=true)
RA_nse_co2_cap_05 = CSV.read(string(inputs_path_plots,sep,"NSE_sum_all_scen_RA_co2_cap.csv"), DataFrame, header=true)
RA_nse_co2_cap_075 =CSV.read(string(inputs_path_plots,sep,"NSE_sum_all_scen_RA_0.75_co2_cap.csv"), DataFrame, header=true)
RA_nse_co2_cap_1 = CSV.read(string(inputs_path_plots,sep,"NSE_sum_all_scen_RA_1_co2_cap.csv"), DataFrame, header=true)

RA_NSE_co2_cap_sum = zeros(5,4)
for s in 1:4
    RA_NSE_co2_cap_sum[1,s] = sum(RA_nse_co2_cap_0[s,2:end])/1000
    RA_NSE_co2_cap_sum[2,s] = sum(RA_nse_co2_cap_025[s,2:end])/1000
    RA_NSE_co2_cap_sum[3,s] = sum(RA_nse_co2_cap_05[s,2:end])/1000
    RA_NSE_co2_cap_sum[4,s] = sum(RA_nse_co2_cap_075[s,2:end])/1000
    RA_NSE_co2_cap_sum[5,s] = sum(RA_nse_co2_cap_1[s,2:end])/1000
end

df_RA_NSE_co2_cap_sum = DataFrame(Scenario = RA_nse_co2_cap_0[:,1], Gamma_0 =RA_NSE_co2_cap_sum[1,:], Gamma_025 =RA_NSE_co2_cap_sum[2,:], Gamma_05 =RA_NSE_co2_cap_sum[3,:], Gamma_075 =RA_NSE_co2_cap_sum[4,:], Gamma_1 =RA_NSE_co2_cap_sum[5,:])
CSV.write(string(results_path_plots,sep,"NSE_dif_RA_co2_cap.csv"), df_RA_NSE_co2_cap_sum)
RA_NSE_co2_cap_PLOT=StatsPlots.groupedbar(RA_NSE_co2_cap_sum[1:end,1:end],bar_position =:stack, bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,RA_stack),title ="Non-Served Energy, CO2 Cap",label = ["Base Demand" "Low Demand" "High Demand" "Very High Demand" ],legendfont =font(8),legendfontsize=6, ylabel = "GWh",xlabel ="Level of Risk Aversion", color =[Palette_bright[1] Palette_bright[9] Palette_bright[3] Palette_bright[4]])
Plots.savefig(RA_NSE_co2_cap_PLOT,"Results/New_Differnt_RA_Plots_with_zones/Grouped_NSE_dif_RA_CO2_Cap_policy.pdf")



# NSE CO2 Tax Policy
RA_nse_co2_tax_0 =  CSV.read(string(inputs_path_plots,sep,"NSE_sum_all_scen_RN_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)
RA_nse_co2_tax_025= CSV.read(string(inputs_path_plots,sep,"NSE_sum_all_scen_RA_0.25_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)
RA_nse_co2_tax_05 = CSV.read(string(inputs_path_plots,sep,"NSE_sum_all_scen_RA_co2_tax_and_Nuclear_const.csv"), DataFrame, header=true)
RA_nse_co2_tax_075 =CSV.read(string(inputs_path_plots,sep,"NSE_sum_all_scen_RA_0.75_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)
RA_nse_co2_tax_1 = CSV.read(string(inputs_path_plots,sep,"NSE_sum_all_scen_RA_1_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)

RA_NSE_co2_tax_sum = zeros(5,4)
for s in 1:4
    RA_NSE_co2_tax_sum[1,s] = sum(RA_nse_co2_tax_0[s,2:end])/1000
    RA_NSE_co2_tax_sum[2,s] = sum(RA_nse_co2_tax_025[s,2:end])/1000
    RA_NSE_co2_tax_sum[3,s] = sum(RA_nse_co2_tax_05[s,2:end])/1000
    RA_NSE_co2_tax_sum[4,s] = sum(RA_nse_co2_tax_075[s,2:end])/1000
    RA_NSE_co2_tax_sum[5,s] = sum(RA_nse_co2_tax_1[s,2:end])/1000
end

df_RA_NSE_co2_tax_sum = DataFrame(Scenario = RA_nse_co2_tax_0[:,1], Gamma_0 =RA_NSE_co2_tax_sum[1,:], Gamma_025 =RA_NSE_co2_tax_sum[2,:], Gamma_05 =RA_NSE_co2_tax_sum[3,:], Gamma_075 =RA_NSE_co2_tax_sum[4,:], Gamma_1 =RA_NSE_co2_tax_sum[5,:])
CSV.write(string(results_path_plots,sep,"NSE_dif_RA_co2_tax.csv"), df_RA_NSE_co2_tax_sum)
RA_NSE_co2_tax_PLOT=StatsPlots.groupedbar(RA_NSE_co2_tax_sum[1:end,1:end],bar_position =:stack, bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,RA_stack),title ="Non-Served Energy, CO2 Tax",label = ["Base Demand" "Low Demand" "High Demand" "Very High Demand" ],legendfont =font(8),legendfontsize=6, ylabel = "GWh",xlabel ="Level of Risk Aversion", color =[Palette_bright[1] Palette_bright[9] Palette_bright[3] Palette_bright[4]])
Plots.savefig(RA_NSE_co2_tax_PLOT,"Results/New_Differnt_RA_Plots_with_zones/Grouped_NSE_dif_RA_CO2_Tax_policy.pdf")


#PLOT Revenues for different RA
#REVENUES NO CLIMATE POLICIES
RA_Revenue_0 =  CSV.read(string(inputs_path_plots,sep,"Revenues_all_scenarios_RN_No_co2_policies.csv"), DataFrame, header=true)
RA_Revenue_025 =  CSV.read(string(inputs_path_plots,sep,"Revenues_all_scenarios_RA_0.25_no_co2_policy.csv"), DataFrame, header=true)
RA_Revenue_05 =  CSV.read(string(inputs_path_plots,sep,"Revenues_all_scenarios_RA_No_co2_policies.csv"), DataFrame, header=true)
RA_Revenue_075 =  CSV.read(string(inputs_path_plots,sep,"Revenues_all_scenarios_RA_0.75_no_co2_policy.csv"), DataFrame, header=true)
RA_Revenue_1 =  CSV.read(string(inputs_path_plots,sep,"Revenues_all_scenarios_RA_1_no_co2_policy.csv"), DataFrame, header=true)
RA_Revenue_sum =zeros(5,S)
for s in 1:4
    RA_Revenue_sum[1,s] = sum(RA_Revenue_0[:,s+1])/1e9
    RA_Revenue_sum[2,s] = sum(RA_Revenue_025[:,s+1])/1e9
    RA_Revenue_sum[3,s] = sum(RA_Revenue_05[:,s+1])/1e9
    RA_Revenue_sum[4,s] = sum(RA_Revenue_075[:,s+1])/1e9
    RA_Revenue_sum[5,s] = sum(RA_Revenue_1[:,s+1])/1e9
end
trans_RA_rev_sum = transpose(RA_Revenue_sum)
df_RA_rev_sum = DataFrame(Scenarios = Scenario_stack, gamma_0 = trans_RA_rev_sum[:,1], gamma_025 = trans_RA_rev_sum[:,2], gamma_05 = trans_RA_rev_sum[:,3], gamma_075 = trans_RA_rev_sum[:,4], gamma_1 = trans_RA_rev_sum[:,5])
CSV.write(string(results_path_plots,sep,"Revenue_sum_dif_RA_no_policy.csv"), df_RA_rev_sum)
RA_Revenue_PLOT=StatsPlots.groupedbar(RA_Revenue_sum[1:end,1:end],bar_position =:stack, bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,RA_stack),title ="Revenues, No Climate Policies",label = ["Base Demand" "Low Demand" "High Demand" "Very High Demand" ],legendfont =font(8),legendfontsize=6, ylabel = "Billion €",xlabel ="Level of Risk Aversion", color =[Palette_bright[1] Palette_bright[9] Palette_bright[3] Palette_bright[4]])
Plots.savefig(RA_Revenue_PLOT,"Results/New_Differnt_RA_Plots_with_zones/Grouped_Revenue_dif_RA_No_climate_policy.pdf")

#REVENUES CO2 CAP
RA_Revenue_Co2_cap_0 =  CSV.read(string(inputs_path_plots,sep,"Revenues_all_scenarios_RN_co2_cap.csv"), DataFrame, header=true)
RA_Revenue_Co2_cap_025 =  CSV.read(string(inputs_path_plots,sep,"Revenues_all_scenarios_RA_0.25_co2_cap.csv"), DataFrame, header=true)
RA_Revenue_Co2_cap_05 =  CSV.read(string(inputs_path_plots,sep,"Revenues_all_scenarios_RA_co2_cap.csv"), DataFrame, header=true)
RA_Revenue_Co2_cap_075 =  CSV.read(string(inputs_path_plots,sep,"Revenues_all_scenarios_RA_0.75_co2_cap.csv"), DataFrame, header=true)
RA_Revenue_Co2_cap_1 =  CSV.read(string(inputs_path_plots,sep,"Revenues_all_scenarios_RA_1_co2_cap.csv"), DataFrame, header=true)
RA_Revenue_co2_cap_sum =zeros(5,S)
for s in 1:4
    RA_Revenue_co2_cap_sum[1,s] = sum(RA_Revenue_Co2_cap_0[:,s+1])/1e9
    RA_Revenue_co2_cap_sum[2,s] = sum(RA_Revenue_Co2_cap_025[:,s+1])/1e9
    RA_Revenue_co2_cap_sum[3,s] = sum(RA_Revenue_Co2_cap_05[:,s+1])/1e9
    RA_Revenue_co2_cap_sum[4,s] = sum(RA_Revenue_Co2_cap_075[:,s+1])/1e9
    RA_Revenue_co2_cap_sum[5,s] = sum(RA_Revenue_Co2_cap_1[:,s+1])/1e9
end
trans_RA_rev_sum_co2_cap = transpose(RA_Revenue_co2_cap_sum)
df_RA_rev_sum_co2_cap = DataFrame(Scenarios = Scenario_stack, gamma_0 = trans_RA_rev_sum_co2_cap[:,1], gamma_025 = trans_RA_rev_sum_co2_cap[:,2], gamma_05 = trans_RA_rev_sum_co2_cap[:,3], gamma_075 = trans_RA_rev_sum_co2_cap[:,4], gamma_1 = trans_RA_rev_sum_co2_cap[:,5])
CSV.write(string(results_path_plots,sep,"Revenue_sum_dif_RA_co2_cap.csv"), df_RA_rev_sum_co2_cap)
RA_Revenue_co2_cap_PLOT=StatsPlots.groupedbar(RA_Revenue_co2_cap_sum[1:end,1:end],bar_position =:stack, bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,RA_stack),title ="Revenues, CO2 Cap",label = ["Base Demand" "Low Demand" "High Demand" "Very High Demand" ],legendfont =font(8),legendfontsize=6, ylabel = "Billion €",xlabel ="Level of Risk Aversion", color =[Palette_bright[1] Palette_bright[9] Palette_bright[3] Palette_bright[4]])
Plots.savefig(RA_Revenue_co2_cap_PLOT,"Results/New_Differnt_RA_Plots_with_zones/Grouped_Revenue_dif_RA_co2_cap_policy.pdf")


#REVENUES CO2 Tax
RA_Revenue_co2_tax_0 =  CSV.read(string(inputs_path_plots,sep,"Revenues_all_scenarios_RN_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)
RA_Revenue_co2_tax_025 =  CSV.read(string(inputs_path_plots,sep,"Revenues_all_scenarios_RA_0.25_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)
RA_Revenue_co2_tax_05 =  CSV.read(string(inputs_path_plots,sep,"Revenues_all_scenarios_RA_co2_tax_and_Nuclear_const.csv"), DataFrame, header=true)
RA_Revenue_co2_tax_075 =  CSV.read(string(inputs_path_plots,sep,"Revenues_all_scenarios_RA_0.75_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)
RA_Revenue_co2_tax_1 =  CSV.read(string(inputs_path_plots,sep,"Revenues_all_scenarios_RA_1_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)
RA_Revenue_co2_tax_sum =zeros(5,S)
for s in 1:4
    RA_Revenue_co2_tax_sum[1,s] = sum(RA_Revenue_co2_tax_0[:,s+1])/1e9
    RA_Revenue_co2_tax_sum[2,s] = sum(RA_Revenue_co2_tax_025[:,s+1])/1e9
    RA_Revenue_co2_tax_sum[3,s] = sum(RA_Revenue_co2_tax_05[:,s+1])/1e9
    RA_Revenue_co2_tax_sum[4,s] = sum(RA_Revenue_co2_tax_075[:,s+1])/1e9
    RA_Revenue_co2_tax_sum[5,s] = sum(RA_Revenue_co2_tax_1[:,s+1])/1e9
end
trans_RA_rev_sum_co2_tax = transpose(RA_Revenue_co2_tax_sum)
df_RA_rev_sum_co2_tax = DataFrame(Scenarios = Scenario_stack, gamma_0 = trans_RA_rev_sum_co2_tax[:,1], gamma_025 = trans_RA_rev_sum_co2_tax[:,2], gamma_05 = trans_RA_rev_sum_co2_tax[:,3], gamma_075 = trans_RA_rev_sum_co2_tax[:,4], gamma_1 = trans_RA_rev_sum_co2_tax[:,5])
CSV.write(string(results_path_plots,sep,"Revenue_sum_dif_RA_co2_tax.csv"), df_RA_rev_sum_co2_tax)
RA_Revenue_co2_tax_PLOT=StatsPlots.groupedbar(RA_Revenue_co2_tax_sum[1:end,1:end],bar_position =:stack, bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,RA_stack),title ="Revenues, CO2 Tax",label = ["Base Demand" "Low Demand" "High Demand" "Very High Demand" ],legendfont =font(8),legendfontsize=6, ylabel = "Billion €",xlabel ="Level of Risk Aversion", color =[Palette_bright[1] Palette_bright[9] Palette_bright[3] Palette_bright[4]])
Plots.savefig(RA_Revenue_co2_tax_PLOT,"Results/New_Differnt_RA_Plots_with_zones/Grouped_Revenue_dif_RA_co2_tax_policy.pdf")


#POWER PRICE FOR DIFFERNT RA
#POWER PRICE NO CLIMATE POLICIES ZONE 1

RA_Price_Z1_0 =  CSV.read(string(inputs_path_plots,sep,"Z1_price_all_scen_RN_No_co2_policies.csv"), DataFrame, header=true)
RA_Price_Z1_025 =  CSV.read(string(inputs_path_plots,sep,"Z1_price_all_scen_RA_0.25_no_co2_policy.csv"), DataFrame, header=true)
RA_Price_Z1_05 =  CSV.read(string(inputs_path_plots,sep,"Z1_price_all_scen_RA_No_co2_policies.csv"), DataFrame, header=true)
RA_Price_Z1_075 =  CSV.read(string(inputs_path_plots,sep,"Z1_price_all_scen_RA_0.75_no_co2_policy.csv"), DataFrame, header=true)
RA_Price_Z1_1 =  CSV.read(string(inputs_path_plots,sep,"Z1_price_all_scen_RA_1_no_co2_policy.csv"), DataFrame, header=true)

RA_Z1_Max_Price = zeros(5,4)
RA_Z1_Mean_Price = zeros(5,4)
for s in 1:4
    RA_Z1_Max_Price[1,s] = maximum(RA_Price_Z1_0[:,s+1])
    RA_Z1_Max_Price[2,s] = maximum(RA_Price_Z1_025[:,s+1])
    RA_Z1_Max_Price[3,s] = maximum(RA_Price_Z1_05[:,s+1])
    RA_Z1_Max_Price[4,s] = maximum(RA_Price_Z1_075[:,s+1])
    RA_Z1_Max_Price[5,s] = maximum(RA_Price_Z1_1[:,s+1])

    RA_Z1_Mean_Price[1,s] = mean(RA_Price_Z1_0[:,s+1])
    RA_Z1_Mean_Price[2,s] = mean(RA_Price_Z1_025[:,s+1])
    RA_Z1_Mean_Price[3,s] = mean(RA_Price_Z1_05[:,s+1])
    RA_Z1_Mean_Price[4,s] = mean(RA_Price_Z1_075[:,s+1])
    RA_Z1_Mean_Price[5,s] = mean(RA_Price_Z1_1[:,s+1])

end

df_max_price_z1 = DataFrame(Level_of_RA = RA_stack[:], Base_Demand = RA_Z1_Max_Price[:,1], Low_Demand = RA_Z1_Max_Price[:,2],High_Demand = RA_Z1_Max_Price[:,3],Very_High_Demand = RA_Z1_Max_Price[:,4])
CSV.write(string(results_path_plots,sep,"Max_price_Z1_dif_RA.csv"), df_max_price_z1)
df_Mean_price_z1 = DataFrame(Level_of_RA = RA_stack[:], Base_Demand = RA_Z1_Mean_Price[:,1], Low_Demand = RA_Z1_Mean_Price[:,2],High_Demand = RA_Z1_Mean_Price[:,3],Very_High_Demand = RA_Z1_Mean_Price[:,4])
CSV.write(string(results_path_plots,sep,"Mean_price_Z1_dif_RA.csv"), df_Mean_price_z1)
plot_max_price_Z1 = StatsPlots.groupedbar(RA_Z1_Max_Price[:,1:end],bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,RA_stack),title ="Maximum Power Price",label = ["Base Demand" "Low Demand" "High Demand" "Very High Demand" ],legendfont =font(8),legendfontsize=6, ylabel = "EUR/MWh",xlabel ="Level of Risk Aversion", color =[Palette_bright[1] Palette_bright[9] Palette_bright[3] Palette_bright[4]])
plot_mean_price_Z1 =StatsPlots.groupedbar(RA_Z1_Mean_Price[:,1:end],bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,RA_stack),title ="Mean Power Price",label = ["Base Demand" "Low Demand" "High Demand" "Very High Demand" ],legendfont =font(8),legendfontsize=6, ylabel = "EUR/MWh",xlabel ="Level of Risk Aversion", color =[Palette_bright[1] Palette_bright[9] Palette_bright[3] Palette_bright[4]])
Plots.savefig(plot_max_price_Z1,"Results/New_Differnt_RA_Plots_with_zones/Grouped_Max_Price_Z1_dif_RA.pdf")
Plots.savefig(plot_mean_price_Z1,"Results/New_Differnt_RA_Plots_with_zones/Grouped_Mean_Price_Z1_dif_RA.pdf")


#POWER PRICE WITH CO2 CAP ZONE 1

RA_Price_Z1_co2_cap_0 =  CSV.read(string(inputs_path_plots,sep,"Z1_price_all_scen_RN_co2_cap.csv"), DataFrame, header=true)
RA_Price_Z1_co2_cap_025 =  CSV.read(string(inputs_path_plots,sep,"Z1_price_all_scen_RA_0.25_co2_cap.csv"), DataFrame, header=true)
RA_Price_Z1_co2_cap_05 =  CSV.read(string(inputs_path_plots,sep,"Z1_price_all_scen_RA_co2_cap.csv"), DataFrame, header=true)
RA_Price_Z1_co2_cap_075 =  CSV.read(string(inputs_path_plots,sep,"Z1_price_all_scen_RA_0.75_co2_cap.csv"), DataFrame, header=true)
RA_Price_Z1_co2_cap_1 =  CSV.read(string(inputs_path_plots,sep,"Z1_price_all_scen_RA_1_co2_cap.csv"), DataFrame, header=true)

RA_Z1_co2_cap_Max_Price = zeros(5,4)
RA_Z1_co2_cap_Mean_Price = zeros(5,4)
for s in 1:4
    RA_Z1_co2_cap_Max_Price[1,s] = maximum(RA_Price_Z1_co2_cap_0[:,s+1])
    RA_Z1_co2_cap_Max_Price[2,s] = maximum(RA_Price_Z1_co2_cap_025[:,s+1])
    RA_Z1_co2_cap_Max_Price[3,s] = maximum(RA_Price_Z1_co2_cap_05[:,s+1])
    RA_Z1_co2_cap_Max_Price[4,s] = maximum(RA_Price_Z1_co2_cap_075[:,s+1])
    RA_Z1_co2_cap_Max_Price[5,s] = maximum(RA_Price_Z1_co2_cap_1[:,s+1])

    RA_Z1_co2_cap_Mean_Price[1,s] = mean(RA_Price_Z1_co2_cap_0[:,s+1])
    RA_Z1_co2_cap_Mean_Price[2,s] = mean(RA_Price_Z1_co2_cap_025[:,s+1])
    RA_Z1_co2_cap_Mean_Price[3,s] = mean(RA_Price_Z1_co2_cap_05[:,s+1])
    RA_Z1_co2_cap_Mean_Price[4,s] = mean(RA_Price_Z1_co2_cap_075[:,s+1])
    RA_Z1_co2_cap_Mean_Price[5,s] = mean(RA_Price_Z1_co2_cap_1[:,s+1])

end

df_max_price_Z1_co2_cap = DataFrame(Level_of_RA = RA_stack[:], Base_Demand = RA_Z1_co2_cap_Max_Price[:,1], Low_Demand = RA_Z1_co2_cap_Max_Price[:,2],High_Demand = RA_Z1_co2_cap_Max_Price[:,3],Very_High_Demand = RA_Z1_co2_cap_Max_Price[:,4])
CSV.write(string(results_path_plots,sep,"Max_price_Z1_co2_cap_dif_RA.csv"), df_max_price_Z1_co2_cap)
df_Mean_price_Z1_co2_cap = DataFrame(Level_of_RA = RA_stack[:], Base_Demand = RA_Z1_co2_cap_Mean_Price[:,1], Low_Demand = RA_Z1_co2_cap_Mean_Price[:,2],High_Demand = RA_Z1_co2_cap_Mean_Price[:,3],Very_High_Demand = RA_Z1_co2_cap_Mean_Price[:,4])
CSV.write(string(results_path_plots,sep,"Mean_price_Z1_co2_cap_dif_RA.csv"), df_Mean_price_Z1_co2_cap)
plot_max_price_Z1_co2_cap = StatsPlots.groupedbar(RA_Z1_co2_cap_Max_Price[:,1:end],bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,RA_stack),title ="Maximum Power Price, CO2 Cap",label = ["Base Demand" "Low Demand" "High Demand" "Very High Demand" ],legendfont =font(8),legendfontsize=6, ylabel = "EUR/MWh",xlabel ="Level of Risk Aversion", color =[Palette_bright[1] Palette_bright[9] Palette_bright[3] Palette_bright[4]])
plot_mean_price_Z1_co2_cap =StatsPlots.groupedbar(RA_Z1_co2_cap_Mean_Price[:,1:end],bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,RA_stack),title ="Mean Power Price, CO2 Cap",label = ["Base Demand" "Low Demand" "High Demand" "Very High Demand" ],legendfont =font(8),legendfontsize=6, ylabel = "EUR/MWh",xlabel ="Level of Risk Aversion", color =[Palette_bright[1] Palette_bright[9] Palette_bright[3] Palette_bright[4]])
Plots.savefig(plot_max_price_Z1_co2_cap,"Results/New_Differnt_RA_Plots_with_zones/Grouped_Max_Price_co2_cap_Z1_dif_RA.pdf")
Plots.savefig(plot_mean_price_Z1_co2_cap,"Results/New_Differnt_RA_Plots_with_zones/Grouped_Mean_Price_co2_cap_Z1_dif_RA.pdf")


#POWER PRICE WITH CO2 TAX ZONE 1 (but easy to switch to other zones)

RA_Price_Z1_co2_Tax_0 =  CSV.read(string(inputs_path_plots,sep,"Z1_price_all_scen_RN_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)
RA_Price_Z1_co2_Tax_025 =  CSV.read(string(inputs_path_plots,sep,"Z1_price_all_scen_RA_0.25_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)
RA_Price_Z1_co2_Tax_05 =  CSV.read(string(inputs_path_plots,sep,"Z1_price_all_scen_RA_co2_tax_and_Nuclear_const.csv"), DataFrame, header=true)
RA_Price_Z1_co2_Tax_075 =  CSV.read(string(inputs_path_plots,sep,"Z1_price_all_scen_RA_0.75_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)
RA_Price_Z1_co2_Tax_1 =  CSV.read(string(inputs_path_plots,sep,"Z1_price_all_scen_RA_1_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)

RA_Z1_co2_Tax_Max_Price = zeros(5,4)
RA_Z1_co2_Tax_Mean_Price = zeros(5,4)
for s in 1:4
    RA_Z1_co2_Tax_Max_Price[1,s] = maximum(RA_Price_Z1_co2_Tax_0[:,s+1])
    RA_Z1_co2_Tax_Max_Price[2,s] = maximum(RA_Price_Z1_co2_Tax_025[:,s+1])
    RA_Z1_co2_Tax_Max_Price[3,s] = maximum(RA_Price_Z1_co2_Tax_05[:,s+1])
    RA_Z1_co2_Tax_Max_Price[4,s] = maximum(RA_Price_Z1_co2_Tax_075[:,s+1])
    RA_Z1_co2_Tax_Max_Price[5,s] = maximum(RA_Price_Z1_co2_Tax_1[:,s+1])

    RA_Z1_co2_Tax_Mean_Price[1,s] = mean(RA_Price_Z1_co2_Tax_0[:,s+1])
    RA_Z1_co2_Tax_Mean_Price[2,s] = mean(RA_Price_Z1_co2_Tax_025[:,s+1])
    RA_Z1_co2_Tax_Mean_Price[3,s] = mean(RA_Price_Z1_co2_Tax_05[:,s+1])
    RA_Z1_co2_Tax_Mean_Price[4,s] = mean(RA_Price_Z1_co2_Tax_075[:,s+1])
    RA_Z1_co2_Tax_Mean_Price[5,s] = mean(RA_Price_Z1_co2_Tax_1[:,s+1])

end

df_max_price_Z1_co2_Tax = DataFrame(Level_of_RA = RA_stack[:], Base_Demand = RA_Z1_co2_Tax_Max_Price[:,1], Low_Demand = RA_Z1_co2_Tax_Max_Price[:,2],High_Demand = RA_Z1_co2_Tax_Max_Price[:,3],Very_High_Demand = RA_Z1_co2_Tax_Max_Price[:,4])
CSV.write(string(results_path_plots,sep,"Max_price_Z1_co2_Tax_dif_RA.csv"), df_max_price_Z1_co2_Tax)
df_Mean_price_Z1_co2_Tax = DataFrame(Level_of_RA = RA_stack[:], Base_Demand = RA_Z1_co2_Tax_Mean_Price[:,1], Low_Demand = RA_Z1_co2_Tax_Mean_Price[:,2],High_Demand = RA_Z1_co2_Tax_Mean_Price[:,3],Very_High_Demand = RA_Z1_co2_Tax_Mean_Price[:,4])
CSV.write(string(results_path_plots,sep,"Mean_price_Z1_co2_Tax_dif_RA.csv"), df_Mean_price_Z1_co2_Tax)
plot_max_price_Z1_co2_Tax = StatsPlots.groupedbar(RA_Z1_co2_Tax_Max_Price[:,1:end],bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,RA_stack),title ="Maximum Power Price, CO2 Tax",label = ["Base Demand" "Low Demand" "High Demand" "Very High Demand" ],legendfont =font(8),legendfontsize=6, ylabel = "EUR/MWh",xlabel ="Level of Risk Aversion", color =[Palette_bright[1] Palette_bright[9] Palette_bright[3] Palette_bright[4]])
plot_mean_price_Z1_co2_Tax =StatsPlots.groupedbar(RA_Z1_co2_Tax_Mean_Price[:,1:end],bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,RA_stack),title ="Mean Power Price, CO2 Tax",label = ["Base Demand" "Low Demand" "High Demand" "Very High Demand" ],legendfont =font(8),legendfontsize=6, ylabel = "EUR/MWh",xlabel ="Level of Risk Aversion", color =[Palette_bright[1] Palette_bright[9] Palette_bright[3] Palette_bright[4]])
Plots.savefig(plot_max_price_Z1_co2_Tax,"Results/New_Differnt_RA_Plots_with_zones/Grouped_Max_Price_co2_Tax_Z1_dif_RA.pdf")
Plots.savefig(plot_mean_price_Z1_co2_Tax,"Results/New_Differnt_RA_Plots_with_zones/Grouped_Mean_Price_co2_Tax_Z1_dif_RA.pdf")



#EMISSIONS FOR DIFFERNT RA
#EMISSIONS NO CLIMATE POLICIES

RA_Emission_0 =  CSV.read(string(inputs_path_plots,sep,"Emission_sum_all_scen_RN_No_co2_policies.csv"), DataFrame, header=true)
RA_Emission_025 =  CSV.read(string(inputs_path_plots,sep,"Emission_sum_all_scen_RA_0.25_no_co2_policy.csv"), DataFrame, header=true)
RA_Emission_05 =  CSV.read(string(inputs_path_plots,sep,"Emission_sum_all_scen_RA_No_co2_policies.csv"), DataFrame, header=true)
RA_Emission_075 =  CSV.read(string(inputs_path_plots,sep,"Emission_sum_all_scen_RA_0.75_no_co2_policy.csv"), DataFrame, header=true)
RA_Emission_1 =  CSV.read(string(inputs_path_plots,sep,"Emission_sum_all_scen_RA_1_no_co2_policy.csv"), DataFrame, header=true)

RA_em_sum = zeros(5,S)
for s in 1:4
    RA_em_sum[1,s] = sum(RA_Emission_0[s,2:end])/1000000
    RA_em_sum[2,s] = sum(RA_Emission_025[s,2:end])/1000000
    RA_em_sum[3,s] = sum(RA_Emission_05[s,2:end])/1000000
    RA_em_sum[4,s] = sum(RA_Emission_075[s,2:end])/1000000
    RA_em_sum[5,s] = sum(RA_Emission_1[s,2:end])/1000000
end
trans_em_sum = transpose(RA_em_sum)
df_em_sum = DataFrame(Scenario = Scenario_stack, gamma_0 = trans_em_sum[:,1], gamma_025 = trans_em_sum[:,2], gamma_05 = trans_em_sum[:,3], gamma_075 = trans_em_sum[:,4], gamma_1 = trans_em_sum[:,5])
CSV.write(string(results_path_plots,sep,"Emissions_dif_RA_no_policy.csv"), df_em_sum)
plot_emission_per_scen = StatsPlots.groupedbar(RA_em_sum[:,1:end],bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,RA_stack),title ="Carbon Emissions",label = ["Base Demand" "Low Demand" "High Demand" "Very High Demand" ],legendfont =font(8),legendfontsize=6, ylabel = "Mt CO2",xlabel ="Level of Risk Aversion", color =[Palette_bright[1] Palette_bright[9] Palette_bright[3] Palette_bright[4]])
Plots.savefig(plot_emission_per_scen,"Results/New_Differnt_RA_Plots_with_zones/Grouped_Emission_dif_RA.pdf")

#EMISSIONS CO2 CAP

RA_Emission_co2_cap_0 =  CSV.read(string(inputs_path_plots,sep,"Emission_sum_all_scen_RN_co2_cap.csv"), DataFrame, header=true)
RA_Emission_co2_cap_025 =  CSV.read(string(inputs_path_plots,sep,"Emission_sum_all_scen_RA_0.25_co2_cap.csv"), DataFrame, header=true)
RA_Emission_co2_cap_05 =  CSV.read(string(inputs_path_plots,sep,"Emission_sum_all_scen_RA_co2_cap.csv"), DataFrame, header=true)
RA_Emission_co2_cap_075 =  CSV.read(string(inputs_path_plots,sep,"Emission_sum_all_scen_RA_0.75_co2_cap.csv"), DataFrame, header=true)
RA_Emission_co2_cap_1 =  CSV.read(string(inputs_path_plots,sep,"Emission_sum_all_scen_RA_1_co2_cap.csv"), DataFrame, header=true)

RA_em_co2_cap_sum = zeros(5,4)
for s in 1:4
    RA_em_co2_cap_sum[1,s] = sum(RA_Emission_co2_cap_0[s,2:end])/1000000
    RA_em_co2_cap_sum[2,s] = sum(RA_Emission_co2_cap_025[s,2:end])/1000000
    RA_em_co2_cap_sum[3,s] = sum(RA_Emission_co2_cap_05[s,2:end])/1000000
    RA_em_co2_cap_sum[4,s] = sum(RA_Emission_co2_cap_075[s,2:end])/1000000
    RA_em_co2_cap_sum[5,s] = sum(RA_Emission_co2_cap_1[s,2:end])/1000000
end
trans_em_co2_cap_sum = transpose(RA_em_co2_cap_sum)
df_em_co2_cap_sum = DataFrame(Scenario = Scenario_stack, gamma_0 = trans_em_co2_cap_sum[:,1], gamma_025 = trans_em_co2_cap_sum[:,2], gamma_05 = trans_em_co2_cap_sum[:,3], gamma_075 = trans_em_co2_cap_sum[:,4], gamma_1 = trans_em_co2_cap_sum[:,5])
CSV.write(string(results_path_plots,sep,"Emissions_dif_RA_Co2_Cap_policy.csv"), df_em_co2_cap_sum)
plot_emission_per_scen_co2_cap = StatsPlots.groupedbar(RA_em_co2_cap_sum[:,1:end],bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,RA_stack),title ="Carbon Emissions, CO2 Cap",label = ["Base Demand" "Low Demand" "High Demand" "Very High Demand" ],legendfont =font(8),legendfontsize=6, ylabel = "Mt CO2",xlabel ="Level of Risk Aversion", color =[Palette_bright[1] Palette_bright[9] Palette_bright[3] Palette_bright[4]])
Plots.savefig(plot_emission_per_scen_co2_cap,"Results/New_Differnt_RA_Plots_with_zones/Grouped_Emission_CO2_Cap_dif_RA.pdf")


#EMISSIONS CO2 TAX

RA_Emission_CO2_Tax_0 =  CSV.read(string(inputs_path_plots,sep,"Emission_sum_all_scen_RN_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)
RA_Emission_CO2_Tax_025 =  CSV.read(string(inputs_path_plots,sep,"Emission_sum_all_scen_RA_0.25_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)
RA_Emission_CO2_Tax_05 =  CSV.read(string(inputs_path_plots,sep,"Emission_sum_all_scen_RA_co2_tax_and_Nuclear_const.csv"), DataFrame, header=true)
RA_Emission_CO2_Tax_075 =  CSV.read(string(inputs_path_plots,sep,"Emission_sum_all_scen_RA_0.75_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)
RA_Emission_CO2_Tax_1 =  CSV.read(string(inputs_path_plots,sep,"Emission_sum_all_scen_RA_1_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)

RA_em_CO2_Tax_sum = zeros(5,4)
for s in 1:4
    RA_em_CO2_Tax_sum[1,s] = sum(RA_Emission_CO2_Tax_0[s,2:end])/1000000
    RA_em_CO2_Tax_sum[2,s] = sum(RA_Emission_CO2_Tax_025[s,2:end])/1000000
    RA_em_CO2_Tax_sum[3,s] = sum(RA_Emission_CO2_Tax_05[s,2:end])/1000000
    RA_em_CO2_Tax_sum[4,s] = sum(RA_Emission_CO2_Tax_075[s,2:end])/1000000
    RA_em_CO2_Tax_sum[5,s] = sum(RA_Emission_CO2_Tax_1[s,2:end])/1000000
end
trans_em_CO2_Tax_sum = transpose(RA_em_CO2_Tax_sum)
df_em_CO2_Tax_sum = DataFrame(Scenario = Scenario_stack, gamma_0 = trans_em_CO2_Tax_sum[:,1], gamma_025 = trans_em_CO2_Tax_sum[:,2], gamma_05 = trans_em_CO2_Tax_sum[:,3], gamma_075 = trans_em_CO2_Tax_sum[:,4], gamma_1 = trans_em_CO2_Tax_sum[:,5])
CSV.write(string(results_path_plots,sep,"Emissions_dif_RA_Co2_Tax_policy.csv"), df_em_CO2_Tax_sum)
plot_emission_per_scen_CO2_Tax = StatsPlots.groupedbar(RA_em_CO2_Tax_sum[:,1:end],bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,RA_stack),title ="Carbon Emissions, CO2 Tax",label = ["Base Demand" "Low Demand" "High Demand" "Very High Demand" ],legendfont =font(8),legendfontsize=6, ylabel = "Mt CO2",xlabel ="Level of Risk Aversion", color =[Palette_bright[1] Palette_bright[9] Palette_bright[3] Palette_bright[4]])
Plots.savefig(plot_emission_per_scen_CO2_Tax,"Results/New_Differnt_RA_Plots_with_zones/Grouped_Emission_CO2_Tax_dif_RA.pdf")

#PLOT GENERATION DUFFERENT LEVELS OF RA

#PLOT GENERATION FOR CO2 TAX
RA_gen_all_scen_co2_tax_0 = CSV.read(string(inputs_path_plots,sep,"Generation_all_scenarios_RN_co2_cap.csv"), DataFrame, header=true)
RA_gen_all_scen_co2_tax_025 = CSV.read(string(inputs_path_plots,sep,"Generation_all_scenarios_RA_0.25_co2_cap.csv"), DataFrame, header=true)
RA_gen_all_scen_co2_tax_05 = CSV.read(string(inputs_path_plots,sep,"Generation_all_scenarios_RA_co2_cap.csv"), DataFrame, header=true)
RA_gen_all_scen_co2_tax_075 = CSV.read(string(inputs_path_plots,sep,"Generation_all_scenarios_RA_0.75_co2_cap.csv"), DataFrame, header=true)
RA_gen_all_scen_co2_tax_1 = CSV.read(string(inputs_path_plots,sep,"Generation_all_scenarios_RA_1_co2_cap.csv"), DataFrame, header=true)

RA_gen_all_in_TWh = zeros(4,R,5)
for s in 1:4
    for r in 1:R
        RA_gen_all_in_TWh[s,r,1] = RA_gen_all_scen_co2_tax_0[s,r]/1000000
        RA_gen_all_in_TWh[s,r,2] = RA_gen_all_scen_co2_tax_025[s,r]/1000000
        RA_gen_all_in_TWh[s,r,3] = RA_gen_all_scen_co2_tax_05[s,r]/1000000
        RA_gen_all_in_TWh[s,r,4] = RA_gen_all_scen_co2_tax_075[s,r]/1000000
        RA_gen_all_in_TWh[s,r,5] = RA_gen_all_scen_co2_tax_1[s,r]/1000000
    end
end

plt_RA_gen_all_scen_co2_tax_0 = StatsPlots.groupedbar(RA_gen_all_in_TWh[:,:,1],bar_position =:stack, bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,demand_stack),title ="Generation, CO2 Cap, γ=0",label =["Nuclear" "Coal" "Gas" "Onshore Wind" "Offshore Wind Fixed" "Offshore Wind Float" "Solar" "Storage"] ,legendfont =font(8),legendfontsize=6, ylabel = "TWh",xlabel ="Demand Scenario", color =[Palette[3] Palette_dark[11] Palette_dark[8] Palette[1] Palette[2] Palette[9] Palette[11] Palette[6]]) 
plt_RA_gen_all_scen_co2_tax_025 = StatsPlots.groupedbar(RA_gen_all_in_TWh[:,:,2],bar_position =:stack, bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,demand_stack),title ="Generation, CO2 Cap, γ=0.25",label =["Nuclear" "Coal" "Gas" "Onshore Wind" "Offshore Wind Fixed" "Offshore Wind Float" "Solar" "Storage"] ,legendfont =font(8),legendfontsize=6, ylabel = "TWh",xlabel ="Demand Scenario", color =[Palette[3] Palette_dark[11] Palette_dark[8] Palette[1] Palette[2] Palette[9] Palette[11] Palette[6]]) 
plt_RA_gen_all_scen_co2_tax_05 = StatsPlots.groupedbar(RA_gen_all_in_TWh[:,:,3],bar_position =:stack, bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,demand_stack),title ="Generation, CO2 Cap, γ=0.5",label =["Nuclear" "Coal" "Gas" "Onshore Wind" "Offshore Wind Fixed" "Offshore Wind Float" "Solar" "Storage"] ,legendfont =font(8),legendfontsize=6, ylabel = "TWh",xlabel ="Demand Scenario", color =[Palette[3] Palette_dark[11] Palette_dark[8] Palette[1] Palette[2] Palette[9] Palette[11] Palette[6]]) 
plt_RA_gen_all_scen_co2_tax_075 = StatsPlots.groupedbar(RA_gen_all_in_TWh[:,:,4],bar_position =:stack, bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,demand_stack),title ="Generation, CO2 Cap, γ=0.75",label =["Nuclear" "Coal" "Gas" "Onshore Wind" "Offshore Wind Fixed" "Offshore Wind Float" "Solar" "Storage"] ,legendfont =font(8),legendfontsize=6, ylabel = "TWh",xlabel ="Demand Scenario", color =[Palette[3] Palette_dark[11] Palette_dark[8] Palette[1] Palette[2] Palette[9] Palette[11] Palette[6]]) 
plt_RA_gen_all_scen_co2_tax_1 = StatsPlots.groupedbar(RA_gen_all_in_TWh[:,:,5],bar_position =:stack, bar_width=0.9, xlim = (0.2,7.8), xticks =(1:8,demand_stack),title ="Generation, CO2 Cap, γ=1",label =["Nuclear" "Coal" "Gas" "Onshore Wind" "Offshore Wind Fixed" "Offshore Wind Float" "Solar" "Storage"] ,legendfont =font(8),legendfontsize=6, ylabel = "TWh",xlabel ="Demand Scenario", color =[Palette[3] Palette_dark[11] Palette_dark[8] Palette[1] Palette[2] Palette[9] Palette[11] Palette[6]]) 



#CO2 Price plot and csv, only fro co2 cap policy
RA_0_co2_price = CSV.read(string(inputs_path_plots,sep,"co2_price_RN_co2_cap.csv"), DataFrame, header=true)
RA_025_co2_price = CSV.read(string(inputs_path_plots,sep,"co2_price_RA_0.25_co2_cap.csv"), DataFrame, header=true)
RA_05_co2_price = CSV.read(string(inputs_path_plots,sep,"co2_price_RA_co2_cap.csv"), DataFrame, header=true)
RA_075_co2_price = CSV.read(string(inputs_path_plots,sep,"co2_price_RA_0.75_co2_cap.csv"), DataFrame, header=true)
RA_1_co2_price = CSV.read(string(inputs_path_plots,sep,"co2_price_RA_1_co2_cap.csv"), DataFrame, header=true)
RA_all_co2_price = zeros(S,5)
for s in 1:S
    RA_all_co2_price[s,1] = RA_0_co2_price[s,2]
    RA_all_co2_price[s,2] = RA_025_co2_price[s,2]
    RA_all_co2_price[s,3] = RA_05_co2_price[s,2]
    RA_all_co2_price[s,4] = RA_075_co2_price[s,2]
    RA_all_co2_price[s,5] = RA_1_co2_price[s,2]
end

df_co2_price_dif_RA = DataFrame(Scenarios = Scenario_stack, gamma_0 = RA_all_co2_price[:,1], gamma_025 = RA_all_co2_price[:,2], gamma_05 = RA_all_co2_price[:,3], gamma_075 = RA_all_co2_price[:,4], gamma_1 = RA_all_co2_price[:,5])
CSV.write(string(results_path_plots,sep,"CO2_PRICE_dif_RA_Co2_Cax_policy.csv"), df_co2_price_dif_RA)



#PLOTTING DURATION CURVE
#FOR RISK NEUTRAL MODEL
demand_z1 = CSV.read(string(inputs_path,sep,"demand_z1.csv"), DataFrame, header=true)
demand_z2 = CSV.read(string(inputs_path,sep,"demand_z2.csv"), DataFrame, header=true)
demand_base = zeros(T)
demand_low = zeros(T)
demand_high = zeros(T)
demand_very_high = zeros(T)
for t in 1:T
    demand_base[t] = (demand_z1[t,2]+demand_z2[t,2])/1000
    demand_low[t] = (demand_z1[t,3]+demand_z2[t,3])/1000
    demand_high[t] = (demand_z1[t,4]+demand_z2[t,4])/1000
    demand_very_high[t] = (demand_z1[t,5]+demand_z2[t,5])/1000
end
generation_z1 = CSV.read(string(inputs_path_plots,sep,"generation_z1_CN_RN_No_co2_policies.csv"), DataFrame, header=true) #For base demand
generation_z2 = CSV.read(string(inputs_path_plots,sep,"generation_z2_SK_RN_No_co2_policies.csv"), DataFrame, header=true)
generation_z3 = CSV.read(string(inputs_path_plots,sep,"generation_z3_GB_RN_No_co2_policies.csv"), DataFrame, header=true)
generation_nuclear = zeros(T)
generation_coal = zeros(T)
generation_gas = zeros(T)
generation_on_wind = zeros(T)
generation_off_wind_fix = zeros(T)
generation_off_wind_float = zeros(T)
generation_solar = zeros(T)
generation_VRES = zeros(T)
for t in 1:T
    generation_nuclear[t] = (generation_z1[t,2] + generation_z2[t,2] + generation_z3[t,2])/1000
    generation_coal[t] = (generation_z1[t,3] + generation_z2[t,3] + generation_z3[t,3])/1000
    generation_gas[t] = (generation_z1[t,4] + generation_z2[t,4] + generation_z3[t,4])/1000
    generation_on_wind[t] = (generation_z1[t,5] + generation_z2[t,5] + generation_z3[t,5])/1000
    generation_off_wind_fix[t] = (generation_z1[t,6] + generation_z2[t,6] + generation_z3[t,6])/1000
    generation_off_wind_float[t] = (generation_z1[t,7] + generation_z2[t,7] + generation_z3[t,7])/1000
    generation_solar[t] = (generation_z1[t,8] + generation_z2[t,8] + generation_z3[t,8])/1000
end
for t in 1:T
    generation_VRES[t] = generation_on_wind[t] + generation_off_wind_fix[t] + generation_off_wind_float[t] + generation_solar[t]
end

dc_demand_base = sort(demand_base, rev= true)    
dc_generation_nuclear = sort(generation_nuclear, rev=true)
dc_generation_coal = sort(generation_coal, rev=true)
dc_generation_gas = sort(generation_gas, rev=true)
dc_generation_on_wind = sort(generation_on_wind, rev=true)
dc_generation_off_wind_fix = sort(generation_off_wind_fix, rev=true)
dc_generation_off_wind_float = sort(generation_off_wind_float, rev=true)
dc_generation_solar = sort(generation_solar, rev=true)
dc_generation_VRES = sort(generation_VRES, rev=true)

net_load = zeros(T)
for t in 1:T
    net_load[t] = demand_base[t] - generation_VRES[t]
end
dc_net_load = sort(net_load, rev=true)

DC_plot_gen = Plots.plot([dc_demand_base,dc_generation_VRES,dc_net_load,dc_generation_coal,dc_generation_gas, dc_generation_on_wind,dc_generation_off_wind_fix,dc_generation_solar],tickfontsize =10,title = "Duration Curve, Risk Neutral, No Policy",label = ["Demand" "VRES" "Net Load" "Coal" "Gas" "Onshore wind" "Offshore wind fixed" "Solar"],xlabel = "Time [h]",ylabel = "GWh",legendfont =font(8),legendfontsize=6, color =[Palette_dark[2] Palette_green[7] Palette[8] Palette_dark[11] Palette_dark[8] Palette[1] Palette[2] palette_clear[6]])
Plots.savefig(DC_plot_gen,"Results/Results_DC_plots/DC_plot_RN_no_policy_base_demand.pdf")

#RISK NEUTRAL GENERATION
generation_base_RN_no_policy = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_Base_demand_scen_RN_No_co2_policies.csv"), DataFrame, header=true) #For base demand
generation_low_RN_no_policy = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_low_demand_scen_RN_No_co2_policies.csv"), DataFrame, header=true)
generation_high_RN_no_policy = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_high_demand_scen_RN_No_co2_policies.csv"), DataFrame, header=true)
generation_Very_high_RN_no_policy = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_Very_high_demand_scen_RN_No_co2_policies.csv"), DataFrame, header=true)

generation_base_RN_co2_cap = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_Base_demand_scen_RN_co2_cap.csv"), DataFrame, header=true) #For base demand
generation_low_RN_co2_cap = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_low_demand_scen_RN_co2_cap.csv"), DataFrame, header=true)
generation_high_RN_co2_cap = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_high_demand_scen_RN_co2_cap.csv"), DataFrame, header=true)
generation_Very_high_RN_co2_cap = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_Very_high_demand_scen_RN_co2_cap.csv"), DataFrame, header=true)

generation_base_RN_co2_tax = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_Base_demand_scen_RN_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true) #For base demand
generation_low_RN_co2_tax = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_low_demand_scen_RN_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)
generation_high_RN_co2_tax = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_high_demand_scen_RN_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)
generation_Very_high_RN_co2_tax = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_Very_high_demand_scen_RN_co2_tax_with_Nuclear_const.csv"), DataFrame, header=true)



#DURATION CURVCE AND STORAGE PLOTS FOR RA GAMMA=0.5 WITH CO2 CAP POLICY
generation_base_RA_05_no_policy = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_Base_demand_scen_RA_No_co2_policies.csv"), DataFrame, header=true) #For base demand
generation_low_RA_05_no_policy = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_low_demand_scen_RA_No_co2_policies.csv"), DataFrame, header=true)
generation_high_RA_05_no_policy = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_high_demand_scen_RA_No_co2_policies.csv"), DataFrame, header=true)
generation_Very_high_RA_05_no_policy = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_Very_high_demand_scen_RA_No_co2_policies.csv"), DataFrame, header=true)

generation_base_RA_05_co2_cap = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_Base_demand_scen_RA_co2_cap.csv"), DataFrame, header=true) #For base demand
generation_low_RA_05_co2_cap = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_low_demand_scen_RA_co2_cap.csv"), DataFrame, header=true)
generation_high_RA_05_co2_cap = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_high_demand_scen_RA_co2_cap.csv"), DataFrame, header=true)
generation_Very_high_RA_05_co2_cap = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_Very_high_demand_scen_RA_co2_cap.csv"), DataFrame, header=true)

generation_base_RA_05_co2_tax = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_Base_demand_scen_RA_co2_tax_and_Nuclear_const.csv"), DataFrame, header=true) #For base demand
generation_low_RA_05_co2_tax = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_low_demand_scen_RA_co2_tax_and_Nuclear_const.csv"), DataFrame, header=true)
generation_high_RA_05_co2_tax = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_high_demand_scen_RA_co2_tax_and_Nuclear_const.csv"), DataFrame, header=true)
generation_Very_high_RA_05_co2_tax = CSV.read(string(inputs_path_plots,sep,"All_Zones_gen_Very_high_demand_scen_RA_co2_tax_and_Nuclear_const.csv"), DataFrame, header=true)


generation_all_RA = zeros(T,R,S)
for t in 1:T, r in 1:R
    generation_all_RA[t,r,1] = generation_base_RA_05_co2_tax[t,r]
    generation_all_RA[t,r,2] = generation_low_RA_05_co2_tax[t,r]
    generation_all_RA[t,r,3] = generation_high_RA_05_co2_tax[t,r]
    generation_all_RA[t,r,4] = generation_Very_high_RA_05_co2_tax[t,r]
end

plt_scen = 4 #choose wich demand scenario to plot for 1 (base), 2(low) 3(high) 4(very high)
generation_nuclear_RA = zeros(T)
generation_coal_RA = zeros(T)
generation_gas_RA = zeros(T)
generation_on_wind_RA = zeros(T)
generation_off_wind_fix_RA = zeros(T)
generation_off_wind_float_RA = zeros(T)
generation_solar_RA = zeros(T)
generation_VRES_RA = zeros(T)
for t in 1:T
    generation_nuclear_RA[t] = (generation_all_RA[t,1,plt_scen])/1000
    generation_coal_RA[t] = (generation_all_RA[t,2,plt_scen])/1000
    generation_gas_RA[t] = (generation_all_RA[t,3,plt_scen])/1000
    generation_on_wind_RA[t] = (generation_all_RA[t,4,plt_scen])/1000
    generation_off_wind_fix_RA[t] = (generation_all_RA[t,5,plt_scen])/1000
    generation_off_wind_float_RA[t] = (generation_all_RA[t,6,plt_scen])/1000
    generation_solar_RA[t] = (generation_all_RA[t,7,plt_scen])/1000
end
for t in 1:T
    generation_VRES_RA[t] = generation_on_wind_RA[t] + generation_off_wind_fix_RA[t] + generation_off_wind_float_RA[t] + generation_solar_RA[t]
end

net_load_RA = zeros(T)
for t in 1:T
    net_load_RA[t] = demand_base[t] - generation_VRES_RA[t]
end
dc_net_load_RA = sort(net_load_RA,rev=true)


dc_demand_base_RA = sort(demand_very_high, rev= true) #Choose demand scenario    
dc_generation_nuclear_RA = sort(generation_nuclear_RA, rev=true)
dc_generation_coal_RA = sort(generation_coal_RA, rev=true)
dc_generation_gas_RA = sort(generation_gas_RA, rev=true)
dc_generation_on_wind_RA = sort(generation_on_wind_RA, rev=true)
dc_generation_off_wind_fix_RA = sort(generation_off_wind_fix_RA, rev=true)
dc_generation_off_wind_float_RA = sort(generation_off_wind_float_RA, rev=true)
dc_generation_solar_RA = sort(generation_solar_RA, rev=true)
dc_generation_VRES_RA = sort(generation_VRES_RA, rev=true)

#NEED TO GET CHARGE AND DISCHARGE FROM OTHER DEMAND SCENARIOS

DC_plot_gen_RA_no_policy = Plots.plot([dc_demand_base_RA,dc_generation_VRES_RA,dc_net_load_RA,dc_generation_coal_RA,dc_generation_gas_RA, dc_generation_on_wind_RA,dc_generation_off_wind_fix_RA,dc_generation_solar_RA],tickfontsize =10,title = "Duration Curve",label = ["Demand" "VRES" "Net Load"  "Coal" "Gas" "Onshore Wind" "Offshore Wind Fixed" "Solar"],xlabel = "Time [h]",ylabel = "GWh",legendfont =font(8),legendfontsize=6, color =[Palette_dark[2] Palette_green[7] Palette[8]  Palette_dark[11] Palette_dark[8] Palette[1] Palette[2] palette_clear[6]])
DC_plot_gen_RA_co2_cap = Plots.plot([dc_demand_base_RA,dc_generation_VRES_RA,dc_net_load_RA,dc_generation_nuclear_RA,dc_generation_gas_RA, dc_generation_on_wind_RA,dc_generation_off_wind_fix_RA,dc_generation_off_wind_float_RA,dc_generation_solar_RA],tickfontsize =10,title = "Duration Curve",label = ["Demand" "VRES" "Net Load"  "Nuclear" "Gas" "Onshore Wind" "Offshore Wind Fixed" "Offshre Wind Float" "Solar"],xlabel = "Time [h]",ylabel = "GWh",legendfont =font(8),legendfontsize=6, color =[Palette_dark[2] Palette_green[7] Palette[8]  Palette[3] Palette_dark[8] Palette[1] Palette[2] Palette[9] palette_clear[6]])
DC_plot_gen_RA_co2_tax = Plots.plot([dc_demand_base_RA,dc_generation_VRES_RA,dc_net_load_RA,dc_generation_nuclear_RA,dc_generation_gas_RA, dc_generation_on_wind_RA,dc_generation_off_wind_fix_RA,dc_generation_solar_RA],tickfontsize =10,title = "Duration Curve",label = ["Demand" "VRES" "Net Load"  "Nuclear" "Gas" "Onshore Wind" "Offshore Wind Fixed" "Solar"],xlabel = "Time [h]",ylabel = "GWh",legendfont =font(8),legendfontsize=6, color =[Palette_dark[2] Palette_green[7] Palette[8]  Palette[3] Palette_dark[8] Palette[1] Palette[2] palette_clear[6]])
Plots.savefig(DC_plot_gen_RA_co2_tax,"Results/Results_DC_plots/DC_plot_RN_co2_tax_very_high_demand.pdf") #,dc_generation_nuclear_RA    "Nuclear" Palette[3]

charge_RA_no_policy = CSV.read(string(inputs_path_plots,sep,"All_Zones_charge_all_demand_scen_RA_No_co2_policies.csv"), DataFrame, header=true)
discharge_RA_no_policy = CSV.read(string(inputs_path_plots,sep,"All_Zones_discharge_all_demand_scen_RA_No_co2_policies.csv"), DataFrame, header=true)
charge_RA_co2_cap = CSV.read(string(inputs_path_plots,sep,"All_Zones_charge_all_demand_scen_RA_co2_cap.csv"), DataFrame, header=true)
discharge_RA_co2_cap = CSV.read(string(inputs_path_plots,sep,"All_Zones_discharge_all_demand_scen_RA_co2_cap.csv"), DataFrame, header=true)
charge_RA_co2_tax = CSV.read(string(inputs_path_plots,sep,"All_Zones_charge_all_demand_scen_RA_co2_tax_and_Nuclear_const.csv"), DataFrame, header=true)
discharge_RA_co2_tax = CSV.read(string(inputs_path_plots,sep,"All_Zones_discharge_all_demand_scen_RA_co2_tax_and_Nuclear_const.csv"), DataFrame, header=true)

storage_sum = zeros(T)

for t in 1:T
    storage_sum[t] = -1*charge_RA_co2_tax[t,plt_scen]/1000 + discharge_RA_co2_tax[t,plt_scen]/1000
end
dc_storage = sort(storage_sum,rev=true)

DC_plot_storage_vs_NetLoad = Plots.plot([dc_net_load_RA,dc_storage],tickfontsize =10,title = "Duration Curve",label = ["Net Load" "Storage"], xlabel = "Time [h]",ylabel = "GWh",legendfont =font(8),legendfontsize=6, color = [Palette[8] Palette[6]])
Plots.savefig(DC_plot_storage_vs_NetLoad,"Results/Results_DC_plots/DC_plot_RN_co2_tax_netLoad_vs_storage.pdf")



#PLOT HOW THE RENEWABLES COVER THE DEMAND OVER 5 DAYS
Demand_1_week = demand_very_high[169:336] #Week number 10
week = zeros(168)
for i in 1:168
    week[i] = 168+i
end
#FOR RA WITH CO2 CAP
on_wind_1_week = generation_on_wind_RA[168:336]
off_wind_fix_1_week = generation_off_wind_fix_RA[168:336]
off_wind_float_1_week = generation_off_wind_float_RA[168:336]
solar_1_week = generation_solar_RA[168:336]
VRES_1_week = zeros(168)
layer_1 = zeros(168)
layer_2 = zeros(168)
layer_3 = zeros(168)
layer_4 = zeros(168)
layer_5 = zeros(168)
for t in 1:168
    VRES_1_week[t] = on_wind_1_week[t] + off_wind_fix_1_week[t] + off_wind_float_1_week[t] + solar_1_week[t]
    layer_1[t] = off_wind_float_1_week[t]
    layer_2[t] = off_wind_float_1_week[t] +off_wind_fix_1_week[t]
    layer_3[t] = on_wind_1_week[t] + off_wind_fix_1_week[t] + off_wind_float_1_week[t]
    layer_4[t] = on_wind_1_week[t] + off_wind_fix_1_week[t] + off_wind_float_1_week[t] + solar_1_week[t]
    layer_5[t] = on_wind_1_week[t] + off_wind_fix_1_week[t] + off_wind_float_1_week[t] + solar_1_week[t] + storage_sum[t+168]
end


plt_all_no_policy = Plots.plot(week,[Demand_1_week,layer_2,layer_3,layer_4],tickfontsize =10,fill=true,fillalpha= [0 0.2 0.2 0.2 0.2],title = "Demand and VRES Production",label = ["Demand" "Offshore Wind Fixed" "Onshore Wind" "Solar"], xlabel = "Time [h]",ylabel = "GWh",legendfont =font(8),legendfontsize=6 ,color=[Palette_dark[2] Palette[4] Palette[9] palette_clear[6]])
plt_all_co2_cap = Plots.plot(week,[Demand_1_week,layer_1,layer_2,layer_3,layer_4,storage_sum[169:336],layer_5],tickfontsize =10,fill=true,fillalpha= [0 0.2 0.2 0.2 0.2 0 0],title = "Demand and VRES Production",label = ["Demand"  "Offshore Wind Float" "Offshore Wind Fixed" "Onshore Wind" "Solar" "Storage" "VRES + Storage"], xlabel = "Time [h]",ylabel = "GWh",legendfont =font(8),legendfontsize=5 ,color=[Palette_dark[2] Palette[1] Palette[2] Palette[9] palette_clear[6] :orange :black])
plt_all_co2_tax = Plots.plot(week,[Demand_1_week,layer_2,layer_3,layer_4,storage_sum[169:336],layer_5],tickfontsize =10,fill=true,fillalpha= [0 0.2 0.2 0.2 0 0],title = "Demand and VRES Production",label = ["Demand" "Offshore Wind Fixed" "Onshore Wind" "Solar" "Storage" "VRES + Storage"], xlabel = "Time [h]",ylabel = "GWh",legendfont =font(8),legendfontsize=5 ,color=[Palette_dark[2] Palette[2] Palette[9] palette_clear[6] :orange :black])
Plots.savefig(plt_all_co2_tax,"Results/Results_DC_plots/Demand_and_VRES_plot_RA_co2_tax_very_high_scenario.pdf")


