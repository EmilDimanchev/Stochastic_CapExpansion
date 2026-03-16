include("gep_stochastic.jl")
using Plots; theme(:bright)
using RDatasets
using StatsPlots
using Distributions
using StatsBase

#using ."gep_stochastic.jl"

Version_1 = 2
x_time = df_emission.Time
x_resources_all = df_cap.Resource
x_resources = df_rev.Resource
if Version_1 == 1
    y_tech1 = df_cap.Capacity
    y_em_coal1 =df_emission.Coal
    y_em_Gas1 = df_emission.Gas
    y_rev1 = df_rev.Revenue
    y_nse1 =df_nse.Non_served_energy
elseif Version_1 == 2
    y_tech2 = df_cap.Capacity
    y_em_coal2 =df_emission.Coal
    y_em_Gas2 = df_emission.Gas
    y_rev2 = df_rev.Revenue 
    y_nse2 = df_nse.Non_served_energy
elseif Version_1 == 3
    y_tech3 = df_cap.Capacity
    y_em_coal3 =df_emission.Coal
    y_em_Gas3 = df_emission.Gas
    y_rev3 = df_rev.Revenue 
    y_nse3 = df_nse.Non_served_energy
elseif Version_1 == 4
    y_tech4 = df_cap.Capacity
    y_em_coal4 =df_emission.Coal
    y_em_Gas4 = df_emission.Gas
    y_rev4 = df_rev.Revenue 
    y_nse4 = df_nse.Non_served_energy
elseif Version_1 == 5
    y_tech5 = df_cap.Capacity
    y_em_coal5 =df_emission.Coal
    y_em_Gas5 = df_emission.Gas
    y_rev5 = df_rev.Revenue 
    y_nse5 = df_nse.Non_served_energy

end


gamma =[0,0.25,0.5,0.75,1]
ra = 5
Staced_tech = zeros(R_all, ra)
Staced_tech[:,1] = y_tech1
Staced_tech[:,2] = y_tech2
Staced_tech[:,3] = y_tech3
Staced_tech[:,4] = y_tech4
Staced_tech[:,5] = y_tech5
st_plt = transpose(Staced_tech)
Y_label_stack =["","0","0.25","0.5","0.75","1"]
Palette = cgrad(:Paired_12)
#Plot Grouped tech mix
#ONLY FOR SIMULATION WITH CO2 CAP!!
groupbar_techmix = StatsPlots.groupedbar([4;st_plt[1:ra,:]],bar_position =:stack, bar_width=0.7, xlim = (1.5,8), xticks =(1:8,Y_label_stack),title ="Different levels of Risk Aversion",label = ["Nuclear" "Coal" "Gas" "Wind" "Solar" "Storage"], ylabel = "MW",xlabel ="Level of Risk Aversion", color =[Palette[3] Palette[12] Palette[10] Palette[1] Palette[11] Palette[8]]) #:green3 :grey5 :red3 :blue3 :lighttest :darkviolet
Plots.savefig(groupbar_techmix,"Plots/grouped_RA_with_CO2_Cap.pdf")
df_different_level_RA_CO2Cap = DataFrame(st_plt,resources[1:R_all,1])
insertcols!(df_different_level_RA_CO2Cap, 1, :gamma => gamma)
CSV.write(string(results_path,sep,"Different_level_RA_with_CO2_cap"), df_different_level_RA_CO2Cap)


#ONLY FOR SIMULATIONS WITH NO CO2 CAP
groupbar_techmix_NoCap = StatsPlots.groupedbar([4;st_plt[1:ra,:]],bar_position =:stack, bar_width=0.7, xlim = (1.5,8), xticks =(1:8,Y_label_stack),title ="Different levels of Risk Aversion",label = ["Nuclear" "Coal" "Gas" "Wind" "Solar" "Storage"], ylabel = "MW",xlabel ="Level of Risk Aversion", color =[Palette[3] Palette[12] Palette[10] Palette[1] Palette[11] Palette[8]]) #:green3 :grey5 :red3 :blue3 :lighttest :darkviolet
Plots.savefig(groupbar_techmix_NoCap,"Plots/grouped_RA_without_CO2_Cap.pdf")
df_different_level_RA_NO_CO2Cap = DataFrame(st_plt,resources[1:R_all,1])
insertcols!(df_different_level_RA_NO_CO2Cap, 1, :gamma => gamma)
CSV.write(string(results_path,sep,"Different_level_RA_with_No_CO2_cap"), df_different_level_RA_NO_CO2Cap)

#Plot Technology mix
p_tecmix_1 = Plots.plot(x_resources_all,y_tech1,label = "Capacity",title = "Technology mix without CO2 cap",xlabel="Technologies",ylabel="Capacity [MW]",seriestype =[:bar])
p_tecmix_2 = Plots.plot(x_resources_all,y_tech2, label = "Capacity",title = "Technology mix with CO2 cap",xlabel="Technologies",ylabel="Capacity [MW]",seriestype = [:bar])
Combined_plots_tech = Plots.plot(p_tecmix_1,p_tecmix_2,layout =(2,1))
Plots.savefig(Combined_plots_tech,"Plots/D_combined_plots_tech.pdf")
Plots.savefig(p_tecmix_1,"Plots/Tecmix_deterministic_no_CO2_cap.pdf")
Plots.savefig(p_tecmix_2,"Plots/Tecmix_deterministi_with_CO2_cap.pdf")

#Plot emissions
p_em_gas_1 = Plots.plot(x_time,y_em_Gas1,tickfontsize =10,titel = "Without CO2 Cap",label = "Emission Gas",xlabel = "Time",ylabel = "Ton CO2")
p_em_gas_2 = Plots.plot(x_time,y_em_Gas2,ylim=(0,35000),tickfontsize =10,titel = "With CO2 Cap",label = "Emission Gas",xlabel = "Time",ylabel = "Ton CO2")
p_em_coal_1 = Plots.plot(x_time,y_em_coal1,tickfontsize =10,ylim=(35000,70000),titel = "Without CO2 Cap",label = "Emission Coal",xlabel = "Time",ylabel = "Ton CO2")
p_em_coal_2 = Plots.plot(x_time,y_em_coal2,tickfontsize =10,titel = "With CO2 Cap",label = "Emission Coal",xlabel = "Time",ylabel = "Ton CO2")
combined_plots_em =Plots.plot(p_em_gas_1,p_em_gas_2,p_em_coal_1,p_em_coal_2, layout =(2,2),plot_title ="Emissions from Coal and Gas")
Plots.savefig(combined_plots_em,"Plots/em_RA.pdf")


#Plot Revenue per technology
p_rev_1 = Plots.plot(x_resources,y_rev1,label = "Revenue",title = "Without CO2 cap",ylabel="USD",xlabel="Resource",seriestype =[:bar])
p_rev_2 = Plots.plot(x_resources,y_rev2,label = "Revenue",title = "With CO2 cap",ylabel="USD",xlabel="Resource",seriestype =[:bar])
Comb_plot_rev = Plots.plot(p_rev_1,p_rev_2,layout = (2,1))
Plots.savefig(Comb_plot_rev,"Plots/rev_RA.pdf")


#Plot NSE
p_nse_1 = Plots.plot(x_time,y_nse1,label= "NSE without CO2 cp", title ="Non served energy " , xlabel = "Time", ylabel = "MWh"  ) 
p_nse_2 = Plots.plot(x_time,y_nse2,label= "NSE with CO2 cap", title ="Non served energy" , xlabel = "Time", ylabel = "MWh"  )
comb_plots_nse = Plots.plot(p_nse_1,p_nse_2, title = "Non served energy")   
Plots.savefig(comb_plots_nse,"Plots/nse_RN.pdf")
