include("gep_stochastic.jl")
#using ."gep_stochastic.jl"

#ploting the recomended output
display(Plots.plot(df_cap.Resource, df_cap.Capacity,title ="Production capacity",ylabel="Capacity [MWh]" ,seriestype =[:bar], palette = cgrad(:greens), fill=0, alpha=0.6))

#plot revenue/ earnings
display(Plots.plot(df_rev.Resource, df_rev.Revenue, title = "Net income per resource", ylabel = "Revenue [USD]" ,seriestype =[:bar], palette = cgrad(:blues), fill=0, alpha=0.6))

#Plot emmissions
display(Plots.plot(collect(1:T), emissions ,title = "Emission per Technology",ylabel = "[Ton CO2]", seriestype =[:bar], palette = cgrad(:reds), fill=0, alpha=0.2, layout = (3,2))) #, label =("Nuclear", "Gas","Wind","Solar","Batteries")


#Plot non served energy
display(Plots.plot(df_nse.Time, df_nse.Non_served_energy, xlabel = "Time [h]", ylabel="[MWh]", title ="Non Served Energy", label = "NSE"))
