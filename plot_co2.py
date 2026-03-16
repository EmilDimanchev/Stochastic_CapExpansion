# %%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
# %%
# Compile data

pol_cases = ["nopol","cap","tax"]
pol_labels = ["No policy","CO2 cap","CO2 tax"]
risk_cases = ["RN","RA"]
scenarios = ["Base demand","Low demand","High demand","Very high demand","Expected"]

data = pd.DataFrame(0, index=pd.MultiIndex.from_product([risk_cases,scenarios], names=["Risk", "Scenario"]), columns=pol_cases)

for risk in risk_cases:
    for pol in pol_cases:
        df_input = pd.read_csv("./Data/Emission_sum_all_scen_"+risk+"_"+pol+".csv",index_col=0)
        all_co2 = list(df_input["Total"])
        exp_co2 = sum(0.25*i for i in all_co2)
        all_co2.append(exp_co2)
        for scenario in scenarios:
            data.loc[(risk,scenario),pol] = all_co2[scenarios.index(scenario)]

data
# %%
# Plot

fig, axes = plt.subplots(1,2, figsize=(14,6), gridspec_kw={'width_ratios': [1, 3]})

font_size = 22
markers = ["*", "D", "s", "o", ">"]
colors = plt.cm.get_cmap("tab20c")(np.linspace(0,1,20))[0:][::4]
z_order = [5,1,2,3,4]

# First plot - no policy
axes[0].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelsize=font_size)
ind = np.arange(1)
to_plot = data["nopol"]
for scenario in scenarios:
    axes[0].plot(np.arange(2)[0], to_plot.loc["RN"].loc[scenario]/1e6, linestyle="", marker=markers[scenarios.index(scenario)], markersize=18, color=colors[scenarios.index(scenario)], zorder=z_order[scenarios.index(scenario)])
    axes[0].plot(np.arange(2)[1], to_plot.loc["RA"].loc[scenario]/1e6, linestyle="", marker=markers[scenarios.index(scenario)], markersize=18, color=colors[scenarios.index(scenario)], zorder=z_order[scenarios.index(scenario)])

axes[0].minorticks_on()
axes[0].spines["top"].set_color("#D4D4D4")      # Top spine to red
axes[0].spines["right"].set_color("white")   # Right spine to blue
axes[0].spines["left"].set_color("#D4D4D4")   # Left spine to green
axes[0].spines["bottom"].set_color("#D4D4D4")# Bottom spine to purple
axes[0].spines["top"].set_linewidth(0.25)     # Top spine thicker
axes[0].spines["right"].set_linewidth(0.25)   # Right spine even thicker
axes[0].spines["left"].set_linewidth(0.25)  # Left spine medium thickness
axes[0].spines["bottom"].set_linewidth(0.25)

axes[0].set_xlim(-0.3,1.3)
axes[0].set_xticks(np.arange(2))
axes[0].set_xticklabels([r"$\gamma=0$",r"$\gamma=0.5$"], fontsize=font_size)
axes[0].set_xlabel("No policy", fontsize=font_size)
axes[0].set_axisbelow(True)
axes[0].grid(True, which='major', linestyle='-', linewidth=0.5, color="#D1D1D1") 
axes[0].grid(True, which='minor', linestyle='-', linewidth=0.15, color="#E2E2E2")

# Second plot
axes[1].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelsize=font_size)
axes[1].minorticks_on()
axes[1].spines["top"].set_color("#D4D4D4")      # Top spine to red
axes[1].spines["right"].set_color("white")   # Right spine to blue
axes[1].spines["left"].set_color("#D4D4D4")   # Left spine to green
axes[1].spines["bottom"].set_color("#D4D4D4")# Bottom spine to purple
axes[1].spines["top"].set_linewidth(0.25)     # Top spine thicker
axes[1].spines["right"].set_linewidth(0.25)   # Right spine even thicker
axes[1].spines["left"].set_linewidth(0.25)  # Left spine medium thickness
axes[1].spines["bottom"].set_linewidth(0.25)





ind = np.arange(4)
to_plot = data[["cap","tax"]]
for scenario in scenarios:
    axes[1].plot(ind[0], to_plot.loc["RN","tax"].loc[scenario]/1e6, linestyle="", marker=markers[scenarios.index(scenario)], markersize=20, label=scenario, color=colors[scenarios.index(scenario)], zorder=z_order[scenarios.index(scenario)])
    axes[1].plot(ind[2], to_plot.loc["RN","cap"].loc[scenario]/1e6, linestyle="", marker=markers[scenarios.index(scenario)], markersize=20, color=colors[scenarios.index(scenario)], zorder=z_order[scenarios.index(scenario)])
    axes[1].plot(ind[1], to_plot.loc["RA","tax"].loc[scenario]/1e6, linestyle="", marker=markers[scenarios.index(scenario)], markersize=20, color=colors[scenarios.index(scenario)], zorder=z_order[scenarios.index(scenario)])
    axes[1].plot(ind[3], to_plot.loc["RA", "cap"].loc[scenario]/1e6, linestyle="", marker=markers[scenarios.index(scenario)], markersize=20, color=colors[scenarios.index(scenario)], zorder=z_order[scenarios.index(scenario)])

axes[1].set_axisbelow(True)
axes[1].set_xlim(-0.3,3.3)
axes[0].set_ylabel("Emissions (MtCO2)", fontsize=font_size)
axes[1].legend(bbox_to_anchor=(0.95, -0.15), ncol = 3, fontsize=font_size)
axes[1].set_xticks(np.arange(4))
axes[1].set_xticklabels([r"$\gamma=0$",r"$\gamma=0.5$",r"$\gamma=0$",r"$\gamma=0.5$"], fontsize=font_size)
axes[1].set_xlabel("Policy cases", fontsize=font_size)
axes[1].set_ylim(0,190)
axes[1].text(2.1,174,"CO\u2082 cap",fontsize=font_size)
axes[1].text(0.1,164,"Base case:\nFixed CO\u2082 Price",fontsize=font_size)
axes[1].grid(True, which='major', linestyle='-', linewidth=0.5, color="#D1D1D1") 
axes[1].grid(True, which='minor', linestyle='-', linewidth=0.15, color="#E2E2E2")

plt.savefig("./emissions_020425.pdf", dpi=600, bbox_inches='tight')