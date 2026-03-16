import matplotlib.pyplot as plt
import numpy as np

# Data
risk_neutral = [0, 0, 0, 872]
risk_averse = [0, 0, 0, 274]
x_labels = ["Low\ndemand", "Base\ndemand", "High\ndemand", "Very\nhigh\ndemand"]

# Create x positions
x = np.arange(len(x_labels))

width = 0.35  # Width of the bars

# Plot
colors = plt.cm.get_cmap("tab20c")(np.linspace(0,1,20))[0:][::4][0:2]
font_size=18

fig, ax = plt.subplots(figsize=(6.5, 5))
# ax.tick_params(labelsize=font_size)
ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelsize=font_size)
ax.minorticks_on()
ax.spines["top"].set_color("#D4D4D4")      # Top spine to red
ax.spines["right"].set_color("white")   # Right spine to blue
ax.spines["left"].set_color("#D4D4D4")   # Left spine to green
ax.spines["bottom"].set_color("#D4D4D4")# Bottom spine to purple
ax.spines["top"].set_linewidth(0.25)     # Top spine thicker
ax.spines["right"].set_linewidth(0.25)   # Right spine even thicker
ax.spines["left"].set_linewidth(0.25)  # Left spine medium thickness
ax.spines["bottom"].set_linewidth(0.25)  
ax.grid(True, which='major', linestyle='-', linewidth=0.5, color="#D1D1D1") 
ax.grid(True, which='minor', linestyle='-', linewidth=0.15, color="#E2E2E2")
bars1 = ax.bar(x - width/2, risk_neutral, width, label='Risk-neutral', color=colors[0], zorder=1, edgecolor="black", linewidth=1)
bars2 = ax.bar(x + width/2, risk_averse, width, label='Risk-averse', color=colors[1],zorder=2, edgecolor="black", linewidth=1)
ax.set_axisbelow(True)
ax.grid(True)
# Labels and titles
ax.set_ylabel('CO$_2$ price (€/tCO$_2$)', fontsize=font_size)
ax.set_xlabel('Scenarios', fontsize=font_size)
ax.set_xticks(x)
ax.set_xticklabels(x_labels)
ax.legend(fontsize=font_size)


# Show plot
plt.savefig("./co2prices_020425.pdf", dpi=600, bbox_inches='tight')
plt.tight_layout()
plt.show()

