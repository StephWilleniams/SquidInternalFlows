"""
Sobol Sensitivity as implemented in SAlib - Squid fluids model
@author: Steve
"""

## Import the needed libraries and functions
# Set the required public imports
import numpy as np
import matplotlib.pyplot as plt
from SALib import ProblemSpec

bounds = [[-1,1],[-1,1],[-1,1],[-1,1]]

sp = ProblemSpec({
        "names": ["x1", "x2", "x3", "x4"],
        "groups": None,
        "bounds": [[-1,1]] * 4,
        "outputs": ["Y"],
    })

params_list = [r'$U_0$',r'$\omega$',r'$U_c$',r'$Y_0$']

# ## Solve the system across the sample set to get the values
param_values = np.loadtxt("sobol_try_4/params_2pow13.txt")
QOI = np.loadtxt("outputs/outputs1.txt")
QOI = np.loadtxt("outputs/outputs2.txt")

# Provide the results to the interface
sp.set_samples(param_values)
sp.set_results(QOI/50)
# sp.set_results(QOI)

## Perform the sobol index calculations
sp.analyze_sobol(calc_second_order=True, print_to_console=True)

bar_width = 0.4
COLOR1 = [0, 0.4470, 0.7410]
COLOR2 = [0.8500, 0.3250, 0.0980]

# Create the plot
fig, ax = plt.subplots()

# Plot bars with error bars
x = np.arange(len(params_list))  # Positions for bars
ax.bar(x - bar_width/2, sp.analysis['S1'], bar_width, label='First Order', color=COLOR1)
ax.bar(x + bar_width/2, sp.analysis['ST'], bar_width, label='Total Order', color=COLOR2)

# Add error bars
yerr_S1 = np.array(sp.analysis['S1_conf']).T  # Transpose for easier indexing
yerr_ST = np.array(sp.analysis['ST_conf']).T

ax.errorbar(x - bar_width/2, sp.analysis['S1'], yerr=yerr_S1, fmt='none', ecolor='black', capsize=7)
ax.errorbar(x + bar_width/2, sp.analysis['ST'], yerr=yerr_ST, fmt='none', ecolor='black', capsize=7)

ax.axhline(0.05, color='y', linestyle='--') 

# Set labels and title
plt.ylim(0, 0.9)
plt.xticks(x, params_list)
plt.xlabel('Parameters', fontsize=16)
plt.ylabel('Sobol Index', fontsize=16)
plt.title('Sobol Indices with Confidence Intervals', fontsize=18)

# Show the plot
plt.tight_layout()
plt.show()
