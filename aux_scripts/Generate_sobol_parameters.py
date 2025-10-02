""": 
Sobol Sensitivity as implemented in SAlib - Squid fluids model
@author: Steve
"""

## Import the needed libraries and functions
import numpy as np
import matplotlib.pyplot as plt
from SALib import ProblemSpec
from SALib.analyze import sobol
from SALib.sample.sobol import sample

# Define the Sobol problem
sp = ProblemSpec({
        "names": ["x1", "x2", "x3", "x4"],
        "groups": None,
        "bounds": [[-1,1]] * 4,
        "outputs": ["Y"],
    })
params_list = [r'$U_0$',r'$\omega$',r'$U_c$',r'$Y_0$']
Nparams = 2**13 # Number of samples to take
param_values = sample(sp, N=Nparams, calc_second_order=True) # Get parameter sets
np.savetxt("inputs/params.txt", param_values) # Save the parameter values to a file
