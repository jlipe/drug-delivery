import numpy as np
from scipy.integrate import *
from matplotlib import pyplot as plt


# define digestion equation - follows model by Li & McClements 2010
# y - concentration of FA (mmol/L), t - digestion time (min)
# p - (c_oil - concentration of ingested TG (mmol/L),
    #D0 - Diameter of oil droplets (cm) - 386nm)
def intestine_digestion(c_FA, t, D0, c_oil):
    # LIPID INPUT PARAMETERS
    Kdig = 3.116e-7              # digestion rate constant (mmol/min*cm^2)
    Vbulk = 15                   # total volume of solution (mL)
    MWoil = 885.4                # molecular weight of ingested TG (g/mol)
    dens_oil = 0.92              # density of ingested oil (g/mL)
    # CORRELATIONS
    m_oil = c_oil*MWoil*Vbulk/1e6      # starting mass of TG (g)
    V_oil = m_oil/dens_oil             # starting volume of oil (mL)
    FA0 = 3*c_oil                      # initial digestible fatty acids, mmol/L
    N = V_oil/(D0**3*3.14/6)           # number of oil droplets

    # PROCESS
    dc_FA = (Kdig/Vbulk)*N*3.14*(D0**2)*((FA0-c_FA)/FA0)**(2/3)*1000

    return dc_FA #mmol/(l*min)

blocksize = 100 # Number of Spaces within each variable

results_array = [] # Initiate an array for the loop results
"""
[c_oil, D0, final_c_FA]
"""


c_oil_array = np.linspace(1, 100, blocksize)  # Sensitivity of c_oil
D0_array = np.linspace(5e-5, 5e-4, blocksize)  # Sensitivity of D0


times = np.linspace(0, 180, 180) #minutes

for c_oil in c_oil_array:
    for D0 in D0_array:
        c_FA = odeint(intestine_digestion, 0, times, (D0, c_oil))
        final_c_FA = float(c_FA[-1])
        results_array.append(final_c_FA)



#         ResultArray[n]= (i, j, np.nanmax(FA))  # Add the c_oil, D0, and max FA into the array
#         n += 1  # Count up

# FAArray = ResultArray[:, 2]  # Splits array into just the max FA
# FAArray = np.reshape(FAArray, (blocksize, blocksize))  # Reshapes from a vector to a square matrix

X, Y = np.meshgrid(c_oil_array, D0_array)
results_array = np.reshape(results_array, (blocksize, blocksize))
print(results_array)
plt.contourf(X, Y, results_array)


plt.title("Sensitivity Analysis")
plt.xlabel("c_oil (mM)")
plt.ylabel("D0 (cm)")
cb = plt.colorbar(plt.contourf(X, Y, results_array))
cb.set_label("Fatty Acid Concentration (mM)")
plt.show()
