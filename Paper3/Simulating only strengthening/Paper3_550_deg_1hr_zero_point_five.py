# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 12:52:36 2024

@author: appee
"""


"""
Attempting to model the change in radius of copper precipitates in 42CrMo4 Quench and Tempering Steel from the following paper:
https://onlinelibrary.wiley.com/doi/full/10.1002/srin.202200623

Then comparing and fitting the following strengthening models to compare the gain in strength:
Orowan
Ashby-Orowan
Jackson-Reed
Russel-Brown
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

r1 = [10e-9]
h = 1
R = 8.314  
pi = np.pi
A_0 = 6.022e23
k_b = 1.380649e-23
wt_cu = 0.5
fe_density = 7800
mol_mass_cu = 63.546 * 1e-3
mol_mass_fe = 55.85e-3
cu_density = 8940
mol_vol_cu = mol_mass_cu / cu_density
mol_vol_fe = 7.09e-6

vol_fe = 100/fe_density
vol_cu = (wt_cu/100)/cu_density
#vol_frac = vol_cu/vol_fe

#vol_frac = 0.012293
vol_frac = (wt_cu/cu_density)/((wt_cu/cu_density)+((100-wt_cu)/fe_density))
print("vol_frac",vol_frac)



#interfacial_energy = [0.43]
#diffusion_coefficient = [2e-21]
T = [550 + 273.15]  # in K
interfacial_energy = [((3/5)*T[0]-3.8899)*(10**-3)]


B = 0.034343885822
A = 34242.66544
#diffusion_coefficient = [B*np.exp((-A)/(T[0]))]
diffusion_coefficient = [2e-21]
print("diffusion", diffusion_coefficient)


solubility = []
for i in T:
    x = 10 ** ((6111850/((i)**2)) - ((16478.2/i)) + 10.3242)
    solubility.append(x)
print("solubility", solubility)

initial_size = (2*interfacial_energy[0])/((R*T[0])/(mol_vol_cu))
print("initial size", initial_size)

nu = 0.25 #can again be used as a fiting parameter (between 0.25 and 0.33)
vol_frac = (wt_cu/cu_density)/(100/fe_density)
Ls = r1[-1]*(np.sqrt(((2*pi)/(3*vol_frac))))
M = 3 #conversion factor between shear and tensile strength
J = 0.8 #assume 0.8 for now (can be used as a fitting parameter, meant to be betwee 0.8 and 1)
G = 48300e6
b = 0.255e-9



"""
OROWAN MODEL 
- Source: https://www.sciencedirect.com/science/article/pii/S0927025614002572
- Fitted the inputs to match the target value of 78MPa for tensile strength
- Target: 44MPa
"""

#shear stress orowan:


gain_tensile_strength_orowan = (((G*b)/Ls) * J * M)*10**-6
print("Orowan", gain_tensile_strength_orowan) 




"""
ASHBY-OROWAN MODEL 
- Source: https://www.sciencedirect.com/science/article/pii/S0927025614002572
- Fitted the inputs to match the target value of 78MPa for tensile strength
- Target: 44MPa
- Assumed Tensile Strength = M * Shear Strength (generally true, but the Taylor Factor for polycrystalline materials is 3)
"""

rs = r1[-1]
ri = 2*b #meant to be between b and 4b 

gain_tensile_strength_Ashby_Orowan = (M*((J*G*b)/(2*pi*np.sqrt(1-nu)*Ls))*np.log(2*rs/ri) ) *10**-6
print("Ashby-Orowan", gain_tensile_strength_Ashby_Orowan)






"""
JACKSON-REED MODEL 
- Source: https://www.sciencedirect.com/science/article/pii/S2589152920300995
- Fitted the inputs to match the target value of 78MPa for tensile strength
- Target: 44MPa
- Assumed Tensile Strength = M * Shear Strength (generally true, but the Taylor Factor for polycrystalline materials is 3)


"""
M = 2
gain_tensile_strength_Jackson_Reed = (((M*G*b)/(r1[-1]))*np.sqrt((1.5*vol_frac)) * ((J)/(pi**(3/2))) * np.sqrt(((2*pi*interfacial_energy[0]*r1[-1])/(J*G*b*b))-1))*10**-6
print("Jackson-Reed", gain_tensile_strength_Jackson_Reed)

"""
Russel-Brown Model 
- Source: https://www.mdpi.com/2075-4701/10/10/1350
- For Yield Stress
- Target: 26MPa
"""

Lx = (1.77*r1[-1])/(np.sqrt(vol_frac))

e1_e2 = 0.6 * (np.log10(r1[-1]/Lx)/np.log10(Lx/ri)) + ( np.log10(Lx/r1[-1])  / np.log10(Lx/ri)  )



gain_yield_strength_Russel_Brown = ((J*G*b)/(Lx))*((1-(e1_e2**2))**(1/2))*10**-6
print("Russel-Brown", gain_yield_strength_Russel_Brown)



#plotting bar chart
categories = ['Orowan', 'Ashby-Orowan', 'Jackson-Reed', 'Russel-Brown',]
values1 = [gain_tensile_strength_orowan,gain_tensile_strength_Ashby_Orowan, gain_tensile_strength_Jackson_Reed, gain_yield_strength_Russel_Brown]
values2 = [44, 44, 44, 26]

# Number of categories
n = len(categories)

# Positions of the bars on the x-axis
r1 = np.arange(n)
r2 = [x + 0.25 for x in r1]

# Create bar chart
plt.bar(r1, values1, color='blue', width=0.25, edgecolor='grey', label='Model')
plt.bar(r2, values2, color='green', width=0.25, edgecolor='grey', label='Experiment')

# Add labels
plt.xlabel('Categories')
plt.ylabel('Values')
plt.title('Side-by-Side Bar Chart')

# Add xticks on the middle of the bars
plt.xticks([r + 0.125 for r in range(n)], categories)

# Add legend
plt.legend()

# Show the plot
plt.show()

