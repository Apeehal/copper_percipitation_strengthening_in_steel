import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

pi = np.pi
e = np.e
# Constants
R = 8.314  # Gas constant in J/(mol*K)
t = np.linspace(0.001, 30*60, 100000)  # Time array from 0 to 30 minutes with 1-second intervals

# Initial conditions
interfacial_energy = [475e-3, 450e-3, 350e-3]  # in J/m^2
#interfacial_energy = [51*10**-3, 123*10**-3, 314*10**-3]  # in J/m^2
#equilibreum_concentration = [0.011, 0.003, 0.0005]  # in mol/m^3
diffusion_coefficient = [2.6e-16, 4.0e-18, 2.0e-21]  # in m^2/s
T = [780 + 273.15, 660 + 273.15, 500 + 273.15]  # in K
molar_volume = [(7.09e-6) * 3 * 16.35 * (T[i] - 25) for i in range(len(T))]

cu_wt = 1.04/100
fe_density = 7800
mol_cu = 63.546*10**-3
cu_density = 8850
#a = 1000*10**-9
a = 348.79*10**-9

k = 1
dvdt = []
v = []
v[0] =( (0.128*10**-9) ** (3) ) * (4/3)
 
for i in len(1, t):
    dvdt[i] = ((k*cu_wt*fe_density)/(mol_cu))  /  ( ( (  (3*e **(np.cbrt(((4*pi)/v)))) - v*e**( (-v/3) + (np.cbrt((4*pi)/3)) ) )    / ( 4* pi  ) )    + (   (k*cu_density*t)/((mol_cu)*(a**3))  +  (k*cu_density*v)/((mol_cu)*(a**3))  )       + (  ((v[-i]) * (e**(-v/3)) * (e**(np.cbrt((4*pi)/3))) )/(4*pi)   )    ) 
