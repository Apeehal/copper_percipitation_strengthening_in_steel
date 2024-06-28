"""
Attempting to model the change in radius of copper precipitates in IF steel from the following paper:
https://www.sciencedirect.com/science/article/pii/S0167577X06012535

Then comparing and fitting the following strengthening models to compare the gain in strength between sample B&C:
Orowan
Ashby-Orowan
Jackson-Reed
Russel-Brown
"""

"""
Assuming precipitates grow from r = 0.128nm when aged at 550deg C
For now assuming Ostwald Ripening, perhaps we need to model the nucleation & growth stages as well. 
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

R = 8.314  # Gas constant in J/(mol*K)
pi = np.pi
interfacial_energy = [0.52]
diffusion_coefficient = [2e-21]
T = [500 + 273.15]  # in K

#a = 348.79e-9

k_b = 1.380649e-23

wt_cu = 1.18 / 100
fe_density = 7800
mol_mass_cu = 63.546 * 1e-3
cu_density = 8940
mol_vol_cu = mol_mass_cu / cu_density
mol_vol_fe = 7.09e-6




solubility = 1514

def radius(t,r1):
    term1 = 8*interfacial_energy[0]*(mol_vol_cu**2)*diffusion_coefficient[0]*solubility* np.exp(  (2*interfacial_energy[0]* 1.182e-29) /  (r1*k_b*T[0])  )   
    term2 = 9*R*T[0]
    numerator = term1/term2
    
    term3 = 3*(r1**2)
    
    term4 = (8*interfacial_energy[0]*(mol_vol_cu**2)*diffusion_coefficient[0]*t)/(9*R*T[0])
    term5 = solubility
    term6 = (2*interfacial_energy[0]*1.182e-29*np.exp((2*interfacial_energy[0]*1.182e-29)/(k_b*T[0]*r1)) )/(k_b*T[0]*(r1**2))
    
    denominator = term3+term4*term5*term6
    
    return numerator/denominator

# Initial condition
y0 = [0.128e-9]

# Time span (start and end times)
t_span = (0, 8*60*60)

# Time points where solution is to be computed
t_eval = np.linspace(t_span[0], t_span[1], 100000)

# Solve the IVP using solve_ivp with the 'RK45' method
sol1 = solve_ivp(radius, t_span, y0, method='RK45', t_eval=t_eval)

r1 = sol1.y[0]
t1 = np.linspace(0,8*60*60,100000)


"""
No precipitate size given
"""



# Plotting the graph
plt.plot(t1, r1, label='Sample A Radius')
plt.xlabel('Time (seconds)')
plt.ylabel('r value')
plt.title('Plot of r vs t')
plt.legend()
plt.grid(True)
plt.show()