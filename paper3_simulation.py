"""
Attempting to model the change in radius of copper precipitates in 42CrMo4 Quench and Tempering Steel from the following paper:
https://onlinelibrary.wiley.com/doi/full/10.1002/srin.202200623

Then comparing and fitting the following strengthening models to compare the gain in strength:
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


"""
0.57% Cu @ 450 deg C
"""

R = 8.314  # Gas constant in J/(mol*K)
pi = np.pi
k_b = 1.380649e-23
wt_cu = 0.57
fe_density = 7800
mol_mass_cu = 63.546 * 1e-3
cu_density = 8940
mol_vol_cu = mol_mass_cu / cu_density
mol_vol_fe = 7.09e-6
h = 0.5
#a = 348.79e-9


interfacial_energy = [0.39]
diffusion_coefficient = [2e-21]
T = [450 + 273.15]  # in K
solubility_wt = 0.1
solubility = (solubility_wt/mol_mass_cu)/(100/fe_density)


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
t_span = (0, h*60*60)

# Time points where solution is to be computed
t_eval = np.linspace(t_span[0], t_span[1], 100000)

# Solve the IVP using solve_ivp with the 'RK45' method
sol1 = solve_ivp(radius, t_span, y0, method='RK45', t_eval=t_eval)

r1 = sol1.y[0]
t1 = np.linspace(0, h*60*60,100000)


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
- Target: 47MPa
"""

#shear stress orowan:


gain_tensile_strength_orowan = (((G*b)/Ls) * J * M)*10**-6
print("Orowan", gain_tensile_strength_orowan) 




"""
ASHBY-OROWAN MODEL 
- Source: https://www.sciencedirect.com/science/article/pii/S0927025614002572
- Fitted the inputs to match the target value of 78MPa for tensile strength
- Target: 47MPa
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
- Target: 365.8MPa
- Assumed Tensile Strength = M * Shear Strength (generally true, but the Taylor Factor for polycrystalline materials is 3)



gain_tensile_strength_Jackson_Reed = (((M*G*b)/(r1[-1]))*np.sqrt((1.5*vol_frac)) * ((J)/(pi**(3/2))) * np.sqrt(((2*pi*interfacial_energy[0]*r1[-1])/(J*G*b*b))-1))*10**-6
print("Jackson-Reed", gain_tensile_strength_Jackson_Reed)




Russel-Brown Model 
- Source: https://www.mdpi.com/2075-4701/10/10/1350
- For Yield Stress
- Target: 16MPa
"""

Lx = (1.77*r1[-1])/(np.sqrt(vol_frac))

e1_e2 = 0.6 * (np.log10(r1[-1]/Lx)/np.log10(Lx/ri)) + ( np.log10(Lx/r1[-1])  / np.log10(Lx/ri)  )



gain_tensile_strength_Russel_Brown = ((J*G*b)/(Lx))*((1-(e1_e2**2))**(3/4))*10**-6
print("Russel-Brown", gain_tensile_strength_Russel_Brown)
