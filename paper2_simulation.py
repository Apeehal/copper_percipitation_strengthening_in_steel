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


vol_frac = (wt_cu/cu_density)/(100/fe_density)
Ls = r1[-1]*(np.sqrt(((2*pi)/(3*vol_frac))))
M = 3#conversion factor between shear and tensile strength
J = 1 #assume 0.8 for now (can be used as a fitting parameter, meant to be betwee 0.8 and 1)
G = 48300e6
b = 0.255e-9


"""
solution strengthening is also simulated as after 8h, a lot of the strength is likely to be in solution
"""


def solution_strengthening(L):
    mass = ( ((wt_cu)*(Ls**3))/(100/fe_density) ) - ((4*(r1[-1]**3)*pi*cu_density)/(3))
    wt_cu_outside = (mass/(fe_density*(L**3)))*100
    strength_gain = (wt_cu_outside/0.1)*3.8
    return strength_gain




"""
OROWAN MODEL 
- Source: https://www.sciencedirect.com/science/article/pii/S0927025614002572
- Fitted the inputs to match the target value of 78MPa for tensile strength
- Target: 365.8MPa
"""

#shear stress orowan:

nu = 0.25 #can again be used as a fiting parameter (between 0.25 and 0.33)

gain_tensile_strength_orowan = (((G*b)/Ls) * J * M)*10**-6 + solution_strengthening(Ls)
print("Orowan", gain_tensile_strength_orowan) 




"""
ASHBY-OROWAN MODEL 
- Source: https://www.sciencedirect.com/science/article/pii/S0927025614002572
- Fitted the inputs to match the target value of 78MPa for tensile strength
- Target: 365.8MPa
- Assumed Tensile Strength = M * Shear Strength (generally true, but the Taylor Factor for polycrystalline materials is 3)
"""

rs = r1[-1]
ri = 4*b #meant to be between b and 4b 

gain_tensile_strength_Ashby_Orowan = (M*((J*G*b)/(2*pi*np.sqrt(1-nu)*Ls))*np.log(2*rs/ri) ) *10**-6 + solution_strengthening(Ls)
print("Ashby-Orowan", gain_tensile_strength_Ashby_Orowan)






"""
JACKSON-REED MODEL 
- Source: https://www.sciencedirect.com/science/article/pii/S2589152920300995
- Fitted the inputs to match the target value of 78MPa for tensile strength
- Target: 365.8MPa
- Assumed Tensile Strength = M * Shear Strength (generally true, but the Taylor Factor for polycrystalline materials is 3)
"""


gain_tensile_strength_Jackson_Reed = (((M*G*b)/(r1[-1]))*np.sqrt((1.5*vol_frac)) * ((J)/(pi**(3/2))) * np.sqrt(((2*pi*interfacial_energy[0]*r1[-1])/(J*G*b*b))-1))*10**-6 + solution_strengthening(Ls)
print("Jackson-Reed", gain_tensile_strength_Jackson_Reed)



"""
Russel-Brown Model 
- Source: https://www.mdpi.com/2075-4701/10/10/1350
- For Yield Stress
- Target: 186.3MPa
"""

Lx = (1.77*r1[-1])/(np.sqrt(vol_frac))

gain_tensile_strength_Russel_Brown = ((0.8*G*b)/(Lx))*(np.sqrt(1-(0.6**2)))*10**-6 + solution_strengthening(Lx)
print("Russel-Brown", gain_tensile_strength_Russel_Brown)