"""
Attempting to model the change in radius of copper precipitates in TRIP steel from the following paper:
https://www.sciencedirect.com/science/article/pii/S0261306914004361

Then comparing and fitting the following strengthening models to compare the gain in strength between sample B&C:
Orowan
Ashby-Orowan
Jackson-Reed
Russel-Brown
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

R = 8.314  # Gas constant in J/(mol*K)
pi = np.pi
interfacial_energy = [0.61, 0.556, 0.52]
diffusion_coefficient = [2.6e-16, 4e-18, 2e-21]
T = [780 + 273.15, 660 + 273.15, 500 + 273.15]  # in K

#a = 348.79e-9

k_b = 1.380649e-23

wt_cu = 1.04 / 100
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
t_span = (0, 30*60)

# Time points where solution is to be computed
t_eval = np.linspace(t_span[0], t_span[1], 100000)

# Solve the IVP using solve_ivp with the 'RK45' method
sol1 = solve_ivp(radius, t_span, y0, method='RK45', t_eval=t_eval)

r1 = sol1.y[0]
t1 = np.linspace(0,30*60,100000)





solubility = 1514

def radius2(t,r2):
    term1 = 8*interfacial_energy[1]*(mol_vol_cu**2)*diffusion_coefficient[1]*solubility* np.exp(  (2*interfacial_energy[1]* 1.182e-29) /  (r2*k_b*T[1])  )   
    term2 = 9*R*T[1]
    numerator = term1/term2
    
    term3 = 3*(r2**2)
    
    term4 = (8*interfacial_energy[1]*(mol_vol_cu**2)*diffusion_coefficient[1]*t)/(9*R*T[1])
    term5 = solubility
    term6 = (2*interfacial_energy[1]*1.182e-29*np.exp((2*interfacial_energy[1]*1.182e-29)/(k_b*T[1]*r2)) )/(k_b*T[1]*(r2**2))
    
    denominator = term3+term4*term5*term6
    
    return numerator/denominator

# Initial condition
y0 = [r1[-1]]

# Time span (start and end times)
t_span = (0, 30*60)

# Time points where solution is to be computed
t_eval = np.linspace(t_span[0], t_span[1], 100000)

# Solve the IVP using solve_ivp with the 'RK45' method
sol3 = solve_ivp(radius2, t_span, y0, method='RK45', t_eval=t_eval)

r2 = sol3.y[0]
t2 = np.linspace(30*60,60*60,100000)





solubility = 1514

def radius3(t,r3):
    term1 = 8*interfacial_energy[2]*(mol_vol_cu**2)*diffusion_coefficient[2]*solubility* np.exp(  (2*interfacial_energy[2]* 1.182e-29) /  (r3*k_b*T[2])  )   
    term2 = 9*R*T[2]
    numerator = term1/term2
    
    term3 = 3*(r3**2)
    
    term4 = (8*interfacial_energy[2]*(mol_vol_cu**2)*diffusion_coefficient[2]*t)/(9*R*T[2])
    term5 = solubility
    term6 = (2*interfacial_energy[2]*1.182e-29*np.exp((2*interfacial_energy[2]*1.182e-29)/(k_b*T[2]*r3)) )/(k_b*T[2]*(r3**2))
    
    denominator = term3+term4*term5*term6
    
    return numerator/denominator

# Initial condition
y0 = [r2[-1]]

# Time span (start and end times)
t_span = (0, 30*60)

# Time points where solution is to be computed
t_eval = np.linspace(t_span[0], t_span[1], 100000)

# Solve the IVP using solve_ivp with the 'RK45' method
sol3 = solve_ivp(radius3, t_span, y0, method='RK45', t_eval=t_eval)

r3 = sol3.y[0]
t3 = np.linspace(60*60,90*60,100000)


plt.figure(figsize=(10, 6))
plt.plot(t1, r1 * 1e9, label='T = {} K'.format(T[0]))
plt.plot(t2, r2 * 1e9, label='T = {} K'.format(T[1]))
plt.plot(t3, r3 * 1e9, label='T = {} K'.format(T[2]))
plt.xlabel('Time (s)')
plt.ylabel('Precipitate Radius (nm)')
plt.title('Particle Radius Evolution over Time')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

"""
OROWAN MODEL 
- Source: https://www.sciencedirect.com/science/article/pii/S0927025614002572
- Fitted the inputs to match the target value of 78MPa for tensile strength
- Target: 78MPa
"""
#vol_frac = 9.11e-3
vol_frac = 0.012293
#shear stress orowan:
J = 0.8 #assume 0.8 for now (can be used as a fitting parameter, meant to be betwee 0.8 and 1)
G = 48300e6
b = 0.255e-9
nu = 0.25 #can again be used as a fiting parameter (between 0.25 and 0.33)
Ls = r3[-1]*(np.sqrt(((2*pi)/(3*vol_frac))))

gain_tensile_strength_orowan = (((G*b)/Ls) * J * 2)*10**-6
print("Orowan", gain_tensile_strength_orowan) 



"""
ASHBY-OROWAN MODEL 
- Source: https://www.sciencedirect.com/science/article/pii/S0927025614002572
- Fitted the inputs to match the target value of 78MPa for tensile strength
- Target: 78MPa
- Assumed Tensile Strength = 2 * Shear Strength (generally true, but the Taylor Factor for polycrystalline materials is 3)
"""

#vol_frac = 9.11e-3
vol_frac = 0.012293
#shear stress orowan:
J = 0.8 #assume 0.8 for now (can be used as a fitting parameter, meant to be betwee 0.8 and 1)
G = 48300e6
b = 0.255e-9
nu = 0.25 #can again be used as a fiting parameter (between 0.25 and 0.33)
Ls = r3[-1]*(np.sqrt(((2*pi)/(3*vol_frac))))
rs = r3[-1]
ri = 4*b #meant to be between b and 4b 

gain_tensile_strength_Ashby_Orowan = (2*((J*G*b)/(2*pi*np.sqrt(1-nu)*Ls))*np.log(2*rs/ri) ) *10**-6
print("Ashby-Orowan", gain_tensile_strength_Ashby_Orowan)







"""
JACKSON-REED MODEL 
- Source: https://www.sciencedirect.com/science/article/pii/S2589152920300995
- Fitted the inputs to match the target value of 78MPa for tensile strength
- Target: 78MPa
- Assumed Tensile Strength = 2 * Shear Strength (generally true, but the Taylor Factor for polycrystalline materials is 3)
"""

#vol_frac = 9.11e-3
vol_frac = 0.012293
#shear stress orowan:
J = 0.8 #assume 0.8 for now (can be used as a fitting parameter, meant to be betwee 0.8 and 1)
G = 48300e6
b = 0.255e-9
M = 2
gain_tensile_strength_Jackson_Reed = (((M*G*b)/(r3[-1]))*np.sqrt((1.5*vol_frac)) * ((J)/(pi**(3/2))) * np.sqrt(((2*pi*interfacial_energy[2]*r3[-1])/(J*G*b*b))-1))*10**-6
print("Jackson-Reed", gain_tensile_strength_Jackson_Reed)



"""
Russel-Brown Model 
- Source: https://www.mdpi.com/2075-4701/10/10/1350
- For Yield Stress
- Target: 15MPa
"""

vol_frac = 0.012293
G = 48300e6
b = 0.255e-9
L = (1.77*r3[-1])/(np.sqrt(vol_frac))

gain_tensile_strength_Russel_Brown = ((0.8*G*b)/(L))*(np.sqrt(1-(0.6**2)))*10**-6
print("Russel-Brown", gain_tensile_strength_Russel_Brown)
import numpy as np
import matplotlib.pyplot as plt

interfacial_energy = 475*10**-3
equilibreum_concentration = 0.02
R = 8.314
t = np.linspace(0,30*60,30*60)

#interfacial_energy = [51*10^-3, 123*10^-3, 314*10^-3]

diffusion_coefficient = 2.6*10**-16
T = 780+273.15
molar_volume = (7.09 * 10**-6) * 3 * 16.35 * (T-25)
r1 = np.cbrt((8*interfacial_energy*(molar_volume**2)*(equilibreum_concentration)*(diffusion_coefficient)*(t))/(9*R*T))


# Plotting the graph
plt.plot(t, r1, label='Sample A Radius')
plt.xlabel('Time (seconds)')
plt.ylabel('r value')
plt.title('Plot of r vs t')
plt.legend()
plt.grid(True)
plt.show()

interfacial_energy = 450*10**-3
diffusion_coefficient = 4.0*10**-18
T = 660+273.15
molar_volume = (7.09 * 10**-6) * 3 * 16.35 * (T-25)
r2 = np.cbrt((8*interfacial_energy*(molar_volume**2)*(equilibreum_concentration)*(diffusion_coefficient)*(t))/(9*R*T) + ((r1[-1])**3))


# Plotting the graph
plt.plot(t, r2, label='Sample B Radius')
plt.xlabel('Time (seconds)')
plt.ylabel('r value')
plt.title('Plot of r vs t')
plt.legend()
plt.grid(True)
plt.show()


interfacial_energy = 350*10**-3
diffusion_coefficient = 2.0*10**-21
T = 500+273.15
molar_volume = (7.09 * 10**-6) * 3 * 16.35 * (T-25)
r3 = np.cbrt((8*interfacial_energy*(molar_volume**2)*(equilibreum_concentration)*(diffusion_coefficient)*(t))/(9*R*T) + ((r2[-1])**3))


# Plotting the graph
plt.plot(t, r3, label='Sample V Radius')
plt.xlabel('Time (seconds)')
plt.ylabel('r value')
plt.title('Plot of r vs t')
plt.legend()
plt.grid(True)
plt.show()

G = 81600*10**6
b = 0.255*10**-9
f = 0.0584
strength_A = ((0.538*G*b*(np.sqrt(f)))/(2*r1[-1]))*(np.log((r1[-1])/(2*b)))
print(strength_A)

strength_B = ((0.538*G*b*(np.sqrt(f)))/(2*r2[-1]))*(np.log((r2[-1])/(2*b)))
print(strength_B)

strength_C = ((0.538*G*b*(np.sqrt(f)))/(2*r3[-1]))*(np.log((r3[-1])/(2*b)))
print(strength_C)
                             
