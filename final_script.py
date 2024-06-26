
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





