
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


pi = np.pi
e = np.exp(1)
R = 8.314e+18  # Gas constant in J/(mol*K)

interfacial_energy = [0.61, 0.556, 0.52 ]
diffusion_coefficient = [260, 4, 0.002]
T = [780 + 273.15, 660 + 273.15, 500 + 273.15]  # in K

a = 348.79

k_b = 1.380649e-5 

wt_cu = 1.04 / 100
fe_density = 7.87e-24
mol_cu = 63.546 * 1e-3
cu_density = 8.94e-24
mol_vol_cu = mol_cu/cu_density
mol_vol_cu = [mol_vol_cu, mol_vol_cu, mol_vol_cu]
mol_vol_fe = [7.09e21, 7.09e21,7.09e21]


r_initial = 0.128*10**-9
def first_iteration(t,r):
    alpha = (2*interfacial_energy[0]*mol_vol_fe[0])/(k_b*T[0])
    beta = (8*interfacial_energy[0]*((mol_vol_cu[0])**2)*diffusion_coefficient[0])/(9*R*T[0])
    term1 = (beta*wt_cu*fe_density)/(mol_cu)
    term2 = (beta*cu_density*4*pi*(r**3))/(3*(a**3)*mol_cu)
        
    term3 = np.exp(alpha/r)  *  (3*r-alpha)

    term4 = np.exp(alpha/r)*(r**3)*( (alpha/(r**2)) - 1 )
            
    term5 = (beta*cu_density*4*pi*(r**2)*t)/((a**3)*mol_cu)
    
    drdt1 =  (term1 + term2) / (term3 + term4 + term5) 

    return drdt1

# Set initial conditions
t0 = 0        # Initial time
r0 = [0.128*10**-9]      # Initial value of y
t_next = 0.0000000001  # The next time step we want to find

# Solve the IVP and evaluate only at the initial and next time step
solution = solve_ivp(first_iteration, (t0, t_next), r0, t_eval=[t0, t_next])

# Extract the results
t = solution.t
y = solution.y[0]

#print(y[1])

r_initial = y[1]
print(r_initial)

