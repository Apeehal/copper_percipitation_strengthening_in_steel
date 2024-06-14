
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


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

cu_wt = 1.4/100
fe_density = 7800
molar_mass_cu = 63.546*10**-3
cu_density = 8850
a = 10000*10**-9



A = (8*interfacial_energy[0]*(molar_volume[0]**2)*diffusion_coefficient[0])/(9*R*T[0])
B = (cu_wt*fe_density)/(molar_mass_cu)
C = (cu_density*(4/3)*(np.pi))/((a**3)*(molar_mass_cu))





r_initial = 0.128*10**-9
def first_iteration(t, r):
    return (1/3)*(((A*B*t+r_initial**3)/(1+A*C*t))**(-2/3))*(((1+A*C*t)*(A*B)-(A*B*t+r_initial**3)*(A*C))/((1+A*C*t)**2))


"""
t0 = 0
r0 = [0.128*10**-9]
t_end = 30


# Solve the IVP with dense output enabled
solution = solve_ivp(first_iteration, (t0, t_end), r0, method='RK45', dense_output=True)

# Evaluate the solution at the next time step
y_next = solution.sol(t_end)

# Print the results
print(f"Value of y at t = {t_end}: {y_next[0]}")

"""

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

def second_iteration(t,r):
    return (1/3)*(((A*B*t+r_initial**3)/(1+A*C*t))**(-2/3))*(((1+A*C*t)*(A*B)-(A*B*t+r_initial**3)*(A*C))/((1+A*C*t)**2))


# Set initial conditions
t0 = 0.0000000001 # Initial time
r0 = [r_initial]      # Initial value of y
t_next = 0.0000000002  # The next time step we want to find

# Solve the IVP and evaluate only at the initial and next time step
solution1 = solve_ivp(first_iteration, (t0, t_next), r0, t_eval=[t0, t_next])

# Extract the results
t1 = solution1.t
y1 = solution1.y[0]

print(y1[1])

