
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# Constants
R = 8.314  # Gas constant in J/(mol*K)
#t = np.linspace(0, 30*60, 100000)  # Time array from 0 to 30 minutes with 1-second intervals

# Initial conditions
interfacial_energy = [475e-3, 450e-3, 350e-3]  # in J/m^2
#equilibreum_concentration = [0.011, 0.003, 0.0005]  # in mol/m^3
diffusion_coefficient = [2.6e-16, 4.0e-18, 2.0e-21]  # in m^2/s
T = [780 + 273.15, 660 + 273.15, 500 + 273.15]  # in K
molar_volume = [(7.09e-6) * 3 * 16.35 * (T[i] - 25) for i in range(len(T))]

cu_wt = 1.4/100
fe_density = 7800
molar_mass_cu = 63.546*10**-3
cu_density = 8850
a = 100*10**-9

r_values = []
#previous_r = 0


total_run_time = 30*60
time_steps = 10000
# Time points for integration
t_eval = np.linspace(0, total_run_time, time_steps)
dt = total_run_time/time_steps
t_span = (0, total_run_time)  # Time span for integration
r0 = [0.0000000000001]  # Initial condition for r(t)



def calculate_r (y,t):
    
    A = (a**3)*(8)*(interfacial_energy[0])*(molar_volume[0]**2)*(diffusion_coefficient[0])*(cu_wt)*(fe_density)
    B = (a**3)*(molar_mass_cu)
    C = (a**3)*(molar_mass_cu)
    D = (8)*(interfacial_energy[0])*(molar_volume[0]**2)*(diffusion_coefficient[0])*(cu_density)*(4/3)*(np.pi)
    

    
    rlag = []
    rlag.append(0.0000000000001)
        
    drdt = (1/3) * ((((A*t) + (B*(rlag[-1]**3))) / (C + D*t)) ** (-2/3)) * (((C + D*t) * A - (A*t + B*(rlag[-1]**3)) * D) / ((C + D*t) ** 2))

    rlag.append(drdt*t)
    
    print(rlag)
    return drdt


sol = solve_ivp(calculate_r, t_span, r0, t_eval=t_eval)

t = sol.t
y = sol.y[0]

# Plot the solution
plt.figure(figsize=(10, 6))
plt.plot(t, y, label='r(t)')
plt.xlabel('Time')
plt.ylabel('r')
plt.title('Solution of Delay Differential Equation using ddeint')
plt.legend()
plt.grid(True)
plt.show()

