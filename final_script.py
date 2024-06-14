import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# Constants
R = 8.314  # Gas constant in J/(mol*K)
#t = np.linspace(0.001, 30*60, 100000)  # Time array from 0 to 30 minutes with 1-second intervals

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
#a = 1000*10**-9
a = 348.79*10**-9

A = (8*interfacial_energy[0]*(molar_volume[0]**2)*diffusion_coefficient[0])/(9*R*T[2])
B = (cu_wt*fe_density)/(molar_mass_cu)
C = (cu_density*(4/3)*(np.pi))/((a**3)*(molar_mass_cu))

r_initial = 0.128 * 10**-9
#r_initial = 1.00782236e-07
t0 = 0.0
step = 0.01
t_next = step

# Define arrays to store results
t_values = [t0]
r_values = [r_initial]

# Define the differential equation function
def ode_function(t, r):
    return (1/3) * (((A * B * t + r_initial**3) / (1 + A * C * t))**(-2/3)) * (((1 + A * C * t) * (A * B) - (A * B * t + r_initial**3) * (A * C)) / ((1 + A * C * t)**2))

# Number of iterations
num_iterations = 10 # Replace with your desired number of iterations

# Iterative solving process
for i in range(num_iterations):
    # Solve the IVP
    solution = solve_ivp(ode_function, (t0, t_next), [r_initial], t_eval=[t_next])
    
    # Extract the results
    t = solution.t[-1]  # Last time point
    r_next = solution.y[0][-1]  # Value of r at the last time point
    
    # Store the results in arrays
    t_values.append(t)
    r_values.append(r_next)
    
    # Print or use the result
    #print(f"Iteration {i+1}: r_next = {r_next}")
    
    # Update initial conditions for the next iteration
    r_initial = r_next
    t0 = t_next
    t_next = t + step  # Update the next time step as desired

# Convert lists to numpy arrays for convenience (optional)
t_values = np.array(t_values)
r_values = np.array(r_values)

# Print or use the final values
#print("Final values:")
#print("t_values:", t_values)
print("r_values:", r_values[-1])

# Plotting the results
plt.figure(figsize=(10, 6))
plt.plot(t_values, r_values, marker='o', linestyle='-', color='b', label='Radius')
plt.title('Radius vs Time')
plt.xlabel('Time')
plt.ylabel('Radius')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()