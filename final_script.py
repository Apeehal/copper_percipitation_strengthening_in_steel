import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


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
molar_mass_cu = 63.546*10**-3
cu_density = 8850
#a = 1000*10**-9
a = 348.79*10**-9

A = (8*interfacial_energy[0]*(molar_volume[0]**2)*diffusion_coefficient[0])/(9*R*T[0])
B = (cu_wt*fe_density)/(molar_mass_cu)
C = (cu_density*(4/3)*(np.pi))/((a**3)*(molar_mass_cu))


num_iterations = 30*60  # Number of iterations
# Initial conditions
dt = (30*60) / 100000
r = np.zeros(num_iterations)
r[0] = 0.128 * 10**-9


for i in range(1, num_iterations):
    r[i] = ((A*B*dt + r[0]**3)/(1 + A*C*dt ))**(1/3)

time_steps = np.arange(num_iterations) * dt  # Time steps corresponding to each iteration

plt.figure(figsize=(10, 6))
plt.plot(time_steps, r, label='r(t)')
plt.title('Plot of r(t) over time')
plt.xlabel('Time')
plt.ylabel('r(t)')
plt.grid(True)
plt.legend()
plt.show()