import numpy as np
import matplotlib.pyplot as plt

# Constants
R = 8.314  # Gas constant in J/(mol*K)
t = np.linspace(0, 30*60, 60*60)  # Time array from 0 to 30 minutes with 1-second intervals

# Initial conditions
interfacial_energy = [475e-3, 450e-3, 350e-3]  # in J/m^2
equilibreum_concentration = [0.011, 0.003, 0.0005]  # in mol/m^3
diffusion_coefficient = [2.6e-16, 4.0e-18, 2.0e-21]  # in m^2/s
T = [780 + 273.15, 660 + 273.15, 500 + 273.15]  # in K
molar_volume = [(7.09e-6) * 3 * 16.35 * (T[i] - 25) for i in range(len(T))]

# Function to calculate r values
def calculate_r(interfacial_energy, molar_volume, equilibreum_concentration, diffusion_coefficient, T, previous_r=0):
    return np.cbrt(
        (8 * interfacial_energy * (molar_volume**2) * equilibreum_concentration * diffusion_coefficient * t) / 
        (9 * R * T) + (previous_r**3)
    )

# Calculating r values for each set of initial conditions
r_values = []
previous_r = 0
for i in range(len(T)):
    r = calculate_r(
        interfacial_energy[i],
        molar_volume[i],
        equilibreum_concentration[i],
        diffusion_coefficient[i],
        T[i],
        previous_r
    )
    r_values.append(r)
    previous_r = r[-1]  # Update previous_r to the last value of the current r

# Plotting the results
for i, r in enumerate(r_values):
    plt.plot(t, r, label=f'Sample {chr(65+i)} Radius')

plt.xlabel('Time (seconds)')
plt.ylabel('r value')
plt.title('Plot of r vs t')
plt.legend()
plt.grid(True)
plt.show()

r1 = r_values[0]
r2 = r_values[1]
r3 = r_values[2]

G = 81600*10**6
b = 0.255*10**-9
f = 0.0584
strength_A = ((0.538*G*b*(np.sqrt(f)))/(2*r1[-1]))*(np.log((r1[-1])/(2*b)))
print(strength_A)

strength_B = ((0.538*G*b*(np.sqrt(f)))/(2*r2[-1]))*(np.log((r2[-1])/(2*b)))
print(strength_B)

strength_C = ((0.538*G*b*(np.sqrt(f)))/(2*r3[-1]))*(np.log((r3[-1])/(2*b)))
print(strength_C)
                     