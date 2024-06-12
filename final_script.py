import numpy as np
import matplotlib.pyplot as plt

interfacial_energy = 475*10**-3
equilibreum_concentration = 0.011
R = 8.314
t = np.linspace(0,30*60,60*60)

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
equilibreum_concentration = 0.003
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
equilibreum_concentration = 0.0005
diffusion_coefficient = 2.0*10**-21
T = 500+273.15
molar_volume = (7.09 * 10**-6) * 3 * 16.35 * (T-25)
r3 = np.cbrt((8*interfacial_energy*(molar_volume**2)*(equilibreum_concentration)*(diffusion_coefficient)*(t))/(9*R*T) + ((r2[-1])**3))
print(r3[-1])

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
                                                    