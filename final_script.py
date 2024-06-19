import numpy as np
import matplotlib.pyplot as plt


# Time parameters
t = np.linspace(0, 30, 100000)
dt = t[1] - t[0]  # Timestep based on the linspace definition
    
# Initial conditions
vA0 = (0.128 * 10) ** 3 * (4 / 3)
    #v0 = 9.6022845
    #vv0 = 
v1 = np.zeros_like(t)
v1[0] = vA0


# Iterate to solve for v and dvdt
dvdt1 = np.zeros_like(t)

for i in range(1, len(t)):
    
    pi = np.pi
    e = np.exp(1)
    R = 8.314e+18  # Gas constant in J/(mol*K)

    interfacial_energy = [0.51*10**-19, 1.23*10**-19, 3.14*10**-19]  # in J/m^2
    diffusion_coefficient = [2.6e+11, 4e+9, 2e+6]  # in m^2/s
    T = [780 + 273.15, 660 + 273.15, 500 + 273.15]  # in K
    molar_volume = [3.57554070675e+26, 3.15822330675e+26, 2.6018001067499996e+26]
    a = 348.79
    k_b = 1.380649e+18

    k = (8*interfacial_energy[0]*((molar_volume[0])**2)*diffusion_coefficient[0])/(9*R*T[0]*(e**((2*interfacial_energy[0]*molar_volume[0])/(k_b*T[0]))))
    cu_wt = 1.04 / 100

    fe_density = 7.8e-24
    mol_cu = 63.546 * 1e-3
    cu_density = 8.85e-24
    a = 348.79

    const1 = (k * cu_wt * fe_density) / mol_cu
    const2 = (4 * pi)
    const3 = (k * cu_density) / ((mol_cu) * (a ** 3))
    const4 = e ** (np.cbrt((4 * pi) / 3))

    
    try:
        term1 = const1

        term2_numerator = (3 * e ** (np.cbrt((4 * pi) / v1[i-1]))) - v1[i-1] * e ** ((-v1[i-1] / 3) + np.cbrt((4 * pi) / 3))
        term2 = term2_numerator / const2

        #term3 = (const3 * t[i]) + (const3 * v[i-1])
        term3 = const3*t[i]
    

        if i > 0:
            term4 = (v1[i-1] * e ** (-v1[i] / 3)) * const4 / const2
        else:
            term4 = (v1[i] * e ** (-v1[i] / 3)) * const4 / const2 # This case is handled implicitly by starting from i=1

        dvdt1[i] = term1 - (const3*v1[i]) / (term2 + term3 + term4)

        # Update v using the computed dvdt
        v1[i] = v1[i-1] + dvdt1[i] * dt

    except ZeroDivisionError as zde:
        print(f"ZeroDivisionError at index {i}: {zde}")
        dvdt1[i] = np.nan  # Handle division by zero
    except Exception as ex:
        print(f"Error at index {i}: {ex}")
        dvdt1[i] = np.nan  # Handle other exceptions

    # Plot the results
plt.figure(figsize=(10, 6))
plt.plot(t, dvdt1, label='dvdt', marker='o')
plt.xlabel('Time t')
plt.ylabel('dvdt1')
plt.title('Plot of dvdt over Time')
plt.legend()
plt.grid(True)
plt.show()

r1 = np.cbrt((3/(4*pi))*(v1))

    # Plot v over time as well
plt.figure(figsize=(10, 6))
plt.plot(t, r1, label='r', marker='o')
plt.xlabel('Time t')
plt.ylabel('r')
plt.title('Plot of r1 over Time')
plt.legend()
plt.grid(True)
plt.show()

print(r1)




v2 = np.zeros_like(t)
v2[0] = v1[-1]


# Iterate to solve for v and dvdt
dvdt2 = np.zeros_like(t)

for i in range(1, len(t)):
    
    pi = np.pi
    e = np.exp(1)
    R = 8.314e+18  # Gas constant in J/(mol*K)

    interfacial_energy = [0.51*10**-19, 1.23*10**-19, 3.14*10**-19]  # in J/m^2
    diffusion_coefficient = [2.6e+11, 4e+9, 2e+6]  # in m^2/s
    T = [780 + 273.15, 660 + 273.15, 500 + 273.15]  # in K
    molar_volume = [3.57554070675e+26, 3.15822330675e+26, 2.6018001067499996e+26]
    a = 348.79
    k_b = 1.380649e+18

    k = (8*interfacial_energy[1]*((molar_volume[1])**2)*diffusion_coefficient[1])/(9*R*T[1]*(e**((2*interfacial_energy[1]*molar_volume[1])/(k_b*T[1]))))
    cu_wt = 1.04 / 100

    fe_density = 7.8e-24
    mol_cu = 63.546 * 1e-3
    cu_density = 8.85e-24
    a = 348.79

    const1 = (k * cu_wt * fe_density) / mol_cu
    const2 = (4 * pi)
    const3 = (k * cu_density) / ((mol_cu) * (a ** 3))
    const4 = e ** (np.cbrt((4 * pi) / 3))

    
    try:
        term1 = const1

        term2_numerator = (3 * e ** (np.cbrt((4 * pi) / v2[i-1]))) - v2[i-1] * e ** ((-v2[i-1] / 3) + np.cbrt((4 * pi) / 3))
        term2 = term2_numerator / const2

        #term3 = (const3 * t[i]) + (const3 * v[i-1])
        term3 = const3*t[i]
    

        if i > 0:
            term4 = (v2[i-1] * e ** (-v2[i] / 3)) * const4 / const2
        else:
            term4 = (v2[i] * e ** (-v2[i] / 3)) * const4 / const2 # This case is handled implicitly by starting from i=1

        dvdt2[i] = term1 - (const3*v2[i]) / (term2 + term3 + term4)

        # Update v using the computed dvdt
        v2[i] = v2[i-1] + dvdt2[i] * dt

    except ZeroDivisionError as zde:
        print(f"ZeroDivisionError at index {i}: {zde}")
        dvdt2[i] = np.nan  # Handle division by zero
    except Exception as ex:
        print(f"Error at index {i}: {ex}")
        dvdt2[i] = np.nan  # Handle other exceptions

    # Plot the results
plt.figure(figsize=(10, 6))
plt.plot(t, dvdt2, label='dvdt', marker='o')
plt.xlabel('Time t')
plt.ylabel('dvdt')
plt.title('Plot of dvdt1 over Time')
plt.legend()
plt.grid(True)
plt.show()


r2 = np.cbrt((3/(4*pi))*(v2))

    # Plot v over time as well
plt.figure(figsize=(10, 6))
plt.plot(t, r2, label='r', marker='o')
plt.xlabel('Time t')
plt.ylabel('r')
plt.title('Plot of r2 over Time')
plt.legend()
plt.grid(True)
plt.show()

print(r2)






v3 = np.zeros_like(t)
v3[0] = v2[-1]


# Iterate to solve for v and dvdt
dvdt3 = np.zeros_like(t)

for i in range(1, len(t)):
    
    pi = np.pi
    e = np.exp(1)
    R = 8.314e+18  # Gas constant in J/(mol*K)

    interfacial_energy = [0.51*10**-19, 1.23*10**-19, 3.14*10**-19]  # in J/m^2
    diffusion_coefficient = [2.6e+11, 4e+9, 2e+6]  # in m^2/s
    T = [780 + 273.15, 660 + 273.15, 500 + 273.15]  # in K
    molar_volume = [3.57554070675e+26, 3.15822330675e+26, 2.6018001067499996e+26]
    a = 348.79
    k_b = 1.380649e+18

    k = (8*interfacial_energy[2]*((molar_volume[2])**2)*diffusion_coefficient[2])/(9*R*T[2]*(e**((2*interfacial_energy[2]*molar_volume[2])/(k_b*T[2]))))
    cu_wt = 1.04 / 100

    fe_density = 7.8e-24
    mol_cu = 63.546 * 1e-3
    cu_density = 8.85e-24
    a = 348.79

    const1 = (k * cu_wt * fe_density) / mol_cu
    const2 = (4 * pi)
    const3 = (k * cu_density) / ((mol_cu) * (a ** 3))
    const4 = e ** (np.cbrt((4 * pi) / 3))

    
    try:
        term1 = const1

        term2_numerator = (3 * e ** (np.cbrt((4 * pi) / v3[i-1]))) - v3[i-1] * e ** ((-v3[i-1] / 3) + np.cbrt((4 * pi) / 3))
        term2 = term2_numerator / const2

        #term3 = (const3 * t[i]) + (const3 * v[i-1])
        term3 = const3*t[i]
    

        if i > 0:
            term4 = (v3[i-1] * e ** (-v3[i] / 3)) * const4 / const2
        else:
            term4 = (v3[i] * e ** (-v3[i] / 3)) * const4 / const2 # This case is handled implicitly by starting from i=1

        dvdt3[i] = term1 - (const3*v3[i]) / (term2 + term3 + term4)

        # Update v using the computed dvdt
        v3[i] = v3[i-1] + dvdt3[i] * dt

    except ZeroDivisionError as zde:
        print(f"ZeroDivisionError at index {i}: {zde}")
        dvdt3[i] = np.nan  # Handle division by zero
    except Exception as ex:
        print(f"Error at index {i}: {ex}")
        dvdt3[i] = np.nan  # Handle other exceptions

    # Plot the results
plt.figure(figsize=(10, 6))
plt.plot(t, dvdt3, label='dvdt3', marker='o')
plt.xlabel('Time t')
plt.ylabel('dvdt3')
plt.title('Plot of dvdt3 over Time')
plt.legend()
plt.grid(True)
plt.show()


r3 = np.cbrt((3/(4*pi))*(v3))
    # Plot v over time as well
plt.figure(figsize=(10, 6))
plt.plot(t, r3, label='r', marker='o')
plt.xlabel('Time t')
plt.ylabel('r')
plt.title('Plot of r over Time')
plt.legend()
plt.grid(True)
plt.show()

print(r3)





all_t = np.linspace(0, 90, 300000)

r = [*r1, *r2, *r3]



plt.figure(figsize=(10, 6))
plt.plot(all_t, r, label='r', marker='o')
plt.xlabel('Time t')
plt.ylabel('r')
plt.title('Plot of r over All Time')
plt.legend()
plt.grid(True)
plt.show()