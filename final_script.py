
import numpy as np
import matplotlib.pyplot as plt

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

# Time parameters
t = np.linspace(0, 30, 100000)
dt = t[1] - t[0]  # Timestep based on the linspace definition
    
# Initial conditions
rA0 = 0.128
r1 = np.zeros_like(t)
r1[0] = rA0


# Iterate to solve for v and dvdt
drdt1 = np.zeros_like(t)

for i in range(1, len(t)):
    
    alpha = (2*interfacial_energy[0]*mol_vol_fe[0])/(k_b*T[0])
    beta = (8*interfacial_energy[0]*((mol_vol_cu[0])**2)*diffusion_coefficient[0])/(9*R*T[0])
    
    try:
        term1 = (beta*wt_cu*fe_density)/(mol_cu)
        term2 = (beta*cu_density*4*pi*(r1[i]**3))/(3*(a**3)*mol_cu)
        
        term3 = np.exp(alpha/r1[i])  *  (3*r1[i]-alpha)
        
        if i > 0:
            term4 = np.exp(alpha/r1[i])*(r1[i-1]**3)*( (alpha/(r1[i]**2)) - 1 )
        else:
            term4 = np.exp(alpha/r1[i])*(r1[i]**3)*( (alpha/(r1[i]**2)) - 1 )
            
        term5 = (beta*cu_density*4*pi*(r1[i]**2)*t[i])/((a**3)*mol_cu)
        
        drdt1[i] = (term1 + term2) / (term3 + term4 + term5)

        # Update v using the computed dvdt
        r1[i] = r1[i-1] + drdt1[i] * dt

    except ZeroDivisionError as zde:
        print(f"ZeroDivisionError at index {i}: {zde}")
        drdt1[i] = np.nan  # Handle division by zero
    except Exception as ex:
        print(f"Error at index {i}: {ex}")
        drdt1[i] = np.nan  # Handle other exceptions

    # Plot the results
plt.figure(figsize=(10, 6))
plt.plot(t, drdt1, label='dvdt', marker='o')
plt.xlabel('Time t')
plt.ylabel('dvdt1')
plt.title('Plot of dvdt over Time')
plt.legend()
plt.grid(True)
plt.show()



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
