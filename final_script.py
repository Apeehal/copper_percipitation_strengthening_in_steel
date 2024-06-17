import numpy as np
import matplotlib.pyplot as plt

# Define constants and initial conditions
k = 10000000000000000000000000000000
cu_wt = 1.04 / 100

fe_density = 7800 / (1e-27)
mol_cu = 63.546 * 1e-3
cu_density = 8850 / (1e-27)
a = 348.79
e = np.exp(1)
pi = np.pi

# Time parameters
t = np.linspace(0, 1, 10000)
dt = t[1] - t[0]  # Timestep based on the linspace definition

# Initial conditions
v0 = (0.128 * 10) ** 3 * (4 / 3)
v = np.zeros_like(t)
v[0] = v0

# Precompute some constants
const1 = (k * cu_wt * fe_density) / mol_cu
const2 = (4 * pi)
const3 = (k * cu_density) / ((mol_cu) * (a ** 3))
const4 = e ** (np.cbrt((4 * pi) / 3))

# Iterate to solve for v and dvdt
dvdt = np.zeros_like(t)

for i in range(1, len(t)):
    try:
        term1 = const1

        term2_numerator = (3 * e ** (np.cbrt((4 * pi) / v[i-1]))) - v[i-1] * e ** ((-v[i-1] / 3) + np.cbrt((4 * pi) / 3))
        term2 = term2_numerator / const2

        term3 = (const3 * t[i]) + (const3 * v[i-1])

        if i > 0:
            term4 = (v[i-1] * e ** (-v[i] / 3)) * const4 / const2
        else:
            term4 = (v[i] * e ** (-v[i] / 3)) * const4 / const2 # This case is handled implicitly by starting from i=1

        dvdt[i] = term1 / (term2 + term3 + term4)

        # Update v using the computed dvdt
        v[i] = v[i-1] + dvdt[i] * dt

    except ZeroDivisionError as zde:
        print(f"ZeroDivisionError at index {i}: {zde}")
        dvdt[i] = np.nan  # Handle division by zero
    except Exception as ex:
        print(f"Error at index {i}: {ex}")
        dvdt[i] = np.nan  # Handle other exceptions

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(t, dvdt, label='dvdt', marker='o')
plt.xlabel('Time t')
plt.ylabel('dvdt')
plt.title('Plot of dvdt over Time')
plt.legend()
plt.grid(True)
plt.show()

r = np.cbrt((3/(4*pi))*(v))

# Plot v over time as well
plt.figure(figsize=(10, 6))
plt.plot(t, r, label='v', marker='o')
plt.xlabel('Time t')
plt.ylabel('v')
plt.title('Plot of v over Time')
plt.legend()
plt.grid(True)
plt.show()

print(r)
