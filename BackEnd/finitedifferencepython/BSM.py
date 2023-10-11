import numpy as np
import matplotlib.pyplot as plt

# Define parameters
T = 1.0  # final time
r = 0.15  # risk-free interest rate
sigma = 0.2  # volatility
K = 5.0  # strike price

# Define grid parameters
N = 50  # Number of price steps
M = 5000  # Number of time steps
S_max = 10.0  # Maximum value of S

# Calculate step sizes
dt = T / M
dS = S_max / N

# Initialize price and time grids
S = np.linspace(0, S_max, N+1)
t = np.linspace(0, T, M+1)

# Initialize option value matrix
V = np.zeros((N+1, M+1))

# Set boundary conditions
V[:, M] = np.maximum(S - K, 0)
V[0, :] = 0
V[N, :] = S_max - K * np.exp(-r * (T - t))

for j in range(M-1, -1, -1):
    for n in range(1, N):
        V[n, j] = V[n, j+1] - dt * (
            r * V[n, j+1] - r * S[n] * (V[n+1, j+1] - V[n-1, j+1]) / (2 * dS) - 
            0.5 * sigma**2 * S[n]**2 * (V[n+1, j+1] - 2 * V[n, j+1] + V[n-1, j+1]) / dS**2
        )

# Generate the 3D graph
S_mesh, T_mesh = np.meshgrid(S, t)
V_mesh = V.T  # Transpose for correct dimensions
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
surf = ax.plot_surface(S_mesh, T_mesh, V_mesh, cmap='viridis')
ax.set_xlabel('S')
ax.set_ylabel('t')
ax.set_zlabel('V')
ax.set_title('Option price 3D Graph')
plt.show()