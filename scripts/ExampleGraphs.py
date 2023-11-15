import matplotlib.pyplot as plt
import numpy as np
import os

# Create a folder named 'graphs' to save the plot
if not os.path.exists('graphs'):
    os.makedirs('graphs')

#threads = 2
size = [1000,1500,2000,2500,3000,3500,4000,4500]

jacobi = [0.759125,1.82483,3.27048,5.11037,7.40779,10.135,13.2603,16.7436]
gauss_elimination = [0.0643407,0.24656,0.777767,1.74874,3.27369,5.46504,8.42173,12.3634]

#Set the figure size and font size
plt.figure(figsize=(8, 8))
plt.rcParams.update({'font.size': 12})

# Plot the data with a solid line and marker
plt.plot(size, jacobi, '-o', label="Jacobi Method", color='red')
plt.plot(size, gauss_elimination, '--s', label="Gauss Eliminaton", color='yellow')
# plt.plot(size, vector, '-.D', label="Vector", color='blue')
# plt.plot(size, noPivot, ':^', label="No Pivot", color='orange')
# plt.plot(size, rayonNoPivot, ':*', label="Rayon+NoPivot", color='purple')
# plt.plot(size, rayonVectorNoPivot, '-.+', label="Rayon+Vector+NoPivot", color='green')

# Add axis labels and title
plt.xlabel('Size of Coefficient Matrix')
plt.ylabel('Execution Time (s)')
plt.suptitle("Size of System Vs. Execution Time")
plt.title("A is tridiagonal, B is first canonical basis vector")

# Add gridlines and legend
plt.grid(False)
plt.legend(loc='upper left')

# Save the plot to the 'graphs' folder with the name 'plot.png'
plt.savefig('graphs/DenseSolvers.png')