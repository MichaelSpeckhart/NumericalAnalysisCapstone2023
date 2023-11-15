import matplotlib.pyplot as plt
import numpy as np
import os

# Create a folder named 'graphs' to save the plot
if not os.path.exists('graphs'):
    os.makedirs('graphs')

#threads = 2
size = [500,1000,1500,2000,2500,3000,3500,4000]

jacobi = [0.184565,0.759483,1.82614,3.27084,5.12254,7.39125,10.0981,13.2249]
ssor10 = [0.373677,1.53101,3.62324,6.71153,10.5258,15.151,20.6353,26.8876]
gauss_sidel = [0.156668,0.69557,1.70365,3.12661,4.96193,7.21092,9.8265,12.9168]
gauss_elimination = [0.00811063,0.0642115,0.245409,0.7746,1.71589,3.13128,5.42308,8.46707]

#Set the figure size and font size
plt.figure(figsize=(8, 8))
plt.rcParams.update({'font.size': 12})

# Plot the data with a solid line and marker
plt.plot(size, ssor10, ':^', label="SSOR w = 1.0", color='orange')
plt.plot(size, jacobi, '-o', label="Jacobi Method", color='red')
plt.plot(size, gauss_sidel, '-.D', label="Gauss Sidel", color='blue')
plt.plot(size, gauss_elimination, '--s', label="Gauss Eliminaton", color='yellow')
#plt.plot(size, ssor15, ':*', label="SSOR w = 1.5", color='purple')
# plt.plot(size, rayonVectorNoPivot, '-.+', label="Rayon+Vector+NoPivot", color='green')

# Add axis labels and title
plt.xlabel('Size of Coefficient Matrix')
plt.ylabel('Execution Time (s)')
plt.suptitle("Size of System Vs. Execution Time")
plt.title("A is Tridiagonal, B is First Canonical Basis Vector, Dense Data Structure")

# Add gridlines and legend
plt.grid(False)
plt.legend(loc='upper left')

# Save the plot to the 'graphs' folder with the name 'plot.png'
plt.savefig('graphs/DenseSolvers.png')