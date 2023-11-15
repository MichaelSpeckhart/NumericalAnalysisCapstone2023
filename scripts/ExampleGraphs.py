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
plt.figure(figsize=(6, 4))
plt.rcParams.update({'font.size': 10})

# Plot the data with a solid line and marker
plt.plot(size, ssor10, ':^', label="DENSE SSOR w = 0.5", color='orange', markersize=8, linewidth=2)
plt.plot(size, jacobi, '-o', label="DENSE Jacobi Method", color='red', markersize=10, linewidth=3)
plt.plot(size, gauss_sidel, '-.D', label="DENSE Gauss Sidel", color='blue', markersize=8, linewidth=2)
plt.plot(size, gauss_elimination, '--s', label="DENSE Gauss Eliminaton", color='green', markersize=8, linewidth=2)
#plt.plot(size, ssor15, ':*', label="SSOR w = 1.5", color='purple', markersize=8, linewidth=2)
# plt.plot(size, rayonVectorNoPivot, '-.+', label="Rayon+Vector+NoPivot", color='yellow', markersize=8, linewidth=2)

# Add axis labels and title
plt.xlabel('Size of Coefficient Matrix')
plt.ylabel('Execution Time (s)')
plt.suptitle("Size of System Vs. Execution Time")
plt.title("A is Tridiagonal, B is First Canonical Basis Vector")

# Add gridlines and legend
plt.grid(False)
plt.legend(loc='upper left')

# Save the plot to the 'graphs' folder with the name 'plot.png'
plt.savefig('graphs/DenseSolvers.png',dpi=600)

size = [500,1000,1500,2000,2500,3000,3500,4000]

jacobi = [0.00183877,0.00351509,0.00519228,0.00718052,0.00909835,0.0105637,0.0124069,0.0140159]
gauss_sidel = [0.0033514,0.00664757,0.0102938,0.0152454,0.021328,0.0246872,0.0281842,0.0314516]
gauss_elimination = [0.00811063,0.0642115,0.245409,0.7746,1.71589,3.13128,5.42308,8.46707]

#Set the figure size and font size
plt.figure(figsize=(6, 4))
plt.rcParams.update({'font.size': 10})

# Plot the data with a solid line and marker
#plt.plot(size, ssor10, ':^', label="SSOR w = 1.0", color='orange')
plt.plot(size, jacobi, '-o', label="SPARSE Jacobi Method", color='red')
plt.plot(size, gauss_sidel, '-.D', label="SPARSE Gauss Sidel", color='blue')
plt.plot(size, gauss_elimination, '--s', label="DENSE Gauss Eliminaton", color='green')
#plt.plot(size, ssor15, ':*', label="SSOR w = 1.5", color='purple')
# plt.plot(size, rayonVectorNoPivot, '-.+', label="Rayon+Vector+NoPivot", color='yellow')

# Add axis labels and title
plt.xlabel('Size of Coefficient Matrix')
plt.ylabel('Execution Time (s)')
plt.suptitle("Sparse Vs. Dense Data Structure")
plt.title("A is Tridiagonal, B is First Canonical Basis Vector")

# Add gridlines and legend
plt.grid(False)
plt.legend(loc='upper left')

# Save the plot to the 'graphs' folder with the name 'plot.png'
plt.savefig('graphs/SparseSolvers.png', dpi=600)