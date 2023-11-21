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
ssor10 = [0.00374438,0.00890029,0.0143356,0.0182472,0.0232362,0.0262019,0.0296359,0.0333711]
gauss_sidel = [0.0033514,0.00664757,0.0102938,0.0152454,0.021328,0.0246872,0.0281842,0.0314516]
gauss_elimination = [0.00811063,0.0642115,0.245409,0.7746,1.71589,3.13128,5.42308,8.46707]

#Set the figure size and font size
plt.figure(figsize=(6, 4))
plt.rcParams.update({'font.size': 10})

# Plot the data with a solid line and marker
plt.plot(size, ssor10, ':^', label="SPARSE SSOR w = 0.5", color='orange', markersize=8, linewidth=2)
plt.plot(size, jacobi, '-o', label="SPARSE Jacobi Method", color='red', markersize=8, linewidth=2)
plt.plot(size, gauss_sidel, '-.D', label="SPARSE Gauss Sidel", color='blue', markersize=8, linewidth=2)
plt.plot(size, gauss_elimination, '--s', label="DENSE Gauss Eliminaton", color='green', markersize=8, linewidth=2)
#plt.plot(size, ssor15, ':*', label="SSOR w = 1.5", color='purple', markersize=8, linewidth=2)
# plt.plot(size, rayonVectorNoPivot, '-.+', label="Rayon+Vector+NoPivot", color='yellow', markersize=8, linewidth=2)

# Add axis labels and title
plt.xlabel('Size of Coefficient Matrix')
plt.ylabel('Log Execution Time (s)')
plt.suptitle("Sparse Vs. Dense Data Structure")
plt.title("A is Tridiagonal, B is First Canonical Basis Vector")

# Add gridlines and legend
plt.grid(False)
plt.legend(loc='upper left')
plt.yscale("log")   

# Save the plot to the 'graphs' folder with the name 'plot.png'
plt.savefig('graphs/SparseVsDenseSolvers.png', dpi=600)

size = [5000,10000,15000,20000,25000,30000,35000,40000]

jacobi = [0.0182609,0.0374146,0.0554051,0.0732703,0.0924255,0.111208,0.130893,0.151109]
ssor10 = [ 0.0415163,0.0789963,0.117733,0.155531,0.194238,0.232461,0.270923,0.309714]
gauss_sidel = [0.0380204,0.0720257,0.105602,0.138193,0.171957,0.205422,0.23823,0.27187]

#Set the figure size and font size
plt.figure(figsize=(6, 4))
plt.rcParams.update({'font.size': 10})

# Plot the data with a solid line and marker
plt.plot(size, ssor10, ':^', label="SPARSE SSOR w = 0.5", color='orange', markersize=8, linewidth=2)
plt.plot(size, jacobi, '-o', label="SPARSE Jacobi Method", color='red', markersize=8, linewidth=2)
plt.plot(size, gauss_sidel, '-.D', label="SPARSE Gauss Sidel", color='blue', markersize=8, linewidth=2)
#plt.plot(size, gauss_elimination, '--s', label="DENSE Gauss Eliminaton", color='green', markersize=8, linewidth=2)
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
plt.savefig('graphs/SparseSolvers.png', dpi=600)
#THREADS
threads = [1,2,3,4,5,6,7,8]

jacobiParallel = [1.90347,1.02416,0.73991,0.640031,0.606254,0.599996,0.599924,0.605556]
jacobiSerial = [1.71814,1.71814,1.71814,1.71814,1.71814,1.71814,1.71814,1.71814]
Gauss_sidel_Parallel = [1.89948,1.01479,0.730396,0.641583,0.601846,0.600159,0.599349,0.60589,]
Gauss_sidel_Serial = [1.70271,1.70271,1.70271,1.70271,1.70271,1.70271,1.70271,1.70271]
SSOR_Parallel = [ 1.71467,0.921314,0.682119,0.617948,0.592276,0.59818,0.597187,0.603722]
SSOR_Serial = [1.78209,1.78209,1.78209,1.78209,1.78209,1.78209,1.78209,1.78209]

#Set the figure size and font size
plt.figure(figsize=(6, 4))
plt.rcParams.update({'font.size': 10})

# Plot the data with a solid line and marker
plt.plot(threads, jacobiParallel, ':^', label="Parallel Sparse Jacobi", color='orange', markersize=8, linewidth=2)
plt.plot(threads, jacobiSerial, '-o', label="Sparse Jacobi", color='red', markersize=8, linewidth=2)
plt.plot(threads, Gauss_sidel_Parallel, '-.D', label="Parallel Sparse Gauss Sidel", color='blue', markersize=8, linewidth=2)
plt.plot(threads, Gauss_sidel_Serial, '--s', label="Sparse Gauss Sidel", color='green', markersize=8, linewidth=2)
plt.plot(threads, SSOR_Parallel, ':*', label="Parallel Sparse SSOR (W = 0.5)", color='purple', markersize=8, linewidth=2)
plt.plot(threads, SSOR_Serial, '-.+', label="Sparse SSOR (W = 0.5)", color='yellow', markersize=8, linewidth=2)

# Add axis labels and title
plt.xlabel('Threads')
plt.ylabel('Execution Time (s)')
plt.suptitle("Number of Threads Vs. Execution Time")
plt.title("1000 iterations on Thread Matrix (NNZ = 4,444,880)")

# Add gridlines and legend
plt.grid(False)
plt.legend(loc='best')

# Save the plot to the 'graphs' folder with the name 'plot.png'
plt.savefig('graphs/ThreadsIterativeSolvers.png', dpi=600)

threads = [1,2,3,4,5,6,7,8]

multParallel = [0.591281,0.310774,0.212397,0.166407,0.134634,0.11685,0.098649,0.0895234]
multSerial = [0.542035,0.542035,0.542035,0.542035,0.542035,0.542035,0.542035,0.542035]

#Set the figure size and font size
plt.figure(figsize=(6, 4))
plt.rcParams.update({'font.size': 10})

# Plot the data with a solid line and marker
plt.plot(threads, multParallel, '-o', label="Parallel Sparse Matrix Multiplication", color='red', markersize=8, linewidth=2)
plt.plot(threads, multSerial, '-.D', label="Sparse Matrix Multiplication", color='blue', markersize=8, linewidth=2)

# Add axis labels and title
plt.xlabel('Threads')
plt.ylabel('Execution Time (s)')
plt.suptitle("Number of Threads Vs. Execution Time")
plt.title("s3rmt3m3 ^2 (NNZ = 207,123)")

# Add gridlines and legend
plt.grid(False)
plt.legend(loc='best')

# Save the plot to the 'graphs' folder with the name 'plot.png'
plt.savefig('graphs/ThreadsSparseMultiplication.png', dpi=600)