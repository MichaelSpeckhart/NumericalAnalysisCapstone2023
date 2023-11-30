import matplotlib.pyplot as plt
import numpy as np
import os

# Create a folder named 'graphs' to save the plot
if not os.path.exists('graphs'):
    os.makedirs('graphs')

#threads = 2
size = [500,1000,1500,2000,2500,3000,3500,4000]

jacobi = [0.00378789,0.0148882,0.0359204,0.0654168,0.102179,0.147571,0.201174,0.263038]
ssor05 = [0.0129413,0.051301,0.122392,0.227018,0.355492,0.511165,0.699199,0.912121]
gauss_sidel = [0.00230138,0.00953875,0.0234736,0.0434253,0.0679959,0.0983122,0.135707,0.177838]
gauss_elimination = [0.00803985,0.0615053,0.239279,0.774963,1.76108,3.29759,5.50983,8.48077]

#Set the figure size and font size
plt.figure(figsize=(6, 4))
plt.rcParams.update({'font.size': 10})

# Plot the data with a solid line and marker
plt.plot(size, gauss_elimination, '--s', label="DENSE Gauss Eliminaton", color='green', markersize=8, linewidth=2)
plt.plot(size, ssor05, ':^', label="DENSE SSOR w = 0.5", color='orange', markersize=8, linewidth=2)
plt.plot(size, jacobi, '-o', label="DENSE Jacobi Method", color='red', markersize=10, linewidth=3)
plt.plot(size, gauss_sidel, '-.D', label="DENSE Gauss Seidel", color='blue', markersize=8, linewidth=2)
#plt.plot(size, ssor15, ':*', label="SSOR w = 1.5", color='purple', markersize=8, linewidth=2)
# plt.plot(size, rayonVectorNoPivot, '-.+', label="Rayon+Vector+NoPivot", color='gray', markersize=8, linewidth=2)

# Add axis labels and title
plt.xlabel('Size of Coefficient Matrix')
plt.ylabel('Execution Time (s)')
plt.suptitle("Size of System Vs. Execution Time")
plt.title("A is Tridiagonal, B is Vector of Ones")

# Add gridlines and legend
plt.grid(False)
plt.legend(loc='upper left')
plt.yscale("log")   


# Save the plot to the 'graphs' folder with the name 'plot.png'
plt.savefig('graphs/DenseSolvers.png',dpi=600)

size = [5000,10000,15000,20000,25000,30000,35000,40000]

jacobi = [0.408604,1.66177,3.7443,6.81123,10.6548,15.3198,20.9042,27.3564]
ssor05 = [1.42598,5.68564,12.7954,22.8108,35.4605,51.0806,70.1463,90.9169]
gauss_sidel = [0.279004,1.12293,2.53999,4.54094,7.10162,10.2421,13.9765,18.1447]
jacobiSparse = [0.000541976,0.00095077,0.00136434,0.00172531,0.00210362,0.00246806,0.00317477,0.00411002]
ssor05Sparse =  [0.00153812,0.00312319,0.00455442,0.00606881,0.00751245,0.00905234,0.0106831,0.0123361]
gauss_sidelSparse = [0.00067122,0.00133102,0.00194759,0.00258418,0.00324621,0.00374712,0.00463873,0.00532041]


#Set the figure size and font size
plt.figure(figsize=(6, 4))
plt.rcParams.update({'font.size': 10})

# Plot the data with a solid line and marker
plt.plot(size, ssor05, ':^', label="DENSE SSOR w = 0.5", color='orange', markersize=8, linewidth=2)
plt.plot(size, jacobi, '-o', label="DENSE Jacobi Method", color='red', markersize=8, linewidth=2)
plt.plot(size, gauss_sidel, '-.D', label="DENSE Gauss Seidel", color='blue', markersize=8, linewidth=2)
plt.plot(size, ssor05Sparse, ':*', label="SPARSE SSOR w = 0.5", color='purple', markersize=8, linewidth=2)
plt.plot(size, gauss_sidelSparse, '-.+', label="SPARSE Gauss Seidel", color='gray', markersize=8, linewidth=2)
plt.plot(size, jacobiSparse, '--s', label="SPARSE Jacobi Method", color='green', markersize=8, linewidth=2)

# Add axis labels and title
plt.xlabel('Size of Coefficient Matrix')
plt.ylabel('Log Execution Time (s)')
plt.suptitle("Sparse Vs. Dense Data Structure")
plt.title("A is Tridiagonal, B is Vector of Ones")

# Add gridlines and legend
plt.grid(False)
plt.legend(loc='best')
plt.yscale("log")   

# Save the plot to the 'graphs' folder with the name 'plot.png'
plt.savefig('graphs/SparseVsDenseSolvers.png', dpi=600)

size = [50000,100000,150000,200000,250000,300000,350000,400000]

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
# plt.plot(size, rayonVectorNoPivot, '-.+', label="Rayon+Vector+NoPivot", color='gray', markersize=8, linewidth=2)

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
plt.plot(threads, SSOR_Serial, '-.+', label="Sparse SSOR (W = 0.5)", color='gray', markersize=8, linewidth=2)
plt.plot(threads, jacobiSerial, '-o', label="Sparse Jacobi", color='red', markersize=8, linewidth=2)
plt.plot(threads, Gauss_sidel_Serial, '--s', label="Sparse Gauss Sidel", color='green', markersize=8, linewidth=2)
plt.plot(threads, jacobiParallel, ':^', label="Parallel Sparse Jacobi", color='orange', markersize=9, linewidth=2)
plt.plot(threads, Gauss_sidel_Parallel, '-.D', label="Parallel Sparse Gauss Sidel", color='blue', markersize=7, linewidth=2)
plt.plot(threads, SSOR_Parallel, ':*', label="Parallel Sparse SSOR (W = 0.5)", color='purple', markersize=8, linewidth=2)

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

size = [500000,1000000,1500000,2000000,2500000,3000000,3500000,4000000]

jacobiSparse =  [0.0502345,0.100394,0.15189,0.20853,0.261648,0.319298,0.373895,0.425974]
ssor05Sparse =  [0.159263,0.319445,0.475897,0.643693,0.807765,0.977107,1.14527,1.31315]
gauss_sidelSparse =  [0.0671957,0.133561,0.195878,0.264329,0.331295,0.396536,0.462269,0.530184]
jacobiParallel = [0.048563,0.0924779,0.139606,0.188805,0.237797,0.287245,0.332312,0.381289]
ssor05Parallel = [0.0798771,0.17022,0.266174,0.360392,0.452413,0.550043,0.636017,0.734703]
gauss_sidelParallel = [0.0322753,0.0630284,0.101977,0.139762,0.175849,0.198169,0.243839,0.280697]


#Set the figure size and font size
plt.figure(figsize=(6, 4))
plt.rcParams.update({'font.size': 10})

# Plot the data with a solid line and marker
plt.plot(size, ssor05Sparse, ':*', label="SSOR w = 0.5", color='purple', markersize=8, linewidth=2)
plt.plot(size, ssor05Parallel, '-o', label="PARALLEL SSOR w = 0.5", color='red', markersize=8, linewidth=2)
plt.plot(size, gauss_sidelSparse, '-.+', label="Gauss Seidel", color='gray', markersize=8, linewidth=2)
plt.plot(size, jacobiSparse, '--s', label="Jacobi Method", color='green', markersize=8, linewidth=2)
plt.plot(size, jacobiParallel, ':^', label="PARALLEL Jacobi Method", color='orange', markersize=8, linewidth=2)
plt.plot(size, gauss_sidelParallel, '-.D', label="PARALLEL Gauss Seidel", color='blue', markersize=8, linewidth=2)


# Add axis labels and title
plt.xlabel('Size of Coefficient Matrix')
plt.ylabel('Log Execution Time (s)')
plt.suptitle("Parallel Vs. Serial for Sparse Data Structure")
plt.title("A is Tridiagonal, B is Vector of Ones")

# Add gridlines and legend
plt.grid(False)
plt.legend(loc='best')
plt.yscale("log")   

# Save the plot to the 'graphs' folder with the name 'plot.png'
plt.savefig('graphs/SparseVsParallelSolvers.png', dpi=600)
