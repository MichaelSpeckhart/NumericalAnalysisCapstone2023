import matplotlib.pyplot as plt
import numpy as np
import os

# Create a folder named 'graphs' to save the plot
if not os.path.exists('graphs'):
    os.makedirs('graphs')

#threads = 2
size = [256,512,1024,2048,4096]
toMS = 1000;

serial = [6.0980,48,395,3.7396*toMS,30.295*toMS]
rayon = [5.9111,41,331,2.79*toMS,22.388*toMS]
vector = [3.5262,27,227,2.5*toMS,23.241*toMS]
noPivot = [6.6467,55,432,3.00*toMS,32.240*toMS]
rayonVectorNoPivot = [4.8619,31,251,2.1669*toMS,17.811*toMS]
rayonNoPivot = [6.5111,46,351,2.81*toMS,21.888*toMS]

#Set the figure size and font size
plt.figure(figsize=(8, 8))
plt.rcParams.update({'font.size': 12})

# Plot the data with a solid line and marker
plt.plot(size, serial, '-o', label="Serial", color='red')
plt.plot(size, rayon, '--s', label="Rayon", color='yellow')
plt.plot(size, vector, '-.D', label="Vector", color='blue')
plt.plot(size, noPivot, ':^', label="No Pivot", color='orange')
plt.plot(size, rayonNoPivot, ':*', label="Rayon+NoPivot", color='purple')
plt.plot(size, rayonVectorNoPivot, '-.+', label="Rayon+Vector+NoPivot", color='green')

# Add axis labels and title
plt.xlabel('Size of Coefficient Matrix')
plt.ylabel('Execution Time (ms)')
plt.suptitle("Size of Coefficient Matrix Vs. Execution Time")
plt.title("2 Threads")

# Add gridlines and legend
plt.grid(False)
plt.legend(loc='upper left')

# Save the plot to the 'graphs' folder with the name 'plot.png'
plt.savefig('graphs/SizeVsExecutionThreadTwo.png')

# Get rid of last element
size = [256,512,1024,2048]
toMS = 1000;

serial = [6.0980,48,395,3.7396*toMS]
rayon = [5.9111,41,331,2.79*toMS]
vector = [3.5262,27,227,2.5*toMS]
noPivot = [6.6467,55,432,3.00*toMS]
rayonVectorNoPivot = [4.8619,31,251,2.1669*toMS]
rayonNoPivot = [6.5111,46,351,2.81*toMS]

#Set the figure size and font size
plt.figure(figsize=(8, 8))
plt.rcParams.update({'font.size': 12})

# Plot the data with a solid line and marker
plt.plot(size, serial, '-o', label="Serial", color='red')
plt.plot(size, rayon, '--s', label="Rayon", color='yellow')
plt.plot(size, vector, '-.D', label="Vector", color='blue')
plt.plot(size, noPivot, ':^', label="No Pivot", color='orange')
plt.plot(size, rayonNoPivot, ':*', label="Rayon+NoPivot", color='purple')
plt.plot(size, rayonVectorNoPivot, '-.+', label="Rayon+Vector+NoPivot", color='green')

# Add axis labels and title
plt.xlabel('Size of Coefficient Matrix')
plt.ylabel('Execution Time (ms)')
plt.suptitle("Size of Coefficient Matrix Vs. Execution Time")
plt.title("2 Threads, no 4096")

# Add gridlines and legend
plt.grid(False)
plt.legend(loc='upper left')

# Save the plot to the 'graphs' folder with the name 'plot.png'
plt.savefig('graphs/SizeVsExecutionThreadTwoNo4096.png')

#threads = 2, log y axis
size = [256,512,1024,2048,4096]
toMS = 1000;

serial = [6.0980,48,395,3.7396*toMS,30.295*toMS]
rayon = [5.9111,41,331,2.79*toMS,22.388*toMS]
vector = [3.5262,27,227,2.5*toMS,23.241*toMS]
noPivot = [6.6467,55,432,3.00*toMS,32.240*toMS]
rayonVectorNoPivot = [4.8619,31,251,2.1669*toMS,17.811*toMS]
rayonNoPivot = [6.5111,46,351,2.81*toMS,21.888*toMS]

#Set the figure size and font size
plt.figure(figsize=(8, 8))
plt.rcParams.update({'font.size': 12})

# Plot the data with a solid line and marker
plt.plot(size, serial, '-o', label="Serial", color='red')
plt.plot(size, rayon, '--s', label="Rayon", color='yellow')
plt.plot(size, vector, '-.D', label="Vector", color='blue')
plt.plot(size, noPivot, ':^', label="No Pivot", color='orange')
plt.plot(size, rayonNoPivot, ':*', label="Rayon+NoPivot", color='purple')
plt.plot(size, rayonVectorNoPivot, '-.+', label="Rayon+Vector+NoPivot", color='green')
plt.yscale('log')

# Add axis labels and title
plt.xlabel('Size of Coefficient Matrix')
plt.ylabel('Execution Time (ms) LOG')
plt.suptitle("Size of Coefficient Matrix Vs. Execution Time")
plt.title("2 Threads")

# Add gridlines and legend
plt.grid(False)
plt.legend(loc='upper left')

# Save the plot to the 'graphs' folder with the name 'plot.png'
plt.savefig('graphs/SizeVsExecutionThreadTwoLogY.png')

# Your data
size = [256, 512, 1024, 2048, 4096]
toMS = 1000
serial = [6.0980, 48, 395, 3.7396 * toMS, 29.902 * toMS]
rayon = [7.1465, 29, 178, 1.5578 * toMS, 14.679 * toMS]
vector = [3.5262, 27, 227, 2.49 * toMS, 21.680 * toMS]
noPivot = [6.6467, 55, 432, 4.0131 * toMS, 32.26 * toMS]
rayonVectorNoPivot = [6.642, 25, 145, 1.5187 * toMS, 14.467 * toMS]
rayonNoPivot = [7.23,28,168,1.54*toMS,14.5*toMS]

# Set the figure size and font size
plt.figure(figsize=(8, 6))
plt.rcParams.update({'font.size': 12})

# Plot the data with different line styles and markers
plt.plot(size, serial, '-o', label="Serial", color='red')
plt.plot(size, rayon, '--s', label="Rayon", color='yellow')
plt.plot(size, vector, '-.D', label="Vector", color='blue')
plt.plot(size, noPivot, ':^', label="No Pivot", color='orange')
plt.plot(size, rayonNoPivot, ':*', label="Rayon+NoPivot", color='purple')
plt.plot(size, rayonVectorNoPivot, '-.+', label="Rayon+Vector+NoPivot", color='green')


# Add axis labels and titles
plt.xlabel('Size of Coefficient Matrix')
plt.ylabel('Execution Time (ms)')
plt.suptitle("Size of Coefficient Matrix Vs. Execution Time")
plt.title("8 Threads")

# Add gridlines and legend
plt.grid(False)
plt.legend(loc='upper left')

# Save the plot to the 'graphs' folder with the name 'plot.png'
plt.savefig('graphs/SizeVsExecutionThreadEightLinearY.png')

plt.ylabel('Execution Time (ms) LOG')
plt.yscale('log')
plt.savefig('graphs/SizeVsExecutionThreadEightLogY.png')


# Your data
thread = [1, 2, 4, 8]
serial = [30, 30,30,30]
rayon = [26.2896, 22.338,15.684,14.679]
vector = [21.751,21.751,21.751,21.751]
noPivot = [32.053, 32.053,32.053,32.053]
rayonVectorNoPivot = [22.342, 17.841,14.200,14.467]
rayonNoPivot = [25.21, 22.56, 15.102, 14.46]

# Set the figure size and font size
plt.figure(figsize=(8, 6))
plt.rcParams.update({'font.size': 12})

# Plot the data with different line styles and markers
plt.plot(thread, serial, '-o', label="Serial", color='red')
plt.plot(thread, rayon, '--s', label="Rayon", color='yellow')
plt.plot(thread, vector, '-.D', label="Vector", color='blue')
plt.plot(thread, noPivot, ':^', label="No Pivot", color='orange')
plt.plot(thread, rayonNoPivot, ':*', label="Rayon+NoPivot", color='purple')
plt.plot(thread, rayonVectorNoPivot, '-.+', label="Rayon+Vector+NoPivot", color='green')


# Add axis labels and titles
plt.xlabel('Number of Threads')
plt.ylabel('Execution Time (Seconds)')
plt.suptitle("Number of Threads Vs. Execution Time")
plt.title("Matrix Size 4096")

# Add gridlines and legend
plt.grid(False)
plt.legend(loc='upper right')

# Save the plot to the 'graphs' folder with the name 'plot.png'
plt.savefig('graphs/NumberOfThreadsVsExecTimeSize4096.png')

# Your data
thread = [1, 2, 4, 8]
serial = [448, 448,448,448]
rayon = [335, 331,222,178]
vector = [238,238,238,238]
noPivot = [434, 434,434,434]
rayonVectorNoPivot = [257, 251,173,145]
rayonNoPivot = [332, 329,220,173]

# Set the figure size and font size
plt.figure(figsize=(8, 6))
plt.rcParams.update({'font.size': 12})

# Plot the data with different line styles and markers
plt.plot(thread, serial, '-o', label="Serial", color='red')
plt.plot(thread, rayon, '--s', label="Rayon", color='yellow')
plt.plot(thread, vector, '-.D', label="Vector", color='blue')
plt.plot(thread, noPivot, ':^', label="No Pivot", color='orange')
plt.plot(thread, rayonNoPivot, ':*', label="Rayon+NoPivot", color='purple')
plt.plot(thread, rayonVectorNoPivot, '-.+', label="Rayon+Vector+NoPivot", color='green')


# Add axis labels and titles
plt.xlabel('Number of Threads')
plt.ylabel('Execution Time (ms)')
plt.suptitle("Number of Threads Vs. Execution Time")
plt.title("Matrix Size 1024")

# Add gridlines and legend
plt.grid(False)
plt.legend(loc='upper right')

# Save the plot to the 'graphs' folder with the name 'plot.png'
plt.savefig('graphs/NumberOfThreadsVsExecTimeSize1024.png')
