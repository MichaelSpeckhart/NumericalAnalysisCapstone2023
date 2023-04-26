import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess

# Create a folder named 'graphs' to save the plot
if not os.path.exists('graphs'):
    os.makedirs('graphs')
#threads
threads = range(1,13)
#1500 by 1500
# parallel = [0.660286,0.349788,0.291875,0.308065,0.283758,0.313248,0.287639,0.311177,0.297052,0.291961,0.287548,0.300015]
# serial = [0.657983,0.637656,0.632463,0.631624,0.648794,0.641446,0.627549,0.627305,0.631666,0.632847,0.628903,0.634367]

# 2000 by 2000
# parallel = [1.87535,1.14586,1.10123,1.10648,1.10977,1.13271,1.13092,1.14381,1.16512,1.18542,1.19614,1.23085]
# serial = [1.83777,1.82945,1.82719,1.84386,1.84174,1.83843,1.83992,1.84405,1.84712,1.85111,1.8325,1.86247]

# 2000 by 2000 no TBB parallel for the devide
parallel = [1.87535,1.14586,1.10123,1.10648,1.10977,1.13271,1.13092,1.14381,1.16512,1.18542,1.19614,1.23085]
serial = [1.83777,1.82945,1.82719,1.84386,1.84174,1.83843,1.83992,1.84405,1.84712,1.85111,1.8325,1.86247]


# Set the figure size and font size
plt.figure(figsize=(8, 6))
plt.rcParams.update({'font.size': 12})

# Plot the data with a solid line and marker
plt.plot(threads, serial, '-o', label="serial", color='blue')
plt.plot(threads, parallel, '-o', label="parallel", color='red')



# Add axis labels and title
plt.xlabel('Threads')
plt.ylabel('Execution Time (seconds)')
plt.suptitle("Threads vs. Execution Time")
plt.title("2000 by 2000 dense matrix")

# Add gridlines and legend
plt.grid(True)
plt.legend(loc='upper right')

# Save the plot to the 'graphs' folder with the name 'plot.png'
plt.savefig('graphs/threadNoTBBDevide.png')

# #matrix size
# threads = [500,1000,1500,2000,2500,3000]
# parallel = [0.0312375,0.0525687,0.29925,1.26904,3.03684,5.59946]
# serial = [0.0157217,0.14666,0.665039,1.9696,3.92505,7.00606]


# # Set the figure size and font size
# plt.figure(figsize=(8, 6))
# plt.rcParams.update({'font.size': 12})

# # Plot the data with a solid line and marker
# plt.plot(threads, serial, '-o', label="serial", color='blue')
# plt.plot(threads, parallel, '-o', label="parallel", color='red')



# # Add axis labels and title
# plt.xlabel('Size of sqaure matrix')
# plt.ylabel('Execution Time (seconds)')
# plt.suptitle("matrix size vs. Execution Time")
# plt.title("")

# # Add gridlines and legend
# plt.grid(True)
# plt.legend(loc='upper right')

# # Save the plot to the 'graphs' folder with the name 'plot.png'
# plt.savefig('graphs/sizeLarge.png')