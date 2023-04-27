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
# parallel = [1.87535,1.14586,1.10123,1.10648,1.10977,1.13271,1.13092,1.14381,1.16512,1.18542,1.19614,1.23085]
# serial = [1.83777,1.82945,1.82719,1.84386,1.84174,1.83843,1.83992,1.84405,1.84712,1.85111,1.8325,1.86247]

# 2000 by 2000 on sunlab
# parallel = [0.938545,0.571836,0.50912,0.47427,0.459708,0.451544,0.444284,0.444054,0.443602,0.445446,0.443874,0.442856]
# serial = [0.912867,0.914606,0.935372,0.92244,0.924902,0.922596,0.92386,0.923549,0.921274,0.920542,0.92076,0.922489]

# 3000 by 3000 on sublab
# parallel = [3.70236,2.72211,2.58541,2.53131,2.48129,2.48064,2.47313,2.4952,2.49963,2.50859,2.51234,2.51945]
# serial = [3.68416,3.6952,3.67579,3.68511,3.68938,3.68659,3.68568,3.69087,3.69132,3.7004,3.69762,3.69059]

# 4000 by 4000 on sublab
# parallel = [9.3059,7.29983,6.98651,6.92029,6.74694,6.78461,6.78477,6.83987,6.8589,6.87411,6.89067,6.90593]
# serial = [9.01141,9.07917,9.14834,9.17323,9.21449,9.19497,9.17806,9.19855,9.20864,9.19725,9.20875,9.23712]

# 3000 by 3000 on sunlab kick out of tbb when less than 200
# parallel = [3.6861,2.71631,2.58569,2.52238,2.47411,2.48103,2.47309,2.4919,2.51586,2.50159,2.5093,2.50831]
# serial = [3.69725,3.68438,3.68933,3.6873,3.70021,3.6943,3.70485,3.68997,3.71045,3.70537,3.69383,3.69587]
# Set the figure size and font size
# plt.figure(figsize=(8, 6))
# plt.rcParams.update({'font.size': 12})

# # Plot the data with a solid line and marker
# plt.plot(threads, serial, '-o', label="serial", color='blue')
# plt.plot(threads, parallel, '-o', label="parallel", color='red')



# # Add axis labels and title
# plt.xlabel('Threads')
# plt.ylabel('Execution Time (seconds)')
# plt.suptitle("Threads vs. Execution Time")
# plt.title("3000 by 3000 dense matrix on SunLab, don't tbb if less than 200 on sub")

# # Add gridlines and legend
# plt.grid(True)
# plt.legend(loc='upper right')

# # Save the plot to the 'graphs' folder with the name 'plot.png'
# plt.savefig('graphs/threadSun3000kick200.png')

#matrix size
# threads = [1000,1500,2000,2500,3000,3500,4000,4500]

# parallel = [0.0312375,0.0525687,0.29925,1.26904,3.03684,5.59946]
# serial = [0.0157217,0.14666,0.665039,1.9696,3.92505,7.00606]

#sunlab
threads = [1000,1500,2000]
parallel = [0.0274525,0.0819047,0.437466]
serial = [0.0708806,0.288978,0.919206]


# Set the figure size and font size
plt.figure(figsize=(8, 6))
plt.rcParams.update({'font.size': 12})

# Plot the data with a solid line and marker
plt.plot(threads, serial, '-o', label="serial", color='blue')
plt.plot(threads, parallel, '-o', label="parallel", color='red')



# Add axis labels and title
plt.xlabel('Size of sqaure matrix')
plt.ylabel('Execution Time (seconds)')
plt.suptitle("matrix size vs. Execution Time on sunlab smaller matrix sizes")
plt.title("")

# Add gridlines and legend
plt.grid(True)
plt.legend(loc='upper right')

# Save the plot to the 'graphs' folder with the name 'plot.png'
plt.savefig('graphs/sizeSunLabSmaller.png')