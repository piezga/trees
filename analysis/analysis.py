import numpy as np
import matplotlib.pyplot as plt
from variables import *
import os

# Do you want the graph?

graph = 0



# Move files to data folder
[os.system(f"mv ../{file}* {simulations_path}/ 2>/dev/null") 
    for file in ["AConf","Abund"]]
"""
# Read the file
data = np.loadtxt(simulations_path + spatial_file, dtype='int')
abundances = np.loadtxt(simulations_path + abundance_file, dtype='int')

# Process spatial data
x_coords = data[:, 0]  
y_coords = data[:, 1]
labels = data[:, 2]  

# Process abundance data
abundances = abundances[abundances[:, 1].argsort()[::-1]]  # Sort by second column

# Determine plot dimensions
max_x, max_y = x_coords.max(), y_coords.max()

# Initialize empty matrix (fill with '' or np.nan for numeric data)
plot = np.full((max_x + 1, max_y + 1), 0)  # +1 for 0-based indexing
most_abundant_plot = np.full((max_x + 1, max_y + 1), 0)  
# All species matrix
for x, y, label in zip(x_coords, y_coords, labels):
    plot[x, y] = label

# Most abundant species matrix
for x, y, label in zip(x_coords, y_coords, labels):
    if label == abundances[0,0]:
        most_abundant_plot[x, y] = 1


# Quantities

tree_number = nx*ny
abundance_ratio = abundances[0,1]/tree_number

# Information printing

print("Analyzed files: ")
print(f"{spatial_file} \n {abundance_file}")
print(f"First species abundance : {abundance_ratio}")


#graph

if graph: 
    plt.pcolormesh(most_abundant_plot, cmap='Greens')
    plt.colorbar()
    plt.title(f'nu = {nu}, K = {kernel}')
    plt.show()

"""
