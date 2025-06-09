import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from functions import load_forest_data, load_senm_data
from variables import *

# Parameters
forest = 'barro'
census = 5
species = 200

# Load forest data
data  = load_forest_data(forest, census, species)  


# Find the most abundant species
species_counts = Counter(data['name'])
most_abundant_species = species_counts.most_common(1)[0][0]
most_abundant_species = 'hybapr'
mean_NN_dist = 1
print(f"Most abundant species: {most_abundant_species}")

# Filter data for the most abundant species
species_data = data[data['name'] == most_abundant_species][['x', 'y']].values

if most_abundant_species == 'all':
    species_data = data[['x', 'y']].values


# Define parameters
min_box_size = 1  # Minimum box size (adjust as needed)
max_box_size = 320  # Maximum box size
step_size = 10 # Step size for increasing box size
N = 100000  # Number of random boxes to throw for each size

# Initialize lists to store results
box_sizes = []
means = []
variances = []
aggregate_counts = []

# Main loop over box sizes
for S in range(min_box_size, int(max_box_size), step_size):
    S = S/mean_NN_dist
    counts = []
    
    
    # Throw N random boxes
    for _ in range(N):
        # Randomly position the box within the plot
        x0 = np.random.uniform(0, data['x'].max() - S)
        y0 = np.random.uniform(0, data['y'].max() - S)
        
        # Count points of the species within the box
        in_box = ((species_data[:, 0] >= x0) & (species_data[:, 0] <= x0 + S)) & \
                 ((species_data[:, 1] >= y0) & (species_data[:, 1] <= y0 + S))
        counts.append(np.sum(in_box))
        aggregate_counts.append(np.sum(in_box))
    # Calculate and store variance
    if len(counts) > 1:
        variance = np.var(counts)
        mean = np.mean(counts)
        rescaled_variance = variance/mean
        box_sizes.append(S)
        variances.append(rescaled_variance)
        print(f"Box size: {S}, Rescaled Variance : {rescaled_variance:.2f}")

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(box_sizes, variances, 'o-')
plt.xlabel('Box size (S)')
plt.ylabel('Variance over mean')
plt.title(f'Density Variance vs Box Size for {most_abundant_species} - census {census}')
plt.grid(True)
plt.show()

