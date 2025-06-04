import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from functions import load_forest_data, load_senm_data
from variables import *



# Parameters for simulated data
NUM_REALIZATIONS = 50  # Number of simulated configurations
NUM_SPECIES = 1         # Analyze top N most abundant species (e.g., 1 for the most abundant)
DOMAIN_SIZE = 500       

# Define parameters for the analysis
min_box_size = 2        # Minimum box size 
max_box_size = 50     # Maximum box size
step_size = 1           # Step size for increasing box size
N = 10000              # Number of random boxes per size

# Initialize storage for aggregated results
all_box_sizes = np.arange(min_box_size, max_box_size, step_size)
all_variances = np.zeros((NUM_REALIZATIONS, len(all_box_sizes)))

# Loop over each realization
for realization in range(1, NUM_REALIZATIONS + 1):  # Realizations are 1-based
    # Load simulated data using the new function
    df = load_senm_data(nx, ny, nu, kernel, realization)
    
    # Find most abundant species
    top_species = df['species_id'].value_counts().head(NUM_SPECIES).index
    most_abundant_species = top_species[0]
    species_data = df[df['species_id'] == most_abundant_species][['x', 'y']].values
    
    print(f"Realization {realization}, Species ID: {most_abundant_species}")

    # Loop over box sizes
    for i, S in enumerate(all_box_sizes):
        counts = []
        
        # Throw N random boxes
        for _ in range(N):
            x0 = np.random.uniform(0, DOMAIN_SIZE - S)
            y0 = np.random.uniform(0, DOMAIN_SIZE - S)
            
            # Count points in the box
            in_box = ((species_data[:, 0] >= x0) & (species_data[:, 0] <= x0 + S)) & \
                     ((species_data[:, 1] >= y0) & (species_data[:, 1] <= y0 + S))
            counts.append(np.sum(in_box))
        
        # Store rescaled variance (σ²/S²)
        if len(counts) > 1:
            all_variances[realization-1, i] = np.var(counts) / np.mean(counts)

# Plot mean variance across realizations
mean_variances = np.mean(all_variances, axis=0)
std_variances = np.std(all_variances, axis=0)

plt.figure(figsize=(10, 6))
plt.plot(all_box_sizes, mean_variances, 'o-', label='Mean ± 1 STD')
plt.fill_between(all_box_sizes, 
                 mean_variances - std_variances, 
                 mean_variances + std_variances, 
                 alpha=0.2)
plt.xlabel('Box size (S)')
plt.ylabel('Variance / Mean')
plt.title(f'Density Variance vs Box Size (Averaged over {NUM_REALIZATIONS} Realizations)\n')
plt.grid(True)
plt.legend()
plt.show()




################################################
# forest

