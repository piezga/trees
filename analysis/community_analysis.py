import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns

from variables import *
from functions import load_forest_data, load_senm_data
from functions import get_top_species
from spectrum_functions import *


forest = 'barro'
censuses = [1,2,3,4,5,6,7,8]
num_species = 250
bin_numbers = range(20,100,2)

senm_communities = []
empirical_communities =  []

for n_bins_unrounded in bin_numbers:
    n_bins = int(n_bins_unrounded)
    print(f"\n{'='*50}\nAnalyzing {n_bins}x{n_bins} bins\n{'='*50}")

    boxlen = GRID_SIZE / n_bins
    boxlen_label = f"{int(round(boxlen))}m"
    
    print(f'Box length is {boxlen}')

    # Compute MP bounds
    lambda_min, lambda_max = marchenko_pastur_bounds(num_species,n_bins, n_bins)
    
    print('Max MP eigenvalue is ')
    print(lambda_max)
    print('Min MP eigenvalues is')
    print(lambda_min)

    # Compute SENM reference spectrum once per bin size
    senm_mean, senm_std = compute_mean_senm_spectrum(num_species,n_bins, n_bins)

    print('SENM spectrum preview:')
    print(senm_mean[:10])


    # Compute forest spectra and average
    forest_spectra = []

    for census in censuses:
        print(f"\nProcessing {forest} census {census}...")
        path = path_template.format(forest=forest)
        os.makedirs(f"{path}plots", exist_ok=True)

        df, names = load_forest_data(forest, census, num_species)
        #randomized_df = shuffle_labels(df)
        forest_spectrum = compute_forest_spectrum(df, names, n_bins, n_bins)
        #random_forest_spectrum = compute_forest_spectrum(randomized_df, names, n_bins, n_bins)

        forest_spectra.append(forest_spectrum)

    mean_forest_spectrum = np.mean(np.array(forest_spectra), axis = 0)

    print('Forest spectrum preview:')
    print(mean_forest_spectrum[:10])

    # Count eigenvalues above MP upper bound for SENM null model
    senm_count = np.sum(senm_mean > lambda_max)
    print(f'Number of SENM communities found: {senm_count}')
    senm_communities.append(senm_count)

    # Count eigenvalues above MP upper bound for empirical forest spectrum
    forest_count = np.sum(mean_forest_spectrum > lambda_max)
    print(f'Number of empirical communities found: {forest_count}')
    empirical_communities.append(forest_count)
        
# Convert to arrays

senm_communities = np.array(senm_communities)
print(f'SENM communities are {senm_communities}')

empirical_communities = np.array(empirical_communities)
print(f'Empirical communities are {empirical_communities}')

# Difference
difference_vector = senm_communities - empirical_communities


# Create the plot
plt.figure(figsize=(10, 6))

# Plot SENM null model results
plt.plot(bin_numbers, senm_communities, 'o-', label='SENM Null Model', color='blue', linewidth=2, markersize=6)

# Plot empirical forest results
plt.plot(bin_numbers, empirical_communities, 's-', label='Empirical Forest', color='red', linewidth=2, markersize=6)

# Add labels and title
plt.xlabel('Number of Bins', fontsize=12)
plt.ylabel('Number of Eigenvalues > λ_max', fontsize=12)
plt.title('Number of Significant Eigenvalues vs. Spatial Scale', fontsize=14)

# Add legend
plt.legend(fontsize=11)

# Add grid for better readability
plt.grid(True, alpha=0.3)

# Adjust layout to prevent label cutting
plt.tight_layout()

# Show the plot
plt.show()



# Difference plot
plt.figure(figsize=(10, 6))


# Plot difference
plt.plot(bin_numbers, difference_vector, 'o-',  color='blue', linewidth=2, markersize=6)


# Add labels and title
plt.xlabel('Number of Bins', fontsize=12)
plt.ylabel('Difference of Eigenvalues > λ_max', fontsize=12)
plt.title('Difference of Significant Eigenvalues ', fontsize=14)

# Add legend
plt.legend(fontsize=11)

# Add grid for better readability
plt.grid(True, alpha=0.3)

# Adjust layout to prevent label cutting
plt.tight_layout()

# Show the plot
plt.show()
