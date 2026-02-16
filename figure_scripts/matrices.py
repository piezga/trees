import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import cmocean
import os
from sklearn.linear_model import LinearRegression

from src.config import load_config
from src.compute import (compute_spectra, compute_average_correlation, 
                         MarchenkoPastur, compute_stripped_correlation_matrix,
                         L_detect_communities)
from src.plotting import plot_correlation_matrices_comparison

config = load_config()
num_species = config['analysis']['num_species']
path_template = config['forests']['templates']['path_template']
names_template = config['forests']['templates']['names_template']

# Parameters
resolution = 20
forest = 'barro'
debug = False


print("Loading species abundance data...")
(senm_mean, senm_std, forest_spectra, 
 bins, senm_abundance, forest_abundances) = compute_spectra(
    resolution, num_species=num_species, calculate=False
)

# Create spatial bins matching species data
n_bins_x = bins[2]
n_bins_y = bins[3]
n_sites = n_bins_x * n_bins_y 

# Compute average forest abundance across censuses
forest_abundance_avg = np.mean(np.array(forest_abundances), axis=0)
print(f"  Species abundance shape: {forest_abundance_avg.shape}")
print(f"  (species={forest_abundance_avg.shape[0]}, sites={forest_abundance_avg.shape[1]})")


print("\nComputing correlation matrix...")
unfiltered_corr = np.corrcoef(forest_abundance_avg)
T = n_sites
filtered_corr = MarchenkoPastur(unfiltered_corr, 100, T)
(reordered_filtered_corr, CM, idx, 
 linkage_matrix_norm, cut_height) = L_detect_communities(filtered_corr, n_communities = 9,
                                                         return_linkage = True)
# Plot 
plot_correlation_matrices_comparison(unfiltered_corr,
                                     reordered_filtered_corr,
                                     data_type = 'forest')
