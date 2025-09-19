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
censuses = [1,2,3,4]
num_species = 200
resolutions = range(5,50,5)
num_replicas = 5
num_removals = 150


for resolution in resolutions:
   
    print(f'Box length is {resolution}')

    n_bins_unrounded = GRID_SIZE/resolution 
    n_bins = int(n_bins_unrounded)
    print(f"\n{'='*50}\nAnalyzing {n_bins}x{n_bins} bins\n{'='*50}") 

    resolution_label = f"{int(round(resolution))}m"
    
    senm_mean, senm_std = senm_bootstrap(num_species, n_bins, n_bins, num_removals, 
                                         num_replicas)

    print('SENM spectrum preview:')
    print(senm_mean[:10])

   
    forest_mean, forest_std = forest_bootstrap(forest, censuses, num_species, n_bins,
                                              n_bins, num_removals, num_replicas)

        
    # Create the plot with nicer colors
    plt.figure(figsize=(10, 6))

    # Create eigenvalue index
    eigenvalues_idx = np.arange(1, len(senm_mean) + 1)

    # Use better colors from seaborn
    colors = sns.color_palette("Set1", 3)

    # Plot SENM spectrum with error bars
    plt.errorbar(eigenvalues_idx, senm_mean, yerr=senm_std,
         fmt='o-', markersize=4, capsize=3, capthick=1.5, elinewidth=1.5,
         label='SENM Null Model', color=colors[1], alpha=0.8)

    # Plot forest spectrum with error bars
    plt.errorbar(eigenvalues_idx, forest_mean, yerr=forest_std,
         fmt='s-', markersize=4, capsize=3, capthick=1.5, elinewidth=1.5,
         label='Empirical Forest', color=colors[2], alpha=0.8)

    # Add MP bounds
    lambda_min, lambda_max = marchenko_pastur_bounds(num_species-num_removals, n_bins, n_bins)
    
    # Add MP bounds as reference lines
    plt.axhline(y=lambda_max, color='red', linestyle='--', alpha=0.7, 
                label=f'MP Upper Bound ({lambda_max:.2f})')
    plt.axhline(y=lambda_min, color='red', linestyle='--', alpha=0.7, 
                label=f'MP Lower Bound ({lambda_min:.2f})')

    # Customize the plot
    plt.xlabel('Eigenvalue Rank', fontsize=12)
    plt.ylabel('Eigenvalue Magnitude', fontsize=12)
    plt.title(f'Eigenvalue Spectra Comparison\n{n_bins}x{n_bins} bins ({resolution_label} resolution)', fontsize=14)
    plt.legend(fontsize=11, frameon=True, fancybox=True, shadow=True)
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.yscale('log')
    plt.xscale('log')

    # Add a nice background
    plt.gca().set_facecolor('#f8f9fa')

    # Save the plot
    plot_dir = f'../data/{forest}/plots/'
    os.makedirs(plot_dir, exist_ok=True)
    filename = f"{plot_dir}spectrum_{resolution}_meters_choose_{num_species}_remove_{num_removals}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Plot saved to: {filename}")

    plt.close()



