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
num_species = 200
bin_numbers = range(10,30,20)
num_replicas = 20
num_removals = 5


for n_bins_unrounded in bin_numbers:
    
    n_bins = int(n_bins_unrounded)
    print(f"\n{'='*50}\nAnalyzing {n_bins}x{n_bins} bins\n{'='*50}") 

    boxlen = GRID_SIZE / n_bins
    boxlen_label = f"{int(round(boxlen))}m"
    print(f'Box length is {boxlen}')
    
    senm_mean, senm_std = senm_bootstrap(num_species, n_bins, n_bins, num_removals, 
                                         num_replicas)

    print('SENM spectrum preview:')
    print(senm_mean[:10])

   
    forest_mean, forest_std = forest_bootstrap(forest, censuses, num_species, n_bins,
                                              n_bins, num_removals, num_replicas)



    # Create the plot
    plt.figure(figsize=(12, 8))
    
    # Create eigenvalue index
    eigenvalues_idx = np.arange(1, len(senm_mean) + 1)
    
    # Plot SENM spectrum with error bars
    plt.errorbar(eigenvalues_idx, senm_mean, yerr=senm_std, 
                 fmt='o-', markersize=4, capsize=3, capthick=1, elinewidth=1,
                 label='SENM Null Model', color='blue', alpha=0.8)
    
    # Plot forest spectrum with error bars
    plt.errorbar(eigenvalues_idx, forest_mean, yerr=forest_std,
                 fmt='s-', markersize=4, capsize=3, capthick=1, elinewidth=1,
                 label='Empirical Forest', color='red', alpha=0.8)
    """   
    # Add MP bounds as reference lines
    plt.axhline(y=lambda_max, color='green', linestyle='--', alpha=0.7, 
                label=f'MP Upper Bound ({lambda_max:.2f})')
    plt.axhline(y=lambda_min, color='purple', linestyle='--', alpha=0.7, 
                label=f'MP Lower Bound ({lambda_min:.2f})')
    """ 

    # Customize the plot
    plt.xlabel('Eigenvalue Rank', fontsize=12)
    plt.ylabel('Eigenvalue Magnitude', fontsize=12)
    plt.title(f'Eigenvalue Spectra Comparison\n{n_bins}x{n_bins} bins ({boxlen_label} resolution)', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.yscale('log')  
    plt.xscale('log')   
    
    # Save the plo
    plot_dir = f'../data/{forest}/plots/'
    filename = f"{plot_dir}spectrum_{n_bins}_bins_choose_{num_species}_remove_{num_removals}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {filename}")
    
    plt.close()  # Close the figure to free memory
