import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns
from scipy.special import rel_entr
from scipy.interpolate import interp1d

from variables import *
from functions import load_forest_data, load_senm_data
from functions import get_top_species
from spectrum_functions import *


forest = 'barro'
censuses = [1, 2, 3, 4, 5, 6, 7, 8]  # Your selected censuses to analyze
bin_numbers = [12] 
num_species = 100




divergence_over_bins = []
square_differences = []

for n_bins_unrounded in bin_numbers:
    n_bins = int(n_bins_unrounded)
    print(f"\n{'='*50}\nAnalyzing {n_bins}x{n_bins} bins\n{'='*50}")
    
    boxlen = GRID_SIZE / n_bins
    boxlen_label = f"{int(round(boxlen))}m"

    # Compute SENM reference spectrum once per bin size
    senm_mean, senm_std = compute_mean_senm_spectrum(num_species,n_bins, n_bins)
   
    # Prepare storage for all census data
    all_spectra = [senm_mean]  # SENM spectrum first
    all_labels = ['SENM']      # Label for SENM
    all_forest_spectra = []
    all_random_spectra = [] 
    
    for census in censuses:
        print(f"\nProcessing {forest} census {census}...")
        path = path_template.format(forest=forest)
        os.makedirs(f"{path}plots", exist_ok=True)
        
        # Compute forest spectrum_to_distribution
        df, names = load_forest_data(forest, census, num_species)
        randomized_df = shuffle_labels(df)
        forest_spectrum = compute_forest_spectrum(df, names, n_bins, n_bins)
        random_forest_spectrum = compute_forest_spectrum(randomized_df, names, n_bins, n_bins)
        
        all_spectra.append(forest_spectrum)  # Add forest spectra to the list
        all_forest_spectra.append(forest_spectrum)
        all_random_spectra.append(random_forest_spectrum)
        all_labels.append(str(census))       # Add census label
    
    mean_forest_spectrum = np.mean(np.array(all_forest_spectra), axis=0)
    
    square_differences.append(np.sum(square_difference(senm_mean[:20], mean_forest_spectrum[:20])))



    # Plot the combined spectra for this bin size (SENM + Forest Spectra)
    plot_combined_spectra(num_species,all_spectra, all_labels, f"Correlation matrix spectrum at L = {boxlen_label}", 
                          f"{path}plots/combined_spectrum_{n_bins}.png", senm_std=senm_std
                          , n_bins_x= n_bins, n_bins_y = n_bins)
    

    """
    # Calculate the KL divergence
    distribution_bins = 20
    forest_distribution, senm_distribution = spectra_to_distributions(forest_spectrum
                                                                      , mean_forest_spectrum, n_bins=distribution_bins)
    low_filter = 4
    kl_divergence = KL_divergence(forest_distribution[low_filter:], senm_distribution[low_filter:])
    
    plot_distributions(forest_distribution, senm_distribution, kl_divergence, 
                       title_prefix=f"Size {boxlen:.2f}, Filter {low_filter} ")

    # Store the divergence for this bin size
    divergence_over_bins.append(kl_divergence)
    
plt.plot(square_differences)
plt.show()


# Print and plot KL divergence over bin numbers
print(f"Divergence over bins is: {divergence_over_bins}")
plot_KL(divergence_over_bins, bin_numbers, GRID_SIZE, prefix=f"Filter {low_filter}") 
plt.show()
"""
