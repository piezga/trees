import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib
matplotlib.use('TkAgg')  # Switch from Qt to Tkinter backend
import matplotlib.pyplot as plt

from src.config import load_config
config = load_config()
from src.functions import *

debug = True

GRID_SIZE = config['grid']['size']
forests = ['barro']
num_species = config['analysis']['num_species']
resolutions = range(3, 5, 1)

# Create the plot
plt.figure(figsize=(10, 6))

for forest in forests:  
    censuses = config['forests']['censuses'][f'{forest}']
    
    if debug:
        print(f"Censuses are {censuses}")

    for num in num_species:

        senm_communities = []
        empirical_communities = []
        square_diffs = []

        # Initialize dict to track individual square diffs per census
        individual_square_differences = {census: [] for census in censuses}

        for resolution in resolutions:

            n_bins_unrounded = GRID_SIZE / resolution
            n_bins = int(n_bins_unrounded)
            boxlen = resolution
            boxlen_label = f"{int(round(resolution))}m"

            # Compute MP bounds
            lambda_min, lambda_max = marchenko_pastur_bounds(num, n_bins, n_bins)
            
            # Compute SENM reference spectrum once per bin size
            senm_mean, senm_std = compute_mean_senm_spectrum(num, n_bins, n_bins)

            forest_spectra = []

            for census in censuses:
                path = path_template.format(forest=forest)
                os.makedirs(f"{path}plots", exist_ok=True)

                df, names = load_forest_data(forest, census, num)

                forest_spectrum = compute_forest_spectrum(df, names, n_bins, n_bins)
                forest_spectra.append(forest_spectrum)

                # Compute square difference for individual census
                individual_square_diff = square_diff_above_MP(
                    senm_mean, forest_spectrum, lambda_max
                )

                # Append to corresponding list in the dict
                individual_square_differences[census].append(individual_square_diff)

            # Mean empirical spectrum over censuses
            mean_forest_spectrum = np.mean(np.array(forest_spectra), axis=0)

            # Count eigenvalues above MP upper bound for SENM null model
            senm_count = np.sum(senm_mean > lambda_max)
            senm_communities.append(senm_count)

            # Count eigenvalues above MP upper bound for empirical forest spectrum
            forest_count = np.sum(mean_forest_spectrum > lambda_max)
            empirical_communities.append(forest_count)

            # Calculate square difference between mean SENM and mean forest
            square_diff = square_diff_above_MP(senm_mean, mean_forest_spectrum, lambda_max)
            square_diffs.append(square_diff)

        # Plot square difference per census (line across resolutions)
        for census in censuses:
            plt.plot(
                list(resolutions),
                individual_square_differences[census],
                'o--',
                label=f"{forest}-{num}-{census}"
            )

# Final plot setup
plt.xlabel('Length', fontsize=12)
plt.ylabel('Square difference > λ_max', fontsize=12)
plt.title('Square difference of Significant Eigenvalues', fontsize=14)
plt.legend(fontsize=9, loc='best')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show(block=True)
plt.close()

