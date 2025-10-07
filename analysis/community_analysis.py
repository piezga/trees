import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib
matplotlib.use('TkAgg')  # Switch from Qt to Tkinter backend
import matplotlib.pyplot as plt

from variables import *
from functions import load_forest_data, load_senm_data
from functions import get_top_species
from spectrum_functions import *


forests = ['wanang']
censuses_barro = [1,2,3,4,5,6,7,8]
censuses_wanang = [1,4]
num_species = [50,100,200,250]
resolutions = range(3,20,1)

# Create the plot
plt.figure(figsize=(10, 6))

for forest in forests:  # ADDED: Outer loop over forests

    if forest == 'barro':
        censuses = censuses_barro
    elif forest == 'wanang':
        censuses = censuses_wanang

    print(f"Censuses are {censuses}")

    for num in num_species:

        senm_communities = []
        empirical_communities =  []
        square_diffs = []

        #  FIX: Reset NN distances for this species count
        census_nn_distances = []

        #  COMPUTE NN DISTANCES ONCE PER num_species (not per resolution)
        for census in censuses:
            print(f"Calculating NN distance for census {census}")
            df, names = load_forest_data(forest, census, num)


            coords = df[['x', 'y']].values
            if len(coords) > 1:
                from scipy.spatial import KDTree
                tree = KDTree(coords)
                distances, indices = tree.query(coords, k=2)
                min_distances = distances[:, 1]
                census_nn = np.mean(min_distances)
                census_nn_distances.append(census_nn)

        for resolution in resolutions:


            n_bins_unrounded = GRID_SIZE/resolution
            n_bins = int(n_bins_unrounded)
            boxlen = resolution
            boxlen_label = f"{int(round(resolution))}m"
            

            # Compute MP bounds
            lambda_min, lambda_max = marchenko_pastur_bounds(num,n_bins, n_bins)
            
            # Compute SENM reference spectrum once per bin size
            senm_mean, senm_std = compute_mean_senm_spectrum(num,n_bins, n_bins)

            forest_spectra = []

            for census in censuses:
                path = path_template.format(forest=forest)
                os.makedirs(f"{path}plots", exist_ok=True)

                df, names = load_forest_data(forest, census, num)

                forest_spectrum = compute_forest_spectrum(df, names, n_bins, n_bins)
                forest_spectra.append(forest_spectrum)

            mean_forest_spectrum = np.mean(np.array(forest_spectra), axis = 0)

            # Count eigenvalues above MP upper bound for SENM null model
            senm_count = np.sum(senm_mean > lambda_max)
            senm_communities.append(senm_count)

            # Count eigenvalues above MP upper bound for empirical forest spectrum
            forest_count = np.sum(mean_forest_spectrum > lambda_max)
            empirical_communities.append(forest_count)

            # Calculate square difference 
            square_diff = square_diff_above_MP(senm_mean, mean_forest_spectrum, lambda_max)
            square_diffs.append(square_diff)

        # CALCULATE mean NN distance after loop
        mean_nn_distance = np.mean(census_nn_distances) if census_nn_distances else 0
        print(f"Mean nearest neighbor distance for {forest} with {num} species: {mean_nn_distance:.6f}m")
                
        # Convert to arrays
        senm_communities = np.array(senm_communities) 

        empirical_communities = np.array(empirical_communities)

        # Difference
        difference_vector = np.abs(senm_communities - empirical_communities)

        # Plot square difference
        plt.plot(resolutions, square_diffs, 'o-', label = f'{forest} - {num} Species', linewidth=2, markersize=6)


# Final plot setup
plt.xlabel('Length', fontsize=12)
plt.ylabel('Square difference  > λ_max', fontsize=12)
plt.title('Square difference of Significant Eigenvalues ', fontsize=14)
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show(block=True)
plt.close()

