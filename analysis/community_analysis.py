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
num_species = [50,100,200]
resolutions = range(5,18,1)

# Create the plot
plt.figure(figsize=(10, 6))

for forest in forests:  # ADDED: Outer loop over forests

    if forest == 'barro':
        censuses = censuses_barro
    elif forest == 'wanang':
        censuses = censuses_wanang

    print(f"Censuses are {censuses}")


    for num in num_species:
        # This is calculated by taking T > 10N
        max_MP_boxlen = 500/(3*np.sqrt(num))
        print(f"Max trustworthy resolution: {max_MP_boxlen}")
        senm_communities = []
        empirical_communities =  []
        square_diffs = []
        for resolution in resolutions:

            print(f'Box length is {resolution}')

            n_bins_unrounded = GRID_SIZE/resolution

            n_bins = int(n_bins_unrounded)
            print(f"\n{'='*50}\nAnalyzing {n_bins}x{n_bins} bins\n{'='*50}")

            boxlen = resolution
            boxlen_label = f"{int(round(boxlen))}m"
            
            print(f'Box length is {boxlen}')

            # Compute MP bounds
            lambda_min, lambda_max = marchenko_pastur_bounds(num,n_bins, n_bins)
            
            print('Max MP eigenvalue is ')
            print(lambda_max)
            print('Min MP eigenvalues is')
            print(lambda_min)

            # Compute SENM reference spectrum once per bin size
            senm_mean, senm_std = compute_mean_senm_spectrum(num,n_bins, n_bins)

            print('SENM spectrum preview:')
            print(senm_mean[:10])


            # Compute forest spectra and average
            forest_spectra = []

            for census in censuses:
                print(f"\nProcessing {forest} census {census}...")
                path = path_template.format(forest=forest)
                os.makedirs(f"{path}plots", exist_ok=True)

                df, names = load_forest_data(forest, census, num)
                #randomized_df = shuffle_labels(df)
                forest_spectrum = compute_forest_spectrum(df, names, n_bins, n_bins)
                #random_forest_spectrum = compute_forest_spectrum(randomized_df, names, n_bins, n_bins)

                forest_spectra.append(forest_spectrum)

            mean_forest_spectrum = np.mean(np.array(forest_spectra), axis = 0)

            print('Forest spectrum preview:')
            print(mean_forest_spectrum[:10])

            # Count eigenvalues above MP upper bound for SENM null model
            senm_count = np.sum(senm_mean > lambda_max)
            print(f'\nNumber of SENM communities found: {senm_count}')
            senm_communities.append(senm_count)

            # Count eigenvalues above MP upper bound for empirical forest spectrum
            forest_count = np.sum(mean_forest_spectrum > lambda_max)
            print(f'\nNumber of empirical communities found: {forest_count}')
            empirical_communities.append(forest_count)

            # Calculate square difference 
            square_diff = square_diff_above_MP(senm_mean, mean_forest_spectrum, lambda_max)
            square_diffs.append(square_diff)
            print(f"\n Square difference is {square_diff}")

                
        # Convert to arrays

        senm_communities = np.array(senm_communities) 
        print(f'\nSENM communities are {senm_communities}')

        empirical_communities = np.array(empirical_communities)
        print(f'\nEmpirical communities are {empirical_communities}')

        # Difference
        difference_vector = np.abs(senm_communities - empirical_communities)
        """
        # Plot SENM null model results
        plt.plot(resolutions, senm_communities, 'o-', label='SENM Null Model', color='blue', linewidth=2, markersize=6)

        # Plot empirical forest results
        plt.plot(resolutions, empirical_communities, 's-', label='Empirical Forest', color='red', linewidth=2, markersize=6)

        # Add labels and title
        plt.xlabel('Length', fontsize=12)
        plt.ylabel('Number of Eigenvalues > λ_max', fontsize=12)
        plt.title('Number of Significant Eigenvalues vs. Spatial Scale', fontsize=14)

        # Add legend
        plt.legend(fontsize=11)

        # Add grid for better readability
        plt.grid(True, alpha=0.3)

        # Adjust layout to prevent label cutting
        plt.tight_layout()


        # Difference plot
        plt.figure(figsize=(10, 6))


        # Plot difference
        plt.plot(resolutions, difference_vector, 'o-',  color='blue', linewidth=2, markersize=6)


        # Add labels and title
        plt.xlabel('Length', fontsize=12)
        plt.ylabel('Difference of Eigenvalues > λ_max', fontsize=12)
        plt.title('Difference of Significant Eigenvalues ', fontsize=14)

        # Add legend
        plt.legend(fontsize=11)

        # Add grid for better readability
        plt.grid(True, alpha=0.3)

        # Adjust layout to prevent label cutting
        plt.tight_layout()
        """


        # Plot square difference
        plt.plot(resolutions, square_diffs, 'o-', label = f'{forest} - {num} Species', 
        linewidth=2, markersize=6)  # MODIFIED: Added forest name to label
        #plt.axvline(x=max_MP_boxlen, color=plt.gca().lines[-1].get_color(), 
        #linestyle='--', alpha=0.7)

# Add labels and title
plt.xlabel('Length', fontsize=12)
plt.ylabel('Square difference  > λ_max', fontsize=12)
plt.title('Square difference of Significant Eigenvalues ', fontsize=14)

# Add legend
plt.legend(fontsize=11)

# Add grid for better readability
plt.grid(True, alpha=0.3)

# Adjust layout to prevent label cutting
plt.tight_layout()

plt.show(block= True)

plt.close()
