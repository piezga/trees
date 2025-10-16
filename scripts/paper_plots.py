import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from src.config import load_config
from src.functions import *

config = load_config()

# SENM parameters

nx = config['senm']['nx']
ny = config['senm']['ny']
nu = config['senm']['nu']         
kernel = config['senm']['kernel']
NUM_REALIZATIONS = config['senm']['num_realizations']

# GRID
forest_grid_width= config['grid']['forest']['width']
forest_grid_height= config['grid']['forest']['height']
senm_grid_width= config['grid']['senm']['width']
senm_grid_height= config['grid']['senm']['height']


# Name templates
senm_spatial_file_template = config['senm_templates']['spatial']
simulations_path = config['senm_templates']['path']
path_template = config['forests']['templates']['path_template']
census_template = config['forests']['templates']['census_template']
names_template = config['forests']['templates']['names_template']

# Other variables
forest = 'barro'
num_species = config['analysis']['num_species']
resolutions = [5,50]
censuses = config['forests']['censuses']
num = 100 # Number of species for plot 1


verbose = True

# Store spectra for plotting later
senm_spectra = []
forest_mean_spectra = []
resolution_labels = []

# Figure 1.a

for resolution in resolutions:

    if verbose: print(f'Resolution is {resolution} m')

    n_bins_x = int(forest_grid_width / resolution)
    n_bins_y = int(forest_grid_height / resolution)
    n_bins_x_senm = int(senm_grid_width/resolution)
    n_bins_y_senm = int(senm_grid_height/resolution)
    total_bins = n_bins_x * n_bins_y
    if verbose: print(f'Number of bins is {n_bins_x} x {n_bins_y} = {total_bins}')
    if verbose: print(f'Number of SENM bins is {n_bins_x_senm} x {n_bins_y_senm} = {n_bins_x_senm * n_bins_y_senm}')


    resolution_label = f'{int(round(resolution))}' 
    resolution_labels.append(f'{resolution_label} m')
    
    # Compute SENM spectrum
    senm_mean, senm_std = compute_mean_senm_spectrum(num, n_bins_x_senm,n_bins_y_senm)
    senm_spectra.append(senm_mean)

    # Compute forest spectra and mean spectrum
    forest_spectra = []

    for census in censuses:


        path = path_template.format(forest=forest)
        os.makedirs(f"{path}plots", exist_ok=True)

        df, names = load_forest_data(forest, census, num)

        forest_spectrum = compute_forest_spectrum(df, names, n_bins_x, n_bins_y)
        forest_spectra.append(forest_spectrum)

    mean_forest_spectrum = np.mean(np.array(forest_spectra), axis=0)
    forest_mean_spectra.append(mean_forest_spectrum)

# Plotting Figure 1.a

sns.set(style="whitegrid", font_scale=1.2, rc={"axes.labelweight": "bold"})

plt.figure(figsize=(10, 6))
x = np.arange(1, num + 1)
colors = sns.color_palette("colorblind", n_colors=len(resolutions))

for i in range(len(resolutions)):
    # SENM spectrum with dashed lines and dots
    plt.plot(x, senm_spectra[i], 'o--', color=colors[i], label=f'SENM ({resolution_labels[i]})', markersize=5)
    # Forest spectrum with solid lines and dots
    plt.plot(x, forest_mean_spectra[i], 'o-', color=colors[i], label=f'Forest ({resolution_labels[i]})', markersize=5)

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Eigenvalue Rank', weight='bold')
plt.ylabel('Eigenvalue Magnitude', weight='bold')
plt.title('Figure 1.a - SENM vs Forest Spectra at Two Resolutions', weight='bold', fontsize=14)
plt.legend(frameon=True)
plt.grid(True, which="both", linestyle='--', alpha=0.6)
plt.tight_layout()

# Save figure
output_path = f"{path}plots/figure_1a.png"
plt.savefig(output_path, dpi=300)
plt.show()

print(f"Saved Figure 1.a to: {output_path}")
