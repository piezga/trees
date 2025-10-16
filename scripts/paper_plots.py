import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from src.config import load_config
from src.functions import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Load config
config = load_config()

# SENM parameters
nx = config['senm']['nx']
ny = config['senm']['ny']
nu = config['senm']['nu']         
kernel = config['senm']['kernel']
NUM_REALIZATIONS = config['senm']['num_realizations']

# GRID
forest_grid_width = config['grid']['forest']['width']
forest_grid_height = config['grid']['forest']['height']
senm_grid_width = config['grid']['senm']['width']
senm_grid_height = config['grid']['senm']['height']

# Templates
senm_spatial_file_template = config['senm_templates']['spatial']
simulations_path = config['senm_templates']['path']
path_template = config['forests']['templates']['path_template']
census_template = config['forests']['templates']['census_template']
names_template = config['forests']['templates']['names_template']

# Other variables
forest = 'barro'
num_species = config['analysis']['num_species']
resolutions = [50, 5]  # Main plot: 50, Inset: 5
censuses = config['forests']['censuses']
num = 100  # Number of species for spectrum

verbose = True

# Setup plot
sns.set(style="whitegrid", font_scale=1.2, rc={"axes.labelweight": "bold"})
fig, ax = plt.subplots(figsize=(10, 6))
x = np.arange(1, num + 1)
colors = sns.color_palette("colorblind", n_colors=2)

# Function to compute both spectra at a given resolution
def compute_spectra(resolution):
    if verbose: print(f'Computing for resolution = {resolution} m')
    
    n_bins_x = int(forest_grid_width / resolution)
    n_bins_y = int(forest_grid_height / resolution)
    n_bins_x_senm = int(senm_grid_width / resolution)
    n_bins_y_senm = int(senm_grid_height / resolution)

    # === SENM spectrum ===
    senm_mean, senm_std = compute_mean_senm_spectrum(num, n_bins_x_senm, n_bins_y_senm)

    # === Forest spectra ===
    forest_spectra = []
    for census in censuses:
        path = path_template.format(forest=forest)
        os.makedirs(f"{path}plots", exist_ok=True)

        df, names = load_forest_data(forest, census, num)
        forest_spectrum = compute_forest_spectrum(df, names, n_bins_x, n_bins_y)
        forest_spectra.append(forest_spectrum)

    forest_array = np.array(forest_spectra)
    mean_forest_spectrum = np.mean(forest_array, axis=0)
    std_forest_spectrum = np.std(forest_array, axis=0)

    return senm_mean, mean_forest_spectrum, senm_std, std_forest_spectrum

# === MAIN PLOT: RESOLUTION 50 ===
res_50_senm, res_50_forest, res_50_senm_std, res_50_forest_std = compute_spectra(50)

ax.errorbar(
    x, res_50_senm, yerr=res_50_senm_std,
    fmt='o--', color=colors[0], label='SENM (50 m)', markersize=5
)
ax.errorbar(
    x, res_50_forest, yerr=res_50_forest_std,
    fmt='o-', color=colors[0], label='Forest (50 m)', markersize=5
)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Eigenvalue Rank', weight='bold')
ax.set_ylabel('Eigenvalue Magnitude', weight='bold')
ax.set_title('Figure 1.a - SENM vs Forest Spectra', weight='bold', fontsize=14)
ax.legend(frameon=True)
ax.grid(True, which="both", linestyle='--', alpha=0.6)

# === INSET PLOT: RESOLUTION 5 ===
res_5_senm, res_5_forest, res_5_senm_std, res_5_forest_std = compute_spectra(5)

axins = inset_axes(ax, width="45%", height="45%", loc='lower left', borderpad=2)

axins.errorbar(
    x, res_5_senm, yerr=res_5_senm_std,
    fmt='o--', color=colors[1], label='SENM (5 m)', markersize=3
)
axins.errorbar(
    x, res_5_forest, yerr=res_5_forest_std,
    fmt='o-', color=colors[1], label='Forest (5 m)', markersize=3
)
axins.set_xscale('log')
axins.set_yscale('log')
axins.set_xlim(1, num)
axins.set_ylim(
    min(min(res_5_senm - res_5_senm_std), min(res_5_forest - res_5_forest_std)),
    max(max(res_5_senm + res_5_senm_std), max(res_5_forest + res_5_forest_std))
)
axins.tick_params(labelsize=8)
axins.grid(True, which="both", linestyle='--', alpha=0.4)
axins.set_title('5 m', fontsize=10, pad=2)

# === Finalize and Save ===
plt.tight_layout()
#output_path = f"{path}plots/figure_1a.png"
#plt.savefig(output_path, dpi=300)
plt.show()

#print(f"Saved Figure 1.a to: {output_path}")

