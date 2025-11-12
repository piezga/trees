import numpy as np
import os
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from src.config import load_config
from src.functions import *
from src.plotting import *
from src.computing_functions import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib
import matplotlib.colors as mcolors
import cmocean
import networkx as nx
from scipy.linalg import expm
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
import matplotlib.patches as patches

matplotlib.use('TkAgg')  # Switch from Qt to Tkinter backend


# === Load config ===
config = load_config()

# === Parameters ===
Nx = config['senm']['nx']
Ny = config['senm']['ny']
nu = config['senm']['nu']
kernel = config['senm']['kernel']
NUM_REALIZATIONS = config['senm']['num_realizations']

forest_grid_width = config['grid']['forest']['width']
forest_grid_height = config['grid']['forest']['height']
senm_grid_width = config['grid']['senm']['width']
senm_grid_height = config['grid']['senm']['height']

senm_spatial_file_template = config['senm_templates']['spatial']
simulations_path = config['senm_templates']['path']
path_template = config['forests']['templates']['path_template']
census_template = config['forests']['templates']['census_template']
names_template = config['forests']['templates']['names_template']
forest = 'barro'
num_species = config['analysis']['num_species']
species_array = config['analysis']['species_array']
censuses = config['forests']['censuses']

# Some parameters
filter_resolution = 32

# Argument parser 

def parse_args():
    parser = argparse.ArgumentParser(description='Run paper plot analyses')
    """ 
    parser.add_argument(
        '--task',
        type=str,
        default='all',
        help='Task(s) to run: spectra, resolution-sweep, correlation, communities, or all (comma-separated)'
    )
    """
    parser.add_argument('-c', '--calculate', action='store_true', help='Calculate instead of loading cached')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print verbose output')
    return parser.parse_args()

args = parse_args()
calculate = args.calculate 
verbose = args.verbose
print(f'Calculate set to {calculate}')


# === Setup style ===
plt.rcParams.update({
    "font.size": 20,
    "axes.labelsize": 24,
    "axes.titlesize": 24,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 18,
    "lines.linewidth": 2,
    "lines.markersize": 6,
    "axes.linewidth": 1.2
})

# === Create subplot mosaic ===
fig, ax_dict = plt.subplot_mosaic("AB", figsize=(14, 6))

# ============================================
# === Panel A: Dual Spectra with Inset ===
# ============================================
ax_a = ax_dict["A"]

# Load or compute data for Panel A
if calculate:
    res_50_senm, res_50_senm_std, res_50_forest_list, _, _, _ = compute_spectra(50)
    np.save(f'quantities/barro_senm_std_{num_species}_50.npy', res_50_senm_std)
else: 
    res_50_senm = np.load(f'quantities/{forest}_senm_spectrum_{num_species}_50.npy')
    res_50_senm_std = np.load(f'quantities/{forest}_senm_std_{num_species}_50.npy')
    res_50_forest_list = np.load(f'quantities/{forest}_forest_spectra_{num_species}_50.npy')

if calculate:
    res_5_senm, res_5_senm_std, res_5_forest_list, _, _, _ = compute_spectra(5)
    np.save(f'quantities/barro_senm_std_{num_species}_5.npy', res_5_senm_std)
else: 
    res_5_senm = np.load(f'quantities/{forest}_senm_spectrum_{num_species}_5.npy')
    res_5_senm_std = np.load(f'quantities/{forest}_senm_std_{num_species}_5.npy')
    res_5_forest_list = np.load(f'quantities/{forest}_forest_spectra_{num_species}_5.npy')

# Plot Panel A
plot_dual_spectra_with_inset(
    res_50_senm, res_50_senm_std, res_50_forest_list,
    res_5_senm, res_5_senm_std, res_5_forest_list,
    num_species=num_species,
    ax=ax_a,
    show=False
)

# ============================
# === Panel B: Size Effect ===
# ============================
ax_b = ax_dict["B"]
resolutions = list(range(4, 34, 1))

# Store results for Panel B and supplementary plots
results_for_panel_b = {}  # Will store: num_species -> list of comm_diffs
results_for_supplementary = []  # Will store all data needed for extra plots

for idx, num_species in enumerate(species_array):
    if verbose:
        print(f'Num species: {num_species} out of {species_array}')
    
    comm_diffs = []
    
    for resolution in resolutions:
        # Load or compute spectra
        if calculate:
            senm_spectrum, _, forest_spectra, bins, _, _ = compute_spectra(resolution)
            np.save(f'quantities/{forest}_senm_spectrum_{num_species}_{resolution}.npy', senm_spectrum)
            np.save(f'quantities/{forest}_forest_spectra_{num_species}_{resolution}.npy', forest_spectra)
            np.save(f'quantities/{forest}_bins_{num_species}_{resolution}.npy', bins)
            np.savetxt(f'quantities/{forest}_senm_spectrum_{num_species}_{resolution}.txt', senm_spectrum)
            np.savetxt(f'quantities/{forest}_forest_spectra_{num_species}_{resolution}.txt', forest_spectra)
            np.savetxt(f'quantities/{forest}_bins_{num_species}_{resolution}.txt', bins)
        else:
            senm_spectrum = np.load(f'quantities/{forest}_senm_spectrum_{num_species}_{resolution}.npy')
            forest_spectra = np.load(f'quantities/{forest}_forest_spectra_{num_species}_{resolution}.npy')
            bins = np.load(f'quantities/{forest}_bins_{num_species}_{resolution}.npy')
        
        if verbose:
            print(f'\nBins array is: {bins}')
            print(f'Resolution is {resolution}')
        
        mean_forest_spectrum = np.mean(np.array(forest_spectra), axis=0)
        
        if verbose:
            print(f'SENM spectrum shape: {np.shape(senm_spectrum)}')
            print(f'Forest spectrum shape: {np.shape(mean_forest_spectrum)}')
        
        # Calculate MP bounds
        _, lambda_max_forest = marchenko_pastur_bounds(num_species, bins[2], bins[3])
        _, lambda_max_senm = marchenko_pastur_bounds(num_species, bins[0], bins[1])
        
        if verbose:
            print(f'Max forest eigenvalue is {lambda_max_forest}')
            print(f'Max senm eigenvalue is {lambda_max_senm}')
        
        if calculate:
            np.savetxt(f'quantities/{forest}_forest_MPmax_{num_species}_{resolution}.txt', [lambda_max_forest])
            np.savetxt(f'quantities/{forest}_senm_MPmax_{num_species}_{resolution}.txt', [lambda_max_senm])
        
        # Calculate differences
        square_diff, comm_diff = square_diff_above_MP(senm_spectrum, mean_forest_spectrum, 
                                                      lambda_max_senm, lambda_max_forest)
        forcom = number_of_communities(mean_forest_spectrum, lambda_max_forest)
        senmcom = number_of_communities(senm_spectrum, lambda_max_senm)
        
        comm_diffs.append(comm_diff)
        
        # Store data for supplementary plots
        results_for_supplementary.append({
            'num_species': num_species,
            'resolution': resolution,
            'forest_spectrum': mean_forest_spectrum,
            'senm_spectrum': senm_spectrum,
            'lambda_max_forest': lambda_max_forest,
            'lambda_max_senm': lambda_max_senm,
            'bins': bins,
            'forest_communities': forcom,
            'senm_communities': senmcom
        })
    
    # Store for Panel B
    results_for_panel_b[num_species] = comm_diffs

# Plot Panel B
plot_size_effect_panel(
    resolutions,
    results_for_panel_b,
    species_array,
    ax=ax_b,
    show=False
)

# === Finalize main figure ===
#plt.tight_layout()
plt.show()

# ============================================
# === Supplementary Plots ===
# ============================================

# Now generate all the supplementary plots using stored results
for result in results_for_supplementary:
    num_species = result['num_species']
    resolution = result['resolution']
    
    # Plot spectra with MP bounds
    plot_spectra_with_mp_bounds(
        result['forest_spectrum'],
        result['senm_spectrum'],
        result['lambda_max_forest'],
        result['lambda_max_senm'],
        resolution,
        num_species,
        filename=f'spectra/{num_species}_species_resolution_{resolution}.png',
        show=False,
        verbose=verbose
    )
    
    # Plot eigenvalue density
    plot_eigenvalue_density_vs_mp(
        result['forest_spectrum'],
        result['senm_spectrum'],
        result['bins'],
        resolution,
        num_species,
        result['lambda_max_forest'],
        result['lambda_max_senm'],
        filename=f'spectra/{num_species}_species_density_resolution_{resolution}.png',
        show=False,
        verbose=verbose
    )

# Plot community count vs resolution for each species count
# Group results by num_species
from collections import defaultdict
community_data = defaultdict(lambda: {'resolutions': [], 'forest': [], 'senm': []})

for result in results_for_supplementary:
    ns = result['num_species']
    community_data[ns]['resolutions'].append(result['resolution'])
    community_data[ns]['forest'].append(result['forest_communities'])
    community_data[ns]['senm'].append(result['senm_communities'])

for num_species, data in community_data.items():
    plot_community_count_vs_resolution(
        data['resolutions'],
        data['forest'],
        data['senm'],
        num_species,
        filename=f'spectra/{num_species}_communities_vs_resolution.png',
        show=False,
        verbose=verbose
    )

# ================================================
# === Correlation Matrix Analysis ===
# ================================================

# Compute spectra and get abundance data
_, _, _, bins, senm_abundance, forest_abundance = compute_spectra(filter_resolution)

# Compute correlation matrices
senm_corr = np.corrcoef(senm_abundance)
forest_corr = np.corrcoef(forest_abundance)

# Apply Marchenko-Pastur filter
filtered_senm_corr = MarchenkoPastur(senm_corr, num_species, bins[0]*bins[1], remove_largest=False)
filtered_forest_corr = MarchenkoPastur(forest_corr, num_species, bins[2]*bins[3], remove_largest=False)

# Detect communities on filtered matrices
tau = 1e-3
Th = 1e-4

senm_reordered, senm_CM, senm_idx = detect_communities(filtered_senm_corr, tau, Th)
forest_reordered, forest_CM, forest_idx = detect_communities(filtered_forest_corr, tau, Th)

# Plot SENM correlation matrices
plot_correlation_matrices_comparison(
    senm_corr,
    senm_reordered,
    data_type='SENM',
    filename='correlation_matrices_senm.png',
    show=True
)

# Plot Forest correlation matrices
plot_correlation_matrices_comparison(
    forest_corr,
    forest_reordered,
    data_type='Forest',
    filename='correlation_matrices_forest.png',
    show=True
)
