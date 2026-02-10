import numpy as np
import argparse
import matplotlib
import matplotlib.pyplot as plt

from src.config import load_config
from src.compute import (
    compute_spectra,
    marchenko_pastur_bounds,
    square_diff_above_MP,
    number_of_communities,
)
from src.plotting import (
    plot_dual_spectra_with_inset,
    plot_size_effect_panel,
    supplementary_plots,
)

matplotlib.use('TkAgg')

# === Config ===
config = load_config()
forest = config['forests']['name']  
num_species = config['analysis']['num_species']
species_array = config['analysis']['species_array']

# === Args ===
def parse_args():
    parser = argparse.ArgumentParser(description='Run spectral analyses')
    parser.add_argument('-c', '--calculate', action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-s', '--supplementary', action='store_true')
    return parser.parse_args()

args = parse_args()
calculate = args.calculate 
verbose = args.verbose
supplementary = args.supplementary

# === Nature style ===
matplotlib.rcParams['font.size'] = 8
matplotlib.rcParams['axes.labelsize'] = 8
matplotlib.rcParams['axes.titlesize'] = 10
matplotlib.rcParams['xtick.labelsize'] = 7
matplotlib.rcParams['ytick.labelsize'] = 7
matplotlib.rcParams['legend.fontsize'] = 7
matplotlib.rcParams['axes.linewidth'] = 0.5
matplotlib.rcParams['xtick.major.width'] = 0.5
matplotlib.rcParams['ytick.major.width'] = 0.5
matplotlib.rcParams['xtick.major.size'] = 2
matplotlib.rcParams['ytick.major.size'] = 2
matplotlib.rcParams['lines.linewidth'] = 1
matplotlib.rcParams['lines.markersize'] = 3

# === Figure: 183mm wide (double column) ===
fig, ax_dict = plt.subplot_mosaic("AB", figsize=(7.08, 3.54), dpi=300, constrained_layout=True)

# === Panel A ===
ax_a = ax_dict["A"]

if calculate:
    res_50_senm, res_50_senm_std, res_50_forest_list, _, _, _ = compute_spectra(50, calculate)
    res_5_senm, res_5_senm_std, res_5_forest_list, _, _, _ = compute_spectra(5, calculate)
else: 
    res_50_senm = np.load(f'quantities/{forest}_senm_spectrum_{num_species}_50.npy')
    res_50_senm_std = np.load(f'quantities/{forest}_senm_std_{num_species}_50.npy')
    res_50_forest_list = np.load(f'quantities/{forest}_forest_spectra_{num_species}_50.npy')
    res_5_senm = np.load(f'quantities/{forest}_senm_spectrum_{num_species}_5.npy')
    res_5_senm_std = np.load(f'quantities/{forest}_senm_std_{num_species}_5.npy')
    res_5_forest_list = np.load(f'quantities/{forest}_forest_spectra_{num_species}_5.npy')

plot_dual_spectra_with_inset(
    res_50_senm, res_50_senm_std, res_50_forest_list,
    res_5_senm, res_5_senm_std, res_5_forest_list,
    num_species=num_species, ax=ax_a, show=False
)
ax_a.set_title('a', loc='left', fontweight='bold', pad=4)

# === Panel B ===
ax_b = ax_dict["B"]
resolutions = list(range(4, 34, 1))
results_for_panel_b = {}
results_for_supplementary = []

for idx, num_species in enumerate(species_array):
    comm_diffs = []
    
    for resolution in resolutions:
        if calculate:
            senm_spectrum, _, forest_spectra, bins, _, _ = compute_spectra(resolution, calculate)
        else:
            senm_spectrum = np.load(f'quantities/{forest}_senm_spectrum_{num_species}_{resolution}.npy')
            forest_spectra = np.load(f'quantities/{forest}_forest_spectra_{num_species}_{resolution}.npy')
            bins = np.load(f'quantities/{forest}_bins_{num_species}_{resolution}.npy')
        
        mean_forest_spectrum = np.mean(np.array(forest_spectra), axis=0)
        _, lambda_max_forest = marchenko_pastur_bounds(num_species, bins[2], bins[3])
        _, lambda_max_senm = marchenko_pastur_bounds(num_species, bins[0], bins[1])
        
        square_diff, comm_diff = square_diff_above_MP(senm_spectrum, mean_forest_spectrum, 
                                                      lambda_max_senm, lambda_max_forest)
        forcom = number_of_communities(mean_forest_spectrum, lambda_max_forest)
        senmcom = number_of_communities(senm_spectrum, lambda_max_senm)
        
        comm_diffs.append(comm_diff)
        results_for_supplementary.append({
            'num_species': num_species, 'resolution': resolution,
            'forest_spectrum': mean_forest_spectrum, 'senm_spectrum': senm_spectrum,
            'lambda_max_forest': lambda_max_forest, 'lambda_max_senm': lambda_max_senm,
            'bins': bins, 'forest_communities': forcom, 'senm_communities': senmcom
        })
    
    results_for_panel_b[num_species] = comm_diffs

plot_size_effect_panel(resolutions, results_for_panel_b, species_array, ax=ax_b, show=False)
ax_b.set_title('b', loc='left', fontweight='bold', pad=4)

# === Save ===
plt.savefig('figures/Fig_spectra.pdf', dpi=300, bbox_inches='tight')
plt.savefig('figures/Fig_spectra.png', dpi=600, bbox_inches='tight')
plt.show()

if supplementary:
    supplementary_plots(results_for_supplementary, verbose)
