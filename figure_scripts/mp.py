import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from src.config import load_config
from src.compute import marchenko_pastur_bounds

# === Config ===
config = load_config()
forest = config['forests']['name']
num_species = config['analysis']['num_species']

# === Nature style ===
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.5
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['lines.markersize'] = 3

# === Load data for Panel A ===
res_28_forest_list = np.load(f'quantities/{forest}_forest_spectra_{num_species}_28.npy')
res_28_bins = np.load(f'quantities/{forest}_bins_{num_species}_28.npy')
res_9_forest_list = np.load(f'quantities/{forest}_forest_spectra_{num_species}_9.npy')
res_9_bins = np.load(f'quantities/{forest}_bins_{num_species}_9.npy')

lambda_min_forest_28, lambda_max_forest_28 = marchenko_pastur_bounds(num_species, res_28_bins[2], res_28_bins[3])
lambda_min_forest_9, lambda_max_forest_9 = marchenko_pastur_bounds(num_species, res_9_bins[2], res_9_bins[3])

# === Load data for Panel B ===
resolutions = np.arange(4, 34, 1)
forest_communities = []
senm_communities = []

for res in resolutions:
    forest_spec = np.load(f'quantities/{forest}_forest_spectra_{num_species}_{res}.npy')
    senm_spec = np.load(f'quantities/{forest}_senm_spectrum_{num_species}_{res}.npy')
    bins = np.load(f'quantities/{forest}_bins_{num_species}_{res}.npy')
    
    _, lambda_max_f = marchenko_pastur_bounds(num_species, bins[2], bins[3])
    _, lambda_max_s = marchenko_pastur_bounds(num_species, bins[0], bins[1])
    
    mean_forest = np.mean(forest_spec, axis=0)
    n_comm_f = np.sum(mean_forest > lambda_max_f)
    n_comm_s = np.sum(senm_spec > lambda_max_s)
    
    forest_communities.append(n_comm_f)
    senm_communities.append(n_comm_s)

# === Figure ===
fig, axes = plt.subplots(1, 2, figsize=(7.08, 3), dpi=300, constrained_layout=True)

# === Panel A: Spectra with MP ===
ax = axes[0]
x = np.arange(1, num_species + 1)
color = '#ff7f0e'
mp_color = '#fdbf6f'

for spectrum in res_9_forest_list:
    ax.plot(x, spectrum, 'o-', color=color, markersize=2, alpha=0.5,
            markerfacecolor='white', markeredgewidth=0.5, linewidth=0.5)

ax.axhline(lambda_max_forest_9, color=mp_color, linestyle='--', linewidth=0.2, alpha=0.8)
ax.axhline(lambda_min_forest_9, color=mp_color, linestyle='--', linewidth=0.2, alpha=0.8)
ax.fill_between(x, lambda_min_forest_9, lambda_max_forest_9, 
                color=mp_color, alpha=0.5, linewidth=0)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1, 100)
ax.set_ylim(5e-1, 3.5)
ax.set_xlabel(r'$\lambda$ rank', fontsize=8)
ax.set_ylabel(r'$\lambda$', fontsize=8)
ax.set_title('(a)', loc='left', fontweight='bold', fontsize=10, pad=4)
ax.grid(True, which="both", linestyle='--', alpha=0.3, linewidth=0.3)
ax.tick_params(labelsize=7)

# Inset
axins = inset_axes(ax, width="35%", height="35%", loc='upper right', borderpad=1)
for spectrum in res_28_forest_list:
    axins.plot(x, spectrum, 'o-', color=color, markersize=1.5, alpha=0.5,
               markerfacecolor='white', markeredgewidth=0.4, linewidth=0.4)

axins.axhline(lambda_max_forest_28, color=mp_color, linestyle='--', linewidth=0.2, alpha=0.8)
axins.axhline(lambda_min_forest_28, color=mp_color, linestyle='--', linewidth=0.2, alpha=0.8)
axins.fill_between(x, lambda_min_forest_28, lambda_max_forest_28, 
                    color=mp_color, alpha=0.5, linewidth=0)

axins.set_xscale('log')
axins.set_yscale('log')
axins.set_xlim(1, 100)
axins.set_ylim(0.15, 10)
axins.set_title('28 m', fontsize=7, pad=2)
axins.grid(True, which="both", linestyle='--', alpha=0.3, linewidth=0.25)
axins.tick_params(labelsize=5)

# === Panel B: Communities vs Resolution ===
ax = axes[1]

ax.plot(resolutions, forest_communities, 'o-', color='#2ca02c', 
        markersize=3, linewidth=1, label='Forest', markerfacecolor='white', markeredgewidth=0.5, alpha=0.7)
ax.plot(resolutions, senm_communities, 's--', color='#1f77b4', 
        markersize=3, linewidth=1, label='SENM', markerfacecolor='white', markeredgewidth=0.5, alpha = 0.7)

ax.set_xlabel('Scale (m)', fontsize=8)
ax.set_ylabel(r'$N_{\mathrm{comm}}$', fontsize=8)
ax.set_title('(b)', loc='left', fontweight='bold', fontsize=10, pad=4)
ax.legend(fontsize=7, frameon=False, loc='best')
ax.grid(True, alpha=0.3, linewidth=0.3)
ax.tick_params(labelsize=7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# === Save ===
plt.savefig('figures/Fig_MP_method.svg', dpi=300, bbox_inches='tight')
plt.show()
