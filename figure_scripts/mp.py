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

res_28_senm = np.load(f'quantities/{forest}_senm_spectrum_{num_species}_28.npy') 
res_9_senm = np.load(f'quantities/{forest}_senm_spectrum_{num_species}_9.npy')

lambda_min_forest_28, lambda_max_forest_28 = marchenko_pastur_bounds(num_species, res_28_bins[2], res_28_bins[3])
lambda_min_forest_9, lambda_max_forest_9 = marchenko_pastur_bounds(num_species, res_9_bins[2], res_9_bins[3])

lambda_min_senm_28, lambda_max_senm_28 = marchenko_pastur_bounds(num_species, res_28_bins[0], res_28_bins[1])
lambda_min_senm_9, lambda_max_senm_9 = marchenko_pastur_bounds(num_species, res_9_bins[0], res_9_bins[1])


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

# === Panel A: Spectra with MP ===

# === Figure ===
fig, ax = plt.subplots(figsize=(3, 3), dpi=300)


x = np.arange(1, num_species + 1)
forest_color = '#2ca02c'
senm_color = '#1f77b4'
forest_mp_color = '#98df8a'  
senm_mp_color = '#aec7e8'   

# Plot mean forest spectrum (main)
mean_forest_9 = np.mean(res_9_forest_list, axis=0)
ax.plot(x, mean_forest_9, 'o-', color=forest_color, markersize=2,
        markerfacecolor='white', markeredgewidth=0.5, linewidth=0.5,
        label='Forest')

# Plot SENM spectrum (main)
ax.plot(x, res_9_senm, 's-', color=senm_color, markersize=2,
        markerfacecolor='white', markeredgewidth=0.5, linewidth=0.5,
        label='SENM')

# MP bounds for forest (main)
ax.axhline(lambda_max_forest_9, color=forest_mp_color, linestyle='--', linewidth=0.8, alpha=0.8)
ax.axhline(lambda_min_forest_9, color=forest_mp_color, linestyle='--', linewidth=0.8, alpha=0.8)
#ax.fill_between(x, lambda_min_forest_9, lambda_max_forest_9, 
#                color=forest_mp_color, alpha=0.5, linewidth=0, label='MP forest')

# MP bounds for SENM (main)
ax.axhline(lambda_max_senm_9, color=senm_mp_color, linestyle='--', linewidth=0.8, alpha=0.8)
ax.axhline(lambda_min_senm_9, color=senm_mp_color, linestyle='--', linewidth=0.8, alpha=0.8)
#ax.fill_between(x, lambda_min_senm_9, lambda_max_senm_9, 
#                color=senm_mp_color, alpha=0.5, linewidth=0, label='MP SENM')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1, 100)
ax.set_ylim(5e-1, 3.5)
ax.set_xlabel(r'$\lambda$ rank', fontsize=8)
ax.set_ylabel(r'$\lambda$', fontsize=8)
ax.grid(False)
ax.tick_params(labelsize=7)
ax.legend(fontsize=6, loc='lower left')

# Inset
axins = inset_axes(ax, width="35%", height="35%", loc='upper right', borderpad=1)

# Plot mean forest spectrum (inset)
mean_forest_28 = np.mean(res_28_forest_list, axis=0)
axins.plot(x, mean_forest_28, 'o-', color=forest_color, markersize=1.5,
           markerfacecolor='white', markeredgewidth=0.4, linewidth=0.4)

# Plot SENM spectrum (inset)
axins.plot(x, res_28_senm, 's-', color=senm_color, markersize=1.5,
           markerfacecolor='white', markeredgewidth=0.4, linewidth=0.4)

# MP bounds for forest (inset)
axins.axhline(lambda_max_forest_28, color=forest_mp_color, linestyle='--', linewidth=0.8, alpha=0.8)
axins.axhline(lambda_min_forest_28, color=forest_mp_color, linestyle='--', linewidth=0.8, alpha=0.8)
#axins.fill_between(x, lambda_min_forest_28, lambda_max_forest_28, 
#                   color=forest_mp_color, alpha=0.5, linewidth=0)

# MP bounds for SENM (inset)
axins.axhline(lambda_max_senm_28, color=senm_mp_color, linestyle='--', linewidth=0.8, alpha=0.8)
axins.axhline(lambda_min_senm_28, color=senm_mp_color, linestyle='--', linewidth=0.8, alpha=0.8)
#axins.fill_between(x, lambda_min_senm_28, lambda_max_senm_28, 
#                  color=senm_mp_color, alpha=0.5, linewidth=0)

axins.set_xscale('log')
axins.set_yscale('log')
axins.set_xlim(1, 100)
axins.set_ylim(0.15, 10)
axins.set_title('28 m', fontsize=7, pad=2)
axins.grid(False)
axins.tick_params(labelsize=5)

plt.savefig('figures/mp/spectra.svg', dpi=300, bbox_inches='tight')


# === Panel B: Eigenvalue Density vs Marchenko-Pastur ===

def marchpast(l, g):
    "Marchenko-Pastur distribution"
    def m0(a):
        return np.maximum(a, np.zeros_like(a))
    gplus  = (1 + g**0.5)**2
    gminus = (1 - g**0.5)**2
    return np.sqrt(m0(gplus - l) * m0(l - gminus)) / (2 * np.pi * g * l)

# g = p/n: num_species / num_sites

chosen_res = 28

if chosen_res == 28:
    n_bins_x = res_28_bins[0]
    n_bins_y = res_28_bins[1]
    all_eigs = np.concatenate(res_28_forest_list)
elif chosen_res == 9:
    n_bins_x = res_9_bins[0]
    n_bins_y = res_9_bins[1]
    all_eigs = np.concatenate(res_9_forest_list)

g = num_species / (n_bins_x * n_bins_y)


bins = np.logspace(np.log10(all_eigs.min()), np.log10(all_eigs.max()), 40)
counts, edges = np.histogram(all_eigs, bins=bins)
centers = np.sqrt(edges[:-1] * edges[1:])          # geometric bin centers
density  = counts / (all_eigs.size * (edges[1:] - edges[:-1]))

l_mp = np.linspace((1 - g**0.5)**2 * 1.001, (1 + g**0.5)**2 * 0.999, 500)
rho_mp = marchpast(l_mp, g)

fig, ax = plt.subplots(figsize=(3, 3), dpi=300)

ax.plot(centers, density, 
       color='#ff7f0e', alpha=0.55, linestyle='--', label='Empirical')
ax.scatter(centers, density, 
       color='#ff7f0e', alpha=1, )
ax.plot(l_mp, rho_mp, color='#1f77b4', linewidth=0.8,
        linestyle='--', label='Marchenko–Pastur')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$\lambda$', fontsize=8)
ax.set_ylabel(r'$\rho(\lambda)$', fontsize=8)
ax.tick_params(labelsize=7, which='both', direction='in',
               top=True, right=True, width=0.5, length=2)
ax.legend(fontsize=6, frameon=False)
for sp in ax.spines.values():
    sp.set_linewidth(0.5)

fig.tight_layout()
fig.savefig('figures/mp/eigenvalue_density.svg', format='svg',
            dpi=300, bbox_inches='tight')

# === Panel C: Communities vs Resolution ===
fig, ax = plt.subplots(figsize=(3, 3), dpi=300)

ax.plot(resolutions, forest_communities, 'o-', color='#2ca02c', 
        markersize=3, linewidth=1, label='Forest', markerfacecolor='white', markeredgewidth=0.5, alpha=0.7)
ax.plot(resolutions, senm_communities, 's--', color='#1f77b4', 
        markersize=3, linewidth=1, label='SENM', markerfacecolor='white', markeredgewidth=0.5, alpha = 0.7)

ax.set_xlabel('Scale (m)', fontsize=8)
ax.set_ylabel(r'$N_{\mathrm{comm}}$', fontsize=8)
ax.legend(fontsize=7, frameon=False, loc='best')
ax.grid(False)
ax.tick_params(labelsize=7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# === Save ===
plt.savefig('figures/mp/communities.svg', dpi=300, bbox_inches='tight')
plt.show()
