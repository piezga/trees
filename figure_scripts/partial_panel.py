import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cmocean
import os
from src.config import load_config

# === Config & Style ===
config = load_config()
forest = config['forests']['name']
num_species = config['analysis']['num_species']
resolution = 20

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 0.5

# === Load species names ===
path_template = config['forests']['templates']['path_template']
names_template = config['forests']['templates']['names_template']
names_file = f"{path_template.format(forest=forest)}{names_template.format(forest=forest, census=4)}"
with open(names_file, encoding='utf-8-sig') as f:
    species_names = [line.strip() for line in f if line.strip()][:num_species]

# Find indices
sp_dict = {name: idx for idx, name in enumerate(species_names)}
pairs_names = [('hirttr', 'guatdu'), ('hirttr', 'licapl'), ('guatdu', 'licapl')]
pairs_full_names = [('H. triandra', 'G. lucens'), 
                    ('H. triandra', 'L. platypus'), 
                    ('G. lucens', 'L. platypus')]
pairs_idx = [(sp_dict[a], sp_dict[b]) for a, b in pairs_names]
# === Load matrices ===
input_dir = "precision_matrix_analysis_shrinkage"
C_stripped = np.load(f'stripped_correlation_analysis/stripped_correlation_20m.npy')
P_partial = np.load(f'{input_dir}/P_partial_correlation.npy')

np.fill_diagonal(C_stripped, np.nan)
np.fill_diagonal(P_partial, np.nan)

# === Figure ===
fig, axes = plt.subplots(1, 2, figsize=(5, 2.5), dpi=300, constrained_layout=True)

cmap = cmocean.cm.balance
vmin, vmax = -0.5, 0.5
norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

# === Panel a: Stripped + filtered ===
im1 = axes[0].imshow(C_stripped, cmap=cmap, norm=norm, rasterized=True)
axes[0].set_xlabel('Species', fontsize=8)
axes[0].set_ylabel('Species', fontsize=8)
axes[0].set_title('a  Soil regression + MP filter', loc='left', fontweight='bold', fontsize=9)
axes[0].tick_params(labelsize=7)
cbar1 = plt.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04)
cbar1.ax.tick_params(labelsize=6)

# === Panel b: Partial correlation (denoised) ===
im2 = axes[1].imshow(P_partial, cmap=cmap, norm=norm, rasterized=True)
axes[1].set_xlabel('Species', fontsize=8)
axes[1].set_ylabel('Species', fontsize=8)
axes[1].set_title('b  Partial correlation', loc='left', fontweight='bold', fontsize=9)
axes[1].tick_params(labelsize=7)
cbar2 = plt.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04)
cbar2.ax.tick_params(labelsize=6)

plt.savefig('figures/partial_correlation_result.svg', dpi=600, bbox_inches='tight')
plt.show()
# === Panel c: Scatter ===
fig, ax = plt.subplots(1, 1, figsize=(3, 3), dpi=300, constrained_layout=True)

# Get upper triangle
triu = np.triu_indices(num_species, k=1)
c_vals = C_stripped[triu]
p_vals = P_partial[triu]
valid = ~(np.isnan(c_vals) | np.isnan(p_vals))

ax.scatter(c_vals[valid], p_vals[valid], s=0.5, alpha=0.2, 
           color='gray', rasterized=True)

# Highlight pairs
colors = ['#d62728', '#2ca02c', '#1f77b4']
for idx, (i, j) in enumerate(pairs_idx):
    c_val = C_stripped[i, j]
    p_val = P_partial[i, j]
    ax.scatter(c_val, p_val, s=30, color=colors[idx], 
               edgecolor='black', linewidth=0.5, zorder=10,
               label=f'{pairs_full_names[idx][0]}-{pairs_full_names[idx][1]}')

# 1:1 line
lim = max(abs(c_vals[valid].min()), abs(c_vals[valid].max()),
          abs(p_vals[valid].min()), abs(p_vals[valid].max()))
ax.plot([-lim, lim], [-lim, lim], 'k-', linewidth=0.5, alpha=0.5)
ax.axhline(0, color='gray', linestyle='-', linewidth=0.3, alpha=0.5)
ax.axvline(0, color='gray', linestyle='-', linewidth=0.3, alpha=0.5)

ax.set_xlabel('Corr. (soil regression)', fontsize=8)
ax.set_ylabel('Partial corr. (denoised)', fontsize=8)
ax.set_title('c', loc='left', fontweight='bold', fontsize=9)
ax.tick_params(labelsize=7)
ax.legend(fontsize=5, frameon=False, loc='best', markerscale=0.8)
ax.set_aspect('equal')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# === Save ===
plt.savefig('figures/partial_correlation_scatter.svg', dpi=600, bbox_inches='tight')
plt.show()
