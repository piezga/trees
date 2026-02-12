import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
from src.config import load_config
import matplotlib.image as mpimg
# === Nature style ===
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 0.5

# === Parameters ===
config = load_config()
forest = 'barro'
census = 4
species_indices = [13, 56, 22]  # Example species
L = 100  # Grid resolution in meters
colors = ['#d62728', '#2ca02c', '#1f77b4']  # Red, green, blue

# === Load data ===
path_template = config['forests']['templates']['path_template']
census_template = config['forests']['templates']['census_template']
names_template = config['forests']['templates']['names_template']

forest_file = f"{path_template.format(forest=forest)}{census_template.format(forest=forest, census=census)}"
df = pd.read_csv(forest_file)

names_file = f"{path_template.format(forest=forest)}{names_template.format(forest=forest, census=4)}"
with open(names_file, encoding='utf-8-sig') as f:
    species_names = [line.strip() for line in f if line.strip()][:100]

# === Grid dimensions ===
n_bins_x = int(1000 / L)
n_bins_y = int(500 / L)
n_total = n_bins_x * n_bins_y

# === Figure ===
fig = plt.figure(figsize=(7.08, 6), dpi=300)
gs = fig.add_gridspec(4, 1, height_ratios=[3, 1, 1, 1], hspace=0.4)

# === Panel A: Spatial map with grid ===
ax_map = fig.add_subplot(gs[0])
plot_names = ['T. panamense', 'T. nervosa', 'C. curvigemmia']
img = mpimg.imread('level.png')
ax_map.imshow(img, extent=[0, 1000, 0, 500], aspect='auto', zorder=0)
for idx, sp_idx in enumerate(species_indices):
    sp_name = species_names[sp_idx]
    sp_data = df[df['name'] == sp_name.strip()]
    ax_map.scatter(sp_data['x'], sp_data['y'], s=1, color=colors[idx], 
                   alpha=0.6, rasterized=True, label=plot_names[idx])

# Grid lines
for i in range(n_bins_x + 1):
    ax_map.axvline(i * L, color='black', linewidth=0.5, alpha=0.8)
for j in range(n_bins_y + 1):
    ax_map.axhline(j * L, color='black', linewidth=0.5, alpha=0.8)

ax_map.set_xlim(0, 1000)
ax_map.set_ylim(0, 500)
ax_map.set_xlabel('X (m)', fontsize=8)
ax_map.set_ylabel('Y (m)', fontsize=8)
ax_map.set_title('a', loc='left', fontweight='bold', fontsize=10)
ax_map.tick_params(labelsize=7)
ax_map.set_aspect('equal')
ax_map.legend(fontsize=6, loc='upper left', bbox_to_anchor = (1,1), frameon=False, markerscale=3)

# === Panels B, C, D: Flattened grids ===
for idx, sp_idx in enumerate(species_indices):
    ax = fig.add_subplot(gs[idx + 1])
    
    sp_name = species_names[sp_idx]
    sp_data = df[df['name'] == sp_name.strip()]
    
    # Bin counts
    abundance = np.zeros(n_total)
    for _, row in sp_data.iterrows():
        x_bin = int(row['x'] / L)
        y_bin = int(row['y'] / L)
        if 0 <= x_bin < n_bins_x and 0 <= y_bin < n_bins_y:
            bin_idx = x_bin * n_bins_y + y_bin
            abundance[bin_idx] += 1
    
    # Plot as bars
    ax.bar(range(n_total), abundance, width=1, color=colors[idx], 
           edgecolor='none', alpha=0.6, rasterized=True)
    
    ax.set_xlim(0, n_total)
    ax.set_ylim(0, abundance.max() * 1.1)
    ax.set_ylabel('Count', fontsize=7)
    ax.tick_params(labelsize=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    panel_label = chr(98 + idx)  # b, c, d
    ax.set_title(panel_label, loc='left', fontweight='bold', fontsize=10)
    
    if idx == len(species_indices) - 1:
        ax.set_xlabel(f'Spatial bin (L = {L} m)', fontsize=7)
    else:
        ax.set_xticklabels([])

plt.savefig('figures/Fig_binning_procedure.svg', dpi=300, bbox_inches='tight')
plt.savefig('figures/Fig_binning_procedure.png', dpi=600, bbox_inches='tight')
plt.show()
