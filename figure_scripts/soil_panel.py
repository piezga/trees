import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import pandas as pd
import cmocean
import os
from sklearn.linear_model import LinearRegression
from src.config import load_config
from src.compute import compute_spectra, compute_average_correlation

# ============================================================
# Config & Style
# ============================================================

config = load_config()
forest = config['forests']['name']
num_species = config['analysis']['num_species']
path_template = config['forests']['templates']['path_template']
census_template = config['forests']['templates']['census_template']
names_template = config['forests']['templates']['names_template']

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
plt.rcParams['font.size'] = 7
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.5
plt.rcParams['xtick.major.size'] = 2
plt.rcParams['ytick.major.size'] = 2

resolution = 20
n_bins_x = int(1000 / resolution)
n_bins_y = int(500 / resolution)
pairs = [(68, 15), (15, 41), (68, 30)]
species_colors = [
    ['#d62728', '#1f77b4'],
    ['#1f77b4', '#8c564b'],
    ['#d62728', '#2ca02c'],
]

os.makedirs('figures/panel_2', exist_ok=True)

# ============================================================
# Load Species Names
# ============================================================

names_file = (f"{path_template.format(forest=forest)}"
              f"{names_template.format(forest=forest, census=4)}")
with open(names_file, encoding='utf-8-sig') as f:
    species_names = [line.strip() for line in f if line.strip()][:num_species]

# ============================================================
# Load & Bin Soil Data
# ============================================================

print("Loading soil data...")
soil_df = pd.read_excel('soil_data/barro_soil_data.xls')

soil_df['x_bin'] = (soil_df['x'] / resolution).astype(int).clip(0, n_bins_x - 1)
soil_df['y_bin'] = (soil_df['y'] / resolution).astype(int).clip(0, n_bins_y - 1)
soil_binned = soil_df.groupby(['x_bin', 'y_bin'])[['N', 'pH']].mean()


def make_nutrient_grid(soil_binned, nutrient):
    grid = np.full((n_bins_y, n_bins_x), np.nan)
    for (x_bin, y_bin), row in soil_binned.iterrows():
        if 0 <= x_bin < n_bins_x and 0 <= y_bin < n_bins_y:
            grid[y_bin, x_bin] = row[nutrient]
    return grid


grid_N = make_nutrient_grid(soil_binned, 'N')
grid_pH = make_nutrient_grid(soil_binned, 'pH')
print(f"  ✓ Nutrient grids: {grid_N.shape}")

# ============================================================
# Load Forest Data & Compute Correlations
# ============================================================

print("Loading species data...")
(_, _, _, bins, _, forest_abundances) = compute_spectra(
    resolution, num_species=num_species, calculate=False
)
forest_abundance_avg = np.mean(np.array(forest_abundances), axis=0)

raw_corr = np.corrcoef(forest_abundance_avg)

nutrient_columns = ['Al', 'B', 'Ca', 'Cu', 'Fe', 'K', 'Mg', 'Mn', 'P', 'Zn']
soil_binned_full = soil_df.groupby(['x_bin', 'y_bin'])[nutrient_columns].mean()

n_sites = n_bins_x * n_bins_y
nutrient_data = np.zeros((len(nutrient_columns), n_sites))
for x_bin in range(n_bins_x):
    for y_bin in range(n_bins_y):
        idx = x_bin * n_bins_y + y_bin
        if (x_bin, y_bin) in soil_binned_full.index:
            nutrient_data[:, idx] = soil_binned_full.loc[(x_bin, y_bin)].values
        else:
            nutrient_data[:, idx] = np.nan

for i in range(nutrient_data.shape[0]):
    col_mean = np.nanmean(nutrient_data[i, :])
    nutrient_data[i, np.isnan(nutrient_data[i, :])] = col_mean

X_species = forest_abundance_avg.T
X_nutrients = nutrient_data.T
residuals = np.zeros_like(X_species)
for i in range(num_species):
    reg = LinearRegression()
    reg.fit(X_nutrients, X_species[:, i])
    residuals[:, i] = X_species[:, i] - reg.predict(X_nutrients)
stripped_corr = np.corrcoef(residuals.T)

delta = stripped_corr - raw_corr
np.fill_diagonal(delta, np.nan)
np.fill_diagonal(raw_corr, np.nan)
np.fill_diagonal(stripped_corr, np.nan)

print("  ✓ Correlations computed")

# ============================================================
# Load Forest Map & Census Data
# ============================================================

print("Loading BCI map and census data...")
bci_img = plt.imread('level.png')

census_file = (f"{path_template.format(forest=forest)}"
               f"{census_template.format(forest=forest, census=4)}")
df_census = pd.read_csv(census_file)

# ============================================================
# Helper: shared axes style
# ============================================================

def style_ax(ax):
    ax.tick_params(labelsize=6)
    ax.grid(True, alpha=0.3, linewidth=0.3)


# ============================================================
# Panel a: N heatmap
# ============================================================

fig, ax = plt.subplots(figsize=(3.5, 2.0), dpi=300)
im = ax.imshow(grid_N, origin='lower', cmap='YlGn',
               extent=[0, 1000, 0, 500], aspect='equal',
               interpolation='nearest')
cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
cbar.ax.tick_params(labelsize=6)
cbar.set_label('N (g/kg)', fontsize=7)
ax.set_xlabel('X (m)', fontsize=7)
ax.set_ylabel('Y (m)', fontsize=7)
style_ax(ax)
fig.tight_layout()
fig.savefig('figures/panel_2/panel_a_N.svg', format='svg', bbox_inches='tight')
plt.close(fig)
print("  ✓ panel_a_N.svg")

# ============================================================
# Panel b: pH heatmap
# ============================================================

fig, ax = plt.subplots(figsize=(3.5, 2.0), dpi=300)
im = ax.imshow(grid_pH, origin='lower', cmap='RdYlBu',
               extent=[0, 1000, 0, 500], aspect='equal',
               interpolation='nearest')
cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
cbar.ax.tick_params(labelsize=6)
cbar.set_label('pH', fontsize=7)
ax.set_xlabel('X (m)', fontsize=7)
ax.set_ylabel('Y (m)', fontsize=7)
style_ax(ax)
fig.tight_layout()
fig.savefig('figures/panel_2/panel_b_pH.svg', format='svg', bbox_inches='tight')
plt.close(fig)
print("  ✓ panel_b_pH.svg")

# ============================================================
# Panel c: Precomputed scatter image
# ============================================================

fig, ax = plt.subplots(figsize=(3.5, 3.5), dpi=300)
img_c = plt.imread('stripped_correlation_analysis/scatter_raw_vs_stripped_20m.png')
ax.imshow(img_c, aspect='auto')
ax.axis('off')
fig.tight_layout()
fig.savefig('figures/panel_2/panel_c_scatter.svg', format='svg', bbox_inches='tight')
plt.close(fig)
print("  ✓ panel_c_scatter.svg")

# ============================================================
# Panel d: Delta heatmap
# ============================================================

fig, ax = plt.subplots(figsize=(3.5, 3.5), dpi=300)
vmax = np.nanpercentile(np.abs(delta), 95)
norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
im_d = ax.imshow(delta, cmap=cmocean.cm.balance, norm=norm, rasterized=True)
cbar = plt.colorbar(im_d, ax=ax, fraction=0.046, pad=0.04)
cbar.ax.tick_params(labelsize=6)
cbar.set_label(r'$\Delta\rho$', fontsize=7)
ax.set_xlabel('Species', fontsize=7)
ax.set_ylabel('Species', fontsize=7)
style_ax(ax)
fig.tight_layout()
fig.savefig('figures/panel_2/panel_d_delta.svg', format='svg', bbox_inches='tight')
plt.close(fig)
print("  ✓ panel_d_delta.svg")

# ============================================================
# Panels e, f, g: Species raster plots
# ============================================================

panel_labels = ['e', 'f', 'g']
panel_names  = ['panel_e_species', 'panel_f_species', 'panel_g_species']

for i, (pair, colors, label, name) in enumerate(
        zip(pairs, species_colors, panel_labels, panel_names)):

    fig, ax = plt.subplots(figsize=(3.5, 2.0), dpi=300)
    ax.imshow(
        bci_img[::-1], origin='lower', alpha=0.9,
        extent=[0, 1000, 0, 500], aspect='equal'
    )
    for idx, sp_idx in enumerate(pair):
        sp_name = species_names[sp_idx]
        sp_data = df_census[df_census['name'] == sp_name.strip()]
        ax.scatter(
            sp_data['x'], sp_data['y'],
            s=0.2, color=colors[idx], alpha=0.9,
            rasterized=True, label=f'{sp_name[:15]}'
        )
    ax.set_xlim(0, 1000)
    ax.set_ylim(0, 500)
    ax.set_xlabel('X (m)', fontsize=7)
    ax.set_ylabel('Y (m)', fontsize=7)
    style_ax(ax)
    ax.legend(fontsize=5, frameon=False, loc='upper right',
              markerscale=4, handletextpad=0.3)
    fig.tight_layout()
    fig.savefig(f'figures/panel_2/{name}.svg', format='svg', bbox_inches='tight')
    plt.close(fig)
    print(f"  ✓ {name}.svg")

print("\n✓ All panels saved to figures/panel_2/")
