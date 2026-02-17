import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import cmocean
import os
from sklearn.linear_model import LinearRegression

from src.config import load_config
from src.compute import (compute_spectra, compute_average_correlation, 
                         MarchenkoPastur, compute_stripped_correlation_matrix)

# ========================================
# Load Configuration
# ========================================

config = load_config()
num_species = config['analysis']['num_species']
path_template = config['forests']['templates']['path_template']
names_template = config['forests']['templates']['names_template']

# Parameters
resolution = 20
forest = 'barro'
debug = True

print(f"\n{'='*80}")
print(f"RAW vs stripped CORRELATION ANALYSIS")
print(f"Resolution: {resolution}m")
print(f"{'='*80}\n")

# ========================================
# Load Species Data
# ========================================

print("Loading species abundance data...")
(senm_mean, senm_std, forest_spectra, 
 bins, senm_abundance, forest_abundances) = compute_spectra(
    resolution, num_species=num_species, calculate=False
)

# Compute average forest abundance across censuses
forest_abundance_avg = np.mean(np.array(forest_abundances), axis=0)
print(f"  Species abundance shape: {forest_abundance_avg.shape}")
print(f"  (species={forest_abundance_avg.shape[0]}, sites={forest_abundance_avg.shape[1]})")

# Load species names
names_file = f"{path_template.format(forest=forest)}{names_template.format(forest=forest, census=4)}"
with open(names_file, 'r', encoding='utf-8-sig') as f:
    species_names = [line.strip() for line in f.readlines()][:num_species]

# ========================================
# Load Nutrient Data
# ========================================

print("\nLoading nutrient data...")

nutrient_file = 'soil_data/barro_soil_data.xls'

# Replace this with your actual data loading code
if os.path.exists(nutrient_file):
    import pandas as pd
    nutrient_df = pd.read_excel(nutrient_file)
    
    # Assuming columns are: x, y, Al,  etc
    nutrient_columns = ['Al', 'B', 'Ca', 'Cu', 'Fe', 'K', 'Mg', 'Mn', 'P', 'Zn', 'N', 'N(min)', 'pH']
    
    # Create spatial bins matching species data
    n_bins_x = bins[2]
    n_bins_y = bins[3]
    
    # Bin the nutrient data to match species resolution
    nutrient_df['x_bin'] = (nutrient_df['x'] / (1000 / n_bins_x)).astype(int).clip(0, n_bins_x - 1)
    nutrient_df['y_bin'] = (nutrient_df['y'] / (500 / n_bins_y)).astype(int).clip(0, n_bins_y - 1)
    
    # Average nutrients within each bin
    nutrient_binned = nutrient_df.groupby(['x_bin', 'y_bin'])[nutrient_columns].mean()
    
    # Convert to array format (n_nutrients, n_sites)
    n_sites = n_bins_x * n_bins_y
    nutrient_data = np.zeros((len(nutrient_columns), n_sites))
    
    for idx, (x_bin, y_bin) in enumerate([(x, y) for x in range(n_bins_x) for y in range(n_bins_y)]):
        if (x_bin, y_bin) in nutrient_binned.index:
            nutrient_data[:, idx] = nutrient_binned.loc[(x_bin, y_bin)].values
        else:
            # Fill missing bins with mean or nearest neighbor
            nutrient_data[:, idx] = np.nan
    
    # Fill NaN values with column means
    for i in range(nutrient_data.shape[0]):
        col_mean = np.nanmean(nutrient_data[i, :])
        nutrient_data[i, np.isnan(nutrient_data[i, :])] = col_mean
    
    print(f"  Nutrient data shape: {nutrient_data.shape}")
    print(f"  Nutrients: {nutrient_columns}")

else:
    print(f"  Warning: {nutrient_file} not found")

# ========================================
# Compute stripped Correlation Matrix
# ========================================


print("\nComputing raw correlation matrix...")
unfiltered_raw_corr = np.corrcoef(forest_abundance_avg)

unfiltered_stripped_corr, stripped_abund = compute_stripped_correlation_matrix(forest_abundance_avg, nutrient_data)


# Filtering scatter_raw_vs_stripped_
T = n_sites
raw_corr = MarchenkoPastur(unfiltered_raw_corr, 100, T)
stripped_corr = MarchenkoPastur(unfiltered_stripped_corr, 100, T)

print(f"  ✓ Raw correlation matrix shape: {raw_corr.shape}")
print(f"  ✓ stripped correlation matrix shape: {stripped_corr.shape}")

# Create the output directory if it doesn't exist
os.makedirs("stripped_correlation_analysis", exist_ok=True)

print(f'Shape of abund is {stripped_abund.shape}')
# ========================================
# Compute Differences 
# ========================================

print("\nComputing differences ...")

# difference
delta = stripped_corr - raw_corr
# Set diagonal to NaN (not meaningful)
np.fill_diagonal(delta, np.nan)
np.fill_diagonal(raw_corr, np.nan)
np.fill_diagonal(stripped_corr, np.nan)
# ========================================
# Identify Interaction Types
# ========================================

print("\nClassifying species pair interactions...")

threshold_high = 0.3
threshold_low = 0.1

direct = []          # High raw, high stripped → direct interaction
indirect = []        # High raw, low stripped → nutrient-mediated
weak = []            # Low raw, low stripped → weak/no interaction
emerging = []        # Low raw, high stripped → suppressed by nutrients

for i in range(num_species):
    for j in range(i+1, num_species):
        r_raw = raw_corr[i, j]
        r_stripped = stripped_corr[i, j]
        
        if abs(r_raw) > threshold_high:
            if abs(r_stripped) > threshold_high:
                direct.append((i, j, r_raw, r_stripped))
            else:
                indirect.append((i, j, r_raw, r_stripped))
        else:
            if abs(r_stripped) > threshold_high:
                emerging.append((i, j, r_raw, r_stripped))
            else:
                weak.append((i, j, r_raw, r_stripped))

total = len(direct) + len(indirect) + len(weak) + len(emerging)

print(f"\nInteraction type summary:")
print(f"  DIRECT: {len(direct)} pairs ({100*len(direct)/total:.1f}%)")
print(f"    → Likely true species interactions")
print(f"  INDIRECT: {len(indirect)} pairs ({100*len(indirect)/total:.1f}%)")
print(f"    → Nutrient-mediated associations")
print(f"  EMERGING: {len(emerging)} pairs ({100*len(emerging)/total:.1f}%)")
print(f"    → Interactions suppressed by nutrients")
print(f"  WEAK: {len(weak)} pairs ({100*len(weak)/total:.1f}%)")
print(f"    → Weak or no association")

# Show examples
print(f"\nTop 5 DIRECT interactions:")
direct_sorted = sorted(direct, key=lambda x: abs(x[3]), reverse=True)[:5]
for i, j, r_raw, r_stripped in direct_sorted:
    print(f"  Species {i:2d} ({species_names[i][:20]:20s}) ↔ "
          f"Species {j:2d} ({species_names[j][:20]:20s})  "
          f"raw={r_raw:+.3f}, stripped={r_stripped:+.3f}")

print(f"\nTop 5 INDIRECT (nutrient-mediated) associations:")
indirect_sorted = sorted(indirect, key=lambda x: abs(x[2]), reverse=True)[:5]
for i, j, r_raw, r_stripped in indirect_sorted:
    print(f"  Species {i:2d} ({species_names[i][:20]:20s}) ↔ "
          f"Species {j:2d} ({species_names[j][:20]:20s})  "
          f"raw={r_raw:+.3f}, stripped={r_stripped:+.3f}")

# ========================================
# Visualization - Nature style
# ========================================
mpl.rcParams['font.size'] = 8
mpl.rcParams['axes.linewidth'] = 0.5
mpl.rcParams['xtick.major.width'] = 0.5
mpl.rcParams['ytick.major.width'] = 0.5
mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['ytick.major.size'] = 2



# Create single figure and axis
fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=300)

# Plot the matrix
im = ax.imshow(delta, cmap='YlOrRd', vmin=0, vmax=np.nanpercentile(delta, 95), rasterized=True)
ax.set_title("c", loc='left', fontweight='bold', fontsize=10)
ax.set_xlabel("Species", fontsize=8)
ax.set_ylabel("Species", fontsize=8)
ax.tick_params(labelsize=6)

# Add colorbar
cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.ax.tick_params(labelsize=6)
cbar.set_label("|Δ|", fontsize=7)

plt.tight_layout()
plt.savefig(f"stripped_correlation_analysis/delta_matrix_{resolution}m.png", dpi=600, bbox_inches='tight')
plt.show()

# Scatter plot
fig, ax = plt.subplots(figsize=(3.35, 3.35), dpi=300)

triu_indices = np.triu_indices_from(raw_corr, k=1)
raw_upper = raw_corr[triu_indices]
stripped_upper = stripped_corr[triu_indices]

valid_mask = ~(np.isnan(raw_upper) | np.isnan(stripped_upper))
raw_upper = raw_upper[valid_mask]
stripped_upper = stripped_upper[valid_mask]

ax.scatter(raw_upper, stripped_upper, alpha=1, s=1, c='#4575b4', edgecolors='none', rasterized=True)

lim = max(abs(raw_upper.min()), abs(raw_upper.max()), 
          abs(stripped_upper.min()), abs(stripped_upper.max()))
ax.plot([-lim, lim], [-lim, lim], 'k-', linewidth=0.5, alpha=1)


ax.set_xlabel('Raw correlation', fontsize=10)
ax.set_ylabel('Correlation after regression', fontsize=10)
ax.tick_params(labelsize=6)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_aspect('equal')

plt.tight_layout()
plt.savefig(f"stripped_correlation_analysis/Fig_scatter_{resolution}m.svg", dpi=300, bbox_inches='tight')
plt.savefig(f"stripped_correlation_analysis/Fig_scatter_{resolution}m.png", dpi=600, bbox_inches='tight')

# ========================================
# Save Results to Files
# ========================================

print("\nSaving results to files...")

# Save matrices
np.save(f"stripped_correlation_analysis/raw_correlation_{resolution}m.npy", raw_corr)
np.save(f"stripped_correlation_analysis/stripped_correlation_{resolution}m.npy", stripped_corr)

# Save interaction classifications
with open(f"stripped_correlation_analysis/interaction_types_{resolution}m.txt", 'w') as f:
    f.write(f"INTERACTION TYPE CLASSIFICATION\n")
    f.write(f"Resolution: {resolution}m\n")
    f.write(f"Threshold for 'high' correlation: {threshold_high}\n")
    f.write(f"="*80 + "\n\n")
    
    f.write(f"DIRECT INTERACTIONS ({len(direct)} pairs)\n")
    f.write(f"High raw correlation, high stripped correlation → likely true species interaction\n")
    f.write("-"*80 + "\n")
    for i, j, r_raw, r_stripped in sorted(direct, key=lambda x: abs(x[3]), reverse=True):
        f.write(f"  {i:3d} ({species_names[i][:30]:30s}) ↔ "
               f"{j:3d} ({species_names[j][:30]:30s})  "
               f"raw={r_raw:+.3f}, stripped={r_stripped:+.3f}\n")
    
    f.write(f"\n\nINDIRECT (NUTRIENT-MEDIATED) ASSOCIATIONS ({len(indirect)} pairs)\n")
    f.write(f"High raw correlation, low stripped correlation → shared nutrient preferences\n")
    f.write("-"*80 + "\n")
    for i, j, r_raw, r_stripped in sorted(indirect, key=lambda x: abs(x[2]), reverse=True):
        f.write(f"  {i:3d} ({species_names[i][:30]:30s}) ↔ "
               f"{j:3d} ({species_names[j][:30]:30s})  "
               f"raw={r_raw:+.3f}, stripped={r_stripped:+.3f}\n")
    
    f.write(f"\n\nEMERGING INTERACTIONS ({len(emerging)} pairs)\n")
    f.write(f"Low raw correlation, high stripped correlation → nutrients suppress interaction\n")
    f.write("-"*80 + "\n")
    for i, j, r_raw, r_stripped in sorted(emerging, key=lambda x: abs(x[3]), reverse=True):
        f.write(f"  {i:3d} ({species_names[i][:30]:30s}) ↔ "
               f"{j:3d} ({species_names[j][:30]:30s})  "
               f"raw={r_raw:+.3f}, stripped={r_stripped:+.3f}\n")

print(f"✓ Saved interaction types to: stripped_correlation_analysis/interaction_types_{resolution}m.txt")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)
print(f"\nAll results saved to: stripped_correlation_analysis/")
print(f"  • Raw and stripped correlation matrices (.npy)")
print(f"  • Comparison visualization (.png)")
print(f"  • Scatter plot (.png)")
print(f"  • Interaction type classifications (.txt)")
print("="*80 + "\n")
