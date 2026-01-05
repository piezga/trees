import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cmocean
import os
from sklearn.linear_model import LinearRegression

from src.config import load_config
from src.compute import compute_spectra, compute_average_correlation, MarchenkoPastur

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

print(f"\n{'='*80}")
print(f"RAW vs PARTIAL CORRELATION ANALYSIS")
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
# Compute Partial Correlation Matrix
# ========================================

def compute_partial_correlation_matrix(species_abundance, nutrient_data):
    """
    Compute partial correlation matrix between species, controlling for nutrients.
    
    Parameters
    ----------
    species_abundance : array, shape (n_species, n_sites)
        Species abundance data (rows = species, cols = spatial bins)
    nutrient_data : array, shape (n_nutrients, n_sites)
        Nutrient abundance data (rows = nutrients, cols = spatial bins)
        
    Returns
    -------
    partial_corr : array, shape (n_species, n_species)
        Partial correlation matrix between species, controlling for nutrients
    """
    n_species = species_abundance.shape[0]
    
    # Transpose so rows are observations (sites) and columns are variables
    X_species = species_abundance.T  # (n_sites, n_species)
    X_nutrients = nutrient_data.T    # (n_sites, n_nutrients)
    
    # Store residuals after regressing out nutrients
    residuals = np.zeros_like(X_species)
    
    print("\nComputing partial correlations...")
    print(f"  Regressing out {X_nutrients.shape[1]} nutrients from {n_species} species...")
    
    # For each species, regress out the effect of all nutrients
    for i in range(n_species):
        if (i + 1) % 20 == 0:
            print(f"    Processed {i+1}/{n_species} species...")
        
        # Fit: species_i = β₀ + β₁*nutrient₁ + ... + βₖ*nutrientₖ + ε
        reg = LinearRegression()
        reg.fit(X_nutrients, X_species[:, i])
        
        # Residuals = what's left after removing nutrient effects
        residuals[:, i] = X_species[:, i] - reg.predict(X_nutrients)
    
    # Compute correlation matrix of residuals
    partial_corr = np.corrcoef(residuals.T)
    print(f"  ✓ Partial correlation matrix computed")
    
    return partial_corr


print("\nComputing raw correlation matrix...")
unfiltered_raw_corr = np.corrcoef(forest_abundance_avg)

unfiltered_partial_corr = compute_partial_correlation_matrix(forest_abundance_avg, nutrient_data)


# Filtering scatter_raw_vs_partial_
T = n_sites
raw_corr = MarchenkoPastur(unfiltered_raw_corr, 100, T)
partial_corr = MarchenkoPastur(unfiltered_partial_corr, 100, T)

print(f"  ✓ Raw correlation matrix shape: {raw_corr.shape}")
print(f"  ✓ Partial correlation matrix shape: {partial_corr.shape}")

# ========================================
# Compute Differences and Mediation
# ========================================

print("\nComputing differences and mediation indices...")

# Absolute difference
delta = np.abs(raw_corr - partial_corr)

# Mediation index: (raw - partial) / raw
# Positive values: correlation reduced after controlling for nutrients (nutrient-mediated)
# Negative values: correlation increased (nutrients were suppressing the relationship)
mediation_index = np.zeros_like(raw_corr)
mask = np.abs(raw_corr) > 1e-10
mediation_index[mask] = (raw_corr[mask] - partial_corr[mask]) / raw_corr[mask]

# Set diagonal to NaN (not meaningful)
np.fill_diagonal(delta, np.nan)
np.fill_diagonal(mediation_index, np.nan)
np.fill_diagonal(raw_corr, np.nan)
np.fill_diagonal(partial_corr, np.nan)

print(f"\nStatistics:")
print(f"  Mean |Δ|: {np.nanmean(delta):.4f}")
print(f"  Max |Δ|: {np.nanmax(delta):.4f}")
print(f"  Mean |mediation index|: {np.nanmean(np.abs(mediation_index)):.4f}")
print(f"  Median |mediation index|: {np.nanmedian(np.abs(mediation_index)):.4f}")

# ========================================
# Identify Interaction Types
# ========================================

print("\nClassifying species pair interactions...")

threshold_high = 0.3
threshold_low = 0.1

direct = []          # High raw, high partial → direct interaction
indirect = []        # High raw, low partial → nutrient-mediated
weak = []            # Low raw, low partial → weak/no interaction
emerging = []        # Low raw, high partial → suppressed by nutrients

for i in range(num_species):
    for j in range(i+1, num_species):
        r_raw = raw_corr[i, j]
        r_partial = partial_corr[i, j]
        
        if abs(r_raw) > threshold_high:
            if abs(r_partial) > threshold_high:
                direct.append((i, j, r_raw, r_partial))
            else:
                indirect.append((i, j, r_raw, r_partial))
        else:
            if abs(r_partial) > threshold_high:
                emerging.append((i, j, r_raw, r_partial))
            else:
                weak.append((i, j, r_raw, r_partial))

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
for i, j, r_raw, r_partial in direct_sorted:
    print(f"  Species {i:2d} ({species_names[i][:20]:20s}) ↔ "
          f"Species {j:2d} ({species_names[j][:20]:20s})  "
          f"raw={r_raw:+.3f}, partial={r_partial:+.3f}")

print(f"\nTop 5 INDIRECT (nutrient-mediated) associations:")
indirect_sorted = sorted(indirect, key=lambda x: abs(x[2]), reverse=True)[:5]
for i, j, r_raw, r_partial in indirect_sorted:
    print(f"  Species {i:2d} ({species_names[i][:20]:20s}) ↔ "
          f"Species {j:2d} ({species_names[j][:20]:20s})  "
          f"raw={r_raw:+.3f}, partial={r_partial:+.3f}")

# ========================================
# Visualization
# ========================================

print("\nGenerating visualization...")

os.makedirs("partial_correlation_analysis", exist_ok=True)

fig, axes = plt.subplots(2, 2, figsize=(16, 14), constrained_layout=True)

# Shared colormap for correlations
cmap = cmocean.cm.balance
vmin, vmax = -0.5, 0.5
norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

# Panel A: Raw correlation
im0 = axes[0, 0].imshow(raw_corr, cmap=cmap, norm=norm)
axes[0, 0].set_title("(A) Raw Species-Species Correlation", fontsize=14, fontweight='bold')
axes[0, 0].set_xlabel("Species index")
axes[0, 0].set_ylabel("Species index")
plt.colorbar(im0, ax=axes[0, 0], label="Correlation")

# Panel B: Partial correlation
im1 = axes[0, 1].imshow(partial_corr, cmap=cmap, norm=norm)
axes[0, 1].set_title("(B) Partial Correlation\n(controlling for nutrients)", 
                    fontsize=14, fontweight='bold')
axes[0, 1].set_xlabel("Species index")
axes[0, 1].set_ylabel("Species index")
plt.colorbar(im1, ax=axes[0, 1], label="Partial Correlation")

# Panel C: Absolute difference
im2 = axes[1, 0].imshow(delta, cmap='YlOrRd', vmin=0, vmax=np.nanpercentile(delta, 95))
axes[1, 0].set_title("(C) Absolute Difference |Raw - Partial|", 
                    fontsize=14, fontweight='bold')
axes[1, 0].set_xlabel("Species index")
axes[1, 0].set_ylabel("Species index")
plt.colorbar(im2, ax=axes[1, 0], label="|Δ|")

# Panel D: Mediation index
med_vmax = np.nanpercentile(np.abs(mediation_index), 95)
im3 = axes[1, 1].imshow(mediation_index, cmap='RdYlGn_r', 
                       vmin=-med_vmax, vmax=med_vmax)
axes[1, 1].set_title("(D) Mediation Index\n(Raw - Partial) / Raw", 
                    fontsize=14, fontweight='bold')
axes[1, 1].set_xlabel("Species index")
axes[1, 1].set_ylabel("Species index")
cbar = plt.colorbar(im3, ax=axes[1, 1], label="Mediation Index")
cbar.ax.text(0.5, 1.05, '+1 = fully nutrient-mediated', 
            transform=cbar.ax.transAxes, ha='center', fontsize=9)
cbar.ax.text(0.5, -0.05, '-1 = nutrients suppress correlation', 
            transform=cbar.ax.transAxes, ha='center', fontsize=9)

fig.suptitle(f"Raw vs. Partial Correlation Analysis (Resolution: {resolution}m)", 
            fontsize=16, fontweight='bold')

filename = f"partial_correlation_analysis/raw_vs_partial_comparison_{resolution}m.png"
plt.savefig(filename, dpi=300, bbox_inches='tight')
print(f"✓ Saved: {filename}")

plt.show()

# ========================================
# Additional Scatter Plot: Raw vs Partial
# ========================================

print("\nGenerating scatter plot...")

fig, ax = plt.subplots(figsize=(10, 10))

# Get upper triangle indices (excluding diagonal)
triu_indices = np.triu_indices_from(raw_corr, k=1)
raw_upper = raw_corr[triu_indices]
partial_upper = partial_corr[triu_indices]

# Remove NaN values
valid_mask = ~(np.isnan(raw_upper) | np.isnan(partial_upper))
raw_upper = raw_upper[valid_mask]
partial_upper = partial_upper[valid_mask]

# Scatter plot
ax.scatter(raw_upper, partial_upper, alpha=0.3, s=10, c='steelblue', edgecolors='none')

# 1:1 line
lim = max(abs(raw_upper.min()), abs(raw_upper.max()), 
          abs(partial_upper.min()), abs(partial_upper.max()))
ax.plot([-lim, lim], [-lim, lim], 'k--', alpha=0.5, label='1:1 line')

# Zero lines
ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
ax.axvline(0, color='gray', linestyle=':', alpha=0.5)

# Quadrant labels
ax.text(0.4, 0.4, 'Direct\ninteraction', transform=ax.transData, 
        ha='center', va='center', fontsize=12, alpha=0.6, 
        bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.3))

ax.text(0.4, -0.1, 'Nutrient-\nmediated', transform=ax.transData,
        ha='center', va='center', fontsize=12, alpha=0.6,
        bbox=dict(boxstyle='round', facecolor='orange', alpha=0.3))

ax.set_xlabel('Raw Correlation', fontsize=14, fontweight='bold')
ax.set_ylabel('Partial Correlation (controlling for nutrients)', fontsize=14, fontweight='bold')
ax.set_title(f'Raw vs. Partial Correlation Scatter\n(Resolution: {resolution}m)', 
            fontsize=16, fontweight='bold')
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
ax.set_aspect('equal')

filename_scatter = f"partial_correlation_analysis/scatter_raw_vs_partial_{resolution}m.png"
plt.savefig(filename_scatter, dpi=300, bbox_inches='tight')
print(f"✓ Saved: {filename_scatter}")

plt.show()

# ========================================
# Save Results to Files
# ========================================

print("\nSaving results to files...")

# Save matrices
np.save(f"partial_correlation_analysis/raw_correlation_{resolution}m.npy", raw_corr)
np.save(f"partial_correlation_analysis/partial_correlation_{resolution}m.npy", partial_corr)
np.save(f"partial_correlation_analysis/mediation_index_{resolution}m.npy", mediation_index)

# Save interaction classifications
with open(f"partial_correlation_analysis/interaction_types_{resolution}m.txt", 'w') as f:
    f.write(f"INTERACTION TYPE CLASSIFICATION\n")
    f.write(f"Resolution: {resolution}m\n")
    f.write(f"Threshold for 'high' correlation: {threshold_high}\n")
    f.write(f"="*80 + "\n\n")
    
    f.write(f"DIRECT INTERACTIONS ({len(direct)} pairs)\n")
    f.write(f"High raw correlation, high partial correlation → likely true species interaction\n")
    f.write("-"*80 + "\n")
    for i, j, r_raw, r_partial in sorted(direct, key=lambda x: abs(x[3]), reverse=True):
        f.write(f"  {i:3d} ({species_names[i][:30]:30s}) ↔ "
               f"{j:3d} ({species_names[j][:30]:30s})  "
               f"raw={r_raw:+.3f}, partial={r_partial:+.3f}\n")
    
    f.write(f"\n\nINDIRECT (NUTRIENT-MEDIATED) ASSOCIATIONS ({len(indirect)} pairs)\n")
    f.write(f"High raw correlation, low partial correlation → shared nutrient preferences\n")
    f.write("-"*80 + "\n")
    for i, j, r_raw, r_partial in sorted(indirect, key=lambda x: abs(x[2]), reverse=True):
        f.write(f"  {i:3d} ({species_names[i][:30]:30s}) ↔ "
               f"{j:3d} ({species_names[j][:30]:30s})  "
               f"raw={r_raw:+.3f}, partial={r_partial:+.3f}\n")
    
    f.write(f"\n\nEMERGING INTERACTIONS ({len(emerging)} pairs)\n")
    f.write(f"Low raw correlation, high partial correlation → nutrients suppress interaction\n")
    f.write("-"*80 + "\n")
    for i, j, r_raw, r_partial in sorted(emerging, key=lambda x: abs(x[3]), reverse=True):
        f.write(f"  {i:3d} ({species_names[i][:30]:30s}) ↔ "
               f"{j:3d} ({species_names[j][:30]:30s})  "
               f"raw={r_raw:+.3f}, partial={r_partial:+.3f}\n")

print(f"✓ Saved interaction types to: partial_correlation_analysis/interaction_types_{resolution}m.txt")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)
print(f"\nAll results saved to: partial_correlation_analysis/")
print(f"  • Raw and partial correlation matrices (.npy)")
print(f"  • Comparison visualization (.png)")
print(f"  • Scatter plot (.png)")
print(f"  • Interaction type classifications (.txt)")
print("="*80 + "\n")
