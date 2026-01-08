import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cmocean
import os
from sklearn.linear_model import LinearRegression
from sklearn.metrics import adjusted_rand_score, confusion_matrix

from src.config import load_config
from src.compute import (
    compute_spectra, compute_average_correlation, 
    MarchenkoPastur, L_detect_communities, compute_partial_correlation_matrix
)
from src.plotting import plot_community_dendrogram, plot_corr_with_communities

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
n_communities = 10
method = 'ward'
laplacian = True
tau = 1e-3

print(f"\n{'='*80}")
print(f"COMMUNITY DETECTION: RAW vs PARTIAL CORRELATION")
print(f"Resolution: {resolution}m")
print(f"Method: {'Laplacian diffusion' if laplacian else method}")
print(f"Target communities: {n_communities}")
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

# Load species names
names_file = f"{path_template.format(forest=forest)}{names_template.format(forest=forest, census=4)}"
with open(names_file, 'r', encoding='utf-8-sig') as f:
    species_names = [line.strip() for line in f.readlines()][:num_species]

# ========================================
# Load Nutrient Data
# ========================================

print("\nLoading nutrient data...")

nutrient_file = 'soil_data/barro_soil_data.xls'

if os.path.exists(nutrient_file):
    import pandas as pd
    nutrient_df = pd.read_excel(nutrient_file)
    nutrient_columns = ['Al', 'B', 'Ca', 'Cu', 'Fe', 'K', 'Mg', 'Mn', 'P', 'Zn', 'N', 'N(min)', 'pH']
    n_bins_x = bins[2]
    n_bins_y = bins[3]
    
    nutrient_df['x_bin'] = (nutrient_df['x'] / (1000 / n_bins_x)).astype(int).clip(0, n_bins_x - 1)
    nutrient_df['y_bin'] = (nutrient_df['y'] / (500 / n_bins_y)).astype(int).clip(0, n_bins_y - 1)
    
    nutrient_binned = nutrient_df.groupby(['x_bin', 'y_bin'])[nutrient_columns].mean()
    
    n_sites = n_bins_x * n_bins_y
    nutrient_data = np.zeros((len(nutrient_columns), n_sites))
    
    for idx, (x_bin, y_bin) in enumerate([(x, y) for x in range(n_bins_x) for y in range(n_bins_y)]):
        if (x_bin, y_bin) in nutrient_binned.index:
            nutrient_data[:, idx] = nutrient_binned.loc[(x_bin, y_bin)].values
        else:
            nutrient_data[:, idx] = np.nan
    
    for i in range(nutrient_data.shape[0]):
        col_mean = np.nanmean(nutrient_data[i, :])
        nutrient_data[i, np.isnan(nutrient_data[i, :])] = col_mean
    
    print(f"  Nutrient data shape: {nutrient_data.shape}")
    print(f"  Nutrients: {nutrient_columns}")

else:
    print(f"  Warning: {nutrient_file} not found")
    print(f"  Creating synthetic nutrient data for demonstration...")
    
    n_sites = forest_abundance_avg.shape[1]
    n_nutrients = 10
    
    np.random.seed(42)
    nutrient_data = np.random.randn(n_nutrients, n_sites)
    
    from scipy.ndimage import gaussian_filter1d
    for i in range(n_nutrients):
        nutrient_data[i, :] = gaussian_filter1d(nutrient_data[i, :], sigma=5)
    
    nutrient_data = (nutrient_data - nutrient_data.mean(axis=1, keepdims=True)) / nutrient_data.std(axis=1, keepdims=True)
    
    nutrient_columns = [f'Nutrient_{i+1}' for i in range(n_nutrients)]
    print(f"  Synthetic nutrient data shape: {nutrient_data.shape}")

# ========================================
# Compute Correlations
# ========================================


print("\nComputing correlation matrices...")
raw_corr = np.corrcoef(forest_abundance_avg)
print(f"  ✓ Raw correlation matrix: {raw_corr.shape}")

partial_corr = compute_partial_correlation_matrix(forest_abundance_avg, nutrient_data)
print(f"  ✓ Partial correlation matrix: {partial_corr.shape}")

# ========================================
# Apply Marchenko-Pastur Filter
# ========================================

print("\nApplying Marchenko-Pastur filter...")
T = bins[2] * bins[3]  # Number of spatial sites

filtered_raw = MarchenkoPastur(
    raw_corr, num_species, T, 
    remove_largest=False, remove_small=False
)
print(f"  ✓ Filtered raw correlation")

filtered_partial = MarchenkoPastur(
    partial_corr, num_species, T,
    remove_largest=False, remove_small=False
)
print(f"  ✓ Filtered partial correlation")

# ========================================
# Community Detection on RAW Correlation
# ========================================

print(f"\n{'─'*80}")
print("COMMUNITY DETECTION ON RAW CORRELATION")
print("(Niche-based communities: co-occurrence for any reason)")
print(f"{'─'*80}")

if laplacian:
    print(f"  Using Laplacian diffusion (τ={tau})...")
    (raw_reordered, raw_CM, raw_idx,
     raw_linkage, raw_cut_height) = L_detect_communities(
        filtered_raw, 
        n_communities=n_communities,
        tau=tau,
        return_linkage=True
    )
else:
    from src.compute import detect_communities_corr
    print(f"  Using {method} linkage...")
    (raw_reordered, raw_CM, raw_idx,
     raw_linkage, raw_cut_height) = detect_communities_corr(
        filtered_raw, 
        n_communities=n_communities,
        return_linkage=True,
        method=method
    )

print(f"  ✓ Communities detected: {len(set(raw_CM))}")
print(f"  ✓ Cut height: {raw_cut_height:.4f}")

# Count species per community
raw_community_sizes = {comm: np.sum(raw_CM == comm) for comm in set(raw_CM)}
print(f"  Community sizes: {dict(sorted(raw_community_sizes.items()))}")

# ========================================
# Community Detection on PARTIAL Correlation
# ========================================

print(f"\n{'─'*80}")
print("COMMUNITY DETECTION ON PARTIAL CORRELATION")
print("(Interaction-based communities: direct associations only)")
print(f"{'─'*80}")

if laplacian:
    print(f"  Using Laplacian diffusion (τ={tau})...")
    (partial_reordered, partial_CM, partial_idx,
     partial_linkage, partial_cut_height) = L_detect_communities(
        filtered_partial, 
        n_communities=n_communities,
        tau=tau,
        return_linkage=True
    )
else:
    print(f"  Using {method} linkage...")
    (partial_reordered, partial_CM, partial_idx,
     partial_linkage, partial_cut_height) = detect_communities_corr(
        filtered_partial, 
        n_communities=n_communities,
        return_linkage=True,
        method=method
    )

print(f"  ✓ Communities detected: {len(set(partial_CM))}")
print(f"  ✓ Cut height: {partial_cut_height:.4f}")

# Count species per community
partial_community_sizes = {comm: np.sum(partial_CM == comm) for comm in set(partial_CM)}
print(f"  Community sizes: {dict(sorted(partial_community_sizes.items()))}")

# ========================================
# Compare Communities
# ========================================

print(f"\n{'─'*80}")
print("COMMUNITY COMPARISON")
print(f"{'─'*80}\n")

# Adjusted Rand Index
ari = adjusted_rand_score(raw_CM, partial_CM)
print(f"Adjusted Rand Index (ARI): {ari:.3f}")
print(f"  1.0 = identical community structure")
print(f"  0.0 = random agreement")
print(f"  <0  = worse than random")

if ari > 0.7:
    print(f"\n  → HIGH AGREEMENT: Communities are mostly driven by DIRECT interactions")
    print(f"     (Nutrients don't strongly shape community structure)")
elif ari < 0.3:
    print(f"\n  → LOW AGREEMENT: Communities are mostly driven by NUTRIENT-MEDIATED associations")
    print(f"     (Shared habitat preferences drive co-occurrence)")
else:
    print(f"\n  → MODERATE AGREEMENT: Communities reflect a mix of both mechanisms")
    print(f"     (Both direct interactions and niche preferences matter)")

# Confusion matrix
conf_matrix = confusion_matrix(raw_CM, partial_CM)
print(f"\nConfusion Matrix (rows=raw communities, cols=partial communities):")
print(conf_matrix)

# ========================================
# Identify Stable vs Reorganized Species
# ========================================

print(f"\n{'─'*80}")
print("SPECIES REORGANIZATION ANALYSIS")
print(f"{'─'*80}\n")

# For each species, see if it's in the "same" community
# (i.e., most of its raw community-mates are also in its partial community)

def get_community_overlap(species_idx, CM1, CM2):
    """
    For a given species, calculate what fraction of its community-mates
    in CM1 are also community-mates in CM2.
    """
    comm1 = CM1[species_idx]
    comm2 = CM2[species_idx]
    
    # Who are the community-mates in CM1?
    mates_in_CM1 = set(np.where(CM1 == comm1)[0])
    mates_in_CM1.discard(species_idx)  # Remove self
    
    # Who are the community-mates in CM2?
    mates_in_CM2 = set(np.where(CM2 == comm2)[0])
    mates_in_CM2.discard(species_idx)
    
    # Overlap
    overlap = mates_in_CM1 & mates_in_CM2
    
    if len(mates_in_CM1) == 0:
        return 1.0  # Singleton in both
    
    return len(overlap) / len(mates_in_CM1)

species_stability = []
for i in range(num_species):
    overlap = get_community_overlap(i, raw_CM, partial_CM)
    species_stability.append((i, overlap, raw_CM[i], partial_CM[i]))

# Sort by stability (most stable first)
species_stability.sort(key=lambda x: x[1], reverse=True)

print("Most STABLE species (retain community-mates after controlling for nutrients):")
for i, overlap, raw_c, partial_c in species_stability[:10]:
    print(f"  Species {i:2d} ({species_names[i][:30]:30s})  "
          f"overlap={overlap:.2f}  raw_comm={raw_c}, partial_comm={partial_c}")

print("\nMost REORGANIZED species (different community-mates after controlling for nutrients):")
for i, overlap, raw_c, partial_c in species_stability[-10:]:
    print(f"  Species {i:2d} ({species_names[i][:30]:30s})  "
          f"overlap={overlap:.2f}  raw_comm={raw_c}, partial_comm={partial_c}")

# ========================================
# Save Results
# ========================================

output_dir = "community_raw_vs_partial"
os.makedirs(output_dir, exist_ok=True)

print(f"\n{'─'*80}")
print("SAVING RESULTS")
print(f"{'─'*80}\n")

# Save community assignments
np.save(f"{output_dir}/raw_communities_{resolution}m.npy", raw_CM)
np.save(f"{output_dir}/partial_communities_{resolution}m.npy", partial_CM)

# Save as text files
for comm in sorted(set(raw_CM)):
    members = np.where(raw_CM == comm)[0]
    with open(f"{output_dir}/raw_community_{comm}.txt", 'w') as f:
        f.write(f"RAW Correlation Community {comm} ({len(members)} species)\n")
        f.write("="*60 + "\n\n")
        for idx in members:
            f.write(f"{idx:3d}  {species_names[idx]}\n")

for comm in sorted(set(partial_CM)):
    members = np.where(partial_CM == comm)[0]
    with open(f"{output_dir}/partial_community_{comm}.txt", 'w') as f:
        f.write(f"PARTIAL Correlation Community {comm} ({len(members)} species)\n")
        f.write("="*60 + "\n\n")
        for idx in members:
            f.write(f"{idx:3d}  {species_names[idx]}\n")

print(f"✓ Saved community assignments to text files")

# ========================================
# Visualize Dendrograms
# ========================================

print("\nGenerating dendrograms...")

plot_community_dendrogram(
    raw_linkage,
    threshold=raw_cut_height,
    title=f'Raw Correlation Communities (Resolution {resolution}m)',
    filename=f'{output_dir}/dendrogram_raw_{resolution}m.png',
    show=False
)
print(f"✓ Saved raw correlation dendrogram")

plot_community_dendrogram(
    partial_linkage,
    threshold=partial_cut_height,
    title=f'Partial Correlation Communities (Resolution {resolution}m)',
    filename=f'{output_dir}/dendrogram_partial_{resolution}m.png',
    show=False
)
print(f"✓ Saved partial correlation dendrogram")

# ========================================
# Visualize Reordered Correlation Matrices
# ========================================

print("\nGenerating reordered correlation matrix plots...")

plot_corr_with_communities(
    filtered_raw,
    raw_CM,
    outfile=f"{output_dir}/raw_corr_with_communities_{resolution}m.png",
    title="Raw Correlation (Niche-based Communities)",
    show=False
)
print(f"✓ Saved raw correlation matrix with communities")

plot_corr_with_communities(
    filtered_partial,
    partial_CM,
    outfile=f"{output_dir}/partial_corr_with_communities_{resolution}m.png",
    title="Partial Correlation (Interaction-based Communities)",
    show=False
)
print(f"✓ Saved partial correlation matrix with communities")

# ========================================
# Side-by-side Comparison
# ========================================

print("\nGenerating side-by-side comparison...")

fig, axes = plt.subplots(1, 2, figsize=(16, 7), constrained_layout=True)

cmap = cmocean.cm.balance
vmin, vmax = -0.5, 0.5
norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

# Raw correlation with communities
im1 = axes[0].imshow(raw_reordered, cmap=cmap, norm=norm)
axes[0].set_title("(A) Raw Correlation\n(Niche-based Communities)", 
                 fontsize=14, fontweight='bold')
axes[0].set_xlabel("Species (reordered)")
axes[0].set_ylabel("Species (reordered)")

# Add community boundaries
def get_community_bounds(CM_sorted):
    bounds = [0]
    current = CM_sorted[0]
    for i, c in enumerate(CM_sorted):
        if c != current:
            bounds.append(i)
            current = c
    bounds.append(len(CM_sorted))
    return bounds

raw_bounds = get_community_bounds(raw_CM[raw_idx])
for i in range(len(raw_bounds)-1):
    start, end = raw_bounds[i], raw_bounds[i+1]
    rect = plt.Rectangle((start-0.5, start-0.5), end-start, end-start,
                         fill=False, edgecolor='black', linewidth=1.5, linestyle='--')
    axes[0].add_patch(rect)

# Partial correlation with communities
im2 = axes[1].imshow(partial_reordered, cmap=cmap, norm=norm)
axes[1].set_title("(B) Partial Correlation\n(Interaction-based Communities)", 
                 fontsize=14, fontweight='bold')
axes[1].set_xlabel("Species (reordered)")
axes[1].set_ylabel("Species (reordered)")

partial_bounds = get_community_bounds(partial_CM[partial_idx])
for i in range(len(partial_bounds)-1):
    start, end = partial_bounds[i], partial_bounds[i+1]
    rect = plt.Rectangle((start-0.5, start-0.5), end-start, end-start,
                         fill=False, edgecolor='black', linewidth=1.5, linestyle='--')
    axes[1].add_patch(rect)

# Shared colorbar
fig.colorbar(im2, ax=axes, orientation='horizontal', 
            fraction=0.046, pad=0.08, label='Correlation')

fig.suptitle(f'Community Structure: Raw vs. Partial Correlation (ARI = {ari:.3f})', 
            fontsize=16, fontweight='bold')

plt.savefig(f'{output_dir}/comparison_raw_vs_partial_{resolution}m.png', 
           dpi=300, bbox_inches='tight')
print(f"✓ Saved side-by-side comparison")

plt.show()

# ========================================
# Confusion Matrix Visualization
# ========================================

print("\nGenerating confusion matrix...")

fig, ax = plt.subplots(figsize=(10, 8))

im = ax.imshow(conf_matrix, cmap='Blues', aspect='auto')

# Add text annotations
for i in range(conf_matrix.shape[0]):
    for j in range(conf_matrix.shape[1]):
        text = ax.text(j, i, conf_matrix[i, j],
                      ha="center", va="center", 
                      color="white" if conf_matrix[i, j] > conf_matrix.max()/2 else "black",
                      fontsize=10)

ax.set_xlabel("Partial Correlation Communities", fontsize=12, fontweight='bold')
ax.set_ylabel("Raw Correlation Communities", fontsize=12, fontweight='bold')
ax.set_title(f"Community Overlap Matrix (ARI = {ari:.3f})", 
            fontsize=14, fontweight='bold')

ax.set_xticks(range(n_communities))
ax.set_yticks(range(n_communities))
ax.set_xticklabels(range(1, n_communities+1))
ax.set_yticklabels(range(1, n_communities+1))

plt.colorbar(im, ax=ax, label="Number of species")
plt.tight_layout()

plt.savefig(f'{output_dir}/confusion_matrix_{resolution}m.png', 
           dpi=300, bbox_inches='tight')
print(f"✓ Saved confusion matrix")

plt.show()

# ========================================
# Summary Report
# ========================================

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)

print(f"\nKey Findings:")
print(f"  • Adjusted Rand Index: {ari:.3f}")
print(f"  • Communities are {'SIMILAR' if ari > 0.5 else 'DIFFERENT'} between raw and partial")
print(f"  • Raw communities: {len(set(raw_CM))} detected")
print(f"  • Partial communities: {len(set(partial_CM))} detected")

print(f"\nAll results saved to: {output_dir}/")
print(f"  • Community assignments (.npy and .txt)")
print(f"  • Dendrograms")
print(f"  • Reordered correlation matrices")
print(f"  • Side-by-side comparison")
print(f"  • Confusion matrix")
print("="*80 + "\n")
