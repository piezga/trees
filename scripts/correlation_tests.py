import numpy as np
import argparse
import matplotlib
import matplotlib.pyplot as plt
import os

from src.config import load_config
from src.compute import (
        compute_spectra, MarchenkoPastur, L_detect_communities,
        marchenko_pastur_pdf, detect_communities_corr, compute_community_overlap,
        compute_average_correlation
        )
from src.plotting import (
    plot_correlation_matrices_comparison, plot_community_dendrogram,
    plot_community_confusion_matrix, plot_correlation_stability
    )
from src.utils import(
        find_consensus_communities
        )
#matplotlib.use('TkAgg')  # Switch from Qt to Tkinter backend

# === Load config ===
config = load_config()
num_species = config['analysis']['num_species']
censuses = config['forests']['censuses']
path_template = config['forests']['templates']['path_template']
census_template = config['forests']['templates']['census_template']
names_template = config['forests']['templates']['names_template']
# === Pick resolution ===
rng = np.random.default_rng()
resolution = int(rng.integers(5,33))
resolution = 10
print(f'Testing for resolution {resolution} m ')

(senm_mean, senm_std, forest_spectra, 
 bins, senm_abundance, forest_abundances) = compute_spectra(resolution, 
                                                           calculate=True)

# Check standardization

print('\nSENM STANDARDIZATION')
print('\nSENM mean')
print(np.mean(senm_abundance, axis = 1))
print('\nSENM std')
print(np.std(senm_abundance, axis = 1))

print('\nFOREST STANDARDIZATION')
print('\nForest mean')
print(np.mean(forest_abundances[2], axis = 1))
print('\nForest std')
print(np.std(forest_abundances[2], axis = 1))

# Is the shape right?
print('\n\n')
print('Testing for correct shape (column)')
print(f'Forest shape is : {np.shape(forest_abundances)}')
print(f'Expecting column: {500*1000/(resolution**2)}')
print(f'SENM shape is   : {np.shape(senm_abundance)}')
print(f'Expecting column: {500*500/(resolution**2)}')

# === Calculating correlation matrices ===

senm_corr= np.corrcoef(senm_abundance)

(forest_corr, 
 forest_corr_std, forest_corrs_all) = compute_average_correlation(forest_abundances)

print('\n\n')
print('The correlation matrices have shape')
print(f'Forest: {forest_corr.shape}')
print(f'SENM: {senm_corr.shape}')
print(f'\nForest preview:\n {forest_corr}')
print(f'\nSENM preview:\n {senm_corr}')

# === MarchenkoPastur Check ===
print('\n\n')
print('Here I am keeping the small eigenvalues:')
filtered_senm_corr = MarchenkoPastur(senm_corr, num_species, 
                                     bins[0]*bins[1], remove_largest=False, 
                                     remove_small=False )
filtered_forest_corr = MarchenkoPastur(forest_corr, num_species, 
                                       bins[2]*bins[3], remove_largest=False,
                                       remove_small=False)
print(f'\nSenm: \n {filtered_senm_corr}')
print(f'\nForest: \n {filtered_forest_corr}')

print('\n\n')
print('Here I am removing the small eigenvalues:')
filtered_senm_corr2 = MarchenkoPastur(senm_corr, num_species, 
                                     bins[0]*bins[1], remove_largest=False, 
                                     remove_small=True )
filtered_forest_corr2 = MarchenkoPastur(forest_corr, num_species, 
                                       bins[2]*bins[3], remove_largest=False,
                                       remove_small=True)
print(f'\nSenm: \n {filtered_senm_corr2}')
print(f'\nForest: \n {filtered_forest_corr2}')

# === Plotting the difference ===

# Create the plot with 2 rows and 2 columns
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
# Making sure I use the same scale
vmin = -1
vmax = 1
# Plot the original SENM correlation matrix
axes[0, 0].imshow(filtered_senm_corr, cmap='coolwarm', aspect='auto',vmin = vmin, vmax = vmax)
axes[0, 0].set_title('Filtered SENM (small eigenvalues included)')
axes[0, 0].set_xlabel('Index')
axes[0, 0].set_ylabel('Index')
plt.colorbar(axes[0, 0].imshow(filtered_senm_corr, cmap='coolwarm', aspect='auto',vmin = vmin, vmax = vmax), ax=axes[0, 0])

# Plot the original Forest correlation matrix
axes[0, 1].imshow(filtered_forest_corr, cmap='coolwarm', aspect='auto',vmin = vmin, vmax = vmax)
axes[0, 1].set_title('Filtered Forest (small eigenvalues included)')
axes[0, 1].set_xlabel('Index')
axes[0, 1].set_ylabel('Index')
plt.colorbar(axes[0, 1].imshow(filtered_forest_corr, cmap='coolwarm', aspect='auto',vmin = vmin, vmax = vmax), ax=axes[0, 1])

# Plot the filtered SENM correlation matrix (no small eigenvalues)
axes[1, 0].imshow(filtered_senm_corr2, cmap='coolwarm', aspect='auto',vmin = vmin, vmax = vmax)
axes[1, 0].set_title('Filtered SENM (no small eigenvalues)')
axes[1, 0].set_xlabel('Index')
axes[1, 0].set_ylabel('Index')
plt.colorbar(axes[1, 0].imshow(filtered_senm_corr2, cmap='coolwarm', aspect='auto',vmin = vmin, vmax = vmax), ax=axes[1, 0])

# Plot the filtered Forest correlation matrix (no small eigenvalues)
axes[1, 1].imshow(filtered_forest_corr2, cmap='coolwarm', aspect='auto',vmin = vmin, vmax = vmax)
axes[1, 1].set_title('Filtered Forest (no small eigenvalues)')
axes[1, 1].set_xlabel('Index')
axes[1, 1].set_ylabel('Index')
plt.colorbar(axes[1, 1].imshow(filtered_forest_corr2, cmap='coolwarm', aspect='auto',vmin = vmin, vmax = vmax), ax=axes[1, 1])

# Adjust layout to make sure the plots are spaced nicely
plt.tight_layout()
#plt.show()


# Ensure output directory exists
plot_dir = "plots/spectral_analysis"
os.makedirs(plot_dir, exist_ok=True)

# Four matrices we will analyze
matrices = {
    "filtered_senm_keep_small": filtered_senm_corr,
    "filtered_forest_keep_small": filtered_forest_corr,
    "filtered_senm_remove_small": filtered_senm_corr2,
    "filtered_forest_remove_small": filtered_forest_corr2,
}

# MP parameters for SENM and Forest
T_senm = bins[0] * bins[1]
T_forest = bins[2] * bins[3]
q_senm = num_species / T_senm
q_forest = num_species / T_forest

# ============================
#   LOOP THROUGH MATRICES
# ============================

for name, M in matrices.items():

    print(f"\n\n=== ANALYZING {name} ===")

    # Eigenvalues & eigenvectors
    eigvals, eigvecs = np.linalg.eigh(M)

    # Basic stats
    trace_eigenvalues = np.sum(eigvals)
    # re-Compute the trace by summing over the diagonal of the matrix
    trace_diagonal = np.sum(np.diag(M))

    # Check that both methods give the same result
    print(f"Trace (from eigenvalues): {trace_eigenvalues}")
    print(f"Trace (from diagonal elements): {trace_diagonal}")

    fr_norm = np.linalg.norm(M, "fro")
    eff_rank = np.exp(-np.sum((eigvals / np.sum(eigvals)) * np.log(eigvals / np.sum(eigvals) + 1e-12)))
    cond_num = eigvals.max() / (eigvals.min() + 1e-12)

    print(f"Frobenius norm: {fr_norm:.4f}")
    print(f"Effective rank: {eff_rank:.4f}")
    print(f"Condition number: {cond_num:.4e}")

    # Determine MP parameters
    if "senm" in name:
        q = q_senm
    else:
        q = q_forest

    # MP support
    lambda_min = (1 - np.sqrt(q))**2
    lambda_max = (1 + np.sqrt(q))**2

    # ====================
    # 1. Scree plot
    # ====================
    plt.figure(figsize=(8, 5))
    plt.plot(np.sort(eigvals)[::-1], marker='.')
    plt.axhline(lambda_min, color='r', linestyle='--', label='MP lower bound')
    plt.axhline(lambda_max, color='g', linestyle='--', label='MP upper bound')
    plt.title(f"Scree Plot: {name}")
    plt.xlabel("Eigenvalue index")
    plt.ylabel("Eigenvalue")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/{name}_scree.png", dpi=200)
    plt.close()

    # ====================
    # 2. Histogram + MP curve
    # ====================
    plt.figure(figsize=(8, 5))
    bins_hist = max(20, int(len(eigvals)/15))
    plt.hist(eigvals, bins=bins_hist, density=True, alpha=0.5, label="Empirical eigvals")

    # Create a grid for eigenvalues
    lgrid = np.linspace(lambda_min, lambda_max, 2000)

    # Call the Marchenko-Pastur PDF function
    mp_density, l_min, l_max = marchenko_pastur_pdf(lgrid, q)

    # Plot the Marchenko-Pastur density
    plt.plot(lgrid, mp_density, label="Marchenko-Pastur PDF", color='b')

    # Plot the MP bounds (just in case you want to show them too)
    plt.axvline(l_min, color='r', linestyle='--', 
                label=f'MP lower bound (l_min={l_min:.3f})')
    plt.axvline(l_max, color='g', linestyle='--', 
                label=f'MP upper bound (l_max={l_max:.3f})')

    plt.title(f"Histogram + Marchenko-Pastur: {name}")
    plt.xlabel("Eigenvalue")
    plt.ylabel("Density")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/{name}_hist_mp.png", dpi=200)
    plt.close()

    # ====================
    # 3. Eigenvalues outside MP bounds
    # ====================
    outliers = eigvals[(eigvals < lambda_min) | (eigvals > lambda_max)]
    above = eigvals[ eigvals > lambda_max]
    print(f"Number of outliers outside MP bounds: {len(outliers)}")
    print("Should be 100!")
    print(f"Number of outliers above MP bounds: {len(above)}")
    # Saving the n_communities for later
    if 'senm_keep' in name:
        n_comms_senm = len(above)
    elif 'forest_keep' in name:
        n_comms_forest = len(above)

    # ====================
    # 4. Leading eigenvector plot
    # ====================
    leading_vec = eigvecs[:, -1]

    plt.figure(figsize=(8, 5))
    plt.plot(leading_vec, marker='.')
    plt.title(f"Leading Eigenvector Components: {name}")
    plt.xlabel("Component index")
    plt.ylabel("Value")
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/{name}_leading_eigenvector.png", dpi=200)
    plt.close()

    # ====================
    # 5. Reconstruction error (after dropping small-eigs)
    # ====================
    # Reconstruct the matrix from eigen decomposition
    M_reconstructed = eigvecs @ np.diag(eigvals) @ eigvecs.T
    reconstruction_error = np.linalg.norm(M - M_reconstructed, "fro")

    print(f"Reconstruction Frobenius error: {reconstruction_error:.4e}")

# ================================
# 6. Compare eigenvalue spectra across all matrices
# ================================

plt.figure(figsize=(10, 6))
for name, M in matrices.items():
    eigvals = np.linalg.eigh(M)[0]
    plt.plot(np.sort(eigvals)[::-1], label=name)

plt.title("Eigenvalue Spectra Comparison (All Filtered Matrices)")
plt.xlabel("Eigenvalue index")
plt.ylabel("Eigenvalue")
plt.legend()
plt.tight_layout()
plt.savefig(f"{plot_dir}/comparison_all_scree.png", dpi=200)
plt.close()

print(f"\nAll spectral plots saved to: {plot_dir}\n")

# === Community analysis | Macmahon&Garlaschelli ===

# The matrices are already filtered (MP)
# We are also considering the small eigenvalues

senm_matrix = filtered_senm_corr
forest_matrix = filtered_forest_corr


# Here we turn the filtered matrix into a distance matrix forcing N communities
print(f'We are looking for {n_comms_senm} communities in SENM')
print(f'We are looking for {n_comms_forest} communities in FOREST')

# === Used methods===
methods = ['ward','average', 'single', 'complete', 'weighted']

senm_communities = {}
forest_communities = {}


for method in methods:
# SENM community detection
    
    print(f'Testing method: {method}')
    (senm_reordered, senm_fcluster, senm_idx,
     senm_linkage_matrix, 
     senm_cut_height) = detect_communities_corr(senm_matrix,  
                                                n_comms_senm, return_linkage = True,
                                                method = method)
    senm_communities[method] = senm_fcluster
    print(f"\n--- SENM Community Detection Results ---")
    print(f"• Communities detected in the SENM matrix: {n_comms_senm}")
    print(f"• Reordered SENM matrix shape: {senm_reordered.shape}")
    print(f"• Number of clusters in SENM: {len(set(senm_fcluster))} unique clusters detected.")
    print(f"• Community detection linkage matrix (SENM) shape: {senm_linkage_matrix.shape}")
    print(f"• Cutting height for the dendrogram: {senm_cut_height:.3f}")
    print(f"--------------------------------------------\n")

# Plot the dendrogram
    plot_community_dendrogram(senm_linkage_matrix,threshold=senm_cut_height,
                              title = 'SENM + ' + method,
                              filename=f'dendrograms/senm_{method}.png')
# Forest community detection

    (forest_reordered, forest_fcluster, forest_idx,
     forest_linkage_matrix, 
     forest_cut_height) = detect_communities_corr(forest_matrix, n_comms_forest,
                                                  return_linkage = True,
                                                  method = method)
    forest_communities[method] = forest_fcluster
    # After community detection for FOREST
    print(f"\n--- FOREST Community Detection Results ---")
    print(f"• Communities detected in the FOREST matrix: {n_comms_forest}")
    print(f"• Reordered FOREST matrix shape: {forest_reordered.shape}")
    print(f"• Number of clusters in FOREST: {len(set(forest_fcluster))} unique clusters detected.")
    print(f"• Community detection linkage matrix (FOREST) shape: {forest_linkage_matrix.shape}")
    print(f"• Cutting height for the dendrogram: {forest_cut_height:.3f}")
    print(f"-------------------------------------------\n")

    plot_community_dendrogram(forest_linkage_matrix,threshold=forest_cut_height,
                              title = 'Forest + ' + method,
                              filename=f'dendrograms/forest_{method}.png')

    # Create output directory
os.makedirs("community_comparison", exist_ok=True)

# SENM: Compare all pairs of methods
print("\n--- SENM: Pairwise Method Comparisons ---")
senm_ari_matrix = np.zeros((len(methods), len(methods)))

for i, method1 in enumerate(methods):
    for j, method2 in enumerate(methods):
        if i <= j:  # Only compute upper triangle + diagonal
            ari, _ = compute_community_overlap(
                senm_communities[method1],
                senm_communities[method2]
            )
            senm_ari_matrix[i, j] = ari
            senm_ari_matrix[j, i] = ari  # Symmetric
            
            if i < j:  # Don't plot diagonal
                print(f"\n{method1} vs {method2}:")
                print(f"  ARI = {ari:.3f}")
                
                plot_community_confusion_matrix(
                    senm_communities[method1],
                    senm_communities[method2],
                    method1_name=method1,
                    method2_name=method2,
                    filename=f"community_comparison/senm_{method1}_vs_{method2}.png",
                    show=False
                )

# FOREST: Compare all pairs of methods
print("\n--- FOREST: Pairwise Method Comparisons ---")
forest_ari_matrix = np.zeros((len(methods), len(methods)))

for i, method1 in enumerate(methods):
    for j, method2 in enumerate(methods):
        if i <= j:
            ari, _ = compute_community_overlap(
                forest_communities[method1],
                forest_communities[method2]
            )
            forest_ari_matrix[i, j] = ari
            forest_ari_matrix[j, i] = ari
            
            if i < j:
                print(f"\n{method1} vs {method2}:")
                print(f"  ARI = {ari:.3f}")
                
                plot_community_confusion_matrix(
                    forest_communities[method1],
                    forest_communities[method2],
                    method1_name=method1,
                    method2_name=method2,
                    filename=f"community_comparison/forest_{method1}_vs_{method2}.png",
                    show=False
                )
# ========================================
# 2. ARI Heatmaps
# ========================================

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# SENM ARI heatmap
im1 = axes[0].imshow(senm_ari_matrix, cmap='RdYlGn', vmin=0, vmax=1)
axes[0].set_xticks(range(len(methods)))
axes[0].set_yticks(range(len(methods)))
axes[0].set_xticklabels(methods, rotation=45, ha='right')
axes[0].set_yticklabels(methods)
axes[0].set_title('SENM: Method Agreement (ARI)')

# Add text annotations
for i in range(len(methods)):
    for j in range(len(methods)):
        text = axes[0].text(j, i, f'{senm_ari_matrix[i, j]:.2f}',
                           ha="center", va="center", 
                           color="white" if senm_ari_matrix[i, j] < 0.5 else "black",
                           fontsize=10)

# FOREST ARI heatmap
im2 = axes[1].imshow(forest_ari_matrix, cmap='RdYlGn', vmin=0, vmax=1)
axes[1].set_xticks(range(len(methods)))
axes[1].set_yticks(range(len(methods)))
axes[1].set_xticklabels(methods, rotation=45, ha='right')
axes[1].set_yticklabels(methods)
axes[1].set_title('FOREST: Method Agreement (ARI)')

# Add text annotations
for i in range(len(methods)):
    for j in range(len(methods)):
        text = axes[1].text(j, i, f'{forest_ari_matrix[i, j]:.2f}',
                           ha="center", va="center",
                           color="white" if forest_ari_matrix[i, j] < 0.5 else "black",
                           fontsize=10)

# Shared colorbar
fig.colorbar(im2, ax=axes, label='Adjusted Rand Index', fraction=0.046, pad=0.04)
plt.tight_layout()
plt.savefig('community_comparison/ari_heatmaps.png', dpi=300, bbox_inches='tight')
plt.close()

# ========================================
# 3. Community Size Distribution
# ========================================

fig, axes = plt.subplots(2, len(methods), figsize=(4*len(methods), 8))

for idx, method in enumerate(methods):
    # SENM
    senm_sizes = np.bincount(senm_communities[method])[1:]  # Skip 0 if present
    axes[0, idx].bar(range(1, len(senm_sizes)+1), senm_sizes)
    axes[0, idx].set_title(f'SENM - {method}')
    axes[0, idx].set_xlabel('Community ID')
    axes[0, idx].set_ylabel('Size (# species)')
    axes[0, idx].grid(alpha=0.3)
    
    # FOREST
    forest_sizes = np.bincount(forest_communities[method])[1:]
    axes[1, idx].bar(range(1, len(forest_sizes)+1), forest_sizes)
    axes[1, idx].set_title(f'FOREST - {method}')
    axes[1, idx].set_xlabel('Community ID')
    axes[1, idx].set_ylabel('Size (# species)')
    axes[1, idx].grid(alpha=0.3)

plt.tight_layout()
plt.savefig('community_comparison/community_sizes.png', dpi=300, bbox_inches='tight')
plt.close()

# ========================================
# 4. Consensus Community Analysis
# ========================================

# Find consensus
senm_consensus, senm_co_cluster = find_consensus_communities(senm_communities, threshold=0.7)
forest_consensus, forest_co_cluster = find_consensus_communities(forest_communities, threshold=0.7)

# Plot consensus matrices
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# SENM co-clustering frequency
im1 = axes[0, 0].imshow(senm_co_cluster, cmap='YlOrRd', vmin=0, vmax=1)
axes[0, 0].set_title('SENM: Co-clustering Frequency')
axes[0, 0].set_xlabel('Species index')
axes[0, 0].set_ylabel('Species index')
plt.colorbar(im1, ax=axes[0, 0], label='Fraction of methods')

# SENM consensus (70% threshold)
im2 = axes[0, 1].imshow(senm_consensus, cmap='binary', vmin=0, vmax=1)
axes[0, 1].set_title('SENM: Consensus Communities (≥70% agreement)')
axes[0, 1].set_xlabel('Species index')
axes[0, 1].set_ylabel('Species index')
plt.colorbar(im2, ax=axes[0, 1])

# FOREST co-clustering frequency
im3 = axes[1, 0].imshow(forest_co_cluster, cmap='YlOrRd', vmin=0, vmax=1)
axes[1, 0].set_title('FOREST: Co-clustering Frequency')
axes[1, 0].set_xlabel('Species index')
axes[1, 0].set_ylabel('Species index')
plt.colorbar(im3, ax=axes[1, 0], label='Fraction of methods')

# FOREST consensus (70% threshold)
im4 = axes[1, 1].imshow(forest_consensus, cmap='binary', vmin=0, vmax=1)
axes[1, 1].set_title('FOREST: Consensus Communities (≥70% agreement)')
axes[1, 1].set_xlabel('Species index')
axes[1, 1].set_ylabel('Species index')
plt.colorbar(im4, ax=axes[1, 1])

plt.tight_layout()
plt.savefig('community_comparison/consensus_communities.png', dpi=300, bbox_inches='tight')
plt.close()

# ========================================
# 5. Summary Statistics
# ========================================

print("\n" + "="*80)
print("SUMMARY STATISTICS")
print("="*80)

print("\n--- SENM ---")
for method in methods:
    CM = senm_communities[method]
    sizes = np.bincount(CM)[1:]

    size_dist = {i: int(s) for i, s in enumerate(sizes, start=1)}

    print(f"\n{method}:")
    print(f"  Communities: {len(sizes)}")
    print(f"  Size range: {sizes.min()} - {sizes.max()}")
    print(f"  Mean size: {sizes.mean():.1f} ± {sizes.std():.1f}")
    print(f"  Size distribution: {size_dist}")

print("\n--- FOREST ---")
for method in methods:
    CM = forest_communities[method]
    sizes = np.bincount(CM)[1:]

    size_dist = {i: int(s) for i, s in enumerate(sizes, start=1)}

    print(f"\n{method}:")
    print(f"  Communities: {len(sizes)}")
    print(f"  Size range: {sizes.min()} - {sizes.max()}")
    print(f"  Mean size: {sizes.mean():.1f} ± {sizes.std():.1f}")
    print(f"  Size distribution: {size_dist}")

print("\n" + "="*80)
print(f"All comparison plots saved to: community_comparison/")
print("="*80 + "\n")


# ========================================
# Ward Method: Detailed Species Analysis
# ========================================

print("\n" + "="*80)
print("WARD CLUSTERING: DETAILED SPECIES COMPOSITION")
print("="*80 + "\n")

# Run Ward clustering
method = 'ward'

# SENM
(senm_reordered, senm_CM, senm_idx,
 senm_linkage, senm_cut_height) = detect_communities_corr(
    senm_matrix, n_comms_senm, 
    return_linkage=True, method=method
)

# FOREST
(forest_reordered, forest_CM, forest_idx,
 forest_linkage, forest_cut_height) = detect_communities_corr(
    forest_matrix, n_comms_forest, 
    return_linkage=True, method=method
)

plot_correlation_stability(
    forest_corrs_all,
    forest_corr,
    forest_corr_std,
    censuses,
    filename='community_comparison/forest_correlations_across_censuses.png'
)

print("\n" + "="*80)
print("CORRELATION ANALYSIS ACROSS CENSUSES")
print("="*80)
print(f"\nComputed correlations for {len(forest_corrs_all)} censuses")
print(f"Using species list from census 4 (reference)")
print(f"\nAverage correlation matrix shape: {forest_corr.shape}")
print(f"Correlation stability (mean std): {np.mean(forest_corr_std):.3f}")
print(f"Max correlation std: {np.max(forest_corr_std):.3f}")
print(f"Species pairs with stable correlations (std < 0.1): {np.sum(forest_corr_std < 0.1)}")
print("="*80 + "\n")
