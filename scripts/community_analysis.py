import numpy as np
import argparse
import matplotlib
import matplotlib.pyplot as plt

from src.config import load_config
from src.compute import (
        compute_spectra
        )

matplotlib.use('TkAgg')  # Switch from Qt to Tkinter backend

# === Load config ===
config = load_config()


# Argument parser

def parse_args():
    parser = argparse.ArgumentParser(description='Run community analysis')
    
    # Positional integer argument
    parser.add_argument(
        "chosen_resolution",
        type=int,
        help="Length of box (integer)"
    )
    
    # Verbose flag
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    
    return parser.parse_args()

args = parse_args()

resolution = args.chosen_resolution

# ================================================
# === Correlation Matrix Analysis ===
# ================================================

# Compute spectra and get abundance data
_, _, _, bins, senm_abundance, forest_abundance = compute_spectra(resolution)

# Compute correlation matrices
senm_corr = np.corrcoef(senm_abundance)
forest_corr = np.corrcoef(forest_abundance)

# Apply Marchenko-Pastur filter
filtered_senm_corr = MarchenkoPastur(senm_corr, num_species, bins[0]*bins[1], remove_largest=False)
filtered_forest_corr = MarchenkoPastur(forest_corr, num_species, bins[2]*bins[3], remove_largest=False)

# Detect communities on filtered matrices (LAPLACIAN)
tau = 1e-3
Th = 1e-4

"""

senm_reordered, senm_CM, senm_idx, senm_linkage =detect_communities(filtered_senm_corr,
                                                                    tau, Th,
                                                                    return_linkage = True)
forest_reordered, forest_CM, forest_idx, forest_linkage  = detect_communities(filtered_forest_corr,tau,Th,return_linkage = True)

# Plot SENM correlation matrices
plot_correlation_matrices_comparison(
    senm_corr,
    senm_reordered,
    data_type='SENM',
    filename='correlation_matrices_senm.png',
    show=True
)

# Plot Forest correlation matrices
plot_correlation_matrices_comparison(
    forest_corr,
    forest_reordered,
    data_type='Forest',
    filename='correlation_matrices_forest.png',
    show=True
)
"""

#_, _, _, senm_linkage, senm_th = detect_communities_corr(filtered_senm_corr, 17, return_linkage = True)

#_, _, _, forest_linkage, forest_th = detect_communities_corr(filtered_forest_corr, 10, return_linkage = True)

# Plot dendrograms
"""
plot_community_dendrogram(
    senm_linkage,
    threshold=senm_th,
    title='SENM',
    filename='dendrogram_senm.png',
    show=True
)


plot_community_dendrogram(
    forest_linkage,
    threshold=forest_th,
    title='Forest',
    filename='dendrogram_forest.png',
    show=True
)

"""
# Run both community detection methods on forest data
# Method 1: Laplacian diffusion (threshold-based)
forest_reordered_L, forest_CM_L, forest_idx_L, forest_linkage_L = L_detect_communities(
    filtered_forest_corr, 
    tau=1e-3, 
    Th=1e-4, 
    return_linkage=True
)

# Method 2: Ward clustering (fixed N communities)
n_communities = 10
forest_reordered_W, forest_CM_W, forest_idx_W, forest_linkage_W, forest_th_W = detect_communities_corr(
    filtered_forest_corr, 
    n_communities=n_communities, 
    return_linkage=True
)

# === Visual Comparison ===
# 1. Side-by-side correlation matrices with community boundaries
plot_community_comparison(
    filtered_forest_corr,
    forest_CM_L,
    forest_CM_W,
    method1_name='Laplacian Diffusion',
    method2_name='Ward Clustering',
    filename='forest_community_comparison.png',
    show=True
)

# 2. Confusion matrix showing overlap
plot_community_confusion_matrix(
    forest_CM_L,
    forest_CM_W,
    method1_name='Laplacian',
    method2_name='Ward',
    filename='forest_community_confusion.png',
    show=True
)

# === Text Comparison ===
# 3. Print membership table (optional: pass species names if you have them)
print_community_membership_comparison(
    forest_CM_L,
    forest_CM_W,
    species_names=names if 'names' in locals() else None,  # Use actual species names if available
    method1_name='Laplacian',
    method2_name='Ward'
)

# 4. Get numerical overlap measure
ari, confusion = compute_community_overlap(forest_CM_L, forest_CM_W)
print(f"\nAdjusted Rand Index: {ari:.3f}")
print("(1.0 = perfect agreement, 0.0 = random, <0 = worse than random)")
