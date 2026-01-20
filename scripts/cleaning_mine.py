import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cmocean
import os
from sklearn.covariance import ShrunkCovariance, LedoitWolf
from sklearn.model_selection import GridSearchCV
from sklearn.base import BaseEstimator

from src.config import load_config
from src.compute import compute_spectra, compute_average_correlation
import src.estimators as es

# X is the data (number of samples T (i.e. times) times number-of-features (N) )

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
print(f" REGULARIZATION FOR PRECISION MATRIX")
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

print(f"  Number of censuses: {len(forest_abundances)}")
print(f"  Species abundance shape per census: {forest_abundances[0].shape}")
print(f"  (species={forest_abundances[0].shape[0]}, sites={forest_abundances[0].shape[1]})")

# Load species names
names_file = f"{path_template.format(forest=forest)}{names_template.format(forest=forest, census=4)}"
with open(names_file, 'r', encoding='utf-8-sig') as f:
    species_names = [line.strip() for line in f.readlines()][:num_species]

# ========================================
# Prepare Data in Required Format
# ========================================

print("\nPreparing data for regularization...")

n_censuses = len(forest_abundances)
n_sites = forest_abundances[0].shape[1]

print(f"  Number of censuses available: {n_censuses}")
print(f"  Number of spatial sites: {n_sites}")

# Concatenate censuses along sites dimension
X = np.concatenate(forest_abundances, axis=1).T  # Shape: (T_train, N)

print(f"\n Data shape: {X.shape} (samples, species)")

T, N = X.shape

# ========================================
# Denoising
# ========================================

print("\n" + "─"*80)
print("DENOISING")
print("─"*80 + "\n")

nsamples = T 

Xtrain = X[:int(0.8*nsamples)]
Xtest = X[int(0.8*nsamples):]

methods = ['shrinkage', 'lasso']

for method in methods:
    
    if method == 'shrinkage':
        res = es.fit_Shrinkage(Xtrain, Xtest, shrinkage=0.1)
    if method == 'lasso':
        res = es.fit_GraphicalLasso(Xtrain, Xtest, alpha=0.1)


    C = X.T @ X / nsamples # empirical covariance
    C_clean = res["Cclean"] # cleaned covariance

# ========================================
# Compute Precision Matrix (Inverse)
# ========================================

    print("\n" + "─"*80)
    print("COMPUTING PRECISION MATRIX")
    print("─"*80 + "\n")

    print("Inverting denoised covariance matrix...")
    J_clean = np.linalg.inv(C_clean)

    print(f"  ✓ Precision matrix shape: {J_clean.shape}")
    print(f"  ✓ Precision matrix range: [{J_clean.min():.3f}, {J_clean.max():.3f}]")

# ========================================
# Standardize Precision Matrix
# ========================================

    def standardize_matrix(J):
        """
        Standardize precision matrix to partial correlation matrix.
        
        Partial correlation: ρ_ij = -J_ij / sqrt(J_ii * J_jj)
        """
        N = J.shape[0]
        P = np.zeros_like(J)
        
        for i in range(N):
            for j in range(N):
                if i == j:
                    P[i, j] = 1.0
                else:
                    P[i, j] = -J[i, j] / np.sqrt(J[i, i] * J[j, j])
        
        return P

    print("Converting precision matrix to partial correlation matrix...")
    P_clean = standardize_matrix(J_clean)

    print(f"  ✓ Partial correlation matrix shape: {P_clean.shape}")
    print(f"  ✓ Diagonal: {np.diag(P_clean)[:5]}...")  # Should be all 1s
    print(f"  ✓ Off-diagonal range: [{P_clean[~np.eye(N, dtype=bool)].min():.3f}, "
          f"{P_clean[~np.eye(N, dtype=bool)].max():.3f}]")

# ========================================
# Compare with Empirical Correlation
# ========================================

    print("\n" + "─"*80)
    print("COMPARISON WITH EMPIRICAL CORRELATION")
    print("─"*80 + "\n")

# Compute empirical correlation from all data
    forest_abundance_avg = np.mean(np.array(forest_abundances), axis=0)
    C_empirical_full = np.corrcoef(forest_abundance_avg)

    print(f"Empirical correlation statistics:")
    print(f"  Mean off-diagonal: {C_empirical_full[~np.eye(N, dtype=bool)].mean():.4f}")
    print(f"  Std off-diagonal: {C_empirical_full[~np.eye(N, dtype=bool)].std():.4f}")

    print(f"\nDenoised covariance (shrunk) statistics:")
    C_clean_corr = C_clean / np.sqrt(np.outer(np.diag(C_clean), np.diag(C_clean)))
    print(f"  Mean off-diagonal: {C_clean_corr[~np.eye(N, dtype=bool)].mean():.4f}")
    print(f"  Std off-diagonal: {C_clean_corr[~np.eye(N, dtype=bool)].std():.4f}")

    print(f"\nPartial correlation (from precision) statistics:")
    print(f"  Mean off-diagonal: {P_clean[~np.eye(N, dtype=bool)].mean():.4f}")
    print(f"  Std off-diagonal: {P_clean[~np.eye(N, dtype=bool)].std():.4f}")

# ========================================
# Visualization
# ========================================

    print("\n" + "─"*80)
    print("GENERATING VISUALIZATIONS")
    print("─"*80 + "\n")

    output_dir = f"precision_matrix_analysis_{method}"
    os.makedirs(output_dir, exist_ok=True)

# Plot 2: Matrix comparison (2x2 grid)
    fig, axes = plt.subplots(2, 2, figsize=(16, 14), constrained_layout=True)

    cmap = cmocean.cm.balance
    vmin, vmax = -0.5, 0.5
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

# Panel A: Empirical correlation
    im0 = axes[0, 0].imshow(C_empirical_full, cmap=cmap, norm=norm)
    axes[0, 0].set_title("(A) Empirical Correlation", fontsize=14, fontweight='bold')
    axes[0, 0].set_xlabel("Species index")
    axes[0, 0].set_ylabel("Species index")
    plt.colorbar(im0, ax=axes[0, 0], label="Correlation", fraction=0.046)

# Panel B: Denoised covariance (as correlation)
    im1 = axes[0, 1].imshow(C_clean_corr, cmap=cmap, norm=norm)
    axes[0, 1].set_title(f"(B) Denoised Correlation", 
                        fontsize=14, fontweight='bold')
    axes[0, 1].set_xlabel("Species index")
    axes[0, 1].set_ylabel("Species index")
    plt.colorbar(im1, ax=axes[0, 1], label="Correlation", fraction=0.046)

# Panel C: Precision matrix
    vmax_prec = np.percentile(np.abs(J_clean[~np.eye(N, dtype=bool)]), 95)
    im2 = axes[1, 0].imshow(J_clean, cmap='RdBu_r', 
                            vmin=-vmax_prec, vmax=vmax_prec)
    axes[1, 0].set_title("(C) Precision Matrix (Inverse)", 
                        fontsize=14, fontweight='bold')
    axes[1, 0].set_xlabel("Species index")
    axes[1, 0].set_ylabel("Species index")
    plt.colorbar(im2, ax=axes[1, 0], label="Precision", fraction=0.046)

# Panel D: Partial correlation
    im3 = axes[1, 1].imshow(P_clean, cmap=cmap, norm=norm)
    axes[1, 1].set_title("(D) Partial Correlation\n(from Precision Matrix)", 
                        fontsize=14, fontweight='bold')
    axes[1, 1].set_xlabel("Species index")
    axes[1, 1].set_ylabel("Species index")
    plt.colorbar(im3, ax=axes[1, 1], label="Partial Correlation", fraction=0.046)

    fig.suptitle(f'{method} Results (Resolution: {resolution}m)', 
                fontsize=16, fontweight='bold')

    plt.savefig(f'{output_dir}/matrix_comparison.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved matrix comparison")
    plt.show()

# Plot 3: Sparsity comparison (histogram)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Empirical correlation distribution
    axes[0].hist(C_empirical_full[~np.eye(N, dtype=bool)].flatten(), 
                bins=50, alpha=0.7, label='Empirical', edgecolor='black')
    axes[0].axvline(0, color='red', linestyle='--', alpha=0.5)
    axes[0].set_xlabel('Correlation', fontsize=12)
    axes[0].set_ylabel('Frequency', fontsize=12)
    axes[0].set_title('Empirical Correlation Distribution', fontsize=14, fontweight='bold')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

# Partial correlation distribution
    axes[1].hist(P_clean[~np.eye(N, dtype=bool)].flatten(), 
                bins=50, alpha=0.7, label='Partial (from precision)', 
                color='orange', edgecolor='black')
    axes[1].axvline(0, color='red', linestyle='--', alpha=0.5)
    axes[1].set_xlabel('Partial Correlation', fontsize=12)
    axes[1].set_ylabel('Frequency', fontsize=12)
    axes[1].set_title('Partial Correlation Distribution', fontsize=14, fontweight='bold')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/correlation_distributions.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved distribution comparison")
    plt.show()

# ========================================
# Save Results
# ========================================

    print("\nSaving matrices to files...")

    np.save(f'{output_dir}/C_empirical_full.npy', C_empirical_full)
    np.save(f'{output_dir}/C_clean_denoised.npy', C_clean)
    np.save(f'{output_dir}/J_precision_matrix.npy', J_clean)
    np.save(f'{output_dir}/P_partial_correlation.npy', P_clean)

# Save metadata
    with open(f'{output_dir}/regularization_info.txt', 'w') as f:
        f.write(f"{method} REGULARIZATION RESULTS\n")
        f.write(f"="*80 + "\n\n")
        f.write(f"Resolution: {resolution}m\n")
        f.write(f"Number of species: {num_species}\n")
        f.write(f"Precision matrix: J = inv(C_clean)\n")
        f.write(f"Partial correlation: P_ij = -J_ij / sqrt(J_ii * J_jj)\n")

    print(f"✓ Saved results to {output_dir}/")

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"\nGenerated files:")
    print(f"  • matrix_comparison.png - All matrices side-by-side")
    print(f"  • correlation_distributions.png - Histogram comparison")
    print(f"  • C_empirical_full.npy - Raw correlation")
    print(f"  • C_clean_denoised.npy - Denoised covariance")
    print(f"  • J_precision_matrix.npy - Precision matrix (inverse)")
    print(f"  • P_partial_correlation.npy - Partial correlations")
    print(f"  • regularization_info.txt - Metadata")
    print("="*80 + "\n")
