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
print(f"LINEAR SHRINKAGE REGULARIZATION FOR PRECISION MATRIX")
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

# We have multiple censuses (8 censuses)
# We'll use some for training and one for testing
# X should be (T, N) where T = spatial sites, N = species

n_censuses = len(forest_abundances)
n_sites = forest_abundances[0].shape[1]

print(f"  Number of censuses available: {n_censuses}")
print(f"  Number of spatial sites: {n_sites}")

# Strategy: Use censuses 1-7 for training, census 8 for testing
# OR: Concatenate all censuses to have more samples

# Option 1: Use last census for testing, rest for training
train_censuses = forest_abundances[:-1]  # Censuses 1-7
test_census = forest_abundances[-1]      # Census 8

# Concatenate training censuses along sites dimension
X_train = np.concatenate(train_censuses, axis=1).T  # Shape: (T_train, N)
X_test = test_census.T  # Shape: (T_test, N)

print(f"\n  Training data shape: {X_train.shape} (samples, species)")
print(f"  Test data shape: {X_test.shape} (samples, species)")

T_train, N = X_train.shape
T_test = X_test.shape[0]

# ========================================
# Standardize Data
# ========================================

def standardize_series(X_train, X_test):
    """
    Standardize features to zero mean and unit variance.
    Fit on training, apply to both train and test.
    """
    mean = np.mean(X_train, axis=0)
    std = np.std(X_train, axis=0)
    std[std == 0] = 1  # Avoid division by zero
    
    X_train_std = (X_train - mean) / std
    X_test_std = (X_test - mean) / std
    
    return X_train_std, X_test_std

print("\nStandardizing data...")
X_train_std, X_test_std = standardize_series(X_train, X_test)

print(f"  Training data: mean={np.mean(X_train_std, axis=0).mean():.4f}, "
      f"std={np.std(X_train_std, axis=0).mean():.4f}")

# ========================================
# Compute Prior Matrix M0
# ========================================

print("\nComputing prior matrix M0...")

# M0 = average correlation matrix across all censuses except the test one
# This is analogous to your: M0 = np.average([np.corrcoef(data[i].T) for i in ...])

M0_list = []
for i, abundance in enumerate(train_censuses):
    # Standardize this census
    abundance_std = (abundance - abundance.mean(axis=1, keepdims=True)) / (abundance.std(axis=1, keepdims=True) + 1e-10)
    corr = np.corrcoef(abundance_std)
    M0_list.append(corr)

M0 = np.mean(M0_list, axis=0)

print(f"  M0 shape: {M0.shape}")
print(f"  M0 diagonal: {np.diag(M0)[:5]}...")  # Should be all 1s
print(f"  M0 off-diagonal range: [{M0[~np.eye(N, dtype=bool)].min():.3f}, "
      f"{M0[~np.eye(N, dtype=bool)].max():.3f}]")

# ========================================
# Linear Shrinkage Estimation
# ========================================

print("\n" + "─"*80)
print("LINEAR SHRINKAGE REGULARIZATION")
print("─"*80 + "\n")

def compute_log_likelihood(X, C):
    """
    Compute log-likelihood of data X under Gaussian with covariance C.
    """
    T = X.shape[0]
    sign, logdet = np.linalg.slogdet(C)
    if sign <= 0:
        return -np.inf
    
    C_inv = np.linalg.inv(C)
    
    # Log-likelihood: -T/2 * log(2π) - T/2 * log|C| - 1/2 * sum_t (x_t^T C^{-1} x_t)
    quadratic_term = 0
    for i in range(T):
        quadratic_term += X[i, :] @ C_inv @ X[i, :]
    
    log_lik = -T/2 * np.log(2*np.pi) - T/2 * logdet - 0.5 * quadratic_term
    
    return log_lik / T  # Normalize by number of samples


def linear_shrinkage_cv(X_train, X_test, M0, alphas=None):
    """
    Linear shrinkage with cross-validation.
    
    C_shrunk(α) = α * M0 + (1-α) * C_empirical
    
    where C_empirical = (1/T) * X^T X (sample covariance)
    
    Parameters
    ----------
    X_train : array, shape (T_train, N)
        Training data
    X_test : array, shape (T_test, N)
        Test data
    M0 : array, shape (N, N)
        Prior/target matrix (e.g., identity or average correlation)
    alphas : array, optional
        Shrinkage parameters to try (0 to 1)
        
    Returns
    -------
    results : dict
        Contains best_alpha, C_clean, test_likelihood
    """
    if alphas is None:
        alphas = np.linspace(0, 1, 21)  # 0.0, 0.05, 0.1, ..., 1.0
    
    # Compute empirical covariance from training data
    C_empirical = np.cov(X_train.T)  # Shape: (N, N)
    
    print(f"  Testing {len(alphas)} alpha values...")
    print(f"  Empirical covariance range: [{C_empirical.min():.3f}, {C_empirical.max():.3f}]")
    
    test_likelihoods = []
    
    for alpha in alphas:
        # Shrinkage: C(α) = α * M0 + (1-α) * C_empirical
        C_shrunk = alpha * M0 + (1 - alpha) * C_empirical
        
        # Evaluate on test set
        test_lik = compute_log_likelihood(X_test, C_shrunk)
        test_likelihoods.append(test_lik)
        
        if alphas[0] == alpha or alphas[-1] == alpha or alpha == alphas[len(alphas)//2]:
            print(f"    α={alpha:.2f}: test log-likelihood={test_lik:.4f}")
    
    # Find best alpha
    test_likelihoods = np.array(test_likelihoods)
    best_idx = np.argmax(test_likelihoods)
    best_alpha = alphas[best_idx]
    best_test_lik = test_likelihoods[best_idx]
    
    # Compute final clean covariance with best alpha
    C_clean = best_alpha * M0 + (1 - best_alpha) * C_empirical
    
    results = {
        'best_alpha': best_alpha,
        'C_clean': C_clean,
        'test_likelihood': best_test_lik,
        'all_alphas': alphas,
        'all_test_likelihoods': test_likelihoods,
        'C_empirical': C_empirical,
        'M0': M0
    }
    
    return results


# Run linear shrinkage
print("Performing linear shrinkage with cross-validation...")
results = linear_shrinkage_cv(X_train_std, X_test_std, M0)

best_alpha = results['best_alpha']
C_clean = results['C_clean']
test_likelihood = results['test_likelihood']

print(f"\n  ✓ Best shrinkage parameter: α = {best_alpha:.3f}")
print(f"  ✓ Test log-likelihood: {test_likelihood:.4f}")

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

output_dir = "precision_matrix_analysis"
os.makedirs(output_dir, exist_ok=True)

# Plot 1: Alpha selection curve
fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(results['all_alphas'], results['all_test_likelihoods'], 
        'o-', linewidth=2, markersize=6)
ax.axvline(best_alpha, color='red', linestyle='--', 
          label=f'Best α = {best_alpha:.3f}')
ax.set_xlabel('Shrinkage parameter α', fontsize=12, fontweight='bold')
ax.set_ylabel('Test log-likelihood', fontsize=12, fontweight='bold')
ax.set_title('Cross-validation: Selecting Optimal Shrinkage Parameter', 
            fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend(fontsize=12)

plt.tight_layout()
plt.savefig(f'{output_dir}/alpha_selection_cv.png', dpi=300, bbox_inches='tight')
print(f"✓ Saved alpha selection curve")
plt.show()

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
axes[0, 1].set_title(f"(B) Denoised Correlation\n(α={best_alpha:.3f})", 
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

fig.suptitle(f'Linear Shrinkage Regularization Results (Resolution: {resolution}m)', 
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
np.save(f'{output_dir}/M0_prior.npy', M0)

# Save metadata
with open(f'{output_dir}/regularization_info.txt', 'w') as f:
    f.write(f"LINEAR SHRINKAGE REGULARIZATION RESULTS\n")
    f.write(f"="*80 + "\n\n")
    f.write(f"Resolution: {resolution}m\n")
    f.write(f"Number of species: {num_species}\n")
    f.write(f"Training samples: {T_train}\n")
    f.write(f"Test samples: {T_test}\n\n")
    f.write(f"Optimal shrinkage parameter: α = {best_alpha:.4f}\n")
    f.write(f"Test log-likelihood: {test_likelihood:.4f}\n\n")
    f.write(f"Formula: C_clean = α * M0 + (1-α) * C_empirical\n")
    f.write(f"  α={best_alpha:.3f}: shrinkage towards prior\n")
    f.write(f"  (1-α)={1-best_alpha:.3f}: weight on empirical covariance\n\n")
    f.write(f"Precision matrix: J = inv(C_clean)\n")
    f.write(f"Partial correlation: P_ij = -J_ij / sqrt(J_ii * J_jj)\n")

print(f"✓ Saved results to {output_dir}/")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)
print(f"\nGenerated files:")
print(f"  • alpha_selection_cv.png - Cross-validation curve")
print(f"  • matrix_comparison.png - All matrices side-by-side")
print(f"  • correlation_distributions.png - Histogram comparison")
print(f"  • C_empirical_full.npy - Raw correlation")
print(f"  • C_clean_denoised.npy - Denoised covariance")
print(f"  • J_precision_matrix.npy - Precision matrix (inverse)")
print(f"  • P_partial_correlation.npy - Partial correlations")
print(f"  • regularization_info.txt - Metadata")
print("="*80 + "\n")
