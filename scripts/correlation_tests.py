import numpy as np
import argparse
import matplotlib
import matplotlib.pyplot as plt
import os

from src.config import load_config
from src.compute import (
        compute_spectra, MarchenkoPastur, L_detect_communities,
        marchenko_pastur_pdf
        )
from src. plotting import (
    plot_correlation_matrices_comparison)
#matplotlib.use('TkAgg')  # Switch from Qt to Tkinter backend

# === Load config ===
config = load_config()
num_species = config['analysis']['num_species']

# === Pick resolution ===
rng = np.random.default_rng()
resolution = int(rng.integers(5,33))
resolution = 6
print(f'Testing for resolution {resolution} m ')

(senm_mean, senm_std, forest_spectra, 
 bins, senm_abundance, forest_abundance) = compute_spectra(resolution, 
                                                           calculate=True)

# Check standardization

print('\nSENM STANDARDIZATION')
print('\nSENM mean')
print(np.mean(senm_abundance, axis = 1))
print('\nSENM std')
print(np.std(senm_abundance, axis = 1))

print('\nFOREST STANDARDIZATION')
print('\nForest mean')
print(np.mean(forest_abundance, axis = 1))
print('\nForest std')
print(np.std(forest_abundance, axis = 1))

# Is the shape right?
print('\n\n')
print('Testing for correct shape (column)')
print(f'Forest shape is : {np.shape(forest_abundance)}')
print(f'Expecting column: {500*1000/(resolution**2)}')
print(f'SENM shape is   : {np.shape(senm_abundance)}')
print(f'Expecting column: {500*500/(resolution**2)}')

# === Calculating correlation matrices ===

senm_corr = np.corrcoef(senm_abundance)
forest_corr = np.corrcoef(forest_abundance)

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
plt.show()


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
    print(f"Number of outliers above MP bounds: {len(above)}")

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

print(f"\nAll spectral plots saved to: {plot_dir}")

"""
# === Laplacian checks === 


tau = 1e-3
Th = 1e-4


(senm_reordered, senm_CM, 
 senm_idx, senm_linkage) = L_detect_communities(filtered_senm_corr,
                                                tau, Th,
                                                return_linkage = True)
(forest_reordered, forest_CM, 
 forest_idx, forest_linkage)  = L_detect_communities(filtered_forest_corr,
                                                     tau,Th,
                                                     return_linkage = True)

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
