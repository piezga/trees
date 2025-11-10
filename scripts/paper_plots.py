import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from src.config import load_config
from src.functions import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib
import matplotlib.colors as mcolors
import cmocean
import networkx as nx
from scipy.linalg import expm
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
import matplotlib.patches as patches

matplotlib.use('TkAgg')  # Switch from Qt to Tkinter backend


# === Load config ===
config = load_config()

# === Parameters ===
Nx = config['senm']['nx']
Ny = config['senm']['ny']
nu = config['senm']['nu']
kernel = config['senm']['kernel']
NUM_REALIZATIONS = config['senm']['num_realizations']

forest_grid_width = config['grid']['forest']['width']
forest_grid_height = config['grid']['forest']['height']
senm_grid_width = config['grid']['senm']['width']
senm_grid_height = config['grid']['senm']['height']

senm_spatial_file_template = config['senm_templates']['spatial']
simulations_path = config['senm_templates']['path']
path_template = config['forests']['templates']['path_template']
census_template = config['forests']['templates']['census_template']
names_template = config['forests']['templates']['names_template']
forest = 'barro'
num_species = config['analysis']['num_species']
censuses = config['forests']['censuses']

# Some parameters
num = 100
verbose = True
calculate = False
print(f'Calculate set to {calculate}')
filter_resolution = 10



# === Compute spectra ===
def compute_spectra(resolution):
    if verbose:
        print(f'Computing for resolution = {resolution} m')
    
    n_bins_x = int(forest_grid_width / resolution)
    n_bins_y = int(forest_grid_height / resolution)
    n_bins_x_senm = int(senm_grid_width / resolution)
    n_bins_y_senm = int(senm_grid_height / resolution)
    
    bins = [n_bins_x_senm, n_bins_y_senm, n_bins_x, n_bins_y]

    senm_mean, senm_std, senm_abundance = compute_mean_senm_spectrum(num, n_bins_x_senm, n_bins_y_senm, standardize=True)

    forest_spectra = []
    for census in censuses:
        path = path_template.format(forest=forest)
        os.makedirs(f"{path}plots", exist_ok=True)
        df, names = load_forest_data(forest, census, num)
        spectrum, forest_abundance = compute_forest_spectrum(df, names, n_bins_x, n_bins_y, standardize = True)
        forest_spectra.append(spectrum)

    return senm_mean, senm_std, forest_spectra, bins, senm_abundance, forest_abundance



# === Setup style ===
plt.rcParams.update({
    "font.size": 20,
    "axes.labelsize": 24,
    "axes.titlesize": 24,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 18,
    "lines.linewidth": 2,
    "lines.markersize": 6,
    "axes.linewidth": 1.2
})

# === Create subplot mosaic ===
fig, ax_dict = plt.subplot_mosaic("AB", figsize=(14, 6))
x = np.arange(1, num + 1)
colors = sns.color_palette("muted", n_colors=2)

# === Panel A: Spectra Plot ===
ax = ax_dict["A"]

if calculate:
    res_50_senm, res_50_senm_std, res_50_forest_list, _ = compute_spectra(50)
    np.save(f'quantities/barro_senm_std_{num}_50.npy', res_50_senm_std)
else: 
    res_50_senm = np.load(f'quantities/{forest}_senm_spectrum_{num}_50.npy')
    res_50_senm_std = np.load(f'quantities/{forest}_senm_std_{num}_50.npy')
    res_50_forest_list= np.load(f'quantities/{forest}_forest_spectra_{num}_50.npy')

ax.plot(x, res_50_senm, 'o--', color=colors[0], label='SENM (50 m)')
ax.set_title('50 m', fontsize = 22, pad = 4)
ax.fill_between(x, res_50_senm - res_50_senm_std, res_50_senm + res_50_senm_std,
                color=colors[0], alpha=0.25)

for spectrum in res_50_forest_list:
    ax.plot(x, spectrum, 'o-', color=colors[0], alpha=0.5, markerfacecolor='white')

if calculate:
    res_5_senm, res_5_senm_std, res_5_forest_list, _ = compute_spectra(5)
    np.save(f'quantities/barro_senm_std_{num}_5.npy', res_5_senm_std)
else: 
    res_5_senm = np.load(f'quantities/{forest}_senm_spectrum_{num}_5.npy')
    res_5_senm_std = np.load(f'quantities/{forest}_senm_std_{num}_5.npy')
    res_5_forest_list= np.load(f'quantities/{forest}_forest_spectra_{num}_5.npy')


axins = inset_axes(ax, width="45%", height="45%", loc='lower left', borderpad=2)
axins.plot(x, res_5_senm, 'o--', color=colors[1], label='SENM (5 m)')
axins.fill_between(x, res_5_senm - res_5_senm_std, res_5_senm + res_5_senm_std,
                   color=colors[1], alpha=0.25)

for spectrum in res_5_forest_list:
    axins.plot(x, spectrum, 'o-', color=colors[1], alpha=0.5, markerfacecolor='white')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1, 90)
ax.set_ylim(1e-3,100)
ax.set_xlabel(r'$\lambda$ rank')
ax.set_ylabel(r'$\lambda$')
ax.grid(True, which="both", linestyle='--', alpha=0.5)
axins.grid(True, linestyle='--', alpha = 0.5)
axins.set_xscale('log')
axins.set_yscale('log')
axins.set_xlim(1, num)
axins.set_ylim(0.5, max(np.max(res_5_senm + res_5_senm_std),
                        max(np.max(spectrum) for spectrum in res_5_forest_list)))
axins.set_title('5 m', fontsize=18, pad=4)
axins.tick_params(labelsize=8)

ax.text(-0.12, 1.05, '(a)', transform=ax.transAxes, fontsize=24, weight='bold')

# === Panel B: Size effect ===
ax_b = ax_dict["B"]
resolutions = list(range(4, 5, 1))
x = resolutions

colors_b = sns.color_palette("colorblind", n_colors=len(num_species))

for idx, num in enumerate(num_species):
    
    if verbose:
        print(f'Num species: {num} out of {num_species}')

    senm_communities = []
    empirical_communities = []
    square_diffs = []

    for resolution in resolutions:
        
        if calculate:
            senm_spectrum, _, forest_spectra, bins = compute_spectra(resolution)
            np.save(f'quantities/{forest}_senm_spectrum_{num}_{resolution}.npy',
                    senm_spectrum)
            np.save(f'quantities/{forest}_forest_spectra_{num}_{resolution}.npy', 
                    forest_spectra)
            np.save(f'quantities/{forest}_bins_{num}_{resolution}.npy', bins)
            
            np.savetxt(f'quantities/{forest}_senm_spectrum_{num}_{resolution}.txt', 
                       senm_spectrum)
            np.savetxt(f'quantities/{forest}_forest_spectra_{num}_{resolution}.txt', 
                       forest_spectra)
            np.savetxt(f'quantities/{forest}_bins_{num}_{resolution}.txt', bins)
        else:
            senm_spectrum = np.load(f'quantities/{forest}_senm_spectrum_{num}_{resolution}.npy')
            forest_spectra = np.load(f'quantities/{forest}_forest_spectra_{num}_{resolution}.npy')
            bins = np.load(f'quantities/{forest}_bins_{num}_{resolution}.npy')
        if verbose:
            print(f'\nBins array is: {bins}')
            print(f'Resolution is {resolution}')

        mean_forest_spectrum = np.mean(np.array(forest_spectra), axis=0)
        
        if verbose:
            print(f'SENM spectrum shape: {np.shape(senm_spectrum)}')
            print(f'Forest spectrum shape: {np.shape(mean_forest_spectrum)}')

        _, lambda_max_forest = marchenko_pastur_bounds(num, bins[2], bins[3])
        _, lambda_max_senm = marchenko_pastur_bounds(num, bins[0], bins[1])
 
        if verbose:
            print(f'Max forest eigenvalue is {lambda_max_forest}')
            print(f'Max senm eigenvalue is {lambda_max_senm}')

        if calculate: 
            np.savetxt(f'quantities/{forest}_forest_MPmax_{num}_{resolution}.txt', 
                       [lambda_max_forest])
            np.savetxt(f'quantities/{forest}_senm_MPmax_{num}_{resolution}.txt', 
                       [lambda_max_senm])
        square_diff = square_diff_above_MP(senm_spectrum, mean_forest_spectrum, 
                                           lambda_max_senm, lambda_max_forest)
        square_diffs.append(square_diff)

        # Plot and save figure
        plt.figure(figsize=(12, 7))

        # Set global font and line styles
        plt.rcParams.update({
            'font.size': 14,
            'font.weight': 'bold',
            'axes.labelweight': 'bold',
            'axes.titlesize': 16,
            'axes.titleweight': 'bold',
            'lines.linewidth': 2.5,
            'lines.markersize': 8
        })
        
        x_spectra = np.arange(1, num + 1)
        plot_spectra = [mean_forest_spectrum, senm_spectrum]
        labels = ['Forest', 'SENM']
        for i in range(len(plot_spectra)):
            plt.loglog(x_spectra, plot_spectra[i], 'o-', label = labels[i])
        
        plt.title(f'{resolution} m')
        plt.grid(True, alpha = 0.5)
        plt.axhline(lambda_max_forest,label = 'MP Forest')
        plt.axhline(lambda_max_senm, color = 'orange', label = 'MP SENM')
        plt.legend()
        plt.ylim(1e-2,1e2)
        filename = f'spectra/{num}_species_resolution_{resolution}.png'
        plt.savefig(filename)
        if verbose: print(f'Saved spectra to {filename}')
        plt.close()

        # === NEW: Plot eigenvalue density vs Marchenko–Pastur ===
        plt.figure(figsize=(8,6))

        # Histogram of eigenvalues
        plt.hist(mean_forest_spectrum, bins=20, density=True, alpha=0.5,
                 label='Forest eigenvalues', color='C0')
        plt.hist(senm_spectrum, bins=20, density=True, alpha=0.5,
                 label='SENM eigenvalues', color='C1')

        # MP theoretical PDF
        q_forest = (bins[2] * bins[3]) / num
        q_senm   = (bins[0] * bins[1]) / num
        l_vals = np.linspace(0, max(lambda_max_forest, lambda_max_senm)*1.2, 400)

        mp_pdf_forest, lmin_forest, lmax_forest = marchenko_pastur_pdf(l_vals, q_forest)
        mp_pdf_senm, lmin_senm, lmax_senm = marchenko_pastur_pdf(l_vals, q_senm)

        plt.plot(l_vals, mp_pdf_forest, 'C0-', lw=2, label='MP Forest')
        plt.plot(l_vals, mp_pdf_senm, 'C1--', lw=2, label='MP SENM')

        # Labels, styling
        plt.xlabel("Eigenvalue")
        plt.ylabel("Density")
        plt.title(f"Eigenvalue Spectrum Density vs MP ({resolution} m, N={num})")
        plt.legend()
        plt.grid(alpha=0.4)
        plt.tight_layout()

        filename_density = f'spectra/{num}_species_density_resolution_{resolution}.png'
        plt.savefig(filename_density)
        if verbose: print(f'Saved MP density plot to {filename_density}')
        plt.close()

       

    # Plot with assigned color and consistent marker style
    ax_b.plot(
        x,
        square_diffs,
        'o--',
        label=f"N = {num}",
        color=colors_b[idx],
        markerfacecolor='white',
        linewidth = '2'
    )

ax_b.set_xlabel("Resolution (m)")
ax_b.set_ylabel(r"$D^2(>\lambda_{\max})$")
ax_b.grid(True, linestyle='--', alpha=0.5)
ax_b.legend(fontsize=22, loc='best', frameon=False)

# Add panel label
ax_b.text(-0.12, 1.05, '(b)', transform=ax_b.transAxes, fontsize=24, weight='bold')

# === Finalize ===
plt.tight_layout()
plt.show()

# === New plot: correlation matrix with filter ===

_, _, _, bins, senm_abundance, forest_abundance = compute_spectra(filter_resolution)

senm_corr = np.corrcoef(senm_abundance)
forest_corr = np.corrcoef(forest_abundance)
filtered_senm_corr = MarchenkoPastur(senm_corr, num, bins[0]*bins[1], remove_largest=False)
filtered_forest_corr = MarchenkoPastur(forest_corr, num, bins[2]*bins[3], remove_largest=False)



# Colormap and normalization centered at zero
cmap = cmocean.cm.balance
vmin, vmax = -0.5, 0.5
norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

# === Plot SENM correlation matrices ===
fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)

im0 = axes[0].imshow(senm_corr, cmap=cmap, norm=norm)
axes[0].set_title("SENM – Unfiltered Correlation")
axes[0].set_xlabel("Species index")
axes[0].set_ylabel("Species index")

im1 = axes[1].imshow(filtered_senm_corr, cmap=cmap, norm=norm)
axes[1].set_title("SENM – Filtered (Marchenko–Pastur)")
axes[1].set_xlabel("Species index")
axes[1].set_ylabel("Species index")

fig.colorbar(im1, ax=axes, orientation="vertical", fraction=0.03, pad=0.04, label="Correlation")
fig.suptitle("SENM Correlation Matrices", fontsize=14, weight="bold")
plt.show()


# ===  Plot Forest correlation matrices ===
fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)

im0 = axes[0].imshow(forest_corr, cmap=cmap, norm=norm)
axes[0].set_title("Forest – Unfiltered Correlation")
axes[0].set_xlabel("Species index")
axes[0].set_ylabel("Species index")

im1 = axes[1].imshow(filtered_forest_corr, cmap=cmap, norm=norm)
axes[1].set_title("Forest – Filtered (Marchenko–Pastur)")
axes[1].set_xlabel("Species index")
axes[1].set_ylabel("Species index")

fig.colorbar(im1, ax=axes, orientation="vertical", fraction=0.03, pad=0.04, label="Correlation")
fig.suptitle("Forest Correlation Matrices", fontsize=14, weight="bold")
plt.show()

# === Helper function for community detection ===
def detect_communities(corr_matrix, tau=1e-3, Th=1e-4):
    """Detect communities from filtered correlation matrix using Laplacian diffusion clustering."""
    corr_pos = np.copy(corr_matrix)
    corr_pos[corr_pos < 0] = 0  # only positive correlations

    # Build graph and Laplacian
    G = nx.from_numpy_array(np.abs(corr_pos))
    G.remove_edges_from(nx.selfloop_edges(G))
    L = nx.laplacian_matrix(G).todense()

    # Diffusion process
    num = expm(-tau * L)
    rho = num / np.trace(num)

    # Symmetric distance matrix
    Trho = np.copy(1.0 / rho)
    Trho = np.tril(Trho) + np.triu(Trho.T, 1)
    np.fill_diagonal(Trho, 0)

    # Hierarchical clustering
    dists = squareform(Trho)
    linkage_matrix = linkage(dists, "complete")
    linkage_matrix = linkage(dists / linkage_matrix[-1, 2], "complete")
    CM = fcluster(linkage_matrix, t=Th, criterion="distance")

    # Reorder matrix by community
    idx = np.argsort(CM)
    reordered = np.array([[corr_matrix[i][j] for j in idx] for i in idx])

    return reordered, CM, idx


# === New plot: correlation matrix with filter ===
_, _, _, bins, senm_abundance, forest_abundance = compute_spectra(filter_resolution)

senm_corr = np.corrcoef(senm_abundance)
forest_corr = np.corrcoef(forest_abundance)
filtered_senm_corr = MarchenkoPastur(senm_corr, num, bins[0]*bins[1], remove_largest=False)
filtered_forest_corr = MarchenkoPastur(forest_corr, num, bins[2]*bins[3], remove_largest=False)

# === Community detection on filtered matrices ===
tau = 1e-3
Th = 1e-4

senm_reordered, senm_CM, senm_idx = detect_communities(filtered_senm_corr, tau, Th)
forest_reordered, forest_CM, forest_idx = detect_communities(filtered_forest_corr, tau, Th)

# === Plot parameters ===
cmap = cmocean.cm.balance
vmin, vmax = -0.5, 0.5
norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

# === Plot SENM correlation matrices ===
fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)

axes[0].imshow(senm_corr, cmap=cmap, norm=norm)
axes[0].set_title("SENM – Unfiltered")
axes[0].set_xlabel("Species index")
axes[0].set_ylabel("Species index")

im1 = axes[1].imshow(senm_reordered, cmap=cmap, norm=norm)
axes[1].set_title("SENM – Filtered & Reordered by Community")
axes[1].set_xlabel("Species index")
axes[1].set_ylabel("Species index")

fig.colorbar(im1, ax=axes, orientation="vertical", fraction=0.03, 
             pad=0.04, label="Correlation")
fig.suptitle("SENM Correlation Matrices", fontsize=14, weight="bold")
plt.show()


# === Plot Forest correlation matrices ===
fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)

axes[0].imshow(forest_corr, cmap=cmap, norm=norm)
axes[0].set_title("Forest – Unfiltered")
axes[0].set_xlabel("Species index")
axes[0].set_ylabel("Species index")

im1 = axes[1].imshow(forest_reordered, cmap=cmap, norm=norm)
axes[1].set_title("Forest – Filtered & Reordered by Community")
axes[1].set_xlabel("Species index")
axes[1].set_ylabel("Species index")

fig.colorbar(im1, ax=axes, orientation="vertical", fraction=0.03, pad=0.04, label="Correlation")
fig.suptitle("Forest Correlation Matrices", fontsize=14, weight="bold")
plt.show()
