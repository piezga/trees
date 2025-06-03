import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from variables import *

# Constants
GRID_SIZE = 512
NUM_SPECIES = 100
NUM_REALIZATIONS = 100


def load_senm_data(realization, n_bins_x, n_bins_y):
    file = senm_spatial_file_template.format(nx=nx, ny=ny, nu=nu, kernel=kernel, realization=realization + 1)
    df = pd.DataFrame(np.loadtxt(f"{simulations_path}{file}"), columns=['x', 'y', 'species_id'])

    top_species = df['species_id'].value_counts().head(NUM_SPECIES).index
    df_top100 = df[df['species_id'].isin(top_species)].copy()
    return df_top100, top_species


def compute_abundance_matrix(df, n_bins_x, n_bins_y):
    df['x_bin'] = (df['x'] / (GRID_SIZE / n_bins_x)).astype(int).clip(0, n_bins_x - 1)
    df['y_bin'] = (df['y'] / (GRID_SIZE / n_bins_y)).astype(int).clip(0, n_bins_y - 1)

    abundance = np.zeros((NUM_SPECIES, n_bins_x * n_bins_y))
    for i, (sp_id, group) in enumerate(df.groupby('species_id')):
        bin_counts = group.groupby(['x_bin', 'y_bin']).size()
        for (x, y), count in bin_counts.items():
            abundance[i, x * n_bins_y + y] = count
    return abundance


def compute_mean_spectrum(n_bins_x, n_bins_y):
    eig_matrix = np.zeros((NUM_REALIZATIONS, NUM_SPECIES))
    for realization in range(NUM_REALIZATIONS):
        df_top100, top_species = load_senm_data(realization, n_bins_x, n_bins_y)
        abundance = compute_abundance_matrix(df_top100, n_bins_x, n_bins_y)
        corr = np.nan_to_num(np.corrcoef(abundance), nan=0)
        eig_matrix[realization] = np.sort(np.linalg.eigvalsh(corr))[::-1]
    return np.mean(eig_matrix, axis=0), np.std(eig_matrix, axis=0)


def load_forest_data(forest, census):
    df = pd.read_csv(f"{path_template.format(forest=forest)}{census_template.format(forest=forest, census=census)}")
    names_file = f"{path_template.format(forest=forest)}{names_template.format(forest=forest, census=census)}"
    
    try:
        names = np.loadtxt(names_file, dtype='str', ndmin=1)[:NUM_SPECIES]
        names = [n[0] if isinstance(n, np.ndarray) else str(n) for n in names]
    except ValueError:
        with open(names_file) as f:
            names = [line.strip() for line in f if line.strip()][:NUM_SPECIES]
    return df, names


def compute_forest_spectrum(df, names, n_bins_x, n_bins_y):
    x_bins = pd.cut(df.x, bins=n_bins_x, labels=False)
    y_bins = pd.cut(df.y, bins=n_bins_y, labels=False)
    species_matrix = np.array([
        np.bincount(x_bins[df.name == name] * n_bins_y + y_bins[df.name == name],
                    minlength=n_bins_x * n_bins_y) for name in names
    ])
    corr = np.nan_to_num(np.corrcoef(species_matrix), nan=0)
    return np.sort(np.linalg.eigvalsh(corr))[::-1]


def plot_spectrum(mean_eig, std_eig, forest_eigenvalues, boxlen_label, forest, n_bins_x, n_bins_y, path):
    plt.figure(figsize=(10, 6))
    x = np.arange(1, NUM_SPECIES + 1)
    plt.errorbar(x, mean_eig, yerr=std_eig, fmt='o-', label='SENM Data', markersize=5, capsize=3, alpha=0.7)
    plt.loglog(x, forest_eigenvalues, 'o-', label=f'{forest} Data', markersize=5)
    plt.xscale('log'); plt.yscale('log')
    plt.xlabel('Eigenvalue Rank'); plt.ylabel('Eigenvalue Magnitude')
    plt.legend(); plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.title(f'{n_bins_x}x{n_bins_y} Bins (100 Realizations)')
    plt.text(0.95, 0.95, f'Box length: {boxlen_label}', transform=plt.gca().transAxes,
             ha='right', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
    filename = f"{path}plots/eigenvalue_spectrum_{boxlen_label}.png"
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"Saved plot: {filename}")


def main():
    forest, census = 'barro', 2
    bin_numbers = [10, 20, 30, 50, 100]
    path = path_template.format(forest=forest, census=census)
    os.makedirs(path + 'plots', exist_ok=True)

    all_senm_spectra, all_senm_stds, all_forest_spectra, box_labels = [], [], [], []

    for n_bins_x in bin_numbers:
        n_bins_y = n_bins_x
        boxlen = GRID_SIZE / n_bins_x
        boxlen_label = f"{int(round(boxlen))}m"
        print(f'\nAnalyzing {forest} census {census}, box length: {boxlen:.1f}m')

        mean_eig, std_eig = compute_mean_spectrum(n_bins_x, n_bins_y)
        df, names = load_forest_data(forest, census)
        forest_eigenvalues = compute_forest_spectrum(df, names, n_bins_x, n_bins_y)

        plot_spectrum(mean_eig, std_eig, forest_eigenvalues, boxlen_label, forest, n_bins_x, n_bins_y, path)

        all_senm_spectra.append(mean_eig)
        all_senm_stds.append(std_eig)
        all_forest_spectra.append(forest_eigenvalues)
        box_labels.append(boxlen_label)

    # Summary plots
    for data, stds, title, fname in [
        (all_senm_spectra, all_senm_stds, f'SENM Spectra for Different Box Sizes_k_{kernel}', 'all_senm_spectra'),
        (all_forest_spectra, None, f'{forest} Spectra for Different Box Sizes - census {census}', f'all_{forest}_spectra_{census}')
    ]:
        plt.figure(figsize=(10, 6))
        for i, (spectrum, label) in enumerate(zip(data, box_labels)):
            if stds:
                plt.errorbar(range(1, NUM_SPECIES + 1), spectrum, yerr=stds[i], fmt='o-', label=label,
                             markersize=4, capsize=2, elinewidth=1, alpha=0.7)
            else:
                plt.loglog(range(1, NUM_SPECIES + 1), spectrum, label=label)
        plt.xlabel('Eigenvalue Rank'); plt.ylabel('Eigenvalue Magnitude')
        plt.title(title); plt.legend(title='Box Size')
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.tight_layout()
        plt.savefig(f"{path}plots/{fname}.png", dpi=300)
        plt.close()
        print(f"Saved {fname}.png")


if __name__ == "__main__":
    main()

