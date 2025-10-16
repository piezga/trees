import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns

from src.variables import *
from src.functions import load_forest_data, load_senm_data
from src.functions import get_top_species

# Constants
GRID_SIZE = 500
NUM_SPECIES = 100
N = NUM_SPECIES
NUM_REALIZATIONS = 100

def shuffle_labels(df):
    """ 
    Shuffles the labels inside of a DF
    """
    df_randomized = df.copy()
    df_randomized['name'] = np.random.permutation(df['name'])  # Shuffle labels
    return df_randomized


def compute_abundance_matrix(df, n_bins_x, n_bins_y):
    """Compute species abundance matrix"""
    
    df = df.copy()
    
    df['x_bin'] = (df['x'] / (GRID_SIZE / n_bins_x)).astype(int).clip(0, n_bins_x - 1)
    df['y_bin'] = (df['y'] / (GRID_SIZE / n_bins_y)).astype(int).clip(0, n_bins_y - 1)
    
    abundance = np.zeros((NUM_SPECIES, n_bins_x * n_bins_y))
    for i, (_, group) in enumerate(df.groupby('species_id')):
        bin_counts = group.groupby(['x_bin', 'y_bin']).size()
        for (x, y), count in bin_counts.items():
            abundance[i, x * n_bins_y + y] = count
    return abundance

def compute_mean_senm_spectrum(n_bins_x, n_bins_y):

    """Compute mean eigenvalue spectrum for SENM"""

    eig_matrix = np.zeros((NUM_REALIZATIONS, NUM_SPECIES))
    for realization in range(NUM_REALIZATIONS):
        df = load_senm_data(nx, ny, nu, kernel, realization + 1)
        df_top_N = get_top_species(df,N) 
        abundance = compute_abundance_matrix(df_top_N, n_bins_x, n_bins_y)
        corr = np.nan_to_num(np.corrcoef(abundance), nan=0)
        eig_matrix[realization] = np.sort(np.linalg.eigvalsh(corr))[::-1]
    return np.mean(eig_matrix, axis=0), np.std(eig_matrix, axis=0)

def compute_forest_spectrum(df, names, n_bins_x, n_bins_y):
    """Compute eigenvalue spectrum for forest data"""
    x_bins = pd.cut(df.x, bins=n_bins_x, labels=False)
    y_bins = pd.cut(df.y, bins=n_bins_y, labels=False)
    species_matrix = np.array([
        np.bincount(x_bins[df.name == name] * n_bins_y + y_bins[df.name == name],
                   minlength=n_bins_x * n_bins_y) for name in names
    ])
    corr = np.nan_to_num(np.corrcoef(species_matrix), nan=0)
    return np.sort(np.linalg.eigvalsh(corr))[::-1]

def plot_combined_spectra(all_spectra, labels, title, filename, senm_std=None):
    """Plot multiple spectra together with elegant styling"""
    plt.figure(figsize=(12, 7))
    x = np.arange(1, NUM_SPECIES + 1)
    
    # Determine if we have SENM reference data
    has_senm = senm_std is not None
    
    # Plot SENM reference first if present
    if has_senm:
        plt.errorbar(x, all_spectra[0], yerr=senm_std, fmt='o-', 
                    color='#1f77b4', markersize=5, capsize=3, alpha=0.7,
                    label='SENM (Reference)')
    
    # Calculate number of forest censuses to plot
    n_forest = len(all_spectra) - (1 if has_senm else 0)
    colors = plt.cm.viridis(np.linspace(0.3, 0.9, n_forest))
    
    # Plot forest censuses with distinct colors
    for i, (spectrum, label) in enumerate(zip(all_spectra[1:] if has_senm else all_spectra, 
                                            labels[1:] if has_senm else labels)):
        plt.loglog(x, spectrum, 'o-', color=colors[i], markersize=4, 
                  label=f'Census {label}', alpha=0.8)
    
    plt.xscale('log'); plt.yscale('log')
    plt.xlabel('Eigenvalue Rank', fontsize=12)
    plt.ylabel('Eigenvalue Magnitude', fontsize=12)
    plt.title(title, fontsize=14)
    plt.legend(fontsize=10, framealpha=0.9)
    plt.grid(True, which="both", ls="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Saved combined plot: {filename}")

def plot_multiple_eigenvalue_distributions(all_eigenvalues, eigenvalues_random_forest, labels,
                                          plot_type = 'hist', title="Eigenvalue Distribution", filename=None):
    """
    Plot multiple eigenvalue distributions P(lambda) vs. lambda with consistent binning.
    
    Parameters:
    -----------
    all_eigenvalues : list of np.ndarray
        A list of eigenvalue arrays, one for each census.
    eigenvalues_random_forest : list of np.ndarray
        Eigenvalues from randomized forests.
    labels : list of str
        The labels for each census.
    title : str
        The title of the plot.
    filename : str or None
        If a filename is provided, save the plot. Otherwise, display the plot.
    """
    plt.figure(figsize=(10, 6))

    # Adding a small number for the log plot 
    eigenvalue_nudge = 0
    threshold = 1e-4
    
    senm_spectrum = np.array(all_eigenvalues)[0, :]  + eigenvalue_nudge
    mean_forest_spectrum = np.mean(np.array(all_eigenvalues)[1:, :], axis=0)  + eigenvalue_nudge
    mean_random_forest_spectrum = np.mean(np.array(eigenvalues_random_forest), axis=0) + eigenvalue_nudge

    # Combine all eigenvalues to determine global bin range
    all_values = np.concatenate([
        senm_spectrum,  
        mean_forest_spectrum,  
        mean_random_forest_spectrum  
    ])

    all_values = all_values[all_values > 0]  # Remove non-positive values
    
    # Create consistent bins based on all data

    bins = np.logspace(np.log10(np.min(all_values)), np.log10(np.max(all_values)), 50) # Logarithmic bins
    bin_centers = np.sqrt(bins[:-1] * bins[1:])

    
    alpha = 0.5

    if plot_type == 'hist':
        # Plot histograms with the same bins
        plt.hist(senm_spectrum, bins=bins, 
                label='SENM - 100 realizations', alpha=alpha, density=True)
        plt.hist(mean_forest_spectrum, axis=0, bins=bins,
                label='BCI - census averaged', alpha=alpha, density=True)
        plt.hist(mean_random_forest_spectrum, axis=0, bins=bins,
                label='Shuffled BCI - census averaged', alpha=alpha, density=True)

    elif plot_type == 'distribution':

        # Keep significant values only        
        senm_spectrum = senm_spectrum[senm_spectrum > threshold]
        mean_forest_spectrum = mean_forest_spectrum[mean_forest_spectrum > threshold]
        mean_random_forest_spectrum = mean_random_forest_spectrum[mean_random_forest_spectrum > threshold]

        # Compute normalized histograms (density=True)
        hist_senm, _ = np.histogram(senm_spectrum, bins=bins, density=True)
        hist_mean_forest, _ = np.histogram(mean_forest_spectrum, bins=bins, density=True)
        hist_mean_random_forest, _ = np.histogram(mean_random_forest_spectrum, bins=bins, density=True)

        # Mask zeros to avoid log(0) issues
        mask_senm = hist_senm > 0
        mask_mean_forest = hist_mean_forest > 0
        mask_mean_random_forest = hist_mean_random_forest > 0

        # Plot log-log lines only where data is non-zero
        plt.loglog(bin_centers[mask_senm], hist_senm[mask_senm], '-o', label='SENM')
        plt.loglog(bin_centers[mask_mean_forest], hist_mean_forest[mask_mean_forest], '-o', label='BCI - census averaged')
        plt.loglog(bin_centers[mask_mean_random_forest], hist_mean_random_forest[mask_mean_random_forest], '-o', label='Shuffled BCI - census averaged')


    else :
        print("You have to define the plot type!")

    plt.xlabel('Eigenvalue (λ)', fontsize=14)
    plt.ylabel('Density P(λ)', fontsize=14)
    plt.title(title, fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.xlim(left=1e-1)  
 
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Saved eigenvalue distribution plot: {filename}")
        plt.close()
    else:
        plt.show()


def main():
    forest = 'barro'
    censuses = [1, 2, 3, 4, 5, 6, 7, 8]  # Your selected censuses to analyze
    bin_numbers = [100]
    
    for n_bins_unrounded in bin_numbers:
        n_bins = int(n_bins_unrounded)
        print(f"\n{'='*50}\nAnalyzing {n_bins}x{n_bins} bins\n{'='*50}")
        
        boxlen = GRID_SIZE / n_bins
        boxlen_label = f"{int(round(boxlen))}m"

        # Compute SENM reference spectrum once per bin size
        senm_mean, senm_std = compute_mean_senm_spectrum(n_bins, n_bins)
       
        # Prepare storage for all census data
        all_spectra = [senm_mean]  # SENM spectrum first
        all_labels = ['SENM']      # Label for SENM
        all_random_spectra = [] 
        for census in censuses:
            print(f"\nProcessing {forest} census {census}...")
            path = path_template.format(forest=forest)
            os.makedirs(f"{path}plots", exist_ok=True)
            
            # Compute forest spectrum
            df, names = load_forest_data(forest, census, N)
            randomized_df = shuffle_labels(df)
            forest_spectrum = compute_forest_spectrum(df, names, n_bins, n_bins)
            random_forest_spectrum = compute_forest_spectrum(randomized_df, names, n_bins, n_bins)
            all_spectra.append(forest_spectrum)  # Add forest spectra to the list
            all_random_spectra.append(random_forest_spectrum)
            all_labels.append(str(census))       # Add census label

        plot_combined_spectra(all_spectra, all_labels,
                      title=f'Species Abundance Spectrum | Box length = {boxlen:.1f} m',
                      filename=f'{path}plots/spectrum_boxlen_{boxlen:.1f}.png')


if __name__ == "__main__":
    main()
