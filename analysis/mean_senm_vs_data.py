import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from variables import *
from functions import load_forest_data, load_senm_data

# Constants
GRID_SIZE = 512
num_species = 100
NUM_REALIZATIONS = 100


def compute_abundance_matrix(df, n_bins_x, n_bins_y):
    """Compute species abundance matrix"""
    df['x_bin'] = (df['x'] / (GRID_SIZE / n_bins_x)).astype(int).clip(0, n_bins_x - 1)
    df['y_bin'] = (df['y'] / (GRID_SIZE / n_bins_y)).astype(int).clip(0, n_bins_y - 1)
    
    abundance = np.zeros((NUM_SPECIES, n_bins_x * n_bins_y))
    for i, (_, group) in enumerate(df.groupby('species_id')):
        bin_counts = group.groupby(['x_bin', 'y_bin']).size()
        for (x, y), count in bin_counts.items():
            abundance[i, x * n_bins_y + y] = count
    return abundance

def compute_mean_spectrum(n_bins_x, n_bins_y):
    """Compute mean eigenvalue spectrum for SENM"""
    eig_matrix = np.zeros((NUM_REALIZATIONS, NUM_SPECIES))
    for realization in range(NUM_REALIZATIONS):
        df_top100, _ = load_senm_data(realization, n_bins_x, n_bins_y)
        abundance = compute_abundance_matrix(df_top100, n_bins_x, n_bins_y)
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
    plt.close()
    print(f"Saved combined plot: {filename}")

def main():
    forest = 'barro'
    censuses = [1,2, 3, 4, 5, 6, 7, 8]  # Your selected censuses to analyze
    bin_numbers = [10, 20, 30, 50, 100]
    
    for n_bins in bin_numbers:
        print(f"\n{'='*50}\nAnalyzing {n_bins}x{n_bins} bins\n{'='*50}")
        
        # Compute SENM reference spectrum once per bin size
        senm_mean, senm_std = compute_mean_spectrum(n_bins, n_bins)
        boxlen = GRID_SIZE / n_bins
        boxlen_label = f"{int(round(boxlen))}m"
        
        # Prepare storage for all census data
        all_spectra = [senm_mean]
        all_labels = ['SENM']
        
        for census in censuses:
            print(f"\nProcessing {forest} census {census}...")
            path = path_template.format(forest=forest)
            os.makedirs(f"{path}plots", exist_ok=True)
            
            # Compute forest spectrum
            df, names = load_forest_data(forest, census)
            forest_spectrum = compute_forest_spectrum(df, names, n_bins, n_bins)
            all_spectra.append(forest_spectrum)
            all_labels.append(str(census))
            
            # Individual plot for this census
            plot_combined_spectra([senm_mean, forest_spectrum], 
                                ['SENM', f'Census {census}'],
                                f'Eigenvalue Spectrum ({n_bins}x{n_bins} bins)',
                                f"{path}plots/spectrum_{n_bins}x{n_bins}_census{census}.png",
                                senm_std=senm_std)
        
        # Combined plot for all censuses at this bin size
        plot_combined_spectra(all_spectra, all_labels,
                            f'Combined Spectra {n_bins}x{n_bins} Bins',
                            f"{path}combined_spectrum_{n_bins}x{n_bins}.png",
                            senm_std=senm_std)

if __name__ == "__main__":
    main()