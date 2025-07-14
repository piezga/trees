
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns
from scipy.special import rel_entr
from scipy.interpolate import interp1d

from variables import *
from functions import load_forest_data, load_senm_data
from functions import get_top_species

# Constants
GRID_SIZE = 500
NUM_SPECIES = 250
N = NUM_SPECIES
NUM_REALIZATIONS = 10



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

def marchenko_pastur_bounds(n_bins_x, n_bins_y, num_species):
    """Calculate the Marchenko-Pastur bounds."""
    ratio = n_bins_x / num_species
    lambda_min = (1 - np.sqrt(ratio))**2
    lambda_max = (1 + np.sqrt(ratio))**2
    return lambda_min, lambda_max

def plot_combined_spectra(all_spectra, labels, title, filename, senm_std=None, n_bins_x=None, n_bins_y=None):
    """Plot multiple spectra together with elegant styling and Marchenko-Pastur bounds."""
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

    # Plot the Marchenko-Pastur bounds
    if n_bins_x is not None and n_bins_y is not None:
        lambda_min, lambda_max = marchenko_pastur_bounds(n_bins_x, n_bins_y, NUM_SPECIES)
        plt.axhline(lambda_min, color='r', linestyle='--', label=f'Marchenko-Pastur Min: {lambda_min:.2f}')
        plt.axhline(lambda_max, color='r', linestyle='--', label=f'Marchenko-Pastur Max: {lambda_max:.2f}')
        plt.fill_between(x, lambda_min, lambda_max, color='gray', alpha=0.2, label='Marchenko-Pastur Range')
    
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


        
        

def spectra_to_distributions(spectrum1, spectrum2, n_bins=20, threshold=0):
    """
    Converts two spectra to probability distributions using log-binning.
    
    Parameters:
    -----------
    spectrum1 : array
        The first spectrum to convert.
    spectrum2 : array
        The second spectrum to convert.
    n_bins : int
        The number of bins to use for the histogram.
    threshold : float
        Values below this threshold will be ignored in the spectra.
    
    Returns:
    --------
    hist1 : array
        The normalized histogram (distribution) for the first spectrum.
    hist2 : array
        The normalized histogram (distribution) for the second spectrum.
    """
    # Filter out small values (below threshold)
    spectrum1 = spectrum1[spectrum1 > threshold]
    spectrum2 = spectrum2[spectrum2 > threshold]

    # Logarithmic binning
    bins = np.logspace(0, np.log10(max(np.max(spectrum1), np.max(spectrum2))), n_bins)
    
    # Compute the histogram for both spectra
    hist1, _ = np.histogram(spectrum1, bins=bins, density=True)
    hist2, _ = np.histogram(spectrum2, bins=bins, density=True)

    # Normalize the histograms to create probability distributions
    hist1 = hist1 / np.sum(hist1)
    hist2 = hist2 / np.sum(hist2)
    
    return hist1, hist2


def KL_divergence(spectrum1, spectrum2):
    """
    Compute the KL divergence between two spectra.

    Parameters:
    -----------
    spectrum1 : array
        The first spectrum (e.g., SENM).
    spectrum2 : array
        The second spectrum (e.g., BCI).

    Returns:
    --------
    divergence : float
        The KL divergence between the two spectra.
    """
    # Convert to numpy arrays
    p = np.asarray(spectrum1, dtype=np.float64)
    q = np.asarray(spectrum2, dtype=np.float64)

    # Normalize both spectra to be probability distributions
    p /= np.sum(p)
    q /= np.sum(q)

    # Avoid division by zero or log(0) by using a small epsilon
    epsilon = 1e-10
    p = np.clip(p, epsilon, 1)
    q = np.clip(q, epsilon, 1)

    # Compute KL divergence
    divergence = np.sum(p * np.log(p / q))
    
    return divergence

        
        
        
def plot_distributions(p, q, divergence, title_prefix=""):
    """
    Plot two distributions side by side with KL divergence annotated.
    MP maximum and minimum eigenvalues are also added
    """
    x = np.arange(len(p))
    width = 0
    
    
    

    plt.figure(figsize=(10, 6))
    plt.plot(x + width/2, q,'o', label='SENM Distribution', alpha=0.7)
    plt.plot(x - width/2, p, 'o',label='Forest Distribution', alpha=0.7)
    plt.xlabel('Bins')
    plt.ylabel('Probability')
    plt.title(f"{title_prefix}KL Divergence: {divergence:.4f}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    

def plot_KL(divergence_values, bin_numbers, grid_size, prefix = ''):
    boxlens = [grid_size / int(b) for b in bin_numbers]
    boxlens_labels = [f"{b:.2f}" for b in boxlens]

    
    
    plt.figure(figsize=(10, 6))
    plt.plot(boxlens, divergence_values, marker='o')
    plt.xticks(boxlens, boxlens_labels, rotation=45)
    plt.xlabel("Box Size (m)")
    plt.ylabel("KL Divergence")
    plt.title(f"{prefix} KL Divergence vs. Box Size")
    plt.grid(True)
    plt.tight_layout()
    