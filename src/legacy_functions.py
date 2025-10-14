import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.special import rel_entr
from scipy.interpolate import interp1d


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


def senm_bootstrap(num_species, n_bins_x, n_bins_y, num_removals, num_iterations):
    """Compute mean eigenvalue spectrum for SENM with bootstrap resampling by removing random species"""
    
    # Initialize matrix to store eigenvalues for all iterations
    eig_matrix = np.zeros((num_iterations, num_species - num_removals))
    
    for iteration in range(num_iterations):
        # Load the base SENM data
        df = load_senm_data(nx, ny, nu, kernel, iteration + 1)
        
        # Get top N species first
        df_top_N = get_top_species(df, num_species)
        
        # Get the unique species IDs from the top N dataframe
        all_species_ids = df_top_N['species_id'].unique()
        
        # Randomly remove 'num_removals' species
        if num_removals > 0:
            # Randomly select species IDs to remove
            species_to_remove = np.random.choice(all_species_ids, size=num_removals, replace=False)
            
            # Remove those species from the dataframe
            df_reduced = df_top_N[~df_top_N['species_id'].isin(species_to_remove)]
            remaining_species = num_species - num_removals
        else:
            df_reduced = df_top_N
            remaining_species = num_species
        
        # Compute abundance matrix with reduced species
        abundance = compute_abundance_matrix(remaining_species, df_reduced, n_bins_x, n_bins_y)
        
        # Compute correlation matrix and eigenvalues
        corr = np.nan_to_num(np.corrcoef(abundance), nan=0)
        eigenvalues = np.sort(np.linalg.eigvalsh(corr))[::-1]
        
        # Store the eigenvalues
        eig_matrix[iteration] = eigenvalues
    
    return np.mean(eig_matrix, axis=0), np.std(eig_matrix, axis=0)

def forest_bootstrap(forest, censuses, num_species, n_bins_x, n_bins_y, num_removals, num_iterations):
    """
    Compute mean eigenvalue spectrum for empirical forest data with bootstrap resampling
    by removing random species and averaging over censuses
    """

    # Initialize matrix to store mean spectra for all iterations
    eig_matrix = np.zeros((num_iterations, num_species - num_removals))

    for iteration in range(num_iterations):
        print(f"\nBootstrap iteration {iteration + 1}/{num_iterations}")

        # List to store spectra for this bootstrap iteration across all censuses
        bootstrap_spectra = []

        for census in censuses:
            print(f"  Processing {forest} census {census}...")

            # Load forest data for this census
            df, names = load_forest_data(forest, census, num_species)

            # Randomly remove 'num_removals' species from both df and names
            if num_removals > 0:
                # Convert names to numpy array for easier indexing
                names_array = np.array(names)

                # Randomly select indices to remove
                indices_to_remove = np.random.choice(len(names_array), size=num_removals,
                                                     replace=False)

                # Get the species names to remove
                species_to_remove = names_array[indices_to_remove]

                # Remove those species from the dataframe
                df_reduced = df[~df['name'].isin(species_to_remove)]

                # Remove those species from the names list
                names_reduced = [name for i,
                                 name in enumerate(names) if i not in indices_to_remove]

                remaining_species = num_species - num_removals
            else:
                df_reduced = df
                names_reduced = names
                remaining_species = num_species

            # Compute forest spectrum with reduced species
            forest_spectrum = compute_forest_spectrum(df_reduced, names_reduced, n_bins_x, 
                                                      n_bins_y)

            bootstrap_spectra.append(forest_spectrum)

        # Average the spectra across all censuses for this bootstrap iteration
        mean_spectrum = np.mean(np.array(bootstrap_spectra), axis=0)

        # Store the mean eigenvalues for this iteration
        eig_matrix[iteration] = mean_spectrum

    return np.mean(eig_matrix, axis=0), np.std(eig_matrix, axis=0)


def load_file_with_padding(filename, N, num_columns):
    # Load available data (will error if file doesn't exist)
    data = np.loadtxt(filename, max_rows=N)

    # Pad with zeros if needed
    if len(data) < N:
        padding = np.zeros((N - len(data), num_columns))
        data = np.vstack((data, padding))

    return data



