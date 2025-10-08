# Import packages
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.special import rel_entr
from scipy.interpolate import interp1d

# Import variables
from variables import *

# Hardcoded constants (delete this!!)
GRID_SIZE = 500
NUM_REALIZATIONS = 10



def load_file_with_padding(filename, N, num_columns):
    # Load available data (will error if file doesn't exist)
    data = np.loadtxt(filename, max_rows=N)
    
    # Pad with zeros if needed
    if len(data) < N:
        padding = np.zeros((N - len(data), num_columns))
        data = np.vstack((data, padding))
    
    return data

def get_forest_file(forest: str, census: int) -> str:
    """
    Returns the formatted filename for the given forest and census.
    Raises FileNotFoundError if the file doesn't exist.
    """
    # Map forest names to their templates
    forest_templates = {
        'wanang': wanang_template,
        'barro': barro_template,
    }
    
    # Check if forest is valid
    if forest not in forest_templates:
        raise ValueError(f"Unknown forest: '{forest}'. Expected one of: {list(forest_templates.keys())}")
    
    # Generate the filename
    file = forest_templates[forest].format(census=census)
    
    # Check if file exists
    if not os.path.exists(file):
        raise FileNotFoundError(f"File not found: '{file}'. Please verify the census or template.")
    
    return file

def load_forest_data(forest, census, num_species):
    # Load main CSV data
    df = pd.read_csv(
        f"{path_template.format(forest=forest)}{census_template.format(forest=forest, census=census)}"
    )

    # Load names file (handling UTF-8 BOM)
    names_file = f"{path_template.format(forest=forest)}{names_template.format(forest=forest, census=census)}"

    try:
        # Safely read lines and strip whitespace
        with open(names_file, encoding='utf-8-sig') as f:
            names = [line.strip() for line in f if line.strip()][:num_species]
    except Exception as e:
        raise RuntimeError(f"Error reading names file '{names_file}': {e}")

    # Filter DataFrame to only top species
    df = df[df['name'].isin(names)]

    return df, names

def load_senm_data(nx, ny, nu, kernel, realization):
    """
    Load SENM spatial data for a single realization.
    
    Parameters:
    -----------
    nx : int
        Number of grid points in x-direction
    ny : int
        Number of grid points in y-direction
    nu : float
        Niche overlap parameter
    kernel : str
        Dispersal kernel type
    realization : int
        Realization number (1-based index)
    simulations_path : str
        Path to simulation data files
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with columns ['x', 'y', 'species_id']
    """
    file_name = senm_spatial_file_template.format(nx = nx,ny = ny,nu = nu ,
                                                  kernel = kernel,realization = realization)
    file_path = f"{simulations_path}{file_name}"
    
    try:
        df = pd.DataFrame(np.loadtxt(file_path), 
                         columns=['x', 'y', 'species_id'])
        return df
    except FileNotFoundError:
        raise FileNotFoundError(f"Simulation file not found at: {file_path}")
    except Exception as e:
        raise Exception(f"Error loading simulation data: {str(e)}")


def get_top_species(df, N):
    """
    Returns the dataframe containing only the N most 
    abundant species
    """

    # 1. Count occurrences of each species
    species_counts = df['species_id'].value_counts()

    # 2. Get the 100 most abundant species
    top_N_species = species_counts.head(N).index

    # 3. Filter the original DataFrame to keep only these species
    filtered_df = df[df['species_id'].isin(top_N_species)]
    
    return filtered_df

def shuffle_labels(df):
    """ 
    Shuffles the labels inside of a DF
    """
    df_randomized = df.copy()
    df_randomized['name'] = np.random.permutation(df['name'])  # Shuffle labels
    return df_randomized

def compute_abundance_matrix(num_species,df, n_bins_x, n_bins_y):
    """Compute species abundance matrix"""
    
    df = df.copy()
    
    df['x_bin'] = (df['x'] / (GRID_SIZE / n_bins_x)).astype(int).clip(0, n_bins_x - 1)
    df['y_bin'] = (df['y'] / (GRID_SIZE / n_bins_y)).astype(int).clip(0, n_bins_y - 1)
    
    abundance = np.zeros((num_species, n_bins_x * n_bins_y))
    for i, (_, group) in enumerate(df.groupby('species_id')):
        bin_counts = group.groupby(['x_bin', 'y_bin']).size()
        for (x, y), count in bin_counts.items():
            abundance[i, x * n_bins_y + y] = count
    return abundance

def compute_mean_senm_spectrum(num_species,n_bins_x, n_bins_y):

    """Compute mean eigenvalue spectrum for SENM"""

    eig_matrix = np.zeros((NUM_REALIZATIONS, num_species))
    for realization in range(NUM_REALIZATIONS):
        df = load_senm_data(nx, ny, nu, kernel, realization + 1)
        df_top_N = get_top_species(df,num_species) 
        abundance = compute_abundance_matrix(num_species,df_top_N, n_bins_x, n_bins_y)
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

def marchenko_pastur_bounds(num_species,n_bins_x, n_bins_y):
    """Calculate the Marchenko-Pastur bounds."""
    ratio = num_species/ (n_bins_x*n_bins_y)

    if ratio > 1 :
        print('Out of MP range')
        lambda_min = 0
        lambda_max = 0
    else : 
        lambda_min = (1 - np.sqrt(ratio))**2
        lambda_max = (1 + np.sqrt(ratio))**2
     
    return lambda_min, lambda_max

def plot_combined_spectra(num_species,all_spectra, labels, title, filename, senm_std=None, n_bins_x=None, n_bins_y=None):
    """Plot multiple spectra together with enhanced styling for presentation."""
    plt.figure(figsize=(12, 7))
    x = np.arange(1, num_species + 1)

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

    has_senm = senm_std is not None

    if has_senm:
        plt.errorbar(x, all_spectra[0], yerr=senm_std, fmt='o-', 
                     color='#1f77b4', capsize=4, alpha=0.8,
                     label='SENM (Reference)', linewidth=3, markersize=9)

    n_forest = len(all_spectra) - (1 if has_senm else 0)
    colors = plt.cm.viridis(np.linspace(0.3, 0.9, n_forest))

    for i, (spectrum, label) in enumerate(zip(all_spectra[1:] if has_senm else all_spectra,
                                              labels[1:] if has_senm else labels)):
        plt.loglog(x, spectrum, 'o-', color=colors[i], 
                   label=f'Census {label}', alpha=0.85)

    # Uncomment if Marchenko-Pastur bounds are to be added
    if n_bins_x is not None and n_bins_y is not None:
        lambda_min, lambda_max = marchenko_pastur_bounds(num_species,n_bins_x, n_bins_y)
        plt.axhline(lambda_min, color='r', linestyle='--', label=f'Marchenko-Pastur Min: {lambda_min:.2f}', linewidth=2)
        plt.axhline(lambda_max, color='r', linestyle='--', label=f'Marchenko-Pastur Max: {lambda_max:.2f}', linewidth=2)
        plt.fill_between(x, lambda_min, lambda_max, color='gray', alpha=0.2, label='Marchenko-Pastur Range')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Eigenvalue Rank', fontsize=15, weight='bold')
    plt.ylabel('Eigenvalue Magnitude', fontsize=15, weight='bold')
    plt.title(title, fontsize=17, weight='bold')
    plt.legend(fontsize=12, framealpha=0.95)
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
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

def square_difference(spectrum1, spectrum2):
    
    difference = spectrum1 - spectrum2
    sq_diff = difference**2
    
    return sq_diff



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


def square_diff_above_MP(spectrum_A, spectrum_B, lambda_max):
    """
    Calculate the squared difference of all eigenvalues above lambda_max
    between two spectra.

    Parameters:
    spectrum_A, spectrum_B: arrays of eigenvalues (sorted in descending order)
    lambda_max: threshold value

    Returns:
    squared_difference: sum of squared differences for eigenvalues > lambda_max
    """
    # Get eigenvalues above lambda_max from both spectra
    eigenvalues_A_above = spectrum_A[spectrum_A > lambda_max]
    eigenvalues_B_above = spectrum_B[spectrum_B > lambda_max]

    # Make sure we have the same number of eigenvalues above threshold
    min_length = min(len(eigenvalues_A_above), len(eigenvalues_B_above))

    if min_length == 0:
        return 0  # No eigenvalues above threshold

    # Take the first min_length eigenvalues from both arrays
    eigenvalues_A_above = eigenvalues_A_above[:min_length]
    eigenvalues_B_above = eigenvalues_B_above[:min_length]

    # Calculate relative squared difference
    squared_diff = np.sum((eigenvalues_A_above - eigenvalues_B_above) ** 2)
    magnitude_A = np.sum(eigenvalues_A_above ** 2)
    magnitude_B = np.sum(eigenvalues_B_above ** 2)
    
    # Normalize by average magnitude
    normalized_diff = squared_diff / ((magnitude_A + magnitude_B) / 2)

    return normalized_diff


