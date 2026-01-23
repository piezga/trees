# Import packages
import numpy as np
import os
import pandas as pd
import networkx as nx

from scipy.interpolate import interp1d
from scipy.linalg import expm
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform

from src.config import load_config
from src.utils import  *

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
species_array = config['analysis']['species_array']
censuses = config['forests']['censuses']



def compute_spectra(resolution, num_species, calculate=True):
    """
    All purpose computing function that returns the basic relevant
    quantities at a given resolution
    """
    if calculate:
        n_bins_x = int(forest_grid_width / resolution)
        n_bins_y = int(forest_grid_height / resolution)
        n_bins_x_senm = int(senm_grid_width / resolution)
        n_bins_y_senm = int(senm_grid_height / resolution)
        
        bins = [n_bins_x_senm, n_bins_y_senm, n_bins_x, n_bins_y]
        
        (senm_mean, 
         senm_std, senm_abundance) = compute_mean_senm_spectrum(num_species, 
                                                                n_bins_x_senm, 
                                                                n_bins_y_senm, 
                                                                standardize=True)
        
        forest_spectra = []
        forest_abundances = []
        
        for census in censuses:
            path = path_template.format(forest=forest)
            os.makedirs(f"{path}plots", exist_ok=True)
            df, names = load_forest_data(forest, census, num_species)

            (spectrum, 
             forest_abundance) = compute_forest_spectrum(df, 
                                                         names, 
                                                         n_bins_x, 
                                                         n_bins_y, 
                                                         standardize=True)
            forest_spectra.append(spectrum)
            forest_abundances.append(forest_abundance)
        
        # Save computed values
        np.save(f'quantities/{forest}_senm_spectrum_{num_species}_{resolution}.npy', senm_mean)
        np.save(f'quantities/{forest}_senm_std_{num_species}_{resolution}.npy', senm_std)
        np.save(f'quantities/{forest}_senm_abundance_{num_species}_{resolution}.npy', senm_abundance)
        np.save(f'quantities/{forest}_forest_spectra_{num_species}_{resolution}.npy', forest_spectra)
        np.save(f'quantities/{forest}_forest_abundances_{num_species}_{resolution}.npy', forest_abundances)
        np.save(f'quantities/{forest}_bins_{num_species}_{resolution}.npy', bins)
        
        return senm_mean, senm_std, forest_spectra, bins, senm_abundance, forest_abundances
    
    else:
        # Load cached values
        senm_mean = np.load(f'quantities/{forest}_senm_spectrum_{num_species}_{resolution}.npy')
        senm_std = np.load(f'quantities/{forest}_senm_std_{num_species}_{resolution}.npy')
        senm_abundance = np.load(f'quantities/{forest}_senm_abundance_{num_species}_{resolution}.npy')
        forest_spectra = np.load(f'quantities/{forest}_forest_spectra_{num_species}_{resolution}.npy')
        forest_abundances = np.load(f'quantities/{forest}_forest_abundances_{num_species}_{resolution}.npy')
        bins = np.load(f'quantities/{forest}_bins_{num_species}_{resolution}.npy')
        
        return senm_mean, senm_std, forest_spectra, bins, senm_abundance, forest_abundances

def compute_average_correlation(forest_abundances):
    """
    Compute average correlation matrix across all census abundances.
    
    Parameters
    ----------
    forest_abundances : list of arrays
        List of abundance matrices, one per census
        Each has shape (num_species, num_bins)
    
    Returns
    -------
    avg_correlation : array
        Average correlation matrix across all censuses
    correlation_std : array
        Standard deviation of correlations across censuses
    all_correlations : list of arrays
        Individual correlation matrices for each census
    """
    all_correlations = []
    
    for abundance in forest_abundances:
        corr = np.corrcoef(abundance)
        all_correlations.append(corr)
    
    # Stack and compute mean and std
    corr_stack = np.stack(all_correlations, axis=0)
    avg_correlation = np.mean(corr_stack, axis=0)
    correlation_std = np.std(corr_stack, axis=0)
    
    return avg_correlation, correlation_std, all_correlations

def L_detect_communities(corr_matrix, n_communities=None, tau=1e-3, Th=1e-4, return_linkage=False):
    """
    Detect communities from filtered correlation matrix using Laplacian diffusion clustering.
    
    Parameters
    ----------
    corr_matrix : array
        Correlation matrix
    n_communities : int, optional
        If provided, cuts dendrogram to get exactly this many communities.
        If None, uses Th threshold with distance criterion.
    tau : float
        Diffusion time parameter
    Th : float
        Distance threshold for cutting (only used if n_communities is None)
    return_linkage : bool
        If True, also return linkage matrix and cut height
        
    Returns
    -------
    reordered : array
        Correlation matrix reordered by community membership
    CM : array
        Community membership labels
    idx : array
        Indices for reordering
    linkage_matrix : array (optional)
        Normalized linkage matrix, returned if return_linkage=True
    cut_height : float (optional)
        Height at which dendrogram was cut, returned if return_linkage=True
    """
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
    linkage_matrix = linkage(dists, "ward")
    
    # Normalize linkage distances by maximum
    linkage_matrix_norm = linkage_matrix.copy()
    linkage_matrix_norm[:, 2] = linkage_matrix[:, 2] / linkage_matrix[-1, 2]
    
    # Cut dendrogram
    if n_communities is not None:
        # Use maxclust criterion to get exactly n_communities
        CM = fcluster(linkage_matrix_norm, t=n_communities, criterion='maxclust')
        
        # Calculate the cut height for n_communities
        n_samples = linkage_matrix_norm.shape[0] + 1
        merge_index = n_samples - n_communities - 1
        cut_height = linkage_matrix_norm[merge_index, 2] if merge_index >= 0 else 0
    else:
        # Use distance threshold
        CM = fcluster(linkage_matrix_norm, t=Th, criterion="distance")
        cut_height = Th
    
    # Reorder matrix by community
    idx = np.argsort(CM)
    reordered = np.array([[corr_matrix[i][j] for j in idx] for i in idx])
    
    if return_linkage:
        return reordered, CM, idx, linkage_matrix_norm, cut_height
    else:
        return reordered, CM, idx

def detect_communities_corr(corr_matrix, n_communities, return_linkage=False,
                            method='ward'):
    """
    Detect communities using Ward hierarchical clustering on correlation-based distances
    """
    distance_matrix = np.sqrt(2 * (1 - np.abs(corr_matrix)))
    distance_matrix = (distance_matrix + distance_matrix.T) / 2
    np.fill_diagonal(distance_matrix, 0)
    
    dists = squareform(distance_matrix)
    linkage_matrix = linkage(dists, method=method)
    
    # Normalize linkage distances by maximum
    linkage_matrix_norm = linkage_matrix.copy()
    linkage_matrix_norm[:, 2] = linkage_matrix[:, 2] / linkage_matrix[-1, 2]
    
    # Cut dendrogram to get n_communities
    CM = fcluster(linkage_matrix_norm, t=n_communities, criterion='maxclust')
    
    # Calculate the cut height for n_communities
    n_samples = linkage_matrix_norm.shape[0] + 1
    merge_index = n_samples - n_communities - 1
    cut_height = linkage_matrix_norm[merge_index, 2] if merge_index >= 0 else 0
    
    # Reorder matrix by community
    idx = np.argsort(CM)
    reordered = np.array([[corr_matrix[i][j] for j in idx] for i in idx])
    
    if return_linkage:
        return reordered, CM, idx, linkage_matrix_norm, cut_height
    else:
        return reordered, CM, idx

def compute_mean_senm_spectrum(num_species,n_bins_x, n_bins_y, standardize):

    """Compute mean eigenvalue spectrum for SENM"""

    eig_matrix = np.zeros((NUM_REALIZATIONS, num_species))
    for realization in range(NUM_REALIZATIONS):
        df = load_senm_data(Nx, Ny, nu, kernel, realization + 1)
        df_top_N = get_top_species(df,num_species) 
        abundance = compute_abundance_matrix(num_species,df_top_N, n_bins_x, n_bins_y, data_type = 'senm', standardize=standardize)
        corr = np.nan_to_num(np.corrcoef(abundance), nan=0)
        eig_matrix[realization] = np.sort(np.linalg.eigvalsh(corr))[::-1]
    return np.mean(eig_matrix, axis=0), np.std(eig_matrix, axis=0), abundance

def compute_forest_spectrum(df, names, n_bins_x, n_bins_y, standardize=False):
    """Compute eigenvalue spectrum for forest data.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing at least 'x', 'y', and 'name' columns.
    names : list-like
        List of species names to include.
    n_bins_x, n_bins_y : int
        Number of bins in the x and y directions.
    standardize : bool, optional
        If True, standardize each species’ abundance vector before correlation.
        
    Returns
    -------
    spectrum : numpy.ndarray
        Sorted eigenvalue spectrum (descending order).
    """
    
    # Compute abundance matrix using your earlier function
    abundance = compute_abundance_matrix(
        num_species=len(names),
        df=df[df['name'].isin(names)],
        n_bins_x=n_bins_x,
        n_bins_y=n_bins_y,
        data_type='forest',
        standardize=standardize
    )
    
    # Compute correlation matrix
    corr = np.nan_to_num(np.corrcoef(abundance), nan=0)
    
    # Compute eigenvalue spectrum, sorted descending
    spectrum = np.sort(np.linalg.eigvalsh(corr))[::-1]
    
    return spectrum, abundance


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


def marchenko_pastur_pdf(lmbda, q, sigma2=1.0):
    """Return the Marchenko–Pastur PDF for given eigenvalue array and aspect ratio q."""
    l_min = sigma2 * (1 - np.sqrt(1.0 / q))**2
    l_max = sigma2 * (1 + np.sqrt(1.0 / q))**2
    pdf = np.zeros_like(lmbda)
    mask = (lmbda >= l_min) & (lmbda <= l_max)
    pdf[mask] = q / (2 * np.pi * sigma2 * lmbda[mask]) * np.sqrt((l_max - lmbda[mask]) * (lmbda[mask] - l_min))
    return pdf, l_min, l_max
        
        
        
def square_difference(spectrum1, spectrum2):
    
    difference = spectrum1 - spectrum2
    sq_diff = difference**2
    
    return sq_diff



def number_of_communities(spectrum, lambda_max ):
    """
    Calculates the number of eigenvalues outside of MP range.

    Parameters:
    spectrum: array of eigenvalues (sorted in descending order)
    lambda_max: threshold value

    Returns:
    (int) num_communities: number of eigenvalues above lambda_max
    """
     
    # Get eigenvalues above lambda_max 
    eigenvalues_above = spectrum[spectrum > lambda_max]
    num_communities = np.count_nonzero(eigenvalues_above) 

    return num_communities


def square_diff_above_MP(spectrum_A, spectrum_B, lambda_max_A, lambda_max_B):
    """
    Calculate the squared difference of all eigenvalues above lambda_max
    between two spectra.

    Parameters:
    spectrum_A, spectrum_B: arrays of eigenvalues (sorted in descending order)
    lambda_max: threshold value

    Returns:
    squared_difference: sum of squared differences for eigenvalues > lambda_max
    """
    debug = True
     
    max_lambda = max(lambda_max_A, lambda_max_B)
    # Get eigenvalues above lambda_max from both spectra
    eigenvalues_A_above = spectrum_A[spectrum_A > max_lambda]
    eigenvalues_B_above = spectrum_B[spectrum_B > max_lambda]
    num_communities_A = np.count_nonzero(eigenvalues_A_above) 
    num_communities_B = np.count_nonzero(eigenvalues_B_above)
    diff_communities = num_communities_A - num_communities_B 
    # Make sure we have the same number of eigenvalues above threshold
    min_length = min(len(eigenvalues_A_above), len(eigenvalues_B_above))
    if min_length == 0:
        print('No eigenvalues above threshold!')
        return 0  # No eigenvalues above threshold

    # Take the first min_length eigenvalues from both arrays and normalize
    eigenvalues_A_above = eigenvalues_A_above[:min_length]
    eigenvalues_B_above = eigenvalues_B_above[:min_length]
    # Calculate relative squared difference
    squared_diff = np.sum((eigenvalues_A_above - eigenvalues_B_above) ** 2)/min_length

    return squared_diff, diff_communities


def MarchenkoPastur(C,N,T,remove_largest=False,remove_small=False):
    """Uses Marchenko-Pastur law to remove noise.
        remove_largest (bool), optional
            If ``False``, all the eigenvectors associated to the
            significant eigenvalues will be used to reconstruct the
            de-noised empirical correlation matrix. If ``True``, the
            eigenvector associated to the largest eigenvalue (normally
            known as the ``market`` mode, [2]) is going to be excluded from
            the recontruction step.  metric_distance (bool), optional: If
            ``False``, a signed graph is obtained.  The weights associated
            to the edges represent the de-noised correlation coefficient
            :math:`\rho_{i,j}` between time series :math:`i` and :math:`j`.
            If ``True``, the correlation is transformed by defining a
            metric distance between each pair of nodes where :math:`d_{i,j}
            = \sqrt{2(1-\rho_{i,j})}` as proposed in [3].  threshold_type
            (str): Which thresholding function to use on the matrix of
            weights. See `netrd.utilities.threshold.py` for
            documentation. Pass additional arguments to the thresholder
            using ``**kwargs``."""

            

            

           
    if N > T:
        raise ValueError("L must be greater or equal than N.")

    Q = T / N
    w, v = np.linalg.eigh(C)  # Spectral decomposition of C
    
    w_min = 1 + 1 / Q - 2 * np.sqrt(1 / Q)
    w_max = 1 + 1 / Q + 2 * np.sqrt(1 / Q)
    
    if remove_small:
        selected = w > w_max
    else:
        selected = (w < w_min) | (w > w_max)
        

    if remove_largest:
        selected[-1] = False

    w_signal = w[selected]
    v_signal = v[:, selected]

    C_new = v_signal.dot(np.diag(w_signal)).dot(v_signal.T)
    return C_new


def get_community_bounds(sorted_CM):
    """Get boundary positions for community boxes."""
    bounds = [0]
    current_comm = sorted_CM[0]
    
    for i, comm in enumerate(sorted_CM):
        if comm != current_comm:
            bounds.append(i)
            current_comm = comm
    
    bounds.append(len(sorted_CM))
    return bounds

def compute_community_overlap(CM_method1, CM_method2):
    """
    Compute overlap/agreement between two community assignments.
    
    Uses Adjusted Rand Index (ARI) to measure similarity.
    ARI = 1 means perfect agreement, 0 means random, negative means worse than random.
    
    Parameters
    ----------
    CM_method1 : array
        Community membership from first method
    CM_method2 : array
        Community membership from second method
        
    Returns
    -------
    ari : float
        Adjusted Rand Index
    confusion : array
        Confusion matrix showing overlap between communities
    """
    from sklearn.metrics import adjusted_rand_score, confusion_matrix
    
    ari = adjusted_rand_score(CM_method1, CM_method2)
    confusion = confusion_matrix(CM_method1, CM_method2)
    
    return ari, confusion

def compute_filtering_variation(resolution, num_species, calculate=True):
    """
    Computes how much each pair of species varies its correlation when
    applying the spectral filter.

    Input:
    --------
    -Resolution
    -Number of species

    Output:
    -Matrix of correlation variations (filtered - original)
    """

    # Get abundances and bins
    _, _, _, bins, _, abundances = compute_spectra(resolution, num_species, calculate)
    T = bins[2]*bins[3]

    # Compute correlation matrix
    unfil_corr_matrix, _, _ = compute_average_correlation(abundances)

    # Filter it 
    fil_corr_matrix = MarchenkoPastur(unfil_corr_matrix, num_species, T)
   
    # Subtract
    difference = fil_corr_matrix - unfil_corr_matrix
    
    # Set diagonal to 0 so it doesn't bother us
    np.fill_diagonal(difference,0)

    return difference, unfil_corr_matrix, fil_corr_matrix

def compute_stripped_correlation_matrix(species_abundance, nutrient_data, debug=False):
    """
    Remove nutrient effects from species abundance, returning cleaned time series.
    
    Parameters
    ----------
    species_abundance : array, shape (n_species, n_sites)
        Species abundance data (rows = species, cols = spatial bins)
    nutrient_data : array, shape (n_nutrients, n_sites)
        Nutrient abundance data (rows = nutrients, cols = spatial bins)
        
    Returns
    -------
    stripped_abundance : array, shape (n_species, n_sites)
        Species abundance with nutrient effects removed (residuals)
    stripped_corr : array, shape (n_species, n_species)
        Correlation matrix of stripped abundances
    """
    from sklearn.linear_model import LinearRegression
    n_species = species_abundance.shape[0]
    
    # Transpose so rows are observations (sites) and columns are variables
    X_species = species_abundance.T  # (n_sites, n_species)
    X_nutrients = nutrient_data.T    # (n_sites, n_nutrients)
    
    if debug:
        print('DEBUG')
        print(f'X_species shape = {X_species.shape}')
        print(f'X_nutrients shape = {X_nutrients.shape}')
    
    # Store residuals after regressing out nutrients
    residuals = np.zeros_like(X_species)
    
    print("\nComputing stripped abundances (removing nutrient effects)...")
    print(f"  Regressing out {X_nutrients.shape[1]} nutrients from {n_species} species...")
    
    # For each species, regress out the effect of all nutrients
    for i in range(n_species):
        if (i + 1) % 20 == 0:
            print(f"    Processed {i+1}/{n_species} species...")
        
        # Fit: species_i = β₀ + β₁*nutrient₁ + ... + βₖ*nutrientₖ + ε
        reg = LinearRegression()
        reg.fit(X_nutrients, X_species[:, i])
        
        # Residuals = what's left after removing nutrient effects
        residuals[:, i] = X_species[:, i] - reg.predict(X_nutrients)
    
    # Transpose back to (n_species, n_sites) format
    stripped_abundance = residuals.T
    
    # Also compute correlation for comparison
    stripped_corr = np.corrcoef(residuals.T)
    
    print(f"  ✓ Stripped abundance shape: {stripped_abundance.shape}")
    print(f"  ✓ Stripped correlation matrix computed")
    
    return stripped_corr, stripped_abundance
