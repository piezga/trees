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
from src.functions import compute_forest_spectrum, compute_mean_senm_spectrum, load_forest_data

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




# === Compute spectra ===
def compute_spectra(resolution):
    """
    All purpose computing function that returns the basic relevant
    quantities at a given resolution
    """

    
    n_bins_x = int(forest_grid_width / resolution)
    n_bins_y = int(forest_grid_height / resolution)
    n_bins_x_senm = int(senm_grid_width / resolution)
    n_bins_y_senm = int(senm_grid_height / resolution)
    
    bins = [n_bins_x_senm, n_bins_y_senm, n_bins_x, n_bins_y]

    senm_mean, senm_std, senm_abundance = compute_mean_senm_spectrum(num_species, n_bins_x_senm, n_bins_y_senm, standardize=True)

    forest_spectra = []
    for census in censuses:
        path = path_template.format(forest=forest)
        os.makedirs(f"{path}plots", exist_ok=True)
        df, names = load_forest_data(forest, census, num_species)
        spectrum, forest_abundance = compute_forest_spectrum(df, names, n_bins_x, n_bins_y, standardize = True)
        forest_spectra.append(spectrum)

    return senm_mean, senm_std, forest_spectra, bins, senm_abundance, forest_abundance

# === Helper function for community detection ===
def detect_communities(corr_matrix, tau=1e-3, Th=1e-4, return_linkage=False):
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
    linkage_matrix = linkage(dists, "ward")
    linkage_matrix = linkage(dists / linkage_matrix[-1, 2], "ward")
    CM = fcluster(linkage_matrix, t=Th, criterion="distance")

    # Reorder matrix by community
    idx = np.argsort(CM)
    reordered = np.array([[corr_matrix[i][j] for j in idx] for i in idx])
    
    if return_linkage:
        return reordered, CM, idx, linkage_matrix
    else:
        return reordered, CM, idx


def detect_communities_corr(corr_matrix, n_communities, return_linkage=False):
    """
    Detect communities using Ward hierarchical clustering on correlation-based distances.
    
    Parameters
    ----------
    corr_matrix : array
        Correlation matrix (NxN)
    n_communities : int
        Number of communities to detect
    return_linkage : bool
        If True, also return the linkage matrix for plotting
        
    Returns
    -------
    reordered : array
        Correlation matrix reordered by community membership
    CM : array
        Community membership labels for each node (1 to n_communities)
    idx : array
        Indices for reordering (sorted by community)
    linkage_matrix : array (optional)
        Normalized linkage matrix, returned if return_linkage=True
    """
    # Convert correlations to distances
    # Using d = sqrt(2(1-r)) which is a proper metric for correlations
    # Alternative: d = 1 - |r| for a simpler transformation
    distance_matrix = np.sqrt(2 * (1 - np.abs(corr_matrix)))
    
    # Ensure symmetry and zero diagonal
    distance_matrix = (distance_matrix + distance_matrix.T) / 2
    np.fill_diagonal(distance_matrix, 0)
    
    # Convert to condensed distance matrix for linkage
    dists = squareform(distance_matrix)
    
    # Perform Ward linkage
    linkage_matrix = linkage(dists, method="ward")
    
    # Normalize linkage distances by maximum
    linkage_matrix_norm = linkage_matrix.copy()
    linkage_matrix_norm[:, 2] = linkage_matrix[:, 2] / linkage_matrix[-1, 2]
    
    # Cut dendrogram to get n_communities
    CM = fcluster(linkage_matrix_norm, t=n_communities, criterion='maxclust')
    
    # Reorder matrix by community
    idx = np.argsort(CM)
    reordered = np.array([[corr_matrix[i][j] for j in idx] for i in idx])
    
    if return_linkage:
        return reordered, CM, idx, linkage_matrix_norm
    else:
        return reordered, CM, idx
