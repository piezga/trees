import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import interp1d
from src.config import load_config


config = load_config()

# Extract SENM parameters
Nx = config['senm']['nx']
Ny = config['senm']['ny']
nu = config['senm']['nu']         
kernel = config['senm']['kernel']
NUM_REALIZATIONS = config['senm']['num_realizations']

# GRID
forest_grid_width= config['grid']['forest']['width']
forest_grid_height= config['grid']['forest']['height']
senm_grid_width= config['grid']['senm']['width']
senm_grid_height= config['grid']['senm']['height']


# Templates
senm_spatial_file_template = config['senm_templates']['spatial']
simulations_path = config['senm_templates']['path']
path_template = config['forests']['templates']['path_template']
census_template = config['forests']['templates']['census_template']
names_template = config['forests']['templates']['names_template']

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

    # Load names file for census 4 (handling UTF-8 BOM)
    names_file = f"{path_template.format(forest=forest)}{names_template.format(forest=forest, census=4)}"

    try:
        # Safely read lines and strip whitespace
        with open(names_file, encoding='utf-8-sig') as f:
            names = [line.strip() for line in f if line.strip()][:num_species]
    except Exception as e:
        raise RuntimeError(f"Error reading names file '{names_file}': {e}")

    # Filter DataFrame to only top species and sort
    df = df[df['name'].isin(names)]
    df = df.sort_values(by='name')
    #print(df['name'].drop_duplicates().reset_index(drop=True))
    return df, names

def load_senm_data(Nx, Ny, nu, kernel, realization):
    """
    Load SENM spatial data for a single realization.
    
    Parameters:
    -----------
    Nx : int
        Number of grid points in x-direction
    Ny : int
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
    file_name = senm_spatial_file_template.format(nx = Nx,ny = Ny,nu = nu ,
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

def compute_abundance_matrix(num_species, df, n_bins_x, n_bins_y, data_type, standardize=False):
    """Compute species abundance matrix.
    
    Parameters
    ----------
    num_species : int
        Number of species to include.
    df : pandas.DataFrame
        Input data with columns 'x', 'y', and a species identifier column.
    n_bins_x, n_bins_y : int
        Number of spatial bins in the x and y directions.
    data_type : str
        Either 'forest' or 'senm', determines grid dimensions and species column name.
    standardize : bool, optional
        If True, standardize each species’ abundance vector (z-score normalization).
    """
    
    if data_type == 'forest':
        width = forest_grid_width
        height = forest_grid_height
    elif data_type == 'senm':
        width = senm_grid_width
        height = senm_grid_height
    else:
        raise ValueError("data_type must be either 'forest' or 'senm'")

    species_col = 'species_id' if data_type == 'senm' else 'name'

    df = df.copy()
    
    # Assign bins
    df['x_bin'] = (df['x'] / (width / n_bins_x)).astype(int).clip(0, n_bins_x - 1)
    df['y_bin'] = (df['y'] / (height / n_bins_y)).astype(int).clip(0, n_bins_y - 1)
    
    # Initialize abundance matrix
    abundance = np.zeros((num_species, n_bins_x * n_bins_y))
    
    for i, (_, group) in enumerate(df.groupby(species_col)):
        if i >= num_species:
            break
        bin_counts = group.groupby(['x_bin', 'y_bin']).size()
        for (x, y), count in bin_counts.items():
            abundance[i, x * n_bins_y + y] = count
    
    # Standardize each row individually if requested
    if standardize:
        row_means = abundance.mean(axis=1, keepdims=True)
        row_stds = abundance.std(axis=1, keepdims=True)
        # Avoid division by zero
        row_stds[row_stds == 0] = 1
        abundance = (abundance - row_means) / row_stds
    
    return abundance

def load_file_with_padding(filename, N, num_columns):
    try:
        data = np.loadtxt(filename, max_rows=N)

        # Handle edge case: single row results in 1D array
        if data.ndim == 1:
            data = np.expand_dims(data, axis=0)

        current_N = data.shape[0]

        # Pad if needed
        if current_N < N:
            padding = np.zeros((N - current_N, num_columns))
            data = np.vstack((data, padding))

        return data

    except Exception as e:
        raise RuntimeError(f"Failed to load or pad file '{filename}': {e}")

def print_community_membership_comparison(
    CM_method1,
    CM_method2,
    species_names=None,
    method1_name='Method 1',
    method2_name='Method 2'
):
    """
    Print a comparison table showing which species belong to which communities.
    
    Parameters
    ----------
    CM_method1 : array
        Community membership from first method
    CM_method2 : array
        Community membership from second method
    species_names : list, optional
        Names of species (if None, uses indices)
    method1_name : str
        Name of first method
    method2_name : str
        Name of second method
    """
    n_species = len(CM_method1)
    
    if species_names is None:
        species_names = [f"Species_{i+1}" for i in range(n_species)]
    
    print(f"\n{'='*80}")
    print(f"Community Membership Comparison")
    print(f"{'='*80}")
    print(f"{'Species':<20} | {method1_name:<20} | {method2_name:<20}")
    print(f"{'-'*80}")
    
    for i in range(n_species):
        print(f"{species_names[i]:<20} | Community {CM_method1[i]:<14} | Community {CM_method2[i]:<14}")
    
    print(f"\n{'='*80}")
    print(f"Summary:")
    print(f"  {method1_name}: {len(np.unique(CM_method1))} communities")
    print(f"  {method2_name}: {len(np.unique(CM_method2))} communities")
    print(f"{'='*80}\n")

def find_consensus_communities(community_dict, threshold=0.7):
    """
    Find species pairs that are consistently grouped together across methods.
    
    Parameters
    ----------
    community_dict : dict
        Dictionary mapping method names to community assignments
    threshold : float
        Fraction of methods that must agree for a pair to be "consensus"
        
    Returns
    -------
    consensus_matrix : array
        Binary matrix where 1 means species are consistently co-clustered
    """
    n_species = len(list(community_dict.values())[0])
    n_methods = len(community_dict)
    
    # Create co-clustering matrix
    co_cluster = np.zeros((n_species, n_species))
    
    for method, CM in community_dict.items():
        # For each pair of species, check if they're in same community
        for i in range(n_species):
            for j in range(i, n_species):
                if CM[i] == CM[j]:
                    co_cluster[i, j] += 1
                    co_cluster[j, i] += 1
    
    # Normalize by number of methods
    co_cluster /= n_methods
    
    # Threshold to get consensus
    consensus = (co_cluster >= threshold).astype(int)
    
    return consensus, co_cluster
