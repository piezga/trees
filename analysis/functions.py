import numpy as np
import os
import pandas as pd
from variables import *


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
    """Load forest census data and top species names"""
    df = pd.read_csv(f"{path_template.format(forest=forest)}{census_template.format(forest=forest, census=census)}")
    names_file = f"{path_template.format(forest=forest)}{names_template.format(forest=forest, census=census)}"
 
    try:
        names = np.loadtxt(names_file, dtype='str', ndmin=1)[:num_species]
        names = [n[0] if isinstance(n, np.ndarray) else str(n) for n in names]
    except ValueError:
        with open(names_file) as f:
            names = [line.strip() for line in f if line.strip()][:num_species]
    
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