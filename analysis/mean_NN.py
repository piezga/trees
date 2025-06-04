import numpy as np
import pandas as pd
from scipy.spatial import KDTree
from variables import *
from functions import load_forest_data

# Parameters for real data
forest = 'barro'
census = 8
species = 100

data  = load_forest_data(forest, census, species) 


def mean_nearest_neighbor_distance(df, species_name):
    """
    Calculate mean nearest neighbor distance for a specific species.
    
    Args:
        df: DataFrame with 'x', 'y', and 'name' columns
        species_name: Target species name (e.g., 'hybapr')
    
    Returns:
        Mean distance to nearest neighbor (or np.nan if <2 points)
    """
    # Convert to DataFrame if not already
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df, columns=['x', 'y', 'name'])
    
    # Filter species points and ensure numeric coordinates
    species_df = df[df['name'] == species_name]
    points = species_df[['x', 'y']].apply(pd.to_numeric).values
    
    if len(points) < 2:
        return np.nan
    
    # Calculate distances
    tree = KDTree(points)
    distances = tree.query(points, k=2)[0][:, 1]  # Get distances to nearest neighbor
    return np.mean(distances)


species = 'hybapr'  
mean_nn_dist = mean_nearest_neighbor_distance(data, species)

if np.isnan(mean_nn_dist):
    print(f"Species '{species}' has fewer than 2 points - cannot compute distance")
else:
    print(f"Mean nearest neighbor distance for {species}: {mean_nn_dist:.2f} units")