import numpy as np
import os
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