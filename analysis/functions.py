import numpy as np


def load_file_with_padding(filename, N, num_columns):
    # Load available data (will error if file doesn't exist)
    data = np.loadtxt(filename, max_rows=N)
    
    # Pad with zeros if needed
    if len(data) < N:
        padding = np.zeros((N - len(data), num_columns))
        data = np.vstack((data, padding))
    
    return data
