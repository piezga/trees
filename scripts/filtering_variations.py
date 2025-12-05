import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

from src.config import load_config
from src.compute import (
        compute_spectra, MarchenkoPastur, detect_communities_corr,
        compute_average_correlation, L_detect_communities,
        compute_filtering_variation
        )
from src.plotting import (
    plot_community_dendrogram,
    plot_community_data,
    plot_corr_with_communities
    )
from src.utils import load_forest_data

# Load config 
config = load_config()
num_species = 100
censuses = config['forests']['censuses']
path_template = config['forests']['templates']['path_template']
census_template = config['forests']['templates']['census_template']
names_template = config['forests']['templates']['names_template']

# Parameters
resolution = 7
print(f'Testing for resolution {resolution} m ')
names_file = f'{path_template.format(forest = "barro")}{names_template.format(forest="barro", census=4)}'

variation = compute_filtering_variation(resolution, num_species)
plt.imshow(variation, cmap='coolwarm')
plt.colorbar()
