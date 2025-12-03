import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

from src.config import load_config
from src.compute import (
        compute_spectra, MarchenkoPastur, detect_communities_corr,
        compute_average_correlation
        )
from src.plotting import (
    plot_community_dendrogram,
    plot_community_data
    )
from src.utils import load_forest_data

# Load config 
config = load_config()
num_species = 200  # Set num_species to 200 as requested
censuses = config['forests']['censuses']
path_template = config['forests']['templates']['path_template']
census_template = config['forests']['templates']['census_template']
names_template = config['forests']['templates']['names_template']

# Parameters
resolution = 10
print(f'Testing for resolution {resolution} m ')
n_comms_forest = 10
method = 'ward'
names_file = f'{path_template.format(forest = "barro")}{names_template.format(forest="barro", census=4)}'
print(names_file)

# Calculating some basic quantities
(senm_mean, senm_std, forest_spectra, 
 bins, senm_abundance, forest_abundances) = compute_spectra(resolution, 
                                                           calculate=True)

# === Calculating correlation matrices and filtering===


(forest_corr, 
 forest_corr_std, forest_corrs_all) = compute_average_correlation(forest_abundances)

filtered_forest_corr = MarchenkoPastur(forest_corr, num_species, 
                                       bins[2]*bins[3], remove_largest=False,
                                       remove_small=False)

# Community detection
(forest_reordered, forest_fcluster, forest_idx,
 forest_linkage_matrix, 
 forest_cut_height) = detect_communities_corr(filtered_forest_corr, n_comms_forest,
                                              return_linkage = True,
                                              method = method)

print(f"\n--- FOREST Community Detection Results ---")
print(f"• Communities detected in the FOREST matrix: {n_comms_forest}")
print(f"• Reordered FOREST matrix shape: {forest_reordered.shape}")
print(f"• Number of clusters in FOREST: {len(set(forest_fcluster))} unique clusters detected.")
print(f"• Community detection linkage matrix (FOREST) shape: {forest_linkage_matrix.shape}")
print(f"• Cutting height for the dendrogram: {forest_cut_height:.3f}")
print(f"-------------------------------------------\n")

plot_community_dendrogram(forest_linkage_matrix,threshold=forest_cut_height,
                          title = 'Forest + ' + method,
                          filename=f'dendrograms/forest_{method}.png')


# Load species names from the names file
with open(names_file, 'r') as f:
    species_names = [line.strip() for line in f.readlines()]

# Ensure the output directory exists
output_dir = 'community_comparison'
os.makedirs(output_dir, exist_ok=True)

# Assuming 'forest_fcluster' contains the community assignments for each species
# We will group species by their community
communities = {}
for idx, community_id in enumerate(forest_fcluster):
    if community_id not in communities:
        communities[community_id] = []
    communities[community_id].append(species_names[idx])

# Write species names for each community to a text file
for community_id, species_list in communities.items():
    # Sort the species names alphabetically
    species_list.sort()

    # Define the output file path
    output_file = os.path.join(output_dir, f'community_{community_id}.txt')

    # Write the sorted species names to the file
    with open(output_file, 'w') as f:
        for species in species_list:
            f.write(f"{species}\n")

    print(f"Community {community_id}: {len(species_list)} species written to {output_file}")


# Loop over each community and generate a plot
for community_id, species_list in communities.items():
    # Get the indices of the species for the current community
    species_indices = [species_names.index(species) for species in species_list]
    
    # Plot and save the community data
    plot_community_data('barro', 4, community_id, species_indices)

