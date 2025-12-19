import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from collections import defaultdict
from itertools import combinations
import shutil
from pathlib import Path

from src.config import load_config
from src.compute import (
    compute_spectra, MarchenkoPastur, detect_communities_corr,
    compute_average_correlation, L_detect_communities
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
method = 'ward'
directory = Path('stability_analysis/')
laplacian = True
tau = 1e-3
resolutions_10 = [7,8, 9, 10, 11, 12]
resolutions_9  = [13, 14, 15, 16, 17, 18, 19, 20]
resolutions = resolutions_10 + resolutions_9
stability_threshold = len(resolutions)  # Appear together in at least ... resolutions

# Clear directory so I don't keep old results
if directory.exists():
    shutil.rmtree(directory)

# Recreate it
directory.mkdir(parents=True, exist_ok=True)

names_file = f'{path_template.format(forest="barro")}{names_template.format(forest="barro", census=4)}'
print(f"Loading species names from: {names_file}")

# Load species names
with open(names_file, 'r') as f:
    species_names = [line.strip() for line in f.readlines()][:num_species]

# Storage for community assignments across resolutions
community_assignments = {}  # resolution -> array of community labels
all_communities = {}  # resolution -> dict mapping community_id to species_list

print("\n" + "="*80)
print("MULTI-RESOLUTION COMMUNITY DETECTION")
print("="*80 + "\n")

# Loop through resolutions
for resolution in resolutions:
    print(f"\n{'─'*80}")
    print(f"Processing resolution: {resolution} m")
    print(f"{'─'*80}")
    
    if resolution in resolutions_9:
        n_communities = 9 

    elif resolution in resolutions_10 :
        n_communities = 10
        
    else:
        raise ValueError(f"ERROR! Wrong resolution: {resolution}")

    print(f'Number of expected communities: {n_communities}') 
    # Compute spectra
    (senm_mean, senm_std, forest_spectra, 
     bins, senm_abundance, forest_abundances) = compute_spectra(
        resolution, num_species=num_species, calculate=False
    )
    
    # Compute average correlation
    (forest_corr, forest_corr_std, 
     forest_corrs_all) = compute_average_correlation(forest_abundances)
    
    # Filter correlation matrix
    filtered_forest_corr = MarchenkoPastur(
        forest_corr, num_species, 
        bins[2]*bins[3], remove_largest=False, remove_small=False
    )
    
    # Community detection
    if laplacian:
        print(f"  Using Laplacian diffusion (τ={tau})")
        (forest_reordered, forest_fcluster, forest_idx,
         forest_linkage_matrix, forest_cut_height) = L_detect_communities(
            filtered_forest_corr, n_communities, 
            tau=tau, return_linkage=True
        )
    else:
        print(f"  Using {method} linkage")
        (forest_reordered, forest_fcluster, forest_idx,
         forest_linkage_matrix, forest_cut_height) = detect_communities_corr(
            filtered_forest_corr, n_communities,
            return_linkage=True, method=method
        )
    
    print(f"  Communities detected: {len(set(forest_fcluster))}")
    print(f"  Cut height: {forest_cut_height:.4f}")
    
    # Store results
    community_assignments[resolution] = forest_fcluster
    
    # Group species by community
    communities = defaultdict(list)
    for idx, community_id in enumerate(forest_fcluster):
        communities[community_id].append(species_names[idx])
    
    all_communities[resolution] = dict(communities)
    
    # Save dendrogram
    os.makedirs('stability_analysis', exist_ok=True)
    plot_community_dendrogram(
        forest_linkage_matrix,
        threshold=forest_cut_height,
        title=f'Forest - Resolution {resolution}m',
        filename=f'stability_analysis/dendrogram_{resolution}m.png',
        show=False
    )

print("\n" + "="*80)
print("STABILITY ANALYSIS")
print("="*80 + "\n")

# ========================================
# Find Stable Species Pairs
# ========================================

def find_stable_pairs(community_assignments, species_names, threshold):
    """
    Find species pairs that co-occur in the same community 
    across multiple resolutions.
    
    Returns
    -------
    stable_pairs : dict
        Maps (species_i, species_j) -> list of resolutions where they co-occur
    """
    n_species = len(species_names)
    pair_cooccurrence = defaultdict(list)
    
    for resolution, assignments in community_assignments.items():
        # Check all species pairs
        for i in range(n_species):
            for j in range(i+1, n_species):
                # If in same community, record this resolution
                if assignments[i] == assignments[j]:
                    pair_cooccurrence[(i, j)].append(resolution)
    
    # Filter for stable pairs (meet threshold)
    stable_pairs = {
        pair: resolutions 
        for pair, resolutions in pair_cooccurrence.items()
        if len(resolutions) >= threshold
    }
    
    return stable_pairs

stable_pairs = find_stable_pairs(
    community_assignments, species_names, stability_threshold
)

print(f"Found {len(stable_pairs)} stable species pairs")
print(f"(co-occurring in ≥{stability_threshold}/{len(resolutions)} resolutions)\n")

# ========================================
# Find Stable Subcommunities (Cliques)
# ========================================

def find_stable_subcommunities(stable_pairs, species_names, min_size=3):
    """
    Find groups of species that all co-occur together (cliques).
    
    Parameters
    ----------
    stable_pairs : dict
        Stable species pairs
    species_names : list
        Species names
    min_size : int
        Minimum subcommunity size
        
    Returns
    -------
    subcommunities : list of sets
        Each set contains species indices forming a stable subcommunity
    """
    # Build adjacency list
    from collections import defaultdict
    graph = defaultdict(set)
    
    for (i, j), _ in stable_pairs.items():
        graph[i].add(j)
        graph[j].add(i)
    
    # Find maximal cliques using greedy approach
    def find_cliques_greedy(graph, min_size):
        cliques = []
        visited = set()
        
        for node in graph:
            if node in visited:
                continue
            
            # Start with this node
            clique = {node}
            candidates = graph[node].copy()
            
            # Greedily add nodes that are connected to all in clique
            while candidates:
                # Pick node with most connections to current clique
                best_node = max(candidates, 
                              key=lambda n: len(graph[n] & clique))
                
                # Check if connected to all in clique
                if graph[best_node] >= clique:
                    clique.add(best_node)
                    candidates &= graph[best_node]
                else:
                    candidates.discard(best_node)
            
            # Add if meets minimum size
            if len(clique) >= min_size:
                cliques.append(clique)
                visited.update(clique)
        
        return cliques
    
    subcommunities = find_cliques_greedy(graph, min_size)
    
    # Remove subsets
    filtered = []
    for c1 in subcommunities:
        is_subset = False
        for c2 in subcommunities:
            if c1 < c2:  # Proper subset
                is_subset = True
                break
        if not is_subset:
            filtered.append(c1)
    
    return filtered

subcommunities = find_stable_subcommunities(stable_pairs, species_names, min_size=3)

print(f"\n{'─'*80}")
print(f"STABLE SUBCOMMUNITIES (≥3 species)")
print(f"{'─'*80}\n")

for idx, subcomm in enumerate(sorted(subcommunities, key=len, reverse=True), 1):
    species_in_subcomm = [species_names[i] for i in sorted(subcomm)]
    print(f"\nSubcommunity {idx} ({len(subcomm)} species):")
    print(f"  Indices: {sorted(subcomm)}")
    print(f"  Species:")
    for sp in species_in_subcomm:
        print(f"    • {sp}")

# ========================================
# Visualize Stability Matrix
# ========================================

def plot_stability_matrix(stable_pairs, species_names, resolutions, filename=None):
    """
    Create a heatmap showing how often each pair co-occurs.
    """
    n_species = len(species_names)
    n_resolutions = len(resolutions)
    
    # Create co-occurrence matrix
    cooccurrence_matrix = np.zeros((n_species, n_species))
    
    for (i, j), res_list in stable_pairs.items():
        cooccurrence_matrix[i, j] = len(res_list)
        cooccurrence_matrix[j, i] = len(res_list)
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 10))
    
    im = ax.imshow(cooccurrence_matrix, cmap='YlOrRd', 
                   vmin=0, vmax=n_resolutions)
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax, label='# Resolutions co-occurring')
    
    # Labels
    ax.set_xlabel('Species Index', fontsize=12)
    ax.set_ylabel('Species Index', fontsize=12)
    ax.set_title(
        f'Species Pair Stability Across Resolutions\n'
        f'({n_resolutions} resolutions tested)', 
        fontsize=14, fontweight='bold'
    )
    
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {filename}")
    
    plt.close()

plot_stability_matrix(
    stable_pairs, species_names, resolutions,
    filename='stability_analysis/stability_matrix.png'
)

# ========================================
# Export Stable Subcommunities
# ========================================

output_dir = 'stability_analysis/stable_subcommunities'
os.makedirs(output_dir, exist_ok=True)

# Save all stable pairs
with open(f'{output_dir}/stable_pairs.txt', 'w') as f:
    f.write(f"Stable Species Pairs (≥{stability_threshold}/{len(resolutions)} resolutions)\n")
    f.write("="*80 + "\n\n")
    
    for (i, j), res_list in sorted(stable_pairs.items(), 
                                   key=lambda x: len(x[1]), reverse=True):
        f.write(f"({i:2d}, {j:2d})  {species_names[i]:<30} + {species_names[j]:<30}  "
               f"[{len(res_list)}/{len(resolutions)} resolutions]\n")
        f.write(f"         Resolutions: {res_list}\n\n")

print(f"✓ Saved stable pairs to: {output_dir}/stable_pairs.txt")

# Save each subcommunity
for idx, subcomm in enumerate(sorted(subcommunities, key=len, reverse=True), 1):
    species_in_subcomm = [species_names[i] for i in sorted(subcomm)]
    
    filename = f'{output_dir}/subcommunity_{idx}.txt'
    with open(filename, 'w') as f:
        f.write(f"Stable Subcommunity {idx}\n")
        f.write("="*80 + "\n")
        f.write(f"Size: {len(subcomm)} species\n")
        f.write(f"Stability: Present in ≥{stability_threshold}/{len(resolutions)} resolutions\n\n")
        f.write("Species:\n")
        for sp_idx, sp_name in zip(sorted(subcomm), species_in_subcomm):
            f.write(f"  {sp_idx:3d}. {sp_name}\n")
    
    print(f"✓ Saved subcommunity {idx} to: {filename}")

# ========================================
# Visualize Subcommunities on Spatial Plot
# ========================================

from src.scatters import plot_species_spatial_distribution

for idx, subcomm in enumerate(sorted(subcommunities, key=len, reverse=True)[:5], 1):
    # Plot spatial distribution of subcommunity members
    species_indices = sorted(list(subcomm))
    print(f'Species indices: {species_indices}') 
    plot_species_spatial_distribution(
        species_indices=species_indices,
        census=4,
        background_image='level.png',
        figsize=(12, 8),
        markersize=3,
        panel_label=f'Subcommunity {idx}',
        filename=f'stability_analysis/stable_subcommunities/subcommunity_{idx}_spatial.png',
        show=True
    )
    
    print(f"✓ Plotted spatial distribution for subcommunity {idx}")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)
print(f"\nResults saved to: stability_analysis/")
print(f"  • Dendrograms for each resolution")
print(f"  • Stability matrix heatmap")
print(f"  • Stable pairs list")
print(f"  • Subcommunity files with species lists")
print(f"  • Spatial distribution plots for top subcommunities")
print("="*80 + "\n")
