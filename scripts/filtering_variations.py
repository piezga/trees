import numpy as np
import matplotlib.pyplot as plt
import os
from src.config import load_config
from src.compute import compute_filtering_variation
from src.scatters import raster_plot
# ========================================
# Helper Plotting Function
# ========================================

def plot_correlation_tracking(
    resolutions,
    top_data,
    bottom_data,
    random_data,
    top_pairs,
    bottom_pairs,
    random_pairs,
    ylabel,
    title,
    filename,
    show=True
):
    """
    Plot correlation tracking across resolutions for three groups of pairs.
    
    Parameters
    ----------
    resolutions : array
        Resolution values for x-axis
    top_data : dict
        Dictionary mapping pairs to values for top pairs
    bottom_data : dict
        Dictionary mapping pairs to values for bottom pairs
    random_data : dict
        Dictionary mapping pairs to values for random pairs
    top_pairs : list
        List of top pair tuples
    bottom_pairs : list
        List of bottom pair tuples
    random_pairs : list
        List of random pair tuples
    ylabel : str
        Y-axis label
    title : str
        Plot title
    filename : str
        Output filename
    show : bool
        Whether to display the plot
    """
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Plot top pairs (solid lines, circles)
    for pair in top_pairs:
        label = f"Species ({pair[0]}, {pair[1]})"
        ax.plot(
            resolutions, top_data[pair], 
            marker='o', linewidth=2, markersize=4, 
            label=label, linestyle='-', alpha=0.8
        )
    
    # Plot bottom pairs (dashed lines, squares)
    for pair in bottom_pairs:
        label = f"Species ({pair[0]}, {pair[1]})"
        ax.plot(
            resolutions, bottom_data[pair], 
            marker='s', linewidth=2, markersize=4,
            label=label, linestyle='--', alpha=0.8
        )
    
    # Plot random pairs (dotted lines, squares)
    for pair in random_pairs:
        label = f"Species ({pair[0]}, {pair[1]})"
        ax.plot(
            resolutions, random_data[pair], 
            marker='s', linewidth=2, markersize=4,
            label=label, linestyle=':', alpha=0.8
        )
    
    # Add horizontal line at zero
    ax.axhline(
        y=0, color='black', linestyle='-.', 
        linewidth=1, alpha=0.5, label='No change'
    )
    
    # Labels and styling
    ax.set_xlabel('Resolution (m)', fontsize=14, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=14, fontweight='bold')
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=10, ncol=2, framealpha=0.9)
    
    # Tight layout
    plt.tight_layout()
    
    # Save plot
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {filename}")
    
    if show:
        plt.show()
    else:
        plt.close()



# Load config 
config = load_config()
num_species = config['analysis']['num_species']

# Parameters
resolutions = np.arange(5, 33, 1)

# Create output directory
os.makedirs("variations", exist_ok=True)

# Hardcoded pairs from resolution 10 + four random pairs
top_pairs = [
    (13, 56),
    (22, 87),
    (30, 62),
    (23, 77),
    (33, 66)
]

bottom_pairs = [
    (14, 72),
    (41, 70),
    (37, 74),
    (35, 86),
    (28, 37)
]

rng = np.random.default_rng()
random_pairs = [
    (rng.integers(1, 99), rng.integers(1, 99)),
    (rng.integers(1, 99), rng.integers(1, 99)),
    (rng.integers(1, 20), rng.integers(1, 20)),
    (rng.integers(70, 99), rng.integers(70, 99))
]

pairs = top_pairs + bottom_pairs + random_pairs

print(f'Extracted random pairs: {random_pairs}')

# Storage for original correlation values and variations across resolutions
top_variations = {pair: [] for pair in top_pairs}
bottom_variations = {pair: [] for pair in bottom_pairs}
top_correlations = {pair: [] for pair in top_pairs}
bottom_correlations = {pair: [] for pair in bottom_pairs}
top_filtered = {pair: [] for pair in top_pairs}
bottom_filtered = {pair: [] for pair in bottom_pairs}

random_variations = {pair: [] for pair in random_pairs}
random_correlations = {pair: [] for pair in random_pairs}
random_filtered = {pair: [] for pair in random_pairs}

# Loop through resolutions and collect correlation values
for resolution in resolutions:
    print(f'Computing variation for resolution {resolution} m')
    
    # Compute filtering variation
    variation, original, filtered = compute_filtering_variation(
        resolution, num_species, calculate=False
    )
    
    # Extract correlation values for tracked pairs
    for pair in top_pairs:
        i, j = pair
        top_variations[pair].append(variation[i, j])
        top_correlations[pair].append(original[i, j])
        top_filtered[pair].append(filtered[i, j])
    
    for pair in bottom_pairs:
        i, j = pair
        bottom_variations[pair].append(variation[i, j])
        bottom_correlations[pair].append(original[i, j])
        bottom_filtered[pair].append(filtered[i, j])

    for pair in random_pairs:
        i, j = pair
        random_variations[pair].append(variation[i, j])
        random_correlations[pair].append(original[i, j])
        random_filtered[pair].append(filtered[i, j])
    
    # Create figure for this resolution
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot with symmetric colormap centered at zero
    vmax = 0.35
    im = ax.imshow(variation, cmap='coolwarm', vmin=-vmax, vmax=vmax)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, label='Correlation Change\n(Filtered - Original)')
    
    # Labels and title
    ax.set_xlabel('Species Index', fontsize=12)
    ax.set_ylabel('Species Index', fontsize=12)
    ax.set_title(
        f'Correlation Matrix Change After MP Filtering\nResolution: {resolution} m', 
        fontsize=14, fontweight='bold'
    )
    
    # Add grid for readability
    ax.grid(False)
    
    # Tight layout
    plt.tight_layout()
    
    # Save figure
    filename = f'variations/correlation_variation_{resolution}m.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f'  → Saved: {filename}')
    
    # Close to free memory
    plt.close()

print(f"\n✓ All variation plots saved to variations/")

# ========================================
# Generate all three plots using helper
# ========================================

# Plot 1: Variation tracking
plot_correlation_tracking(
    resolutions=resolutions,
    top_data=top_variations,
    bottom_data=bottom_variations,
    random_data=random_variations,
    top_pairs=top_pairs,
    bottom_pairs=bottom_pairs,
    random_pairs=random_pairs,
    ylabel='Correlation Change (Filtered - Original)',
    title='Tracking Correlation Changes Across Resolutions\nTop 5 & Bottom 5 Species Pairs from Resolution 10m',
    filename='variations/variation_tracking_across_resolutions.png',
    show=True
)

# Plot 2: Original correlations
plot_correlation_tracking(
    resolutions=resolutions,
    top_data=top_correlations,
    bottom_data=bottom_correlations,
    random_data=random_correlations,
    top_pairs=top_pairs,
    bottom_pairs=bottom_pairs,
    random_pairs=random_pairs,
    ylabel='Unfiltered correlation',
    title='Original Correlations\nTop 5 & Bottom 5 Species Pairs from Resolution 10m',
    filename='variations/original_correlation_across_resolutions.png',
    show=True
)

# Plot 3: Filtered correlations
plot_correlation_tracking(
    resolutions=resolutions,
    top_data=top_filtered,
    bottom_data=bottom_filtered,
    random_data=random_filtered,
    top_pairs=top_pairs,
    bottom_pairs=bottom_pairs,
    random_pairs=random_pairs,
    ylabel='Filtered correlation',
    title='Filtered Correlations\nTop 5 & Bottom 5 Species Pairs from Resolution 10m',
    filename='variations/filtered_correlation_across_resolutions.png',
    show=True
)


#===================
#=== Raster plots===
#===================
for pair in pairs:
    print(pair)
    raster_plot(*pair)  # unpack tuple to separate arguments


