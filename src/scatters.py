import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from src.config import load_config
import os
import glob

# Delete previous raster plots
for f in glob.glob("variations/raster_plot_pair*.png"):
    os.remove(f)
print('Removed previous raster plots')

def plot_species_spatial_distribution(
    species_indices,
    census=8,
    forest='barro',
    background_image=None,
    ax=None,
    colors=None,
    markersize=2,
    marker_scale=10,
    figsize=(10, 6),
    xlim=(0, 1000),
    ylim=(0, 500),
    show_legend=True,
    panel_label=None,
    show=True,
    filename=None
):
    """
    Plot spatial distribution of species on the forest plot.
    
    Parameters
    ----------
    species_indices : list or int
        Species index/indices to plot (0-99 for top 100 species)
        If int, converts to single-item list
    census : int
        Census number (1-8)
    forest : str
        Forest name (e.g., 'barro')
    background_image : str, optional
        Path to background image (e.g., elevation map)
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, creates new figure
    colors : list, optional
        List of colors for each species. If None, uses default palette
    markersize : float
        Size of scatter plot markers
    marker_scale : float
        Scale factor for legend markers
    figsize : tuple
        Figure size if ax is None
    xlim : tuple
        X-axis limits (min, max)
    ylim : tuple
        Y-axis limits (min, max)
    show_legend : bool
        Whether to show legend
    panel_label : str, optional
        Label for panel (e.g., 'A', 'B')
    show : bool
        Whether to display the plot
    filename : str, optional
        Path to save figure
        
    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    # Convert single index to list
    if isinstance(species_indices, int):
        species_indices = [species_indices]
    
    # Load config
    config = load_config()
    num_species = config['analysis']['num_species']
    path_template = config['forests']['templates']['path_template']
    census_template = config['forests']['templates']['census_template']
    names_template = config['forests']['templates']['names_template']
    
    # Load forest data
    forest_file = f"{path_template.format(forest=forest)}{census_template.format(forest=forest, census=census)}"
    print(forest_file)
    df = pd.read_csv(forest_file)
    
    # Load species names (using census 4 as reference)
    names_file = f"{path_template.format(forest=forest)}{names_template.format(forest=forest, census=4)}"
    with open(names_file, encoding='utf-8-sig') as f:
        species_names = [line.strip() for line in f if line.strip()][:num_species]
    
    # Create axes if not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    
    # Plot background image if provided
    if background_image is not None:
        try:
            im = plt.imread(background_image)
            ax.imshow(im[::-1], origin='lower', alpha=0.8, extent=[xlim[0], xlim[1], ylim[0], ylim[1]])
        except Exception as e:
            print(f"Warning: Could not load background image: {e}")
    
    # Default colors if not provided
    if colors is None:

        colors = ['blue', 'red', 'green', 'orange', 
                      'purple', 'cyan', 'magenta', 'yellow', 'black', 'brown', 'fuchsia']
        colors = colors[:len(species_indices)]
    
    # Plot each species
    for idx, species_idx in enumerate(species_indices):
        if species_idx >= len(species_names):
            print(f"Warning: Species index {species_idx} out of range (max: {len(species_names)-1})")
            continue
        
        # Get species name
        species_name = species_names[species_idx]
        # Filter data for this species
        species_data = df[df['name'] == species_name.strip()]
        
        if len(species_data) == 0:
            print(f"Warning: No data found for species {species_idx} ({species_name})")
            continue
        
        # Extract coordinates (assuming columns 0=x, 1=y)
        x = species_data['x'].values
        y = species_data['y'].values
        
        # Plot
        ax.scatter(
            x, y, 
            s=markersize, 
            color=colors[idx],
            rasterized=True,
            label=f"{species_idx}: {species_name}",
            alpha=0.7
        )
    
    # Set limits
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    # Labels
    ax.set_xlabel('X (m)', fontsize=12)
    ax.set_ylabel('Y (m)', fontsize=12)
    ax.set_title(f'Spatial Distribution - Census {census}', fontsize=14, fontweight='bold')
    
    # Legend
    if show_legend and len(species_indices) > 0:
        lgd = ax.legend(
            ncol=min(2, len(species_indices)),
            loc='upper center',
            markerscale=marker_scale,
            columnspacing=0.4,
            handlelength=0.5,
            handletextpad=0.4,
            bbox_to_anchor=(0.5, 1.15),
            fontsize=10,
            framealpha=0.9
        )
        # Make legend markers more visible
        for legend_handle in lgd.legend_handles:
            legend_handle.set_sizes([20])
    
    # Panel label
    if panel_label is not None:
        ax.text(-0.05, 1.03, panel_label, transform=ax.transAxes, 
                size=16, weight='bold')
    
    # Save
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {filename}")
    
    # Show
    if show:
        plt.show()
    else:
        plt.close()
    
    return ax


# ========================================
# Simple wrapper for quick plotting
# ========================================

def raster_plot(*species_indices, census=4, background=None, filename=None, show=True):
    """
    Quick plotting function for species spatial distribution.
    
    Parameters
    ----------
    *species_indices : int
        One or more species indices to plot
    census : int
        Census number
    background : str, optional
        Path to background image
    filename : str, optional
        Save filename
        
    Example
    -------
    raster_plot(12, 37)
    raster_plot(12, 37, 56, census=5)
    """
    plot_species_spatial_distribution(
        species_indices=list(species_indices),
        census=census,
        background_image=background if background else None,
        filename=filename,
        show=show
    )
