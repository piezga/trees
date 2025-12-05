import numpy as np
import matplotlib.pyplot as plt
import os
from src.config import load_config
from src.compute import compute_filtering_variation

# Load config 
config = load_config()
num_species = config['analysis']['num_species']

# Parameters
resolutions = np.arange(5, 33, 1)

# Create output directory
os.makedirs("variations", exist_ok=True)

# Loop through resolutions
for resolution in resolutions:
    print(f'Computing variation for resolution {resolution} m')
    
    # Compute filtering variation
    variation = compute_filtering_variation(resolution, num_species, calculate=False)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot with symmetric colormap centered at zero
    vmax = 0.35
    im = ax.imshow(variation, cmap='coolwarm', vmin=-vmax, vmax=vmax)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, label='Correlation Change\n(Filtered - Original)')
    
    # Labels and title
    ax.set_xlabel('Species Index', fontsize=12)
    ax.set_ylabel('Species Index', fontsize=12)
    ax.set_title(f'Correlation Matrix Change After MP Filtering\nResolution: {resolution} m', 
                 fontsize=14, fontweight='bold')
    
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
