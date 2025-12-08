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

# Hardcoded pairs from resolution 10
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

# Storage for original correlation values and variations across resolutions
top_variations = {pair: [] for pair in top_pairs}
bottom_variations = {pair: [] for pair in bottom_pairs}
top_correlations = {pair: [] for pair in top_pairs}
bottom_correlations = {pair: [] for pair in bottom_pairs}
top_filtered = {pair: [] for pair in top_pairs}
bottom_filtered = {pair: [] for pair in bottom_pairs}

# Loop through resolutions and collect correlation values
for resolution in resolutions:
    print(f'Computing variation for resolution {resolution} m')
    
    # Compute filtering variation
    variation, original, filtered = compute_filtering_variation(resolution, 
                                                      num_species, calculate=False)
    
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

# ========================================
# Plot correlation tracking across resolutions
# ========================================

fig, ax = plt.subplots(figsize=(14, 8))

# Plot top 5 pairs (positive correlation changes)
for pair in top_pairs:
    label = f"Species ({pair[0]}, {pair[1]})"
    ax.plot(resolutions, top_variations[pair], 
            marker='o', linewidth=2, markersize=4, 
            label=label, linestyle='-', alpha=0.8)

# Plot bottom 5 pairs (negative correlation changes)
for pair in bottom_pairs:
    label = f"Species ({pair[0]}, {pair[1]})"
    ax.plot(resolutions, bottom_variations[pair], 
            marker='s', linewidth=2, markersize=4,
            label=label, linestyle='--', alpha=0.8)

# Add horizontal line at zero
ax.axhline(y=0, color='black', linestyle=':', linewidth=1, alpha=0.5, label='No change')

# Labels and styling
ax.set_xlabel('Resolution (m)', fontsize=14, fontweight='bold')
ax.set_ylabel('Correlation Change (Filtered - Original)', fontsize=14, fontweight='bold')
ax.set_title('Tracking Correlation Changes Across Resolutions\nTop 5 & Bottom 5 Species Pairs from Resolution 10m', 
             fontsize=16, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend(loc='best', fontsize=10, ncol=2, framealpha=0.9)

# Tight layout
plt.tight_layout()

# Save tracking plot
tracking_filename = 'variations/variation_tracking_across_resolutions.png'
plt.savefig(tracking_filename, dpi=300, bbox_inches='tight')
print(f"\n✓ Variation tracking plot saved: {tracking_filename}")

plt.show()

#============================
#=== Plot original values ===
#============================

fig, ax = plt.subplots(figsize=(14, 8))

for pair in top_pairs:
    label = f"Species ({pair[0]}, {pair[1]})"
    ax.plot(resolutions, top_correlations[pair], 
            marker='o', linewidth=2, markersize=4, 
            label=label, linestyle='-', alpha=0.8)

# Plot bottom 5 pairs (negative correlation changes)
for pair in bottom_pairs:
    label = f"Species ({pair[0]}, {pair[1]})"
    ax.plot(resolutions, bottom_correlations[pair], 
            marker='s', linewidth=2, markersize=4,
            label=label, linestyle='--', alpha=0.8)

# Add horizontal line at zero
ax.axhline(y=0, color='black', linestyle=':', linewidth=1, alpha=0.5, label='No change')

# Labels and styling
ax.set_xlabel('Resolution (m)', fontsize=14, fontweight='bold')
ax.set_ylabel('Unfiltered correlation', fontsize=14, fontweight='bold')
ax.set_title('Original Correlations \nTop 5 & Bottom 5 Species Pairs from Resolution 10m', 
             fontsize=16, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend(loc='best', fontsize=10, ncol=2, framealpha=0.9)

# Tight layout
plt.tight_layout()

# Save tracking plot
tracking_filename = 'variations/original_correlation_across_resolutions.png'
plt.savefig(tracking_filename, dpi=300, bbox_inches='tight')
print(f"\n✓ Correlation tracking plot saved: {tracking_filename}")

plt.show()

#============================
#=== Plot filtered values ===
#============================

fig, ax = plt.subplots(figsize=(14, 8))

for pair in top_pairs:
    label = f"Species ({pair[0]}, {pair[1]})"
    ax.plot(resolutions, top_filtered[pair], 
            marker='o', linewidth=2, markersize=4, 
            label=label, linestyle='-', alpha=0.8)

# Plot bottom 5 pairs (negative correlation changes)
for pair in bottom_pairs:
    label = f"Species ({pair[0]}, {pair[1]})"
    ax.plot(resolutions, bottom_filtered[pair], 
            marker='s', linewidth=2, markersize=4,
            label=label, linestyle='--', alpha=0.8)

# Add horizontal line at zero
ax.axhline(y=0, color='black', linestyle=':', linewidth=1, alpha=0.5, label='No change')

# Labels and styling
ax.set_xlabel('Resolution (m)', fontsize=14, fontweight='bold')
ax.set_ylabel('Filtered correlation', fontsize=14, fontweight='bold')
ax.set_title('Filtered Correlations \nTop 5 & Bottom 5 Species Pairs from Resolution 10m', 
             fontsize=16, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend(loc='best', fontsize=10, ncol=2, framealpha=0.9)

# Tight layout
plt.tight_layout()

# Save filtered plot
filtered_filename = 'variations/filtered_correlation_across_resolutions.png'
plt.savefig(filtered_filename, dpi=300, bbox_inches='tight')
print(f"\n✓ Filtered correlation tracking plot saved: {filtered_filename}")

plt.show()

