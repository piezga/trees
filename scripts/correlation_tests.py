import numpy as np
import argparse
import matplotlib
import matplotlib.pyplot as plt

from src.config import load_config
from src.compute import (
        compute_spectra, MarchenkoPastur
        )

matplotlib.use('TkAgg')  # Switch from Qt to Tkinter backend

# === Load config ===
config = load_config()
num_species = config['analysis']['num_species']

# === Pick resolution ===
resolution = 25
print(f'Testing for resolution {resolution} m ')

(senm_mean, senm_std, forest_spectra, 
 bins, senm_abundance, forest_abundance) = compute_spectra(resolution, 
                                                           calculate=True)

# Check standardization

print('\nSENM STANDARDIZATION')
print('\nSENM mean')
print(np.mean(senm_abundance, axis = 1))
print('\nSENM std')
print(np.std(senm_abundance, axis = 1))

print('\nFOREST STANDARDIZATION')
print('\nForest mean')
print(np.mean(forest_abundance, axis = 1))
print('\nForest std')
print(np.std(forest_abundance, axis = 1))

# Is the shape right?
print('\n\n')
print('Testing for correct shape (column)')
print(f'Forest shape is : {np.shape(forest_abundance)}')
print(f'Expecting column: {500*1000/(resolution**2)}')
print(f'SENM shape is   : {np.shape(senm_abundance)}')
print(f'Expecting column: {500*500/(resolution**2)}')

# === Calculating correlation matrices ===

senm_corr = np.corrcoef(senm_abundance)
forest_corr = np.corrcoef(forest_abundance)

print('\n\n')
print('The correlation matrices have shape')
print(f'Forest: {forest_corr.shape}')
print(f'SENM: {senm_corr.shape}')
print(f'\nForest preview:\n {forest_corr}')
print(f'\nSENM preview:\n {senm_corr}')

# === MarchenkoPastur Check ===
print('\n\n')
print('Here I am keeping the small eigenvalues:')
filtered_senm_corr = MarchenkoPastur(senm_corr, num_species, 
                                     bins[0]*bins[1], remove_largest=False, 
                                     remove_small=False )
filtered_forest_corr = MarchenkoPastur(forest_corr, num_species, 
                                       bins[2]*bins[3], remove_largest=False,
                                       remove_small=False)
print(f'\nSenm: \n {filtered_senm_corr}')
print(f'\nForest: \n {filtered_forest_corr}')

print('\n\n')
print('Here I am removing the small eigenvalues:')
filtered_senm_corr2 = MarchenkoPastur(senm_corr, num_species, 
                                     bins[0]*bins[1], remove_largest=False, 
                                     remove_small=True )
filtered_forest_corr2 = MarchenkoPastur(forest_corr, num_species, 
                                       bins[2]*bins[3], remove_largest=False,
                                       remove_small=True)
print(f'\nSenm: \n {filtered_senm_corr2}')
print(f'\nForest: \n {filtered_forest_corr2}')

# === Plotting the difference ===

# Create the plot with 2 rows and 2 columns
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot the original SENM correlation matrix
axes[0, 0].imshow(senm_corr, cmap='coolwarm', aspect='auto')
axes[0, 0].set_title('Filtered SENM (small eigenvalues included)')
axes[0, 0].set_xlabel('Index')
axes[0, 0].set_ylabel('Index')
plt.colorbar(axes[0, 0].imshow(filtered_senm_corr2, cmap='coolwarm', aspect='auto'), ax=axes[0, 0])

# Plot the original Forest correlation matrix
axes[0, 1].imshow(forest_corr, cmap='coolwarm', aspect='auto')
axes[0, 1].set_title('Filtered Forest (small eigenvalues included)')
axes[0, 1].set_xlabel('Index')
axes[0, 1].set_ylabel('Index')
plt.colorbar(axes[0, 1].imshow(filtered_forest_corr2, cmap='coolwarm', aspect='auto'), ax=axes[0, 1])

# Plot the filtered SENM correlation matrix (no small eigenvalues)
axes[1, 0].imshow(filtered_senm_corr, cmap='coolwarm', aspect='auto')
axes[1, 0].set_title('Filtered SENM (no small eigenvalues)')
axes[1, 0].set_xlabel('Index')
axes[1, 0].set_ylabel('Index')
plt.colorbar(axes[1, 0].imshow(filtered_senm_corr, cmap='coolwarm', aspect='auto'), ax=axes[1, 0])

# Plot the filtered Forest correlation matrix (no small eigenvalues)
axes[1, 1].imshow(filtered_forest_corr, cmap='coolwarm', aspect='auto')
axes[1, 1].set_title('Filtered Forest (no small eigenvalues)')
axes[1, 1].set_xlabel('Index')
axes[1, 1].set_ylabel('Index')
plt.colorbar(axes[1, 1].imshow(filtered_forest_corr, cmap='coolwarm', aspect='auto'), ax=axes[1, 1])

# Adjust layout to make sure the plots are spaced nicely
plt.tight_layout()
plt.show()
