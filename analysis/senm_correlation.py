import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from variables import *


chosen_file = spatial_file_template.format(
            nx=nx,  
            ny=ny,
            nu=nu,
            kernel=kernel,
            realization=realization
)

# Load data (space/tab separated: x y species_id)
data = np.loadtxt(simulations_path + chosen_file)  
df = pd.DataFrame(data, columns=['x','y','species_id'])

# --- Get top 100 species ---
species_counts = df['species_id'].value_counts()
top_100_species = species_counts.head(100).index.values
df_top100 = df[df['species_id'].isin(top_100_species)].copy()

# --- Binning and Matrix Creation ---
n_bins_x, n_bins_y = 15, 15  # 225 total boxes
grid_size = 512

df_top100['x_bin'] = (df_top100['x'] / (grid_size/n_bins_x)).astype(int).clip(0, n_bins_x-1)
df_top100['y_bin'] = (df_top100['y'] / (grid_size/n_bins_y)).astype(int).clip(0, n_bins_y-1)

# Create abundance matrix (100 species x 225 boxes)
abundance = np.zeros((100, n_bins_x * n_bins_y))
for i, sp in enumerate(top_100_species):
    box_counts = df_top100[df_top100['species_id'] == sp].groupby(['x_bin','y_bin']).size()
    for (x,y), count in box_counts.items():
        abundance[i, x*n_bins_y + y] = count

# --- REAL Correlation Matrix ---
real_corr = np.corrcoef(abundance)
real_corr = np.nan_to_num(real_corr, nan=0)


# --- RANDOMIZED Matrix ---
df_rand = df_top100.copy()
df_rand['species_id'] = np.random.permutation(df_top100['species_id'])



rand_abundance = np.zeros((100, n_bins_x * n_bins_y))
for i, sp in enumerate(top_100_species):
    box_counts = df_rand[df_rand['species_id'] == sp].groupby(['x_bin','y_bin']).size()
    for (x,y), count in box_counts.items():
        rand_abundance[i, x*n_bins_y + y] = count

rand_corr = np.corrcoef(rand_abundance)
rand_corr = np.nan_to_num(rand_corr, nan=0)

# Set diagonal equal to 0
for i in range(100):
    real_corr[i][i] = 0
    rand_corr[i][i] = 0



# --- Visualization ---
fig, axes = plt.subplots(1, 3, figsize=(18,5))

# Heatmaps
axes[0].imshow(real_corr, vmin=-1, vmax=1, cmap='coolwarm')
axes[0].set_title('Real Correlations (Top 100)')

axes[1].imshow(rand_corr, vmin=-1, vmax=1, cmap='coolwarm')
axes[1].set_title('Randomized Correlations')

# Eigenvalue spectrum
real_eig = np.linalg.eigvalsh(real_corr)
rand_eig = np.linalg.eigvalsh(rand_corr)
axes[2].plot(np.sort(real_eig)[::-1], 'o-', label='Real')
axes[2].plot(np.sort(rand_eig)[::-1], 'o-', label='Randomized')
axes[2].set_yscale('log')
axes[2].set_title('Eigenvalue Spectrum')
axes[2].legend()

plt.tight_layout()
plt.show()

# --- Key Statistics ---
print(f"Top 100 Species Analysis:")
print(f"Mean |Real Corr|: {np.abs(real_corr[np.triu_indices(100, k=1)]).mean():.3f}")
print(f"Mean |Random Corr|: {np.abs(rand_corr[np.triu_indices(100, k=1)]).mean():.3f}")
print(f"Max Real Eigenvalue: {np.max(real_eig):.2f} vs Random: {np.max(rand_eig):.2f}")