import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.stats import gaussian_kde
from variables import *
from functions import get_forest_file

# Array of bin sizes to test
bin_numbers = [10, 20, 40, 80, 160,320]  

forest = 'wanang'
census = 1

path = path_template.format(forest = forest, census = census)

all_senm_spectra = []
all_forest_spectra = []
box_labels = []


for n_bins_x in bin_numbers:
    n_bins_y = n_bins_x  # Keep square bins
    boxlen = 512/n_bins_x  
    
    print(f'Considering {forest} census {census}')
    print(f'\nAnalyzing with box length: {boxlen:.1f}m)')
   

    # Initialize eigenvalue matrix properly
    num_realizations = 100
    eig_matrix = np.zeros((num_realizations, 100))

    for realization in range(num_realizations):
        chosen_file = spatial_file_template.format(
                    nx=nx,  
                    ny=ny,
                    nu=nu,
                    kernel=kernel,
                    realization=realization+1
        )

        # Load data
        data = np.loadtxt(simulations_path + chosen_file)  
        df = pd.DataFrame(data, columns=['x','y','species_id'])

        # Get top 100 species
        species_counts = df['species_id'].value_counts()
        top_100_species = species_counts.head(100).index.values
        df_top100 = df[df['species_id'].isin(top_100_species)].copy()

        # Binning
        grid_size = 512
        df_top100['x_bin'] = (df_top100['x'] / (grid_size/n_bins_x)).astype(int).clip(0, n_bins_x-1)
        df_top100['y_bin'] = (df_top100['y'] / (grid_size/n_bins_y)).astype(int).clip(0, n_bins_y-1)

        # Create abundance matrix
        abundance = np.zeros((100, n_bins_x * n_bins_y))
        for i, sp in enumerate(top_100_species):
            box_counts = df_top100[df_top100['species_id'] == sp].groupby(['x_bin','y_bin']).size()
            for (x,y), count in box_counts.items():
                abundance[i, x*n_bins_y + y] = count

        # Correlation matrix
        real_corr = np.corrcoef(abundance)
        real_corr = np.nan_to_num(real_corr, nan=0)
        
        # Get sorted eigenvalues
        eig_vals = np.linalg.eigvalsh(real_corr)
        eig_vals_sorted = np.sort(eig_vals)[::-1]
        
        # Store in matrix
        eig_matrix[realization, :] = eig_vals_sorted

    # Calculate mean spectrum
    mean_eig = np.mean(eig_matrix, axis=0)
    std_eig = np.std(eig_matrix, axis=0)

    # FOREST ANALYSIS
    forest_path = path_template.format(forest=forest)
    census_file = forest_path + census_template.format(forest=forest, census=census)
    df = pd.DataFrame(pd.read_csv(census_file))

    num_boxes = n_bins_x * n_bins_y
    num_species = 100

    # Create bins
    x_bins = pd.cut(df.x, bins=n_bins_x, labels=False)
    y_bins = pd.cut(df.y, bins=n_bins_y, labels=False)

    # Load names
    names_file = forest_path + names_template.format(forest=forest, census=census)
    try:
        names = np.loadtxt(names_file, dtype='str', ndmin=1)[:num_species]
        names = [name[0] if isinstance(name, np.ndarray) else str(name) for name in names]
    except ValueError:
        with open(names_file, 'r') as f:
            names = [line.strip() for line in f if line.strip()][:num_species]

    # Create species matrix
    species_matrix = np.empty([num_species, num_boxes])
    for i, name in enumerate(names):
        try:
            mask = df['name'] == name
            grid_counts = np.zeros((n_bins_x, n_bins_y), dtype=int)
            for x, y in zip(x_bins[mask], y_bins[mask]):
                if not np.isnan(x) and not np.isnan(y):
                    grid_counts[int(x), int(y)] += 1
            species_matrix[i, :] = grid_counts.ravel()
        except Exception as e:
            print(f"Error processing species {name}: {e}")
            species_matrix[i, :] = 0

    # Calculate forest eigenvalues
    raw_correlation_matrix = np.corrcoef(species_matrix)
    forest_eigenvalues = np.linalg.eigh(raw_correlation_matrix)[0][::-1]

    # Plotting
    plt.figure(figsize=(10,6))
    plt.loglog(np.arange(1,101), mean_eig, 'o-', label='SENM Data', markersize=5)
    plt.loglog(np.arange(1,101), forest_eigenvalues, 'o-', label=f'{forest} Data', markersize=5)
    
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Eigenvalue Rank')
    plt.ylabel('Eigenvalue Magnitude')
    plt.legend()
    plt.title(f'{n_bins_x}x{n_bins_y} Bins ({num_realizations} Realizations)')
    
    # Add text box with box length information
    textstr = f'Box length: {boxlen:.1f}m'
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    plt.text(0.95, 0.95, textstr, transform=plt.gca().transAxes, 
             fontsize=12, verticalalignment='top', horizontalalignment='right',
             bbox=props)
    
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.tight_layout()
    
    # Save each plot
    plot_filename = path + f'plots/eigenvalue_spectrum_{boxlen}.png'
    plt.savefig(plot_filename, dpi=300)
    print(f'Saved plot as {plot_filename}')
    
    # Append
    all_senm_spectra.append(mean_eig)
    all_forest_spectra.append(forest_eigenvalues)
    box_labels.append(f'{boxlen}')

    plt.close()  # Close the figure to free memory

    
    # Print statistics
    print(f"\nStatistics for {n_bins_x}x{n_bins_y} bins:")
    print(f"Box length: {boxlen:.1f}m")
    print(f"Mean largest SENM eigenvalue: {mean_eig[0]:.3f} ± {std_eig[0]:.3f}")
    print(f"Mean smallest SENM non-zero eigenvalue: {mean_eig[-2]:.3f} ± {std_eig[-2]:.3f}")
    

# Plot all SENM spectra
plt.figure(figsize=(10, 6))
for spectrum, label in zip(all_senm_spectra, box_labels):
    plt.loglog(range(1, 101), spectrum, label=label)
plt.xlabel('Eigenvalue Rank')
plt.ylabel('Eigenvalue Magnitude')
plt.title('SENM Spectra for Different box Sizes')
plt.legend(title='Box Size')
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout()
plt.savefig(path + 'plots/all_senm_spectra.png', dpi=300)
plt.close()
print("Saved combined SENM spectra plot as 'all_senm_spectra.png'")


# Plot all Wanang spectra
plt.figure(figsize=(10, 6))
for spectrum, label in zip(all_forest_spectra, box_labels):
    plt.loglog(range(1, 101), spectrum, label=label)
plt.xlabel('Eigenvalue Rank')
plt.ylabel('Eigenvalue Magnitude')
plt.title('Wanang Spectra for Different Box Sizes')
plt.legend(title='Box Size')
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout()
plt.savefig(path + 'plots/all_wanang_spectra.png', dpi=300)
plt.close()
print("Saved combined Wanang spectra plot as 'all_wanang_spectra.png'")

