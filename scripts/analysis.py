import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from src.config import load_config
from src.functions import load_file_with_padding

# Load config
config = load_config()

# Extract from config
forest = config['forests']['name'][0]
censuses = config['forests']['censuses']
num_species = config['analysis']['num_species'][0]  # Only one value in list
chosen_kernels = config['senm']['chosen_kernels']
num_realizations = config['senm']['num_realizations']
nu = config['senm']['nu']
nx = config['senm']['nx']
ny = config['senm']['ny']

# Templates
path_template = config['forests']['templates']['path_template']
abundance_template = config['forests']['templates']['abundance_template']
senm_abundance_template = config['senm_templates']['abundance']
simulations_path = config['senm_templates']['path']

# =======================
# Load empirical dataset
# =======================

dataset = pd.DataFrame()
path = path_template.format(forest=forest)

for census in censuses:
    abundance_file = path + abundance_template.format(forest=forest, census=census)
    df = pd.read_csv(abundance_file)

    # Relative abundance
    rel_values = df['count'] / df['count'].sum()

    # Build column per census
    temp_df = pd.DataFrame({f'Census {census}': rel_values})

    # Merge
    dataset = temp_df if dataset.empty else dataset.join(temp_df, how='outer')

# Fill missing species with 0s
dataset = dataset.fillna(0)

# =======================
# Load SENM simulations
# =======================

senm = pd.DataFrame()

for kernel in chosen_kernels:
    print(f'Processing kernel {kernel}...')

    mean_senm_abundances = np.zeros(num_species)

    for realization in range(1, num_realizations + 1):
        current_file = senm_abundance_template.format(
            nx=nx, ny=ny, nu=nu, kernel=kernel, realization=realization
        )
        file_path = simulations_path + current_file

        # Load abundance data with padding
        abundances = load_file_with_padding(file_path, num_species, 2)[:, 1]

        # Sort in descending order
        abundances = np.sort(abundances)[::-1]

        mean_senm_abundances += abundances

    # Average and normalize
    mean_senm_abundances /= num_realizations
    relative_abundances = mean_senm_abundances / mean_senm_abundances.sum()

    senm[f'Kernel {kernel}'] = relative_abundances

# =======================
# Plotting
# =======================

fig, ax = plt.subplots(figsize=(12, 7))

# Plot a couple of censuses from dataset
selected_censuses = ['Census 1', 'Census 4']  # Adjust if needed
dataset[selected_censuses].plot(logy=True, style='o', markersize=5,
                                ax=ax, label=selected_censuses,
                                color=['blue', 'green'])

# Plot SENM
senm.plot(logy=True, style='.', markersize=4,
          ax=ax, label=senm.columns,
          color=['red', 'purple', 'orange'])

# Customize
ax.set_ylabel('Relative abundance (log scale)')
ax.set_xlabel('Species rank')
ax.set_title('Empirical vs Modeled Species Abundance Distributions')
ax.legend(fontsize=12)
plt.tight_layout()
plt.show()

