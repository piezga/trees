import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from src.config import load_config
import os

# Load config
config = load_config()

# Parameters
N = 100  # Number of top species to consider
forest = 'barro'
censi = config['forests']['censuses']

# Paths and templates
path_template = config['forests']['templates']['path_template']
census_template = config['forests']['templates']['census_template']
forest_path = path_template.format(forest=forest)

senm_path = config['senm_templates']['path']
senm_abundance_template = config['senm_templates']['abundance']
nx = config['senm']['nx']
ny = config['senm']['ny']
nu = config['senm']['nu']
kernel = config['senm']['kernel']
num_realizations = config['senm']['num_realizations']

# --- Load empirical abundances ---
empirical_curves = []
for census in censi:
    file_path = os.path.join(forest_path, census_template.format(forest=forest, census=census))
    df = pd.read_csv(file_path)
    counts = df['name'].value_counts().values
    top_counts = np.sort(counts)[::-1][:N]
    top_rel = top_counts / top_counts.sum()
    empirical_curves.append(top_rel)

# --- Load SENM abundances and average ---
all_senm_abund = []
for realization in range(1, num_realizations + 1):
    file_name = senm_abundance_template.format(
        nx=nx, ny=ny, nu=nu, kernel=kernel, realization=realization)
    full_path = os.path.join(senm_path, file_name)
    try:
        data = np.loadtxt(full_path)[:, 1]  # column 1 = abundance
    except Exception as e:
        print(f"Error reading {full_path}: {e}")
        continue
    sorted_abund = np.sort(data)[::-1][:N]
    padded = np.pad(sorted_abund, (0, max(0, N - len(sorted_abund))), constant_values=0)
    all_senm_abund.append(padded)

senm_mean = np.mean(all_senm_abund, axis=0)
senm_curve = senm_mean / senm_mean.sum()

# --- Plotting ---
plt.figure(figsize=(10, 6))
x = np.arange(1, N + 1)

# Empirical
for i, curve in enumerate(empirical_curves):
    plt.plot(x[:len(curve)], curve, label=f'Census {censi[i]}', marker='o', linestyle='-')

# SENM
plt.plot(x, senm_curve, label='SENM (mean)', color='black', linewidth=2.5, linestyle='--', marker='x')

# Plot formatting
plt.xlabel('Species Rank', fontsize=13, weight='bold')
plt.ylabel('Relative Abundance', fontsize=13, weight='bold')
plt.title(f'Relative Abundance Curves (Top {N} Species)', fontsize=14, weight='bold')
plt.yscale('log')
plt.legend()
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()

