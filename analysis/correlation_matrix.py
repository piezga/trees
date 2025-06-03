import pandas as pd
from variables import *
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# File to analyse


forest = 'barro'
census = 8

path = path_template.format(forest = forest, census = census)

barro_names_file = path + 'names_barro_8.txt'

forest_path = path_template.format(forest=forest)
census_file = forest_path + census_template.format(forest=forest, census=census)
df = pd.DataFrame(pd.read_csv(census_file))

# Parameters
n_bins_x = 20 
n_bins_y = 10

num_boxes = n_bins_x * n_bins_y
num_species = 100

# Create bins

x_bins = pd.cut(df.x, bins=n_bins_x, labels=False)
y_bins = pd.cut(df.y, bins=n_bins_y, labels=False)


# Select N most abundant names

names = np.loadtxt(barro_names_file, dtype = 'str').tolist()[:num_species]


# Preallocation

species_matrix = np.empty([num_species, num_boxes])


# Loop over each species and count grid occurrences
for i, name in enumerate(names):
    # Filter data for the current species
    mask = df['name'] == name
    # Initialize a temporary grid counter
    grid_counts = np.zeros((n_bins_x, n_bins_y), dtype=int)
    
    # Count specimens in each grid box
    for x, y in zip(x_bins[mask], y_bins[mask]):
        if not np.isnan(x) and not np.isnan(y):  # Skip NaN bins (outliers)
            grid_counts[int(x), int(y)] += 1
    
    # Unfurl grid into a 200-length vector (row-major order)
    species_matrix[i, :] = grid_counts.ravel()

# Result: species_matrix is (n_species x 200)
print("Species Matrix Shape:", species_matrix.shape)
print(species_matrix)

# Calculate correlation


raw_correlation_matrix = np.corrcoef(species_matrix)

for i in range(100):
    raw_correlation_matrix[i][i] = 0
    

mean = 0
for i in range(100):
    mean += np.sum(raw_correlation_matrix[i][i+1:])
    
mean /= 50*99

print(f'mean correlation is : {mean}')

# Check matrix properties
print("Shape:", raw_correlation_matrix.shape)  # Should be (n_species, n_species)
print("Is symmetric?", np.allclose(raw_correlation_matrix, 
                                   raw_correlation_matrix.T))


plt.figure(figsize=(10, 8))
sns.heatmap(
    raw_correlation_matrix,      # Direct NumPy array input
    annot=False,              # Show values
    fmt=".2f",              # Round to 2 decimal places
    cmap="coolwarm",        # Red (positive) to Blue (negative)
    vmin=-1,                # Fix color scale limits
    vmax=1,
    linewidths=0.5,         # Add grid lines
    square=True             # Keep aspect ratio square
)



plt.title("Raw Species Correlation Matrix", fontsize=14)
plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels
plt.tight_layout()
plt.show()


# Preferred method for symmetric matrices (preferred for correlation matrices)
eigenvalues = np.linalg.eigh(raw_correlation_matrix)[0]  # Returns sorted eigenvalues


# Randomize species labels (shuffle the 'name' column)
df_randomized = df.copy()
df_randomized['name'] = np.random.permutation(df['name'])  # Shuffle labels


# Verify randomization
print("Original species counts:\n", df['name'].value_counts())
print("\nRandomized species counts:\n", df_randomized['name'].value_counts())


random_species_matrix = np.empty([num_species, num_boxes])


# Loop over each species and count grid occurrences
for i, name in enumerate(names):
    # Filter data for the current species
    mask = df_randomized['name'] == name
    # Initialize a temporary grid counter
    grid_counts = np.zeros((n_bins_x, n_bins_y), dtype=int)
    
    # Count specimens in each grid box
    for x, y in zip(x_bins[mask], y_bins[mask]):
        if not np.isnan(x) and not np.isnan(y):  # Skip NaN bins (outliers)
            grid_counts[int(x), int(y)] += 1
    
    # Unfurl grid into a 200-length vector (row-major order)
    random_species_matrix[i, :] = grid_counts.ravel()

# Result: random_species_matrix is (n_species x 200)
print("Species Matrix Shape:", random_species_matrix.shape)
print(random_species_matrix)

# Calculate correlation


randomized_correlation_matrix = np.corrcoef(random_species_matrix)

randomized_eigenvalues = np.linalg.eigh(randomized_correlation_matrix)[0]  # Returns sorted eigenvalues

# Sort eigenvalues in descending order
sorted_eigs = np.sort(eigenvalues)[::-1]
sorted_rand_eigs = np.sort(randomized_eigenvalues)[::-1]

# Create rank-order x-axis
ranks = np.arange(1, len(sorted_eigs)+1)

# Create the log-log plot
plt.figure(figsize=(10, 6))
plt.loglog(ranks, sorted_eigs, 'o-', label='Real Data', markersize=5)
plt.loglog(ranks, sorted_rand_eigs, 's-', label='Randomized', markersize=5, alpha=0.7)

# Add labels and legend
plt.xlabel('Rank (log scale)')
plt.ylabel('Eigenvalue (log scale)')
plt.title('Eigenvalue Distribution (Log-Log Scale)')
plt.legend()
plt.grid(True, which="both", ls="--")

# Add theoretical random matrix expectation (Marchenko-Pastur law)
# Only if you want to compare to null expectation
if len(eigenvalues) > 10:  # Only for large matrices
    Q = len(eigenvalues)/num_boxes  # Matrix shape ratio
    lambda_plus = (1 + np.sqrt(1/Q))**2  # MP upper bound
    lambda_minus = (1 - np.sqrt(1/Q))**2  # MP lower bound
    plt.axhline(y=lambda_plus, color='r', linestyle='--', label='MP upper bound')
    plt.axhline(y=lambda_minus, color='g', linestyle='--', label='MP lower bound')

plt.tight_layout()
plt.show()
