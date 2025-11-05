import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- Function under test ---
def compute_forest_spectrum(df, names, n_bins_x, n_bins_y):
    """Compute eigenvalue spectrum for forest data"""
    x_bins = pd.cut(df.x, bins=n_bins_x, labels=False)
    y_bins = pd.cut(df.y, bins=n_bins_y, labels=False)
    species_matrix = np.array([
        np.bincount(x_bins[df.name == name] * n_bins_y + y_bins[df.name == name],
                    minlength=n_bins_x * n_bins_y)
        for name in names
    ])
    corr = np.nan_to_num(np.corrcoef(species_matrix), nan=0)
    return np.sort(np.linalg.eigvalsh(corr))[::-1]

# --- Generate random test data ---
def generate_random_forest_data(num_species=10, num_points=1000, 
                                x_range=(0, 1000), y_range=(0, 500)):
    """Generate a random forest-like dataset with random species and coordinates."""
    np.random.seed(42)
    names = [f"species_{i+1}" for i in range(num_species)]
    df = pd.DataFrame({
        "x": np.random.uniform(x_range[0], x_range[1], num_points),
        "y": np.random.uniform(y_range[0], y_range[1], num_points),
        "name": np.random.choice(names, num_points)
    })
    return df, names

# --- Marchenko–Pastur distribution ---
def marchenko_pastur_pdf(lmbda, q, sigma2=1.0):
    """
    Compute the Marchenko–Pastur probability density function.
    Parameters:
      - lmbda: array of eigenvalue positions
      - q: aspect ratio (N_samples / N_features)
      - sigma2: variance scale
    """
    l_min = sigma2 * (1 - np.sqrt(1.0 / q))**2
    l_max = sigma2 * (1 + np.sqrt(1.0 / q))**2
    pdf = np.zeros_like(lmbda)
    mask = (lmbda >= l_min) & (lmbda <= l_max)
    pdf[mask] = q / (2 * np.pi * sigma2 * lmbda[mask]) * np.sqrt((l_max - lmbda[mask]) * (lmbda[mask] - l_min))
    return pdf, l_min, l_max

# --- Main test ---
if __name__ == "__main__":
    # Generate random data
    num_species = 100
    num_points = 200*10**3
    n_bins_x, n_bins_y = 20, 10

    df, names = generate_random_forest_data(num_species=num_species, num_points=num_points)

    # Compute eigenvalue spectrum
    spectrum = compute_forest_spectrum(df, names, n_bins_x, n_bins_y)

    # Compute MP distribution
    num_bins = n_bins_x * n_bins_y
    q = num_bins / num_species  # aspect ratio (samples/features)
    l_vals = np.linspace(0, 4, 400)
    mp_pdf, l_min, l_max = marchenko_pastur_pdf(l_vals, q)

    # --- Plot ---
    plt.figure(figsize=(7,5))
    plt.plot(spectrum, 'o-', label='Empirical eigenvalues (forest data)')
    plt.axhline(1, color='gray', linestyle='--', alpha=0.6, label='Mean correlation=1')
    plt.title("Eigenvalue Spectrum vs. Marchenko–Pastur Prediction")
    plt.xlabel("Eigenvalue Rank")
    plt.ylabel("Eigenvalue")
    plt.grid(alpha=0.3)
    plt.legend()

    plt.figure(figsize=(7,5))
    plt.hist(spectrum, bins=20, density=True, alpha=0.6, label="Empirical spectrum")
    plt.plot(l_vals, mp_pdf, 'r-', lw=2, label="Marchenko–Pastur PDF")
    plt.axvline(l_min, color='k', linestyle='--', alpha=0.6)
    plt.axvline(l_max, color='k', linestyle='--', alpha=0.6)
    plt.title("Marchenko–Pastur vs. Empirical Eigenvalue Distribution")
    plt.xlabel("Eigenvalue")
    plt.ylabel("Density")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

    print(f"Marchenko–Pastur theoretical range: λ₋={l_min:.3f}, λ₊={l_max:.3f}")

