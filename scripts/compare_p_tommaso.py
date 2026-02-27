import numpy as np
import matplotlib.pyplot as plt
from src.config import load_config

# ========================================
# Load Configuration
# ========================================

config = load_config()
num_species = config['analysis']['num_species']
path_template = config['forests']['templates']['path_template']
names_template = config['forests']['templates']['names_template']
forest = config['forests']['name']

# --- Load matrices ---
P = np.load('precision_matrix_analysis_shrinkage/P_partial_correlation.npy')
T = np.loadtxt('correlations_tommaso/correlations.txt')

# Load species names
names_file = f"{path_template.format(forest=forest)}{names_template.format(forest=forest, census=4)}"
with open(names_file, 'r', encoding='utf-8-sig') as f:
    species_names = [line.strip() for line in f.readlines()][:num_species]

assert P.shape == T.shape, f"Shape mismatch: P={P.shape}, T={T.shape}"
assert P.shape[0] == len(species_names), "Matrix size != num_species"

np.fill_diagonal(P, np.nan)
np.fill_diagonal(T, np.nan)

# --- Top 15 pairs helper ---
def top_pairs(mat, names, n=15):
    mask = np.triu(np.ones_like(mat, dtype=bool), k=1)
    vals = mat[mask]
    idxs = np.argwhere(mask)
    order = np.argsort(vals)[-n:][::-1]
    return [(names[i], names[j], vals[k]) for k, (i, j) in enumerate(idxs) if k in set(order)]

def top_pairs(mat, names, n=15):
    mask = np.triu(np.ones_like(mat, dtype=bool), k=1)
    vals = mat[mask]
    idxs = np.argwhere(mask)
    pos_mask = vals > 0
    pos_vals = vals[pos_mask]
    pos_idxs = idxs[pos_mask]
    order = np.argsort(pos_vals)[-n:][::-1]
    return [(names[i], names[j], pos_vals[k]) for k, (i, j) in enumerate(pos_idxs) if k in set(order)]

top_P = top_pairs(P, species_names)
top_T = top_pairs(T, species_names)

pairs_P = {(a, b) for a, b, _ in top_P}
pairs_T = {(a, b) for a, b, _ in top_T}
common  = pairs_P & pairs_T

# --- Print top 15 side by side ---
print(f"\n{'RANK':<5} {'P matrix':<65} {'Tommaso':<65}")
print("-" * 135)
for rank, (p, t) in enumerate(zip(top_P, top_T), 1):
    p_str = f"{p[0][:20]:<22} <-> {p[1][:20]:<22} {p[2]:>8.5f}"
    t_str = f"{t[0][:20]:<22} <-> {t[1][:20]:<22} {t[2]:>8.5f}"
    flag  = "  ***" if (p[0], p[1]) in common else ""
    print(f"{rank:<5} {p_str:<65} {t_str:<65}{flag}")

print(f"\nCommon pairs in top 15: {len(common)}/{15}")
for a, b in sorted(common):
    print(f"  {a} <-> {b}")

# --- Scatter plot of all shared pairs ---
mask = np.triu(np.ones_like(P, dtype=bool), k=1)
p_vals = P[mask]
t_vals = T[mask]

# remove nans
valid = ~(np.isnan(p_vals) | np.isnan(t_vals))
p_vals = p_vals[valid]
t_vals = t_vals[valid]

fig, ax = plt.subplots(figsize=(3.5, 3.5), dpi=300)
ax.scatter(t_vals, p_vals, s=1.5, alpha=0.3, color='#1f77b4',
           linewidths=0, rasterized=True, label='All pairs')

# highlight common top-15 pairs
for a, b in common:
    i, j = species_names.index(a), species_names.index(b)
    ax.scatter(T[i, j], P[i, j], s=20, color='#d62728',
               zorder=5, linewidths=0)

lims = [min(t_vals.min(), p_vals.min()), max(t_vals.max(), p_vals.max())]
ax.plot(lims, lims, 'k--', linewidth=0.6, alpha=0.6)
ax.set_xlabel("Tommaso's matrix", fontsize=8)
ax.set_ylabel('P (partial correlation)', fontsize=8)
ax.tick_params(labelsize=7, which='both', direction='in',
               top=True, right=True, width=0.5, length=2)
for sp in ax.spines.values():
    sp.set_linewidth(0.5)

# pearson r
r = np.corrcoef(p_vals, t_vals)[0, 1]
ax.text(0.05, 0.92, f'r = {r:.3f}', transform=ax.transAxes,
        fontsize=7, va='top')
ax.text(0.05, 0.84, '● common top-15', transform=ax.transAxes,
        fontsize=6, color='#d62728', va='top')

fig.tight_layout()
fig.savefig('figures/P_vs_tommaso_scatter.svg',
            format='svg', dpi=300, bbox_inches='tight')
plt.close(fig)
print("\n✓ Saved figures/P_vs_tommaso_scatter.svg")
