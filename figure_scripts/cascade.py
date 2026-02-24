import numpy as np
import matplotlib.pyplot as plt
from src.config import load_config

# === Config & Style ===
config = load_config()
forest = config['forests']['name']
num_species = config['analysis']['num_species']
resolution = 20

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 0.5

# === Load species names ===
path_template = config['forests']['templates']['path_template']
names_template = config['forests']['templates']['names_template']
names_file = f"{path_template.format(forest=forest)}{names_template.format(forest=forest, census=4)}"
with open(names_file, encoding='utf-8-sig') as f:
    species_names = [line.strip() for line in f if line.strip()][:num_species]

sp_dict = {name: idx for idx, name in enumerate(species_names)}

# === Load matrices ===
input_dir = "stripped_correlation_analysis"
C_raw = np.load(f'{input_dir}/unfiltered_raw_correlation_{resolution}m.npy')
C_raw_MP = np.load(f'{input_dir}/raw_correlation_{resolution}m.npy')
C_stripped = np.load(f'{input_dir}/stripped_correlation_{resolution}m.npy')
P_partial = np.load(f'precision_matrix_analysis_shrinkage/P_partial_correlation.npy')

# === Pairs ===
pairs_names = [('hirttr', 'guatdu'), ('hirttr', 'licapl'), ('guatdu', 'licapl')]
pairs_full_names = [('H. triandra', 'G. lucens'), 
                    ('H. triandra', 'L. platypus'), 
                    ('G. lucens', 'L. platypus')]
pairs_idx = [(sp_dict[a], sp_dict[b]) for a, b in pairs_names]

# === Extract correlations for each pair at each stage ===
stages = ['Raw', 'MP filtered', 'Nutrient\nregression', 'Partial corr.']
matrices = [C_raw, C_raw_MP, C_stripped, P_partial]

pair_trajectories = []
for i, j in pairs_idx:
    traj = [M[i, j] for M in matrices]
    pair_trajectories.append(traj)

# === Figure ===
fig, ax = plt.subplots(figsize=(3, 3), dpi=300)

colors = ['#d62728', '#2ca02c', '#1f77b4']

for idx, (traj, full_name) in enumerate(zip(pair_trajectories, pairs_full_names)):
    ax.plot(range(len(stages)), traj, 'o-', color=colors[idx], 
            linewidth=1.5, markersize=5, markerfacecolor='white', 
            markeredgewidth=1, label=f'{full_name[0]} – {full_name[1]}')

# Horizontal line at zero
ax.axhline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

ax.set_xticks(range(len(stages)))
ax.set_xticklabels(stages, fontsize=8, rotation=45, ha='right')
ax.set_ylabel('Correlation', fontsize=8)
ax.tick_params(labelsize=7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='y', alpha=0.3, linewidth=0.3)
ax.legend(fontsize=6, frameon=False, loc='best')

# === Print summary ===
print("\n" + "="*80)
print("CORRELATION TRAJECTORIES THROUGH FILTERING PIPELINE")
print("="*80)
for idx, (pair_name, traj) in enumerate(zip(pairs_full_names, pair_trajectories)):
    print(f"\n{pair_name[0]} – {pair_name[1]}:")
    for stage, val in zip(stages, traj):
        print(f"  {stage:20s}: {val:+.4f}")
    total_change = traj[-1] - traj[0]
    pct_total = 100 * total_change / (abs(traj[0]) + 1e-10)
    print(f"  {'Total change':20s}: {total_change:+.4f} ({pct_total:+.1f}%)")
print("="*80 + "\n")

# === Save ===
plt.tight_layout()
plt.savefig('figures/filtering_trajectories.svg', dpi=600, bbox_inches='tight')
plt.show()
