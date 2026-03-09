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
C_raw      = np.load(f'{input_dir}/unfiltered_raw_correlation_{resolution}m.npy')
C_raw_MP   = np.load(f'{input_dir}/raw_correlation_{resolution}m.npy')
C_stripped = np.load(f'{input_dir}/stripped_correlation_{resolution}m.npy')
P_partial  = np.load(f'precision_matrix_analysis_shrinkage/P_partial_correlation.npy')

matrices = [C_raw, C_raw_MP, C_stripped, P_partial]

# === FIND PAIRS ===
mask = np.triu(np.ones(num_species, dtype=bool), k=1)
ii, jj = np.where(mask)


# Pair that drops sharply at MP step (C_raw > 0, largest drop C_raw -> C_raw_MP)
drop_MP = C_raw[ii, jj] - C_raw_MP[ii, jj]          # positive = dropped
started_positive_MP = C_raw[ii, jj] > 0

# Filter to valid pairs and sort drops
valid_drops = drop_MP[started_positive_MP]
valid_indices = np.where(started_positive_MP)[0]

# Get one of the largest drops
largest_idx = valid_indices[np.argsort(valid_drops)[-1]]  # -2 = second largest
mp_pair_idx = (ii[largest_idx], jj[largest_idx])
mp_pair_names = (species_names[mp_pair_idx[0]], species_names[mp_pair_idx[1]])


# Pair that drops sharply at nutrient regression step (C_raw_MP > 0, largest drop C_raw_MP -> C_stripped)
drop_nutr = C_raw_MP[ii, jj] - C_stripped[ii, jj]
started_positive_nutr = C_raw_MP[ii, jj] > 0
best_nutr = np.argmax(drop_nutr * started_positive_nutr)
nutr_pair_idx = (ii[best_nutr], jj[best_nutr])
nutr_pair_names = (species_names[nutr_pair_idx[0]], species_names[nutr_pair_idx[1]])

# Pair that drops sharply at partial step 
drop_partial = C_stripped[ii, jj] - P_partial[ii, jj]          # positive = dropped

started_positive_partial = C_stripped[ii, jj] > 0

# Filter to valid pairs and sort drops
valid_drops = drop_partial[started_positive_partial]
valid_indices = np.where(started_positive_partial)[0]

# Get one of the largest drops
largest_idx = valid_indices[np.argsort(valid_drops)[-1]]  # -2 = second largest
partial_pair_idx = (ii[largest_idx], jj[largest_idx])
partial_pair_names = (species_names[partial_pair_idx[0]], species_names[partial_pair_idx[1]])

print(f"Auto-selected | drops at MP step      : {mp_pair_names[0]} – {mp_pair_names[1]}")
print(f"  trajectory: {[M[mp_pair_idx] for M in matrices]}")
print(f"Auto-selected | drops at nutrient step : {nutr_pair_names[0]} – {nutr_pair_names[1]}")
print(f"  trajectory: {[M[nutr_pair_idx] for M in matrices]}")
print(f"Auto-selected | drops at partial step : {partial_pair_names[0]} – {partial_pair_names[1]}")
print(f"  trajectory: {[M[partial_pair_idx] for M in matrices]}")

mp_pair_names = ['S. simplex', 'H. purpurea']
nutr_pair_names = ['O. maporta', 'M. guianensis']
partial_pair_names = ['C. curvigemmia', 'P. marginata']

# === Find mediating species for (cou2cu, psycma) ===
i, j = partial_pair_idx

# For every other species k, take the minimum of its partial correlation with i and j
# The best mediator maximises this joint signal
scores = np.array([
    min(P_partial[i, k], P_partial[j, k])
    for k in range(num_species)
    if k != i and k != j
])
candidate_indices = [k for k in range(num_species) if k != i and k != j]

top_n = 5
top_idx = np.argsort(scores)[-top_n:][::-1]

print(f"\nTop mediator candidates for {partial_pair_names[0]} – {partial_pair_names[1]}:")
print(f"{'Species':<30}  P(i,k)    P(j,k)    min score")
print("-" * 65)
for idx in top_idx:
    k = candidate_indices[idx]
    p_ik = P_partial[i, k]
    p_jk = P_partial[j, k]
    print(f"{species_names[k]:<30}  {p_ik:+.4f}    {p_jk:+.4f}    {scores[idx]:+.4f}")


# === Pairs ===
hardcoded_names = [('hirttr', 'licapl')]
hardcoded_full  = [('H. triandra', 'L. platypus')]
hardcoded_idx   = [(sp_dict[a], sp_dict[b]) for a, b in hardcoded_names]

all_pairs_idx   = hardcoded_idx   + [mp_pair_idx,   nutr_pair_idx, partial_pair_idx]
all_pairs_label = hardcoded_full  + [mp_pair_names, nutr_pair_names, partial_pair_names]

# === Extract trajectories ===
stages = ['Raw', 'MP filtered', 'Nutrient\nregression + MP']
pair_trajectories = [[M[i, j] for M in matrices[:-1]] for i, j in all_pairs_idx]
# === Figure ===
fig, ax = plt.subplots(figsize=(3.5, 3.2), dpi=300)

colors = ['#0077BB', '#EE7733', '#33BBEE', '#CC3311']

for idx, (traj, full_name) in enumerate(
        zip(pair_trajectories, all_pairs_label)):
    ax.plot(range(len(stages)), traj, '--o',
            color=colors[idx], 
            linewidth=0.4, markersize=4,
            markerfacecolor='white', markeredgewidth=0.8,
            label=f'{full_name[0]} – {full_name[1]}')

ax.axhline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
ax.set_xticks(range(len(stages)))
ax.set_xticklabels(stages, fontsize=8, rotation=45, ha='right')
ax.set_ylabel('Correlation', fontsize=8)
ax.tick_params(labelsize=7, direction='in', top=True, right=False,
               width=0.5, length=2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(fontsize=5.5, frameon=False, loc='best')

plt.tight_layout()
plt.savefig('figures/filtering_trajectories.svg', dpi=600, bbox_inches='tight')
plt.savefig('figures/filtering_trajectories.pdf', dpi=600, bbox_inches='tight')
plt.show()
print("✓ Saved figures/filtering_trajectories.svg")
print("✓ Saved figures/filtering_trajectories.pdf")
