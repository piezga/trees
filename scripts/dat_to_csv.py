import os
import pandas as pd
from src.config import load_config

config = load_config()

forest = config['forests']['name'][0]  # Assuming single forest
censuses = config['forests']['censuses']
path_template = config['forests']['templates']['path_template']
raw_data_template = config['forests']['templates']['raw_data_template']
census_template = config['forests']['templates']['census_template']
names_template = config['forests']['templates']['names_template']

for census in censuses:
    # Resolve full paths
    path = path_template.format(forest=forest)
    input_file = os.path.join(path, raw_data_template.format(census=census))
    output_file = os.path.join(path, census_template.format(forest=forest, census=census))
    names_file = os.path.join(path, names_template.format(forest=forest, census=census))

    # Read raw data
    data = []
    with open(input_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 4:
                x, y, _, name = parts
                data.append((float(x), float(y), name))

    # Create DataFrame
    df = pd.DataFrame(data, columns=['x', 'y', 'name'])

    # Save cleaned CSV
    df.to_csv(output_file, index=False)
    print(f"Saved cleaned CSV to: {output_file}")

    # Count species and sort by descending abundance
    species_counts = df['name'].value_counts()
    sorted_species = species_counts.index.tolist()

    # Save species names to text file
    with open(names_file, 'w') as f:
        for name in sorted_species:
            f.write(f"{name}\n")

    print(f"Saved species name list to: {names_file}")

