import pandas as pd
from variables import *

forest = 'wanang'
census = '1'

path = path_template.format(forest=forest)

# Input and output file paths
input_file = path + census_template.format(forest=forest, census=census)  
coord_output_file = input_file
abundance_output_file = path + f"names_{forest}_{census}.txt"

# Read the data
df = pd.read_csv(input_file, header=0, skip_blank_lines=True)

# Process coordinates
coord_df = df[['mnemonic', 'px', 'py']].rename(columns={
    'mnemonic': 'name',
    'px': 'x',
    'py': 'y'
}).dropna()
coord_df = coord_df[(coord_df['x'] != 0) & (coord_df['y'] != 0)]

# Convert to absolute coordinates starting from bottom-left (0,0)
min_x, min_y = coord_df['x'].min(), coord_df['y'].min()
coord_df['x'] = (coord_df['x'] - min_x).round(2)
coord_df['y'] = (coord_df['y'] - min_y).round(2)

# Filter to keep only points within 512x512 grid
coord_df = coord_df[(coord_df['x'] < 512) & (coord_df['y'] < 512)]

# Reset index after filtering
coord_df.reset_index(drop=True, inplace=True)

# Save coordinate data
coord_df.to_csv(coord_output_file, index=False)

# Create and save abundance-sorted mnemonics (using original data)
abundance_df = df['mnemonic'].value_counts().reset_index()
abundance_df.columns = ['mnemonic', 'count']
abundance_df = abundance_df.sort_values('count', ascending=False)

# Save just the mnemonics in abundance order
abundance_df['mnemonic'].to_csv(abundance_output_file, 
                               index=False, 
                               header=False)

print(f"Filtered coordinate data saved to {coord_output_file}")
print(f"Mnemonics sorted by abundance saved to {abundance_output_file}")
print(f"Top 5 most abundant mnemonics:")
print(abundance_df.head())
print(f"\nData bounds after filtering:")
print(f"X range: {coord_df['x'].min():.2f} - {coord_df['x'].max():.2f}")
print(f"Y range: {coord_df['y'].min():.2f} - {coord_df['y'].max():.2f}")
print(f"Total points within 512x512 grid: {len(coord_df)}")