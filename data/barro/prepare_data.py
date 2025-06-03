import os
import pandas as pd
from collections import Counter

def process_dat_files(input_folder='.', output_folder='output', max_x=500, max_y=500):
    """Process all bcitree*.dat files with 500x500m filtering and create 3 output files."""
    os.makedirs(output_folder, exist_ok=True)
    
    dat_files = [f for f in os.listdir(input_folder) 
               if f.startswith('bcitree') and f.endswith('.dat')]
    
    for dat_file in dat_files:
        file_num = dat_file[7:-4]  # Gets X from bcitreeX.dat
        
        barro_data = []
        species_counts = Counter()
        
        with open(os.path.join(input_folder, dat_file), 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 4:
                    x, y = float(parts[0]), float(parts[1])
                    if x <= max_x and y <= max_y:
                        count, name = int(parts[2]), parts[3]
                        barro_data.append({'x': x, 'y': y, 'name': name})
                        species_counts[name] += count
        
        # Sort species by abundance (descending)
        sorted_species = [name for name, count in 
                         sorted(species_counts.items(), 
                         key=lambda x: x[1], 
                         reverse=True)]
        
        # 1. Save barro_X.csv (coordinates)
        pd.DataFrame(barro_data).to_csv(
            os.path.join(output_folder, f'barro_{file_num}.csv'), 
            index=False
        )
        
        # 2. Save barro_abundances_X.csv
        pd.DataFrame(
            [(name, count) for name, count in 
             sorted(species_counts.items(), 
                    key=lambda x: x[1], 
                    reverse=True)],
            columns=['name', 'count']
        ).to_csv(
            os.path.join(output_folder, f'barro_abundances_{file_num}.csv'), 
            index=False
        )
        
        # 3. Save names_barro_X.txt (sorted by abundance)
        with open(os.path.join(output_folder, f'names_barro_{file_num}.txt'), 'w') as f:
            f.write("\n".join(sorted_species))
        
        print(f"Processed {dat_file}:")
        print(f"  • {len(barro_data)} trees in 500x500m area")
        print(f"  • {len(species_counts)} species total")

if __name__ == '__main__':
    process_dat_files(max_x=500, max_y=500)
