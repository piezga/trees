import pandas as pd
import numpy as np
from variables import *
import matplotlib.pyplot as plt
from functions import load_file_with_padding

censi = ['1', '4']
wanang = pd.DataFrame()


# Wanang data

for census in censi:
    wanang_abundance_file = wanang_path + f"wanang_abundances_{census}.csv"
    df = pd.read_csv(wanang_abundance_file)
    
    # Calculate relative abundance
    rel_values = df['count'] / df['count'].sum()
    
    # Create temporary DataFrame with species as index
    temp_df = pd.DataFrame({
        f'Wanang census {census}': rel_values})
    
    # Merge with main DataFrame
    if wanang.empty:
        wanang = temp_df
    else:
        wanang = wanang.join(temp_df, how='outer')

# Fill NA values with 0 (for species missing in one census)
wanang = wanang.head(100)


# Barro data
barro_abundance_file = barro_path + f"barro_abundances.csv"
df = pd.read_csv(barro_abundance_file)
rel_values = df['count'] / df['count'].sum()
barro = pd.DataFrame({'Barro' : rel_values})
barro = barro.head(100)

# Create figure and axis
fig, ax = plt.subplots(figsize=(12, 7))

# Plot empirical data (wanang)
wanang.plot(logy=True, style='o', 
                  markersize=5, 
                  color=['blue', 'green'],  # Different colors for each census
                  ax=ax, 
                  label=['Wanang Census 1', 'Wanang Census 4'])

# Barro
barro.plot(logy=True, style='x', 
                  markersize=5, 
                  color=['black'],  
                  ax=ax, 
                  label=['Barro Census 1'])



# Customize plot
ax.set_ylabel('Relative abundance (log scale)')
ax.set_xlabel('Species rank')
ax.set_title('Comparison between Barro Colorado and wanang (First 100 species)')
ax.legend(fontsize = '20')

plt.tight_layout()
plt.show()

