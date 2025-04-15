import pandas as pd
import numpy as np
from variables import *
import matplotlib.pyplot as plt
from functions import load_file_with_padding

censi = ['1', '4']
wanang = pd.DataFrame()
senm = pd.DataFrame()


# Wanang data

for census in censi:
    abundance_file = wanang_path + f"wanang_abundances_{census}.csv"
    df = pd.read_csv(abundance_file)
    
    # Calculate relative abundance
    rel_values = df['count'] / df['count'].sum()
    
    # Create temporary DataFrame with species as index
    temp_df = pd.DataFrame({
        f'census {census}': rel_values})
    
    # Merge with main DataFrame
    if wanang.empty:
        wanang = temp_df
    else:
        wanang = wanang.join(temp_df, how='outer')

# Fill NA values with 0 (for species missing in one census)
wanang = wanang.fillna(0)

# Barro data
"""
abundance_file = barro_path + f"barro_abundances.csv"
df = pd.read_csv(abundance_file)
rel_values = df['count'] / df['count'].sum()
barro = pd.DataFrame(rel_values)
"""

for kernel in chosen_kernels:


    for realization in range(1,num_realizations+1):



        current_file = abundance_file_template.format(
            nx=nx,  
            ny=ny,
            nu=nu,
            kernel=kernel,
            realization=realization
        )
        
        senm_abundances = np.full(len(df['count']),0.1)
        abundances = load_file_with_padding(simulations_path + current_file,
                                            len(rel_values),2)[:,1]
        abundances = abundances[abundances.argsort()[::-1]]  # Sort
        senm_abundances += abundances
                
        mean_senm_abundances = senm_abundances/num_realizations
        relative_mean_senm_abundances = mean_senm_abundances/mean_senm_abundances.sum()
                
        senm[f'k = {kernel}'] = relative_mean_senm_abundances.tolist()




# Create figure and axis
fig, ax = plt.subplots(figsize=(12, 7))

# Plot empirical data (wanang)
wanang.plot(logy=True, style='o', 
                  markersize=5, 
                  color=['blue', 'green'],  # Different colors for each census
                  ax=ax, 
                  label=['Wanang Census 1', 'Wanang Census 4'])

# Barro
"""
barro.plot(logy=True, style='o', 
                  markersize=5, 
                  color=['black'],  
                  ax=ax, 
                  label=['Barro Census 1'])
"""




# Plot model results (senm)
senm.plot(logy=True, style='.', 
         markersize=2, 
         color=['red', 'purple', 'orange'],
         ax=ax,
         label=[f'SENM (k={k})' for k in chosen_kernels])

# Customize plot
ax.set_ylabel('Relative abundance (log scale)')
ax.set_xlabel('Species rank')
ax.set_title('Empirical vs Modeled Species Abundance Distributions')
ax.legend(fontsize = '20')

plt.tight_layout()
plt.show()