import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from variables import *


forest = 'barro'
path = path_template.format(forest = forest, census = census)
censi = ['1','2','3','4','5','6','7']


for census in censi:
    #Load wanang data
    
    wanang_abundances_file = wanang_path + f"wanang{census}.txt"
    
    print(f"Analysing file {wanang_abundances_file}")
    
    # Dropping NaNs
    df = pd.DataFrame(pd.read_csv(wanang_abundances_file,
                                  low_memory=False)).dropna()
    
    #Filter indets
    
    filter_ = df['latin'] != 'indet indet'
    df = df[filter_]
    
    
    
    # Species count 
    
    # Create and sort the counts (highest to lowest)
    species_counts = (df['latin'].value_counts()
                      .reset_index()
                      .sort_values('count', ascending=False))
    
    # Rename columns
    species_counts.columns = ['species', 'count']
    
    # Show top 10 most abundant species
    print(species_counts.head(10))
    
    species_counts.to_csv(wanang_path + f"wanang_abundances_{census}.csv", 
                          index= False)


