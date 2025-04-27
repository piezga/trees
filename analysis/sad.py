import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from variables import *
from functions import load_file_with_padding


# Load Barro data
barro_abundances_file = barro_path + 'barro_abundances.csv'
df = pd.DataFrame(pd.read_csv(barro_abundances_file))
df = df.drop(columns =  ['name'])
df['count'] = df['count']/df['count'].sum()


for kernel in chosen_kernels:

    senm_abundances = np.full(len(df.count),0.1)

    for realization in range(1,num_realizations+1):



        print(realization)
        current_file = abundance_file_template.format(
            nx=nx,  
            ny=ny,
            nu=nu,
            kernel=kernel,
            realization=realization
        )
        
        print(current_file)
        
        abundances = load_file_with_padding(simulations_path + current_file,
                                            len(df.count),2)[:,1]
        abundances = abundances[abundances.argsort()[::-1]]  # Sort
        senm_abundances += abundances
                
        mean_senm_abundances = senm_abundances/num_realizations
        relative_mean_senm_abundances = mean_senm_abundances/mean_senm_abundances.sum()
                
        df[f'k = {kernel}'] = relative_mean_senm_abundances.tolist()


df.plot(logy = True)
plt.legend(fontsize = '16')        







