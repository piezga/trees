
#senm parameters

nx = 512
ny = 512
kernel = '5'
nu = '9.000000e-05'
num_realizations = 100

chosen_kernels= ['1','2','5','10']

#other senm variables

spatial_file_template = 'AConf_NX{nx}_NY{ny}_NU{nu}_Square_Kernsize{kernel}_rel.{realization}'
abundance_file_template = 'Abund_NX{nx}_NY{ny}_NU{nu}_Square_Kernsize{kernel}_rel.'


simulations_path = '../senm/'

#real data 

forests = ['barro', 'wanang']

path_template = '../data/{forest}/'
census_template = '{forest}_{census}.csv'
names_template = 'names_{forest}_{census}.txt'
