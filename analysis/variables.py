
#senm parameters

nx = 500
ny = 500
kernel = '5'
nu = '9.000000e-05'
num_realizations = 100

chosen_kernels= ['2','3','5']

#other senm variables

senm_spatial_file_template = 'AConf_NX{nx}_NY{ny}_NU{nu}_Square_Kernsize{kernel}_rel.{realization}'
senm_abundance_file_template = 'Abund_NX{nx}_NY{ny}_NU{nu}_Square_Kernsize{kernel}_rel.{realization}'


simulations_path = '../senm/'

#real data 

forests = ['barro', 'wanang']

path_template = '../data/{forest}/'
census_template = '{forest}_{census}.csv'
abundance_template = '{forest}_abundances_{census}.csv'
names_template = 'names_{forest}_{census}.txt'
