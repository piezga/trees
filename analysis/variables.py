
#senm parameters

nx = 512
ny = 512
kernel = '10'
nu = '9.000000e-05'
num_realizations = 50

chosen_kernels= ['1','2','5','10']

#other senm variables

spatial_file_template = 'AConf_NX{nx}_NY{ny}_NU{nu}_Square_Kernsize{kernel}_rel.'
abundance_file_template = 'Abund_NX{nx}_NY{ny}_NU{nu}_Square_Kernsize{kernel}_rel.'


simulations_path = '../senm/'

#real data 

forests = ['barro', 'wanang']

barro_path = '../data/barro/'
wanang_path = '../data/wanang/'

barro_spatial_file = barro_path + "bcitree8.csv"
barro_names_file = barro_path + "names_bci_sorted.txt"