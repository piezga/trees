import numpy as np
from variables import *
import os


# Move files to data folder
[os.system(f"mv ../{file}* {simulations_path}/ 2>/dev/null") 
    for file in ["AConf","Abund"]]

