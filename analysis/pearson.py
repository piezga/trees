import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

x,y = [12,5,2,65], [43,2,6,1]

r, _ = stats.pearsonr(x, y)

plt.scatter(x,y)