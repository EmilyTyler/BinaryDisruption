import csv
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.constants import au, parsec, giga, year
plt.rc('font', family='serif')

average_dm_density = np.array([0.0113, 0.00903, 0.00682, 0.000763, 0.0131, 0.000917])
r_max = np.array([9.13, 10.3, 25.5, 186, 8.41, 66.0])
plt.xlabel(r'$r_\mathrm{max}$, kpc')
plt.ylabel(r'Time-averaged dark matter density, $M_\odot$pc$^{-3}$')
plt.scatter(r_max, average_dm_density)
plt.show()