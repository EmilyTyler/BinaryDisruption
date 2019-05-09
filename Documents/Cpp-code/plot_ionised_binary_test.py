import matplotlib.pyplot as plt
import csv
import numpy as np
from scipy.constants import au
plt.rc('font', family='serif')

N_bin = 10**5
M_p = np.array([10.0**0.0, 10.0**1.0, 10.0**2.0, 10.0**3.0, 10.0**4.0])*2.0*10.0**30.0
a = np.array([10.0**3.0, 10.0**4.0, 10.0**5.0, 10.0**6.0]) * au

#Number of binaries bound that were once unbound
number_rebound = np.array([	[0, 0, ],
							[],
							[],
							[],
							[]])

#Number of binaries unbound and within 100pc
number_close = np.array([	[33, 342, ],
							[],
							[],
							[],
							[]])