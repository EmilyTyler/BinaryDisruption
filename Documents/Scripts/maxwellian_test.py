#Generate maxwellian distribution

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from matplotlib import pyplot as plt
from monte_carlo import draw_maxwellian

v_rms = 10.0**5.0
v_max = 10.0**6.0
v_min = 10.0**(-3.0)
N = 10**5

max_dist = draw_maxwellian(v_rms, v_min, v_max, N)

#Create velocity bins
N_v = 1000
dlogv = (np.log(v_max) - np.log(v_min))/N_v
v_bins = np.array([v_min*np.exp(dlogv*i) for i in range(N_v)])

#Number of points in each bin
p = np.zeros(N_v, dtype = int)
for x in max_dist:
        i = int(np.floor(np.log(x/v_min)/dlogv))
        p[i] += 1
        
     
plt.plot(v_bins/1000.0, p)
plt.xlabel('v, km/s')
plt.ylabel('pdf')
plt.show()