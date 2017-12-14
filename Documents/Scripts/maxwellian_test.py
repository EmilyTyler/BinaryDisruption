#Generate maxwellian distribution

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from matplotlib import pyplot as plt
from monte_carlo import draw_maxwellian
from monte_carlo import maxwellianPdf
from monte_carlo import maxwellianComparison

v_rms = 10.0**5.0
v_max = 10.0**1.0 * v_rms
v_min = 10.0**(-1.0) * v_rms
N = 10**6

max_dist = draw_maxwellian(v_rms, v_min, v_max, N)

#Create velocity bins
N_v = 1000
#dlogv = (np.log(v_max) - np.log(v_min))/N_v
#v_bins = np.array([v_min*np.exp(dlogv*i) for i in range(N_v)])
dv = (v_max - v_min)/N_v
v_bins = np.array([v_min + i*dv for i in range(N_v)])

#Number of points in each bin
p = np.zeros(N_v, dtype = int)
for x in max_dist:
        #i = int(np.floor(np.log(x/v_min)/dlogv))
        i = int(np.floor((x-v_min)/dv))
        p[i] += 1

#Normalise
p = p /(np.max(p)) * 4.0*(2.0*np.pi)**(-0.5)/v_rms*np.exp(-1.0)


plt.plot(v_bins/1000.0, p)
pdf = np.array([maxwellianPdf(v, v_rms) for v in v_bins])
plt.plot(v_bins/1000.0, pdf)
#comparison = [maxwellianComparison(v, v_rms, v_min, v_max) for v in v_bins]
#plt.plot(v_bins/1000.0, comparison)
plt.xlabel('v, km/s')
plt.ylabel('pdf')
plt.show()