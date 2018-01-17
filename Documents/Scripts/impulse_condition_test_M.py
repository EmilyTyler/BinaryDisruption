#Test for my impulse condition

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from matplotlib import pyplot as plt
from scipy.constants import au
from impulse_test_encounter import impulseTestEncounter, bImpulseValid, MImpulseValid

#Initialise variables
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 2.0*10.0**30.0
m2 = 2.0*10.0**30.0
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/((3.086*10.0**16.0)**3.0)
#RMS of Maxwellian velocity distribution, m/s
v_rms = 220.0 * 1000.0


#Semi-major axes
a_min = 10.0**3.0 * au
a_max = 10.0**12.0 * au
#Number of a's to test
N_a = 10
#a array
dloga = (np.log(a_max)-np.log(a_min))/(N_a-1.0)
a_bins = np.array([a_min*np.exp(dloga*i) for i in range(N_a)])

#Number of encounters for each a value
N_enc = 10

#Average fractional error
a_frac_avg = np.zeros(N_a, dtype=float)
for i in range(N_a):
        M_p = MImpulseValid(a_bins[i])
        n_p = rho/M_p
        b = bImpulseValid(a_bins[i], n_p, v_rms, M_p, m1, m2)
        for j in range(N_enc):
                (notBound_imp, a_imp, e_imp, a_frac, e_diff) = impulseTestEncounter(m1, m2, v_rms, b, a_bins[i], e, M_p) 
                a_frac_avg[i] += a_frac
                    
#Normalise
a_frac_avg /= N_enc
print(np.absolute(a_frac_avg))
#Plot
plt.loglog(a_bins/au, np.absolute(a_frac_avg))
plt.ylabel('Average fractional error in semi-major axis')
plt.xlabel('Semi-major axis, au')
plt.show()