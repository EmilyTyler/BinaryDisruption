#To test the impulse approximation

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from matplotlib import pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors

from impulse_test_encounter import encounterGrid_M

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
a_min = 10.0**3.0 * 1.5*10.0**11.0
a_max = 10.0**12.0 * 1.5*10.0**11.0
#Number of a's to test
N_a = 20
#PBH masses to test
M_p_min = 10.0**1.0 * 2.0*10.0**30.0
M_p_max = 10.0**3.0 * 2.0*10.0**30.0
#Number of M_p's to test
N_M = 10

#Number of encounters per each pair of values
#TAKES 5.5 HOURS TO RUN for 20, 20, 100
N_enc = 20

a_frac_avg, a_bins, M_p_bins = encounterGrid_M(m1, m2, v_rms, rho, e, M_p_min, M_p_max, a_min, a_max, N_M, N_a, N_enc)

#Log absolute value
plt.title('Absolute average fractional error in semi-major axis due to impulse approximation')
ax = plt.gca()
cs = ax.contourf(a_bins/(1.5*10.0**11.0), M_p_bins/(2.0*10.0**30.0), np.transpose(np.absolute(a_frac_avg)), locator=ticker.LogLocator())
plt.colorbar(cs)
plt.ylabel('Perturber mass, solar masses')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()