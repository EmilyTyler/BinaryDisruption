#To evolve a distribution of binaries using MC method

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from monte_carlo import MCEncounters
from random_power import randomPower
import matplotlib.pyplot as plt
from scipy.constants import au
from frequency import calcFrequency

#Mass of binary stars
m1 = 0.5 * 2.0*10.0**30.0
m2 = 0.5 * 2.0*10.0**30.0
#Mass of perturbers
M_p = 3.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 100.0 * 1000.0
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/((3.086*10.0**16.0)**3.0)
#Number density of perturbers
n_p = rho/M_p
#Time to run simulation for
T = 10.0**10.0*365.25*24.0*60.0*60.0

#Number of binary pairs
N_bin = 10000

#Semi-major axis array
#Draw a from log distribution
#x = log(a/AU)
x_min = 1.0
x_max = 5.5
x = np.random.random(N_bin) * (x_max - x_min) + x_min
a = 10.0**x * au

a_bins, N_a = calcFrequency(a, 100)
plt.loglog(a_bins, N_a)
plt.show()

#Eccentricity array
#Draw e from distribution uniform in e^2 between 0 and 1
e = (np.random.random(N_bin))**(1.0/3.0)

e_bins, N_e = calcFrequency(e, 100)
plt.plot(e_bins, N_e)
plt.show()

































