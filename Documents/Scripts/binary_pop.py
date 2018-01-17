#To evolve a distribution of binaries using MC method

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from monte_carlo import MCEncounters
import matplotlib.pyplot as plt
from scipy.constants import au, year, mega, giga
from frequency import calcFrequency

#Mass of binary stars
m1 = 0.5 * 2.0*10.0**30.0
m2 = 0.5 * 2.0*10.0**30.0
#Mass of perturbers
M_p = 10.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 220.0 * 1000.0
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/((3.086*10.0**16.0)**3.0)
#Number density of perturbers
n_p = rho/M_p
#Time to run simulation for
T = 10.0 * giga*year

#Number of binary pairs
#TAKES 1.5 HOURS TO RUN 1000
N_bin = 4000

print('Generating initial population')
#Semi-major axis array
#Draw a from distribution dN/dloga\propto a^{1-alpha}
alpha = 2.0
a_min = 10.0**1.0
a_max = 10.0**5.5
if alpha == 2.0:
        c = np.log(a_min)/np.log(a_max/a_min)
        a = (a_max/a_min)**(np.random.random(N_bin) + c) * au
else:
        a = (np.random.random(N_bin)*(a_max**(2.0-alpha) - a_min**(2.0-alpha)) + a_min**(2.0-alpha))**(1.0/(2.0-alpha)) * au

#a_bins, N_a = calcFrequency(a/au, 1000)
#plt.loglog(a_bins, N_a)
#plt.xlim(10.0**3.0, 3.0*10.0**5.0)
#plt.ylim(5.0*10.0**(-1.0), 3.0*10.0**3.0)
#plt.show()

#Eccentricity array
#Draw e from distribution uniform in e^2 between 0 and 1
e = (np.random.random(N_bin))**(1.0/3.0)

#e_bins, N_e = calcFrequency(e, 1000)
#plt.plot(e_bins, N_e)
#plt.show()

#Evolve distribution in time
print('Evolving population')
a_end, e_end = MCEncounters(v_rms, n_p, T, m1, m2, M_p, a, e, N_bin)
a_end = a_end[np.nonzero(a_end)]
e_end = e_end[np.nonzero(a_end)]


print('Plotting')
#Plot final and initial distributions
#Number of bins
N_bins = 10
a_bins_old, N_a_old = calcFrequency(a, N_bins, log=True)
a_bins_new, N_a_new = calcFrequency(a_end, N_bins, log=True)
e_bins_old, N_e_old = calcFrequency(e, N_bins)
e_bins_new, N_e_new = calcFrequency(e_end, N_bins)
#Semi-major axis
plt.loglog(a_bins_old/au, N_a_old)
plt.loglog(a_bins_new/au, N_a_new)
#plt.xlim(10.0**3.0, 3.0*10.0**5.0)
#plt.ylim(5.0*10.0**(-1.0), 3.0*10.0**3.0)
plt.xlabel('Semi-major axis, au')
plt.ylabel('Number of binaries')
plt.show()
#Eccentricity
plt.plot(e_bins_old, N_e_old)
plt.plot(e_bins_new, N_e_new)
plt.xlim(0.0, 1.0)
plt.xlabel('Eccentricity')
plt.ylabel('Number of binaries')
plt.show()

print('Finished')



























