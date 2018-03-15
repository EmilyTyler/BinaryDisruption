#To evolve a distribution of binaries using MC method

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from monte_carlo import MCEncounters
import matplotlib.pyplot as plt
from scipy.constants import au, year, mega, giga
from frequency import calcFrequency

#Mass of binary stars
m1 = 2.0*10.0**30.0
m2 = 2.0*10.0**30.0
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
T = 10.0 * giga * year

#Number of binary pairs
#TAKES 25-40 minutes TO RUN 1000 for T=10Gyr
N_bin = 10000

print('Initial number of binaries =', N_bin)

print('Generating initial population')
#Semi-major axis array
#Draw a from distribution dN/dloga\propto a^{1-alpha}
alpha = 2.0
a_min = 10.0**3.0*au
a_max = 10.0**6.0*au
if alpha == 2.0:
        c = np.log(a_min)/np.log(a_max/a_min)
        a_ini = (a_max/a_min)**(np.random.random(N_bin) + c)
else:
        a_ini = (np.random.random(N_bin)*(a_max**(2.0-alpha) - a_min**(2.0-alpha)) + a_min**(2.0-alpha))**(1.0/(2.0-alpha))

#a_bins, N_a = calcFrequency(a/au, 1000)
#plt.loglog(a_bins, N_a)
#plt.xlim(10.0**3.0, 3.0*10.0**5.0)
#plt.ylim(5.0*10.0**(-1.0), 3.0*10.0**3.0)
#plt.show()

#Eccentricity array
#Draw e from distribution uniform in e^2 between 0 and 1
e_ini = (np.random.random(N_bin))**(1.0/3.0)

#e_bins, N_e = calcFrequency(e, 1000)
#plt.plot(e_bins, N_e)
#plt.show()

#Evolve distribution in time
print('Evolving population')
a_fin, e_fin, N_broken = MCEncounters(v_rms, n_p, T, m1, m2, M_p, a_ini, e_ini, N_bin)
a_fin = a_fin[np.where(a_fin!=-1.0)]
e_fin = e_fin[np.where(e_fin!=-1.0)]

print('Number of binaries broken =', N_broken)

print('Plotting')
#Plot final and initial distributions
#Number of bins
N_bins = 50
a_bins_ini, N_a_ini, a_ini_binwidth = calcFrequency(a_ini, N_bins, log=True)
a_bins_fin, N_a_fin, a_fin_binwidth = calcFrequency(a_fin, N_bins, log=True)
#Actual bin widths
da_ini = a_bins_ini *(np.exp(a_ini_binwidth)-1.0)
da_fin = a_bins_fin *(np.exp(a_fin_binwidth)-1.0)

e_bins_ini, N_e_ini, e_ini_binwidth = calcFrequency(e_ini, N_bins)
e_bins_fin, N_e_fin, e_fin_binwidth = calcFrequency(e_fin, N_bins)
#Semi-major axis
plt.loglog(a_bins_ini/au, N_a_ini/(N_bin-N_broken)/(da_ini/au))
plt.loglog(a_bins_fin/au, N_a_fin/(N_bin-N_broken)/(da_fin/au))
#plt.xlim(10.0**3.0, 3.0*10.0**5.0)
#plt.ylim(5.0*10.0**(-1.0), 3.0*10.0**3.0)
plt.xlabel('Semi-major axis, au')
plt.ylabel(r'Probability density, au$^{-1}$')
plt.show()
#Eccentricity
plt.plot(e_bins_ini, N_e_ini/(N_bin-N_broken)/e_ini_binwidth)
plt.plot(e_bins_fin, N_e_fin/(N_bin-N_broken)/e_fin_binwidth)
plt.xlim(0.0, 1.0)
plt.xlabel('Eccentricity')
plt.ylabel('Probability density')
plt.show()

#Plot absolute numbers rather than densities
#Set up bins
a_min = np.min([np.min(a_ini), np.min(a_fin)])
a_max = np.max([np.max(a_ini), np.max(a_fin)])
dloga = (np.log(a_max)-np.log(a_min))/N_bins
a_bins = np.array([a_min*np.exp(i*dloga) for i in range(N_bins)])
# e
e_min = np.min([np.min(e_ini), np.min(e_fin)])
e_max = np.max([np.max(e_ini), np.max(e_fin)])
dloge = (np.log(e_max)-np.log(e_min))/N_bins
e_bins = np.array([e_min*np.exp(i*dloge) for i in range(N_bins)])
#Count frequency
N_a_ini = np.zeros(N_bins, dtype=int)
N_a_fin = np.zeros(N_bins, dtype=int)
for val in a_ini:
        if val == a_max:
                i = N_bins - 1
        else:
                i = int(np.floor(np.log(val/a_min)/dloga))
                N_a_ini[i] += 1
for val in a_fin:
        if val == a_max:
                i = N_bins - 1
        else:
                i = int(np.floor(np.log(val/a_min)/dloga))
                N_a_fin[i] += 1
for val in e_ini:
        if val == e_max:
                i = N_bins - 1
        else:
                i = int(np.floor(np.log(val/e_min)/dloge))
                N_e_ini[i] += 1
for val in e_fin:
        if val == e_max:
                i = N_bins - 1
        else:
                i = int(np.floor(np.log(val/e_min)/dloge))
                N_e_fin[i] += 1
plt.loglog(a_bins/au, N_a_ini, label='Initial')
plt.loglog(a_bins/au, N_a_fin, label='Final')
plt.legend()
plt.xlabel('Semi-major axis, au')
plt.ylabel('Number of binaries')
plt.show()
plt.loglog(e_bins, N_e_ini, label='Initial')
plt.loglog(e_bins, N_e_fin, label='Final')
plt.legend()
plt.xlabel('Eccentricity')
plt.ylabel('Number of binaries')
plt.show()

print('Finished')



























