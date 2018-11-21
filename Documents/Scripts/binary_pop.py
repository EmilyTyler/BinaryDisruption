#To evolve a distribution of binaries using MC method

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from monte_carlo import MCEncounters_new, JTEncounters, WSWEncounters, ClosestEncounters
import matplotlib.pyplot as plt
from scipy.constants import au, year, mega, giga, parsec, kilo
from frequency import calcFrequency
from internal_units import *

#Mass of binary stars
m1 = 1.0*10.0**30.0 / mass_scale()
m2 = 1.0*10.0**30.0 / mass_scale()
#Mass of perturbers
M_p = 10.0 * 2.0*10.0**30.0 / mass_scale()
#RMS of Maxwellian velocity distribution, m/s
v_rms = 220.0 * 1000.0 /length_scale() * time_scale()
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/(parsec**3.0)
#Internal
rho = rho /mass_scale() * length_scale()**3.0
#Number density of perturbers
n_p = rho/M_p
#Time to run simulation for
T = 10.0 * giga * year / time_scale()

#Number of binary pairs
#TAKES 25-40 minutes TO RUN 1000 for T=10Gyr
N_bin = 10**3

print('Initial number of binaries =', N_bin)

print('Generating initial population')
#Semi-major axis array
#Draw a from distribution dN/dloga\propto a^{1-alpha}
alpha = 1.0
a_min = 10.0**3.0*au / length_scale()
a_max = 10.0**6.0*au / length_scale()
if alpha == 2.0:
        c = np.log(a_min)/np.log(a_max/a_min)
        a_ini = (a_max/a_min)**(np.random.random(N_bin) + c)
else:
        a_ini = (np.random.random(N_bin)*(a_max**(2.0-alpha) - a_min**(2.0-alpha)) + a_min**(2.0-alpha))**(1.0/(2.0-alpha))
        
#a_ini = np.array([4.55*10.0**5.0*au]*N_bin)

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
a_fin, e_fin, N_broken = MCEncounters_new(v_rms, n_p, T, m1, m2, M_p, a_ini, e_ini, N_bin, a_T=1.0*parsec/length_scale(), prefactor=0.001)
a_fin = a_fin[np.where(a_fin!=-1.0)]
e_fin = e_fin[np.where(e_fin!=-1.0)]

#Save data
print('Saving data')
np.savez('simulation_test_{}Msol_{}e{}.npz'.format(int(M_p*mass_scale()/(2.0*10.0**30.0)), int(N_bin/10**int(np.floor(np.log10(N_bin)))), int(np.floor(np.log10(N_bin)))), a_fin=a_fin, e_fin=e_fin, N_broken=N_broken)

#Load data
#print('Loading data')
#data = np.load('binary_pop.npz')
#a_fin, e_fin, N_broken = data['a_fin'], data['e_fin'], data['N_broken']

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
plt.loglog(a_bins_ini*length_scale()/au, N_a_ini/(da_ini*length_scale()/au))
plt.loglog(a_bins_fin*length_scale()/au, N_a_fin/(da_fin*length_scale()/au))
#plt.xlim(10.0**3.0, 3.0*10.0**5.0)
#plt.ylim(5.0*10.0**(-1.0), 3.0*10.0**3.0)
plt.xlabel('Semi-major axis, au')
plt.ylabel(r'Probability density, au$^{-1}$')
plt.show()
#Eccentricity
plt.plot(e_bins_ini, N_e_ini/N_bin/e_ini_binwidth)
plt.plot(e_bins_fin, N_e_fin/(N_bin-N_broken)/e_fin_binwidth)
plt.xlim(0.0, 1.0)
plt.xlabel('Eccentricity')
plt.ylabel('Probability density')
plt.show()

print('Plotting cumulative')
N_bins = 1000
a_min = np.min([np.min(a_fin), np.min(a_ini)])
a_max = np.max([np.max(a_fin), np.max(a_ini)])
dloga = (np.log(a_max) - np.log(a_min))/N_bins
a_bins = np.array([a_min * np.exp(i*dloga) for i in range(N_bins)])
N_a_ini = np.zeros(N_bins)
N_a_fin = np.zeros(N_bins)
for i in range(1, N_bins):
        N_a_fin[i] = np.size(np.where(a_fin <= a_bins[i]))
        N_a_ini[i] = np.size(np.where(a_ini <= a_bins[i]))
#Normalise
N_a_ini /= N_bins
N_a_fin /= np.size(a_fin)
plt.semilogx(a_bins*length_scale()/au, N_a_ini, label='Initial')
plt.semilogx(a_bins*length_scale()/au, N_a_fin, label='Final')
plt.xlabel(r'Semimajor axis $a$, au')
plt.ylabel(r'Fraction of binaries with semimajor axis less than $a$')
plt.legend()
plt.show()

print('Finished')



























