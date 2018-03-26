import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from encounters import calc_b_max, integrateEncounter, impulseEncounter
from matplotlib import pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
from scipy.constants import giga, year, au, parsec

#Initialise variables
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 2.0*10.0**30.0
m2 = 2.0*10.0**30.0
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/(parsec**3.0)
#RMS of Maxwellian velocity distribution, m/s
v_rms = 220.0 * 1000.0

#Semi-major axes
a_min = 10.0**3.0 * au
a_max = 10.0**6.0 * au
#Number of a's to test
N_a = 20
#PBH masses to test
M_p_min = 10.0**1.0 * 2.0*10.0**30.0
M_p_max = 10.0**3.0 * 2.0*10.0**30.0
#Number of M_p's to test
N_M = 20

#Number of encounters per each pair of values
N_enc = 1000

#Set up logarithmic a bins
dloga = (np.log(a_max)-np.log(a_min))/N_a
a_bins = np.array([a_min*np.exp(dloga*i) for i in range(N_a)])
#Set up logarithmic M_p bins
dlogM = (np.log(M_p_max)-np.log(M_p_min))/N_M
M_p_bins = np.array([M_p_min*np.exp(dlogM*i) for i in range(N_M)])
#Average fractional difference in a
a_frac_avg = np.zeros((N_a, N_M), dtype=float)

#b=parsec
for j in range(N_M):
        n_p = rho/M_p_bins[j]
        for i in range(N_a):
                b = calc_b_max(M_p_bins[j], v_rms, a_bins[i], m1, m2)
                for k in range(N_enc):
                        (notBound, a_new, e_new) = impulseEncounter(m1, m2, v_rms, b, a_bins[i], e, M_p_bins[j])
                        a_frac_avg[i,j] += (a_new-a_bins[i])/a_bins[i]
#Normalise a_frac_avg
a_frac_avg /= N_enc

#Plot
plt.title(r'Absolute average fractional change in semi-major axis due to single encounter at $b=b_{\mathrm{max}}$')
ax = plt.gca()
cs = ax.contourf(a_bins/au, M_p_bins/(2.0*10.0**30.0), np.transpose(np.absolute(a_frac_avg)), locator=ticker.LogLocator())
plt.colorbar(cs)
plt.ylabel(r'Perturber mass, $M_\odot$')
plt.xlabel('Initial semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()
        
        
        
        
        
        
        
        
        
        