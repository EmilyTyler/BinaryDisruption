#Plot the fraction of binaries that break over a two dimensional parameter space of a and M_p

import numpy as np
import os
os.system("python setup.py build_ext --inplace")
import matplotlib.pyplot as plt

from encounters import impulseEncounter
from scipy.constants import parsec, au, giga, year

#Initialise variables
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 1.0*10.0**30.0
m2 = 1.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 220.0 * 1000.0
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/(parsec**3.0)


#Total number of encounters per mass and semi-major axis bin
N_enc = 1000

#Perturber mass
#Minimum mass
M_min = 2.0*10.0**30.0
M_max = 100000000.0 * 2.0*10.0**30.0
#Number of masses to test
N_M = 100
#Set up logarithmic mass bins
dlogM = (np.log(M_max)-np.log(M_min))/N_M
M_bins = np.array([M_min*np.exp(dlogM*i) for i in range(N_M)])

#Semi-major axis
#Minimum a
a_min = 1000.0 * au
a_max = 100000000.0 * au
#Number of a's to test
N_a = 100
#Set up logarithmic a bins
dloga = (np.log(a_max)-np.log(a_min))/N_a
a_bins = np.array([a_min*np.exp(dloga*i) for i in range(N_a)])

#Fraction of encounters that break binary
F_enc = np.zeros((N_M, N_a))
#Number density of perturbers
n_p = rho/M_bins
#Minimum impact parameter
b_min = (np.pi*n_p*v_rms*(10.0*giga*year))**(-0.5)

for k in range(N_a):
        for j in range(N_M):      
                for i in range(N_enc):
                        (notBound, a_new, e_new) = impulseEncounter(m1, m2, v_rms, b_min[j], a_bins[k], e, M_bins[j])
                        if notBound:
                                F_enc[j,k] += 1

#Normalise F_enc
F_enc /= N_enc

#Contour plot
plt.title(r'Fraction of binaries broken due to a single encounter at $b=b_\mathrm{min}$')
plt.contourf(M_bins/(2.0*10.0**30.0), a_bins/au, np.transpose(F_enc), levels=np.linspace(0.0, 1.0, 11))
plt.xlabel('Perturber mass, solar masses')
plt.ylabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.colorbar();
plt.show()
