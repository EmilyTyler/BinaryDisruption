#Plot the fraction of binaries that break over a two dimensional parameter space of a and b

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

import matplotlib.pyplot as plt
from encounters import impulseEncounter
from scipy.constants import parsec, au, giga, year, G

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
#Perturber mass
M_p = 10.0*2.0*10.0**30.0
#Number density of perturbers
n_p = rho/M_p

#Total number of encounters per a and b bin
N_enc = 100

#Semi-major axis
#Minimum a
a_min = 1000.0 * au
a_max = 1000000.0 * au
#Number of a's to test
N_a = 100
#Set up logarithmic a bins
dloga = (np.log(a_max)-np.log(a_min))/N_a
a_bins = np.array([a_min*np.exp(dloga*i) for i in range(N_a)])

#Impact parameter
#Minimum b
b_min = 0.01 * 1.5*(G*M_p**2.0*a_min**3.0/((m1+m2)*v_rms**2.0))**0.25
b_max = 100.0 * 1.5*(G*M_p**2.0*a_max**3.0/((m1+m2)*v_rms**2.0))**0.25
#Number of masses to test
N_b = 100
#Set up logarithmic mass bins
dlogb = (np.log(b_max)-np.log(b_min))/N_b
b_bins = np.array([b_min*np.exp(dlogb*i) for i in range(N_b)])


#Fraction of encounters that break binary
F_enc = np.zeros((N_b, N_a))


for k in range(N_a):
        for j in range(N_b):      
                for i in range(N_enc):
                        (notBound, a_new, e_new) = impulseEncounter(m1, m2, v_rms, b_bins[j], a_bins[k], e, M_p)
                        if notBound:
                                F_enc[j,k] += 1

#Normalise F_enc
F_enc /= N_enc

#Contour plot
plt.title(r'Fraction of binaries broken due to a single encounter with $M_p=10M_\odot$')
plt.contourf(b_bins/parsec, a_bins/parsec, np.transpose(F_enc), levels=np.linspace(0.0, 1.0, 11))
plt.xlabel('Impact parameter, pc')
plt.ylabel('Semi-major axis, pc')
plt.xscale('log')
plt.yscale('log')
plt.colorbar();
plt.show()