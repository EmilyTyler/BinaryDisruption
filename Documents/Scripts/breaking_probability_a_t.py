#Fraction of binaries that break when a is changed
import numpy as np
import matplotlib.pyplot as plt

from encounters import encounter

#Initialise variables
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 1.0*10.0**30.0
m2 = 1.0*10.0**30.0
#Perturber mass
M_p = 10000.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 100.0 * 1000.0
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/((3.086*10.0**16.0)**3.0)
#Number density of perturbers
n_p = rho/M_p
#Minimum impact parameter
b_min = (np.pi*n_p*v_rms*(10.0**10.0*365.25*24.0*60.0*60.0))**(-0.5)

#Total number of encounters per semi-major axis
N_enc = 1000
#Minimum a
a_min = 1000.0 * 1.5*10.0**11.0
a_max = 1000000.0 * 1.5*10.0**11.0
#Number of masses to test
N_a = 100
#Set up logarithmic mass bins
dloga = (np.log(a_max)-np.log(a_min))/N_a
a_bins = np.array([a_min*np.exp(dloga*i) for i in range(N_a)])

v = v_rms
b = b_min

#Fraction of encounters that break binary
F_enc = np.zeros(N_a)
for j in range(N_a):
        for i in range(N_enc):
                (notBound, a_new, e_new) = encounter(m1, m2, v, b, a_bins[j], e, M_p)
                if notBound:
                        F_enc[j] += 1

#Normalise F_enc
F_enc /= N_enc

#Plot
plt.semilogx(a_bins/(1.5*10.0**11.0), F_enc)
plt.xlabel('Semi-major axis, au')
plt.ylabel('Fraction of binaries broken')
plt.title(r'Fraction of binaries broken due to single encounter at $b=b_\mathrm{min}$, $M_p = 10000M_\odot$')
plt.show()