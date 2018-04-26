#Plot probability of binary breaking against perturber mass
import numpy as np
import matplotlib.pyplot as plt

from encounters import impulseEncounter
from scipy.constants import parsec, giga, year

#Initialise variables
#Semi-major axis, m
a = 0.1 * parsec
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 1.0*10.0**30.0
m2 = 1.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 100.0 * 1000.0
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/(parsec**3.0)


#Total number of encounters per mass
N_enc = 1000
#Minimum mass
M_min = 2.0*10.0**30.0
M_max = 1000000.0 * 2.0*10.0**30.0
#Number of masses to test
N_M = 100
#Set up logarithmic mass bins
dlogM = (np.log(M_max)-np.log(M_min))/N_M
M_bins = np.array([M_min*np.exp(dlogM*i) for i in range(N_M)])

v = v_rms
#Fraction of encounters that break binary
F_enc = np.zeros(N_M)
for j in range(N_M):
        #Number density of perturbers
        n_p = rho/M_bins[j]
        #Minimum impact parameter
        b_min = (np.pi*n_p*v_rms*(10.0*giga*year))**(-0.5)
        b = b_min
        #Implement encounters
        for i in range(N_enc):
                (notBound, a_new, e_new) = impulseEncounter(m1, m2, v, b, a, e, M_bins[j])
                if notBound:
                        F_enc[j] += 1

#Normalise F_enc
F_enc /= N_enc

#Plot
plt.semilogx(M_bins/(2.0*10.0**30.0), F_enc)
plt.xlabel('Perturber mass, solar masses')
plt.ylabel('Fraction of binaries broken')
plt.title(r'Fraction of binaries broken due to single encounter at $b=b_\mathrm{min}$, $a = 0.1$pc')
plt.show()
