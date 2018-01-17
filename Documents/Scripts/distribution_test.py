#Test b and v distributions and N_enc for binning and monte carlo methods
import numpy as np
import os
os.system("python setup.py build_ext --inplace")

import matplotlib.pyplot as plt
from scipy.constants import giga, year, parsec
from scipy.stats import maxwell
from monte_carlo import draw_b, maxwellianPdf
from encounters import calc_b_max, encounterRate
from frequency import calcFrequency

#Initialise variables
#Semi-major axis, m
a = 0.1 * parsec
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 2.0*10.0**30.0
m2 = 2.0*10.0**30.0
#Density of dark matter halo solar masses/pc**3
rho = 0.008
#Convert to SI
rho = rho * 2.0*10.0**30.0/(parsec**3.0)
#Simulation time
T = 10.0*giga*year
#Mass of perturbers
M_p = 3.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 220.0 * 1000.0
#Number density of perturbers
n_p = rho/M_p

#Minimum impact parameter
b_min = (np.pi*n_p*v_rms*T)**(-0.5)
#Maximum impact parameter
b_max = calc_b_max(M_p, v_rms, a, m1, m2)
#Minimum velocity
v_min = 10.0**(-2.0) * v_rms
#Maximum velocity
v_max = 10.0**2.0 * v_rms

#Monte Carlo
#Number of encounters
N_enc_MC = int(T*encounterRate(n_p, v_rms, b_min, b_max, v_min, v_max))
#b values
b_MC = draw_b(b_max, N_enc_MC)
#v values
v_MC = maxwell.rvs(scale=v_rms, size=N_enc_MC)
#Bin b values
b_bins_MC, N_b_MC = calcFrequency(b_MC, 100)
#Bin v values
v_bins_MC, N_v_MC = calcFrequency(v_MC, 100)

#Binning
#b bins for encounter rate
N_b = 100
#dlogb = (np.log(b_max)-np.log(b_min))/N_b
#b = np.array([b_min*np.exp(i*dlogb) for i in range(N_b)])
db = (b_max - b_min)/(N_b-1)
b = np.array([b_min + i*db for i in range(N_b)])
#v bins for encounter rate
N_v = 100
#dlogv = (np.log(v_max)-np.log(v_min))/N_v
#v = np.array([v_min*np.exp(i*dlogv) for i in range(N_v)])
dv = (v_max - v_min)/(N_v-1)
v = np.array([v_min + i*dv for i in range(N_v)])
#R[i,j] is the encounter rate for objects with impact parameter b[i] and relative velocity v[j]
R = np.zeros([N_b,N_v], dtype=float)
for i in range(N_b):
        for j in range(N_v):
                #R[i,j] = encounterRate(n_p, v_rms, b[i], b[i]*np.exp(dlogb), v[j], v[j]*np.exp(dlogv))
                R[i,j] = encounterRate(n_p, v_rms, b[i], b[i]+db, v[j], v[j]+dv)
#Time step
dt = 1.0/np.amax(R)
#Number of timesteps
N_t = int(np.ceil(T/dt + 1))
#Sum of mean values
print('Sum of mean values = ', np.sum(R*dt))
#Number of encounters matrix
N = np.rollaxis(np.array([[np.random.poisson(R[i,j]*dt, size=N_t) for j in range(N_v)] for i in range(N_b)]), 2)
#Total number of encounters
N_enc_B = np.sum(N)
#Array of indices where encounters happen
i_enc = np.transpose(np.array(np.nonzero(N)))
'''
#b values
b_B = np.array([])
#v values
v_B = np.array([])
for (i,j,k) in i_enc:
        b_B = np.append(b_B, b[j])
        v_B = np.append(v_B, v[k])
#Re-bin b and v in linear bins
b_bins_B, N_b_B = calcFrequency(b_B, 100)
v_bins_B, N_v_B = calcFrequency(v_B, 100)
'''
#b distribution
N_b_B = np.zeros(N_b, dtype=int)
#v distribution
N_v_B = np.zeros(N_v, dtype=int)
for (i,j,k) in i_enc:
        N_b_B[j] += N[i,j,k]
        N_v_B[k] += N[i,j,k]
b_bins_B = b
v_bins_B = v

#Print number of encounters
print('N_enc_MC = ', N_enc_MC)
print('N_enc_B = ', N_enc_B)

#Plot distributions
#Plot b distributions
plt.plot(b_bins_MC/parsec, N_b_MC/N_enc_MC, label='Monte Carlo')
plt.plot(b_bins_B/parsec, N_b_B/N_enc_B, label='Binning')
plt.xlabel('Impact parameter, pc')
plt.ylabel('Fraction of encounters at impact parameter b')
plt.legend()
plt.show()

#Plot v distributions
plt.plot(v_bins_MC/1000.0, N_v_MC/N_enc_MC, label='Monte Carlo')
plt.plot(v_bins_B/1000.0, N_v_B/N_enc_B, label='Binning')
pdf = np.array([maxwellianPdf(x, v_rms) for x in v_bins_MC])
plt.plot(v_bins_MC/1000.0, pdf*N_enc_MC, label='PDF MC')
pdf = np.array([maxwellianPdf(x, v_rms) for x in v_bins_B])
plt.plot(v_bins_B/1000.0, pdf*N_enc_B, label='PDF Binning')
plt.xlabel('v, km/s')
plt.ylabel('Fraction of encounters at relative velocity v')
plt.legend()
plt.show()















































