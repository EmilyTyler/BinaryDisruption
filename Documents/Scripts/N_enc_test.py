#To test the distributions of number of encounters

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

N_v = 100
N_b = 100
#Number of numbers of encounters to generate
N_N_enc = 1000000

#Monte Carlo
#Number of encounters
MC_mean = T*encounterRate(n_p, v_rms, b_min, b_max, v_min, v_max) 
N_enc_MC = np.random.poisson(MC_mean, size=N_N_enc)
print('MC mean = ', MC_mean)

#Binning
#b bins for encounter rate
dlogb = (np.log(b_max)-np.log(b_min))/N_b
b = np.array([b_min*np.exp(i*dlogb) for i in range(N_b)])
db_B = b * (np.exp(dlogb) - 1.0)

#db_B = (b_max - b_min)/(N_b)
#b = np.array([b_min + i*db_B for i in range(N_b)])

#v bins for encounter rate
dlogv = (np.log(v_max)-np.log(v_min))/N_v
v = np.array([v_min*np.exp(i*dlogv) for i in range(N_v)])
dv_B = v * (np.exp(dlogv) - 1.0)

#dv_B = (v_max - v_min)/(N_v)
#v = np.array([v_min + i*dv_B for i in range(N_v)])

#R[i,j] is the encounter rate for objects with impact parameter b[i] and relative velocity v[j]
R = np.zeros([N_b,N_v], dtype=float)
for i in range(N_b):
        for j in range(N_v):
                R[i,j] = encounterRate(n_p, v_rms, b[i], b[i]+db_B[i], v[j], v[j]+dv_B[j])
#Time step
dt = 1.0/np.amax(R)
#Number of timesteps
N_t = int(np.ceil(T/dt + 1))
#Sum of mean values
B_mean = np.sum(R*T)
print('Binning mean = ', B_mean)
N_enc_B = np.zeros(N_N_enc)
for k in range(N_N_enc):
        #Number of encounters matrix
        N = np.rollaxis(np.array([[np.random.poisson(R[i,j]*dt, size=N_t) for j in range(N_v)] for i in range(N_b)]), 2)
        #Total number of encounters
        N_enc_B[k] = np.sum(N)
        
#Bin N_enc values
N_enc_bins_MC, N_N_enc_MC, dN_enc_MC = calcFrequency(N_enc_MC, 200)
N_enc_bins_B, N_N_enc_B, dN_enc_B = calcFrequency(N_enc_B, 200)


#Plot N_enc distribution
plt.plot(N_enc_bins_MC, N_N_enc_MC/N_N_enc/dN_enc_MC, label='Monte Carlo')
plt.plot(N_enc_bins_B, N_N_enc_B/N_N_enc/dN_enc_B, label='Binning')
plt.plot([MC_mean]*np.size(N_N_enc_MC), N_N_enc_MC/N_N_enc/dN_enc_MC, label='MC mean')
plt.plot([B_mean]*np.size(N_N_enc_B), N_N_enc_B/N_N_enc/dN_enc_B, label='Binning mean')
plt.xlabel('Number of encounters')
plt.ylabel('Probability density of number of encounters')
plt.legend()
plt.show()

