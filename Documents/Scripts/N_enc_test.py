#To test the distributions of number of encounters

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

import random
import matplotlib.pyplot as plt
from scipy.constants import giga, year, parsec, mega
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
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/(parsec**3.0)
#Simulation time
T = 10.0*giga*year
print('T =', T)
#Mass of perturbers
M_p = 1.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 220.0 * 1000.0
#Number density of perturbers
n_p = rho/M_p

#Minimum impact parameter
b_min = (np.pi*n_p*v_rms*T)**(-0.5)
#b_min = 0.0
#Maximum impact parameter
b_max = calc_b_max(M_p, v_rms, a, m1, m2)
#Minimum velocity
v_min = 10.0**(-2.0) * v_rms
#Maximum velocity
v_max = 10.0**2.0 * v_rms

#Number of numbers of encounters to generate
N_N_enc = 10000

#Monte Carlo
#Number of encounters
#print('b_min =', b_min)
#print('b_max =', b_max)
#print('v_min =', v_min)
#print('v_max =', v_max)
MC_mean = T*encounterRate(n_p, v_rms, b_min, b_max, v_min, v_max) 
N_enc_MC = np.random.poisson(MC_mean, size=N_N_enc)
print('MC mean = ', MC_mean)
#print('N_enc_MC =', N_enc_MC)        

#Filter method
b_max_filter = 10.0 * b_max
F_mean = T*encounterRate(n_p, v_rms, b_min, b_max_filter, v_min, v_max)
N_enc_F = np.zeros(N_N_enc)
for i in range(N_N_enc):
        n = np.random.poisson(F_mean)
        b_F = draw_b(b_max_filter, n)
        b_F = b_F[np.where(b_F<=b_max)]
        N_enc_F[i] = np.size(b_F)
        
#Mean time method
N_enc_T = np.zeros(N_N_enc)
for i in range(N_N_enc):
        #time passed
        t = 0.0
        while t <= T:
                rate = encounterRate(n_p, v_rms, b_min, b_max, v_min, v_max)
                t += np.random.exponential(1.0/rate)
                N_enc_T[i] += 1
#print('N_enc_T =', N_enc_T)        
                
               
#Bin N_enc values
#print('N_enc_B =', N_enc_B)
N_enc_bins_MC, N_N_enc_MC, dN_enc_MC = calcFrequency(N_enc_MC, 100)
N_enc_bins_T, N_N_enc_T, dN_enc_T = calcFrequency(N_enc_T, 100)
N_enc_bins_F, N_N_enc_F, dN_enc_F = calcFrequency(N_enc_F, 100)

#Plot N_enc distribution
plt.plot(N_enc_bins_MC, N_N_enc_MC/N_N_enc/dN_enc_MC, label='Monte Carlo')
plt.plot(N_enc_bins_T, N_N_enc_T/N_N_enc/dN_enc_T, label='Mean time')
plt.plot(N_enc_bins_F, N_N_enc_F/N_N_enc/dN_enc_F, label='Filter')
plt.plot([MC_mean]*np.size(N_N_enc_MC), N_N_enc_MC/N_N_enc/dN_enc_MC, label='MC mean')
plt.xlabel('Number of encounters')
plt.ylabel('Probability density of number of encounters')
plt.legend()
plt.show()






















