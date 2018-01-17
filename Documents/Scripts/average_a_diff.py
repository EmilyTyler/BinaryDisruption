#Plot average change in a over time for multiple systems along with YCG prediction for binning and Monte Carlo

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from monte_carlo import MCEncounters
from encounters import binning
from random_binary import setupRandomBinary
import matplotlib.pyplot as plt
from scipy.constants import G, year, mega, giga

#Initial semi-major axis
a_0 = 0.1 * 3.086*10.0**16.0
#Initial eccentricity
e_0 = 0.7
#Mass of binary stars
m1 = 1.0*10.0**30.0
m2 = 1.0*10.0**30.0
#Mass of perturbers
M_p = 3.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 220.0 * 1000.0
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/((3.086*10.0**16.0)**3.0)
#Number density of perturbers
n_p = rho/M_p
#Time to run simulation for
t_end = 10.0 * giga * year
#Number of systems to run
# TAKES 10 MINUTES FOR 10Gyr and 100 systems
N_bin = 1000

#Binning
t_B, A_B_new, es_B_new, N_broken_new = binning(v_rms, n_p, t_end, a_0, e_0, m1, m2, M_p)
N_t = np.size(t_B)
#
A_B = np.zeros((N_bin, N_t), dtype=float)
A_B[:,0] = np.array([a_0]*N_bin)
A_B[1] = A_B_new
#
es_B = np.zeros((N_bin, N_t), dtype=float)
es_B[:,0] = np.array([e_0]*N_bin)
es_B[1] = es_B_new
#
N_broken_B_ind = np.zeros((N_bin, N_t), dtype=int)
N_broken_B_ind[:,0] = np.array([0]*N_bin)
N_broken_B_ind[1] = N_broken_new

for i in range(2, N_bin):
        t_B, A_B[i], es_B[i], N_broken_B_ind[i] = binning(v_rms, n_p, t_end, a_0, e_0, m1, m2, M_p)        
#Number of binaries broken
N_broken_B = np.sum(N_broken_B_ind, axis=0)
#Make cumulative
for i in range(1, N_t):
        N_broken_B[i] += N_broken_B[i-1]
#Average
A_B_avg = np.zeros(N_t, dtype=float)
for i in range(N_t):
        A_B_avg[i] = np.sum(A_B[:,i])/np.size(np.nonzero(A_B[:,i]))

#Monte Carlo
#Number of timesteps
N_t = 50
#Time array
t_MC = np.linspace(0.0, t_end, N_t)
dt = t_MC[1] - t_MC[0]
#Semi-major axis array
A_MC = np.zeros((N_t, N_bin), dtype=float)
A_MC[0] = np.array([a_0]*N_bin)
#Eccentricity array
es_MC = np.zeros((N_t, N_bin), dtype=float)
es_MC[0] = np.array([e_0]*N_bin)
N_broken_MC = np.zeros(N_t, dtype=int)
for i in range(1, N_t):
        A_MC[i], es_MC[i], N_broken_MC[i] = MCEncounters(v_rms, n_p, dt, m1, m2, M_p, A_MC[i-1], es_MC[i-1], N_bin)
#Average
A_MC_avg = np.zeros(N_t, dtype=float)
for i in range(N_t):
        A_MC_avg[i] = np.sum(A_MC[i])/np.size(np.nonzero(A_MC[i]))
                

#Plots

#Plot binning and average binning
plt.title('Binning')
for i in range(N_bin):
        plt.plot(t_B/(10.0**6.0*365.25*24.0*60.0*60.0), A_B[i]/(3.086*10.0**16.0), color='darkgrey')
plt.plot(t_B/(10.0**6.0*365.25*24.0*60.0*60.0), A_B_avg/(3.086*10.0**16.0), color='forestgreen')
plt.xlabel('Time/Myr')
plt.ylabel('Semimajor axis/pc')
plt.show()

#Plot MC and average MC
plt.title('Monte Carlo')
for i in range(N_bin):
        plt.plot(t_MC/(10.0**6.0*365.25*24.0*60.0*60.0), A_MC[:,i]/(3.086*10.0**16.0), color='darkgrey')
plt.plot(t_MC/(10.0**6.0*365.25*24.0*60.0*60.0), A_MC_avg/(3.086*10.0**16.0), color='forestgreen')
plt.xlabel('Time/Myr')
plt.ylabel('Semimajor axis/pc')
plt.show()

#Plot YCG and averages
b_min = (np.pi*n_p*v_rms*(10.0**10.0*365.25*24.0*60.0*60.0))**(-0.5)
A_YCG = np.zeros(N_t)
A_YCG[0] = a_0
for i in range(1,N_t):
        X = setupRandomBinary(A_YCG[i-1], e_0, m1, m2)
        V = np.linalg.norm(X[2]-X[3])
        delta_v = np.sqrt(16.0*np.pi*G**2.0*rho*M_p*(t_MC[i]-t_MC[i-1])*np.log(A_YCG[i-1]/b_min)/v_rms)
        delta_a = 2.0*A_YCG[i-1]**2.0*V*delta_v/(G*m1)
        A_YCG[i] = A_YCG[i-1] + delta_a

#Original YCG calc        
delta_v = np.sqrt(16.0*np.pi*G**2.0*rho*M_p*t_MC*np.log(a_0/b_min)/v_rms)
X = setupRandomBinary(a_0, e_0, m1, m2)
R = np.linalg.norm(X[0]-X[1])
V = np.linalg.norm(X[2]-X[3])
delta_a = 2.0*a_0**2.0*V*delta_v/(G*m1)
delta_a += a_0
A_YCG2 = delta_a

plt.plot(t_MC/(10.0**6.0*365.25*24.0*60.0*60.0), A_YCG/(3.086*10.0**16.0), label='YCG')
plt.plot(t_MC/(10.0**6.0*365.25*24.0*60.0*60.0), A_YCG2/(3.086*10.0**16.0), label='Original YCG')
plt.plot(t_B/(10.0**6.0*365.25*24.0*60.0*60.0), A_B_avg/(3.086*10.0**16.0), label='Binning')
plt.plot(t_MC/(10.0**6.0*365.25*24.0*60.0*60.0), A_MC_avg/(3.086*10.0**16.0), label='MC')
plt.xlabel('Time/Myr')
plt.ylabel('Semimajor axis/pc')
plt.ylim(np.min([np.min(A_MC_avg), np.min(A_B_avg)])/(3.086*10.0**16.0), np.max([np.max(A_MC_avg), np.max(A_B_avg)])/(3.086*10.0**16.0))
plt.legend()
plt.show()

#Plot number of binaries broken against time
plt.plot(t_B/(mega*year), N_broken_B, label='Binning')
plt.plot(t_MC/(mega*year), N_broken_MC, label='Monte Carlo')
plt.xlabel('Time/Myr')
plt.ylabel('Number of binaries broken')
plt.legend()
plt.show()



























