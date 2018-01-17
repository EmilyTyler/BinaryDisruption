import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from monte_carlo import MCEncounters
import matplotlib.pyplot as plt
from scipy.constants import year, mega, giga, parsec

#Initial semi-major axis
a_0 = 0.1 * 3.086*10.0**16.0
#Initial eccentricity
e_0 = 0.7
#Mass of binary stars
m1 = 1.0*10.0**30.0
m2 = 1.0*10.0**30.0
#Mass of perturbers
M_p = 3.0 * 2.0*10.0**30.0
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/((3.086*10.0**16.0)**3.0)
#Number density of perturbers
n_p = rho/M_p
#Time to run simulation for
t_end = 10.0 * giga * year
#Number of systems to run
# TAKES 2h10m MINUTES FOR 10Gyr and 1000 systems and N_v=5
N_bin = 1000

#Velocity bins
N_v = 5
v_min = 200.0 * 1000.0
v_max = 240.0 * 1000.0
dv = (v_max - v_min)/(N_v-1)
v_bins = np.array([v_min + i*dv for i in range(N_v)])

#Number of timesteps
N_t = 100
#Time array
t = np.linspace(0.0, t_end, N_t)
dt = t[1] - t[0]
#Semi-major axis array
A = np.zeros((N_t, N_v, N_bin), dtype=float)
A[0] = np.array([a_0]*N_bin)
#Eccentricity array
es = np.zeros((N_t, N_v, N_bin), dtype=float)
es[0] = np.array([e_0]*N_bin)
N_broken = np.zeros((N_t, N_v), dtype=int)
for j in range(N_v):
        for i in range(1, N_t):
                A[i,j], es[i,j], N_broken[i,j] = MCEncounters(v_bins[j], n_p, dt, m1, m2, M_p, A[i-1,j], es[i-1,j], N_bin)
#Average
A_avg = np.zeros((N_t, N_v), dtype=float)
for j in range(N_v):
        for i in range(N_t):
                A_avg[i,j] = np.sum(A[i,j])/np.size(np.nonzero(A[i,j]))
        
#Plot average A against time
for i in range(N_v):
        plt.plot(t/(mega*year), A_avg[:,i]/(parsec), label=r'$\sigma = {}$km/s'.format(v_bins[i]/1000.0))
plt.xlabel('Time, Myr')
plt.ylabel('Average semi-major axis, pc')
plt.title('Average semi-major axis over time of {} binaries evolved using the Monte Carlo method for different PBH velocity dispersions'.format(N_bin))
plt.legend()
plt.show()

#Plot number of binaries broken against time
for i in range(N_v):
        plt.plot(t/(mega*year), N_broken[:,i], label=r'$\sigma = {}$km/s'.format(v_bins[i]/1000.0))
plt.xlabel('Time, Myr')
plt.ylabel('Number of binaries broken')
plt.title('Number of binaries broken over time for {} binaries evolved using the Monte Carlo method for different PBH velocity dispersions'.format(N_bin))
plt.legend()
plt.show()



















