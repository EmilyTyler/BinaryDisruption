#To solve the PDE:
#dN/dt = x d^2N/dx^2 + 7/2 dN/dx

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import au, G, parsec, giga, year
from mpl_toolkits import mplot3d

#Total mass of binaries
M = 4.0*10.0**30.0
#Initial distribution parameter
alpha = 3.0
#Semi-major axes
a_min = 10.0**3.0*au
a_max = 10.0**6.0*au
#Mass density of perturbers
rho = 0.009 * 2.0*10.0**30.0/(parsec**3.0)
#Mass of perturbers
M_p = 2.0*10.0**30.0
#Number density of perturbers
n_p = rho / M_p
#Relative velcity dispersion
v_rel = 2.2*10.0**5.0
#Coulomb logarithm
#coul_log = ln(b_max/b_FP)
coul_log = 1.0
#End time of simulation
t_end = 10.0 * giga * year
#Number of binaries
N_bin = 10**6



#First order diffusion coefficient
#e_1/P in Weinberg et al
epsilon = 4.0*(2.0*np.pi)**0.5*n_p*G**2.0*M_p**2.0/v_rel * coul_log
#Minimum energy
#E_0 = abs(E(a_max))
E_0 = G*M/(2.0*a_max)
#Dimensionless time steps
N_tau = 10
tau_min = 0.0
tau_max = 2.0*epsilon*t_end/(3.0*E_0)
d_tau = (tau_max - tau_min)/N_tau
tau_bins = np.array([tau_min + i*d_tau for i in range(N_tau)])
#Dimensionless energy steps
N_x = 10
x_min = 1.0
x_max = G*M/(2.0*a_min*E_0)
dx = (x_max - x_min)/N_x
x_bins = np.array([x_min + i*dx for i in range(N_x)])
#Initial distribution
if alpha == 2.0:
        c = np.log(a_min)/np.log(a_max/a_min)
        a_ini = (a_max/a_min)**(np.random.random(N_bin) + c)
else:
        a_ini = (np.random.random(N_bin)*(a_max**(2.0-alpha) - a_min**(2.0-alpha)) + a_min**(2.0-alpha))**(1.0/(2.0-alpha))
x_ini = G*M/(2.0*a_ini*E_0)
#Matrix of N values
N = np.zeros((N_tau, N_x), dtype=float)
#Initialize
for val in x_ini:
        if val == x_max:
                i = N_x - 1
        else:
                i = int(np.floor((val - x_min)/dx))
                N[0,i] += 1


#
dN_dx = np.zeros(N_x)
d2N_dx2 = np.zeros(N_x)
for i in range(1,N_tau):
        for j in range(N_x-1):
                dN_dx[j] = (N[i-1,j+1] - N[i-1,j])/dx
        dN_dx[N_x-1] = - N[i-1,N_x-1]/dx
        for j in range(N_x-1):
                d2N_dx2[j] = (dN_dx[j+1] - dN_dx[j])/dx
        d2N_dx2[N_x-1] = - dN_dx[N_x-1]/dx
        N[i] = d_tau * (x_bins*d2N_dx2 + 7.0/2.0 * dN_dx)
        
#Plot N
#Timesteps in giga years
t_bins = 3.0*E_0*tau_bins/(2.0*epsilon) /(giga*year)
#Semi-major axis in au
a_bins = G*M/(2.0*x_bins*E_0)/au
for i in range(0,N_tau):
        plt.loglog(a_bins, N[i], label='t={}Gyr'.format(t_bins[i]))
plt.xlabel('Semi-major axis, au')
plt.ylabel('Number of binaries')
plt.legend()
plt.show()


























