#To solve the PDE:
#dN/dt = x d^2N/dx^2 + 7/2 dN/dx

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import au, G, parsec, giga, year
from mpl_toolkits import mplot3d
from numpy.linalg import inv

print('Initialising')

#Total mass of binaries
M_b = 2.0*10.0**30.0
#Initial distribution parameter
alpha = 1.0
#Semi-major axes
a_min = 10.0**2.0*au
a_max = 10.0**6.0*au
#Mass density of perturbers
rho = 0.009 * 2.0*10.0**30.0/(parsec**3.0)
#Mass of perturbers
M_p = np.array([10.0, 100.0, 1000.0]) * 2.0*10.0**30.0
N_M = np.size(M_p)
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
N_bin = 10**9



#First order diffusion coefficient
#e_1/P in Weinberg et al
epsilon = 4.0*(2.0*np.pi)**0.5*n_p*G**2.0*M_p**2.0/v_rel * coul_log
#Minimum energy
#E_0 = abs(E(a_max))
E_0 = G*M_b/(2.0*a_max)
#Dimensionless time steps
N_tau = 10001
tau_min = 0.0
tau_max = 2.0*epsilon*t_end/(3.0*E_0)
d_tau = (tau_max - tau_min)/(N_tau - 1)
tau_bins = np.array([[tau_min + i*d_tau[j] for i in range(N_tau)] for j in range(N_M)])
#Dimensionless energy steps
N_x = 10000
x_min = 1.0
x_max = G*M_b/(2.0*a_min*E_0)
dx = (x_max - x_min)/N_x
x_bins = np.array([x_min + i*dx for i in range(N_x)])
#
r = d_tau/dx**2.0
print('r =', r)
print('Generating initial distribution')
#Initial distribution

if alpha == 2.0:
        c = np.log(a_min)/np.log(a_max/a_min)
        #a_ini = (a_max/a_min)**(np.random.random(N_bin) + c)
#else:
        #a_ini = (np.random.random(N_bin)*(a_max**(2.0-alpha) - a_min**(2.0-alpha)) + a_min**(2.0-alpha))**(1.0/(2.0-alpha))
#x_ini = G*M_b/(2.0*a_ini*E_0)

#Matrix of N values
N = np.zeros((N_M, N_tau, N_x), dtype=float)
#Initialize
for i in range(N_bin):
        if alpha == 2.0:
                a_ini = (a_max/a_min)**(np.random.random() + c)
        else:
                a_ini = (np.random.random()*(a_max**(2.0-alpha) - a_min**(2.0-alpha)) + a_min**(2.0-alpha))**(1.0/(2.0-alpha))
        val = G*M_b/(2.0*a_ini*E_0)
        if val == x_max:
                i = N_x - 1
        else:
                i = int(np.floor((val - x_min)/dx))
        N[:,0,i] += 1

print(N[0,0])

print('Generating evolution matrix')
#Crank-Nicolson Method
M = np.zeros((N_M,N_x,N_x))
k = np.zeros(N_x)
M[:,0,0] = 1.0/d_tau + x_bins[0]/dx**2.0
M[:,0,1] = - 7.0/(8.0*dx) - x_bins[0]/(2.0*dx**2.0)
M[:,N_x-1,N_x-2] = 7.0/(8.0*dx) - x_bins[N_x-1]/(2.0*dx**2.0)
M[:,N_x-1,N_x-1] = 1.0/d_tau + x_bins[N_x-1]/dx**2.0
for j in range(1,N_x-1):
                M[:,j,j-1] = 7.0/(8.0*dx) - x_bins[j]/(2.0*dx**2.0)
                M[:,j,j] = 1.0/d_tau + x_bins[j]/dx**2.0
                M[:,j,j+1] = - 7.0/(8.0*dx) - x_bins[j]/(2.0*dx**2.0)
M_inv = np.array([inv(M[i]) for i in range(N_M)])

print('Evolving in time')
for l in range(N_M):
        for i in range(1,N_tau):
                k[0] = N[l,i-1,0]/d_tau[l] + 7.0/(8.0*dx)*N[l,i-1,1] + x_bins[0]/(2.0*dx**2.0)*(N[l,i-1,1] - 2.0*N[l,i-1,0])
                k[N_x-1] = N[l,i-1,N_x-1]/d_tau[l] - 7.0/(8.0*dx)*N[l,i-1,N_x-2] + x_bins[N_x-1]/(2.0*dx**2.0)*(-2.0*N[l,i-1,N_x-1] + N[l,i-1,N_x-2])
                for j in range(1,N_x-1):
                        k[j] = N[l,i-1,j]/d_tau[l] + 7.0/(8.0*dx)*(N[l,i-1,j+1] - N[l,i-1,j-1]) + x_bins[j]/(2.0*dx**2.0)*(N[l,i-1,j+1] - 2.0*N[l,i-1,j] + N[l,i-1,j-1])
                N[l,i] = np.dot(M_inv[l],k)
#
print('Plotting')
#Plot N
#Timesteps in giga years
t_bins = np.array([3.0*E_0*tau_bins[i]/(2.0*epsilon[i]) /(giga*year) for i in range(N_M)])
#Semi-major axis in au
a_bins = G*M_b/(2.0*x_bins*E_0)/au
da = np.zeros(N_x)
for i in range(N_x-1):
        da[i] = a_bins[i-1] - a_bins[i]
da[N_x-1] = a_max - a_bins[N_x-1]
#
'''
for i in range(0,N_tau,N_tau//5):
        plt.loglog(a_bins, N[i]/da, label='t={}Gyr'.format(t_bins[i]))
plt.xlabel('Semi-major axis, au')
plt.ylabel(r'Number of binaries, au$^{-1}$')
plt.gca().set_ylim(bottom=0.1)
plt.legend()
plt.show()
'''
plt.loglog(a_bins, N[0,0]/da, label='Initial')
for i in range(N_M):
        plt.loglog(a_bins, N[i,N_tau-1]/da, label=r'$M={}M_\odot$'.format(int(M_p[i]/(2.0*10.0**30.0))))
plt.xlabel('Semi-major axis, au')
plt.ylabel(r'Number of binaries, au$^{-1}$')
plt.title('Semi-major axis distribution after {}Gyr found by numerically solving the Fokker-Planck equation (Weinberg et al B11)'.format(int(t_end/(giga*year))), wrap=True)
plt.axis([1000.0, 300000.0, 0.1, 600.0])
plt.legend()
plt.show()

print('Finished')






















