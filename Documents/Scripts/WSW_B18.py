#To compare my simulations with WSW B18

import numpy as np
import os
os.system("python setup.py build_ext --inplace")
import matplotlib.pyplot as plt

from scipy.constants import parsec, au, giga, year, G
from scipy.special import jv, iv
from scipy.integrate import quad, quadrature, romberg

#Total mass of binaries
M_b = 4.0*10.0**30.0
#Initial distribution parameter
alpha = 1.0
#Mass density of perturbers
rho = 0.009 * 2.0*10.0**30.0/(parsec**3.0)
#Mass of perturbers
M_p = 1.0 * 2.0*10.0**30.0
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
N_bin = 10**4

#Semi-major axis bins
N_a = 10
a_min = 10.0**3.0*au
a_max = 10.0**6.0*au
da = (a_max - a_min)/N_a
a_bins = np.array([a_min + i*da for i in range(N_a)])

#Time
N_t = 11
dt = t_end/(N_t-1.0)
t = np.array([i*dt for i in range(N_t)])

#Initial distribution
if alpha == 2.0:
        c = np.log(a_min)/np.log(a_max/a_min)
        a_ini = (a_max/a_min)**(np.random.random(N_bin) + c)
else:
        a_ini = (np.random.random(N_bin)*(a_max**(2.0-alpha) - a_min**(2.0-alpha)) + a_min**(2.0-alpha))**(1.0/(2.0-alpha))
        
#Number of binaries with a=a_bins[i] at time t[j]: N[j,i]
N = np.zeros((N_t, N_a), dtype=float)

#First order diffusion coefficient
#e_1/P in Weinberg et al
epsilon = 8.0*(2.0*np.pi)**0.5*n_p*G**2.0*M_p**2.0/v_rel * coul_log
#Minimum energy
#E_0 = abs(E(a_max))
E_0 = G*M_b/(2.0*a_max)

def F_k(k,x):
        return (2.0*(k*x)**0.5)**(-5.0/2.0) * (jv(5.0/2.0, 2.0*(k*x)**0.5)*jv(-5.0/2.0, 2.0*(k)**0.5) - jv(5.0/2.0, 2.0*(k)**0.5)*jv(-5.0/2.0, 2.0*(k*x)**0.5))        

def integrand_B13(k, x_0, tau, x):
        return 512.0*np.pi*x_0**(5.0/2.0)*(np.exp(-k*tau)*k**5.0*F_k(k,x_0)*F_k(k,x))/(9.0 + 12.0*k + 16.0*k**2.0)      
                
def B13(x, x_0, tau):
        value, error = quad(integrand_B13, 0.001, 10.0**10.0, args=(x_0, tau, x))
        #print('B13:', value, error)
        return value     

#Initialize
for val in a_ini:
        if val == a_max:
                i = N_a - 1
        else:
                i = int(np.floor((val - a_min)/da))
                N[0,i] += 1

tau = 2.0*epsilon*t/(3.0*E_0)
x_bins = G*M_b/(2.0*a_bins*E_0)
print('x_bins =', x_bins)
dx = np.zeros(N_a)
for i in range(N_a-1):
        dx[i] = x_bins[i+1] - x_bins[i]
dx[N_a-1] = G*M_b/(2.0*a_max*E_0) - x_bins[N_a-1]

N_initial = np.zeros((N_t,N_a))
N_initial[0] = np.array([N[0,i] for i in range(N_a)]) 

'''
print('Here')
N_k = 100
k_min = 0.001
k_max = 10.0**10.0
dk = (k_max - k_min)/N_k
k_bins = np.array([k_min + i*dk for i in range(N_k)])
integ1 = np.zeros(N_k)
for i in range(N_k):
        integ1[i] = integrand_B13(k_bins[i], x_bins[0], tau[0], x_bins[2])
plt.plot(k_bins, integ1)
plt.show()
print('integ1 =', integ1)

print('Here 2')

integ2 = np.zeros(N_a)
for i in range(N_a):
        integ2[i] = B13(x_bins[i], x_bins[0], tau[0])
plt.plot(x_bins, integ2)
plt.show()
print('integ2 =', integ2)

print('Here 3')

for i in range(1,N_t):       
        for j in range(N_a):               
                for k in range(N_a-1):
                        N[i,k] += N[0,j] * quad(B13, x_bins[k+1], x_bins[k], args=(x_bins[j], tau[i]))[0]
                        #print('this =', N[0,j] * quad(B13, x_bins[k], x_bins[k+1], args=(x_bins[j], tau[i]))[0])
                        #print('this =', N[0,j] * quad(B13, x_bins[k], x_bins[k+1], args=(x_bins[j], tau[i]))[1])
                N[i,N_a-1] += N[0,j] * quad(integrand2, x_bins[N_a-1], 1.0, args=(x_bins[j], tau[i]))[0]

print('N =', N)

for i in range(0,N_t,N_t//5):
        plt.loglog(a_bins/au, N[i]/(da/au), label='t={}Gyr'.format(t[i]/(giga*year)))
plt.xlabel('Semi-major axis, au')
plt.ylabel(r'Number of binaries, au$^{-1}$')
plt.gca().set_ylim(bottom=0.1)
plt.legend()
plt.show()

print('Here 4')
'''
N = N_initial

def B18(x_0, x, t):
        print('B18 calculation')
        print('iv(5.0/2.0, 2.0*np.sqrt(x*x_0)/t) =', iv(5.0/2.0, 2.0*np.sqrt(x*x_0)/t))
        print('1.0/t * (x_0/x)**(5.0/4.0) * np.exp(-(x + x_0)/t) =', 1.0/t * (x_0/x)**(5.0/4.0) * np.exp(-(x + x_0)/t))
        print('np.exp(-(x + x_0)/t) * iv(5.0/2.0, 2.0*np.sqrt(x*x_0)/t) =', np.exp(-(x + x_0)/t) * iv(5.0/2.0, 2.0*np.sqrt(x*x_0)/t))
        result = 1.0/t * (x_0/x)**(5.0/4.0) * np.exp(-(x + x_0)/t) * iv(5.0/2.0, 2.0*np.sqrt(x*x_0)/t) 
        print('result =', result)
        return result

integ3 = np.zeros(N_a)
for i in range(N_a):
        integ3[i] = B18(x_bins[0], x_bins[i], tau[1]) 
print('B18 =', integ3)
#print(2.0*np.sqrt(x_bins*x_bins[0])/tau[10])
plt.plot(x_bins, integ3)
plt.show()


for i in range(1,N_t):
        for j in range(N_a):
                for k in range(N_a):
                        N[i,k] += N[0,j] * quad(B18, x_bins[j]+dx[j], x_bins[j], args=(x_bins[k], tau[i]))[0]
print('N =', N)
for i in range(0,N_t,N_t//5):
        plt.loglog(a_bins/au, N[i]/(da/au), label='t={}Gyr'.format(t[i]/(giga*year)))
plt.xlabel('Semi-major axis, au')
plt.ylabel(r'Number of binaries, au$^{-1}$')
plt.gca().set_ylim(bottom=0.1)
plt.legend()
plt.show()































