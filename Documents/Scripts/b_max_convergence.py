#Convergence tests for b_max

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from scipy.constants import au, parsec, giga, year, mega, G
from monte_carlo import draw_b, draw_vmaxwellian, MCEncounters
from encounters import encounterRate, encounter, calc_b_max
import matplotlib.pyplot as plt

#Initialise variables
#Semi-major axis, m
a = 10.0**4.0 * au
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 2.0*10.0**30.0
m2 = 2.0*10.0**30.0
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/(parsec**3.0)
#Mass of perturbers
M_p = 10.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rel = 220.0 * 1000.0
#Number density of perturbers
n_p = rho/M_p
#Simulation run time
T = 10.0 * giga * year

#b_max
b_max = calc_b_max(M_p, v_rel, a, m1, m2)
#b_max values
b_maxs = np.array([0.01*b_max, 0.1*b_max, 1.0*b_max])
b_maxs_size = np.size(b_maxs)
b_max_max = np.max(b_maxs)

#Plotting colors list
colors = ['red', 'forestgreen', 'blue', 'darkorange', 'darkorchid']
# Function to make custom legend
def make_legend(b_maxes, colours):
        ls = np.array(np.round(b_maxes/parsec, decimals=2), dtype=str)
        lc = np.array(colours)
        lc = np.delete(lc, range(np.size(ls),np.size(lc)))
        plt.figtext(0.27,0.85,r'$b_{\mathrm{max}}$/pc =')
        i = 0
        for s, c in zip(ls,lc):
                if i==(np.size(ls)-1):
                        plt.figtext(0.46+i*0.055,0.85,s,color=c)
                else:
                        plt.figtext(0.46+i*0.055,0.85,s+',',color=c)
                i += 1
  
'''
#Calculate average energy change for a large number of encounters for different values of b_max
print('Calculating average energy changes')
#Number of encounters for each value of b_max
#10**6 takes 15 minutes
#10**7 takes ~2.5 hours
N_enc = 10**7
#Average energy change
delta_E_avg = np.zeros(b_maxs_size)
#Initial energy
E_ini = -G*(m1+m2)/(2.0*a)
#Generate b values
b = draw_b(b_max_max, N_enc)
#Generate v values
v = draw_vmaxwellian(v_rel, 0.01*v_rel, 100.0*v_rel, N_enc)
#Implement encounters for different values of b_max
for i in range(b_maxs_size):
        b_subset = b[np.where(b <= b_maxs[i])]
        v_subset = v[np.where(b <= b_maxs[i])]
        for j in range(np.size(b_subset)):
                notBound, a_new, e_new = encounter(m1, m2, v_subset[j], b_subset[j], a, e, M_p)
                E_new = -G*(m1+m2)/(2.0*a_new)
                delta_E_avg[i] += (E_new - E_ini)
#Normalise
delta_E_avg /= N_enc
#Print results
for i in range(b_maxs_size):
        print('b_max =', b_maxs[i]/parsec, 'pc    ', 'dE_avg =', delta_E_avg[i])


#Plot delta a against t for a small population
print('Evolving small population with same initial semi-major axis')
#Number of binaries
#10 takes 10 minutes
N_bin = 10
#Number of time steps
N_t = 1000
#Time array
t = np.linspace(0.0, T, N_t)
dt = t[1] - t[0]
#Array of a's: A[i,j,k] is the semi-major axis for binary i at time t[j] with b_max=b_maxs[k]
A = np.zeros((N_bin, N_t, b_maxs_size))
#Initialise A
A[:,0,:] += a
#Array of es
es = np.zeros((N_bin, N_t, b_maxs_size))
#Initialise es
es[:,0,:] += e
for i in range(N_bin):
        for j in range(1,N_t):
                #Calculate number of encounters
                N_enc = np.random.poisson(dt*encounterRate(n_p, v_rel, 0.0, b_max_max, 0.01*v_rel, 100.0*v_rel))
                #Generate b values
                b = draw_b(b_max_max, N_enc)
                #Generate v values
                v = draw_vmaxwellian(v_rel, 0.01*v_rel, 100.0*v_rel, N_enc)
                #Implement encounters for different values of b_max
                for k in range(b_maxs_size):
                        b_subset = b[np.where(b <= b_maxs[k])]
                        v_subset = v[np.where(b <= b_maxs[k])]
                        A[i,j,k] = A[i,j-1,k]
                        if A[i,j,k] < 0.0:
                                break
                        for l in range(np.size(b_subset)):
                                notBound, A[i,j,k], es[i,j,k] = encounter(m1, m2, v_subset[l], b_subset[l], A[i,j,k], es[i,j,k], M_p)
                                if notBound:
                                        A[i,j,k] = -1.0
                                        es[i,j,k] = -1.0
                                        break
#Find change in semi-major axis
delta_a = np.zeros((N_bin, N_t, b_maxs_size))
delta_a[0] = A[0] - a
for j in range(1,N_t):
        delta_a[:,j,:] = A[:,j,:] - A[:,j-1,:]
#Set delta_a to zero for broken binaries
delta_a = np.where(A < 0.0, 0.0, delta_a)
#Plot delta_a against time for different b_maxs
for i in range(N_bin):
        for k in range(b_maxs_size):
                plt.plot(t/(mega*year), delta_a[i,:,k]/au, color=colors[k])                                
plt.title(r'Change in semi-major axis over time for different values of $b_{\mathrm{max}}$')
plt.xlabel('Time, Myr')
plt.ylabel(r'$\Delta a$, au')
make_legend(b_maxs, colors)
plt.show()
#Plot a against time for different b_maxs
for i in range(N_bin):
        for k in range(b_maxs_size):
                plt.plot(t[np.where(A[i,:,k] > 0.0)]/(mega*year), A[i,:,k][np.where(A[i,:,k] > 0.0)]/au, color=colors[k])                                
plt.title(r'Semi-major axis over time for different values of $b_{\mathrm{max}}$')
plt.xlabel('Time, Myr')
plt.ylabel('Semi-major axis, au')
make_legend(b_maxs, colors)
plt.show()

'''


#Plot evolved semi-major axis distribution
print('Evolving a distribution')
#Number of binary pairs for each b_max
N_bin = 10**3
print('Number of binaries =', N_bin)
#Semi-major axis array
#Draw a from distribution dN/dloga\propto a^{1-alpha}
alpha = 2.0
a_min = 10.0**3.0*au
a_max = 10.0**6.0*au
if alpha == 2.0:
        c = np.log(a_min)/np.log(a_max/a_min)
        a_ini = (a_max/a_min)**(np.random.random(N_bin) + c)
else:
        a_ini = (np.random.random(N_bin)*(a_max**(2.0-alpha) - a_min**(2.0-alpha)) + a_min**(2.0-alpha))**(1.0/(2.0-alpha))
#Eccentricity array
#Draw e from distribution uniform in e^2 between 0 and 1
e_ini = (np.random.random(N_bin))**(1.0/3.0)
#Final distributions
a_fin = np.zeros((b_maxs_size, N_bin))
e_fin = np.zeros((b_maxs_size, N_bin))
for i in range(b_maxs_size):
        #Evolve distribution in time
        a_fin[i], e_fin[i], N_broken = MCEncounters(v_rel, n_p, T, m1, m2, M_p, a_ini, e_ini, N_bin, prefactor=b_maxs[i]/b_max)
        print('b_max =', b_maxs[i]/parsec, 'pc    ', 'N_broken =', N_broken)
#Plot final and initial distributions
#Set up bins
#Number of bins
N_bins = 50
a_min = np.min([np.min(a_ini), np.min(a_fin[np.where(a_fin>0.0)])])
a_max = np.max([np.max(a_ini), np.max(a_fin)])
dloga = (np.log(a_max)-np.log(a_min))/N_bins
a_bins = np.array([a_min*np.exp(i*dloga) for i in range(N_bins)])
# e
e_min = np.min([np.min(e_ini), np.min(e_fin[np.where(e_fin>=0.0)])])
e_max = np.max([np.max(e_ini), np.max(e_fin)])
dloge = (np.log(e_max)-np.log(e_min))/N_bins
e_bins = np.array([e_min*np.exp(i*dloge) for i in range(N_bins)])
#Count frequency
N_a_ini = np.zeros(N_bins, dtype=int)
N_e_ini = np.zeros(N_bins, dtype=int)
N_a_fin = np.zeros((b_maxs_size, N_bins), dtype=int)
N_e_fin = np.zeros((b_maxs_size, N_bins), dtype=int)
#Frequency of initial distributions
for val in a_ini:
        if val == a_max:
                i = N_bins - 1
        else:
                i = int(np.floor(np.log(val/a_min)/dloga))
                N_a_ini[i] += 1
for val in e_ini:
        if val == e_max:
                i = N_bins - 1
        else:
                i = int(np.floor(np.log(val/e_min)/dloge))
                N_e_ini[i] += 1
#Frequency of final distributions
for k in range(b_maxs_size):     
        for val in a_fin[k]:
                if val > 0.0:
                        if val == a_max:
                                i = N_bins - 1
                        else:
                                i = int(np.floor(np.log(val/a_min)/dloga))
                                N_a_fin[k,i] += 1       
        for val in e_fin[k]:
                if val > 0.0:
                        if val == e_max:
                                i = N_bins - 1
                        else:
                                i = int(np.floor(np.log(val/e_min)/dloge))
                                N_e_fin[k,i] += 1
#Plot
plt.loglog(a_bins/au, N_a_ini, color = 'black')
for i in range(b_maxs_size):
        plt.loglog(a_bins/au, N_a_fin[i], color=colors[i])
plt.xlabel('Semi-major axis, au')
plt.ylabel('Number of binaries')
make_legend(b_maxs, colors)
plt.show()

plt.loglog(e_bins, N_e_ini, color='black')
for i in range(b_maxs_size):
        plt.loglog(e_bins, N_e_fin[i], color=colors[i])
plt.xlabel('Eccentricity')
plt.ylabel('Number of binaries')
make_legend(b_maxs, colors)
plt.show()

print('Finished')

                
                
                
                
                
                
                
                
                
                
                
                
                
                
        
        

