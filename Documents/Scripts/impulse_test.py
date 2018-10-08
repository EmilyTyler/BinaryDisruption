#To test the impulse approximation

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from matplotlib import pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
from scipy.constants import G, au, parsec, giga, year

from impulse_test_encounter import encounterGrid
from encounters import calc_b_max
from random_direction import randomDirection
from random_binary import setupRandomBinary

'''
#Initialise variables
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
M_p = 3.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 220.0 * 1000.0
#Number density of perturbers
n_p = rho/M_p

#Semi-major axes
a_min = 10.0**3.0 * au
a_max = 10.0**5.0 * au
#Number of a's to test
N_a = 20
#Impact parameters
#b_min = (np.pi*n_p*v_rms*(10.0*giga*year))**(-0.5)
#b_max = calc_b_max(M_p, v_rms, a_max, m1, m2)
b_min = a_min
b_max = a_max
#Number of b's to test
N_b = 20

#Number of encounters per each pair of values
#TAKES 2000 MINUTES TO RUN for 20, 20, 1000
N_enc = 1000



dloga = (np.log(a_max)-np.log(a_min))/N_a
a_bins = np.array([a_min*np.exp(dloga*i) for i in range(N_a)])
dlogb = (np.log(b_max)-np.log(b_min))/N_b
b_bins = np.array([b_min*np.exp(dlogb*i) for i in range(N_b)])
a_frac_avg, E_frac_avg, a_bins, b_bins = encounterGrid(m1, m2, v_rms, e, M_p, a_min, a_max, N_a, b_min, b_max, N_b, N_enc)

#Save
np.savez('impulsecontour_2BHTenc.npz', a_frac_avg=a_frac_avg, E_frac_avg=E_frac_avg, a_bins=a_bins, b_bins=b_bins)
'''
#Load
a_frac_avg = np.load('impulsecontour_hypenc.npz')['a_frac_avg']
E_frac_avg = np.load('impulsecontour_hypenc.npz')['E_frac_avg']
a_bins = np.load('impulsecontour_hypenc.npz')['a_bins']
b_bins = np.load('impulsecontour_hypenc.npz')['b_bins']

#Contour plot
'''
#Linear
plt.title('Average fractional error in semi-major axis due to impulse approximation')
plt.contourf(a_bins/(1.5*10.0**11.0), b_bins/(3.086*10.0**16.0), np.transpose(a_frac_avg))
plt.ylabel('Impact parameter, pc')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.colorbar()
plt.show()
'''

'''
#Symlog
plt.title('Average fractional error in semi-major axis due to impulse approximation')
ax = plt.gca()
lt = np.min(np.absolute(a_frac_avg[np.nonzero(a_frac_avg)]))
pcm = ax.pcolormesh(a_bins/au, b_bins/au, np.transpose(a_frac_avg), norm=colors.SymLogNorm(linthresh=lt))
plt.colorbar(pcm)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()

'''

#Log absolute value
plt.title('Absolute average fractional error in semi-major axis due to impulse approximation (double Bahcall et al. equations)', wrap=True)
ax = plt.gca()
cs = ax.contourf(a_bins/au, b_bins/au, np.transpose(np.absolute(a_frac_avg)), locator=ticker.LogLocator())
plt.colorbar(cs)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()

plt.title('Absolute average fractional error in semi-major axis due to impulse approximation (double Bahcall et al. equations)', wrap=True)
ax = plt.gca()
cs = ax.pcolormesh(a_bins/au, b_bins/au, np.transpose(np.absolute(a_frac_avg)), norm=colors.LogNorm())
plt.colorbar(cs)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()

#Log absolute value
plt.title('Absolute average fractional error in energy due to impulse approximation (double Bahcall et al. equations)', wrap=True)
ax = plt.gca()
norm= colors.LogNorm(vmin=10.0**(-7.0),vmax=10.0**2.0)
cs = ax.contourf(a_bins/au, b_bins/au, np.transpose(np.absolute(E_frac_avg)), norm=norm)
plt.colorbar(cs)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()

plt.title('Absolute average fractional error in energy due to impulse approximation (hyperbolic equations)', wrap=True)
ax = plt.gca()
cs = ax.pcolormesh(a_bins/au, b_bins/au, np.transpose(np.absolute(E_frac_avg)), norm=colors.LogNorm(), vmin=10.0**(-6.0), vmax=10.0**1.0)
plt.colorbar(cs)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()

'''
plt.title('Sign of average fractional error in semi-major axis due to impulse approximation')
plt.contourf(a_bins/au, b_bins/au, np.transpose(np.sign(a_frac_avg)))
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.colorbar()
plt.show()



#Plot time ratio test for impulse approximation
t_T = np.zeros((N_a, N_b))
for i in range(N_a):
        t_T[i] = np.array([(1.0 + e)/(2.0*np.pi*v_rms)*np.sqrt(G*(m1 + m2)/a_bins[i])]*N_b)               
plt.title('Time taken for PBH to travel between points of closest approach divided by orbital period')
ax = plt.gca()
cs = ax.contourf(a_bins/au, b_bins/au, np.transpose(t_T), locator=ticker.LogLocator())
plt.colorbar(cs)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()

#Plot better time ratio test 
t_T2 = np.zeros((N_a, N_b))
for i in range(N_a):
        for j in range(N_b):
                if M_p*a_bins[i]**2.0/(10.0**(-6.0)*np.min([m1,m2])) - b_bins[j]**2.0 < 0.0:
                        t_T2[i,j] = 0.0
                else:
                        t_T2[i,j] = 2.0/v_rms*np.sqrt(M_p*a_bins[i]**2.0/(10.0**(-6.0)*np.min([m1,m2])) - b_bins[j]**2.0)/(2.0*np.pi*np.sqrt(a_bins[i]**3.0/(G*(m1+m2))))
plt.title('Duration of PBH encounter with binary divided by orbital period')
ax = plt.gca()
cs = ax.contourf(a_bins/au, b_bins/au, np.transpose(t_T2), locator=ticker.LogLocator())
plt.colorbar(cs)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()

#Plot third time ratio test
t_T3 = np.zeros((N_a, N_b))
for i in range(N_a):
        P = 2.0 * np.pi * np.sqrt(a_bins[i]**3.0/(G*(m1+m2)))
        b_max = calc_b_max(M_p, v_rms, a_bins[i], m1, m2)
        for j in range(N_b):  
                if b_bins[j] < b_max:
                        t_T3[i,j] = 2.0 * b_bins[j] / (v_rms * P)
plt.title('Crossing time of encounter divided by orbital period')
ax = plt.gca()
cs = ax.contourf(a_bins/au, b_bins/au, np.transpose(t_T3), locator=ticker.LogLocator())
plt.colorbar(cs)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()


#Plot deltaV/V for PBH velocity
V_frac = np.zeros((N_b, N_a))
for i in range(N_b):
        V_frac[i] = np.array([2.0*m1*G/(b_bins[i]*v_rms**2.0)]*N_a)
plt.title('PBH fractional velocity difference after an encounter with one star')
ax = plt.gca()
cs = ax.contourf(a_bins/au, b_bins/au, V_frac, locator=ticker.LogLocator())
plt.colorbar(cs)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()
'''
'''
#Crossing time over period
t_crossing_P = np.zeros((N_a, N_b))
b_star = np.zeros(2)
#Plot maximum crossing time
for i in range(N_a):
        #Orbital period
        P = 2.0 * np.pi * np.sqrt(a_bins[i]**3.0/(G*(m1+m2)))
        b_max = calc_b_max(M_p, v_rms, a_bins[i], m1, m2)
        for j in range(N_b):
                #if b_bins[j] < b_max:
                t_crossing_P[i,j] = 2.0*b_bins[j]/(v_rms*P)
                for k in range(N_enc):
                        #Perturber velocity vector
                        v_vec = v_rms * randomDirection()
                        #Setup random binary
                        X = setupRandomBinary(a_bins[i], e, m1, m2)
                        #Centre of mass vector
                        R = (m1*X[0] + m2*X[1])/(m1 + m2)
                        #Find impact parameter vector
                        b_vec = np.dot(R,v_vec)/v_rms**2.0*v_vec - R
                        b_vec_norm = np.sqrt(b_vec[0]**2.0+b_vec[1]**2.0+b_vec[2]**2.0)
                        if b_vec_norm > 0.0:
                                b_vec = b_bins[j] * b_vec/b_vec_norm
                        #Impact parameters for individual stars
                        for l in range(2):
                                b_star_vec = (np.dot(X[l],v_vec) - np.dot(b_vec,v_vec))/v_rms**2.0 * v_vec + b_vec - X[l]
                                b_star[l] = np.linalg.norm(b_star_vec)
                        t_crossing_P[i,j] = np.amax([t_crossing_P[i,j], 2.0*b_star[0]/(v_rms*P), 2.0*b_star[1]/(v_rms*P)])                        
plt.title('Maximum crossing time of encounter divided by orbital period')
ax = plt.gca()
cs = ax.contourf(a_bins/au, b_bins/au, np.transpose(t_crossing_P), locator=ticker.LogLocator())
plt.colorbar(cs)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show() 

plt.title('Maximum crossing time of encounter divided by orbital period')
ax = plt.gca()
cs = ax.pcolormesh(a_bins/au, b_bins/au, np.transpose(t_crossing_P), norm=colors.LogNorm())
plt.colorbar(cs)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()   
        
'''
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        


