#To test the BHT encounter equations

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from matplotlib import pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
from scipy.constants import G, au, parsec, giga, year

from encounters import calc_b_max, impactAndVelocityVectors, dot_3d
from random_binary import setupRandomBinary
from orbital_elements import orbitalElements

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
v_rel = 220.0 * 1000.0
#Number density of perturbers
n_p = rho/M_p

#Semi-major axes
a_min = 10.0**3.0 * au
a_max = 10.0**5.0 * au
#Number of a's to test
N_a = 20
#Impact parameters
b_min = (np.pi*n_p*v_rel*(10.0*giga*year))**(-0.5)
b_max = calc_b_max(M_p, v_rel, a_max, m1, m2)
#Number of b's to test
N_b = 20

#Number of encounters per each pair of values
#TAKES 2000 MINUTES TO RUN for 20, 20, 1000
N_enc = 1000

#Semi-major axis bins
dloga = (np.log(a_max)-np.log(a_min))/N_a
a_bins = np.array([a_min*np.exp(dloga*i) for i in range(N_a)])
#Impact parameter bins
dlogb = (np.log(b_max)-np.log(b_min))/N_b
b_bins = np.array([b_min*np.exp(dlogb*i) for i in range(N_b)])

#Average fractional difference in a
a_frac_avg = np.zeros((N_a, N_b), dtype=float)
#Average fractional difference in energy
E_frac_avg = np.zeros((N_a, N_b), dtype=float)

#Star masses
m = np.array([m1, m2])
#90 degree deflection radius
b_90 = G*(M_p+m)/v_rel**2.0

for i in range(N_a):
        for j in range(N_b):
                for k in range(N_enc):
                        #Open binary
                        X = setupRandomBinary(a_bins[i], e, m1, m2)
                        #Find impact parameter vector and velocity vector
                        b_vec, v_vec = impactAndVelocityVectors(b_bins[j], v_rel)
                        #Velocity change from my equations:
                        v_my = np.zeros((2,3))
                        #Velocity change from BHT equations
                        v_BHT = np.zeros((2,3))
                        for l in range(2):
                                #Calculate impact parameter for this star
                                b_star = dot_3d(X[l],v_vec)/v_rel**2.0 * v_vec + b_vec - X[l]
                                b_star_norm = np.sqrt(b_star[0]**2.0 + b_star[1]**2.0 + b_star[2]**2.0)
                                #Calculate velocity change from my equations
                                v_perp = 2.0*M_p*v_rel/(m[l]+M_p) * (b_star_norm/b_90[l])/(1.0 + b_star_norm**2.0/b_90[l]**2.0) * (b_star/b_star_norm)
                                #Calculate velocity change in -v direction
                                v_parr = 2.0*M_p*v_rel/(m[l]+M_p) * 1.0/(1.0 + b_star_norm**2.0/b_90[l]**2.0) * (-v_vec/v_rel)
                                v_my[l] = v_parr + v_perp
                                #Calculate velocity change from BHT equations
                                if b_bins[j] > a_bins[i]:
                                        v_BHT[l] = G*M_p*a_bins[i]/(b_star_norm**2.0*v_rel) * (b_star/b_star_norm)
                                else:
                                        v_BHT[l] = G*M_p/(b_star_norm*v_rel) * (b_star/b_star_norm)
                        #Calculate new semi-major axes
                        notBound, a_my, e_my = orbitalElements(np.array([X[0], X[1], X[2]+v_my[0], X[3]+v_my[1]]), m1, m2)
                        notBound, a_BHT, e_BHT = orbitalElements(np.array([X[0], X[1], X[2]+v_BHT[0], X[3]+v_BHT[1]]), m1, m2)
                        #Calculate a_frac_avg
                        a_frac_avg[i,j] += (a_my-a_BHT)/a_BHT
                        #Calculate E_frac_avg
                        E_frac_avg[i,j] += a_BHT/a_my - 1.0

#Normalise a_frac_avg
a_frac_avg /= N_enc
#Normalise E_frac_avg
E_frac_avg /= N_enc

#Log absolute value
plt.title('Absolute average fractional difference in semi-major axis', wrap=True)
ax = plt.gca()
cs = ax.contourf(a_bins/au, b_bins/au, np.transpose(np.absolute(a_frac_avg)), locator=ticker.LogLocator())
plt.colorbar(cs)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()

plt.title('Absolute average fractional difference in semi-major axis', wrap=True)
ax = plt.gca()
cs = ax.pcolormesh(a_bins/au, b_bins/au, np.transpose(np.absolute(a_frac_avg)), norm=colors.LogNorm())
plt.colorbar(cs)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()

#Log absolute value
plt.title('Absolute average fractional difference in energy')
ax = plt.gca()
cs = ax.contourf(a_bins/au, b_bins/au, np.transpose(np.absolute(E_frac_avg)), locator=ticker.LogLocator(), wrap=True)
plt.colorbar(cs)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()

plt.title('Absolute average fractional difference in energy', wrap=True)
ax = plt.gca()
cs = ax.pcolormesh(a_bins/au, b_bins/au, np.transpose(np.absolute(E_frac_avg)), norm=colors.LogNorm())
plt.colorbar(cs)
plt.ylabel('Impact parameter, au')
plt.xlabel('Semi-major axis, au')
plt.xscale('log')
plt.yscale('log')
plt.show()