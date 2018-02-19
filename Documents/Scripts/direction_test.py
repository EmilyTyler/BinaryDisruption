#To test the randomness of the randomly generated b and v vectors

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

import matplotlib.pyplot as plt
from encounters import impactAndVelocityVectors
from frequency import calcFrequency

#Find theta in spherical polars
def theta(x):
        return np.arccos(x[2]/np.linalg.norm(x))

#Find phi
def phi(x):
        return np.arctan(x[1]/x[0])

#Impact parameter
b = 1.0
#Perturber speed
v = 1.0
#Number of vectors to generate
N_vec = 100000

v_vectors = np.zeros((N_vec,3))
b_vectors = np.zeros((N_vec,3))
v_phis = np.zeros(N_vec)
b_phis = np.zeros(N_vec)
for i in range(N_vec):
        b_vectors[i], v_vectors[i] = impactAndVelocityVectors(b, v)
        b_phis[i] = phi(b_vectors[i])
        v_phis[i] = phi(v_vectors[i])
        
plt.title('z test for v')        
vz_bins, N_vz, d_vz = calcFrequency(v_vectors[:,2], 100)
plt.plot(vz_bins, N_vz/N_vec/d_vz)
plt.show()

plt.title('z test for b') 
bz_bins, N_bz, d_bz = calcFrequency(b_vectors[:,2], 100)
plt.plot(bz_bins, N_bz/N_vec/d_bz)
plt.plot()
plt.show()

plt.title('phi test for v')
vphi_bins, N_vphi, d_vphi = calcFrequency(v_phis, 100)
plt.plot(vphi_bins, N_vphi/N_vec/d_vphi)
plt.show()

plt.title('phi test for b')
bphi_bins, N_bphi, d_bphi = calcFrequency(b_phis, 100)
plt.plot(bphi_bins, N_bphi/N_vec/d_bphi)
plt.show()

plt.title('x and y test for v')
u = np.sqrt(1.0 - v_vectors[:,0]**2.0 - v_vectors[:,1]**2.0)
vu_bins, N_vu, d_vu = calcFrequency(u, 100)
plt.plot(vu_bins, N_vu/N_vec/d_vu)
plt.show()

plt.title('x and y test for b')
u = np.sqrt(1.0 - b_vectors[:,0]**2.0 - b_vectors[:,1]**2.0)
bu_bins, N_bu, d_bu = calcFrequency(u, 100)
plt.plot(bu_bins, N_bu/N_vec/d_bu)
plt.show()













