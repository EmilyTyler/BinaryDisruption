#To plot V, dV and cos(theta)
import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

data = np.zeros((10**6,4), dtype=float)

with open('WSW_encounters_V_dV_theta_e0.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	row_number = 0
	for row in reader:
		#data = np.append(data, np.array([[float(row[0]), float(row[1]), float(row[2])]]), axis=0)
		data[row_number] = np.array([float(row[0]), float(row[1]), float(row[2]), float(row[3])])
		row_number += 1

N_bins = 200
'''
v_min = np.min(data[:,0])
v_max = np.max(data[:,0])
dv = (v_max - v_min)/(N_bins-1)
v_bins = np.array([v_min + i*dv for i in range(N_bins)])
N_v = np.zeros(N_bins)
for i in range(np.size(data[:,0])):
	j = int(np.floor((data[i,0]-v_min)/dv))
	N_v[j] += 1
N_v /= np.size(data[:,0])
plt.plot(v_bins, N_v/dv)
#plt.plot(d_bins, 2.0*(d_bins+0.5*dd)/(b_max**2.0)
plt.xlabel(r'$|\mathbf{V}|$, ms$^{-1}$')
plt.ylabel(r'Probability density, m$^{-1}$s')
plt.show()
'''
#Energy distribution
'''
dE = np.zeros(10**6)
for i in range(10**6):
	dE[i] = data[i,0] * data[i,1] * data[i,2]
'''
dE_min = np.min(data[:,3])
dE_max = np.max(data[:,3])
ddE = (dE_max - dE_min)/(N_bins-1)
dE_bins = np.array([dE_min + i*ddE for i in range(N_bins)])
N_dE = np.zeros(N_bins)
for i in range(np.size(data[:,3])):
	j = int(np.floor((data[i,3]-dE_min)/ddE))
	N_dE[j] += 1
plt.plot(dE_bins, N_dE)
plt.xlabel(r'$\mathbf{V}\cdot\Delta\mathbf{V}$ term J')
plt.ylabel('Number of encounters')
plt.show()

dv_min = np.min(data[:,1])
dv_max = np.max(data[:,1])
ddv = (dv_max - dv_min)/(N_bins-1)
dv_bins = np.array([dv_min + i*ddv for i in range(N_bins)])
N_dv = np.zeros(N_bins)
for i in range(np.size(data[:,1])):
	j = int(np.floor((data[i,1]-dv_min)/ddv))
	N_dv[j] += 1
plt.plot(dv_bins, N_dv)
plt.xlabel(r'$|\Delta\mathbf{V}|$, ms$^{-1}$')
plt.ylabel('Number of encounters')
plt.show()

ct_min = np.min(data[:,2])
ct_max = np.max(data[:,2])
dct = (ct_max - ct_min)/(N_bins-1)
ct_bins = np.array([ct_min + i*dct for i in range(N_bins)])
N_ct = np.zeros(N_bins)
for i in range(np.size(data[:,2])):
	j = int(np.floor((data[i,2]-ct_min)/dct))
	N_ct[j] += 1
plt.plot(ct_bins, N_ct)
plt.xlabel(r'$\cos(\theta)$')
plt.ylabel('Number of encounters')
plt.show()


#Contour plot
N_dvct = np.zeros((N_bins, N_bins))
for i in range(np.size(data[:,1])):
	j = int(np.floor((data[i,1]-dv_min)/ddv))
	k = int(np.floor((data[i,2]-ct_min)/dct))
	N_dvct[j,k] += 1
levels = np.linspace(np.min(N_dvct), np.max(N_dvct), num=25)
#N_dv /= np.size(data[:,1])
cs = plt.contourf(dv_bins, ct_bins, np.transpose(N_dvct), levels=levels)
cbar = plt.colorbar(cs)
cbar.ax.set_ylabel('Number of encounters')
plt.xlabel(r'$|\Delta\mathbf{V}|$, ms$^{-1}$')
plt.ylabel(r'$\cos(\theta)$')
plt.show()

