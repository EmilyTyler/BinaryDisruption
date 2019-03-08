#To plot V, dV and cos(theta)
import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

e = np.array([0.0, 0.3, 0.5, 0.7, 0.9])
N_e = np.size(e)
N_data=5*10**6
data = np.zeros((N_data, N_e, 4), dtype=float)
G = 6.67*10.0**(-11.0)
M_p = 6.0*10.0**30.0
a = 10.0**5.0*1.5*10.0**11.0
v_rel = 2.2*10.0**5.0
b = 10.0**5.0*1.5*10.0**11.0

def load(filename, index):
	with open(filename) as csvfile:
		reader = csv.reader(csvfile, delimiter=',')
		row_number = 0
		for row in reader:
			data[row_number, index] = np.array([float(row[0]), float(row[1]), float(row[2]), float(row[3])])
			row_number += 1
	return

load('WSW_encounters_V_dV_theta_e0.csv', 0)
load('WSW_encounters_V_dV_theta_e0_3.csv', 1)
load('WSW_encounters_V_dV_theta_e0_5.csv', 2)
load('WSW_encounters_V_dV_theta_e0_7.csv', 3)
load('WSW_encounters_V_dV_theta_e0_9.csv', 4)

'''
N_bins = 500

v_min = 0.0
v_max = np.max(data[:,:,0])
dv = (v_max - v_min)/(N_bins-1)
v_bins = np.array([v_min + i*dv for i in range(N_bins)])
N_v = np.zeros((N_e, N_bins))
for k in range(N_e):
	for i in range(N_data):
		j = int(np.floor((data[i,k,0]-v_min)/dv))
		N_v[k,j] += 1
	plt.plot(v_bins, N_v[k], label=r'$e={}$'.format(e[k]))
plt.xlabel(r'$|\mathbf{V}|$, ms$^{-1}$')
plt.ylabel('Number of encounters')
plt.legend()
plt.show()


N_bins = 1000000

#Energy distribution
dE_min = np.min(data[:,:,3])
print('Min vdv = ', dE_min)
dE_max = np.max(data[:,:,3])
print('Max vdv = ', dE_max)
ddE = (dE_max - dE_min)/(N_bins-1)
dE_bins = np.array([dE_min + i*ddE for i in range(N_bins)])
N_dE = np.zeros((N_e, N_bins))
for k in range(N_e):
	for i in range(N_data):
		j = int(np.floor((data[i,k,3]-dE_min)/ddE))
		N_dE[k,j] += 1
	plt.plot(dE_bins, N_dE[k], label=r'$e={}$'.format(e[k]))
plt.xlabel(r'$\mathbf{V}\cdot\Delta\mathbf{V}$ term, J')
plt.ylabel('Number of encounters')
plt.legend()
plt.show()

'''
N_bins = 1000000

dv_min = 0.0
dv_max = np.max(data[:,:,1])
ddv = (dv_max - dv_min)/(N_bins-1)
dv_bins = np.array([dv_min + i*ddv for i in range(N_bins)])
N_dv = np.zeros((N_e,N_bins))
for k in range(N_e):
	for i in range(N_data):
		j = int(np.floor((data[i,k,1]-dv_min)/ddv))
		N_dv[k,j] += 1
	plt.plot(dv_bins, N_dv[k], label=r'$e={}$'.format(e[k]))
	plt.plot([2*G*M_p*a/(v_rel*b**2)*(1+e[k])]*N_bins, N_dv[k], color = 'black')
plt.xlabel(r'$|\Delta\mathbf{V}|$, ms$^{-1}$')
plt.ylabel('Number of encounters')
plt.legend()
plt.show()
'''

N_bins = 25

ct_min = -1.0
ct_max = 1.0
dct = (ct_max - ct_min)/(N_bins)
ct_bins = np.array([ct_min + i*dct for i in range(N_bins)])
N_ct = np.zeros((N_e,N_bins))
for k in range(N_e):
	for i in range(N_data):
		j = int(np.floor((data[i,k,2]-ct_min)/dct))
		N_ct[k,j] += 1
	plt.plot(ct_bins, N_ct[k], label=r'$e={}$'.format(e[k]))
plt.xlabel(r'$\cos(\theta)$')
plt.ylabel('Number of encounters')
plt.legend()
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
'''
