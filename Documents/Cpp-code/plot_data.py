#To plot test data
import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.constants import au

'''
data = np.zeros(10**7, dtype=float)

with open('test_data.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	row_number = 0
	for row in reader:
		data[row_number] = float(row[0])
		row_number += 1

N_bins = 1000

d_min = np.min(data)
d_max = np.max(data)
print('min =', d_min)
print('max =', d_max)
dd = (d_max - d_min)/(N_bins)
d_bins = np.array([d_min + i*dd for i in range(N_bins)])
N_d = np.zeros(N_bins)
for i in range(np.size(data)):
	j = int(np.floor((data[i]-d_min)/dd))
	if (data[i] == d_max):
		j = N_bins-1
	N_d[j] += 1
N_d /= np.size(data)
#Move bins into centre for plotting and calculations
d_bins += 0.5*dd
plt.plot(d_bins/au, N_d*au/dd, label='Simulation')
e=0.3
a = 10.0**5.0 * au
plt.plot(d_bins/au, au*1.0/np.pi*d_bins/a*(a**2.0*e**2.0-(d_bins-a)**2.0)**(-1/2), label='Analytic')
#plt.plot(d_bins, 1.0/(2.0*np.pi)*(1.0-e)**(3.0/2.0)*(1.0+e*np.cos(d_bins))**(-2.0))
plt.xlabel(r'$r$, au')
plt.ylabel(r'Probability density, au$^{-1}$')
plt.legend()
plt.show()




#Vector direction tests
data = np.zeros((5*10**6, 3), dtype=float)

with open('test_data.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	row_number = 0
	for row in reader:
		data[row_number, 0] = float(row[0])
		data[row_number, 1] = float(row[1])
		data[row_number, 2] = float(row[2])
		row_number += 1

#Plot 3D
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
for i in range(1000):
	ax1.scatter(data[i,0], data[i,1], data[i,2])
plt.show()

#Test theta = arctan(y/x), should be uniform [0,2pi)
theta = np.arctan2(data[:,1], data[:,0])
N_bins = 100
t_min = 0.0
t_max = 2.0*np.pi
dt = (t_max - t_min)/(N_bins)
t_bins = np.array([t_min + i*dt for i in range(N_bins)])
N_t = np.zeros(N_bins)
for i in range(np.size(theta)):
	j = int(np.floor((theta[i]-t_min)/dt))
	N_t[j] += 1
plt.plot(t_bins, N_t/dt/np.size(theta))
plt.plot(t_bins, [1.0/N_bins/dt]*N_bins)
plt.show()
#Test z, should be uniform [-1,1]
z = data[:,2]
N_bins = 100
z_min = -1.0
z_max = 1.1
dz = (z_max - z_min)/(N_bins-1)
z_bins = np.array([z_min + i*dz for i in range(N_bins)])
print('z_bins =', z_bins)
N_z = np.zeros(N_bins)
for i in range(np.size(z)):
	j = int(np.floor((z[i]-z_min)/dz))
	N_z[j] += 1
plt.plot(z_bins, N_z/dz/np.size(z))
plt.plot(z_bins, [1.0/N_bins/dz if (z_bins[i]<=1.0) else 0.0 for i in range(N_bins)])
plt.show()
'''

#2D vdv and i distribution
N_data=10**6
data = np.zeros((N_data, 2), dtype=float)

def load(filename):
	with open(filename) as csvfile:
		reader = csv.reader(csvfile, delimiter=',')
		row_number = 0
		for row in reader:
			data[row_number] = np.array([float(row[0]), float(row[1])])
			row_number += 1
	return

load('vdv_i_b10e5au.csv')

N_bins = 100
vdv_min = np.min(data[:,0])
vdv_max = np.max(data[:,0])
dvdv = (vdv_max - vdv_min)/(N_bins-1)
vdv_bins = np.array([vdv_min + i*dvdv for i in range(N_bins)])
#Move to centre of bins for plotting
vdv_bins += 0.5*dvdv
N_vdv = np.zeros(N_bins)
for i in range(N_data):
	j = int(np.floor((data[i,0]-vdv_min)/dvdv))
	N_vdv[j] += 1
plt.plot(vdv_bins, N_vdv)
plt.xlabel(r'$\mathbf{V}\cdot\Delta\mathbf{V}$, J')
plt.ylabel('Number of encounters')
plt.legend()
plt.show()

N_bins = 100
i_min = 0.0
i_max = np.pi/2.0
di = (i_max - i_min)/(N_bins-1)
i_bins = np.array([i_min + i*di for i in range(N_bins)])
#Move to centre of bins for plotting
i_bins += 0.5*di
N_i = np.zeros(N_bins)
for i in range(N_data):
	j = int(np.floor((data[i,1]-i_min)/di))
	N_i[j] += 1
plt.plot(i_bins, N_i)
plt.xlabel('Perturber inclination, rad')
plt.ylabel('Number of encounters')
plt.legend()
plt.show()

#Contour plot
N_bins = 100
N_vdvi = np.zeros((N_bins, N_bins))
for i in range(np.size(data[:,0])):
	j = int(np.floor((data[i,0]-vdv_min)/dvdv))
	k = int(np.floor((data[i,1]-i_min)/di))
	N_vdvi[j,k] += 1
levels = np.linspace(np.min(N_vdvi), np.max(N_vdvi), num=25)
#N_dv /= np.size(data[:,1])
cs = plt.contourf(vdv_bins, i_bins, np.transpose(N_vdvi), levels=levels)
cbar = plt.colorbar(cs)
cbar.ax.set_ylabel('Number of encounters')
ax = plt.gca()
ax.set_yticks([0., .25*np.pi, .5*np.pi])
ax.set_yticklabels(["$0$", r"$\frac{1}{4}\pi$", r"$\frac{1}{2}\pi$"])
plt.xlabel(r'$\mathbf{V}\cdot\Delta\mathbf{V}$, J')
plt.ylabel(r'Perturber inclination $i$, rad')
plt.show()


