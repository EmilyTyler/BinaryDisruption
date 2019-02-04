#To plot test data
import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

b_max = 3.086*10.0**16.0


#Vector direction tests
data = np.zeros((10**3, 6), dtype=float)

with open('WSW_encounters_V_dV_vectors.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	row_number = 0
	for row in reader:
		data[row_number, 0] = float(row[0])
		data[row_number, 1] = float(row[1])
		data[row_number, 2] = float(row[2])
		data[row_number, 3] = float(row[3])
		data[row_number, 4] = float(row[4])
		data[row_number, 5] = float(row[5])
		row_number += 1


#Normalise vectors
for i in range(np.size(data[:,0])):
	factor = np.sqrt(data[i,0]**2.0 + data[i,1]**2.0+data[i,2]**2.0)
	data[i,:3] /= factor
	factor = np.sqrt(data[i,3]**2.0 + data[i,4]**2.0+data[i,5]**2.0)
	data[i,3:] /= factor


#Plot 3D
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
for i in range(1000):
	ax1.scatter(data[i,0], data[i,1], data[i,2])
plt.show()

fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
for i in range(1000):
	ax1.scatter(data[i,3], data[i,4], data[i,5])
plt.show()



#Test theta = arctan(y/x), should be uniform [0,2pi)
theta = np.arctan2(data[:,4], data[:,3])
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
z = data[:,5]
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
