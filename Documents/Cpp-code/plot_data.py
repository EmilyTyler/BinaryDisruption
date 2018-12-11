#To plot test data
import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


v_rel = 50.0
data = np.zeros(1000000, dtype=float)

with open('test_data.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	row_number = 0
	for row in reader:
		data[row_number] = float(row[0])
		row_number += 1

N_bins = 100

d_min = np.min(data)
d_max = np.max(data)
dd = (d_max - d_min)/(N_bins-1.0)
d_bins = np.array([d_min + i*dd for i in range(N_bins)])
N_d = np.zeros(N_bins)
for i in range(np.size(data)):
	j = int(np.floor((data[i]-d_min)/dd))
	N_d[j] += 1
N_d /= np.size(data)
plt.plot(d_bins, N_d/dd)
plt.plot(d_bins, (d_bins+0.5*dd)**3.0/(2.0*v_rel**4.0)*np.exp(-(d_bins+0.5*dd)**2.0/(2.0*v_rel**2.0)))
plt.show()










'''
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')

data = np.zeros((1000, 3), dtype=float)

with open('test_data.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	row_number = 0
	for row in reader:
		data[row_number, 0] = float(row[0])
		data[row_number, 1] = float(row[1])
		data[row_number, 2] = float(row[2])
		row_number += 1

for i in range(1000):
	ax1.scatter(data[i,0], data[i,1], data[i,2])
plt.show()
'''