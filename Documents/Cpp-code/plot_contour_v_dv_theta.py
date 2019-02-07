import numpy as np
import csv
from matplotlib import pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors

N_data = 10**7

data = np.zeros((N_data,3), dtype=float)

with open('WSW_encounters_V_dV_theta.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	row_number = 0
	for row in reader:
		data[row_number,0] = float(row[0])
		data[row_number,1] = float(row[1])
		data[row_number,2] = float(row[2])
		row_number += 1

#V bins
N_v = 100
v_min = np.min(data[:,0])
v_max = np.max(data[:,0])
dv = (v_max-v_min)/(N_v-1)
v_bins = np.array([v_min + i*dv for i in range(N_v)])

#dV bins
N_dv = 100
dv_min = np.min(data[:,1])
dv_max = np.max(data[:,1])
ddv = (dv_max-dv_min)/(N_dv-1)
dv_bins = np.array([dv_min + i*ddv for i in range(N_dv)])

#Sort data
Avg_costh = np.zeros((N_v, N_dv), dtype=float)
N_avg_costh = np.zeros((N_v, N_dv), dtype=int)

for i in range(N_data):
	#v index
	j = int(np.floor((data[i,0]-v_min)/dv))
	#dv index
	k = int(np.floor((data[i,1]-dv_min)/ddv))
	#
	Avg_costh[j,k] += data[i,2]
	N_avg_costh[j,k] += 1

#Normalise
Avg_costh /= N_avg_costh

print('N_avg_costh =', N_avg_costh)

#Contour plot
ax = plt.gca()
#cs = ax.contourf(dv_bins, v_bins, Avg_costh)
cs = ax.contourf(dv_bins, v_bins, np.abs(Avg_costh), locator=ticker.LogLocator())
clb = plt.colorbar(cs)
clb.set_label(r'$\log($abs$(\cos(\theta)))$')
plt.xlabel(r'$|\Delta\mathbf{V}|$, ms$^{-1}$')
plt.ylabel(r'$|\mathbf{V}|$, ms$^{-1}$')
#plt.xscale('log')
#plt.yscale('log')
plt.show()