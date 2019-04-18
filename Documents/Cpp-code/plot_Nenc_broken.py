import matplotlib.pyplot as plt
import csv
import numpy as np
from scipy.constants import au
plt.rc('font', family='serif')

def loadData(filename, plot_label, plot_color, marker_style='x'):
	N_enc = np.zeros(0)
	a_ini = np.zeros(0)
	with open(filename) as csvfile:
	        reader = csv.reader(csvfile, delimiter=',')
	        for row in reader:
	        	N_enc = np.append(N_enc, float(row[0]))
	        	a_ini = np.append(a_ini, float(row[1]))
	plt.scatter(a_ini/au, N_enc, label = plot_label, color=plot_color, marker=marker_style)

loadData('N_enc_broken_dist_10Msol.csv', plot_label=r'$M_p=10M_\odot$', plot_color='red')
loadData('N_enc_broken_dist_100Msol.csv', plot_label=r'$M_p=100M_\odot$', plot_color='forestgreen')
loadData('N_enc_broken_dist_1000Msol.csv', plot_label=r'$M_p=1000M_\odot$', plot_color='dodgerblue')

plt.xlabel('Initial semi-major axis, au')
plt.ylabel('Number of encounters before binary disruption')
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlim([10.0**3.0, 4*10.0**5.0])
plt.legend()
plt.show()