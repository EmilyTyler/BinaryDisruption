import matplotlib.pyplot as plt
import csv
import numpy as np
from scipy.constants import au
plt.rc('font', family='serif')

a_min = 10.0**3.0 * au
a_max = 2.0*10.0**5.0 * au
N_bins = 40

dloga = (np.log(a_max)-np.log(a_min))/(N_bins)
a_bins = np.array([a_min*np.exp(i*dloga) for i in range(N_bins)])

N_bin = 10**5

def loadData(filename, plot_label, plot_initial, plot_color, plot_linestyle):
	a_ini = np.zeros(N_bin)
	e_ini = np.zeros(N_bin)
	a_fin = np.zeros(N_bin)
	e_fin = np.zeros(N_bin)
	with open(filename) as csvfile:
	        reader = csv.reader(csvfile, delimiter=',')
	        i = 0
	        for row in reader:
	        	a_ini[i] = float(row[0])
	        	e_ini[i] = float(row[1])
	        	a_fin[i] = float(row[2])
	        	e_fin[i] = float(row[3])
	        	i += 1
	N_a_ini = np.zeros(N_bins)
	N_a_fin = np.zeros(N_bins)
	for i in range(N_bin):
		if (a_ini[i] > a_min) and (a_ini[i] < a_max):
			j = int(np.floor(np.log(a_ini[i]/a_min)/dloga))
			N_a_ini[j] += 1
		if (a_fin[i] > a_min) and (a_fin[i] < a_max):
			j = int(np.floor(np.log(a_fin[i]/a_min)/dloga))
			N_a_fin[j] += 1
	if (plot_initial):
		plt.plot(a_bins/au, N_a_ini, label = 'Simulation, Initial', color='red', linestyle=plot_linestyle)
	plt.plot(a_bins/au, N_a_fin, label = plot_label, color=plot_color, linestyle=plot_linestyle)

def loadYCGData(filename, plot_label, plot_color, plot_linestyle, y_offset=90/N_bins):
	a = np.zeros(0)
	N = np.zeros(0)
	with open(filename) as csvfile:
	        reader = csv.reader(csvfile, delimiter=',')
	        for row in reader:
	        	a = np.append(a, float(row[0]))
	        	N = np.append(N, y_offset*float(row[1]))
	plt.plot(a, N, label = plot_label, color=plot_color, linestyle=plot_linestyle)

loadYCGData('YCGfig2_initial.csv', 'Yoo et al., Initial', plot_color='red', plot_linestyle='-')
loadYCGData('YCGfig2_10Msol.csv', r'Yoo et al., $10M_\odot$', plot_color='dodgerblue', plot_linestyle='-')
loadYCGData('YCGfig2_100Msol.csv', r'Yoo et al., $100M_\odot$', plot_color='forestgreen', plot_linestyle='-')
loadYCGData('YCGfig2_1000Msol.csv', r'Yoo et al., $1000M_\odot$', plot_color='darkorange', plot_linestyle='-')

loadData('binary_pop_YCG10Msol.csv', r'Simulation, $10M_\odot$', True, plot_color='dodgerblue', plot_linestyle='--')
loadData('binary_pop_YCG100Msol.csv', r'Simulation, $100M_\odot$', False, plot_color='forestgreen', plot_linestyle='--')
loadData('binary_pop_YCG1000Msol.csv', r'Simulation, $1000M_\odot$', False, plot_color='darkorange', plot_linestyle='--')

loadData('binary_pop_YCG10Msol_100closest.csv', r'Simulation 100 closest, $10M_\odot$', False, plot_color='dodgerblue', plot_linestyle='-.')
loadData('binary_pop_YCG100Msol_100closest.csv', r'Simulation 100 closest, $100M_\odot$', False, plot_color='forestgreen', plot_linestyle='-.')
loadData('binary_pop_YCG1000Msol_100closest.csv', r'Simulation 100 closest, $1000M_\odot$', False, plot_color='darkorange', plot_linestyle='-.')

loadData('binary_pop_YCG10Msol_100closest_maxwellian.csv', r'Simulation 100 closest Maxwellian, $10M_\odot$', False, plot_color='dodgerblue', plot_linestyle=':')
loadData('binary_pop_YCG100Msol_100closest_maxwellian.csv', r'Simulation 100 closest Maxwellian, $100M_\odot$', False, plot_color='forestgreen', plot_linestyle=':')
loadData('binary_pop_YCG1000Msol_100closest_maxwellian.csv', r'Simulation 100 closest Maxwellian, $1000M_\odot$', False, plot_color='darkorange', plot_linestyle=':')

#Move to centre of bins for plotting
a_bins *= np.exp(0.5*dloga)

ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('Semi-major axis, au')
plt.ylabel('Number of binaries')
plt.legend()
plt.show()                