import matplotlib.pyplot as plt
import csv
import numpy as np
from scipy.constants import au, year, giga
plt.rc('font', family='serif')

'''
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
'''

'''
N_enc_min = 10**0
N_enc_max = 10**4
N_bins = 20

#dN_enc = (N_enc_max-N_enc_min)/(N_bins - 1)
#N_enc_bins = np.array([N_enc_min + i*dN_enc for i in range(N_bins)])
dN_enc = np.log(N_enc_max/N_enc_min)/(N_bins - 1)
N_enc_bins = np.array([N_enc_min*np.exp(i*dN_enc) for i in range(N_bins)])

N_N_enc = np.zeros(N_bins, dtype=int)

def loadData(filename, plot_label, plot_color):
	N_bin = 0
	with open(filename) as csvfile:
	        reader = csv.reader(csvfile, delimiter=',')
	        for row in reader:
	        	N_bin += 1
	        	N_enc = float(row[0])
	        	#j = int(np.floor((N_enc - N_enc_min)/dN_enc))
	        	j = int(np.floor(np.log(N_enc/N_enc_min)/dN_enc))
	        	N_N_enc[j] += 1
	plt.plot(N_enc_bins + 0.5*dN_enc, N_N_enc/N_bin, label = plot_label, color=plot_color)

loadData('N_enc_broken_dist_10Msol_a10e4au.csv', plot_label=r'$M_p=10M_\odot$', plot_color='red')
loadData('N_enc_broken_dist_100Msol_a10e4au.csv', plot_label=r'$M_p=100M_\odot$', plot_color='forestgreen')
loadData('N_enc_broken_dist_1000Msol_a10e4au.csv', plot_label=r'$M_p=1000M_\odot$', plot_color='dodgerblue')

plt.xlabel('Number of encounters required to break binary')
plt.ylabel('Fraction of broken binaries')
ax = plt.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
plt.legend()
plt.show()
'''
'''
t_min = 0.0
t_max = 10.0**10.0 * year
N_bins = 20

dt = (t_max-t_min)/(N_bins)
t_bins = np.array([t_min + i*dt for i in range(N_bins)])

N_t = np.zeros(N_bins, dtype=int)

def loadData(filename, plot_label, plot_color):
	N_bin = 0
	with open(filename) as csvfile:
	        reader = csv.reader(csvfile, delimiter=',')
	        for row in reader:
	        	t = float(row[1])
	        	if (t<=t_max):
	        		N_bin += 1
	        		j = int(np.floor((t - t_min)/dt))
	        		N_t[j] += 1
	plt.plot((t_bins + 0.5*dt)/(10.0**9*year), N_t, label = plot_label, color=plot_color)

loadData('t_broken_dist_10Msol_a10e4au.csv', plot_label=r'$M_p=10M_\odot$', plot_color='red')
loadData('t_broken_dist_100Msol_a10e4au.csv', plot_label=r'$M_p=100M_\odot$', plot_color='forestgreen')
loadData('t_broken_dist_1000Msol_a10e4au.csv', plot_label=r'$M_p=1000M_\odot$', plot_color='dodgerblue')

plt.xlabel('Disruption time, Gyr')
plt.ylabel(r'Number of binaries (out of $10^5$)')
ax = plt.gca()
plt.xlim([0,10])
plt.legend()
plt.show()
'''

'''
t_min = 0.0
t_max = 10.0**10.0 * year
N_bins = 200

dt = (t_max-t_min)/(N_bins)
t_bins = np.array([t_min + i*dt for i in range(N_bins)])

def loadData(filename, plot_label, plot_color):
	N_bin = 0
	N_t = np.zeros(N_bins, dtype=int)
	with open(filename) as csvfile:
	        reader = csv.reader(csvfile, delimiter=',')
	        for row in reader:
	        	t = float(row[1])
	        	if (t<=t_max):
	        		N_bin += 1
	        		j = int(np.floor((t - t_min)/dt))
	        		N_t[j] += 1
	#Make it cumulative
	for i in range(1, N_bins):
		N_t[i] += N_t[i-1]
	plt.plot((t_bins + 0.5*dt)/(10.0**9*year), N_t, label = plot_label, color=plot_color)

loadData('t_broken_dist_10Msol_a10e4au.csv', plot_label=r'$M_p=10M_\odot$', plot_color='red')
loadData('t_broken_dist_100Msol_a10e4au.csv', plot_label=r'$M_p=100M_\odot$', plot_color='forestgreen')
loadData('t_broken_dist_1000Msol_a10e4au.csv', plot_label=r'$M_p=1000M_\odot$', plot_color='dodgerblue')

plt.xlabel('Time, Gyr')
plt.ylabel(r'Number of binaries disrupted (out of $10^5$)')
ax = plt.gca()
#ax.set_yscale('log')
plt.xlim([0,10])
plt.legend()
plt.show()
'''


'''
t_min = 0.0
t_max = 2.0*10.0**10.0 * year
N_bins = 200

dt = (t_max-t_min)/(N_bins)
t_bins = np.array([t_min + i*dt for i in range(N_bins)])

def loadData(filename, plot_label, plot_color):
	N_bin = 0
	N_t = np.zeros(N_bins, dtype=int)
	with open(filename) as csvfile:
	        reader = csv.reader(csvfile, delimiter=',')
	        for row in reader:
	        	t = float(row[1])
	        	if (t<=t_max):
	        		N_bin += 1
	        		j = int(np.floor((t - t_min)/dt))
	        		N_t[j] += 1
	#Make it cumulative
	for i in range(1, N_bins):
		N_t[i] += N_t[i-1]
	plt.plot((t_bins + 0.5*dt)/(10.0**9*year), N_t, label = plot_label, color=plot_color)

loadData('t_broken_dist_10Msol_a10e4au_t20Gyr.csv', plot_label=r'$M_p=10M_\odot$', plot_color='red')
loadData('t_broken_dist_100Msol_a10e4au_t20Gyr.csv', plot_label=r'$M_p=100M_\odot$', plot_color='forestgreen')
loadData('t_broken_dist_1000Msol_a10e4au_t20Gyr.csv', plot_label=r'$M_p=1000M_\odot$', plot_color='dodgerblue')

plt.xlabel('Time, Gyr')
plt.ylabel(r'Number of binaries disrupted (out of $10^5$)')
ax = plt.gca()
#ax.set_yscale('log')
plt.xlim([0,20])
plt.legend()
plt.show()
'''


'''
t_min = 0.0
t_max = 2.0*10.0**10.0 * year
N_bins = 20

dt = (t_max-t_min)/(N_bins)
t_bins = np.array([t_min + i*dt for i in range(N_bins)])

def loadData(filename, plot_label, plot_color):
	N_bin = 0
	N_t = np.zeros(N_bins, dtype=float)
	with open(filename) as csvfile:
	        reader = csv.reader(csvfile, delimiter=',')
	        for row in reader:
	        	t = float(row[1])
	        	if (t<=t_max):
	        		N_bin += 1
	        		j = int(np.floor((t - t_min)/dt))
	        		N_t[j] += 1
	print(N_t)
	N_tc = np.zeros(N_bins)
	N_tc[0] = N_t[0]
	for i in range(1, N_bins):
		N_tc[i] = N_t[i] + N_tc[i-1]
	print(N_tc)
	N_t[0] = N_t[0]/(10.0**5.0)
	for i in range(1, N_bins):
		N_t[i] /= (10.0**5.0 - N_tc[i-1])
	print(N_t)
	plt.plot((t_bins + 0.5*dt)/(10.0**9*year), N_t, label = plot_label, color=plot_color)

loadData('t_broken_dist_10Msol_a10e4au_t20Gyr.csv', plot_label=r'$M_p=10M_\odot$', plot_color='red')
loadData('t_broken_dist_100Msol_a10e4au_t20Gyr.csv', plot_label=r'$M_p=100M_\odot$', plot_color='forestgreen')
loadData('t_broken_dist_1000Msol_a10e4au_t20Gyr.csv', plot_label=r'$M_p=1000M_\odot$', plot_color='dodgerblue')

plt.xlabel('Time, Gyr')
plt.ylabel(r'Fraction of surviving binaries disrupted')
ax = plt.gca()
#ax.set_yscale('log')
plt.xlim([0,20])
plt.legend()
plt.show()
'''

t_min = 0.0
t_max = 2.0*10.0**10.0 * year
N_bins = 200

dt = (t_max-t_min)/(N_bins)
t_bins = np.array([t_min + i*dt for i in range(N_bins)])

def loadData(filename, plot_label, plot_color):
	N_bin = 0
	N_t = np.zeros(N_bins, dtype=int)
	with open(filename) as csvfile:
	        reader = csv.reader(csvfile, delimiter=',')
	        for row in reader:
	        	t = float(row[1])
	        	if (t<=t_max):
	        		N_bin += 1
	        		j = int(np.floor((t - t_min)/dt))
	        		N_t[j] += 1
	#Make it cumulative
	for i in range(1, N_bins):
		N_t[i] += N_t[i-1]
	plt.plot((t_bins + 0.5*dt)/(10.0**9*year), N_t, label = plot_label, color=plot_color)

loadData('t_broken_dist_10Msol_a10e4au_t20Gyr.csv', plot_label=r'$M_p=10M_\odot$', plot_color='red')
loadData('t_broken_dist_100Msol_a10e4au_t20Gyr.csv', plot_label=r'$M_p=100M_\odot$', plot_color='forestgreen')
loadData('t_broken_dist_1000Msol_a10e4au_t20Gyr.csv', plot_label=r'$M_p=1000M_\odot$', plot_color='dodgerblue')
loadData('t_broken_dist_10Msol_a10e4au.csv', plot_label=r'$M_p=10M_\odot$', plot_color='darkorange')
loadData('t_broken_dist_100Msol_a10e4au.csv', plot_label=r'$M_p=100M_\odot$', plot_color='darkviolet')
loadData('t_broken_dist_1000Msol_a10e4au.csv', plot_label=r'$M_p=1000M_\odot$', plot_color='saddlebrown')

plt.xlabel('Time, Gyr')
plt.ylabel(r'Number of binaries disrupted (out of $10^5$)')
ax = plt.gca()
#ax.set_yscale('log')
plt.xlim([0,20])
plt.legend()
plt.grid(linestyle='--')
plt.show()