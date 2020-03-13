import matplotlib.pyplot as plt
import csv
import numpy as np
from scipy.constants import au
plt.rc('font', family='serif')

a_min = 10.0**(-2.0) * au
a_max = 10.0**15.0 * au
N_bins = 40

dloga = (np.log(a_max)-np.log(a_min))/(N_bins)
a_bins = np.array([a_min*np.exp(i*dloga) for i in range(N_bins)])

N_bin = 10**6

def loadData(filename, plot_label, plot_initial, plot_color, plot_linestyle, yoffset=1.0):
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
	print(np.sum(N_a_fin))
	if (plot_initial):
		plt.plot(a_bins/au, N_a_ini*yoffset, label = 'Initial distribution', color=plot_color, linestyle=':')
	plt.plot(a_bins/au, N_a_fin*yoffset, label = plot_label, color=plot_color, linestyle=plot_linestyle)

def loadYCGData(filename, plot_label, plot_color, plot_linestyle, y_offset=90/N_bins):
	a = np.zeros(0)
	N = np.zeros(0)
	with open(filename) as csvfile:
	        reader = csv.reader(csvfile, delimiter=',')
	        for row in reader:
	        	a = np.append(a, float(row[0]))
	        	N = np.append(N, y_offset*float(row[1]))
	plt.plot(a, N, label = plot_label, color=plot_color, linestyle=plot_linestyle)

def loadMRAData(filename, plot_label, plot_color, plot_linestyle, y_offset=1720/N_bins):
	a = np.zeros(0)
	N = np.zeros(0)
	with open(filename) as csvfile:
	        reader = csv.reader(csvfile, delimiter=',')
	        for row in reader:
	        	a = np.append(a, float(row[0]))
	        	N = np.append(N, y_offset*float(row[1]))
	plt.plot(a, N, label = plot_label, color=plot_color, linestyle=plot_linestyle)

#loadYCGData('YCGfig2_initial.csv', 'Yoo et al., Initial', plot_color='black', plot_linestyle=':')
#loadYCGData('YCGfig2_10Msol.csv', r'Yoo et al., $10M_\odot$, fitted', plot_color='black', plot_linestyle='--')
#loadYCGData('YCGfig2_100Msol.csv', r'Yoo et al., $100M_\odot$, fitted', plot_color='black', plot_linestyle='-.')
#loadYCGData('YCGfig2_1000Msol.csv', r'Yoo et al., $1000M_\odot$, fitted', plot_color='black', plot_linestyle='-')

#loadYCGData('YCGfig2_10Msol_points.csv', r'Yoo et al., $10M_\odot$', plot_color='grey', plot_linestyle='--')
#loadYCGData('YCGfig2_100Msol_points.csv', r'Yoo et al., $100M_\odot$', plot_color='grey', plot_linestyle='-.')
#loadYCGData('YCGfig2_1000Msol_points.csv', r'Yoo et al., $1000M_\odot$', plot_color='grey', plot_linestyle='-')

#loadMRAData('MRAfig1_10Msol.csv', r'Monroy-Rodr$\mathrm{\'{\i}}$guez & Allen, $10M_\odot$', plot_color='darkorange', plot_linestyle='--')
#loadMRAData('MRAfig1_100Msol.csv', r'Monroy-Rodr$\mathrm{\'{\i}}$guez & Allen, $100M_\odot$', plot_color='darkorange', plot_linestyle='-.')
#loadMRAData('MRAfig1_1000Msol.csv', r'Monroy-Rodr$\mathrm{\'{\i}}$guez & Allen, $1000M_\odot$', plot_color='darkorange', plot_linestyle='-')

#loadData('binary_pop_YCG10Msol.csv', r'Our simulation, $10M_\odot$', True, plot_color='dodgerblue', plot_linestyle='--')
#loadData('binary_pop_YCG100Msol.csv', r'Our simulation, $100M_\odot$', False, plot_color='dodgerblue', plot_linestyle='-.')
#loadData('binary_pop_YCG1000Msol.csv', r'Our simulation, $1000M_\odot$', False, plot_color='dodgerblue', plot_linestyle='-')

#loadData('binary_pop_YCG10Msol_YCGbmax_Nencclosest.csv', r'Simulation, All encounters, Yoo et al. $b_\mathrm{max}$, $10M_\odot$', False, plot_color='forestgreen', plot_linestyle='--')
#loadData('binary_pop_YCG100Msol_YCGbmax_Nencclosest.csv', r'Simulation, All encounters, Yoo et al. $b_\mathrm{max}$, $100M_\odot$', False, plot_color='forestgreen', plot_linestyle='-.')
#loadData('binary_pop_YCG1000Msol_YCGbmax_Nencclosest.csv', r'Full simulation with closest encounters code, Yoo et al. $b_\mathrm{max}$, $1000M_\odot$', False, plot_color='forestgreen', plot_linestyle='-')

#loadData('binary_pop_YCG10Msol_YCGbmax_100closest.csv', r'Simulation 100 closest, $10M_\odot$', False, plot_color='forestgreen', plot_linestyle='--')
#loadData('binary_pop_YCG100Msol_YCGbmax_100closest.csv', r'Simulation 100 closest, $100M_\odot$', False, plot_color='forestgreen', plot_linestyle='-.')
#loadData('binary_pop_YCG1000Msol_YCGbmax_100closest.csv', r'Simulation, 100 closest, Yoo et al. $b_\mathrm{max}$, $1000M_\odot$', False, plot_color='forestgreen', plot_linestyle='-')

#loadData('binary_pop_YCG10Msol_YCGbmax_100closest_maxwellian.csv', r'Simulation 100 closest Maxwellian, $10M_\odot$', False, plot_color='darkorange', plot_linestyle='--')
#loadData('binary_pop_YCG100Msol_YCGbmax_100closest_maxwellian.csv', r'Simulation 100 closest Maxwellian, $100M_\odot$', False, plot_color='darkorange', plot_linestyle='-.')
#loadData('binary_pop_YCG1000Msol_YCGbmax_100closest_maxwellian.csv', r'Simulation 100 closest Maxwellian, $1000M_\odot$', False, plot_color='darkorange', plot_linestyle='-')

#loadData('binary_pop_YCG10Msol_YCGbmax_1000closest.csv', r'Simulation 1000 closest, $10M_\odot$', False, plot_color='m', plot_linestyle='--')
#loadData('binary_pop_YCG100Msol_YCGbmax_1000closest.csv', r'Simulation 1000 closest, $100M_\odot$', False, plot_color='m', plot_linestyle='-.')
#loadData('binary_pop_YCG1000Msol_YCGbmax_1000closest.csv', r'Simulation 1000 closest, $1000M_\odot$', False, plot_color='m', plot_linestyle='-')

#loadData('binary_pop_YCG10Msol_YCGbmax_1000closest_maxwellian.csv', r'Simulation 1000 closest Maxwellian, $10M_\odot$', False, plot_color='red', plot_linestyle='--')
#loadData('binary_pop_YCG100Msol_YCGbmax_1000closest_maxwellian.csv', r'Simulation 1000 closest Maxwellian, $100M_\odot$', False, plot_color='red', plot_linestyle='-.')
#loadData('binary_pop_YCG1000Msol_YCGbmax_1000closest_maxwellian.csv', r'Simulation 1000 closest Maxwellian, $1000M_\odot$', False, plot_color='red', plot_linestyle='-')

#loadData('binary_pop_YCG10Msol_YCGbmax_10000closest.csv', r'Simulation 10000 closest, $10M_\odot$', False, plot_color='saddlebrown', plot_linestyle='--')
#loadData('binary_pop_YCG100Msol_YCGbmax_10000closest.csv', r'Simulation 10000 closest, $100M_\odot$', False, plot_color='saddlebrown', plot_linestyle='-.')
#loadData('binary_pop_YCG1000Msol_YCGbmax_10000closest.csv', r'Simulation 10000 closest, $1000M_\odot$', False, plot_color='saddlebrown', plot_linestyle='-')

#loadData('binary_pop_YCG10Msol_YCGbmax_10000closest_maxwellian.csv', r'Simulation 10000 closest Maxwellian, $10M_\odot$', False, plot_color='mediumblue', plot_linestyle='--')
#loadData('binary_pop_YCG100Msol_YCGbmax_10000closest_maxwellian.csv', r'Simulation 10000 closest Maxwellian, $100M_\odot$', False, plot_color='mediumblue', plot_linestyle='-.')
#loadData('binary_pop_YCG1000Msol_YCGbmax_10000closest_maxwellian.csv', r'Simulation 10000 closest Maxwellian, $1000M_\odot$', False, plot_color='mediumblue', plot_linestyle='-')

#loadData('binary_pop_YCG10Msol_mybmax_Nencclosest.csv', r'Full simulation with closest encounters code, my $b_\mathrm{max}$, $10M_\odot$', False, plot_color='darkorange', plot_linestyle='--')
#loadData('binary_pop_YCG100Msol_mybmax_Nencclosest.csv', r'Full simulation with closest encounters code, my $b_\mathrm{max}$, $100M_\odot$', False, plot_color='darkorange', plot_linestyle='-.')
#loadData('binary_pop_YCG1000Msol_mybmax_Nencclosest.csv', r'Full simulation with closest encounters code, my $b_\mathrm{max}$, $1000M_\odot$', False, plot_color='darkorange', plot_linestyle='-')

#172 km/s
'''
loadData('binary_pop_YCG10Msol_172kms.csv', r'My simulation $\sigma_\mathrm{rel} = 172\mathrm{kms}^{-1}$, $10M_\odot$', True, plot_color='dodgerblue', plot_linestyle='--')
loadData('binary_pop_YCG100Msol_172kms.csv', r'My simulation $\sigma_\mathrm{rel} = 172\mathrm{kms}^{-1}$, $100M_\odot$', False, plot_color='dodgerblue', plot_linestyle='-.')
loadData('binary_pop_YCG1000Msol_172kms.csv', r'My simulation $\sigma_\mathrm{rel} = 172\mathrm{kms}^{-1}$, $1000M_\odot$', False, plot_color='dodgerblue', plot_linestyle='-')

loadData('binary_pop_YCG10Msol_YCGbmax_100closest_172kms.csv', r'Simulation 100 closest $\sigma_\mathrm{rel} = 172\mathrm{kms}^{-1}$, $10M_\odot$', False, plot_color='forestgreen', plot_linestyle='--')
loadData('binary_pop_YCG100Msol_YCGbmax_100closest_172kms.csv', r'Simulation 100 closest $\sigma_\mathrm{rel} = 172\mathrm{kms}^{-1}$, $100M_\odot$', False, plot_color='forestgreen', plot_linestyle='-.')
loadData('binary_pop_YCG1000Msol_YCGbmax_100closest_172kms.csv', r'Simulation 100 closest $\sigma_\mathrm{rel} = 172\mathrm{kms}^{-1}$, $1000M_\odot$', False, plot_color='forestgreen', plot_linestyle='-')

loadData('binary_pop_YCG10Msol_YCGbmax_100closest_maxwellian_172kms.csv', r'Simulation 100 closest $\sigma_\mathrm{rel} = 172\mathrm{kms}^{-1}$, Maxwellian, $10M_\odot$', False, plot_color='darkorange', plot_linestyle='--')
loadData('binary_pop_YCG100Msol_YCGbmax_100closest_maxwellian_172kms.csv', r'Simulation 100 closest $\sigma_\mathrm{rel} = 172\mathrm{kms}^{-1}$, Maxwellian, $100M_\odot$', False, plot_color='darkorange', plot_linestyle='-.')
loadData('binary_pop_YCG1000Msol_YCGbmax_100closest_maxwellian_172kms.csv', r'Simulation 100 closest $\sigma_\mathrm{rel} = 172\mathrm{kms}^{-1}$, Maxwellian, $1000M_\odot$', False, plot_color='darkorange', plot_linestyle='-')
'''
'''
#all encs expect 1000msol, maxwellian, 200kms, YCG bmax
loadData('binary_pop_YCG10Msol_YCGbmax_allencs_maxwellian.csv', r'Simulation, All encounters, Yoo et al. $b_\mathrm{max}$, Maxwellian, $10M_\odot$', False, plot_color='darkorange', plot_linestyle='--')
loadData('binary_pop_YCG100Msol_YCGbmax_allencs_maxwellian.csv', r'Simulation, All encounters, Yoo et al. $b_\mathrm{max}$, Maxwellian, $100M_\odot$', False, plot_color='darkorange', plot_linestyle='-.')
loadData('binary_pop_YCG1000Msol_YCGbmax_100closest_maxwellian.csv', r'Simulation, 100 closest, Yoo et al. $b_\mathrm{max}$, Maxwellian, $1000M_\odot$', False, plot_color='darkorange', plot_linestyle='-')
'''

#bmax convergence testing
'''
loadData('binary_pop_10Msol_0_2bmax.csv', r'0.2 $b_\mathrm{max}$', True, plot_color='darkorange', plot_linestyle='-')
loadData('binary_pop_10Msol_0_5bmax.csv', r'0.5 $b_\mathrm{max}$', False, plot_color='forestgreen', plot_linestyle='-')
loadData('binary_pop_10Msol_1bmax.csv', r'$b_\mathrm{max}$', False, plot_color='darkviolet', plot_linestyle='-')
loadData('binary_pop_10Msol_2bmax.csv', r'2 $b_\mathrm{max}$', False, plot_color='red', plot_linestyle='-')
loadData('binary_pop_10Msol_5bmax.csv', r'5 $b_\mathrm{max}$', False, plot_color='mediumblue', plot_linestyle='-')
'''

#New 100 closest plots
#loadData('binary_pop_100closest_1000Msol.csv', r'Our simulation, 100 closest, $1000M_\odot$', False, plot_color='darkorange', plot_linestyle='-')

#New plots, including rebound binaries
loadData('final_semi_major_axes_distribution_1Msol_initial_log_dist_rebound_included_Nbin10e6.csv', r'Rebound included, $1M_\odot$', True, plot_color='dodgerblue', plot_linestyle='--')
loadData('final_semi_major_axes_distribution_10Msol_initial_log_dist_rebound_included_Nbin10e6.csv', r'Rebound included, $10M_\odot$', False, plot_color='dodgerblue', plot_linestyle='-.')
loadData('final_semi_major_axes_distribution_100Msol_initial_log_dist_rebound_included_Nbin10e6.csv', r'Rebound included, $100M_\odot$', False, plot_color='dodgerblue', plot_linestyle='-')
loadData('final_semi_major_axes_distribution_1Msol_initial_log_dist_rebound_not_included_Nbin10e6.csv', r'Rebound not included, $1M_\odot$', False, plot_color='darkorange', plot_linestyle='--')
loadData('final_semi_major_axes_distribution_10Msol_initial_log_dist_rebound_not_included_Nbin10e6.csv', r'Rebound not included, $10M_\odot$', False, plot_color='darkorange', plot_linestyle='-.')
loadData('final_semi_major_axes_distribution_100Msol_initial_log_dist_rebound_not_included_Nbin10e6.csv', r'Rebound not included, $100M_\odot$', False, plot_color='darkorange', plot_linestyle='-')



#Move to centre of bins for plotting
#a_bins *= np.exp(0.5*dloga)

ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('Semi-major axis, au')
plt.xlim(a_min/au, a_max/au)
plt.ylim(6.0*10.0**(-1.0), 3.0*10.0**5.0)
plt.ylabel('Number of binaries')
plt.legend()
plt.show()

#Plot separations
loadData('final_separation_distribution_1Msol_initial_log_dist_rebound_included_Nbin10e6.csv', r'Rebound included, $1M_\odot$', False, plot_color='dodgerblue', plot_linestyle='--')
loadData('final_separation_distribution_10Msol_initial_log_dist_rebound_included_Nbin10e6.csv', r'Rebound included, $10M_\odot$', False, plot_color='dodgerblue', plot_linestyle='-.')
loadData('final_separation_distribution_100Msol_initial_log_dist_rebound_included_Nbin10e6.csv', r'Rebound included, $100M_\odot$', False, plot_color='dodgerblue', plot_linestyle='-')
loadData('final_separation_distribution_1Msol_initial_log_dist_rebound_not_included_Nbin10e6.csv', r'Rebound not included, $1M_\odot$', False, plot_color='darkorange', plot_linestyle='--')
loadData('final_separation_distribution_10Msol_initial_log_dist_rebound_not_included_Nbin10e6.csv', r'Rebound not included, $10M_\odot$', False, plot_color='darkorange', plot_linestyle='-.')
loadData('final_separation_distribution_100Msol_initial_log_dist_rebound_not_included_Nbin10e6.csv', r'Rebound not included, $100M_\odot$', False, plot_color='darkorange', plot_linestyle='-')
#With J+T params
loadData('final_separation_distribution_J+Tparams_initial_J+Tlog_dist_rebound_included_Nbin10e6.csv', r'Rebound included, J+T params', False, plot_color='forestgreen', plot_linestyle='-')
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('Separation, au')
plt.xlim(a_min/au, a_max/au)
plt.ylim(6.0*10.0**(-1.0), 3.0*10.0**5.0)
plt.ylabel('Number of binaries')
plt.legend()
plt.show()                