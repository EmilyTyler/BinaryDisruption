#To plot binary_pop results to look like YCG fig 2

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

import matplotlib.pyplot as plt
from scipy.constants import au, year, mega, giga, parsec, kilo
from frequency import calcFrequency
import csv

plt.rc('font', family='serif')

#Number of bins
N_bins = 20

def read_plot_csv(filename, legendlabel, linestyle, color, y_offset=1.0):
        data = np.array([[]])
        with open(filename) as csvfile:
                reader = csv.reader(csvfile, delimiter=',')
                for row in reader:
                        data = np.append(data, row)                
        data = data.astype(np.float)
        data = np.reshape(data, (-1,2))
        plt.loglog(data[:,0], data[:,1]*y_offset, label=legendlabel, linestyle=linestyle, color=color)
        return None

def read_bin_plot_npz(filename, legendlabel, linestyle, color, y_offset=1.0):
        data = np.load(filename)
        a, e, N_broken = data['a_fin'], data['e_fin'], data['N_broken']
        a_bins, N_a, a_binwidth = calcFrequency(a, N_bins, log=True)
        da = a_bins *(np.exp(a_binwidth)-1.0)
        plt.loglog(a_bins/au, N_a/(da/au)*y_offset, label=legendlabel, linestyle=linestyle, color=color)
     

#Plot YCG
read_plot_csv('YCGfig2_initial.csv', 'Yoo et al: Initial', 'dashed', 'red', y_offset=1.0/600.0)
read_plot_csv('YCGfig2_10Msol.csv', r'Yoo et al: $M_p = 10M_\odot$', 'dashed', 'royalblue', y_offset=1.0/600.0)
read_plot_csv('YCGfig2_100Msol.csv', r'Yoo et al: $M_p = 100M_\odot$', 'dashed', 'darkorange', y_offset=1.0/600.0)
read_plot_csv('YCGfig2_1000Msol.csv', r'Yoo et al: $M_p = 1000M_\odot$', 'dashed', 'forestgreen', y_offset=1.0/600.0)

#Plot filter star frame
read_bin_plot_npz('simulation_filter_10Msol_1e6.npz', r'$M_p = 10M_\odot$', '-', 'royalblue', y_offset=1.0)
read_bin_plot_npz('simulation_filter_100Msol_1e6.npz', r'$M_p = 100M_\odot$', '-', 'darkorange', y_offset=1.0)
read_bin_plot_npz('simulation_filter_1000Msol_1e6.npz', r'$M_p = 1000M_\odot$', '-', 'forestgreen', y_offset=1.0)
'''
#Plot filter binary frame
read_bin_plot_npz('simulation_filter_binaryframe_10Msol_1e6.npz', r'Binary frame: $M_p = 10M_\odot$', 'dashdot', 'royalblue', y_offset=1.0)
read_bin_plot_npz('simulation_filter_binaryframe_100Msol_1e6.npz', r'Binary frame: $M_p = 100M_\odot$', 'dashdot', 'darkorange', y_offset=1.0)
read_bin_plot_npz('simulation_filter_binaryframe_1000Msol_1e6.npz', r'Binary frame: $M_p = 1000M_\odot$', 'dashdot', 'forestgreen', y_offset=1.0)
'''
'''
#Plot filter larger tidal radius
read_bin_plot_npz('simulation_filter_a_T8e5au_10Msol_1e5.npz', r'$a_T=8e5\mathrm{au}: M_p = 10M_\odot$', 'dotted', 'royalblue', y_offset=10.0)
read_bin_plot_npz('simulation_filter_a_T8e5au_100Msol_1e5.npz', r'$a_T=8e5\mathrm{au}: M_p = 100M_\odot$', 'dotted', 'darkorange', y_offset=10.0)
read_bin_plot_npz('simulation_filter_a_T8e5au_1000Msol_1e5.npz', r'$a_T=8e5\mathrm{au}: M_p = 1000M_\odot$', 'dotted', 'forestgreen', y_offset=10.0)
'''
'''
#Plot filter wrong random angle
read_bin_plot_npz('simulation_filter_wrongangle_10Msol_1e6.npz', r'Random $\theta$: $M_p = 10M_\odot$', 'dotted', 'royalblue', y_offset=1.0)
read_bin_plot_npz('simulation_filter_wrongangle_100Msol_1e6.npz', r'Random $\theta$: $M_p = 100M_\odot$', 'dotted', 'darkorange', y_offset=1.0)
read_bin_plot_npz('simulation_filter_wrongangle_1000Msol_1e6.npz', r'Random $\theta$: $M_p = 1000M_\odot$', 'dotted', 'forestgreen', y_offset=1.0)
'''
'''
#Plot 100 closest encounters
read_bin_plot_npz('simulation_100closest_10Msol_1e6.npz', r'100 closest: $M_p = 10M_\odot$', 'dotted', 'royalblue', y_offset=1.0)
read_bin_plot_npz('simulation_100closest_100Msol_1e6.npz', r'100 closest: $M_p = 100M_\odot$', 'dotted', 'darkorange', y_offset=1.0)
read_bin_plot_npz('simulation_100closest_1000Msol_1e6.npz', r'100 closest: $M_p = 1000M_\odot$', 'dotted', 'forestgreen', y_offset=1.0)
'''
#Plot GPU results
read_bin_plot_npz('simulation_GPU_10Msol_1e6.npz', r'GPU: $M_p = 10M_\odot$', 'dotted', 'royalblue', y_offset=1.0)
read_bin_plot_npz('simulation_GPU_100Msol_1e6.npz', r'GPU: $M_p = 100M_\odot$', 'dotted', 'darkorange', y_offset=1.0)
read_bin_plot_npz('simulation_GPU_1000Msol_1e6.npz', r'GPU: $M_p = 1000M_\odot$', 'dotted', 'forestgreen', y_offset=1.0)

plt.legend()
plt.xlim(10.0**3.0, 3.0*10.0**5.0)
plt.ylim(10.0**(-3.0), 10.0**1.0)
plt.xlabel('Semi-major axis, au')
plt.ylabel(r'Probability density, au$^{-1}$')
#plt.title('Simulations: filter method')
plt.show()

'''
#Jiang and Tremaine, no diffusive regime
#Load data
data_JTnd_10Msol = np.load('JT_nodiff_10Msol_1e6.npz')
a_JTnd_10Msol, e_JTnd_10Msol, N_broken_JTnd_10Msol = data_JTnd_10Msol['a_fin'], data_JTnd_10Msol['e_fin'], data_JTnd_10Msol['N_broken']
data_JTnd_100Msol = np.load('JT_nodiff_100Msol_1e6.npz')
a_JTnd_100Msol, e_JTnd_100Msol, N_broken_JTnd_100Msol = data_JTnd_100Msol['a_fin'], data_JTnd_100Msol['e_fin'], data_JTnd_100Msol['N_broken']
data_JTnd_1000Msol = np.load('JT_nodiff_1000Msol_1e6.npz')
a_JTnd_1000Msol, e_JTnd_1000Msol, N_broken_JTnd_1000Msol = data_JTnd_1000Msol['a_fin'], data_JTnd_1000Msol['e_fin'], data_JTnd_1000Msol['N_broken']
#Bin data
a_bins_JTnd_10Msol, N_a_JTnd_10Msol, a_binwidth_JTnd_10Msol = calcFrequency(a_JTnd_10Msol, N_bins, log=True)
da_JTnd_10Msol = a_bins_JTnd_10Msol *(np.exp(a_binwidth_JTnd_10Msol)-1.0)
a_bins_JTnd_100Msol, N_a_JTnd_100Msol, a_binwidth_JTnd_100Msol = calcFrequency(a_JTnd_100Msol, N_bins, log=True)
da_JTnd_100Msol = a_bins_JTnd_100Msol *(np.exp(a_binwidth_JTnd_100Msol)-1.0)
a_bins_JTnd_1000Msol, N_a_JTnd_1000Msol, a_binwidth_JTnd_1000Msol = calcFrequency(a_JTnd_1000Msol, N_bins, log=True)
da_JTnd_1000Msol = a_bins_JTnd_1000Msol *(np.exp(a_binwidth_JTnd_1000Msol)-1.0)
#Plot
plt.loglog(a_bins_JTnd_10Msol/au, N_a_JTnd_10Msol/(da_JTnd_10Msol/au), label=r'$M_p = 10M_\odot$')
plt.loglog(a_bins_JTnd_100Msol/au, N_a_JTnd_100Msol/(da_JTnd_100Msol/au), label=r'$M_p = 100M_\odot$')
plt.loglog(a_bins_JTnd_1000Msol/au, N_a_JTnd_1000Msol/(da_JTnd_1000Msol/au), label=r'$M_p = 1000M_\odot$')
plt.legend()
plt.xlim(10.0**3.0, 3.0*10.0**5.0)
plt.ylim(10.0**(-3.0), 10.0**1.0)
plt.xlabel('Semi-major axis, au')
plt.ylabel(r'Number density, au$^{-1}$')
plt.title('Jiang and Tremaine without diffusive regime')
plt.show()

#Numerical solution to Fokker-Planck
#Load data
data_numerical = np.load('numerical.npz')
N_num, a_bins_num, da_num = data_numerical['N'], data_numerical['a_bins'], data_numerical['da']
M_p = np.array([10.0, 100.0, 1000.0])*2.0*10.0**30.0
#Plot
plt.loglog(a_bins_num/au, N_num[0,0]/(da_num/au), label='Initial')
N_tau = np.size(N_num[0,:,0])
for i in range(3):
        plt.loglog(a_bins_num/au, N_num[i,N_tau-1]/(da_num/au), label=r'$M={}M_\odot$'.format(int(M_p[i]/(2.0*10.0**30.0))))
plt.xlabel('Semi-major axis, au')
plt.ylabel(r'Number of binaries, au$^{-1}$')
plt.title('Numerical solution to Fokker-Planck', wrap=True)
plt.axis([1000.0, 300000.0, 10.0**(-3.0), 10.0])
plt.legend()
plt.show()

#Weinberg et al
#Load data
data_WSW_10Msol = np.load('WSW_10Msol_1e6.npz')
a_WSW_10Msol, e_WSW_10Msol, N_broken_WSW_10Msol = data_WSW_10Msol['a_fin'], data_WSW_10Msol['e_fin'], data_WSW_10Msol['N_broken']
data_WSW_100Msol = np.load('WSW_100Msol_1e6.npz')
a_WSW_100Msol, e_WSW_100Msol, N_broken_WSW_100Msol = data_WSW_100Msol['a_fin'], data_WSW_100Msol['e_fin'], data_WSW_100Msol['N_broken']
data_WSW_1000Msol = np.load('WSW_1000Msol_1e6.npz')
a_WSW_1000Msol, e_WSW_1000Msol, N_broken_WSW_1000Msol = data_WSW_1000Msol['a_fin'], data_WSW_1000Msol['e_fin'], data_WSW_1000Msol['N_broken']
#Bin data
a_bins_WSW_10Msol, N_a_WSW_10Msol, a_binwidth_WSW_10Msol = calcFrequency(a_WSW_10Msol, N_bins, log=True)
da_WSW_10Msol = a_bins_WSW_10Msol *(np.exp(a_binwidth_WSW_10Msol)-1.0)
a_bins_WSW_100Msol, N_a_WSW_100Msol, a_binwidth_WSW_100Msol = calcFrequency(a_WSW_100Msol, N_bins, log=True)
da_WSW_100Msol = a_bins_WSW_100Msol *(np.exp(a_binwidth_WSW_100Msol)-1.0)
a_bins_WSW_1000Msol, N_a_WSW_1000Msol, a_binwidth_WSW_1000Msol = calcFrequency(a_WSW_1000Msol, N_bins, log=True)
da_WSW_1000Msol = a_bins_WSW_1000Msol *(np.exp(a_binwidth_WSW_1000Msol)-1.0)
#Plot
plt.loglog(a_bins_WSW_10Msol/au, N_a_WSW_10Msol/(da_WSW_10Msol/au), label=r'$M_p = 10M_\odot$')
plt.loglog(a_bins_WSW_100Msol/au, N_a_WSW_100Msol/(da_WSW_100Msol/au), label=r'$M_p = 100M_\odot$')
plt.loglog(a_bins_WSW_1000Msol/au, N_a_WSW_1000Msol/(da_WSW_1000Msol/au), label=r'$M_p = 1000M_\odot$')
plt.legend()
plt.xlim(10.0**3.0, 3.0*10.0**5.0)
plt.ylim(10.0**(-3.0), 10.0**1.0)
plt.xlabel('Semi-major axis, au')
plt.ylabel(r'Number density, au$^{-1}$')
plt.title('Weinberg et al')
plt.show()
'''





















