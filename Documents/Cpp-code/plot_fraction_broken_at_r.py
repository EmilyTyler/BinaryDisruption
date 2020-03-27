import csv
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.constants import au, parsec, giga, year

#Plot fraction of binaries that broke at each separation
r_min = 10**(-1.0)*parsec
r_max = 10.0**(5.0)*parsec
N_r_bins = 100
dr = (np.log(r_max) - np.log(r_min))/(N_r_bins)
r_bins = np.array([r_min*np.exp(i*dr) for i in range(N_r_bins)])
dr_log = np.zeros(N_r_bins)
for i in range(N_r_bins-1):
	dr_log[i] = r_bins[i+1]-r_bins[i]
dr_log[N_r_bins-1] = r_max - r_bins[N_r_bins-1]

total_at_r = np.zeros(N_r_bins)
broke_at_r = np.zeros(N_r_bins)
bound_at_r = np.zeros(N_r_bins)


with open("separation_over_time_and_report_broken_binaries_Mp1Msol_Nbin10e4_a_i1pc.csv") as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	for row in reader:
		r = float(row[0])
		broke = int(row[1])
		try:
			j = int(np.floor(np.log(r/r_min)/dr))
		except:
			j=0
		total_at_r[j] += 1
		if (broke==1):
			broke_at_r[j] += 1
		else:
			bound_at_r[j] += 1


#Find fraction broken at each r value
fraction_broken_at_r = np.zeros(N_r_bins)
for i in range(N_r_bins):
	if (total_at_r[i] != 0.0):
		fraction_broken_at_r[i] = broke_at_r[i]/total_at_r[i]
	else:
		print(broke_at_r[i])

plt.semilogx(r_bins/parsec, fraction_broken_at_r, color='dodgerblue')
plt.ylabel(r'Fraction of times binaries with separation $r$ broke with separation $r$', wrap=True, color='dodgerblue')
ax1 = plt.gca()
ax2 = ax1.twinx()
ax2.semilogx(r_bins/parsec, bound_at_r, color='darkorange')
plt.xlabel(r'$r$/pc')
ax2.set_ylabel(r'Number of times a bound binary has separation $r$', wrap=True, color='darkorange')
plt.show()