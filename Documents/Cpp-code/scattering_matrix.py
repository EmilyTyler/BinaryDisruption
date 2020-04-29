import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import au, parsec
plt.rc('font', family='serif')

#Set up scattering matrix
a_min = 10.*au
a_max = 500.*au
N_a_bins = 10
da = (np.log(a_max) - np.log(a_min))/(N_a_bins)
a_bins = np.array([a_min*np.exp(i*da) for i in range(N_a_bins)])
#Probability that a binary with initial semi-major axis r_bins[i] will have final separation r_bins[j] is prob[i,j]
prob = np.zeros((N_a_bins, N_a_bins))
N_binaries = 0
with open("final_semimajoraxis_distribution_MRAparams_Mp10_Nbin10e6.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        a_ini = float(row[0])
        a_fin = float(row[2])
        if (a_fin > 0.0):
            try:
                i = int(np.floor(np.log(a_ini/a_min)/da))
                j = int(np.floor(np.log(a_fin/a_min)/da))
                prob[i,j] += 1
            except:
                if ((a_fin > 0.0) and (a_fin < 10000.*parsec)):
                    print('Not added: a_ini/au =', a_ini/au, ', a_fin/au=', a_fin/au)
                continue
            N_binaries += 1
#Normalise to make it a probability
print(N_binaries)
for i in range(N_a_bins):
    N_binaries_with_a_ini_a_bins_i = np.sum(prob[i])
    prob[i, np.nonzero(prob[i])] /= N_binaries_with_a_ini_a_bins_i

print('a_bins/au =', a_bins/au)
print(prob)
