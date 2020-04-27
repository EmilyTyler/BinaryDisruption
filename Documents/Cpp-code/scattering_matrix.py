import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import au, parsec
plt.rc('font', family='serif')

#Set up scattering matrix
r_min = 10.*au
r_max = 10000.*parsec
N_r_bins = 10
dr = (np.log(r_max) - np.log(r_min))/(N_r_bins)
r_bins = np.array([r_min*np.exp(i*dr) for i in range(N_r_bins)])
#Probability that a binary with initial semi-major axis r_bins[i] will have final separation r_bins[j] is prob[i,j]
prob = np.zeros((N_r_bins, N_r_bins))
N_binaries = 0
with open("final_semimajoraxis_distribution_MRAparams_Mp1000_Nbin10e6.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        a_ini = float(row[0])
        r_fin = float(row[2])
        try:
            i = int(np.floor(np.log(a_ini/r_min)/dr))
            j = int(np.floor(np.log(r_fin/r_min)/dr))
            prob[i,j] += 1
        except:
            if ((r_fin > 0.0) and (r_fin < 10000.*parsec)):
                print('Not added: a_ini/au =', a_ini/au, ', r_fin/au=', r_fin/au)
            continue
        N_binaries += 1
#Normalise to make it a probability
print(N_binaries)
prob[np.nonzero(prob)] /= N_binaries

print(prob)
