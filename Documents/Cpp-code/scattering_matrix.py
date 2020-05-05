import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import au, parsec
import random
plt.rc('font', family='serif')

#MRA observational data
#The semi-major axis of the 25 most halo like wide binaries
a_MRA = au*np.array([340309, 54283, 5139, 22320, 685, 2805, 12155, 9135, 213, 537, 55410, 79139, 125073, 2198, 12380, 1176, 140, 468, 33416, 9439, 443, 672, 2986, 19173, 8565])
a_min = 10.0*au
a_max = np.max(a_MRA)
N_a_bins = 100
da = (np.log(a_max) - np.log(a_min))/(N_a_bins)
a_bins = np.array([a_min*np.exp(i*da) for i in range(N_a_bins)])

#Set up scattering matrix
#Probability that a binary with initial semi-major axis r_bins[i] will have final separation r_bins[j] is prob[i,j]
prob_a = np.zeros((N_a_bins, N_a_bins))
prob_r = np.zeros((N_a_bins, N_a_bins))
N_binaries = 0
with open("final_r_and_a_distributions_MRAparams_Mp10_Nbin10e5_format_ai_ri_ei_af_rf_ef.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        a_ini = float(row[0])
        r_ini = float(row[1])
        a_fin = float(row[3])
        r_fin = float(row[4])
        if (a_fin > 0.0):
            try:
                i = int(np.floor(np.log(a_ini/a_min)/da))
                j = int(np.floor(np.log(a_fin/a_min)/da))
                prob_a[i,j] += 1
                #i = int(np.floor(np.log(r_ini/a_min)/da))
                #j = int(np.floor(np.log(r_fin/a_min)/da))
                #prob_r[i,j] += 1
            except:
                #if ((r_fin < 10000.*parsec) and (a_fin < 10000.*parsec)):
                    #print('Not added: a_ini/pc =', a_ini/parsec, ', a_fin/pc=', a_fin/parsec, ", r_fin/pc=", r_fin/parsec)
                continue
            N_binaries += 1
#Normalise to make it a probability
#print(N_binaries)
for i in range(N_a_bins):
    N_binaries_with_a_ini_a_bins_i = np.sum(prob_a[i])
    prob_a[i, np.nonzero(prob_a[i])] /= N_binaries_with_a_ini_a_bins_i
    N_binaries_with_r_ini_a_bins_i = np.sum(prob_r[i])
    prob_r[i, np.nonzero(prob_r[i])] /= N_binaries_with_r_ini_a_bins_i

print('Scattering matrix generated with', N_binaries, 'binaries')

def generate_evolved_binaries_with_a_ini_i(i, N_bins=1):
    random_a_fins = np.zeros(N_bins)
    for j in range(N_bins):
        random_number = random.random()
        for k in range(N_a_bins):
            if np.sum(prob_a[i,0:k+1]) > random_number:
                random_a_fins[j] = a_bins[k]
                break
    return random_a_fins

def generate_evolved_binaries_from_initial_dist(a_0, a_1, alpha, N_bins):
    a_inis = np.zeros(N_bins)
    if (alpha == 1.0):
        c = np.log(a_0)/np.log(a_1/a_0)
        for i in range(N_bins):
            a_inis[i] = (a_1/a_0)**(random.random() + c)
    else:
        for i in range(N_bins):
            a_inis[i] = (random.random()*(a_1**(1.0-alpha) - a_0**(1.0-alpha)) + a_0**(1.0-alpha))**(1.0/(1.0 - alpha))
    a_fins = np.zeros(N_bins)
    for j in range(N_bins):
        i = int(np.floor(np.log(a_inis[j]/a_min)/da))
        a_fins[j] = generate_evolved_binaries_with_a_ini_i(i, 1)[0]
    return a_fins

#Generate population of evolved binaries
a_virtual_evolved = generate_evolved_binaries_from_initial_dist(10.0*au, 3.0*10.0**5.0*au, 1.0, 1000000)

#Bin virtual binaries
N_virtual_evolved = np.zeros(N_a_bins)
for i in range(np.size(a_virtual_evolved)):
    try:
        j = int(np.floor(np.log(a_virtual_evolved[i]/a_min)/da))
        N_virtual_evolved[j] += 1
    except:
        print(a_virtual_evolved[i])

plt.loglog(a_bins/au, N_virtual_evolved)
plt.xlim([10.0**3.0, 3.0*10.0**5.5])
plt.show()


        