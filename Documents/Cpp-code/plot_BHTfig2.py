import matplotlib.pyplot as plt
import csv
import numpy as np
from scipy.constants import giga, year
plt.rc('font', family='serif')

#Plot simulations of 25 binaries
for i in range(26):
        BHTsim = np.array([[]])
        with open('BHTfig2_mysim_25bin_{}.csv'.format(i)) as csvfile:
                reader = csv.reader(csvfile, delimiter=',')
                for row in reader:
                        BHTsim = np.append(BHTsim, row)              
        BHTsim = BHTsim.astype(np.float)
        BHTsim = np.reshape(BHTsim, (-1,2))
        plt.plot(BHTsim[:,0]/(giga*year), BHTsim[:,1], color='darkgrey')

#Simulation of 1000 binaries
BHTsim = np.array([[]])
with open('BHTfig2_mysim_1000bin_0.csv') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
                BHTsim = np.append(BHTsim, row)              
BHTsim = BHTsim.astype(np.float)
BHTsim = np.reshape(BHTsim, (-1,2))
plt.plot(BHTsim[:,0]/(giga*year), BHTsim[:,1], label='My simulation', color='red')


#Plot BHT data
BHTsim = np.array([[]])
with open('BHTfig2simdata.csv') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
                BHTsim = np.append(BHTsim, row)              
BHTsim = BHTsim.astype(np.float)
BHTsim = np.reshape(BHTsim, (-1,2))
BHTsim = BHTsim[BHTsim[:,0].argsort()]
plt.plot(BHTsim[:,0], BHTsim[:,1], label='BHT sim.', color='forestgreen')
BHTexp = np.array([[]])
with open('BHTfig2expdecay.csv') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
                BHTexp = np.append(BHTexp, row)              
BHTexp = BHTexp.astype(np.float)
BHTexp = np.reshape(BHTexp, (-1,2))
BHTexp = BHTexp[BHTexp[:,0].argsort()]
plt.plot(BHTexp[:,0], BHTexp[:,1], label='BHT exp.', color='dodgerblue')

#Setup axes
plt.legend()
plt.xlabel('Time, Gyr')
plt.ylabel('Survival probability')
plt.show()