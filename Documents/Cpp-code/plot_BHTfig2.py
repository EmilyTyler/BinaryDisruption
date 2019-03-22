import matplotlib.pyplot as plt
import csv
import numpy as np
from scipy.constants import giga, year, G, parsec
plt.rc('font', family='serif')

'''
#Plot simulations of 25 binaries
for i in range(59):
        BHTsim = np.array([[]])
        with open('BHTfig2_mysim_25bin_{}.csv'.format(i)) as csvfile:
                reader = csv.reader(csvfile, delimiter=',')
                for row in reader:
                        BHTsim = np.append(BHTsim, row)              
        BHTsim = BHTsim.astype(np.float)
        BHTsim = np.reshape(BHTsim, (-1,2))
        plt.plot(BHTsim[:,0]/(giga*year), BHTsim[:,1], color='darkgrey')
'''
#Simulation of 1000 binaries
BHTsim = np.array([[]])
with open('BHTfig2_mysim_1000bin_0.csv') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
                BHTsim = np.append(BHTsim, row)              
BHTsim = BHTsim.astype(np.float)
BHTsim = np.reshape(BHTsim, (-1,2))
plt.plot(BHTsim[:,0]/(giga*year), BHTsim[:,1], label=r'My simulation, V$_\mathrm{rel}\times$ Maxwellian', color='red')

#With 2*b_max
'''
#Plot simulations of 25 binaries
for i in range(28):
        BHTsim = np.array([[]])
        with open('BHTfig2_mysim_2bmax_25bin_{}.csv'.format(i)) as csvfile:
                reader = csv.reader(csvfile, delimiter=',')
                for row in reader:
                        BHTsim = np.append(BHTsim, row)              
        BHTsim = BHTsim.astype(np.float)
        BHTsim = np.reshape(BHTsim, (-1,2))
        plt.plot(BHTsim[:,0]/(giga*year), BHTsim[:,1], color='darkgrey')
'''
#Simulation of 1000 binaries
BHTsim = np.array([[]])
with open('BHTfig2_mysim_2bmax_1000bin_0.csv') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
                BHTsim = np.append(BHTsim, row)              
BHTsim = BHTsim.astype(np.float)
BHTsim = np.reshape(BHTsim, (-1,2))
plt.plot(BHTsim[:,0]/(giga*year), BHTsim[:,1], label=r'My simulation, up to $2\times b_\mathrm{max}$', color='darkorange')


#With 5*b_max
'''
#Plot simulations of 25 binaries
for i in range(28):
        BHTsim = np.array([[]])
        with open('BHTfig2_mysim_5bmax_25bin_{}.csv'.format(i)) as csvfile:
                reader = csv.reader(csvfile, delimiter=',')
                for row in reader:
                        BHTsim = np.append(BHTsim, row)              
        BHTsim = BHTsim.astype(np.float)
        BHTsim = np.reshape(BHTsim, (-1,2))
        plt.plot(BHTsim[:,0]/(giga*year), BHTsim[:,1], color='darkgrey')
'''
#Simulation of 1000 binaries
BHTsim = np.array([[]])
with open('BHTfig2_mysim_5bmax_1000bin_0.csv') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
                BHTsim = np.append(BHTsim, row)              
BHTsim = BHTsim.astype(np.float)
BHTsim = np.reshape(BHTsim, (-1,2))
plt.plot(BHTsim[:,0]/(giga*year), BHTsim[:,1], label=r'My simulation, up to $5\times b_\mathrm{max}$', color='darkviolet')


#Maxwellian rather than vmaxwellian
'''
#Plot simulations of 25 binaries
for i in range(4):
        BHTsim = np.array([[]])
        with open('BHTfig2_mysim_maxwellian_25bin_{}.csv'.format(i)) as csvfile:
                reader = csv.reader(csvfile, delimiter=',')
                for row in reader:
                        BHTsim = np.append(BHTsim, row)              
        BHTsim = BHTsim.astype(np.float)
        BHTsim = np.reshape(BHTsim, (-1,2))
        plt.plot(BHTsim[:,0]/(giga*year), BHTsim[:,1], color='darkgrey')

BHTsim = np.array([[]])
with open('BHTfig2_mysim_maxwellian_1000bin_0.csv') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
                BHTsim = np.append(BHTsim, row)              
BHTsim = BHTsim.astype(np.float)
BHTsim = np.reshape(BHTsim, (-1,2))
plt.plot(BHTsim[:,0]/(giga*year), BHTsim[:,1], label=r'My simulation, Maxwellian', color='darkviolet')
'''
'''
#Plot exponential decay
M_b = 2.0 * 2.0*10.0**30.0
V_rms = 10.0**5.0
rho = 0.1 *2.0*10.0**30.0/(parsec**3.0)
M_p = 3.0 * 2.0*10.0**30.0
a = 0.1 * parsec
#t_half = 0.194*10.0**10.0*year
#t_half = np.log(2.0)/40.0*(3.0/(2.0*np.pi))**0.5*M_b*V_rms/(G*rho*M_p*a)
#t_half = 5.14*10.0**9.0*year*2.0*10.0**30.0/M_p + 3.3*10.0**8.0*year
t_half = 3.52*10.0**9.0*year
print('t_half, yr =', t_half/year)
rate = np.log(2.0)/t_half
t = np.linspace(0, 10.0*giga*year, num=100)
y = np.exp(-rate*t)
plt.plot(t/(giga*year), y, label='Exponential decay', color='darkorange')
'''

#Plot BHT data
BHTsim = np.array([[]])
with open('BHTfig2simdata.csv') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
                BHTsim = np.append(BHTsim, row)              
BHTsim = BHTsim.astype(np.float)
BHTsim = np.reshape(BHTsim, (-1,2))
BHTsim = BHTsim[BHTsim[:,0].argsort()]
plt.plot(BHTsim[:,0], BHTsim[:,1], label='Bahcall et al. simulation', color='forestgreen')
BHTexp = np.array([[]])
with open('BHTfig2expdecay.csv') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
                BHTexp = np.append(BHTexp, row)              
BHTexp = BHTexp.astype(np.float)
BHTexp = np.reshape(BHTexp, (-1,2))
BHTexp = BHTexp[BHTexp[:,0].argsort()]
plt.plot(BHTexp[:,0], BHTexp[:,1], label='Bahcall et al. exponential', color='dodgerblue')

#Setup axes
plt.legend()
plt.xlabel('Time, Gyr')
plt.ylabel('Survival probability')
plt.show()