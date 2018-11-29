import matplotlib.pyplot as plt
import csv
import numpy as np
from scipy.constants import giga, year
plt.rc('font', family='serif')
'''
#Simulations of 25 binaries
for i in range(0,100):
        N_broken_sim_25 = np.load('BHTfig2_25bin_mysim_{}.npz'.format(i))['N_broken']
        t_sim_25 = np.load('BHTfig2_25bin_mysim_{}.npz'.format(i))['t']
        plt.plot(t_sim_25/(giga*year), 1-N_broken_sim_25, color='darkgrey')
      
'''
'''
#Encounters with dv given by equations in BHT
#Simulations of 25 binaries
for i in range(0,100):
        N_broken_sim_25 = np.load('BHTfig2_25bin_mysim_BHTenc_{}.npz'.format(i))['N_broken']
        t_sim_25 = np.load('BHTfig2_25bin_mysim_BHTenc_{}.npz'.format(i))['t']
        plt.plot(t_sim_25/(giga*year), 1-N_broken_sim_25, color='darkgrey')
 '''     
        
#My simulation data
N_broken_sim = np.load('BHTfig2_1000bin_mysim_0.npz')['N_broken']
t_sim = np.load('BHTfig2_1000bin_mysim_0.npz')['t']
plt.plot(t_sim/(giga*year), 1-N_broken_sim, label='My sim, 1000 binaries', color='red')

#My simulation t method
N_broken_sim = np.load('BHTfig2_1000bin_mysim_tmethod_0.npz')['N_broken']
t_sim = np.load('BHTfig2_1000bin_my sim_tmethod_0.npz')['t']
plt.plot(t_sim/(giga*year), 1-N_broken_sim, label='My sim, t method, 1000 binaries', color='gold')


#Encounters with dv given by equations in BHT   
#My simulation data
N_broken_sim = np.load('BHTfig2_1000bin_mysim_BHTenc_0.npz')['N_broken']
t_sim = np.load('BHTfig2_1000bin_mysim_BHTenc_0.npz')['t']
plt.plot(t_sim/(giga*year), 1-N_broken_sim, label='BHT encounters, 1000 binaries', color='darkorange')

#Encounters with dv given by double equations in BHT   
#My simulation data
N_broken_sim = np.load('BHTfig2_1000bin_2BHT_0.npz')['N_broken']
t_sim = np.load('BHTfig2_1000bin_2BHT_0.npz')['t']
plt.plot(t_sim/(giga*year), 1-N_broken_sim, label='2x BHT encounters, 1000 binaries', color='violet')
'''
#Encounters with change in relative star speed given by equations in BHT   
N_broken_sim = np.load('BHTfig2_1000bin_BHTenc_changev_tmethod_0.npz')['N_broken']
t_sim = np.load('BHTfig2_1000bin_BHTenc_changev_tmethod_0.npz')['t']
plt.plot(t_sim/(giga*year), 1-N_broken_sim, label='BHT new, 1000 binaries', color='gold')
'''
'''
#Different separations
#Average separation
N_broken_sim = np.load('BHTfig2_100bin_mysim_avgsep_0.npz')['N_broken']
t_sim = np.load('BHTfig2_100bin_mysim_avgsep_0.npz')['t']
plt.plot(t_sim/(giga*year), 1-N_broken_sim, label='Avg. sep., 100 binaries', color='darkorange')

#Maximum separation
N_broken_sim = np.load('BHTfig2_100bin_mysim_maxsep_0.npz')['N_broken']
t_sim = np.load('BHTfig2_100bin_mysim_maxsep_0.npz')['t']
plt.plot(t_sim/(giga*year), 1-N_broken_sim, label='Max. sep., 100 binaries', color='brown')

#Minimum separation
N_broken_sim = np.load('BHTfig2_100bin_mysim_minsep_0.npz')['N_broken']
t_sim = np.load('BHTfig2_100bin_mysim_minsep_0.npz')['t']
plt.plot(t_sim/(giga*year), 1-N_broken_sim, label='Min. sep., 100 binaries', color='violet')
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
plt.plot(BHTsim[:,0],BHTsim[:,1],label='BHT sim.', color='forestgreen')
BHTexp = np.array([[]])
with open('BHTfig2expdecay.csv') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
                BHTexp = np.append(BHTexp, row)              
BHTexp = BHTexp.astype(np.float)
BHTexp = np.reshape(BHTexp, (-1,2))
BHTexp = BHTexp[BHTexp[:,0].argsort()]
plt.plot(BHTexp[:,0],BHTexp[:,1],label='BHT exp.', color='dodgerblue')

#Setup axes
plt.legend()
plt.xlabel('Time, Gyr')
plt.ylabel('Survival probability')
plt.show()