#Compare with WSW

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from matplotlib import pyplot as plt
from scipy.constants import au, giga, year, parsec, G

from encounters import calc_b_max, impulseEncounter, integrateEncounter
from frequency import calcFrequency


#Initialise variables
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 2.0*10.0**30.0
m2 = 2.0*10.0**30.0
#Density of dark matter halo solar masses/pc**3
rho = 0.009
#Convert to SI
rho = rho * 2.0*10.0**30.0/(parsec**3.0)
#Mass of perturbers
M_p = 1.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 220.0 * 1000.0
#Number density of perturbers
n_p = rho/M_p

#Semi-major axis
a = 10.0**6.0 * au
#Initial energy
E_old = -G*(m1+m2)/(2.0*a)
#Impact parameter
b_min = (np.pi*n_p*v_rms*(10.0*giga*year))**(-0.5)
b_max = calc_b_max(M_p, v_rms, a, m1, m2)
b = b_min

#Number of encounters
N_enc = 500

#Average energy change from impulse
dE_imp = np.zeros(N_enc)
for i in range(N_enc):
        notBound_new, a_new, e_new = impulseEncounter(m1, m2, v_rms, b, a, e, M_p)
        E_new = -G*(m1+m2)/(2.0*a_new)
        dE_imp[i] = (E_new-E_old)
#Put into frequency bins
E_bins_imp, N_imp, binwidth_imp = calcFrequency(dE_imp, N_enc//10)

#Average energy change from three body
dE_thr = np.zeros(N_enc)
for i in range(N_enc):
        notBound_new, a_new, e_new = integrateEncounter(m1, m2, v_rms, b, a, e, M_p)
        E_new = -G*(m1+m2)/(2.0*a_new)
        dE_thr[i] = (E_new-E_old)
#Put into frequency bins
E_bins_thr, N_thr, binwidth_thr = calcFrequency(dE_thr, N_enc//10)

#Calculate mean and dispersion from WSW
if a > b:
        #Single kick regime
        dE_wsw_mean = 2.0*(G*M_p/(b*v_rms))**2.0
        wsw_sigma_sq = 4.0/3.0*(G*(m1+m2)/a)*(G*M_p/(b*v_rms))**2.0
        wsw_sigma = np.sqrt(wsw_sigma_sq)
else:
        #Tidal regime
        dE_wsw_mean = 4.0/3.0*(G*M_p/(b*v_rms))**2.0 * (a/b)**2.0 * (1.0+3.0/2.0*e**2.0)
        wsw_sigma_sq = 4.0/5.0*(G*(m1+m2)/a)*(G*M_p/(b*v_rms))**2.0*(a/b)**2.0*(1.0 - e**2.0/3.0) + 16.0/45.0*(G*M_p/(b*v_rms))**4.0*(a/b)**4.0*(1.0 + 15.0*e**2.0)
        wsw_sigma = np.sqrt(wsw_sigma_sq)
print(dE_wsw_mean)
dE_wsw = np.random.normal(dE_wsw_mean, wsw_sigma, size=N_enc)
#Put into frequency bins
E_bins_wsw, N_wsw, binwidth_wsw = calcFrequency(dE_wsw, N_enc//10)

print(dE_imp[:100])
print(dE_thr[:100])
print(dE_wsw[:100])
print(N_imp)
print(E_bins_imp)
print(N_thr)
print(E_bins_thr)
print(N_wsw)
print(E_bins_wsw)
#Plot distributions
plt.title('Comparison of the average energy change during encounters with that predicted by Weinberg et al 1987')
plt.plot(E_bins_imp, N_imp/binwidth_imp/N_enc, label='Impulse code')
plt.plot(E_bins_thr, N_thr/binwidth_thr/N_enc, label='Three body code')
plt.plot(E_bins_wsw, N_wsw/binwidth_wsw/N_enc, label='WSW')
plt.xlabel('Average change in orbital energy per reduced mass after one encounter, J/kg')
plt.ylabel('Probability density of energy change, kg/J')
plt.legend()
plt.show()

























