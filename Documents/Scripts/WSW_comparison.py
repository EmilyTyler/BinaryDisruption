#Compare with WSW

import numpy as np
import os
os.system("python setup.py build_ext --inplace")

from matplotlib import pyplot as plt
from scipy.constants import au, giga, year, parsec, G

from encounters import calc_b_max, impulseEncounter, integrateEncounter, multipleImpulseEncounter
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
a = 10.0**4.0 * au
print('a = ', a)
#Orbital period       
P = 2.0 * np.pi * np.sqrt(a**3.0/(G*(m1+m2)))
#Initial energy
E_old = -G*(m1+m2)/(2.0*a)
print('E_old = ', E_old)
#Impact parameter
b_min = (np.pi*n_p*v_rms*(10.0*giga*year))**(-0.5)
b_max = calc_b_max(M_p, v_rms, a, m1, m2)
#WSW b_max
b_max_wsw = v_rms*P/(2.0*np.pi)
b = 0.05 * b_max_wsw
print('b = ', b)

#Number of encounters
N_enc = 10**4

#Average energy change from impulse
dE_imp_mean = 0.0
dE_imp_meansq = 0.0
#Number of negative energy changes
#N_neg = 0

for i in range(N_enc):
        notBound_new, a_new, e_new = impulseEncounter(m1, m2, v_rms, b, a, e, M_p)
        #print('a_new =', a_new)
        #if notBound_new:
                #print('BINARY BROKEN!')
        E_new = -G*(m1+m2)/(2.0*a_new)
        #print('E_new =', E_new)
        #print('E_new-E_old =', E_new-E_old)
        #if (E_new-E_old) < 0.0:
                #N_neg += 1
        dE_imp_mean += (E_new-E_old)
        dE_imp_meansq += (E_new-E_old)**2.0

#Normalise
dE_imp_mean /= N_enc
dE_imp_meansq /= N_enc
#Calculate variance
dE_imp_var = dE_imp_meansq - dE_imp_mean**2.0

'''
print('Impulse approximation code:')
print('Mean = ', dE_imp_mean)
print('Variance = ', dE_imp_var)



a_new = multipleImpulseEncounter(m1, m2, v_rms, b, a, e, M_p, N_enc)
E_new = -G*(m1+m2)/(2.0*a_new)
dE_imp_mean = np.sum(E_new - E_old)/N_enc
dE_imp_meansq = np.sum((E_new-E_old)**2.0)/N_enc
#Calculate variance
dE_imp_var = dE_imp_meansq - dE_imp_mean**2.0
'''
#print('Fraction of energy changes < 0 =', N_neg/N_enc)

'''
#Average energy change from three body
dE_thr_mean = 0.0
dE_thr_meansq = 0.0
for i in range(N_enc):
        notBound_new, a_new, e_new = integrateEncounter(m1, m2, v_rms, b, a, e, M_p)
        E_new = -G*(m1+m2)/(2.0*a_new)
        dE_thr_mean += (E_new-E_old)
        dE_thr_meansq += (E_new-E_old)**2.0
#Normalise
dE_thr_mean /= N_enc
dE_thr_meansq /= N_enc
#Calculate variance
dE_thr_var = dE_thr_meansq - dE_thr_mean**2.0
'''

#Calculate mean and variance from WSW
if a > b:
        #Single kick regime
        print('Single kick regime')
        dE_wsw_mean = 2.0*(G*M_p/(b*v_rms))**2.0
        dE_wsw_var = 4.0/3.0*(G*(m1+m2)/a)*(G*M_p/(b*v_rms))**2.0
else:
        #Tidal regime
        print('Tidal regime')
        dE_wsw_mean = 4.0/3.0*(G*M_p/(b*v_rms))**2.0 * (a/b)**2.0 * (1.0+3.0/2.0*e**2.0)
        dE_wsw_var = 4.0/5.0*(G*(m1+m2)/a)*(G*M_p/(b*v_rms))**2.0*(a/b)**2.0*(1.0 - e**2.0/3.0) + 16.0/45.0*(G*M_p/(b*v_rms))**4.0*(a/b)**4.0*(1.0 + 15.0*e**2.0)

#Crossing time over orbital period
#t_crossing_P = 2.0*b/(v_rms*P)
#print('Crossing time over orbital period = ', t_crossing_P)

print('Impulse approximation code:')
print('Mean = ', dE_imp_mean)
print('Variance = ', dE_imp_var)
'''
print('Integration code:')
print('Mean = ', dE_thr_mean)
print('Variance = ', dE_thr_var)
'''
print('Weinberg et al equations:')
print('Mean = ', dE_wsw_mean)
print('Variance = ', dE_wsw_var)

print('Fractional difference between impulse and Weinberg means:')
print((dE_imp_mean-dE_wsw_mean)/dE_wsw_mean)
print('Fractional difference between impulse and Weinberg variances:')
print((dE_imp_var-dE_wsw_var)/dE_wsw_var)












