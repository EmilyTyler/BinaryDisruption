#Test b and v distributions and N_enc for binning and monte carlo methods
import numpy as np
import os
os.system("python setup.py build_ext --inplace")

import matplotlib.pyplot as plt
from scipy.constants import mega, giga, year, parsec
from scipy.stats import maxwell
from monte_carlo import draw_b, vmaxwellianPdf, draw_maxwellian, draw_vmaxwellian, vmaxwellianComparison
from encounters import calc_b_max, encounterRate
from frequency import calcFrequency

#Initialise variables
#Semi-major axis, m
a = 0.1 * parsec
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 2.0*10.0**30.0
m2 = 2.0*10.0**30.0
#Density of dark matter halo solar masses/pc**3
rho = 0.008
#Convert to SI
rho = rho * 2.0*10.0**30.0/(parsec**3.0)
#Simulation time
T = 10.0**7.0 *year
#print('T =', T)
#Mass of perturbers
M_p = 3.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 220.0 * 1000.0
#Number density of perturbers
n_p = rho/M_p

#Minimum impact parameter
b_min = (np.pi*n_p*v_rms*T)**(-0.5)
#Maximum impact parameter
b_max = calc_b_max(M_p, v_rms, a, m1, m2)
#Minimum velocity
v_min = 10.0**(-2.0) * v_rms *0.0
#Maximum velocity
v_max = 10.0**2.0 * v_rms

N_v = 1000
N_b = 1000

#Monte Carlo
#Number of encounters
print('M Mean =', T*encounterRate(n_p, v_rms, b_min, b_max, v_min, v_max))
N_enc_MC = np.random.poisson(T*encounterRate(n_p, v_rms, b_min, b_max, v_min, v_max))
#b values
#b_MC = draw_b(b_max, N_enc_MC)
#v values
#v_MC = maxwell.rvs(scale=v_rms, size=N_enc_MC)
v_MC = draw_vmaxwellian(v_rms, v_min, v_max, N_enc_MC)
#Bin b values
#b_bins_MC, N_b_MC, db_MC = calcFrequency(b_MC, N_b)
#Bin v values
v_bins_MC, N_v_MC, dv_MC = calcFrequency(v_MC, N_v)

'''
#Binning
#b bins for encounter rate
dlogb = (np.log(b_max)-np.log(b_min))/N_b
b = np.array([b_min*np.exp(i*dlogb) for i in range(N_b)])
db_B = b * (np.exp(dlogb) - 1.0)
#print('b =', b)
#print('db_B =', db_B)


#db_B = (b_max - b_min)/(N_b)
#b = np.array([b_min + i*db_B for i in range(N_b)])

#v bins for encounter rate
dlogv = (np.log(v_max)-np.log(v_min))/N_v
v = np.array([v_min*np.exp(i*dlogv) for i in range(N_v)])
dv_B = v * (np.exp(dlogv) - 1.0)
#print('v =', v)
#print('dv_B =', dv_B)

#dv_B = (v_max - v_min)/(N_v)
#v = np.array([v_min + i*dv_B for i in range(N_v)])


#Bin MC values

N_b_MC = np.zeros(N_b, dtype=int)
for val in b_MC:
        if val == b_max:
                i = N_b - 1
        else:
                i = int(np.floor(np.log(val/b_min)/dlogb))
                N_b_MC[i] += 1

N_v_MC = np.zeros(N_v, dtype=int)
for val in v_MC:
        if val == v_max:
                i = N_v - 1
        else:
                i = int(np.floor(np.log(val/v_min)/dlogv))
                N_v_MC[i] += 1

'''
'''
#R[i,j] is the encounter rate for objects with impact parameter b[i] and relative velocity v[j]
R = np.zeros([N_b,N_v], dtype=float)
for i in range(N_b):
        for j in range(N_v):
                R[i,j] = encounterRate(n_p, v_rms, b[i], b[i]+db_B[i], v[j], v[j]+dv_B[j])
#print('R =', R)
#Time step
dt = 1.0/np.amax(R)
#print('dt =', dt)
#Number of timesteps
N_t = int(np.around(T/dt))
#Adjust dt to take rounding into account
dt = T/N_t
#print('dt =', dt)
#Sum of mean values
print('B Mean =', np.sum(R*T))
#Number of encounters matrix
N = np.rollaxis(np.array([[np.random.poisson(R[i,j]*dt, size=N_t) for j in range(N_v)] for i in range(N_b)]), 2)
#Total number of encounters
N_enc_B = np.sum(N)
#Array of indices where encounters happen
i_enc = np.transpose(np.array(np.nonzero(N)))
#b distribution
N_b_B = np.zeros(N_b, dtype=int)
#v distribution
N_v_B = np.zeros(N_v, dtype=int)
for (i,j,k) in i_enc:
        N_b_B[j] += N[i,j,k]
        N_v_B[k] += N[i,j,k]
b_bins_B = b
v_bins_B = v

#Print number of encounters
print('N_enc_MC = ', N_enc_MC)
print('N_enc_B = ', N_enc_B)

#Calculate binning means to plot
#b
RT_b = np.zeros(N_b)
#v
RT_v = np.zeros(N_v)
for i in range(N_b):
        for j in range(N_v):
                RT_b[i] += R[i,j]*T
                RT_v[j] += R[i,j]*T

#print(np.sum(RT_b))
#print(np.sum(RT_v))
'''
#Plot distributions
'''
#Plot b distributions
plt.plot(b_bins_MC/parsec, N_b_MC/N_enc_MC/(db_MC/parsec), label='Monte Carlo')
plt.plot(b_bins_B/parsec, N_b_B/N_enc_B/(db_B/parsec), label='Binning')
plt.plot(b_bins_B/parsec, RT_b/N_enc_B/(db_B/parsec), label='Binning means')
plt.xlabel('Impact parameter $b$, pc')
plt.ylabel(r'Probability density of encounters at impact parameter $b$, pc$^{-1}$')
plt.legend()
plt.show()
'''
#Plot v distributions
#Move to centre of bins
v_bins_MC += 0.5*dv_MC
plt.plot(v_bins_MC/1000.0, N_v_MC/N_enc_MC/(dv_MC/1000.0), label='Monte Carlo')
#plt.plot(v_bins_B/1000.0, N_v_B/N_enc_B/(dv_B/1000.0), label='Binning')
#plt.plot(v_bins_B/1000.0, RT_v/N_enc_B/(dv_B/1000.0), label='Binning means')
pdf = np.array([vmaxwellianPdf(x, v_rms) for x in v_bins_MC])
plt.plot(v_bins_MC/1000.0, pdf*1000.0, label='PDF MC')
comp = np.array([vmaxwellianComparison(x, v_rms, v_max, 3.0**0.5*v_rms) for x in v_bins_MC])
plt.plot(v_bins_MC/1000.0, comp*1000.0, label='Comp. MC')
#pdf = np.array([maxwellianPdf(x, v_rms) for x in v_bins_B])
#plt.plot(v_bins_B/1000.0, pdf*1000.0, label='PDF Binning')
plt.xlabel(r'Relative velocity $v$, kms$^{-1}$')
plt.ylabel(r'Probability density of encounters at relative velocity $v$, km$^{-1}$s')
plt.legend()
plt.show()
'''
#Plot b distributions
plt.plot(b_bins_B/parsec, N_b_MC/N_enc_MC, label='Monte Carlo')
plt.plot(b_bins_B/parsec, N_b_B/N_enc_B, label='Binning')
plt.plot(b_bins_B/parsec, RT_b/N_enc_B, label='Binning means')
plt.xlabel('Impact parameter $b$, pc')
plt.ylabel(r'Fraction of encounters at impact parameter $b$')
plt.legend()
plt.show()

#Plot v distributions
plt.plot(v_bins_B/1000.0, N_v_MC/N_enc_MC, label='Monte Carlo')
plt.plot(v_bins_B/1000.0, N_v_B/N_enc_B, label='Binning')
plt.plot(v_bins_B/1000.0, RT_v/N_enc_B, label='Binning means')
plt.xlabel(r'Relative velocity $v$, kms$^{-1}$')
plt.ylabel(r'Fraction of encounters at relative velocity $v$')
plt.legend()
plt.show()

'''













































