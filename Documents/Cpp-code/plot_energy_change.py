#To plot the energy change against WSW average energy
import numpy as np
from matplotlib import pyplot as plt
import csv
from scipy.constants import au, G

plt.rc('font', family='serif')

#Parameters
a = 10.0**5.0 * au
e = 0.7
M_p = 3.0 * 2.0*10.0**30.0
v_rel = 220.0 * 1000.0
m1 = 2.0*10.0**30.0
m2 = 2.0*10.0**30.0


'''
#Impact parameter bins
N_b = 10
b_min = 10.0**(0.0)*au
print('b_min =', b_min)
b_max = 10.0**8.5*au
print('b_max =', b_max)
dlogb = (np.log(b_max)-np.log(b_min))/(N_b)
b_bins = np.array([b_min*np.exp(dlogb*i) for i in range(N_b)])
print('b_bins =', b_bins)
N_b_bins = np.zeros(N_b, dtype=int)

dE_sum = np.zeros(N_b)
dE2_sum = np.zeros(N_b)

#Import data function
def loadData(filename):
        with open(filename) as csvfile:
                reader = csv.reader(csvfile, delimiter=',')
                row_number = 0
                for row in reader:
                        E_ini = float(row[0])
                        E_fin = float(row[1])
                        b_star = float(row[2])
                        row_number += 1
                        j = int(np.floor(np.log(b_star/b_min)/dlogb))
                        if (b_min > b_star) or (b_star > b_max):
                                print('Out of range:', b_star)
                                print('b_min =', b_min)
                                print('b_max =', b_max)
                        print('Added part = ', E_fin-E_ini)
                        print('Average before addition = ', dE_sum[j])
                        dE_sum[j] += (E_fin - E_ini)
                        print('Average after addition  = ', dE_sum[j])
                        dE2_sum[j] += (E_fin - E_ini)**2.0
                        N_b_bins[j] += 1


#Choose which files to load
loadData('WSW_encounters_10e5.csv')


loadData('WSW_encounters_10e7.csv')
loadData('WSW_encounters_10e7_1.csv')
loadData('WSW_encounters_3x10e7.csv')
loadData('WSW_encounters_5x10e7.csv')


print('N_b_bins =', N_b_bins)
#Normalise
dE_sum /= N_b_bins
dE2_sum /= N_b_bins
#Mean
dE_mean = dE_sum
#Standard deviation
dE_mean_error = np.sqrt(dE2_sum - dE_mean**2.0)



#Move to the centre of the bins for calculations and plotting
b_bins *= np.exp(0.5*dlogb)

#Analytical prediction for average energy change (Weinberg et al.)
dE_avg_analytic = np.zeros(N_b)
for i in range(N_b):
        if b_bins[i] < a:
                dE_avg_analytic[i] = 2.0*(G*M_p/(b_bins[i]*v_rel))**2.0
        else:
                dE_avg_analytic[i] = 4.0/3.0 * (G*M_p/(b_bins[i]*v_rel))**2.0 * (a/b_bins[i])**2.0 * (1.0 + 3.0*e**2.0/2.0)
#Errors
dE_error_analytic = np.zeros(N_b)
for i in range(N_b):
        if b_bins[i] < a:
                variance = 4.0/3.0*G*(m1+m2)/a*(G*M_p/(b_bins[i]*v_rel))**2.0
        else:
                variance = 4.0/5.0*G*(m1+m2)/a*(G*M_p/(b_bins[i]*v_rel))**2.0*(a/b_bins[i])**2.0*(1.0 - e**2.0/3.0) + 16.0/45.0*(G*M_p/(b_bins[i]*v_rel))**4.0*(a/b_bins[i])**4.0*(1.0 +15.0*e**2.0)
        dE_error_analytic[i] = np.sqrt(variance)
#Change from reduced energy to total energy
dE_avg_analytic *= m1*m2/(m1+m2)
dE_error_analytic *= m1*m2/(m1+m2)


#Plot
plt.xlabel('Impact parameter, au')
ax = plt.gca()
ax.set_xscale('log')
plt.ylabel('Mean energy change due to an encounter, J')
ax.set_yscale('symlog')
plt.scatter(b_bins[np.where(dE_mean!=0.0)]/au, (dE_mean[np.where(dE_mean!=0.0)]), label='Impulse', marker='^')
plt.scatter(b_bins/au, dE_avg_analytic, label='Analytic', marker='x')
plt.axis([10.0**0.0, 10.0**8.5, -10.0**(42.0), 10.0**42.0])
plt.legend()
plt.title(r'$10^5$ encounters')
plt.show()

plt.xlabel('Impact parameter, au')
ax = plt.gca()
ax.set_xscale('log')
plt.ylabel('Standard deviation of energy change due to an encounter, J')
ax.set_yscale('symlog')
plt.scatter(b_bins[np.where(dE_mean_error!=0.0)]/au, (dE_mean_error[np.where(dE_mean_error!=0.0)]), label='Impulse', marker='^')
plt.scatter(b_bins/au, dE_error_analytic, label='Analytic', marker='x')
plt.axis([10.0**0.0, 10.0**8.5, 10.0**23.0, 10.0**39.0])
plt.legend()
plt.title(r'$10^5$ encounters')
plt.show()
'''



#Plot average against number of encounters
N_enc_min = 10**0
N_enc_max = 10**7
b = 10.0**4.0*au
with open("WSW_encounters_N_enc_log_b10e4au.csv") as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        row_number = 0
        for row in reader:
                dE_mean = float(row[0])
                std_dev = float(row[1])
                N_enc = float(row[2])
                row_number += 1
                plt.scatter(N_enc, dE_mean, marker='x', color='dodgerblue')
N_encs = np.linspace(N_enc_min, N_enc_max, num = 2)
dE_avg_analytic = np.zeros(2, dtype=float)
b_min = 0.9*b
b_max = 1.1*b
if b_min < a:
        dE_avg_analytic[0] = 2.0*(G*M_p/(b_min*v_rel))**2.0
else:
        dE_avg_analytic[0] = 4.0/3.0 * (G*M_p/(b_min*v_rel))**2.0 * (a/b_min)**2.0 * (1.0 + 3.0*e**2.0/2.0)
if b_max < a:
        dE_avg_analytic[1] = 2.0*(G*M_p/(b_max*v_rel))**2.0
else:
        dE_avg_analytic[1] = 4.0/3.0 * (G*M_p/(b_max*v_rel))**2.0 * (a/b_max)**2.0 * (1.0 + 3.0*e**2.0/2.0)
        #Change from reduced energy to total energy
dE_avg_analytic *= m1*m2/(m1+m2)
plt.plot(N_encs, [dE_avg_analytic[0]]*2, color='darkorange')
plt.plot(N_encs, [dE_avg_analytic[1]]*2, color='darkorange')
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('symlog')
plt.title(r'b = $10^{}$ au $\pm10\%$'.format(int(np.floor(np.log10(b/au)))))
plt.ylabel('Average energy change, J')
plt.xlabel('Number of encounters')
plt.show()


'''
#Plot distribution of energy changes
#Energy change bins
N_bins = 600
dE_min = 10.0**(22.0)
print('dE_min =', dE_min)
dE_max = 10.0**(33.0)
print('dE_max =', dE_max)
dlogdE = (np.log(dE_max)-np.log(dE_min))/(N_bins)
dE_bins = np.array([dE_min*np.exp(dlogdE*i) for i in range(N_bins)])
#print('dE_bins =', dE_bins)
N_dE = np.zeros((2,N_bins), dtype=float)
N_dE_v_dv = np.zeros((2,N_bins), dtype=float)
N_dE_dv_dv = np.zeros((2,N_bins), dtype=float)

filename = 'WSW_encounters_dE_b10e5au.csv'
with open(filename) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        row_number = 0
        for row in reader:
                dE = float(row[0])
                dE_v_dv = float(row[1])
                dE_dv_dv = float(row[2])
                row_number += 1

                j = int(np.floor(np.log(abs(dE)/dE_min)/dlogdE))
                if (dE_min > abs(dE)) or (abs(dE) > dE_max):
                        print('Out of range:', dE)
                        print('dE_min =', dE_min)
                        print('dE_max =', dE_max)
                if (dE > 0.0):
                        N_dE[0,j] += 1
                else:
                        N_dE[1,j] += 1

                j = int(np.floor(np.log(abs(dE_v_dv)/dE_min)/dlogdE))
                if (dE_min > abs(dE_v_dv)) or (abs(dE_v_dv) > dE_max):
                        print('Out of range:', dE_v_dv)
                        print('dE_min =', dE_min)
                        print('dE_max =', dE_max)
                if (dE_v_dv > 0.0):
                        N_dE_v_dv[0,j] += 1
                else:
                        N_dE_v_dv[1,j] += 1

                j = int(np.floor(np.log(abs(dE_dv_dv)/dE_min)/dlogdE))
                if (dE_min > abs(dE_dv_dv)) or (abs(dE_dv_dv) > dE_max):
                        print('Out of range:', dE_dv_dv)
                        print('dE_min =', dE_min)
                        print('dE_max =', dE_max)
                if (dE_dv_dv > 0.0):
                        N_dE_dv_dv[0,j] += 1
                else:
                        N_dE_dv_dv[1,j] += 1

#Move to the centre of the bins for calculations and plotting
dE_bins *= np.exp(0.5*dlogdE)

plt.plot(dE_bins, N_dE[0])
plt.plot(dE_bins, N_dE[1])
ax = plt.gca()
#ax.set_yscale('log')
ax.set_xscale('symlog', linthreshx=dE_min)
plt.xlabel('Energy change due to encounter, J')
plt.ylabel('Number of encounters')
plt.show()

plt.plot(dE_bins, N_dE_v_dv[0])
plt.plot(dE_bins, N_dE_v_dv[1])
ax = plt.gca()
#ax.set_yscale('log')
ax.set_xscale('symlog', linthreshx=dE_min)
plt.xlabel('Energy change due to encounter, vdv, J')
plt.ylabel('Number of encounters')
plt.show()

plt.plot(dE_bins, N_dE_dv_dv[0])
plt.plot(-dE_bins, N_dE_dv_dv[1])
ax = plt.gca()
#ax.set_yscale('log')
ax.set_xscale('symlog', linthreshx=dE_min)
plt.xlabel('Energy change due to encounter, dvdv, J')
plt.ylabel('Number of encounters')
plt.show()
'''