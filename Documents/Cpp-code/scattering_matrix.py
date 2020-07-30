import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import au, parsec
import random
from scipy import stats
plt.rc('font', family='serif')
random.seed(10)

a_min = 30.0*au
a_max = 20000.0*parsec
N_a_bins = 40
da = (np.log(a_max) - np.log(a_min))/(N_a_bins)
a_bins = np.array([a_min*np.exp(i*da) for i in range(N_a_bins)])

r_cutoff = 10000.*parsec

#Initial distribution
alpha = 1.0
print('alpha =', alpha)

M_ps = np.array([10, 100, 1000])
N_Mps = np.size(M_ps)
linestyles = np.array(['dotted', 'dashed', 'dashdot'])

#Set up scattering matrix
print('Importing simulation data and generating scattering matrices')
#Probability that a binary with initial semi-major axis r_bins[i] will have final separation r_bins[j] is prob[i,j] given that it never exceeds the cut-off semi-major axis
prob_a = np.zeros((N_Mps, N_a_bins, N_a_bins))
prob_a_cutoff = np.zeros((N_Mps, N_a_bins))
N_input_dist = np.zeros((N_Mps, N_a_bins))
N_input_dist_r_proj = np.zeros((N_Mps, N_a_bins))
N_input_initial = np.zeros((N_Mps, N_a_bins))
prob_r = np.zeros((N_Mps, N_a_bins, N_a_bins))
prob_r_projected = np.zeros((N_Mps, N_a_bins, N_a_bins))
N_broken_cutoff = np.zeros(N_Mps, dtype=int)
N_binaries = np.zeros(N_Mps, dtype=int)
for k in range(N_Mps):
    with open("final_r_and_a_distributions_MRAparams_Mp{}_Nbin10e6_format_ai_ri_ei_af_rf_ef.csv".format(M_ps[k])) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            a_ini = float(row[0])
            r_ini = float(row[1])
            a_fin = float(row[3])
            r_fin = float(row[4])
            i = int(np.floor(np.log(a_ini/a_min)/da))
            N_input_initial[k,i] += 1
            if (a_fin > 0.0):
            #if True:
                try:
                    j = int(np.floor(np.log(a_fin/a_min)/da))
                    prob_a[k,i,j] += 1
                    N_input_dist[k,j] += 1
                    if (a_ini==a_fin):
                        r_fin = r_ini
                    #Find final r projected
                    u = random.uniform(-1.0, 1.0)
                    theta = random.uniform(0.0, 2.0*np.pi)
                    rand_vec = np.array([np.sqrt(1.0-u**2.0)*np.cos(theta), np.sqrt(1.0-u**2.0)*np.sin(theta), u])
                    rand_vec /= np.linalg.norm(rand_vec)
                    r_fin_projected = r_fin*np.sqrt(1.0-rand_vec[0]**2.0)
                    #
                    j_r = int(np.floor(np.log(r_fin/a_min)/da))
                    j_r_proj = int(np.floor(np.log(r_fin_projected/a_min)/da))
                    prob_r[k,i,j_r] += 1
                    prob_r_projected[k,i,j_r_proj] += 1
                    N_input_dist_r_proj[k,j_r_proj] += 1
                except:
                    if ((r_fin < r_cutoff) and (a_fin < r_cutoff)):
                        print(r_fin_projected/parsec)
                        print(r_ini/parsec, r_fin/parsec)
                        print(a_ini/parsec, a_fin/parsec)
                        input()
                    else:
                        N_broken_cutoff[k] += 1
                        prob_a_cutoff[k,i] += 1
                        continue
                N_binaries[k] += 1
            else:
                N_broken_cutoff[k] += 1
                prob_a_cutoff[k,i] += 1
    #Normalise to make it a probability
    for i in range(N_a_bins):
        prob_a[k,i, np.nonzero(prob_a[k,i])] /= np.sum(prob_a[k,i])
        prob_r[k,i, np.nonzero(prob_r[k,i])] /= np.sum(prob_r[k,i])
        prob_r_projected[k,i, np.nonzero(prob_r_projected[k,i])] /= np.sum(prob_r_projected[k,i])
        if prob_a_cutoff[k,i] > 0.0:
            prob_a_cutoff[k,i] /= N_input_initial[k,i]
    print('Scattering matrix generated with', N_binaries[k], 'binaries for pertuber mass', M_ps[k], 'M_sol')
    print('Number of binaries removed due to cut-off:', N_broken_cutoff[k])

def eccentricAnomaly(es, Ms):
    Es = (Ms + np.copysign(1.0, np.sin(Ms))*0.85*es) % (2.0*np.pi)
    for i in range(np.size(es)):
        count = 0
        while (abs(Es[i] - es[i]*np.sin(Es[i]) - Ms[i]) > 10.0**(-8.0)):
            f = Es[i] - es[i]*np.sin(Es[i]) - Ms[i]
            f_p = 1.0 - es[i]*np.cos(Es[i])
            f_pp = es[i]*np.sin(Es[i])
            f_ppp = es[i]*np.cos(Es[i])

            d_1 = -f/f_p
            d_2 = -f/(f_p + 0.5*d_1*f_pp)
            d_3 = -f/(f_p + 0.5*d_2*f_pp + d_2**2.0*f_ppp/6.0)

            Es[i] += d_3

            count += 1
            if (count > 100):
                print("eccentricAnomaly did not converge, E =", Es[i], ', e =', es[i], ', M =', Ms[i])
                input()
                break
    return Es % (2.0*np.pi)

def generate_evolved_binaries_with_a_ini_i(mass_index, i, N_bins=1, r_proj_not_a=False):
    random_a_fins = np.zeros(N_bins)
    for j in range(N_bins):
        #Check if it's cutoff
        random_number1 = random.random()
        if (random_number1 <= prob_a_cutoff[mass_index, i]):
            random_a_fins[j] = -1.0
            continue
        #If it's not cut-off find the final semi-major axis
        random_number2 = random.random()
        for k in range(N_a_bins):
            if (r_proj_not_a):
                comparison_number = np.sum(prob_r_projected[mass_index, i,0:k+1])
            else:
                comparison_number = np.sum(prob_a[mass_index, i,0:k+1])
            if (random_number2 <= comparison_number):
                random_a_fins[j] = a_bins[k]
                break
        if (random_a_fins[j] == 0.0):
            print("No data for binaries with a_ini/au =", a_bins[i]/au)
            print("Continue?")
            input()
    return random_a_fins

def generate_evolved_binaries_from_initial_dist(mass_index, a_0, a_1, alpha, N_bins, replace_cutoffs=False, r_proj_not_a=False):
    a_inis = np.zeros(N_bins)
    if (alpha == 1.0):
        c = np.log(a_0)/np.log(a_1/a_0)
        for i in range(N_bins):
            a_inis[i] = (a_1/a_0)**(random.random() + c)
    else:
        for i in range(N_bins):
            a_inis[i] = (random.random()*(a_1**(1.0-alpha) - a_0**(1.0-alpha)) + a_0**(1.0-alpha))**(1.0/(1.0 - alpha))
    #Convert a_inis to r_proj_inis for finding L(0,0) 
    e_inis = np.random.random(size=N_bins)**(1.0/3.0)
    M = np.random.uniform(0.0, 2.0*np.pi, size=N_bins)
    E = eccentricAnomaly(e_inis, M)
    f = 2.0*np.arctan(np.sqrt((1.0+e_inis)/(1.0-e_inis))*np.tan(E/2.0))
    r_proj_inis = a_inis*(1.0 - e_inis**2.0)/(1.0 + e_inis*np.cos(f))
    #
    a_fins = np.zeros(N_bins)
    for j in range(N_bins):
        i = int(np.floor(np.log(a_inis[j]/a_min)/da))
        a_fins[j] = generate_evolved_binaries_with_a_ini_i(mass_index, i, 1, r_proj_not_a=r_proj_not_a)[0]
        if replace_cutoffs:
            while a_fins[j]<0.0:
                a_fins[j] = generate_evolved_binaries_with_a_ini_i(mass_index, i, 1, r_proj_not_a=r_proj_not_a)[0]
    return r_proj_inis

print('Importing real binary data')
#Binaries from AMR (1406.5164) catalog
s_MRA = np.zeros(0)
a_MRA = np.zeros(0)
t_dt_MRA = np.zeros(0)
d_MRA = np.zeros(0)
with open("AMR_data_s_a_t_dt_d.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        s_MRA = np.append(s_MRA, float(row[0]))
        a_MRA = np.append(a_MRA, float(row[1])*au)
        d_MRA = np.append(d_MRA, float(row[3])*parsec)
        if row[2]=='':
            t_dt_MRA = np.append(t_dt_MRA, -1.0)
        else:
            t_dt_MRA = np.append(t_dt_MRA, float(row[2]))
r_projected_MRA = s_MRA*d_MRA/parsec*au

plt.scatter(s_MRA, r_projected_MRA/au)
plt.xlabel('Angular separation, arcsec')
plt.ylabel('Projected separation, au')
plt.show()

plt.scatter(d_MRA/au, r_projected_MRA/au)
plt.xlabel('Distance, au')
plt.ylabel('Projected separation, au')
plt.show()
#Maximum time in disk
max_td = np.sort(t_dt_MRA[np.nonzero(t_dt_MRA+1.0)])[149]
#Remove real binaries with angular separation less than real_binary_cutoff
real_binary_cutoff = 1000.0*au
for i in range(N_a_bins):
    if real_binary_cutoff < a_bins[i]:
        real_binary_cutoff = a_bins[i]
        break
print('Projected separation lower cut-off, au =', real_binary_cutoff/au)
#Only keep those binaries with -1 < t_dt_MRA <= max_td and s > real_binary_cutoff
indices_to_keep = np.zeros(0, dtype=int)
for i in range(np.size(a_MRA)):
    #if ((t_dt_MRA[i]> -1.0) and (t_dt_MRA[i]<max_td) and (r_projected_MRA[i] > real_binary_cutoff)):
    if (r_projected_MRA[i] > real_binary_cutoff):
        indices_to_keep = np.append(indices_to_keep, i)
a_MRA = a_MRA[indices_to_keep]
s_MRA = s_MRA[indices_to_keep]
d_MRA = d_MRA[indices_to_keep]
t_dt_MRA = t_dt_MRA[indices_to_keep]
r_projected_MRA = r_projected_MRA[indices_to_keep]
#print(s_MRA)
#print(d_MRA/parsec)
#print(r_projected_MRA/au)
print('Keeping the', np.size(a_MRA), 'most halo-like binaries')

#Find i_real_min and i_real_max
i_real_min = int(np.floor(np.log(np.min(r_projected_MRA)/a_min)/da))
i_real_max = int(np.floor(np.log(np.max(r_projected_MRA)/a_min)/da))
print(i_real_min, i_real_max)

print('Binning real binaries')
#Bin real binaries
N_real_r_projected = np.zeros(N_a_bins)
for i in range(np.size(a_MRA)):
    try:
        j = int(np.floor(np.log(r_projected_MRA[i]/a_min)/da))
        N_real_r_projected[j] += 1
    except:
        print('Real binary semi-major axis/au:', a_MRA[i]/au)
        input()
'''
#Bin and plot in terms of angular separation
s_min = np.min(s_MRA)
print('s_min =', s_min)
s_max = np.max(s_MRA)
print('s_max =', s_max)
N_s_bins = 10
ds = (np.log(s_max) - np.log(s_min))/(N_s_bins-1.0)
s_bins = np.array([s_min*np.exp(i*ds) for i in range(N_s_bins)])
N_real_s = np.zeros(N_s_bins)
for i in range(np.size(s_MRA)):
    try:
        j = int(np.floor(np.log(s_MRA[i]/s_min)/ds))
        N_real_s[j] += 1
    except:
        print('Real binary angular separation/arcsec:', s_MRA[i])
        input()
plt.loglog(s_bins, N_real_s)
#plt.loglog(s_bins*np.exp(0.5*ds), N_real_s)
plt.xlabel('Angular separation, arsec')
plt.ylabel('Number of binaries')
plt.show()
'''
'''
#Plot distance distribution
d_min = np.min(d_MRA)
d_max = np.max(d_MRA)
N_d_bins = 20
dd = (np.log(d_max) - np.log(d_min))/(N_d_bins-1.0)
d_bins = np.array([d_min*np.exp(i*dd) for i in range(N_d_bins)])
N_real_d = np.zeros(N_d_bins)
for i in range(np.size(d_MRA)):
    try:
        j = int(np.floor(np.log(d_MRA[i]/d_min)/dd))
        N_real_d[j] += 1
    except:
        print('Real binary distance/pc:', d_MRA[i]/pc)
        input()
plt.loglog(d_bins/parsec, N_real_d)
plt.loglog(d_bins/parsec, 3.0*np.size(d_MRA)*(d_bins/parsec)**2.0/((d_max/parsec)**3.0-(d_min/parsec)**3.0)*(d_bins/parsec*(np.exp(dd)-1.0)))
#plt.loglog(s_bins*np.exp(0.5*ds), N_real_s)
plt.xlabel('Distance, pc')
plt.ylabel('Number of binaries')
plt.show()
'''
'''
#Just gonna make a cumulative distribution plot here
N_a_bins_cumulative = 500
N_cumulative = np.zeros(N_a_bins_cumulative)
da_cumulative = (np.log(s_max) - np.log(s_min))/(N_a_bins_cumulative)
a_bins_cumulative = np.array([s_min*np.exp(i*da_cumulative) for i in range(N_a_bins_cumulative)])
for i in range(np.size(a_MRA)):
    j = int(np.floor(np.log(s_MRA[i]/s_min)/da_cumulative))
    N_cumulative[j:] += 1
plt.semilogx(a_bins_cumulative, N_cumulative)
plt.show()
'''


print('Generating population of binaries from scattering matrix')
#Generate population of evolved binaries
a_virtual_evolved = np.zeros((N_Mps, np.size(a_MRA)))
r_virtual_projected = np.zeros((N_Mps, np.size(a_MRA)))
for k in range(N_Mps):
    #a_virtual_evolved[k] = generate_evolved_binaries_from_initial_dist(k, a_bins[i_real_min], a_bins[i_real_max], alpha, np.size(a_MRA), replace_cutoffs=False)
    r_virtual_projected[k] = generate_evolved_binaries_from_initial_dist(k, 10.0*au, 3.0*10.0**5.0*au, alpha, np.size(a_MRA), replace_cutoffs=True, r_proj_not_a=True)


'''
print('Calculating goodness of fit for this model: M_p=10Msol')
#Goodness of fit
#Generate a binary sample with the same number of binaries as the data
a_virtual_evolved_averaged = np.zeros(np.size(a_MRA))
for i in range(np.size(a_MRA)):
    a_virtual_evolved_averaged[i] = np.median(a_virtual_evolved[(np.size(a_virtual_evolved)//np.size(a_MRA))*i:(np.size(a_virtual_evolved)//np.size(a_MRA))*(i+1)])
#Sort the semi-major axes
a_virtual_evolved_averaged = np.sort(a_virtual_evolved_averaged)
a_MRA = np.sort(a_MRA)
print(a_virtual_evolved_averaged/au)
print(a_MRA/au)
#Compare the semi-major axes
sigma = 0.0
for i in range(np.size(a_MRA)):
    sigma += (a_MRA[i]/au - a_virtual_evolved_averaged[i]/au)**2.0
sigma = np.sqrt(sigma)/np.size(a_MRA)
sigma_0 = 1561.724988768935
print('sigma/sigma_0 =', sigma/sigma_0)
'''

print('Binning virtual binaries')
#Bin virtual binaries
N_virtual_evolved = np.zeros((N_Mps, N_a_bins))
for k in range(N_Mps):
    for i in range(np.size(r_virtual_projected[k])):
        if r_virtual_projected[k,i] >= 0.0:
            try:
                j = int(np.floor(np.log(r_virtual_projected[k,i]/a_min)/da))
                N_virtual_evolved[k,j] += 1
            except:
                print('Virtual evolved binary projected separation/au:', r_virtual_projected[k,i]/au)
                input()


print('Attempt at Yoo et al.s analysis')
P = np.zeros((N_Mps, N_a_bins))
for k in range(N_Mps):
    print('M_p =', M_ps[k], 'M_sol')
    #Find P (equation 15 Yoo et al.)
    #Take into account cutoff binaries in prob_a:
    for i in range(N_a_bins):
        prob_r_projected[k,i] *= 1.0 - prob_a_cutoff[k,i]
    #Set-up a_min and a_max limits
    #prob_r_projected[k,:,0:i_real_min] *= 0.0
    #prob_r_projected[k,:,i_real_max+1:] *= 0.0
    P[k] = np.matmul(np.transpose(prob_r_projected[k]), a_bins**(1.0-alpha))
    #Normalise
    #P /= np.sum(a_bins**(1.0-alpha))
    P[k] /= np.sum(P[k, i_real_min:i_real_max+1])
    P[k] *= np.size(a_MRA)
    #Let's plot this to see how it looks
    plt.plot(a_bins*np.exp(0.5*da)/au, P[k], label='Distribution from P_N, M_p = {}M_sol'.format(M_ps[k]), linestyle=linestyles[k])
    #Calculate likelihood function
    L = 1.0
    for i in range(N_a_bins):
        L *= P[k,i]**(N_real_r_projected[i])
    print('L =', L)
    print('lnL =', np.log(L))
    lnL_0 = 456.3159808642587
    print('sigma =', np.sqrt(2.0*(lnL_0-np.log(L)))) 


print('Plotting binned distributions')
#Edit real distribution to have 10 times fewer bins
reduction_factor = 1
N_real_condensed = np.zeros(np.size(N_real_r_projected)//reduction_factor)
for i in range(np.size(N_real_condensed)):
    N_real_condensed[i] = np.sum(N_real_r_projected[(i*reduction_factor):((i+1)*reduction_factor)])
a_bins_condensed = a_bins[0::reduction_factor]
if np.size(a_bins_condensed) > np.size(N_real_condensed):
    a_bins_condensed = a_bins_condensed[0:np.size(a_bins_condensed)-1]
#plt.loglog(a_bins/au, prob_a[40]*1000000, linestyle=':')
plt.loglog(a_bins/au, N_input_dist_r_proj[0]/np.sum(N_input_dist_r_proj[0,i_real_min:i_real_max+1])*np.size(a_MRA), label='Distribution from simulation')
#for k in range(N_Mps):
    #plt.plot(a_bins/au, N_virtual_evolved[k], label='Distribution from scattering matrix, M_p = {}M_sol'.format(M_ps[k]))
#plt.xlim([10.0**1.0, 3.0*10.0**5.5])
plt.scatter(a_bins_condensed*np.exp(0.5*(np.log(a_max/a_min)/np.size(a_bins_condensed)))/au, N_real_condensed/reduction_factor, label='Real binary distribution')
#plt.scatter(a_bins*np.exp(0.5*da)/au, N_real, label='Real binary distribution')
#plt.plot(a_bins*np.exp(0.5*da)/au, a_bins**(1.0-alpha)/(np.sum(a_bins[i_real_min:i_real_max+1]**(1.0-alpha)))*np.size(a_MRA), label=r'$\alpha =${}'.format(alpha))
#plt.xlim([a_min/au, np.min([a_bins[i_real_max]/au, a_max/au, r_cutoff/au])])
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.legend()
plt.xlabel('Projected separation, au')
plt.ylabel('Number of binaries')
plt.show()

print('Finished')
        