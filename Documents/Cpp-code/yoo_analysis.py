import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import au, parsec
import random
from scipy import stats
from scipy.special import factorial
plt.rc('font', family='serif')
random.seed(10)

a_min = 30.0*au
a_max = 20000.0*parsec
N_a_bins = 6
da = (np.log(a_max) - np.log(a_min))/(N_a_bins)
a_bins = np.array([a_min*np.exp(i*da) for i in range(N_a_bins)])

s_min = 3.5
s_max = 10000.0
N_s_bins = N_a_bins
ds = (np.log(s_max) - np.log(s_min))/(N_s_bins)
s_bins = np.array([s_min*np.exp(i*ds) for i in range(N_s_bins)])
print('s_bins =', s_bins)

r_cutoff = 10000.*parsec

M_ps = np.array([10, 100, 1000])
N_Mps = np.size(M_ps)
linestyles = np.array(['dotted', 'dashed', 'dashdot'])

rhos = np.array([0.007])
N_rhos = np.size(rhos)


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
#Maximum time in disk
max_td = np.sort(t_dt_MRA[np.nonzero(t_dt_MRA+1.0)])[105]
#Remove real binaries with projected separation less than real_binary_cutoff (to the nearest bin boundary)
real_binary_cutoff = 21.5 
for i in range(N_a_bins):
    if real_binary_cutoff < s_bins[i]:
        real_binary_cutoff = s_bins[i]
        break
print('Projected separation lower cut-off, au =', real_binary_cutoff/au)
#Only keep those binaries with -1 < t_dt_MRA <= max_td and r_proj > real_binary_cutoff
indices_to_keep = np.zeros(0, dtype=int)
for i in range(np.size(a_MRA)):
    #if ((t_dt_MRA[i]> -1.0) and (t_dt_MRA[i]<max_td)):
    #if (r_projected_MRA[i] > real_binary_cutoff):
    if(s_MRA[i] > real_binary_cutoff):
        indices_to_keep = np.append(indices_to_keep, i)
#indices_to_keep = np.random.randint(low=0, high=np.size(a_MRA), size=79)
a_MRA = a_MRA[indices_to_keep]
s_MRA = s_MRA[indices_to_keep]
d_MRA = d_MRA[indices_to_keep]
t_dt_MRA = t_dt_MRA[indices_to_keep]
r_projected_MRA = r_projected_MRA[indices_to_keep]
print('Keeping', np.size(s_MRA), 'real binaries')

#Find i_real_min and i_real_max
i_real_min = int(np.floor(np.log(np.min(r_projected_MRA)/a_min)/da))
i_real_max = int(np.floor(np.log(np.max(r_projected_MRA)/a_min)/da))
i_real_min = int(np.floor(np.log(np.min(s_MRA)/s_min)/ds))
i_real_max = int(np.floor(np.log(np.max(s_MRA)/s_min)/ds))
print(i_real_min, i_real_max)

print('Binning real binaries')
#Bin real binaries
N_real_r_projected = np.zeros(N_a_bins)
for i in range(np.size(a_MRA)):
    try:
        #j = int(np.floor(np.log(r_projected_MRA[i]/a_min)/da))
        j = int(np.floor(np.log(s_MRA[i]/s_min)/ds))
        N_real_r_projected[j] += 1
    except:
        print('Real binary semi-major axis/au:', a_MRA[i]/au)
        print('s =', s_MRA[i])
        input()

print('a_bins/au =', a_bins/au)
print('Real binaries =', N_real_r_projected)
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

#Set up scattering matrices
print('Importing simulation data and generating scattering matrices')
#Probability that a binary with initial semi-major axis a_bins[i] will have final projected separation a_bins[j] is prob[?,?,i,j] given that it never exceeds the cut-off semi-major axis
prob = np.zeros((N_Mps, N_rhos, N_a_bins, N_a_bins))
prob_cutoff = np.zeros((N_Mps, N_rhos, N_a_bins))
prob_cutoff_ini = np.zeros((N_Mps, N_rhos, N_a_bins))
prob_Mp0 = np.zeros((N_Mps, N_rhos, N_a_bins, N_a_bins))
N_input_dist_a = np.zeros((N_Mps, N_rhos, N_a_bins))
N_input_dist_r_proj = np.zeros((N_Mps, N_rhos, N_a_bins))
N_input_initial_a = np.zeros((N_Mps, N_rhos, N_a_bins))
N_input_initial_r_proj = np.zeros((N_Mps, N_rhos, N_a_bins))
N_broken_cutoff = np.zeros((N_Mps, N_rhos), dtype=int)
N_binaries = np.zeros((N_Mps, N_rhos), dtype=int)
for k in range(N_Mps):
    for l in range(N_rhos):
        with open("final_r_and_a_distributions_MRAparams_Mp{}_Nbin10e6_format_ai_ri_ei_af_rf_ef.csv".format(M_ps[k])) as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            for row in reader:
                a_ini = float(row[0])
                r_ini = float(row[1])
                a_fin = float(row[3])
                r_fin = float(row[4])
                if (a_ini==a_fin):
                    r_fin = r_ini
                i = int(np.floor(np.log(a_ini/a_min)/da))
                N_input_initial_a[k,l,i] += 1
                #Find initial r projected
                u = random.uniform(-1.0, 1.0)
                theta = random.uniform(0.0, 2.0*np.pi)
                rand_vec = np.array([np.sqrt(1.0-u**2.0)*np.cos(theta), np.sqrt(1.0-u**2.0)*np.sin(theta), u])
                rand_vec /= np.linalg.norm(rand_vec)
                r_ini_projected = r_ini*np.sqrt(1.0-rand_vec[0]**2.0)
                j_ri_proj = int(np.floor(np.log(r_ini_projected/a_min)/da))
                distance = np.random.choice(d_MRA)
                s_ini = r_ini_projected/au/(distance/parsec)
                s_ini = a_ini/au/1.4
                j_si = int(np.floor(np.log(s_ini/s_min)/ds))
                try:
                    prob_Mp0[k,l,i,j_si] += 1
                    N_input_initial_r_proj[k,l,j_si] += 1
                except:
                    N_broken_cutoff[k,l] += 1
                    prob_cutoff[k,l,i] += 1
                    prob_cutoff_ini[k,l,i] += 1
                    continue
                if (a_fin > 0.0):
                    try:
                        j = int(np.floor(np.log(a_fin/a_min)/da))
                        N_input_dist_a[k,l,j] += 1
                        #Find final r projected
                        u = random.uniform(-1.0, 1.0)
                        theta = random.uniform(0.0, 2.0*np.pi)
                        rand_vec = np.array([np.sqrt(1.0-u**2.0)*np.cos(theta), np.sqrt(1.0-u**2.0)*np.sin(theta), u])
                        rand_vec /= np.linalg.norm(rand_vec)
                        r_fin_projected = r_fin*np.sqrt(1.0-rand_vec[0]**2.0)
                        #
                        j_r_proj = int(np.floor(np.log(r_fin_projected/a_min)/da))
                        s_fin = r_fin_projected/au/(distance/parsec)
                        s_fin = a_fin/au/1.4
                        j_s = int(np.floor(np.log(s_fin/s_min)/ds))
                        prob[k,l,i,j_s] += 1
                        N_input_dist_r_proj[k,l,j_r_proj] += 1
                    except:
                        if (s_fin > s_max) or (s_fin < s_min):
                            N_broken_cutoff[k,l] += 1
                            prob_cutoff[k,l,i] += 1
                            continue
                        if ((r_fin < r_cutoff) and (a_fin < r_cutoff)):
                            print(j_s)
                            print(r_fin_projected/au)
                            print(r_ini/au, r_fin/au)
                            print(a_ini/au, a_fin/au)
                            input()
                        else:
                            N_broken_cutoff[k,l] += 1
                            prob_cutoff[k,l,i] += 1
                            continue
                    N_binaries[k,l] += 1
                else:
                    N_broken_cutoff[k,l] += 1
                    prob_cutoff[k,l,i] += 1
        #Normalise to make it a probability
        for i in range(N_a_bins):
            prob[k,l,i, np.nonzero(prob[k,l,i])] /= np.sum(prob[k,l,i])
            prob_Mp0[k,l,i, np.nonzero(prob_Mp0[k,l,i])] /= np.sum(prob_Mp0[k,l,i])
            if prob_cutoff[k,l,i] > 0.0:
                prob_cutoff[k,l,i] /= N_input_initial_a[k,l,i]
            if prob_cutoff_ini[k,l,i] > 0.0:
                prob_cutoff_ini[k,l,i] /= N_input_initial_a[k,l,i]
        print('Scattering matrix generated with', N_binaries[k,l], 'binaries for pertuber mass', M_ps[k], 'M_sol and dark matter density', rhos[l], 'M_sol/pc^3')
        print('Number of binaries removed due to cut-off:', N_broken_cutoff[k,l])

print('Initial simulation distribution =', N_input_initial_r_proj[0,0]/np.sum(N_input_initial_r_proj[0,0,i_real_min:i_real_max+1])*np.size(s_MRA))

print('Starting Yoo et al.s analysis')
#Alphas to test
alphs = np.linspace(1.0, 3.0, num=2000)
lnLs = np.zeros(np.size(alphs))
#Results
sigmas = np.zeros((N_Mps, N_rhos))
alphas = np.zeros((N_Mps, N_rhos))
P = np.zeros((N_Mps, N_rhos, N_a_bins))
for k in range(N_Mps):
    for l in range(N_rhos):
        #Find P (equation 15 Yoo et al.)
        #Take into account cutoff binaries in prob:
        for i in range(N_a_bins):
            prob[k,l,i] *= 1.0 - prob_cutoff[k,l,i]
            prob_Mp0[k,l,i] *= 1.0 - prob_cutoff_ini[k,l,i]
        #Set-up a_min and a_max limits
        #prob[k,l,0:i_real_min] *= 0.0
        #prob[k,l,i_real_max+1:] *= 0.0
        #Find P_Mp0
        #Set-up a_min and a_max limits
        #prob_Mp0[k,l,0:i_real_min] *= 0.0
        #prob_Mp0[k,l,i_real_max+1:] *= 0.0
        #Find L(0,0) (L_0) and corresponding alpha (alpha_0)
        for j in range(np.size(alphs)):
            P[k,l] = np.matmul(np.transpose(prob_Mp0[k,l]), a_bins**(1.0-alphs[j]))
            P[k,l] /= np.sum(P[k,l,i_real_min:i_real_max+1])
            P[k,l] *= np.size(a_MRA)
            #likelihood function
            lnLs[j] = 0.0
            for i in range(i_real_min, i_real_max+1):
                if (P[k,l,i]>0.0):
                    #lnLs[j] += N_real_r_projected[i]*np.log(P[k,l,i])
                    lnLs[j] += N_real_r_projected[i]*np.log(P[k,l,i])-np.log(float(factorial(N_real_r_projected[i], exact=True)))-P[k,l,i]
        lnL_0 = np.max(lnLs)
        print('lnL_0 =', lnL_0)
        alpha_0 = alphs[np.argmax(lnLs)]
        alpha_0 = 1.0
        P[k,l] = np.matmul(np.transpose(prob_Mp0[k,l]), a_bins**(1.0-alpha_0))
        P[k,l] /= np.sum(P[k,l,i_real_min:i_real_max+1])
        P[k,l] *= np.size(a_MRA)
        plt.plot(s_bins*np.exp(0.5*ds), P[k,l], label='Distribution from P_N, M_p = 0M_sol')
        print('P =', P[k,l])
        #print('Contributions by bin (Mp=0) =', [N_real_r_projected[i]*np.log(P[k,l,i])-np.log(float(factorial(N_real_r_projected[i], exact=True)))-P[k,l,i] if ((P[k,l,i]>0.0) and (i_real_min<=i<i_real_max+1)) else 0.0 for i in range(N_a_bins)])
        #Find sigma and corresponding alpha
        for j in range(np.size(alphs)):
            P[k,l] = np.matmul(np.transpose(prob[k,l]), a_bins**(1.0-alphs[j]))
            P[k,l] /= np.sum(P[k,l,i_real_min:i_real_max+1])
            P[k,l] *= np.size(a_MRA)
            #likelihood function
            lnLs[j] = 0.0
            for i in range(i_real_min, i_real_max+1):
                if (P[k,l,i]>0.0):
                    #lnLs[j] += N_real_r_projected[i]*np.log(P[k,l,i])
                    lnLs[j] += N_real_r_projected[i]*np.log(P[k,l,i])-np.log(float(factorial(N_real_r_projected[i], exact=True)))-P[k,l,i]
        print('lnL =', np.max(lnLs))
        #input()
        sigmas[k,l] = np.sqrt(2.0*(lnL_0-np.max(lnLs)))
        alphas[k,l] = alphs[np.argmax(lnLs)]
print('Value of alpha that corresponds with L(0,0) =', alpha_0)
print('Results for sigma')
print('Perturber masses, M_sol =', M_ps)
print('Dark matter densities, M_sol/pc^3 =', rhos)
print('Sigma =', sigmas)
print('Corresponding alphas =', alphas)

print('Plotting Yoo et al. fig 5')
#Plot Yoo et al. fig 5
for k in np.array([0,1,2]):
    for l in np.array([0]):
        P[k,l] = np.matmul(np.transpose(prob[k,l]), a_bins**(1.0-alphas[k,l]))
        #Normalise (necessary for plotting, not analysis)
        #P /= np.sum(a_bins**(1.0-alpha))
        P[k,l] /= np.sum(P[k,l,i_real_min:i_real_max+1])
        P[k,l] *= np.size(a_MRA)
        print('P =', P[k,l])
        print('Contributions by bin =', [N_real_r_projected[i]*np.log(P[k,l,i])-np.log(float(factorial(N_real_r_projected[i], exact=True)))-P[k,l,i] if ((P[k,l,i]>0.0) and (i_real_min<=i<i_real_max+1)) else 0.0 for i in range(N_a_bins)])
        #plt.plot(a_bins*np.exp(0.5*da)/au, P[k,l], label='Distribution from P_N, M_p = {}M_sol'.format(M_ps[k]), linestyle=linestyles[k])
        plt.plot(s_bins*np.exp(0.5*ds), P[k,l], label='Distribution from P_N, M_p = {}M_sol'.format(M_ps[k]), linestyle=linestyles[k])
#plt.plot(a_bins*np.exp(0.5*da)/au, np.size(a_MRA)*np.matmul(np.transpose(prob_Mp0[0,0]), a_bins**(1.0-alpha_0))/np.sum(np.matmul(np.transpose(prob_Mp0[0,0]), a_bins**(1.0-alpha_0))[i_real_min:i_real_max+1]), label='Distribution from P_N, M_p = 0M_sol'.format(M_ps[k]), linestyle=linestyles[k])
#plt.plot(s_bins*np.exp(0.5*ds), np.size(s_MRA)*np.matmul(np.transpose(prob_Mp0[0,0]), a_bins**(1.0-alpha_0))/np.sum(np.matmul(np.transpose(prob_Mp0[0,0]), a_bins**(1.0-alpha_0))[i_real_min:i_real_max+1]), label='Distribution from P_N, M_p = 0M_sol'.format(M_ps[k]), linestyle=linestyles[k])
#plt.plot(a_bins*np.exp(0.5*da)/au, N_input_initial_r_proj[0,0], label='Initial distribution from simulation')
plt.plot(s_bins*np.exp(0.5*ds), N_input_initial_r_proj[0,0]/np.sum(N_input_initial_r_proj[0,0,i_real_min:i_real_max+1])*np.size(a_MRA), label='Input s dist')
#Edit real distribution to have reduction_factor times fewer bins
reduction_factor = 1
N_real_condensed = np.zeros(np.size(N_real_r_projected)//reduction_factor)
for i in range(np.size(N_real_condensed)):
    N_real_condensed[i] = np.sum(N_real_r_projected[(i*reduction_factor):((i+1)*reduction_factor)])
a_bins_condensed = a_bins[0::reduction_factor]
s_bins_condensed = s_bins[0::reduction_factor]
if np.size(a_bins_condensed) > np.size(N_real_condensed):
    a_bins_condensed = a_bins_condensed[0:np.size(a_bins_condensed)-1]
if np.size(s_bins_condensed) > np.size(N_real_condensed):
    s_bins_condensed = s_bins_condensed[0:np.size(s_bins_condensed)-1]
#plt.scatter(a_bins_condensed*np.exp(0.5*(np.log(a_max/a_min)/np.size(a_bins_condensed)))/au, N_real_condensed/reduction_factor, label='Real binary distribution')
plt.scatter(s_bins_condensed*np.exp(0.5*(np.log(s_max/s_min)/np.size(s_bins_condensed))), N_real_condensed/reduction_factor, label='Real binary distribution')
#plt.xlim([a_bins[i_real_min]/au, a_bins[i_real_max+1]/au])
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.legend()
plt.xlabel('Projected separation, au')
plt.ylabel('Number of binaries')
plt.show()

print('Finished')