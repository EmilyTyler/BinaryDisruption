import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import au, parsec
from scipy.stats import chisquare
from scipy import stats
from scipy.optimize import curve_fit
from mc_draw_from_dist import draw_from_poly, seed
from numpy.polynomial import Polynomial
plt.rc('font', family='serif', size=14)
np.random.seed(seed)

#Settings
s_analysis_not_r_proj = False
#Parameters to test
alphs = np.linspace(0.7, 2.0, num=200)
As = np.linspace(0.0, 1.0, num=100)
#Corresponding number of degrees of freedom
ddof = 2

#Semi-major axis/projected separation bins
a_min = 30.0*au
a_max = 20000.0*parsec
N_a_bins = 25
da = (np.log(a_max) - np.log(a_min))/(N_a_bins)
a_bins = np.array([a_min*np.exp(i*da) for i in range(N_a_bins)])
das = a_bins*(np.exp(da)-1.0)
print('a_bins/au =', a_bins/au)

#Angular separation bins
s_min = 0.1
s_max = a_max/au/1.4
N_s_bins = N_a_bins
ds = (np.log(s_max) - np.log(s_min))/(N_s_bins)
s_bins = np.array([s_min*np.exp(i*ds) for i in range(N_s_bins)])
#print('s_bins =', s_bins)

r_cutoff = 10000.*parsec

#Perturber masses
M_ps = np.array([10, 30, 50, 70, 79, 89, 100, 173, 300, 500, 595, 707, 840, 1000, 2000, 3000])
N_Mps = np.size(M_ps)
linestyles = np.array(['dotted', 'dashed', 'dashdot', 'dotted', 'dashed', 'dashdot', 'dotted', 'dashed', 'dashdot', 'dotted', 'dashed', 'dashdot', 'dotted', 'dashed', 'dashdot', 'dotted'])

#Perturber densities
rhos = np.array([0.0005, 0.0009, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.009, 0.012])
N_rhos = np.size(rhos)

#Import projected separations of observed binaries
print('Importing real binary data')
#Binaries from AMR (1406.5164) catalog
s_MRA = np.zeros(0)
a_MRA = np.zeros(0)
t_dt_MRA = np.zeros(0)
d_MRA = np.zeros(0)
with open("data/AMR_data_s_a_t_dt_d.csv") as csvfile:
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

AMR_data_names = np.zeros(0, dtype=str)
with open("data/AMR_data_names.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        AMR_data_names = np.append(AMR_data_names, row[0])

#Import galpy time in disk
t_dt_galpy = np.full(np.size(t_dt_MRA), -1.0)
with open("data/AMR_data_names_galpyt_d.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        label = row[0]
        i = np.argwhere(AMR_data_names == label)
        t_dt_galpy[i] = float(row[1])

#Maximum time in disk
#max_td = np.sort(t_dt_MRA[np.nonzero(t_dt_MRA+1.0)])[25]
max_td = np.sort(t_dt_galpy[np.nonzero(t_dt_galpy+1.0)])[38]
#Remove real binaries with projected separation less than real_binary_cutoff (to the nearest bin boundary)
indices_to_keep = np.zeros(0, dtype=int)
if (s_analysis_not_r_proj):
    real_binary_cutoff = 3.5
    for i in range(N_a_bins):
        if real_binary_cutoff <= s_bins[i]:
            real_binary_cutoff = s_bins[i]
            break
    print('Angular separation lower cut-off, arcsec =', real_binary_cutoff)
    for i in range(np.size(a_MRA)):
        #if ((t_dt_MRA[i]> -1.0) and (t_dt_MRA[i]<max_td)):
        if(s_MRA[i] > real_binary_cutoff):
            #if ((t_dt_galpy[i]> -1.0) and (t_dt_galpy[i]<max_td)):
                indices_to_keep = np.append(indices_to_keep, i)
else:
    real_binary_cutoff = 1000.0*au
    for i in range(N_a_bins):
        if real_binary_cutoff <= a_bins[i]:
            real_binary_cutoff = a_bins[i]
            break
    print('Projected separation lower cut-off, au =', real_binary_cutoff/au)
    for i in range(np.size(a_MRA)):
    #if ((t_dt_MRA[i]> -1.0) and (t_dt_MRA[i]<max_td)):
        #if ((t_dt_galpy[i]> -1.0) and (t_dt_galpy[i]<max_td)):
        if (r_projected_MRA[i] > real_binary_cutoff):
            indices_to_keep = np.append(indices_to_keep, i)


#Plot like SUPERWIDE fig 28
plt.figure(figsize=(7.5,5))
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.scatter(r_projected_MRA/au, d_MRA/parsec, color='tab:blue', marker='x')
plt.plot([1273]*np.size(d_MRA), d_MRA/parsec, label='1273au', color='tab:orange')
plt.xlabel('Projected separation / au')
plt.ylim([np.min(d_MRA/parsec), np.max(d_MRA/parsec)])
plt.ylabel('Distance / pc')
plt.legend()
plt.tight_layout()
plt.show()


a_MRA = a_MRA[indices_to_keep]
s_MRA = s_MRA[indices_to_keep]
d_MRA = d_MRA[indices_to_keep]
t_dt_MRA = t_dt_MRA[indices_to_keep]
r_projected_MRA = r_projected_MRA[indices_to_keep]

print('Keeping', np.size(s_MRA), 'real binaries')

#Find i_real_min and i_real_max
if (s_analysis_not_r_proj):
    i_real_min = int(np.floor(np.log(np.min(s_MRA)/s_min)/ds))
    i_real_max = int(np.floor(np.log(np.max(s_MRA)/s_min)/ds))
else:
    i_real_min = int(np.floor(np.log(np.min(r_projected_MRA)/a_min)/da))
    i_real_max = int(np.floor(np.log(np.max(r_projected_MRA)/a_min)/da))

print('Binning real binaries')
#Bin real binaries
N_real_r_projected = np.zeros(N_a_bins)
for i in range(np.size(a_MRA)):
    try:
        if (s_analysis_not_r_proj):
            j = int(np.floor(np.log(s_MRA[i]/s_min)/ds))
        else:
            j = int(np.floor(np.log(r_projected_MRA[i]/a_min)/da))
        N_real_r_projected[j] += 1
    except:
        print('Real binary semi-major axis/au:', a_MRA[i]/au)
        print('s =', s_MRA[i])
        input()

print('Real binaries =', N_real_r_projected)
'''
#Plot real distribution
plt.figure(figsize=(7.5,5))
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
if (s_analysis_not_r_proj):
    plt.scatter(s_bins*np.exp(0.5*ds), N_real_r_projected, marker='x')
    plt.xlabel('Angular separation / arcsec')
else:
    plt.scatter(a_bins*np.exp(0.5*da)/au, N_real_r_projected, marker='x')
    plt.xlabel('Projected separation / au')
plt.ylabel('Number of binaries')
plt.tight_layout()
plt.show()
'''
'''
#Bin and plot in terms of angular separation
plt.figure(figsize=(7.5,5))
s_min = np.min(s_MRA)
#print('s_min =', s_min)
s_max = np.max(s_MRA)
#print('s_max =', s_max)
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
#plt.loglog(s_bins, N_real_s)
plt.scatter(s_bins*np.exp(0.5*ds), N_real_s, marker='x')
plt.xlabel('Angular separation / arcsec')
plt.ylabel('Number of binaries')
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.tight_layout()
plt.show()
'''

def power_law_a(x, a):
    return a*x**(a-1.0)

#Plot distance distribution
#This is all only necessary for angular separation analysis where we have to generate a distance distribtion to draw from
d_min = np.min(d_MRA)
d_max = np.max(d_MRA)+1.0*parsec
N_d_bins = 1000
dd = (np.log(d_max) - np.log(d_min))/(N_d_bins)
d_bins = np.array([d_min*np.exp(i*dd) for i in range(N_d_bins)])
dds = d_bins*(np.exp(dd)-1.0)
N_real_d = np.zeros(N_d_bins)
for i in range(np.size(d_MRA)):
    try:
        j = int(np.floor(np.log(d_MRA[i]/d_min)/dd))
        N_real_d[j] += 1
    except:
        print('Real binary distance/pc:', d_MRA[i]/parsec)
        input()
#Find best fit
coeffs = np.polyfit(d_bins/parsec, N_real_d/np.sum(N_real_d)/(dds/parsec), deg=11)
terms = np.array([coeffs[i]*(d_bins/parsec)**(np.size(coeffs)-1-i) for i in range(np.size(coeffs))])
#Plot this cumulatively
def make_cumulative(arr):
    cumulative_arr = np.zeros(np.size(arr))
    for i in range(np.size(arr)):
        cumulative_arr[i] = np.sum(arr[:(i+1)])
    return cumulative_arr
plt.plot(d_bins/parsec*np.exp(0.5*dd), make_cumulative(N_real_d)/np.sum(N_real_d), label='Real binaries')
plt.plot(d_bins/parsec*np.exp(0.5*dd), make_cumulative(np.sum(terms, axis=0)*dds)/np.sum(np.sum(terms, axis=0)*dds), label='Polynomial fit')
ax = plt.gca()
ax.set_xscale('log')
plt.xlabel('Distance, pc')
plt.ylabel('Cumulative distribution')
plt.legend()
plt.show()

#Find power law fit to initial conditions from Daniel's thesis
Fdims = [1.6,2.0,2.6,3.0]
plotcolors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
a2_min = 10.0*au
a2_max = 21000.0*au
N_a2_bins = 1000
da2 = (np.log(a2_max) - np.log(a2_min))/(N_a2_bins)
a2_bins = np.array([a2_min*np.exp(i*da2) for i in range(N_a2_bins)])
da2s = a2_bins*(np.exp(da2)-1.0)

def power_law(x, b, alpha, c):
    return b*x**(1.0-alpha) + c

def differentiated_power_law(x, b, alpha, c):
    return abs(b*(1.0-alpha)*x**(-alpha))

plt.figure(figsize=(7.5,5))
for k in range(np.size(Fdims)):  
    semimajor_axis = np.zeros(0)
    cumulative_dist = np.zeros(0)
    with open("data/DanielGriffithsThesisFig5_13Fdim{}_{}.csv".format(str(Fdims[k])[0], str(Fdims[k])[2])) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            semimajor_axis = np.append(semimajor_axis, float(row[0])*au)
            cumulative_dist = np.append(cumulative_dist, float(row[1]))
    #Fit a power law to this distribution
    pars, cov = curve_fit(f=power_law, xdata=semimajor_axis[np.nonzero(semimajor_axis)]/au, ydata=cumulative_dist[np.nonzero(semimajor_axis)])
    hybrid_alpha = pars[1]
    hybrid_c = pars[2]
    print('Fitted alpha =', hybrid_alpha)
    print('Fitted c =', hybrid_c)
    plt.plot(semimajor_axis/au, cumulative_dist, color=plotcolors[k], label=Fdims[k])
    plt.plot(semimajor_axis[np.nonzero(semimajor_axis)]/au, power_law(semimajor_axis[np.nonzero(semimajor_axis)]/au, *pars), color=plotcolors[k], linestyle='--')
    #Differentiate to get prob_density and convert to our bins
    lognormalmean = np.log(100.0*au/au)
    lognormalsigma = 1.5
    log_normal_addition = (np.exp(-(np.log(a2_bins/au) - lognormalmean)**2 / (2 * lognormalsigma**2))/(a2_bins/au * lognormalsigma * np.sqrt(2 * np.pi)))
    prob_density = differentiated_power_law(a2_bins/au, *pars)/(np.sum(differentiated_power_law(a2_bins/au, *pars)*(da2s/au)))
    #Add log normal addition
    prob_density += log_normal_addition
    #Renormalise
    prob_density /= np.sum(prob_density*(da2s/au))
plt.xlabel('Semi-major axis / au')
plt.ylabel('Cumulative distribution')
plt.legend()
plt.tight_layout()
plt.show()

def initial_semimajor_axis_dist(x, alpha, A):
    if (alpha == 1.0):
        power_law_part = A * (x/au)**(-alpha)/(np.log(a_max/a_min))
    else:
        power_law_part = A * (1.0-alpha)  * (x/au)**(-alpha) /((a_max/au)**(1.0-alpha) - (a_min/au)**(1.0-alpha))
    log_normal_part = (1.0 - A) / (np.sqrt(2.0*np.pi)*lognormalsigma*(x/au)) * np.exp(-(np.log(x/au) - lognormalmean)**2.0/(2.0*lognormalsigma**2.0))
    return power_law_part + log_normal_part

#Set up scattering matrices
print('Importing simulation data and generating scattering matrices')
#Probability that a binary with initial semi-major axis a_bins[i] will have final projected separation a_bins[j] is prob[?,?,i,j] given that it never exceeds the cut-off semi-major axis
prob = np.zeros((N_Mps, N_rhos, N_a_bins, N_a_bins))
prob_Mp0 = np.zeros((N_Mps, N_rhos, N_a_bins, N_a_bins))
N_input_dist_a = np.zeros((N_Mps, N_rhos, N_a_bins))
N_input_dist_r_proj = np.zeros((N_Mps, N_rhos, N_a_bins))
N_input_initial_a = np.zeros((N_Mps, N_rhos, N_a_bins))
N_input_initial_r_proj = np.zeros((N_Mps, N_rhos, N_a_bins))
N_broken_cutoff = np.zeros((N_Mps, N_rhos), dtype=int)
N_binaries = np.zeros((N_Mps, N_rhos), dtype=int)
distances = draw_from_poly(coeffs, d_min/parsec, d_max/parsec, size=100000)*parsec
#distances = np.array([100.0*parsec]*100000)
for k in range(N_Mps):
    for l in range(N_rhos):
        with open("output/final_r_and_a_distributions_rho0_{}_Mp{}_vrel_220_Nbin10e5_format_ai_ri_ei_af_rf_ef.csv".format(str(rhos[l])[2:], M_ps[k])) as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            row_number = 0
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
                u = np.random.uniform(-1.0, 1.0)
                theta = np.random.uniform(0.0, 2.0*np.pi)
                rand_vec = np.array([np.sqrt(1.0-u**2.0)*np.cos(theta), np.sqrt(1.0-u**2.0)*np.sin(theta), u])
                rand_vec /= np.linalg.norm(rand_vec)
                r_ini_projected = r_ini*np.sqrt(1.0-rand_vec[0]**2.0)
                #Index corresponding to this initial projected separation
                j_ri_proj = int(np.floor(np.log(r_ini_projected/a_min)/da))
                distance = distances[row_number]
                s_ini = r_ini_projected/au/(distance/parsec)
                #Index corresponding to this initial angular separation
                j_si = int(np.floor(np.log(s_ini/s_min)/ds))
                try:
                    if (s_analysis_not_r_proj):
                        prob_Mp0[k,l,i,j_si] += 1
                        N_input_initial_r_proj[k,l,j_si] += 1
                    else:
                        prob_Mp0[k,l,i,j_ri_proj] += 1
                        N_input_initial_r_proj[k,l,j_ri_proj] += 1
                except:
                    N_broken_cutoff[k,l] += 1
                    continue
                #If the binary didn't break
                if (a_fin > 0.0):
                    try:
                        j = int(np.floor(np.log(a_fin/a_min)/da))
                        N_input_dist_a[k,l,j] += 1
                        #Find final r projected
                        u = np.random.uniform(-1.0, 1.0)
                        theta = np.random.uniform(0.0, 2.0*np.pi)
                        rand_vec = np.array([np.sqrt(1.0-u**2.0)*np.cos(theta), np.sqrt(1.0-u**2.0)*np.sin(theta), u])
                        rand_vec /= np.linalg.norm(rand_vec)
                        r_fin_projected = r_fin*np.sqrt(1.0-rand_vec[0]**2.0)
                        #Index corresponding to final projected separation
                        j_r_proj = int(np.floor(np.log(r_fin_projected/a_min)/da))
                        s_fin = r_fin_projected/au/(distance/parsec)
                        #Index corresponding to final angular separation
                        j_s = int(np.floor(np.log(s_fin/s_min)/ds))
                        if (s_analysis_not_r_proj):
                            prob[k,l,i,j_s] += 1
                            N_input_dist_r_proj[k,l,j_s] += 1
                        else:
                            prob[k,l,i,j_r_proj] += 1
                            N_input_dist_r_proj[k,l,j_r_proj] += 1
                    except:
                        if ((r_fin < r_cutoff) and (a_fin < r_cutoff)):
                            print(j_s)
                            print(r_fin_projected/au)
                            print(r_ini/au, r_fin/au)
                            print(a_ini/au, a_fin/au)
                            input()
                        else:
                            N_broken_cutoff[k,l] += 1
                            continue
                    N_binaries[k,l] += 1
                else:
                    N_broken_cutoff[k,l] += 1
                row_number += 1

print('Starting analysis')
#Results
alphas = np.zeros((N_Mps, N_rhos))
alpha_0 = np.zeros((N_Mps, N_rhos))
As_Mp0 = np.zeros((N_Mps, N_rhos))
As_Mp = np.zeros((N_Mps, N_rhos))
P = np.zeros((N_Mps, N_rhos, N_a_bins))
P_minchi2 = np.zeros((N_Mps, N_rhos, N_a_bins))
chi2_Mp0 = np.zeros((N_Mps, N_rhos))
chi2_Mps = np.zeros((N_Mps, N_rhos))
chi2_Mp0 += np.inf
chi2_Mps += np.inf
for k in range(N_Mps):
    for l in range(N_rhos):
        #Find P (equation 15 Yoo et al.)
        #Start with Mp=0
        for j in range(np.size(alphs)):
            for m in range(np.size(As)):
                P[k,l] = np.matmul(np.transpose(prob_Mp0[k,l]), a_bins*initial_semimajor_axis_dist(a_bins, alphs[j], As[m]))
                if np.sum(P[k,l,i_real_min:i_real_max+1]) == 0.0:
                    break
                P[k,l] /= np.sum(P[k,l,i_real_min:i_real_max+1])
                P[k,l] *= np.size(a_MRA)
                #Calculate chi^2
                chi2, p_value = chisquare(f_obs=N_real_r_projected[np.nonzero(N_real_r_projected)], f_exp=P[k,l][np.nonzero(N_real_r_projected)], ddof=ddof)
                #Update chi^2 if it's smaller
                if (chi2_Mp0[k,l] > chi2):
                    chi2_Mp0[k,l] = chi2
                    alpha_0[k,l] = alphs[j]
                    As_Mp0[k,l] = As[m]
        #Same process but for non-zero pertuber masses
        for j in range(np.size(alphs)):
            for m in range(np.size(As)):
                P[k,l] = np.matmul(np.transpose(prob[k,l]), a_bins*initial_semimajor_axis_dist(a_bins, alphs[j], As[m]))
                if np.sum(P[k,l,i_real_min:i_real_max+1]) == 0.0:
                    break
                P[k,l] /= np.sum(P[k,l,i_real_min:i_real_max+1])
                P[k,l] *= np.size(a_MRA)
                #Calculate chi^2
                chi2, p_value = chisquare(f_obs=N_real_r_projected[np.nonzero(N_real_r_projected)], f_exp=P[k,l][np.nonzero(N_real_r_projected)], ddof=ddof)
                #Update chi^2 if it's smaller
                if (chi2_Mps[k,l] > chi2):
                    chi2_Mps[k,l] = chi2
                    alphas[k,l] = alphs[j]
                    As_Mp[k,l] = As[m]
                    P_minchi2[k,l] = P[k,l]  

#Calculate Y^2 and p-values
print('Plotting constraints from Lucy multinomial Y^2')
Y2 = np.zeros((N_Mps, N_rhos))
p_values_Y2 = np.zeros((N_Mps, N_rhos))
I = i_real_max - i_real_min + 1
nu = I - 1 - ddof
for k in range(N_Mps):
    V_X2 = 2.0*nu + 1.0/np.size(a_MRA)*(np.sum(1/(P_minchi2[k,:,i_real_min:i_real_max+1]/np.size(a_MRA)), axis=1) - I**2.0 - 2.0*I + 2.0)
    xi = np.sqrt(V_X2/(2.0*nu))
    Y2[k] = nu + 1/xi*(chi2_Mps[k] - nu)
    #Calculate corresponding p-values
    for l in range(N_rhos):
        p_values_Y2[k,l] = 1.0 - stats.chi2.cdf(Y2[k,l], 2)


#Constraints from Y^2 parameter estimation (page 815 numerical recipes)
#Number of degrees of freedom
nu_pe = 2.0
delta_pe_2sigma = stats.chi2.isf(0.045500264, nu_pe)
delta_pe_3sigma = stats.chi2.isf(0.002699796, nu_pe)
#Difference in Y^2 compared to smallest Y^2 value
deltaY2 = Y2 - np.min(Y2)
sigmas2_Y2_pe = np.full(N_Mps, rhos[N_rhos-1])
sigmas3_Y2_pe = np.full(N_Mps, rhos[N_rhos-1])
for k in range(N_Mps):
    for l in range(N_rhos):
        if ((deltaY2[k,l] > delta_pe_2sigma) and (rhos[l] < sigmas2_Y2_pe[k])):
            sigmas2_Y2_pe[k] = rhos[l]
        if ((deltaY2[k,l] > delta_pe_3sigma) and (rhos[l] < sigmas3_Y2_pe[k])):
            sigmas3_Y2_pe[k] = rhos[l]

#Smooth through interpolation
sigmas2_Y2_pe_interp = np.full(N_Mps, 10000.0)
sigmas3_Y2_pe_interp = np.full(N_Mps, 10000.0)
for k in range(N_Mps):
    poly_interp = Polynomial.fit(rhos, deltaY2[k], deg=3)
    x_interp, y_interp = poly_interp.linspace(n=200)
    sigma2_interp_density = (poly_interp - delta_pe_2sigma).roots()
    sigma3_interp_density = (poly_interp - delta_pe_3sigma).roots()
    for num in sigma2_interp_density:
        if ((np.isreal(num)) and (num > 0.0) and (num <= 10.0*rhos[N_rhos-1])):
            sigmas2_Y2_pe_interp[k] = num
            break
    if ((sigmas2_Y2_pe_interp[k] > 9999.0) and (np.min(y_interp) > delta_pe_2sigma)):
        sigmas2_Y2_pe_interp[k] = 0.0
    for num in sigma3_interp_density:
        if ((np.isreal(num)) and (num > 0.0)):
            sigmas3_Y2_pe_interp[k] = num
            break
    if ((sigmas3_Y2_pe_interp[k] > 9999.0) and (np.min(y_interp) > delta_pe_3sigma)):
        sigmas3_Y2_pe_interp[k] = 0.0

#Import YCG fig 6 and MRA fig 9 and Quinn fig 3
YCGfig6_2sigma_data = np.zeros((2,0))
with open("YCGfig6_2sigma.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        YCGfig6_2sigma_data = np.append(YCGfig6_2sigma_data, [[float(row[0])], [float(row[1])]], axis=1)
YCGfig6_3sigma_data = np.zeros((2,0))
with open("YCGfig6_3sigma.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        YCGfig6_3sigma_data = np.append(YCGfig6_3sigma_data, [[float(row[0])], [float(row[1])]], axis=1)
Quinn_fig3_updated2sigma_data = np.zeros((2,0))
with open("Quinn_fig3_updated2sigma.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        Quinn_fig3_updated2sigma_data = np.append(Quinn_fig3_updated2sigma_data, [[float(row[0])], [float(row[1])]], axis=1)
Quinn_fig3_old2sigma_data = np.zeros((2,0))
with open("Quinn_fig3_old2sigma.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        Quinn_fig3_old2sigma_data = np.append(Quinn_fig3_old2sigma_data, [[float(row[0])], [float(row[1])]], axis=1)
Quinn_fig3_Yoo2sigma_data = np.zeros((2,0))
with open("Quinn_fig3_Yoo2sigma.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        Quinn_fig3_Yoo2sigma_data = np.append(Quinn_fig3_Yoo2sigma_data, [[float(row[0])], [float(row[1])]], axis=1)
MRAfig9_25_data = np.zeros((2,0))
with open("MRAfig9_25.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        MRAfig9_25_data = np.append(MRAfig9_25_data, [[float(row[0])], [float(row[1])]], axis=1)
MRAfig9_50_data = np.zeros((2,0))
with open("MRAfig9_50.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        MRAfig9_50_data = np.append(MRAfig9_50_data, [[float(row[0])], [float(row[1])]], axis=1)
MRAfig9_100_data = np.zeros((2,0))
with open("MRAfig9_100.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        MRAfig9_100_data = np.append(MRAfig9_100_data, [[float(row[0])], [float(row[1])]], axis=1)

#Plot constraints
print('Minimum Y^2 =', np.min(Y2))
plt.plot(M_ps, sigmas2_Y2_pe_interp/0.009, label=r'$2\sigma$ interpolation', linestyle = 'solid', color='dodgerblue')
plt.plot(M_ps, sigmas3_Y2_pe_interp/0.009, label=r'$3\sigma$ interpolation', linestyle = 'dotted', color='dodgerblue')
plt.plot(YCGfig6_2sigma_data[0], YCGfig6_2sigma_data[1], label=r'Yoo et al. $2\sigma$', linestyle = 'solid', color='darkorange')
plt.plot(YCGfig6_3sigma_data[0], YCGfig6_3sigma_data[1], label=r'Yoo et al. $3\sigma$', linestyle = 'dotted', color='darkorange')
plt.plot(Quinn_fig3_updated2sigma_data[0], Quinn_fig3_updated2sigma_data[1], label=r'Quinn et al. $2\sigma$', linestyle = 'solid', color='blueviolet')
plt.plot(Quinn_fig3_old2sigma_data[0], Quinn_fig3_old2sigma_data[1], label=r'Quinn et al. old $2\sigma$', linestyle = 'dashed', color='blueviolet')
plt.plot(Quinn_fig3_Yoo2sigma_data[0], Quinn_fig3_Yoo2sigma_data[1], label=r'Quinn et al. Yoo $2\sigma$', linestyle = 'dashdot', color='blueviolet')
plt.plot(MRAfig9_25_data[0], MRAfig9_25_data[1], label=r'M-R&A 25 m.h-l.', linestyle = 'dashdot', color='forestgreen')
plt.plot(MRAfig9_50_data[0], MRAfig9_50_data[1], label=r'M-R&A 50 m.h-l.', linestyle = 'dashed', color='forestgreen')
plt.plot(MRAfig9_100_data[0], MRAfig9_100_data[1], label=r'M-R&A 100 m.h-l.', linestyle = 'solid', color='forestgreen')
plt.legend()
plt.xlabel(r'Perturber mass / $M_\odot$')
plt.ylabel(r'Halo density / $0.009M_\odot$pc$^{-3}$')
plt.xlim([1,10000])
#plt.ylim([rhos[0]*0.9/0.009, rhos[N_rhos-1]*1.1/0.009])
plt.ylim([0.08, 1.0])
plt.xscale('log')
plt.yscale('log')
plt.show()

#Plot best fit final distribution
for k,l in zip([np.argmin(chi2_Mps)//N_rhos], [np.argmin(chi2_Mps)%N_rhos]):
    print(r'Best fit for M_p = {}$M_\odot$, $\rho$={}$M_\odot$pc$^{{-3}}$: $\alpha$={}, $A$={}'.format(M_ps[k], rhos[l], str(alphas[k,l])[:4], str(As_Mp[k,l])[:4]))
    plt.plot(a_bins*np.exp(0.5*da)/au, 100000/(np.sum(N_input_dist_r_proj[k,l]))*initial_semimajor_axis_dist(a_bins, alphas[k,l], As_Mp[k,l])*das*np.size(a_MRA)/np.sum((initial_semimajor_axis_dist(a_bins, alphas[k,l], As_Mp[k,l])*das)[i_real_min:i_real_max+1]), label='Initial', linestyle='-', color='tab:orange')
    plt.plot(a_bins*np.exp(0.5*da)/au, P_minchi2[k,l], label='Final', linestyle='-', color='tab:green')
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
plt.scatter(a_bins_condensed*np.exp(0.5*(np.log(a_max/a_min)/np.size(a_bins_condensed)))/au, N_real_condensed/reduction_factor, label='Observed', marker='x', color='tab:blue')
#plt.scatter(s_bins_condensed*np.exp(0.5*(np.log(s_max/s_min)/np.size(s_bins_condensed))), N_real_condensed/reduction_factor, label='Real binary distribution')
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.legend()
plt.xlabel('Projected separation / au')
#plt.xlabel('Angular separation, arcsec')
plt.ylabel('Number of binaries')
plt.tight_layout()
plt.show()

print('Plotting p-value contour plot')
plt.figure(figsize=(7.5,5))
ax = plt.gca()
cs = ax.contourf(M_ps, rhos/0.009, np.transpose(p_values_Y2), cmap='Blues')
plt.scatter([M_ps[np.argmin(chi2_Mps)//N_rhos]], [rhos[np.argmin(chi2_Mps)%N_rhos]/0.009], marker='x', color='white')
plt.plot(M_ps, sigmas2_Y2_pe_interp/0.009, linestyle = 'solid', color='black')
plt.plot(M_ps, sigmas3_Y2_pe_interp/0.009, linestyle = 'dotted', color='black')
clb = plt.colorbar(cs)
plt.xlabel(r'Perturber mass / $M_\odot$')
plt.ylabel(r'Halo density / $0.009M_\odot$pc$^{-3}$')
#plt.title(r'Best fit p-values, according to $Y^2$')
plt.xscale('log')
plt.yscale('log')
plt.ylim([0.056, 1.34])
plt.tight_layout()
plt.show()

print('Plotting best fit alpha contour plot')
plt.figure(figsize=(7.5,5))
ax = plt.gca()
cs = ax.contourf(M_ps, rhos/0.009, np.transpose(np.where(chi2_Mps < chi2_Mp0, alphas, alpha_0)), cmap='Blues')
plt.scatter([M_ps[np.argmin(chi2_Mps)//N_rhos]], [rhos[np.argmin(chi2_Mps)%N_rhos]/0.009], marker='x', color='white')
clb = plt.colorbar(cs)
plt.xlabel(r'Perturber mass / $M_\odot$')
plt.ylabel(r'Halo density / $0.009M_\odot$pc$^{-3}$')
#plt.title(r'Best fit $\alpha$, according to $\chi^2$')
plt.xscale('log')
plt.yscale('log')
plt.ylim([0.056, 1.34])
plt.tight_layout()
plt.show()

print('Plotting best fit A contour plot')
plt.figure(figsize=(7.5,5))
ax = plt.gca()
cs = ax.contourf(M_ps, rhos/0.009, np.transpose(np.where(chi2_Mps < chi2_Mp0, As_Mp, As_Mp0)), cmap='Blues')
plt.scatter([M_ps[np.argmin(chi2_Mps)//N_rhos]], [rhos[np.argmin(chi2_Mps)%N_rhos]/0.009], marker='x', color='white')
clb = plt.colorbar(cs)
plt.xlabel(r'Perturber mass / $M_\odot$')
plt.ylabel(r'Halo density / $0.009M_\odot$pc$^{-3}$')
#plt.title(r'Best fit A, according to $\chi^2$')
plt.xscale('log')
plt.yscale('log')
plt.ylim([0.056, 1.34])
plt.tight_layout()
plt.show()

#print('N_binaries =', np.size(s_MRA), ', Seed =', seed, ', s_analysis_not_r_proj?', s_analysis_not_r_proj, ', Best fit M_p =', M_ps[np.argmin(chi2_Mps)//N_rhos], ', Best fit rho =', rhos[np.argmin(chi2_Mps)%N_rhos])

print('Finished')