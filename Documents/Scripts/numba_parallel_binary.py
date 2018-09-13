import os
os.system("python setup.py build_ext --inplace")
import datetime
import math
import numpy as np
from numba import cuda, config
from numba import autojit, jit, int64, float32, int32
from numba.cuda.random import create_xoroshiro128p_states, xoroshiro128p_uniform_float64
from encounters import encounterRate
from scipy.constants import au, year, mega, giga, parsec, kilo, G
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

print('GPU:', cuda.gpus)
print('Simulator enabled =', config.ENABLE_CUDASIM)
print('Initialising')

#INTERNAL UNITS
length_unit = 1.0 * au
mass_unit = 2.0*10.0**30.0
#Make new G=1
#https://gandalfcode.github.io/gandalf-school/Units.pdf
G_new = 1.0
time_unit = np.sqrt(G_new*length_unit**3.0/(G*mass_unit))
#print('time_unit =', time_unit)
G_new = G*mass_unit*time_unit**2.0/(length_unit**3.0)
#print('G_new =', G_new)

#REWRITTEN PYTHON FUNCTIONS
#Draw random number between 0.0 and 1.0
@cuda.jit(device=True)
def rand(rng_states):
        return xoroshiro128p_uniform_float64(rng_states, cuda.grid(1))

#Draw impact parameter
@cuda.jit(device=True)
def draw_b(b_max, rng_states):   
        return b_max * math.sqrt(rand(rng_states))

#Calculate b_max
@cuda.jit(device=True)
def calc_b_max(M_p, v_rel, a, m1, m2, delta, G_new):
        return (64.0*G_new*M_p**2.0*a**3.0/((m1+m2)*v_rel**2.0*delta**2.0))**0.25

#Functions for draw_vmaxwellian
@cuda.jit(device=True)
def draw_vmaxwellian(v_rms, v_min, v_max, rng_states):
        #Total area under comparison function
        area = (v_max - v_min)*3.0**(3.0/2.0)/(2.0*v_rms)*math.exp(-3.0/2.0)
        while True:
                u = area*rand(rng_states)
                x = vmaxwellianX_from_area(u, v_rms, v_min)
                y = vmaxwellianComparison(x, v_rms, v_min, v_max)*rand(rng_states)
                if y < vmaxwellianPdf(x, v_rms):
                        return x

@cuda.jit(device=True)
def vmaxwellianPdf(x, v_rms):
        return x**3.0/(2.0*v_rms**4.0)*math.exp(-x**2.0/(2.0*v_rms**2.0))

@cuda.jit(device=True)
def vmaxwellianComparison(x, v_rms, v_min, v_max):
        answer = 0.0
        if v_min<x<v_max:
                answer = 3.0**(3.0/2.0)/(2.0*v_rms)*math.exp(-3.0/2.0)
        return answer
        
@cuda.jit(device=True)
def vmaxwellianX_from_area(A, v_rms, v_min):
        return v_min + 2.0*v_rms*A*math.exp(3.0/2.0)/3.0**(3.0/2.0)

#Functions for impulse encounter
@cuda.jit(device=True)
def findEccentricAnomaly(e, M):
        #Solves Kepler's equation E-esinE=M to find eccentric anomaly E given eccentricity e and mean anomaly M.
        #From page 36 of Solar System Dynamics, or Danby 1988
        E = M + math.copysign(1.0, math.sin(M))*0.85*e
        count = 0    
        while abs(E - e*math.sin(E) - M) > 10.0**(-8.0):
                f = E - e*math.sin(E) - M
                f_p = 1.0 - e*math.cos(E)
                f_pp = e*math.sin(E)
                f_ppp = e*math.cos(E)

                d_1 = - f/f_p
                d_2 = -f/(f_p + 0.5*d_1*f_pp)
                d_3 = -f/(f_p + 0.5*d_2*f_pp + d_2**2.0*f_ppp/6.0)
        
                E += d_3
                
                count += 1
        return E

@cuda.jit(device=True)
def randomDirection(rng_states):        
        u = 2.0*rand(rng_states) - 1.0
        theta = 2.0*math.pi*rand(rng_states)
        v1 = math.sqrt(1.0 - u**2.0)*math.cos(theta)
        v2 = math.sqrt(1.0 - u**2.0)*math.sin(theta)
        v3 = u    
        v_norm = math.sqrt(v1**2.0 + v2**2.0 + v3**2.0)
        v1 /= v_norm
        v2 /= v_norm
        v3 /= v_norm
        return v1, v2, v3
        
@cuda.jit(device=True)
def setupRandomBinary(a, e, m1, m2, G_new, rng_states):
        #pos = cuda.grid(1)
        #Randomise mean anomaly
        M = 2.0*math.pi*rand(rng_states)
        #print('M =', pos, M)
        #Find eccentric anomaly from Kepler's equation
        E = findEccentricAnomaly(e, M)
        #print('E =', pos, E)
        #print('E-esinE=M =', pos, E-e*math.sin(E)-M)
        #Find true anomaly
        f = 2.0*math.atan(((1.0+e)/(1.0-e))**0.5 * math.tan(0.5*E))
        #print('f =', pos, f)
        
        #Randomise true anomaly
        #f = np.random.uniform(0.0, 2.0*np.pi)
        
        #Initial separation
        r = a*(1.0 - e**2.0)/(1.0 + e*math.cos(f))    
        #print('r/pc =', pos, r)    
        #Mean motion
        #print('G_new =', pos, G_new)
        n = math.sqrt(G_new*(m1+m2)/(a**3.0))
        #print('n =', pos, n)
        #Initial coordinates of first star (cartesian)
        X11 = 0.0
        X12 = 0.0
        X13 = 0.0
        #Initial velocity of first star
        X31 = 0.0
        X32 = 0.0
        X33 = 0.0
        #Initial coordinates of second star (cartesian)
        X21 = r*math.cos(f)
        X22 = r*math.sin(f)
        X23 = 0.0
        #Initial velocity of second star
        X41 = - n*a/(math.sqrt(1.0-e**2.0)) * math.sin(f)
        X42 = n*a/(math.sqrt(1.0-e**2.0)) * (e + math.cos(f))
        X43 = 0.0
        #Move into centre of mass rest frame
        #x
        R1 = (m1*X11 + m2*X21)/(m1 + m2)
        V1 = (m1*X31 + m2*X41)/(m1 + m2)
        #X11 -= R1
        #X21 -= R1
        X31 -= V1
        X41 -= V1
        #y
        R2 = (m1*X12 + m2*X22)/(m1 + m2)
        V2 = (m1*X32 + m2*X42)/(m1 + m2)
        #X12 -= R2
        #X22 -= R2
        X32 -= V2
        X42 -= V2
        #z
        R3 = (m1*X13 + m2*X23)/(m1 + m2)
        V3 = (m1*X33 + m2*X43)/(m1 + m2)
        #X13 -= R3
        #X23 -= R3
        X33 -= V3
        X43 -= V3
        #print('R =', pos, R1, R2, R3)
        #print('V =', pos, V1, V2, V3)
        return X11, X12, X13, X21, X22, X23, X31, X32, X33, X41, X42, X43

@cuda.jit(device=True)
def impactAndVelocityVectors(b, v, rng_states):
        #Velocity vector
        v1, v2, v3 = randomDirection(rng_states)
        v1 *= v
        v2 *= v
        v3 *= v
        #Random vector
        n1, n2, n3 = randomDirection(rng_states)
        #Unnormalised impact parameter vector
        b1 = v2*n3 - v3*n2
        b2 = v3*n1 - v1*n3
        b3 = v1*n2 - v2*n1
        #Normalise b
        b_norm = (b1**2.0 + b2**2.0 + b3**2.0)**0.5
        b1 *= b/b_norm
        b2 *= b/b_norm
        b3 *= b/b_norm
        return b1, b2, b3, v1, v2, v3

@cuda.jit(device=True)
def orbitalElements(X11, X12, X13, X21, X22, X23, X31, X32, X33, X41, X42, X43, m1, m2, G_new):
        #pos = cuda.grid(1)
        R = math.sqrt((X11 - X21)**2.0 + (X12 - X22)**2.0 + (X13 - X23)**2.0)
        V = math.sqrt((X31 - X41)**2.0 + (X32 - X42)**2.0 + (X33 - X43)**2.0)
        #print('R =', pos, R)
        #print('V =', pos, V)
        #Total energy
        E = m1*m2*(V**2.0/(2.0*(m1+m2)) - G_new/R) 
        #print('E =', pos, E)      
        #Total angular momentum
        L_squared = (m1*m2/(m1+m2))**2.0*(((X12 - X22)*(X33 - X43)-(X13 - X23)*(X32 - X42))**2.0 + ((X13 - X23)*(X31 - X41)-(X11 - X21)*(X33 - X43))**2.0 + ((X11 - X21)*(X32 - X42)-(X12 - X22)*(X31 - X41))**2.0)
        #print('L_squared =', pos, L_squared)
        #Semi-major axis
        a = G_new*m1*m2/(2.0*abs(E))  
        #print('a =', pos, a)    
        #Eccentricity
        e = math.sqrt(1.0 + 2.0*(m1+m2)*L_squared*E/(G_new**2.0*(m1*m2)**3.0))     
        #print('e =', pos, e)
        return((E >= 0.0), a, e)
                
@cuda.jit(device=True)
def impulseEncounter(m1, m2, v, b, a, e, M_p, G_new, rng_states):
        #print('m1 =', m1)
        pos = cuda.grid(1)
        #print('v =', pos, v)
        #print('b =', pos, b)
        #print('a =', pos, a)
        #print('e =', pos, e)
        #90 degree deflection radius
        b_90_1 = G_new*(M_p+m1)/v**2.0   
        b_90_2 = G_new*(M_p+m2)/v**2.0
        #print('b_90_1/pc =', pos, b_90_1)
        #print('b_90_2/pc =', pos, b_90_2)
        #Open binary
        #print('G_new =', G_new)
        X11, X12, X13, X21, X22, X23, X31, X32, X33, X41, X42, X43 = setupRandomBinary(a, e, m1, m2, G_new, rng_states)
        #print('Before encounter x1 =', pos, X11, X12, X13) 
        #print('Before encounter x2 =', pos, X21, X22, X23)
        #print('Before encounter v1 =', pos, X31, X32, X33)
        #print('Before encounter v2 =', pos, X41, X42, X43)
        #Find impact parameter vector and velocity vector
        b1, b2, b3, v1, v2, v3 = impactAndVelocityVectors(b, v, rng_states)
        #print('b_vec =', pos, b1, b2, b3)
        #print('v_vec =', pos, v1, v2, v3)
        
        #Implement encounter for star 1
        #Calculate impact parameter for this star
        #x
        b_star1 = (X11*v1 + X12*v2 + X13*v3)/v**2.0 * v1 + b1 - X11
        #y
        b_star2 = (X11*v1 + X12*v2 + X13*v3)/v**2.0 * v2 + b2 - X12
        #z
        b_star3 = (X11*v1 + X12*v2 + X13*v3)/v**2.0 * v3 + b3 - X13
        #print('b_star =', pos, b_star1, b_star2, b_star3)
        b_star_norm = math.sqrt(b_star1**2.0 + b_star2**2.0 + b_star3**2.0)
        #Velocity changes
        #x
        #b direction
        X31 += 2.0*M_p*v/(m1+M_p) * (b_star_norm/b_90_1)/(1.0 + b_star_norm**2.0/b_90_1**2.0) * (b_star1/b_star_norm)
        #-v direction
        X31 += 2.0*M_p*v/(m1+M_p) * 1.0/(1.0 + b_star_norm**2.0/b_90_1**2.0) * (-v1/v)
        #y
        #b direction
        X32 += 2.0*M_p*v/(m1+M_p) * (b_star_norm/b_90_1)/(1.0 + b_star_norm**2.0/b_90_1**2.0) * (b_star2/b_star_norm)
        #-v direction
        X32 += 2.0*M_p*v/(m1+M_p) * 1.0/(1.0 + b_star_norm**2.0/b_90_1**2.0) * (-v2/v)
        #z
        #b direction
        X33 += 2.0*M_p*v/(m1+M_p) * (b_star_norm/b_90_1)/(1.0 + b_star_norm**2.0/b_90_1**2.0) * (b_star3/b_star_norm)
        #-v direction
        X33 += 2.0*M_p*v/(m1+M_p) * 1.0/(1.0 + b_star_norm**2.0/b_90_1**2.0) * (-v3/v)
        
        #Implement encounter for star 2
        #Calculate impact parameter for this star
        #x
        b_star1 = (X21*v1 + X22*v2 + X23*v3)/v**2.0 * v1 + b1 - X21
        #y
        b_star2 = (X21*v1 + X22*v2 + X23*v3)/v**2.0 * v2 + b2 - X22
        #z
        b_star3 = (X21*v1 + X22*v2 + X23*v3)/v**2.0 * v3 + b3 - X23
        #print('b_star =', pos, b_star1, b_star2, b_star3)
        b_star_norm = math.sqrt(b_star1**2.0 + b_star2**2.0 + b_star3**2.0)
        #Velocity changes
        #x
        #b direction
        X41 += 2.0*M_p*v/(m1+M_p) * (b_star_norm/b_90_2)/(1.0 + b_star_norm**2.0/b_90_2**2.0) * (b_star1/b_star_norm)
        #-v direction
        X41 += 2.0*M_p*v/(m1+M_p) * 1.0/(1.0 + b_star_norm**2.0/b_90_2**2.0) * (-v1/v)
        #y
        #b direction
        X42 += 2.0*M_p*v/(m1+M_p) * (b_star_norm/b_90_2)/(1.0 + b_star_norm**2.0/b_90_2**2.0) * (b_star2/b_star_norm)
        #-v direction
        X42 += 2.0*M_p*v/(m1+M_p) * 1.0/(1.0 + b_star_norm**2.0/b_90_2**2.0) * (-v2/v)
        #z
        #b direction
        X43 += 2.0*M_p*v/(m1+M_p) * (b_star_norm/b_90_2)/(1.0 + b_star_norm**2.0/b_90_2**2.0) * (b_star3/b_star_norm)
        #-v direction
        X43 += 2.0*M_p*v/(m1+M_p) * 1.0/(1.0 + b_star_norm**2.0/b_90_2**2.0) * (-v3/v)

        #print('HERE')
        #print(str.format('{0:.15f}', X11))
        #print('After encounter x1 =', pos, X11, X12, X13) 
        #print('After encounter x2 =', pos, X21, X22, X23)
        #print('After encounter v1 =', pos, X31, X32, X33)
        #print('After encounter v2 =', pos, X41, X42, X43)        
        #Close binary
        notBound, a_new, e_new = orbitalElements(X11, X12, X13, X21, X22, X23, X31, X32, X33, X41, X42, X43, m1, m2, G_new)       
        return notBound, a_new, e_new


#SETUP ENCOUNTERS
#Mass of binary stars IN SOLAR MASSES
m1 = 0.5 * 2.0*10.0**30.0 / mass_unit
m2 = 0.5 * 2.0*10.0**30.0 / mass_unit
#Mass of perturbers IN SOLAR MASSES
M_p = 1000.0 * 2.0*10.0**30.0 / mass_unit
#Velocity dispersion of Maxwellian velocity distribution in internal units
v_rel = 220.0 * 1000.0 * time_unit/length_unit
#Density of dark matter halo solar masses/length_unit**3
rho = 0.009 * (length_unit/parsec)**3.0
#Number density of perturbers
n_p = rho/M_p
#Time to run simulation for
T = 10.0**10.0 * year /time_unit
#Number of binary pairs
N_bin = 10**6
#Tidal radius in au
a_T = 1.0 * parsec/length_unit
#print('a_T,au =', a_T)
#Delta
delta = 10.0**(-3.0)
#Semi-major axis array
#Draw a from distribution dN/dloga\propto a^{1-alpha}
alpha = 1.0
a_min = 10.0**3.0*au/length_unit
a_max = np.min([10.0**6.0*au/length_unit, a_T])
if alpha == 2.0:
        c = np.log(a_min)/np.log(a_max/a_min)
        a = (a_max/a_min)**(np.random.random(N_bin) + c)
else:
        a = (np.random.random(N_bin)*(a_max**(2.0-alpha) - a_min**(2.0-alpha)) + a_min**(2.0-alpha))**(1.0/(2.0-alpha))
#Draw e from distribution uniform in e^2 between 0 and 1
e = (np.random.random(N_bin))**(1.0/3.0)
#Number of binaries broken
N_broken = 0

#SETUP
#Minimum impact parameter
b_min = 0.0
#Maximum maximum impact parameter
b_max_max = (64.0*G_new*M_p**2.0*a_T**3.0/((m1+m2)*v_rel**2.0*delta**2.0))**0.25
#print('b_max_max/pc =', b_max_max)
#Minimum velocity
v_min = 10.0**(-2.0)*v_rel
#Maximum velocity
v_max = 10.0**2.0*v_rel

@cuda.jit
def run_encounters(a, e, m1, m2, M_p, v_rel, a_T, b_max_max, v_min, v_max, N_enc, delta, G_new, N_broken, rng_states):
    #Position in a and e arrays to update
    pos = cuda.grid(1)
    if pos < len(a):
        #print('N_enc[pos] =', N_enc[pos])
        #print('G_new =', G_new)
        #print('m1 =', m1, 'm2 =', m2)
        #print('delta =', delta)
        #print('M_p =', M_p)
        #print('v_rel =', v_rel)
        #print('a[pos] =', a[pos])
        #print('e[pos] =', e[pos])
        #Maximum impact parameter
        b_max = calc_b_max(M_p, v_rel, a[pos], m1, m2, delta, G_new)
        #print('b_max/pc =', pos, b_max)
        #Implement encounters
        #print('Number of encounters =', pos,  N_enc[pos])
        for j in range(N_enc[pos]):
            #print('N_enc[pos] =', pos, N_enc[pos])
            #Impact parameter of encounter
            b = draw_b(b_max, rng_states)
            #print(b/(3.086*10.0**16.0))
            #print(b/(3.086*10.0**16.0))
            if (b<=b_max):
                #print('Encounter', pos)
                #print('b =', pos, b)
                #Draw relative velocity of encounter
                v = draw_vmaxwellian(v_rel, v_min, v_max, rng_states)
                #print('v =', pos, v)
                #Implement encounter
                notBound, a[pos], e[pos] = impulseEncounter(m1, m2, v, b, a[pos], e[pos], M_p, G_new, rng_states)
                if (notBound or a[pos]>=a_T):
                    N_broken += 1
                    a[pos] = -1.0
                    e[pos] = -1.0
                #Update maximum impact parameter
                b_max = calc_b_max(M_p, v_rel, a[pos], m1, m2, delta, G_new)
                #print('a/pc =', pos, a[pos])
                #print('e =', pos, e[pos])
                #print('b_max/pc =', pos, b_max)
            N_enc[pos] -= 1
        cuda.syncthreads()

def plotCumulative(x):
	x_min = np.min(x)
	x_max = np.max(x)
	N_bin = 1000
	dx = (x_max - x_min)/(N_bin)
	x_bins = np.array([x_min + i*dx for i in range(N_bin)])
	N = np.zeros(N_bin)
	for i in range(N_bin):
		N[i] = np.size(np.where(x<=x_bins[i]))
	#Normalise
	N /= np.size(x)
	#Plot
	plt.plot(x_bins, N)
	plt.show()
	return

'''
#Test distributions
N=1000
x=np.zeros((N,3))
@cuda.jit
def test_distributions(x, rng_states):
	pos = cuda.grid(1)
	if pos < len(x[:,0]):
		x[pos,0], x[pos,1], x[pos,2], p, q, r = impactAndVelocityVectors(10.0, 50.0, rng_states)
	cuda.syncthreads()
	return

threadsperblock = 32
blockspergrid = math.ceil(x[:,0].size / threadsperblock)
rng_states = create_xoroshiro128p_states(threadsperblock * blockspergrid, seed=1)
test_distributions[blockspergrid, threadsperblock](x, rng_states)
print('x =', x)
#plotCumulative(x)
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
ax1.scatter(x[:,0], x[:,1], x[:,2])
plt.show()
'''



#print('a_initial/au =', a)
#print('e_initial =', e)
threadsperblock = 128
blockspergrid = math.ceil(a.size / threadsperblock)
#Call multiple times to avoid timeout
#Time per run:
t_run = 10.0**4 * year / time_unit
if T > t_run:
    #Number of runs
    N_runs = int(T // t_run)
    #Adjust t_run:
    t_run = T / N_runs
else:
    N_runs = 1
    t_run = T
#print('N_runs =', N_runs)
#Mean number of encounters in time t_run
N_mean = t_run*encounterRate(n_p, v_rel, b_min, b_max_max, v_min, v_max)
print('Evolving population, start time:', datetime.datetime.now())
for i in range(N_runs):
    #Number of encounters
    N_enc = np.random.poisson(N_mean, size=N_bin)
    rng_states = create_xoroshiro128p_states(threadsperblock * blockspergrid, seed=1)
    with cuda.profiling():
        run_encounters[blockspergrid, threadsperblock](a, e, m1, m2, M_p, v_rel, a_T, b_max_max, v_min, v_max, N_enc, delta, G_new, N_broken, rng_states)
    e=e[np.where(a>0.0)]
    a=a[np.where(a>0.0)]
#N_broken = N_bin - np.size(a)
#print('a_final/au =', a)
#print('e_final =', e)
print('Finish time:', datetime.datetime.now())
print('N_broken =', N_broken)

print('Saving')
np.savez('simulation_GPU_{}Msol_{}e{}.npz'.format(int(M_p*mass_unit/(2.0*10.0**30.0)), int(N_bin/10**int(np.floor(np.log10(N_bin)))), int(np.floor(np.log10(N_bin)))), a_fin=a*length_unit, e_fin=e, N_broken=N_broken)
print('Finished')

'''
an_array = np.arange(10).reshape(2,5).astype(np.float32)
print(an_array)

increment = 2
cuda.to_device(increment)

@cuda.jit
def increment_a_2D_array(an_array):
        pos = cuda.grid(2)
        if pos[0]<an_array.shape[0] and pos[1]<an_array.shape[1]:
                an_array[pos] += increment
        
threadsperblock = (16, 16)
blockspergrid_x = math.ceil(an_array.shape[0] / threadsperblock[0])
blockspergrid_y = math.ceil(an_array.shape[1] / threadsperblock[1])
blockspergrid = (blockspergrid_x, blockspergrid_y)
increment_a_2D_array[blockspergrid, threadsperblock](an_array)

print(an_array)
'''





























