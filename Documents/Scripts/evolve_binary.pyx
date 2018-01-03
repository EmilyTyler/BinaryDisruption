import math
import numpy as np
cimport numpy as np
import itertools as it
from scipy.constants import G

from orbital_elements import orbitalElements
from random_binary import setupRandomBinary

#Evolution of binary through integration
def integrateBinary(int N, np.ndarray[double, ndim=2] X, np.ndarray[double, ndim=1] m):
        # dt is step size
        # X = [x1, x2, x3,..., v1, v2, v3,...]
        # N is number of objects
        # m = [m1, m2, m3,...]
        X_0 = X[:N]
        V_0 = X[N:]      
        #Hermite scheme (Dehnen and Read 2011) with Aarseth criterion timesteps
        #Find initial acceleration and jerk and timestep
        (A_0, J_0, dt) = acc_jerk_and_timestep(N, X_0, V_0, m)
        #Predict positions and velocities
        X_p = X_0 + V_0*dt + A_0*dt**2.0/2.0 + J_0*dt**3.0/6.0
        V_p = V_0 + A_0*dt + J_0*dt**2.0/2.0                   
        #Estimate acceleration and jerk
        (A_1, J_1) = acc_and_jerk(N, X_p, V_p, m)     
        #Obtain corrected positions and velocities
        X_1 = X_0 + (V_p+V_0)*dt/2.0 + (A_0-A_1)*dt**2.0/12.0
        V_1 = V_0 + (A_1+A_0)*dt/2.0 + (J_0-J_1)*dt**2.0/12.0      
        return (np.concatenate([X_1, V_1]), dt)

#See Nitadori and Makino 2008 for higher derivatives of acceleration
def calc_alpha(int N, np.ndarray[double, ndim=3] x, np.ndarray[double, ndim=3] v):
        cdef np.ndarray alpha = np.zeros((N, N), dtype=float)
        for i in range(N):
                for j in range(N):
                        alpha[i,j] = np.dot(x[i,j], v[i,j])/np.dot(x[i,j], x[i,j])
        return alpha

def calc_beta(int N, np.ndarray[double, ndim=3] x, np.ndarray[double, ndim=3] v, np.ndarray[double, ndim=3] a, np.ndarray[double, ndim=2] alpha):
        cdef np.ndarray beta = np.zeros((N, N), dtype=float)
        for i in range(N):
                for j in range(N):
                        beta[i,j] = (np.dot(v[i,j], v[i,j]) + np.dot(x[i,j], a[i,j]))/np.dot(x[i,j], x[i,j]) + alpha[i,j]**2.0
        return beta

def calc_gamma(int N, np.ndarray[double, ndim=3] x, np.ndarray[double, ndim=3] v, np.ndarray[double, ndim=3] a, np.ndarray[double, ndim=3] j, np.ndarray[double, ndim=2] alpha, np.ndarray[double, ndim=2] beta):
        cdef np.ndarray gamma = np.zeros((N,N), dtype=float)
        for i in range(N):
                for k in range(N):
                        gamma[i,k] = (3.0*np.dot(v[i,k], a[i,k]) + np.dot(x[i,k], j[i,k]))/np.dot(x[i,k], x[i,k]) + alpha[i,k]*(3.0*beta[i,k] - 4.0*alpha[i,k]**2.0)
        return gamma
        
def calc_acc(int N, np.ndarray[double, ndim=2] X, np.ndarray[double, ndim=1] m):
        cdef np.ndarray acc = np.zeros((N,3), dtype=float)
        cdef np.ndarray x = np.zeros((N,N,3), dtype=float)
        cdef np.ndarray x_norm = np.zeros((N,N), dtype=float)
        cdef np.ndarray A = np.zeros((N,N,3), dtype=float)
        cdef double xij_norm
        for i in range(N):
                for j in range(N):
                        x[i,j] = X[j] - X[i]
                        x_norm[i,j] = np.sqrt(x[i,j,0]**2.0+x[i,j,1]**2.0+x[i,j,2]**2.0)
                        A[i,j] = m[j]*x[i,j]/(x_norm[i,j])**3.0
        for i in range(N):
                for j in it.chain(range(i), range(i+1, N)):
                        acc[i] += A[i,j]
        return (G*acc, x, x_norm, A)

def calc_jerk(int N, np.ndarray[double, ndim=3] x, np.ndarray[double, ndim=2] x_norm, np.ndarray[double, ndim=2] V, np.ndarray[double, ndim=1] m, np.ndarray[double, ndim=3] A):
        cdef np.ndarray jerk = np.zeros((N,3), dtype=float)
        cdef np.ndarray v = np.zeros((N,N,3), dtype=float)
        cdef np.ndarray J = np.zeros((N,N,3), dtype=float)
        for i in range(N):
                for j in range(N):
                        v[i,j] = V[j] - V[i]
        cdef np.ndarray alpha = calc_alpha(N, x, v)
        for i in range(N):
                for j in range(N):
                        J[i,j] = m[j]*v[i,j]/(x_norm[i,j])**3.0 - 3.0*alpha[i,j]*A[i,j]     
        for i in range(N):
                for j in it.chain(range(i), range(i+1, N)):
                        jerk[i] += J[i,j]
        return (G*jerk, v, J, alpha)
        
def calc_snap(int N, np.ndarray[double, ndim=3] x, np.ndarray[double, ndim=2] x_norm, np.ndarray[double, ndim=3] v,  np.ndarray[double, ndim=1] m, np.ndarray[double, ndim=3] A, np.ndarray[double, ndim=3] J, np.ndarray[double, ndim=2] alpha, np.ndarray[double, ndim=2] acc):
        cdef np.ndarray snap = np.zeros((N,3), dtype=float)
        cdef np.ndarray a = np.zeros((N,N,3), dtype=float)
        cdef np.ndarray S = np.zeros((N,N,3), dtype=float)
        for i in range(N):
                for j in range(N):
                        a[i,j] = acc[j] - acc[i]
        cdef np.ndarray beta = calc_beta(N, x, v, a, alpha)
        for i in range(N):
                for j in range(N):
                        S[i,j] = m[j]*a[i,j]/(x_norm[i,j])**3.0 - 6.0*alpha[i,j]*J[i,j] - 3.0*beta[i,j]*A[i,j]
        for i in range(N):
                for j in it.chain(range(i), range(i+1, N)):
                        snap[i] += S[i,j]
        return (G*snap, a, S, beta)

def calc_crackle(int N, np.ndarray[double, ndim=3] x, np.ndarray[double, ndim=2] x_norm, np.ndarray[double, ndim=3] v, np.ndarray[double, ndim=3] a, np.ndarray[double, ndim=1] m, np.ndarray[double, ndim=3] A, np.ndarray[double, ndim=3] J, np.ndarray[double, ndim=3] S, np.ndarray[double, ndim=2] alpha, np.ndarray[double, ndim=2] beta, np.ndarray[double, ndim=2] jerk):
        cdef np.ndarray crackle = np.zeros((N,3), dtype=float)
        cdef np.ndarray j = np.zeros((N,N,3), dtype=float)
        cdef np.ndarray C = np.zeros((N,N,3), dtype=float)
        for i in range(N):
                for k in range(N):
                        j[i,k] = jerk[k] - jerk[i]
        cdef np.ndarray gamma = calc_gamma(N, x, v, a, j, alpha, beta)
        for i in range(N):
                for k in range(N):
                        C[i,k] = m[k]*j[i,k]/(x_norm[i,k])**3.0 - 9.0*alpha[i,k]*S[i,k] - 9.0*beta[i,k]*J[i,k] - 3.0*gamma[i,k]*A[i,k]
        for i in range(N):
                for k in it.chain(range(i), range(i+1, N)):
                        crackle[i] += C[i,k]
        return (G*crackle)
        
#Aarseth criterion
def timestep(double eta, int N, np.ndarray[double, ndim=2] acc, np.ndarray[double, ndim=2] jerk, np.ndarray[double, ndim=2] snap, np.ndarray[double, ndim=2] crackle):
        cdef np.ndarray timesteps = np.zeros(N, dtype=float)
        cdef double acc_norm, jerk_norm, snap_norm, crackle_norm
        for i in range(N):
                acc_norm = np.sqrt(acc[i,0]**2.0+acc[i,1]**2.0+acc[i,2]**2.0)
                jerk_norm = np.sqrt(jerk[i,0]**2.0+jerk[i,1]**2.0+jerk[1,2]**2.0)
                snap_norm = np.sqrt(snap[i,0]**2.0+snap[i,1]**2.0+snap[i,2]**2.0)
                crackle_norm = np.sqrt(crackle[i,0]**2.0+crackle[i,1]**2.0+crackle[i,2]**2.0)
                timesteps[i] = np.sqrt(eta*(acc_norm*snap_norm + jerk_norm**2.0)/(jerk_norm*crackle_norm + snap_norm**2.0))
        return np.min(timesteps)

def acc_jerk_and_timestep(int N, np.ndarray[double, ndim=2] X, np.ndarray[double, ndim=2] V, np.ndarray[double, ndim=1] m):
        cdef np.ndarray acc, x, x_norm, A, jerk, v, J, alpha, snap, a, S, beta, crackle
        #Find acceleration
        (acc, x, x_norm, A) = calc_acc(N, X, m)
        #Find jerk
        (jerk, v, J, alpha) = calc_jerk(N, x, x_norm, V, m, A)
        #Find snap
        (snap, a, S, beta) = calc_snap(N, x, x_norm, v, m, A, J, alpha, acc)
        #Find crackle
        crackle = calc_crackle(N, x, x_norm, v, a, m, A, J, S, alpha, beta, jerk)
        #Find timestep
        cdef double dt = timestep(0.02, N, acc, jerk, snap, crackle)       
        return (acc, jerk, dt)

def acc_and_jerk(int N, np.ndarray[double, ndim=2] X, np.ndarray[double, ndim=2] V, np.ndarray[double, ndim=1] m):
        cdef np.ndarray acc, x, x_norm, A, jerk, v, J, alpha
        #Find acceleration
        (acc, x, x_norm, A) = calc_acc(N, X, m)
        #Find jerk
        (jerk, v, J, alpha) = calc_jerk(N, x, x_norm, V, m, A)
        return (acc, jerk)
             
                
