import math
import numpy as np
cimport numpy as np
import itertools as it

from orbital_elements import orbitalElements
from random_binary import setupRandomBinary

G = 6.67 * 10.0**(-11.0)



#Evolution of binary through analytics
def analyticBinary(np.ndarray[double, ndim=2] X, double m1, double m2):
        #Find orbital elements
        (a, e) = orbitalElements(X, m1, m2)
        #Rebuild binary
        return setupRandomBinary(a, e, m1, m2)

cdef double a1, a2, j1, j2, x1_p, v1_p, x2_p, v2_p, a1_1, a2_1, j1_1, j2_1, x1_1, x2_1, v1_1, v2_1


#Evolution of binary through integration
def integrateBinary(int N, np.ndarray[double, ndim=2] X, np.ndarray[double, ndim=1] m, double dt):
        # dt is step size
        # X = [x1, x2, x3,..., v1, v2, v3,...]
        # N is number of objects
        # m = [m1, m2, m3,...]
        X_0 = X[:N]
        V_0 = X[N:]
        
        #Hermite scheme (Dehnen and Read 2011)
      
        #Predict positions and velocities
        A_0 = acc(N, X_0, m)
        J_0 = jerk(N, X_0, V_0, m)
        X_p = X_0 + V_0*dt + A_0*dt**2.0/2.0 + J_0*dt**3.0/6.0
        V_p = V_0 + A_0*dt + J_0*dt**2.0/2.0            
        
        #Estimate acceleration and jerk
        A_1 = acc(N, X_p, m)
        J_1 = jerk(N, X_p, V_p, m)
        
        #Obtain corrected positions and velocities
        X_1 = X_0 + (V_p+V_0)*dt/2.0 + (A_0-A_1)*dt**2.0/12.0
        V_1 = V_0 + (A_1+A_0)*dt/2.0 + (J_0-J_1)*dt**2.0/12.0
        
        return np.concatenate([X_1, V_1])
        
        
def acc(int N, np.ndarray[double, ndim=2] X, np.ndarray[double, ndim=1] m):
        
        cdef np.ndarray A = np.zeros((N,3), dtype=float)
        for i in range(N):
                for j in it.chain(range(i), range(i+1, N)):
                        x_ij = X[i] - X[j]
                        A[i] += m[j]*x_ij/(np.linalg.norm(x_ij))**3.0
        return -G*A

def jerk(int N, np.ndarray[double, ndim=2] X, np.ndarray[double, ndim=2] V, np.ndarray[double, ndim=1] m):
        
        cdef np.ndarray J = np.zeros((N,3), dtype=float)
        for i in range(N):
                for j in it.chain(range(i), range(i+1, N)):
                        x_ij = X[i] - X[j]
                        v_ij = V[i] - V[j]
                        J[i] += m[j]*(np.dot(x_ij,x_ij)*v_ij - 3.0*x_ij*np.dot(x_ij,v_ij))/(np.linalg.norm(x_ij))**5.0
        return -G*J
        
        
        
        
        
        
        
        
        
        
        
        
        
        
                