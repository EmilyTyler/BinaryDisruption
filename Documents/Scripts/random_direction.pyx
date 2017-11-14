import random
import numpy as np
cimport numpy as np

#To find a random 3D direction



def randomDirection():
        
        cdef double u, theta
        cdef np.ndarray vec = np.zeros(3, dtype=float)
        
        u = random.uniform(-1.0, 1.0)
        theta = random.uniform(0.0, 2.0*np.pi)
        vec = np.array([np.sqrt(1.0 - u**2.0)*np.cos(theta), np.sqrt(1.0 - u**2.0)*np.sin(theta), u])
        
        return vec/np.linalg.norm(vec)