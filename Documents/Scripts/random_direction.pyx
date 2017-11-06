import random
import numpy as np

#To find a random 3D direction

def randomDirection():

        u = random.uniform(0.0, 1.0)
        theta = random.uniform(0.0, 2.0*np.pi)
        vec = np.array([np.sqrt(1.0 - u**2.0)*np.cos(theta), np.sqrt(1.0 - u**2.0)*np.sin(theta), u])
        
        return vec/np.linalg.norm(vec)