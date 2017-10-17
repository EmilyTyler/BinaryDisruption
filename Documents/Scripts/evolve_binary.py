import math
import numpy as np

G = 6.67 * 10.0**(-11.0)

def evolveBinary(x1, x2, v1, v2, m1, m2, dt):
        # dt is step size
        
        #Taylor expansion for now
        a1 = acc1(x1, x2, m1, m2)
        a2 = acc2(x1, x2, m1, m2)
        
        v1 = v1 + dt*a1
        v2 = v2 + dt*a2
        
        x1 = x1 + dt*v1 + 0.5*dt**2.0*a1
        x2 = x2 + dt*v2 + 0.5*dt**2.0*a2
        
        return(x1, x2, v1, v2)
        
             
        
        
def acc1(x1, x2, m1, m2):
        r = np.linalg.norm(x2-x1)
        a1 = G*m1*m2/(m1*r**3.0) * (x2 - x1)
        return a1

def acc2(x1, x2, m1, m2):
        r = np.linalg.norm(x2-x1)
        a2 = - G*m1*m2/(m2*r**3.0) * (x2 - x1)
        return a2

        