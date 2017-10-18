import math
import numpy as np

G = 6.67 * 10.0**(-11.0)

def evolveBinary(x1, x2, v1, v2, m1, m2, dt):
        # dt is step size
        
        #Hermite scheme (Dehnen and Read 2011)
        
        a1 = acc1(x1, x2, m1, m2)
        a2 = acc2(x1, x2, m1, m2)
        j1 = jerk1(x1, x2, v1, v2, m1, m2)
        j2 = jerk2(x1, x2, v1, v2, m1, m2)
        
        #Predict positions and velocities
        x1_p = x1 + v1*dt + 0.5*a1*dt**2.0 + (j1*dt**3.0)/6.0
        v1_p = v1 + a1*dt + 0.5*j1*dt**2.0
        x2_p = x2 + v2*dt + 0.5*a2*dt**2.0 + (j2*dt**3.0)/6.0
        v2_p = v2 + a2*dt + 0.5*j2*dt**2.0
        
        #Estimate acceleration and jerk
        a1_1 = acc1(x1_p, x2_p, m1, m2)
        a2_1 = acc2(x1_p, x2_p, m1, m2)
        j1_1 = jerk1(x1_p, x2_p, v1_p, v2_p, m1, m2)
        j2_1 = jerk2(x1_p, x2_p, v1_p, v2_p, m1, m2)
        
        #Obtain corrected positions and velocities
        x1_1 = x1 + (v1_p+v1)*dt/2.0 +((a1-a1_1)*dt**2.0)/12.0
        x2_1 = x2 + (v2_p+v2)*dt/2.0 +((a2-a2_1)*dt**2.0)/12.0
        v1_1 = v1 + (a1_1+a1)*dt/2.0 +((j1-j1_1)*dt**2.0)/12.0
        v2_1 = v2 + (a2_1+a2)*dt/2.0 +((j2-j2_1)*dt**2.0)/12.0        
        
        return(x1_1, x2_1, v1_1, v2_1)
        
             
        
        
def acc1(x1, x2, m1, m2):
        r = np.linalg.norm(x2-x1)
        a1 = G*m2/(r**3.0) * (x2 - x1)
        return a1

def acc2(x1, x2, m1, m2):
        r = np.linalg.norm(x2-x1)
        a2 = - G*m1/(r**3.0) * (x2 - x1)
        return a2


def jerk1(x1, x2, v1, v2, m1, m2):
        r = np.linalg.norm(x2-x1)
        x12 = x1-x2
        v12 = v1-v2
        j1 = -G*m2*(np.dot(x12,x12)*v12-3.0*x12*np.dot(x12,v12))/r**5.0
        return j1

def jerk2(x1, x2, v1, v2, m1, m2):
        r = np.linalg.norm(x2-x1)
        x21 = x2-x1
        v21 = v2-v1
        j2 = -G*m1*(np.dot(x21,x21)*v21-3.0*x21*np.dot(x21,v21))/r**5.0
        return j2