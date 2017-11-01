import math
import numpy as np
from orbital_elements import orbitalElements
from eccentric_anomaly import findEccentricAnomaly

G = 6.67 * 10.0**(-11.0)

#Evolution of binary through analytics
def evolveBinary(X, m1, m2, t):
        #Find orbital elements
        (a, e, f, I, Omega, omega) = orbitalElements(X, m1, m2)
        #Find eccentric anomaly
        E = np.arccos((1.0 - (np.linalg.norm(X[0]-X[1]))/a)/e)
        #Mean motion
        n = np.sqrt(G*(m1+m2)/a**3.0)
        #Find eccentric anomaly after time t
        E_new = findEccentricAnomaly(e, (E-e*np.sin(E)+n*t))
        #Find new true anomaly
        f_new = 2.0*np.arctan(((1.0+e)/(1.0-e))**0.5 * np.tan(0.5*E_new))
        #Find x and y positions in orbital plane
        r_new = a*(1.0 - e**2.0)/(1.0 + e*np.cos(f_new))
        #New coordinates of first star (cartesian)
        X[0] = np.array([0.0, 0.0, 0.0])
        #New velocity of first star
        X[2] = np.array([0.0, 0.0, 0.0])
        #New coordinates of second star (cartesian)
        X[1] = np.array([r_new*math.cos(f_new), r_new*math.sin(f_new), 0.0])
        #New velocity of second star
        X[3] = np.array([- n*a/(math.sqrt(1.0-e**2.0)) * math.sin(f_new), n*a/(math.sqrt(1.0-e**2.0)) * (e + math.cos(f_new)), 0.0])
        
        #Rotate into reference frame
        R = np.array([[np.cos(Omega)*np.cos(omega) - np.sin(Omega)*np.cos(I)*np.sin(omega), -np.cos(Omega)*np.sin(omega)-np.sin(Omega)*np.cos(I)*np.cos(omega), np.sin(I)*np.sin(Omega)],
             [np.sin(Omega)*np.cos(omega) + np.cos(Omega)*np.cos(I)*np.sin(omega), -np.sin(Omega)*np.sin(omega)+np.cos(Omega)*np.cos(I)*np.cos(omega),-np.sin(I)*np.cos(Omega)],
             [np.sin(I)*np.sin(omega), np.sin(I)*np.cos(omega), np.cos(I)]])
        X[1] = np.transpose(np.dot(R, np.transpose(X[1])))
        X[3] = np.transpose(np.dot(R, np.transpose(X[3])))
        
        return X


#Evolution of binary through integration
def integrateBinary(x1, x2, v1, v2, m1, m2, dt):
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