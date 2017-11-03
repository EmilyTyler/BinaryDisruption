#Solve Kepler's equation for eccentric anomaly
import numpy as np

def findEccentricAnomaly(e, M):

        #Solves Kepler's equation E-esinE=M to find eccentric anomaly E given eccentricity e and mean anomaly M.
        #From page 36 of Solar System Dynamics, or Danby 1988
        E = M + np.sign(np.sin(M))*0.85*e
        count = 0
        while abs(E - e*np.sin(E) - M) > 10.0**(-8.0):
                f = E - e*np.sin(E) - M
                f_p = 1.0 - e*np.cos(E)
                f_pp = e*np.sin(E)
                f_ppp = e*np.cos(E)

                d_1 = - f/f_p
                d_2 = -f/(f_p + 0.5*d_1*f_pp)
                d_3 = -f/(f_p + 0.5*d_2*f_pp + d_2**2.0*f_ppp/6.0)
        
                E += d_3
        
                if count > 100:
                        print ('findEccentricAnomaly did not converge')
                        break
         
        return E 