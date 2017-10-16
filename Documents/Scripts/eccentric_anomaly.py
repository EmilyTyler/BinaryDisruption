#Solve Kepler's equation for eccentric anomaly
import math

def findEccentricAnomaly(e, M):

        #Solves Kepler's equation E-esinE=M to find eccentric anomaly E given eccentricity e and mean anomaly M.
        #From: http://alpheratz.net/dynamics/twobody/KeplerIterations_summary.pdf
        tol = 10.0**(-14.0)
        Mnorm = M % (2.0*math.pi)
        E0 = startingValue(e, Mnorm)
        dE = tol + 1.0     
        count = 0
        while (dE > tol):
                E = E0 - eps3(e,Mnorm,E0)
                dE = abs(E-E0)
                E0 = E
                count = count + 1
                if (count==100):
                        print('findEccentricAnomaly failed to converge!')
                        break
        print 'Mean anomaly = ', M
        print 'Mean anomaly calculated from E = ', (E-e*math.sin(E))
        return E
        
        
        
def startingValue(e, M):
        t34 = e**2.0
        t35 = e*t34
        t33 = math.cos(M)
        return M+(-1.0/2.0*t35+e+(t34+3.0/2.0*t33*t35)*t33)*math.sin(M)


def eps3(e,M,x):
        t1 = math.cos(x)
        t2 = -1.0+e*t1
        t3 = math.sin(x)
        t4 = e*t3
        t5 = -x+t4+M
        t6 = t5/(1.0/2.0*t5*t4/t2+t2)
        return t5/((1.0/2.0*t3 - 1.0/6.0*t1*t6)*e*t6+t2)


        
        
        