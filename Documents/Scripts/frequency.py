import numpy as np

def calcFrequency(x, N_bins, normalise=False, log=False):
        x_min = np.min(x)
        x_max = np.max([np.max(x), x_min+1.0])
        
        #x bins
        if log:
                dlogx = (np.log(x_max)-np.log(x_min))/(N_bins-1)
                x_bins = np.array([x_min * np.exp(i*dlogx) for i in range(N_bins)])
        else:
                dx = (x_max - x_min)/(N_bins-1)
                x_bins = np.array([x_min + i*dx for i in range(N_bins)])
        
        #Count frequency
        N = np.zeros(N_bins, dtype=int)
        if log:
                for val in x:
                        i = int(np.floor(np.log(val/x_min)/dlogx))
                        N[i] += 1
        else:
                for val in x:
                        i = int(np.floor((val - x_min)/dx))
                        N[i] += 1
        
        #Normalise
        if normalise:
                N /= np.size(x)
                
        return(x_bins, N)
        
        