import numpy as np

def calcFrequency(x, N_bins, normalise=False):
        x_min = np.min(x)
        x_max = np.max(x)
        
        #x bins
        dx = (x_max - x_min)/(N_bins-1)
        x_bins = np.array([x_min + i*dx for i in range(N_bins)])
        
        #Count frequency
        N = np.zeros(N_bins, dtype=int)
        for val in x:
                i = int(np.floor((val - x_min)/dx))
                N[i] += 1
        
        #Normalise
        if normalise:
                N /= np.size(x)
                
        return(x_bins, N)
        
        