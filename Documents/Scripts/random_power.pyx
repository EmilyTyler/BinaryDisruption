#Draws numbers from a power law distribution

import numpy as np
cimport numpy as np

def randomPower(a, b, g, size=1):
        #Draws from power law for pdf(x)\propto x^{g-1} for a<=x<=b
        r = np.random.random(size=size)
        ag, bg = a**g, b**g
        return (ag + (bg - ag)*r)**(1.0/g)