from scipy.constants import G as G_SI
from scipy.constants import giga, year
 
 

def time_scale():
        return giga*year

def mass_scale():
        return 2.0*10.0**30.0

def G():
        return 1.0

def length_scale():
        m_scale = mass_scale()
        t_scale = time_scale()
        return (G_SI * m_scale * t_scale**2.0)**(1.0/3.0)