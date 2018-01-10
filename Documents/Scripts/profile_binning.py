import os
os.system("python setup.py build_ext --inplace")
#os.system("python setup_profiling.py build_ext --inplace")

import pstats, cProfile
#import line_profiler_master.line_profiler
import numpy as np
from encounters import binning
from monte_carlo import MCEncounters

from scipy.constants import G

#Initialise variables
#Semi-major axis, m
a = 0.1 * 3.086*10.0**16.0
#Eccentricity
e = 0.7
#Mass of binary stars
m1 = 2.0*10.0**30.0
m2 = 2.0*10.0**30.0
#Density of dark matter halo solar masses/pc**3
rho = 0.008
#Convert to SI
rho = rho * 2.0*10.0**30.0/((3.086*10.0**16.0)**3.0)
#Mass of perturbers
M_p = 3.0 * 2.0*10.0**30.0
#RMS of Maxwellian velocity distribution, m/s
v_rms = 100.0 * 1000.0
#Time to run simulation for
t_end = 10.0**10.0*365.25*24.0*60.0*60.0
#Number density of perturbers
n_p = rho/M_p

cProfile.runctx("binning(v_rms, n_p, t_end, a, e, m1, m2, M_p)", globals(), locals(), "binning.prof")

#s = pstats.Stats("binning.prof")
#s.strip_dirs().sort_stats("time").print_stats()

cProfile.runctx("MCEncounters(v_rms, n_p, t_end, m1, m2, M_p, np.array([a]), np.array([e]), 1)", globals(), locals(), "MCEncounters.prof")



'''
#Line profiler
#From https://github.com/cython/cython/blob/master/tests/run/line_profile_test.srctree

def assert_stats(profile, name):
        profile.print_stats()
        stats = profile.get_stats()
        assert len(stats.timing) > 0, "No profile stats."
        for key, timings in stats.timings.items():
                if key[-1] == name:
                        assert len(timings) > 0
                        break
        else:
                raise ValueError("No stats for %s." % name)

func = binning
profile = line_profiler.LineProfiler(func)
profile.runcall(binning, v_rms, n_p, t_end, a, e, m1, m2, M_p)
assert_stats(profile, func.__name__)

func = MCEncounters
profile = line_profiler.LineProfiler(func)
profile.runcall(MCEncounters, v_rms, n_p, t_end, m1, m2, M_p, np.array([a]), np.array([e]), 1)
assert_stats(profile, func.__name__)
'''









