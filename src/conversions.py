import numpy as np
from astropy.cosmology import Planck13
from funcs import two_lines

def s2h(ms, mv0=11.88, ms0=10.13, a1=0.556, a2=1.5):
    """
    Return an estimate of the virial mass, given stellar mass. A fit to the
    abundance matching is used.

    :param ms: linear, 10^10 solMass unit
    """
    return np.power(10, two_lines(10+np.log10(ms), ms0, mv0, a1, a2))

def mvir2c(mvir, z=0):
    """
    Return the concentration parameter based on the virial mass of a halo, and z
    using the fits from Gutton (2014) -- all haloes

    :param mvir: linear, 1 solMass unit
    """
    a = 0.520 + (0.905-0.520)*np.exp(-0.617*z**1.21)
    b = -0.101 + 0.026*z

    return np.power(10, a+b*np.log10(mvir/1e12))