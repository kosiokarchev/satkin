import numpy as np
from astropy import units as u
from astropy.cosmology import Planck13
from astropy.constants import G
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

def mvir2rvir(mvir, z):
    rho = Planck13.clone(H0=100).critical_density(z).to('solMass/Mpc3').value
    return (np.array(mvir, copy=False) / ((4*np.pi/3)*200*rho))**(1/3)

def mvir2vvir(mvir, z):
    rvir = mvir2rvir(mvir, z) * u.Mpc
    return np.sqrt(G * mvir * u.solMass / rvir).to('km/s').value

def mv2s2(mv):
    """
    Lokas (2001) with c <- Dutton (2014) at z=1 and beta=const=0.25

    :param mv: in logarithmic units
    :return: sigma^2 in (km/s)^2
    """
    return np.power(10, 0.66172139*mv - 3.74132182)