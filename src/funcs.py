from __future__ import print_function, division
import json
import numpy as np
from scipy.io import readsav
from scipy.special import erf
from astropy.stats import sigma_clip
from astropy.table import Table, join, unique
from astropy.cosmology import Planck13
import lmfit
from lmfit import models
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from print_function import print
from settings import *


loaded_tables = {}
def load(fname, **kwargs):
    if fname not in loaded_tables:
        print('Loading', fname)
        t = Table.read(fname, **kwargs)
        loaded_tables[fname] = t
    return loaded_tables[fname]
def unload(fname):
    if fname in loaded_tables:
        del loaded_tables[fname]
def write(t, fname, **kwargs):
    print('Writing to', fname)
    t.write(fname, overwrite=True, **kwargs)

    loaded_tables[fname] = t


def get_sparse(length, sparse=333):
    return np.random.randint(0, length, int(length / sparse))

def load_table(fname, **kwargs):
    print('Loading', fname)
    return Table.read(fname, **kwargs)
def write_table(t, fname, **kwargs):
    print('Writing to', fname)
    t.write(fname, overwrite=True, **kwargs)

def read_idl(fname, prefix=None):
    s = readsav(fname)
    names = s.keys()
    return Table([s[key] for key in names],
                 names=[k.replace(prefix+'_', '') for k in names] if prefix is not None else names)

def fixsw(fname):
    from astropy import units as u
    from astropy.cosmology import Planck13
    cosmo = Planck13.clone(H0=100*u.km/u.s/u.Mpc)

    t = load(fname)

    h = 0.7

    t['galaxyId'] = np.arange(len(t), dtype=np.int64)
    t = t['galaxyId', 'ra', 'dec', 'zbest', 'sed_lmass'][t['zbest_type'] == 1]
    t.rename_column('zbest', 'z_app')
    t.rename_column('sed_lmass', 'stellarMass')
    t = t[t['z_app'] > 0]
    t['stellarMass'] = np.power(10, t['stellarMass']) / 1e10 * h
    t['ra'] = np.deg2rad(t['ra'])
    t['dec'] = np.deg2rad(t['dec'])
    t['d_comoving'] = cosmo.comoving_distance(t['z_app'])

    return t

def deal(t, minms=6):
    t = t[np.where(t['stellarMass'] > np.power(10., minms-10))]
    t['mv'] = 10 + np.log10(t['mvir'])
    t['ms'] = 10 + np.log10(t['stellarMass'])
    return t


def percentile_wrapper(q):
    return lambda a: np.percentile(a, q)
def percentile(x, y, p, d=None):
    if d is None:
        d = np.full_like(p, np.nan)
    return np.array([x[np.searchsorted(np.cumsum(a), p)]
                    if not np.all(np.isnan(a)) else d
                    for a in y]).T


def fit_or(model, data, outlier_func=sigma_clip, niter=3, outlier_kwargs=None, **model_kwargs):
    if outlier_kwargs is None:
        outlier_kwargs = {}
        if outlier_func is sigma_clip:
            outlier_kwargs['iters'] = 1
            outlier_kwargs['sigma'] = 3

    fit_res = model.fit(data, **model_kwargs)
    indices = np.arange(len(data))
    for i in range(niter):
        filtered_data = outlier_func(fit_res.residual, **outlier_kwargs)

        mask = ~filtered_data.mask
        indices = indices[mask]
        if mask.all():
            break

        data = data[mask]
        for iv in model.independent_vars:
            niv = model_kwargs[iv][mask]
            del model_kwargs[iv]
            model_kwargs[iv] = niv

        fit_res = model.fit(data, **model_kwargs)

    return indices, fit_res

def powerlaw(x, amplitude, exponent, scale=1):
    return amplitude*np.power(x/scale, exponent)
def logpowerlaw(x, amplitude, exponent, scale=1):
    return np.log10(amplitude) + exponent * (x - np.log10(scale))
def line(x, intercept, slope, offset=0):
    return intercept + slope * (x-offset)

def two_lines_1(x, x0, y0, a1, a2):
    res = np.full_like(x, y0)
    lower = x<x0
    res[lower] += a1*(x[lower]-x0)
    res[~lower] += a2*(x[~lower]-x0)

    return res


def two_lines_2(x, x0, y0, a1, a2, scale=1):
    da = a2-a1
    dx = (x-x0) / scale
    logpart = np.where(dx < 300,
                       np.log10((1 + np.power(10, dx)) / 2),
                       dx - np.log10(2))
    return y0 + a1*dx + da*logpart


two_lines = two_lines_2


def ecurve(x, a=1, x0=1, b=0):
    return np.power(10, a*(x-x0)) + b


def cumgauss(x, sigma, A, B):
    return A * erf(x / (np.sqrt(2) * sigma)) + B * x

def schechter_func(x, a, x0, b, c):
    dx = x - x0
    return c + (a+1) * dx - np.power(10, b*dx) / np.log(10)
class SchechterModel(lmfit.model.Model):
    def __init__(self, **kws):
        super(SchechterModel, self).__init__(schechter_func, **kws)
        self.name = 'Schechter'

from conversions import s2h, mvir2c, mvir2rvir, mvir2vvir, mv2s2

# from binner import bin2d, binX
# from observer import Observer, Regressor
# from data import Data
# from proclib import *
