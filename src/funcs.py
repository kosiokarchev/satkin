from __future__ import print_function, division
import json
import numpy as np
from astropy.stats import sigma_clip
from astropy.table import Table, join, unique
import lmfit
from lmfit import models
from matplotlib import pyplot as plt

from print_function import print
from settings import *


loaded_tables = {}
def load(fname, **kwargs):
    if fname not in loaded_tables:
        print('Loading', fname)
        t = Table.read(fname, **kwargs)
        loaded_tables[fname] = t
    return loaded_tables[fname]
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

def deal(t):
    t = t[np.where(t['stellarMass'] > 1e-4)]
    t['mv'] = 10 + np.log10(t['mvir'])
    t['ms'] = 10 + np.log10(t['stellarMass'])
    return t


def percentile_wrapper(q):
    return lambda a: np.percentile(a, q)


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


# from binner import bin2d, binX
# from observer import Observer, Regressor
# from data import Data
# from proclib import *
