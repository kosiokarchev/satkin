from __future__ import print_function, division
import json
import numpy as np
from astropy.stats import sigma_clip
from astropy.table import Table, join as join
from lmfit import models
from matplotlib import pyplot as plt


SNDATA = {
    34: {'z': 1.5, 'clr': 'red'},
    38: {'z': 1.0, 'clr': 'green'},
    45: {'z': 0.5, 'clr': 'blue'}
}

BINWIDTH = 0.1
PLOTDATA = {
    'bins': {
        'mvir': np.arange(10, 14 + BINWIDTH, BINWIDTH),
        'stellarMass': np.arange(7, 12 + BINWIDTH, BINWIDTH)
    }
}

FILES = {
    'cube': lambda sn, sats=False: 'data/cube{}{}0.fits'.format(sn, '>' if sats else '='),
    'masses': lambda sn: 'data/masses{}.fits'.format(sn),
    'nums': lambda sn: 'data/nums{}.fits'.format(sn),
    '3dbins': lambda sn: 'data/3dbins{}.json'.format(sn),
    'sats': lambda sn: 'data/sats{}.fits'.format(sn),

    'sw-reg-mvir': lambda sn: 'data/sw-reg-mvir{}.csv'.format(sn),
    'sw-sigmas': lambda sn, mname='stellarMass', observe=False: 'data/{}sigmas-{}{}.csv'.format('observe/' if observe else '', mname, sn),
    'plot-sw-sigmas': lambda sn, mname='stellarMass', observe=False: 'plots/{}sigmas-{}{}.pdf'.format('observe/' if observe else '', mname, sn),

    'hw-reg-mvir': lambda sn: 'data/hw-reg-mvir{}.json'.format(sn),
    'hw-rms': lambda sn, observe=False: 'data/{}hw-rms{}.fits'.format('observe/' if observe else '', sn),
    'hwp-sigmas': lambda sn, observe=False: 'data/{}hwp-sigmas{}.csv'.format('observe/' if observe else '', sn)
}


def get_sparse(length, sparse=333):
    return np.random.randint(0, length, int(length / sparse))


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


from binner import bin2d, binX
from observer import Observer, Regressor
from data import Data
from proclib import *
