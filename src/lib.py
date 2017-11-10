from __future__ import print_function, division
import numpy as np
from astropy.table import Table, join as join
from matplotlib import pyplot as plt


SNDATA = {
    34: {'z': 1.5, 'clr': 'red'},
    38: {'z': 1.0, 'clr': 'green'},
    45: {'z': 0.5, 'clr': 'blue'}
}

BINWIDTH = 0.1
PLOTDATA = {
    'bins': {
        'mvir': np.arange(10, 15 + BINWIDTH, BINWIDTH),
        'stellarMass': np.arange(7, 12 + BINWIDTH, BINWIDTH)
    }
}

FILES = {
    'cube': lambda sn, sats=False: 'data/cube{}{}0.fits'.format(sn, '>' if sats else '='),
    'masses': lambda sn: 'data/masses{}.fits'.format(sn),
    '3dbins': lambda sn: 'data/3dbins{}.pickle'.format(sn),
    'sats': lambda sn: 'data/sats{}.fits'.format(sn),
    'sigmas': lambda sn, mname='stellarMass', observe=False: 'data/{}sigmas-{}{}.csv'.format('observe/' if observe else '', mname, sn),
    'plotsigmas': lambda sn, mname='stellarMass', observe=False: 'plots/{}sigmas-{}{}.pdf'.format('observe/' if observe else '', mname, sn)
}


def deal(t):
    t = t[np.where(t['stellarMass'] > 1e6)]
    t['mv'] = 10 + np.log10(t['mvir'])
    t['ms'] = 10 + np.log10(t['stellarMass'])
    return t


def percentile_wrapper(q):
    return lambda a: np.percentile(a, q)


from binner import bin2d
from observer import Observer, Regressor
from data import Data
