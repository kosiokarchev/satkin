from __future__ import division, print_function

import numpy as np
from astropy.table import Table
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from mysql_getter import mysql_load

def bin2D(t, x, y, nbins=20, bins=None, funcsx=None, funcsy=None):
    if funcsx is None:
        funcsx = (np.mean, np.std)
    elif not isinstance(funcsx, tuple):
        funcsx = tuple(funcsx)

    if funcsy is None:
        funcsy = (np.mean, np.std)
    elif funcsy is True:
        funcsy = funcsx
    elif not isinstance(funcsy, tuple):
        funcsy = tuple(funcsy)


    if bins is not None:
        width = bins[1] - bins[0]
        left = bins[0]
        nbins = len(bins) - 1
    else:
        left = np.min(t[x])
        width = (np.max(t[x]) - left) / nbins

    t['bin'] = np.floor_divide(t[x] - left, width)
    g = t.group_by('bin')

    xstats = []
    ystats = []
    for key, group in zip(g.groups.keys, g.groups):
        if 0 <= key['bin'] < nbins:
            xstats.append([f(group[x]) for f in funcsx])
            ystats.append([f(group[y]) for f in funcsy])

    return np.array(xstats).transpose(), np.array(ystats).transpose()


def binner(t, x, y, nbinsx=20, nbinsy=20, binsx=None, binsy=None, funcsx=None, funcsy=None):
    binx = bin2D(t, x, y, nbinsx, binsx, funcsx, funcsy)
    biny = bin2D(t, y, x, nbinsy, binsy, funcsx, funcsy)
    biny = biny[1], biny[0]

    return binx, biny


def percentile_wrapper(q):
    return lambda a: np.percentile(a, q)


def plot_both_binnings(bin1, bin2, label1, label2):
    plt.fill_between(bin1[0][0], bin1[1][3], bin1[1][4],
                     color='blue', alpha=0.5, label=label1)
    plt.plot(bin1[0][0], bin1[1][2], 'o', color='blue', markeredgecolor='red')
    plt.errorbar(bin1[0][0], bin1[1][0], bin1[1][1], bin1[0][1], color='black')

    plt.fill_betweenx(bin2[1][0], bin2[0][3], bin2[0][4],
                      color='yellow', alpha=0.5, label=label2)
    plt.plot(bin2[0][2], bin2[1][0], 'o', color='yellow', markeredgecolor='red')
    plt.errorbar(bin2[0][0], bin2[1][0], bin2[1][1], bin2[0][1], color='black')


def mstar_mhalo(sparse=333, ):
    t = Table.read('data/central-mvir-mste.fits')
    t = t[np.where(t['redshift'] < 1)]


    t['mv'] = 10+np.log10(t['mvir'])
    t['ms'] = 10+np.log10(t['stellarMass'])

    xbounds = np.linspace(11, 14, 16)
    ybounds = np.linspace(9, 11.6, 14)

    bhalo, bstar = binner(t, 'mv', 'ms', binsx=xbounds, binsy=ybounds,
                          funcsy=(np.mean, np.std, np.median, percentile_wrapper(16), percentile_wrapper(84)))

    x, y = t['mv'].copy(), t['ms'].copy()
    x.sort(), y.sort()

    fig = plt.figure()
    ax = fig.gca()

    r = (np.random.random(int(len(t) / sparse))*len(t)).astype(int)

    plt.plot(t['mv'][r], t['ms'][r], 'r,', zorder=-1, label='_nolegend_')
    plot_both_binnings(bhalo, bstar, '$M_{halo}$ binning', '$M_{stars}$ binning')
    plt.plot(x, y, 'r-', label='abundance matching')

    plt.legend(loc='upper left')

    ax.set_xlabel('$\log_{10}(M_{halo} / M_{\odot})$')
    ax.set_ylabel('$\log_{10}(M_{stars} / M_{\odot})$')

    # ax.set_xbound(xbounds[0], xbounds[-1])
    # ax.set_ybound(ybounds[0], ybounds[-1])
    ax.set_aspect('equal', adjustable='box')
    fig.tight_layout()


def mhalo_sigmav():
    t = Table.read('data/central-stats.csv')

    t['mv'] = 10 + np.log10(t['mvir'])
    t['sv'] = np.log10(t['sigma_vpec'])

    mbounds = np.linspace(12, 15, 16)
    sbounds = np.linspace(0, 3, 16)

    bsigma, bmass = binner(t, 'sv', 'mv', binsx=sbounds, binsy=mbounds,
                           funcsy=(
                           np.mean, np.std, np.median, percentile_wrapper(16),
                           percentile_wrapper(84)))

    fig = plt.figure()
    ax = fig.gca()

    plt.plot(t['sv'][::3], t['mv'][::3], 'r,', zorder=-1, label='_nolabel_')
    plot_both_binnings(bsigma, bmass, '$\sigma_v$ binning',
                       '$M_{halo}$ binning')

    plt.legend(loc='upper left')

    ax.set_xlabel('$\log_{10}(\sigma_v\ [km/s])$')
    ax.set_ylabel('$\log_{10}(M_{halo} / M_{\odot})$')

    ax.set_xbound(sbounds[0], sbounds[-1])
    ax.set_ybound(mbounds[0], mbounds[-1])
    ax.set_aspect('equal', adjustable='box')
    fig.tight_layout()


if __name__=='__main__':
    pass