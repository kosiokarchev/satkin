import numpy as np
from astropy.table import Table, join as tjoin
from matplotlib import pyplot as plt

def percentile_wrapper(q):
    return lambda a: np.percentile(a, q)


def binX(t, x, y, nbins=20, bins=None, funcsx=None, funcsy=None):
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

    t['bin'] = np.floor_divide(t[x] - left, width).astype(int)
    g = t.group_by('bin')

    xstats = []
    ystats = []
    for i in range(len(g.groups.keys)):
        if 0 <= g.groups.keys[i]['bin'] < nbins:
            xstats.append([f(g[x].groups[i]) for f in funcsx])
            ystats.append([f(g[y].groups[i]) for f in funcsy])

    return np.array(xstats).transpose(), np.array(ystats).transpose()
def bin2D(t, x, y, nbinsx=20, nbinsy=20, binsx=None, binsy=None, funcsx=None, funcsy=None):
    binx = binX(t, x, y, nbinsx, binsx, funcsx, funcsy)
    biny = binX(t, y, x, nbinsy, binsy, funcsx, funcsy)
    biny = biny[1], biny[0]

    return binx, biny

def plot_both_binnings(bin1, bin2, label1=None, label2=None, clr1=None, clr2=None):
    if label1 is None:
        label1 = 'x binning'
    if label2 is None:
        label2 = 'y binning'
    if clr1 is None:
        clr1 = 'blue'
    if clr2 is None:
        clr2 = 'yellow'

    plt.fill_between(bin1[0][0], bin1[1][5], bin1[1][6],
                     color=clr1, alpha=0.25, label=label1)
    plt.fill_between(bin1[0][0], bin1[1][3], bin1[1][4],
                     color=clr1, alpha=0.5, label='_nolegend_')
    plt.plot(bin1[0][0], bin1[1][2], 'o', color=clr1, markeredgecolor='red')
    plt.errorbar(bin1[0][0], bin1[1][0], bin1[1][1], bin1[0][1], color='black')

    plt.fill_betweenx(bin2[1][0], bin2[0][5], bin2[0][6],
                      color=clr2, alpha=0.25, label=label2)
    plt.fill_betweenx(bin2[1][0], bin2[0][3], bin2[0][4],
                      color=clr2, alpha=0.5, label='_nolegend_')
    plt.plot(bin2[0][2], bin2[1][0], 'o', color=clr2, markeredgecolor='red')
    plt.errorbar(bin2[0][0], bin2[1][0], bin2[1][1], bin2[0][1], color='black')