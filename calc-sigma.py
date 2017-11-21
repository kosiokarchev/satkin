from __future__ import print_function
import numpy as np
from astropy.table import Table, join
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def bin(n, mname, left, width):
    print('Loading table...')
    sats = Table.read('data/sats' + str(n) + '.fits', format='fits')
    sats = sats[np.where(sats['stellarMass'] > 10 ** -6)]

    sats['masslog'] = 10 + np.log10(sats[mname])

    sats['bin'] = np.floor_divide(sats['masslog'] - left, width)
    satsbinned = sats.group_by('bin')

    stdev = satsbinned[
        'bin', 'vpec', 'vpecx', 'vpecy', 'vpecz'].groups.aggregate(np.std)

    stdev.rename_column('vpec', 'sigma')
    stdev.rename_column('vpecx', 'sigmax')
    stdev.rename_column('vpecy', 'sigmay')
    stdev.rename_column('vpecz', 'sigmaz')

    stdev['mean'] = satsbinned['vpec'].groups.aggregate(np.mean)
    stdev['median'] = satsbinned['vpec'].groups.aggregate(np.median)
    stdev['mass'] = left + (stdev['bin'] + 0.5) * width

    return sats, stdev

def sigma(n, mname, xlim, xlabel='$\log_{10}(M / M_\odot)$', width=0.1, sparse=333, write=False):
    sats, stdev = bin(n, mname, xlim[0], width)
    if write:
        stdev.write('data/mstar-stdev'+str(n)+'.fits', format = 'fits')

    print('Plotting...')

    fig, ax = plt.subplots(2, 1, sharex='col', squeeze=True,
                           gridspec_kw={'height_ratios': [3, 1], 'hspace': 0})

    plt.sca(ax[0])
    r = (np.random.random(int(len(sats) / sparse))*len(sats)).astype(int)
    plt.plot(sats['masslog'][r], sats['vpec'][r], ',', label='_nolegend_')
    plt.plot(stdev['mass'], stdev['sigma'], '-r', linewidth = 3, label ='$\sigma_{3D}$')
    plt.plot(stdev['mass'], stdev['sigmax'], '-y', linewidth = 2, label = '$\sigma_{los}$')
    plt.plot(stdev['mass'], stdev['sigmay'], '-y', linewidth = 2, label = '_nolegend_')
    plt.plot(stdev['mass'], stdev['sigmaz'], '-y', linewidth = 2, label = '_nolegend_')
    plt.plot(stdev['mass'], stdev['median'], '-k', linewidth = 3, label = '_nolegend_')
    plt.plot(stdev['mass'], stdev['mean'], '-g', linewidth = 3, label = '$\langle |v| \\rangle$')

    plt.legend(loc = 'lower right')

    plt.ylabel('Speed, km/s')

    plt.semilogy()

    print('Plotting...')
    plt.sca(ax[1])
    plt.plot(stdev['mass'], stdev['mean'] / stdev['sigmax'])
    plt.plot(xlim, [np.sqrt(3), np.sqrt(3)], ':')
    plt.plot(xlim, [np.pi/2, np.pi/2], ':')

    plt.xlabel(xlabel)
    plt.ylabel('$v_{3D}$ / $\sigma_{los}$')

    plt.ylim(1, 2)
    plt.xlim(*xlim)

def plot_sigma(sats, stdev, xlim, xlabel='$\log_{10}(M / M_\odot)$', sparse=333):
    print('Plotting...')
    fig, ax = plt.subplots(2, 1, sharex='col', squeeze=True,
                           gridspec_kw={'height_ratios': [3, 1], 'hspace': 0})

    plt.sca(ax[0])
    r = (np.random.random(int(len(sats) / sparse)) * len(sats)).astype(int)
    plt.plot(sats['masslog'][r], sats['vpec'][r], ',', label='_nolegend_')
    plt.plot(stdev['mass'], stdev['sigma'], '-r', linewidth=3, label='$\sigma_{3D}$')
    plt.plot(stdev['mass'], stdev['sigmax'], '-y', linewidth=2, label='$\sigma_{los}$')
    plt.plot(stdev['mass'], stdev['sigmay'], '-y', linewidth=2, label='_nolegend_')
    plt.plot(stdev['mass'], stdev['sigmaz'], '-y', linewidth=2, label='_nolegend_')
    plt.plot(stdev['mass'], stdev['median'], '-k', linewidth=3, label='_nolegend_')
    plt.plot(stdev['mass'], stdev['mean'], '-g', linewidth=3, label='$\langle |v| \\rangle$')

    plt.legend(loc='lower right')

    plt.ylabel('Speed, km/s')

    plt.semilogy()

    print('Plotting...')
    plt.sca(ax[1])
    plt.plot(stdev['mass'], stdev['mean'] / stdev['sigmax'])
    plt.plot(xlim, [np.sqrt(3), np.sqrt(3)], ':')
    plt.plot(xlim, [np.pi / 2, np.pi / 2], ':')

    plt.xlabel(xlabel)
    plt.ylabel('$v_{3D}$ / $\sigma_{los}$')

    plt.ylim(1, 2)
    plt.xlim(*xlim)

def dobin(t, left, width=0.1):
    t['bin'] = np.floor_divide(t['masslog'] - left, width)
    binned = t.group_by('bin')

    res = binned['bin', 'vpec', 'vpecx', 'vpecy', 'vpecz'].groups.aggregate(np.std)

    res.rename_column('vpec', 'sigma')
    res.rename_column('vpecx', 'sigmax')
    res.rename_column('vpecy', 'sigmay')
    res.rename_column('vpecz', 'sigmaz')

    res['mean'] = binned['vpec'].groups.aggregate(np.mean)
    res['median'] = binned['vpec'].groups.aggregate(np.median)
    res['mass'] = left + (res['bin'] + 0.5) * width

    return res

def doplots(n, write=False):
    print('Loading table...')
    sats = Table.read('data/sats' + str(n) + '.fits', format='fits')
    sats = sats[np.where(sats['stellarMass'] > 1e-6)]

    for mname, xlim in zip(('mvir', 'stellarMass'),
                           ((9, 15), (7, 12))):
        sats['masslog'] = 10 + np.log10(sats[mname])
        stdev = dobin(sats, xlim[0])

        if write:
            stdev.write('data/sigmas-'+mname+str(n)+'.fits', overwrite=True)

        plot_sigma(sats, stdev, xlim)
        plt.savefig('plots/sigmas-'+mname+str(n)+'.pdf')
        print('Figure ' + mname + str(n) + ' saved.')
        plt.close()

def fit(func,x_data,y_data):
    popt, pcov = curve_fit(func, x_data, y_data)
    return popt

for n in [34,38,45]:
    doplots(n)
    # sigma(n, 'mvir', (9, 15))
    # plt.savefig('plots/sigmas-mstar'+str(n)+'.pdf')
    # print('Figure '+str(n)+ ' saved.')
    # plt.close()
    fit(linear, stdev['mass'], np.log10(stdev['sigmax']))