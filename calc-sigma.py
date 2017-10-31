from __future__ import print_function
import numpy as np
from astropy.table import Table, join
import matplotlib.pyplot as plt

def rms(x):
    return np.sqrt(np.mean(x*x))


# sats = Table.read('data/cube38>0.fits', format = 'fits')
# centrals = Table.read('data/cube38=0.fits', format = 'fits')

cube = Table.read('cube.fits', format = 'fits')
sats = cube['velX', 'velY', 'velZ', 'stellarMass', 'centralMVir', 'fofCentralId'][np.where(cube['type']>0)]
centrals = cube['velX', 'velY', 'velZ', 'stellarMass', 'galaxyId'][np.where(cube['type']==0)]
del cube

centrals.rename_column('galaxyId', 'fofCentralId')
sats = join(sats,centrals,'fofCentralId',table_names=['s', 'c']) #column names like name_s/name_c

sats['vpecx'] = sats['velX_s'] - sats['velX_c']
sats['vpecy'] = sats['velY_s'] - sats['velY_c']
sats['vpecz'] = sats['velZ_s'] - sats['velZ_c']
sats['vpec'] = np.sqrt(sats['vpecx']**2+sats['vpecy']**2+sats['vpecz']**2) #magn of peculiar velocity

sats['vpeclog'] = np.log10(sats['vpec'])
sats['centralMVirlog'] = np.log10(sats['centralMVir'])

# plt.figure(1)
# plt.plot(sats['vpeclog'], sats['centralMVirlog'], ',')
# plt.xlabel('log(vpec/(km/s))')
# plt.ylabel('log(centralMVir/(10^10Msun))')

# bins = np.linspace(0,6,51)
left = 1
width = 0.1
sats['bin'] = np.floor_divide(sats['centralMVirlog']-left, width)
satsbinned = sats.group_by('bin')
stdev = satsbinned['bin', 'vpec', 'vpecx', 'vpecy', 'vpecz'].groups.aggregate(np.std)
stdev.rename_column('vpec', 'sigma')
stdev.rename_column('vpecx', 'sigmax')
stdev.rename_column('vpecy', 'sigmay')
stdev.rename_column('vpecz', 'sigmaz')
stdev['mean'] = satsbinned['vpec'].groups.aggregate(np.mean)
stdev['median'] = satsbinned['vpec'].groups.aggregate(np.median)
# stdev['rms'] = satsbinned['vpec'].groups.aggregate(rms)
# stdev['rmsx'] = satsbinned['vpecx'].groups.aggregate(rms)
# stdev['rmsy'] = satsbinned['vpecy'].groups.aggregate(rms)
# stdev['rmsz'] = satsbinned['vpecz'].groups.aggregate(rms)
stdev['mass'] = left + (stdev['bin'] + 0.5) * width

plt.figure(2)
stdev['mass'] += 10
sats['centralMVirlog'] += 10
plt.plot(sats['centralMVirlog'][::10], sats['vpec'][::10], ',')
# plt.plot(sats['centralMVirlog'][::10], sats['vpecx'][::10], ',', alpha = 0.3)
# plt.plot(sats['centralMVirlog'][::10], sats['vpecy'][::10], ',', alpha = 0.3)
# plt.plot(sats['centralMVirlog'][::10], sats['vpecz'][::10], ',', alpha = 0.3)
plt.plot(stdev['mass'], stdev['sigma'], '-r', linewidth = 3, label ='sigma')
plt.plot(stdev['mass'], stdev['sigmax'], '-y', linewidth = 2, label = 'sigmax')
plt.plot(stdev['mass'], stdev['sigmay'], '-y', linewidth = 2, label = 'sigmay')
plt.plot(stdev['mass'], stdev['sigmaz'], '-y', linewidth = 2, label = 'sigmaz')
# plt.plot(stdev['mass'], stdev['rms'], '-y', linewidth = 3, label = 'rms')
# plt.plot(stdev['mass'], stdev['rmsx'], '-y', linewidth = 2, label = 'rmsx')
# plt.plot(stdev['mass'], stdev['rmsy'], '-y', linewidth = 2, label = 'rmsy')
# plt.plot(stdev['mass'], stdev['rmsz'], '-y', linewidth = 2, label = 'rmsz')
plt.plot(stdev['mass'], stdev['mean'], '-g', linewidth = 3, label = 'mean')
plt.plot(stdev['mass'], stdev['median'], '-k', linewidth = 3, label = 'median')
plt.xlabel('log($M_{halo}$ / $M_\odot$)')
plt.ylabel('Speed, km/s')
plt.legend()
plt.semilogy()
# plt.savefig('plots/sigmas.pdf')

r = stdev['mean'] / stdev['sigmax']
plt.figure(4)
plt.xlabel('log($M_{halo}$ / $M_\odot$)')
plt.ylabel('$v_{3D}$ / $\sigma_{los}$')
plt.plot(stdev['mass'], r)
plt.plot([9,15], [np.pi/2, np.pi/2], '--')
# plt.savefig('plots/pi.pdf')


# stdevx = satsbinned['bin', 'vpecx'].groups.aggregate(np.std)
# stdevx.rename_column('vpecx', 'sigma_x')
# stdevx['mean'] = satsbinned['vpecx'].groups.aggregate(np.mean)
# stdevx['median'] = satsbinned['vpecx'].groups.aggregate(np.median)
# stdevx['mass'] = left + (stdevx['bin'])
#
# plt.figure(3)
# plt.plot(sats['centralMVirlog'][::10], sats['vpecx'][::10], ',')

# plt.figure(3)
# plt.hist(satsbinned['vpec'].groups[25])