from __future__ import print_function
import numpy as np
from astropy.table import Table, join

def combine(n):
    print('Loading tables ' + str(n) + ' ...')
    centrals = Table.read('data/cube'+str(n)+'=0.fits', format = 'fits')
    sats = Table.read('data/cube'+str(n)+'>0.fits', format = 'fits')

    print('Joining ' + str(n) + ' ...')
    centrals.rename_column('galaxyId', 'fofCentralId')
    sats = join(sats,centrals,'fofCentralId',table_names=['s', 'c']) #column names like name_s/name_c

    print('Calculating velocities and distances ' + str(n) + ' ...')
    sats['vpecx'] = sats['velX_s'] - sats['velX_c']
    sats['vpecy'] = sats['velY_s'] - sats['velY_c']
    sats['vpecz'] = sats['velZ_s'] - sats['velZ_c']
    sats['vpec'] = np.sqrt(sats['vpecx']**2+sats['vpecy']**2+sats['vpecz']**2) #magn of peculiar velocity

    satellites = sats['mvir', 'stellarMass', 'vpecx', 'vpecy', 'vpecz', 'vpec', 'distanceToCentralGalX', 'distanceToCentralGalY', 'distanceToCentralGalZ']
    satellites.rename_column('distanceToCentralGalX', 'rx')
    satellites.rename_column('distanceToCentralGalY', 'ry')
    satellites.rename_column('distanceToCentralGalZ', 'rz')
    satellites['r'] = np.sqrt(satellites['rx']**2+satellites['ry']**2+satellites['rz']**2)

    print('Writing tables ' + str(n) + ' ...')
    satellites.write('data/sats'+str(n)+'.fits', format = 'fits')

combine(34)
combine(38)
combine(45)