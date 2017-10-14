#!/usr/bin/env python3
import subprocess


def download(sql, fname):
    with open('credentials.txt', encoding='ascii') as f:
        usr, pwd = f.read().replace('\n', '').split(':')
    url = 'http://gavo.mpa-garching.mpg.de/MyMillennium/?action=doQuery&SQL='+sql

    wget = ('wget', '--http-user=' + usr, '--http-passwd=' + pwd, '-O', fname, url)
    subprocess.run(wget, check=True)


def gen_select(cols, table, where='', top=0):
    return 'SELECT {top} {cols} FROM {table} {where}'.format(
        top=('TOP '+str(top)) if top else '',
        cols=', '.join(cols), table=table,
        where=('WHERE '+where) if where else '')


TABLES = {'H15-cube': 'Henriques2015a..MRscPlanck1',
          'H15-cone': lambda n=1: 'Henriques2015a.cones.MRscPlanck1_BC03_{:0>3}'.format(n)}
COLUMNS = {
    'H15-cube': [
        'galaxyId',
        'phKey', 'x', 'y', 'z', 'velX', 'velY', 'velZ',
        'redshift', 'lookBackTime',
        'rvir', 'bulgeSize', 'stellarDiskRadius',
        'mvir', 'coldGas', 'stellarMass', 'bulgeMass', 'diskMass', 'hotGas', 'ejectedMass', 'blackHoleMass', 'icmStellarMass',
        'sfr', 'sfrBulge',
        'type', 'centralMVir', 'centralRvir', 'distanceToCentralGalX', 'distanceToCentralGalY', 'distanceToCentralGalZ',
        'treeId', 'descendantId', 'mainLeafId', 'treeRootId', 'firstProgenitorId', 'nextProgenitorId', 'lastProgenitorId', 'haloId', 'subHaloId', 'fofCentralId', 'fofSubhaloId'],
    'H15-cone': [
        'galaxyId',
        'ra', 'dec', 'inclination', 'PA',
        'z_geo', 'z_app', 'vpec',
        'd_comoving', 'dlum', 'ang_dist'
    ]
}


if __name__ == '__main__':
    download(gen_select(COLUMNS['H15-cube'], TABLES['H15-cube'], top=10000), 'test.csv')
