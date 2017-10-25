#!/usr/bin/env python3
import os
from time import time
import requests


class Query:
    def __init__(self, cols, table, where=None, order=None):
        self.cols = list(cols)
        self.table = table
        self.where = where if where is not None else ''
        self.order = list(order) if order is not None else []

    def __call__(self):
        return 'SELECT {cols}\nFROM {table}{where}{order}'.format(
            cols = ','.join(self.cols),
            table = self.table,
            where = '\nWHERE '+self.where if self.where else '',
            order = '\nORDER BY ' + ','.join(self.order) if self.order else ''
        )


class Downloader:
    def __init__(self, query, out, index='galaxyId'):
        self.query = query
        self.out = os.path.abspath(out)
        self.index = index

    @staticmethod
    def load_pass():
        pf = os.path.join(os.path.dirname(__file__), 'credentials.txt')
        with open(pf, encoding='ascii') as f:
            usr, pwd = f.read().replace('\n', '').split(':')
        return usr, pwd

    def go(self):
        baseUrl = 'http://gavo.mpa-garching.mpg.de/MyMillennium/'

        if self.index not in self.query.cols:
            self.query.cols.append(self.index)
        if not (self.query.order and self.query.order[0] == self.index):
            self.query.order = [self.index+' OPTION (FAST 100000)'] + self.query.order

        indexindex = self.query.cols.index(self.index)
        maxindex = 0
        basewhere = self.query.where

        with requests.Session() as sess:
            sess.auth = Downloader.load_pass()
            sess.stream = True

            print('='*80)
            print('='*80)
            print(self)
            print('='*80)
            print('='*80)

            while True:
                self.query.where = '{}>{} AND ({})'.format(self.index, maxindex, basewhere)

                q = self.query()

                start_time = time()

                r = sess.get(baseUrl, params={'action': 'doQuery', 'SQL': q})
                r.raise_for_status()

                lines = r.iter_lines(decode_unicode=True)

                line = next(lines)
                if not line == '#OK':
                    raise RuntimeError('Query not OK')

                nlines = 1
                nchars = len(line)
                with open(self.out, 'a') as f:
                    for line in lines:
                        if line[0] == '#':
                            pass
                            # print('Ignoring as comment:', line)
                        else:
                            f.write(line + '\n')
                            maxindex = line.split(',')[indexindex]

                        nlines += 1
                        nchars += len(line)
                        dt = time() - start_time


                        if nlines % 10000 == 0:
                            msg = '({:,.0f} lines)'.format(nlines)
                            msg += ' [{:.0f}s]'.format(dt)

                            if nchars < 1000000:
                                msg += ' {:.0f} kc'.format(nchars / 1000)
                            else:
                                msg += ' {:.0f} Mc'.format(nchars / 1000000)

                            print(msg+'\u001b[0K\r', end='', flush=True)

                print()
                print('Last line is', line)
                print('maxindex is', maxindex)
                print('=' * 80)

                if line == '#OK':
                    break

        print('DOWNLOAD COMPLETE!!')
        print('Data saved in', self.out)


    def __repr__(self):
        return self.query() + '\n--> ' + self.out


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
    cols_central = [
        'galaxyId',
        'x', 'y', 'z', 'velX', 'velY', 'velZ', 'stellarSpinX', 'stellarSpinY',
        'stellarSpinZ',
        'mvir', 'rvir', 'stellarMass', 'sfr',
    ]
    cols_satellites = [
        'galaxyId', 'fofCentralId',
        'velX', 'velY', 'velZ',
        'distanceToCentralGalX', 'distanceToCentralGalY', 'distanceToCentralGalZ'
    ]

    snapnums = (38, 34, 45)

    for snapnum in snapnums:
        for typewhere, cols in zip(('=0', '>0'), (cols_central, cols_satellites)):
            where = 'snapnum={} AND type{}'.format(snapnum, typewhere)
            q = Query(cols, TABLES['H15-cube'], where)

            fname = 'data/cube{}{}.csv'.format(snapnum, typewhere)

            d = Downloader(q, fname)
            d.go()
