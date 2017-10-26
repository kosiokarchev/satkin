#!/usr/bin/env python3
import os
from time import time, sleep
import threading
import math
import requests

SI_PREFIXES = ('n', 'u', 'm', '', 'k', 'M', 'G', 'T')
def SI(n):
    if n==0:
        return '0 '
    p = int(math.floor(math.log(abs(n), 1000)))
    m = n / 1000**p
    pr = SI_PREFIXES[p+3] if (-3 <= p <= 4) else '10^'+str(3*p)
    return '{:.0f} {}'.format(m, pr)


class Query:
    def __init__(self, cols, table, where=None, order=None):
        self.cols = list(cols)
        self.table = table
        self.where = where if where is not None else ''
        self.order = list(order) if order is not None else []

    def __call__(self):
        return 'SELECT TOP 100000 {cols}\nFROM {table}{where}{order}'.format(
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

        if self.index not in self.query.cols:
            self.query.cols.append(self.index)
        if not (self.query.order and self.query.order[0] == self.index):
            self.query.order = [self.index+' OPTION (FAST 100000)'] + self.query.order

        self.indexindex = self.query.cols.index(self.index)
        self.maxindex = 0
        self.basewhere = self.query.where

        self.baseUrl = 'http://gavo.mpa-garching.mpg.de/MyMillennium/'
        self.baseQuery = self.query()

        self.file_start_time = None
        self.nfiles = self.nlines = self.nchars = 0

        self.sess = None
        self.overwritten = False
        self.ok = False

        self.writeout_count = 10000

    @staticmethod
    def load_pass():
        pf = os.path.join(os.path.dirname(__file__), 'credentials.txt')
        with open(pf, encoding='ascii') as f:
            usr, pwd = f.read().replace('\n', '').split(':')
        return usr, pwd

    def download_one(self):
        self.query.where = '{}>{} AND ({})'.format(self.index, self.maxindex, self.basewhere)

        self.file_start_time = time()
        r = self.sess.get(self.baseUrl, params={'action': 'doQuery',
                                                'SQL': self.query()})
        r.raise_for_status()

        lines = r.iter_lines(decode_unicode=True)

        line = next(lines)
        if not line == '#OK':
            raise RuntimeError('Query not OK')

        self.nfiles += 1
        self.nlines += 1
        self.nchars += len(line)

        with open(self.out, 'a', ) as f:
            for line in lines:
                if line[0] != '#':
                    f.write(line+'\n')
                    self.maxindex = line.split(',')[self.indexindex]

                self.nlines += 1
                self.nchars += len(line)

        if line == '#OK':
            self.ok = True

    def report(self):
        start_time = time()
        prevchars = self.nchars
        prevtime = start_time

        print('\u001b[31m'+'Overwriting'+'\u001b[m' if self.overwritten else 'Saving to',
              '\u001b[1m'+self.out+'\u001b[m')
        print('-'*80)
        print('\u001b[3m'+self.baseQuery+'\u001b[m')
        print('-'*80)
        print('\u001b[?25l', end='', flush=True)

        while True:
            sleep(1)

            dchars = self.nchars - prevchars
            prevchars = self.nchars
            dt = time() - prevtime
            prevtime = time()

            v = dchars / dt

            elapsed_file = prevtime - self.file_start_time
            elapsed = prevtime - start_time
            vavg = self.nchars / elapsed

            msg = '[{:.0f}s]'.format(elapsed)
            msg += ' ({}: [{:.0f}s])'.format(self.nfiles, elapsed_file)
            msg += ' {:,.0f} lines, {}c'.format(self.nlines, SI(self.nchars))
            msg += ' ({}c/s, avg: {}c/s)'.format(SI(v), SI(vavg))
            msg += ' --> {}'.format(self.maxindex)

            print('\r\u001b[0K' + msg, end='', flush=True)

            if self.ok:
                break

        print()
        print('\u001b[?25h', end='', flush=True)
        print('\u001b[32m'+'Download complete!'+'\u001b[m')
        print('({} lines) [{}s] {}'.format(self.nlines, time()-start_time, SI(self.nchars)))
        print('='*80)

    def go(self):
        try:
            os.remove(self.out)
            self.overwritten = True
        except OSError as e:
            if not e.errno == 2:
                raise e

        rep_thread = threading.Thread(target=self.report, daemon=True)
        rep_thread.start()

        with requests.Session() as self.sess:
            self.sess.auth = Downloader.load_pass()
            self.sess.stream = True

            while not self.ok:
                self.download_one()

        rep_thread.join()


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
    #
    # snapnums = (38, 34, 45)
    #
    # for snapnum in snapnums:
    #     for typewhere, cols in zip(('=0', '>0'), (cols_central, cols_satellites)):
    #         where = 'snapnum={} AND type{}'.format(snapnum, typewhere)
    #         q = Query(cols, TABLES['H15-cube'], where)
    #
    #         fname = 'data/cube{}{}.csv'.format(snapnum, typewhere)
    #
    #         d = Downloader(q, fname)
    #         d.go()
    q = Query(cols_central, TABLES['H15-cube'], 'snapnum=34 AND type=0')
    d = Downloader(q, 'data/cube34=0.csv')
    d.go()
