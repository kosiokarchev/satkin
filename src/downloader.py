from __future__ import print_function
import os
from time import time, sleep
import threading
import math
import requests

SI_PREFIXES = ('n', 'u', 'm', '', 'k', 'M', 'G', 'T')
def SI(n, r=3):
    if n==0:
        return '0 '
    p10 = int(math.floor(math.log10(abs(n))))
    p = int(math.floor(p10 / 3))
    m = round(n, r-1-p10) / 1000**p
    if m.is_integer():
        m = int(m)
    pr = SI_PREFIXES[p+3] if (-3 <= p <= 4) else '10^' + str(3*p)

    return str(m) + ' ' + pr


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
    def __init__(self, query, out, index='galaxyId', isline=None):
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
        self.nfiles = self.nlines = self.nchars = self.datalines = self.datachars = 0

        self.isline = isline
        if self.isline is None:
            self.isline = lambda line: not line[0] == '#'

        self.sess = None
        self.overwritten = False
        self.ok = self.error = False

        self.writeout_count = 10000

    @staticmethod
    def load_pass():
        pf = os.path.join(os.path.dirname(__file__), 'credentials.txt')
        with open(pf, encoding='ascii') as f:
            usr, pwd = f.read().replace('\n', '').split(':')
        return usr, pwd

    def report(self):
        bigsep = '='*80
        smallsep = '-'*80

        start_time = time()
        prevchars = self.nchars
        prevtime = start_time

        print(self)
        print('\u001b[?25l', end='', flush=True)

        while True:
            sleep(1)

            dt = time() - prevtime
            prevtime = time()

            dchars = self.nchars - prevchars
            prevchars = self.nchars
            v = dchars / dt

            elapsed_file = prevtime - self.file_start_time
            elapsed = prevtime - start_time
            vavg = self.nchars / elapsed

            msg = '[{:.0f}s] ({}: [{:.0f}s]) {}c in {:,} lines ({}c/s, avg: {}c/s) --> {}'.format(
                elapsed,
                self.nfiles, elapsed_file,
                SI(self.datachars), self.datalines,
                SI(v), SI(vavg),
                self.maxindex
            )

            print('\r\u001b[0K' + msg, end='', flush=True)

            if self.ok:
                break
            if self.error:
                print('\u001b[?25h')
                return

        elapsed = time() - start_time

        print('\u001b[?25h')
        print('\u001b[32m'+'Download complete!'+'\u001b[m')
        print('{}c ({:,} lines) in [{:.0f}s] ({}c/s)'.format(
            SI(self.datachars), self.datalines, elapsed, SI(self.nchars/elapsed)
        ))
        print(bigsep)

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
                if self.isline(line):
                    f.write(line+'\n')
                    self.maxindex = line.split(',')[self.indexindex]

                    self.datalines += 1
                    self.datachars += len(line)

                self.nlines += 1
                self.nchars += len(line)

        if line == '#OK':
            self.ok = True

    def go(self):
        if os.path.exists(self.out):
            open(self.out, 'w').close() # Just empty the file, no delete
            self.overwritten = True

        rep_thread = threading.Thread(target=self.report)
        rep_thread.start()

        try:
            with requests.Session() as self.sess:
                self.sess.auth = Downloader.load_pass()
                self.sess.stream = True

                while not self.ok:
                    self.download_one()

                rep_thread.join()
        except BaseException as e:
            self.error = True
            rep_thread.join()
            raise e


    def __repr__(self):
        s = '\u001b[31m' + 'Overwriting' + '\u001b[m' if self.overwritten else 'Saving to'
        s += ' \u001b[1m' + self.out + '\u001b[m' + '\n'
        s += '-'*80 + '\n'
        s += '\u001b[3m' + self.baseQuery + '\u001b[m' +'\n'
        s += '-'*80
        return s


def csv2fits(fname, outname, cols):
    from astropy.io import ascii
    print('Loading', fname)
    t = ascii.read(fname, format='csv', names=cols, data_start=0,
                   fast_reader={'parallel': True})
    print('Loaded', fname)
    for c in t.colnames:
        if 'Id' in c:
            t[c] = t[c].astype(int)
    t.write(outname, format='fits', overwrite=True)
    print('Written', outname)
