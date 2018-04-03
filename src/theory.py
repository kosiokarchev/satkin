from multiprocessing import Pool

import numpy as np
from scipy.stats import poisson
from astropy.table import join, vstack


class Poissonizer:
    def __init__(self, t, bincols, axis):
        self.t = t
        self.m = None
        self.bincols = bincols
        self.axis = axis

        self.calced = self.collapsed = self.predicted = self.meaned = False

    def calc_points(self):
        print('Getting counts')
        self.t = self.t.group_by(self.bincols + [self.axis]).groups
        self.t.keys['n'] = self.t.indices[1:] - self.t.indices[:-1]
        self.t = self.t.keys

        self.calced = True
        return self

    def collapse(self):
        if not self.calced:
            self.calc_points()

        print('Getting totals')
        N = self.t[self.bincols + ['n']].group_by(self.bincols).groups.aggregate(np.sum)
        N.rename_column('n', 'N')

        self.t = join(self.t, N, self.bincols, 'left')

        self.collapsed = True
        return self

    def set_collapsed(self, t):
        self.t = t
        self.collapsed = True
        return self

    def predict(self):
        if not self.collapsed:
            self.collapse()

        self.t['f'] = self.t['n'] / self.t['N']

        self.t = self.t.group_by(self.bincols)

        pois = []
        for keys, vals in zip(self.t.groups.keys, self.t.groups):
            pois.extend(poisson.pmf(vals[self.axis],
                                    np.sum(vals['f'] * vals[self.axis])))
        self.t['poisson'] = pois

        self.predicted = True
        return self

    def mean(self):
        if not self.collapsed:
            self.collapse()

        self.m = self.t.group_by(self.bincols)
        self.m.groups.keys['mean'] = self.m[self.axis].groups.aggregate(np.mean)
        self.m.groups.keys['stdev'] = self.m[self.axis].groups.aggregate(np.std)
        self.m = self.m.groups.keys

        self.meaned = True
        return self


class BootstrapPoissonizer(Poissonizer):
    def __init__(self, cen, cenp, bincols, axis, nproc=None, ntrials=10):
        super(BootstrapPoissonizer, self).__init__(None, bincols, axis)
        self.cen = cen
        self.cenp = cenp

        self.nproc = nproc
        self.ntrials = ntrials

        self.observed = False

    def observe(self):
        sats = self.cenp[self.cenp['p'] > np.random.random(len(self.cenp))]
        sats = sats[sats['galaxyId'] != sats['fofCentralId']]
        sats = sats.group_by('fofCentralId')['fofCentralId', 'galaxyId'].groups.aggregate(len)
        sats.rename_column('galaxyId', self.axis)

        t = join(self.cen[['fofCentralId'] + self.bincols], sats, 'fofCentralId', 'left')
        t[self.axis] = t[self.axis].filled(0)
        self.t = t

        self.observed = True

    def _calc_points_one(self, i, seed):
        np.random.seed(seed)
        print('=' * 32, 'BOOTSTRAP {:0>4}'.format(i), '=' * 32)
        super(BootstrapPoissonizer, self).calc_points()
        return self.t

    def calc_points(self):
        if not self.observed:
            self.observe()

        print('Bootstrap collapse: nproc={}, ntrials={}'.format(self.nproc, self.ntrials))

        with Pool(self.nproc) as p:
            g = p.starmap(self._calc_points_one,
                          zip(range(self.ntrials),
                              np.uint32(np.random.random(self.ntrials) * 2**32)))

        self.t = (vstack(g)
                  .group_by(self.bincols + [self.axis])
                  .groups.aggregate(np.sum))

        self.calced = True
        return self
