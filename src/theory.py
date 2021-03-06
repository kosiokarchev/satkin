from multiprocessing import Pool

import numpy as np
from scipy.stats import poisson
from astropy.table import join, vstack


class Poissonizer:
    def __init__(self, t, bincols, axis, verbose=True):
        self.t = t
        self.m = None
        self.bincols = bincols
        self.axis = axis

        self.calced = self.collapsed = self.predicted = self.meaned = False

        self.verbose = verbose

    def calc_points(self):
        self.print('Getting counts')
        self.t = self.t.group_by(self.bincols + [self.axis]).groups
        self.t.keys['n'] = self.t.indices[1:] - self.t.indices[:-1]
        self.t = self.t.keys

        self.calced = True
        return self

    def collapse(self):
        if not self.calced:
            self.calc_points()

        self.print('Getting totals')
        N = self.t[self.bincols + ['n']].group_by(self.bincols).groups.aggregate(np.sum)
        N.rename_column('n', 'N')

        self.t = join(self.t, N, self.bincols, 'left')
        self.t['f'] = self.t['n'] / self.t['N']

        self.collapsed = True
        return self

    def set_collapsed(self, t):
        self.t = t
        self.collapsed = True
        return self

    def predict(self):
        if not self.collapsed:
            self.collapse()

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

        self.m = self.t.group_by(self.bincols).groups
        means = []
        medians = []
        stdevs = []
        for g in self.m:
            mean = np.sum(g['f']*g[self.axis])
            means.append(mean)
            medians.append(g[self.axis][np.searchsorted(np.cumsum(g['f']), 0.5)])
            stdevs.append(np.sqrt(np.sum(g['f'] * (g[self.axis] - mean)**2)))
        self.m = self.m.keys
        self.m['mean'] = means
        self.m['median'] = medians
        self.m['stdev'] = stdevs

        self.meaned = True
        return self

    def print(self, msg):
        if self.verbose:
            print(msg)


class BootstrapPoissonizer(Poissonizer):
    def __init__(self, observe, bincols, axis, nproc=None, ntrials=10):
        super(BootstrapPoissonizer, self).__init__(None, bincols, axis, False)

        self.observe = observe
        self.nproc = nproc
        self.ntrials = ntrials

    def _calc_points_one(self, i, seed):
        np.random.seed(seed)
        print('=' * 32, 'BOOTSTRAP {:0>4}'.format(i), '=' * 32)

        self.t = self.observe()
        super(BootstrapPoissonizer, self).calc_points()
        return self.t

    def calc_points(self):
        print('Bootstrap: nproc={}, ntrials={}'.format(self.nproc, self.ntrials))

        with Pool(self.nproc) as p:
            g = p.starmap(self._calc_points_one,
                          tuple(zip(
                              range(self.ntrials),
                              np.uint32(np.random.random(self.ntrials) * 2**32))
                          ))

        print('Bootstrap complete. Stacking', len(g), 'tables...')
        self.t = (vstack(g)
                  .group_by(self.bincols + [self.axis])
                  .groups.aggregate(np.sum))

        self.calced = True
        return self
