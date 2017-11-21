from funcs import *
from scipy.special import erf
import numpy.random as random


def p(x):
    y = 0.5*(1+erf(-(x-26)/2))
    return y

class Observer:
    _observer = None

    @classmethod
    def get(cls):
        if Observer._observer is None:
            Observer._observer = Observer()
        return Observer._observer

    def __init__(self):
        print('Preparing observer:')
        print('Loading magnitudes...')
        self.mags = Table.read('data/cones/cone.fits')
        print('Observer prepared.')

    def observe(self, t):
        print('Observing')

        self.mags['rand'] = random.random(len(self.mags))
        self.mags['func'] = p(self.mags['ACS775'])
        sample = self.mags[np.where(self.mags['func'] > self.mags['rand'])]
        t = join(sample, t, 'galaxyId')

        return t