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
        self.mags = load_table('data/cones/cone.fits')

        rand = random.random(len(self.mags))
        func = p(self.mags['ACS775'])
        self.sample = self.mags[np.where(func > rand)]

        print('Observer prepared.')

    def observe(self, t):
        print('Observing')
        return join(self.sample, t, 'galaxyId')