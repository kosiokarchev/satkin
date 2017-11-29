from funcs import *
from scipy.special import erf
import numpy.random as random


def prob_erf(x, m, s, A, H):
    return H * 0.5*(1+erf((x-m) / s)) + A
prob_erf_params = (25.718585789217027,
                   -1.7715191963066255,
                   0.044643065344693235,
                   0.9564014441048382)
def prob_fd(x, m, s, A, H):
    return H / (1+np.exp((x-m) / s)) + A
prob_fd_params = (25.708477496243546,
                  0.7560238533784606,
                  0.04035120415731956,
                  0.9670816213373163)
prob_func = lambda x: prob_erf(x, *prob_erf_params)

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
        prob = prob_func(self.mags['ACS775'])
        self.sample = self.mags[np.where(prob > rand)]

        print('Observer prepared.')

    def observe(self, t):
        print('Observing')
        return join(self.sample, t, 'galaxyId')