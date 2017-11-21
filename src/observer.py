from funcs import *
import numpy.random as random
# from astropy.table import Table

# class Regressor():
#     _regressor = None
#
#     @classmethod
#     def extract(cls, sns=SNDATA, params=('amp', 'err_amp', 'exp', 'err_exp')):
#         cols = {'snapnum': sns}
#         for p in params:
#             cols[p] = []
#         for sn in sns:
#             r = Table.read(FILES['sw-reg-mvir'](sn), format='ascii.ecsv')
#             for p in params:
#                 cols[p].append(r.meta[p])
#
#         t = Table(data=cols,
#                   names=['snapnum'] + list(params))
#
#         return t
#
#     fname = 'data/reg-mvir.csv'
#     format = 'ascii.ecsv'
#
#     @classmethod
#     def get(cls):
#         if Regressor._regressor is None:
#             try:
#                 t = Table.read(Regressor.fname, format=Regressor.format)
#             except FileNotFoundError:
#                 t = Regressor.extract()
#                 t.write(Regressor.fname, format=Regressor.format)
#             Regressor._regressor = cls(t)
#
#         return Regressor._regressor
#
#     def __init__(self, t):
#         self.t = t
#         self.sns = self.t['snapnum'].tolist()
#
#     def __getitem__(self, item):
#         return self.t[self.sns.index(item)]

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
        self.mags['rand'] = random.random(len(self.mags))
        self.mags['func'] = p(self.mags['ACS775'])
        sample = self.mags[np.where(self.mags['func'] > self.mags['rand'])]
        t = join(sample, t, 'galaxyId')
        return t