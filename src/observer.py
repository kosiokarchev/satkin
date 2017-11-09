from lib import *


class Regressor():
    _regressor = None

    @classmethod
    def extract(cls, sns=SNDATA, params=('amp', 'err_amp', 'exp', 'err_exp')):
        cols = {'snapnum': sns}
        for p in params:
            cols[p] = []
        for sn in sns:
            r = Table.read('data/reg-mvir{}.csv'.format(sn),
                           format='ascii.ecsv')
            for p in params:
                cols[p].append(r.meta[p])

        t = Table(data=cols,
                  names=['snapnum'] + list(params))

        return t

    fname = 'data/reg-mvir.csv'
    format = 'ascii.ecsv'

    @classmethod
    def get(cls):
        if Regressor._regressor is None:
            try:
                t = Table.read(Regressor.fname, format=Regressor.format)
            except FileNotFoundError:
                t = Regressor.extract()
                t.write(Regressor.fname, format=Regressor.format)
            Regressor._regressor = cls(t)

        return Regressor._regressor

    def __init__(self, t):
        self.t = t
        self.sns = self.t['snapnum'].tolist()

    def __getitem__(self, item):
        return self.t[self.sns.index(item)]


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
        t = join(self.mags, t, 'galaxyId')
        return t[np.where(t['ACS775'] < 27)]