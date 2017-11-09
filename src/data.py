from lib import *


class Data:
    @classmethod
    def obtain_sw(cls, sn, observe=True):
        params = Regressor.get()[sn]
        sn, amp, err_amp, exp, err_exp = params

        sname = 'data{}/sigmas-stellarMass{}.fits'.format('/observe' if observe else '', sn)
        data = Table.read(sname)
        data['mass'] = np.power(10, data['mass'])
        data.rename_column('mass', 'mstar')

        data['mvir'] = np.power((data['sigmax'] / amp), 1 / exp)

        data['relerr_mvir'] = np.sqrt(
            ((1/exp) * (err_amp/amp))**2 + (np.log(data['mvir']) * (err_exp/exp))**2
        )

        return cls(mvir=data['mvir'], mstar=data['mstar'],
                   relerr_mvir=data['relerr_mvir'],
                   N=data['N'])

    @classmethod
    def obtain_hw(cls):
        mvir = [1e11, 2.3e12]
        mstar = [1.8e9, 3e10]

        return cls(mvir=mvir, mstar=mstar)

    def __init__(self,
                 mvir=None, mstar=None,
                 err_mvir=None, err_mstar=None,
                 relerr_mvir=None, relerr_mstar=None,
                 N=None,
                 label=None):
        self.mvir = mvir
        self.mstar = mstar
        self.err_mvir = err_mvir
        self.err_mstar = err_mstar
        self.relerr_mvir = relerr_mvir
        self.relerr_mstar = relerr_mstar
        self.N = N
        self.label = label