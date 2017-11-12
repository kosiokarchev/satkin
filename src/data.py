from lib import *


class Data:
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