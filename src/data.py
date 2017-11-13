from lib import *


class Data:
    def __init__(self,
                 mvir=None, mstar=None,
                 err_mvir=None, err_mstar=None,
                 p16=None, p84=None,
                 N=None,
                 label=None):
        self.mvir = mvir
        self.mstar = mstar
        self.err_mvir = err_mvir
        self.err_mstar = err_mstar
        self.mvir_p16 = p16
        self.mvir_p84 = p84
        self.N = N
        self.label = label

    def get_ebars(self):
        if self.err_mvir is None:
            return (self.mvir-self.mvir_p16, self.mvir_p84-self.mvir)
        else:
            return self.err_mvir

    def plot(self, **kwargs):
        kw = dict(errorevery=7, marker='.', markersize=5, elinewidth=2, capsize=2)
        kw.update(kwargs)
        plt.errorbar(self.mvir, self.mstar, xerr=self.get_ebars(), **kw)