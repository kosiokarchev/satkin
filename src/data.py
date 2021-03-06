from funcs import *


class Data:
    def __init__(self, mvir=None, mstar=None, err_mvir=None, err_mstar=None,
                 mvir_p16=None, mvir_p84=None, N=None, label=None, **kwargs):
        self.mvir = np.array(mvir)
        self.mstar = np.array(mstar)
        self.err_mvir = np.array(err_mvir) if err_mvir is not None else None
        self.err_mstar = np.array(err_mstar) if err_mstar is not None else None
        self.mvir_p16 = np.array(mvir_p16) if mvir_p16 is not None else None
        self.mvir_p84 = np.array(mvir_p84) if mvir_p84 is not None else None
        self.N = np.array(N) if N is not None else None
        self.label = label

    def save(self, fname):
        d = {
            key: self.__dict__[key].tolist()
            if isinstance(self.__dict__[key], np.ndarray)
            else self.__dict__[key]
            for key in self.__dict__
        }
        with open(fname, 'w') as f:
            json.dump(d, f)

    @classmethod
    def load(cls, *args):
        if len(args)==1:
            fname = args[0]
        else:
            fname = FILES['data'](*args)

        with open(fname) as f:
            d = json.load(f)
        return cls(**d)


    def get_ebars(self):
        if self.err_mvir is None:
            return self.mvir-self.mvir_p16, self.mvir_p84-self.mvir
        else:
            return self.err_mvir

    def plot(self, **kwargs):
        kw = dict(errorevery=7, elinewidth=1, capsize=2, linewidth=2)
        kw.update(kwargs)
        plt.errorbar(self.mvir, self.mstar, xerr=self.get_ebars(), **kw)