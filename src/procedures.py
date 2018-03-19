from funcs import *
from binner import binX
from observer import Observer
from data import Data


__all__ = (
    'load_sats',
    'procs', 'procnames',
    'SWProcedure',
    'HWAProcedure', 'HWMProcedure', 'HWNProcedure', 'HWPProcedure'
)


def load_sats(sn):
    t = load(FILES['sats'](sn))
    t = t[np.where(t['stellarMass'] > 1e-4)]
    t['mvir'] = 10 + np.log10(t['mvir'])
    t['stellarMass'] = 10 + np.log10(t['stellarMass'])
    return t


class Procedure:
    stdev_file = None
    reg_file = None

    def __init__(self, sn, observe):
        self.sn = sn
        self.observe = observe

        self.sats = None
        self.stdev = None
        self.regression = None
        self.predictions = None

        if self.stdev_file is not None:
            self.stdev_file = self.stdev_file.__func__(sn, observe) # Python 2....

    def load_sats(self):
        if self.sats is None:
            t = load_sats(self.sn)

            if self.observe:
                t = Observer.get().observe(t)

            self.sats = t

    def load_stdev(self):
        if self.stdev is None:
            try:
                self.stdev = load(self.stdev_file)
            except IOError:
                print('Table does not exist! Calculating it.')
                self.calculate_stdev(True)

    def write_stdev(self):
        print('Writing to', self.stdev_file)
        self.stdev.write(self.stdev_file, overwrite=True)

    def get_subset(self):
        raise NotImplementedError

    def regress(self):
        raise NotImplementedError

    def calculate_stdev(self, write=True):
        raise NotImplementedError

    def predict(self):
        self.load_stdev()
        p = self.get_regression()

        mvir = np.log10(self.stdev['sigmax'] / p['amp']) / p['exp']

        relerr_mvir = np.sqrt(
            ((1 / p['exp']) * (p['err_amp'] / p['amp'])) ** 2
            + (np.log(mvir) * (p['err_exp'] / p['exp'])) ** 2
        )

        self.predictions = Data(mvir=mvir, mstar=self.stdev['stellarMass'],
                                err_mvir=relerr_mvir / np.log(10),
                                N=self.stdev['N'])

    def get_predictions(self):
        if self.predictions is None:
            self.predict()
        return self.predictions

    def write_regression(self):
        try:
            with open(self.reg_file) as f:
                reg = json.load(f)
        except IOError:
            reg = dict()
        reg[str(self.sn)] = self.regression
        with open(self.reg_file, 'w') as f:
            f.write(json.dumps(reg))

    def get_regression(self):
        if self.regression is None:
            try:
                with open(self.reg_file) as f:
                    self.regression = json.load(f)[str(self.sn)]
            except IOError:
                if self.observe:
                    self.regression = type(self)(self.sn,
                                                 False).get_regression()
                else:
                    print('The regression does not exist! Calculating it...')
                    self.regress()
        return self.regression


class SWProcedure(Procedure):
    reg_file = FILES['sw-reg']
    stdev_file = FILES['sw-sigmas']

    def __init__(self, sn, observe):
        Procedure.__init__(self, sn, observe)

        self.mvir = None
        self.mvir_file = FILES['sw-sigmas'](sn, observe, 'mvir')

    def load_mvir(self):
        if self.mvir is None:
            try:
                self.mvir = load(self.mvir_file)
            except IOError:
                print('Table does not exist! Calculating it...')
                self.bin(PLOTDATA['bins']['mvir'][0], BINWIDTH, 'mvir',
                         self._store_mvir, self.write_mvir)

    def write_mvir(self):
        write(self.mvir, self.mvir_file)

    def get_subset(self):
        t = load(FILES['nums'](self.sn))
        return t[np.where(t['num_sat'] > 2)]

    def _store_stdev(self, stdev):
        self.stdev = stdev

    def _store_mvir(self, mvir):
        self.mvir = mvir

    def bin(self, left, width, mname, store, write=None):
        self.load_sats()

        self.sats['bin'] = np.floor_divide(self.sats[mname] - left, width)
        binned = self.sats.group_by('bin')
        self.sats.remove_column('bin')

        res = binned['bin', 'vpec', 'vpecx', 'vpecy', 'vpecz'].groups.aggregate(
            np.std)

        res.rename_column('vpec', 'sigma')
        res.rename_column('vpecx', 'sigmax')
        res.rename_column('vpecy', 'sigmay')
        res.rename_column('vpecz', 'sigmaz')

        res['mean'] = binned['vpec'].groups.aggregate(np.mean)
        res['median'] = binned['vpec'].groups.aggregate(np.median)
        res[mname] = left + (res['bin'] + 0.5) * width
        res['N'] = binned['fofCentralId'].groups.aggregate(
            lambda x: len(np.unique(x)))

        store(res)

        if callable(write):
            write()

    def regress(self, write=True):
        print('Regressing...')
        self.load_mvir()

        print('Fitting to {} points'.format(len(self.mvir)))
        indices, fit = fit_or(models.PowerLawModel(), self.mvir['sigmax'],
                              x=np.power(10, self.mvir['mvir']))

        self.regression = {
            'exp': fit.params['exponent'].value,
            'err_exp': fit.params['exponent'].stderr,
            'amp': fit.params['amplitude'].value,
            'err_amp': fit.params['amplitude'].stderr,
            'indices': indices.tolist()
        }

        if write:
            self.write_regression()

    def calculate_stdev(self, write=True):
        print('Calculating stdev')
        self.bin(PLOTDATA['bins']['stellarMass'][0], BINWIDTH, 'stellarMass',
                 self._store_stdev, self.write_stdev if write else False)

    def plot(self, mname, xlim, xlabel='$\log_{10}(M / M_\odot)$', sparse=333):
        print('Plotting...')
        fig, ax = plt.subplots(2, 1, sharex='col', squeeze=True,
                               gridspec_kw={'height_ratios': [3, 1],
                                            'hspace': 0})

        plt.sca(ax[0])
        r = np.random.randint(0, len(self.sats), int(len(self.sats) / sparse))
        plt.plot(self.sats[mname][r], self.sats['vpec'][r], ',',
                 label='_nolegend_')

        X = self.stdev[mname]
        Y = self.stdev
        plt.plot(X, Y['sigma'], '-r', linewidth=3, label='$\sigma_{3D}$')
        plt.plot(X, Y['sigmax'], '-y', linewidth=2, label='$\sigma_{los}$')
        plt.plot(X, Y['sigmay'], '-y', linewidth=2, label='_nolegend_')
        plt.plot(X, Y['sigmaz'], '-y', linewidth=2, label='_nolegend_')
        plt.plot(X, Y['mean'], '-g', linewidth=3,
                 label='$\langle |v| \\rangle$')
        plt.semilogy()

        plt.legend(loc='lower right')

        plt.ylabel('Speed, km/s')

        print('Plotting...')
        plt.sca(ax[1])
        plt.plot(X, Y['mean'] / Y['sigmax'])
        plt.plot(xlim, [np.sqrt(3)] * 2, ':')
        plt.plot(xlim, [np.sqrt(8 / np.pi)] * 2, ':')

        plt.xlabel(xlabel)
        plt.ylabel('$v_{3D}$ / $\sigma_{los}$')

        plt.ylim(1, 2)
        plt.xlim(xlim)


class HWProcedure(Procedure):
    reg_file = FILES['hw-reg']

    def __init__(self, sn, observe):
        Procedure.__init__(self, sn, observe)

        self.rms = None
        self.rms_file = FILES['hw-rms'](sn, observe)

    def bin(self, left, width, write=True):
        raise NotImplementedError

    def load_rms(self, write=True):
        if self.rms is None:
            try:
                self.rms = load(self.rms_file)
            except IOError:
                print('The table does not exist! Calculating it...')
                self.calculate_rms(write)

    def write_rms(self):
        write(self.rms, self.rms_file)

    def get_subset(self):
        self.load_rms()
        return self.rms

    def calculate_rms(self, write=True):
        print('Calculating rms')
        self.load_sats()

        print('Grouping by fofCentralId...')
        b = self.sats.group_by('fofCentralId')

        cols = ['vpec', 'vpecx', 'vpecy', 'vpecz']
        for col in cols:
            b[col] = b[col] ** 2

        print('Calculating mean square...')
        rms = b[['fofCentralId'] + cols].groups.aggregate(np.mean)

        for col in cols:
            rms[col] = np.sqrt(rms[col])
            rms.rename_column(col, col.replace('vpec', 'rms'))

        print('Joining with centrals data...')
        nums = load(FILES['nums'](self.sn))
        rms = join(rms, nums, keys='fofCentralId')
        self.rms = rms

        if write:
            self.write_rms()

    def regress(self, write=True):
        print('Regressing...')
        self.load_rms()

        t = deal(self.rms)
        print('Fitting to {} points'.format(len(t)))
        indices, fit = fit_or(models.PowerLawModel(), t['rmsx'],
                              x=np.power(10, t['mv']))

        self.regression = {
            'exp': fit.params['exponent'].value,
            'err_exp': fit.params['exponent'].stderr,
            'amp': fit.params['amplitude'].value,
            'err_amp': fit.params['amplitude'].stderr
        }

        if write:
            self.write_regression()

    def calculate_stdev(self, write=True):
        print('Calculating stdev')
        self.load_rms(write)
        try:
            self.bin(PLOTDATA['bins']['stellarMass'][0], BINWIDTH, write)
        except NotImplementedError:
            pass


class HWSigmaProcedure(HWProcedure):
    @staticmethod
    def collapse(binned, cols):
        raise NotImplementedError

    def bin(self, left, width, write=True):
        self.load_rms()

        t = deal(self.rms)
        t['bin'] = np.floor_divide(t['ms'] - left, width)
        binned = t.group_by('bin')

        cols = ['rms', 'rmsx', 'rmsy', 'rmsz']
        self.stdev = self.collapse(binned, cols)
        for col in cols:
            self.stdev.rename_column(col, col.replace('rms', 'sigma'))
        self.stdev['stellarMass'] = left + (self.stdev['bin'] + 0.5) * width
        self.stdev['N'] = binned.groups.indices[1:] - binned.groups.indices[:-1]

        if write:
            self.write_stdev()


class HWPProcedure(HWSigmaProcedure):
    stdev_file = FILES['hwp-sigmas']

    @staticmethod
    def collapse(binned, cols):
        for col in cols:
            binned[col] = binned[col] ** 2

        res = binned[['bin'] + cols].groups.aggregate(np.mean)

        for col in cols:
            res[col] = np.sqrt(res[col])

        return res


class HWAProcedure(HWSigmaProcedure):
    stdev_file = FILES['hwa-sigmas']

    @staticmethod
    def collapse(binned, cols):
        return binned[['bin'] + cols].groups.aggregate(np.mean)


class HWMProcedure(HWSigmaProcedure):
    stdev_file = FILES['hwm-sigmas']

    @staticmethod
    def collapse(binned, cols):
        return binned[['bin'] + cols].groups.aggregate(np.median)


class HWNProcedure(HWProcedure):
    stdev_file = None

    def bin(self, *args, **kwargs):
        raise NotImplementedError

    def predict(self):
        print('Predicting')
        p = self.get_regression()
        self.load_rms()

        t = deal(self.rms['stellarMass', 'mvir', 'rmsx'])
        t['mv'] = np.log10(t['rmsx'] / p['amp']) / p['exp']

        b = binX(t, 'ms', 'mv', bins=PLOTDATA['bins']['stellarMass'],
                 funcsx={'mean': np.mean, 'N': len},
                 funcsy={
                     'mean': np.mean, 'median': np.median, 'std': np.std,
                     'p16': percentile_wrapper(16),
                     'p84': percentile_wrapper(84)})

        self.predictions = {
            'median': Data(mvir=b[1]['median'], mstar=b[0]['mean'],
                           mvir_p16=b[1]['p16'], mvir_p84=b[1]['p84'],
                           N=b[0]['N']),
            'mean': Data(mvir=b[1]['mean'], mstar=b[0]['mean'],
                         err_mvir=b[1]['std'], N=b[0]['N'])
        }


procs = SWProcedure, HWAProcedure, HWMProcedure, HWNProcedure, HWPProcedure
procnames = 'sw', 'hwa', 'hwm', 'hwn', 'hwp'