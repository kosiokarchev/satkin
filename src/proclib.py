from lib import *


def combine(sn):
    print('Loading tables ' + str(sn) + ' ...')
    centrals = Table.read('data/cube' + str(sn) + '=0.fits', format ='fits')
    sats = Table.read('data/cube' + str(sn) + '>0.fits', format ='fits')

    print('Joining ' + str(sn) + ' ...')
    centrals.rename_column('galaxyId', 'fofCentralId')
    sats = join(sats, centrals, 'fofCentralId', table_names=['s', 'c']) #column names like name_s/name_c

    print('Calculating velocities and distances ' + str(sn) + ' ...')
    sats['vpecx'] = sats['velX_s'] - sats['velX_c']
    sats['vpecy'] = sats['velY_s'] - sats['velY_c']
    sats['vpecz'] = sats['velZ_s'] - sats['velZ_c']
    sats['vpec'] = np.sqrt(sats['vpecx']**2+sats['vpecy']**2+sats['vpecz']**2) #magn of peculiar velocity

    sats = sats['fofCentralId', 'galaxyId', 'mvir', 'stellarMass', 'vpecx', 'vpecy', 'vpecz', 'vpec', 'distanceToCentralGalX', 'distanceToCentralGalY', 'distanceToCentralGalZ']
    sats.rename_column('distanceToCentralGalX', 'rx')
    sats.rename_column('distanceToCentralGalY', 'ry')
    sats.rename_column('distanceToCentralGalZ', 'rz')
    sats['r'] = np.sqrt(sats['rx']**2+sats['ry']**2+sats['rz']**2)

    print('Writing tables ' + str(sn) + ' ...')
    sats.write('data/sats' + str(sn) + '.fits', format ='fits', overwrite=True)


def get_nums(sn):
    print(sn)
    print('Loading...')
    s = Table.read('data/cube' + str(sn) + '>0.fits')

    print('Grouping...')
    g = s.group_by('fofCentralId').groups
    s = Table(names=('fofCentralId', 'num_sat'),
              data=(
              g.keys['fofCentralId'].data, g.indices[1:] - g.indices[:-1]))

    print('Loading...')
    c = Table.read('data/cube' + str(sn) + '=0.fits')
    c.rename_column('galaxyId', 'fofCentralId')

    print('Joining...')
    s = join(s, c, 'fofCentralId')['fofCentralId', 'stellarMass', 'mvir', 'num_sat']

    print('Writing...')
    s.write('data/nums' + str(sn) + '.fits', overwrite = True)


def joinsats(sn, clean=False, outname='data/satsnum{}.fits', **kwargs):
    print(sn)
    print('Loading sats')
    sats = Table.read('data/sats{}.fits'.format(sn))
    print('Loading nums')
    nums = Table.read('data/nums{}.fits'.format(sn))['fofCentralId', 'num_sat']

    if clean:
        nums = nums[np.where(nums['num_sat'] > 2)]

    print('Joining...')
    sats = join(nums, sats, 'fofCentralId')

    print('Writing...')
    sats.write(outname.format(sn), **kwargs)


def satkin_combine(sn):
    centrals = load_table(FILES['cube'](sn, sats=False))
    sats = load_table(FILES['cube'](sn, sats=True))

    print('Joining ...')
    centrals.rename_column('galaxyId', 'fofCentralId')
    sats = join(sats, centrals, 'fofCentralId',
                table_names=['s', 'c'])  # column names like name_s/name_c

    print('Calculating velocities and distances ' + str(sn) + ' ...')
    sats['vpecx'] = sats['velX_s'] - sats['velX_c']
    sats['vpecy'] = sats['velY_s'] - sats['velY_c']
    sats['vpecz'] = sats['velZ_s'] - sats['velZ_c']
    sats['vpec'] = np.sqrt(sats['vpecx'] ** 2
                           + sats['vpecy'] ** 2
                           + sats['vpecz'] ** 2)  # magn of peculiar velocity

    sats = sats['fofCentralId', 'galaxyId',
                'mvir', 'stellarMass',
                'vpecx', 'vpecy', 'vpecz', 'vpec',
                'distanceToCentralGalX', 'distanceToCentralGalY', 'distanceToCentralGalZ']
    sats.rename_column('distanceToCentralGalX', 'rx')
    sats.rename_column('distanceToCentralGalY', 'ry')
    sats.rename_column('distanceToCentralGalZ', 'rz')
    sats['r'] = np.sqrt(sats['rx'] ** 2
                        + sats['ry'] ** 2
                        + sats['rz'] ** 2)

    write_table(sats, FILES['sats'](sn))

    return centrals, sats

def satkin_get_nums(centrals, sats):
    g = sats.group_by('fofCentralId').groups
    nums = Table(names=('fofCentralId', 'num_sat'),
                 data=(g.keys['fofCentralId'].data, g.indices[1:] - g.indices[:-1]))

    nums = join(nums, centrals, 'fofCentralId')['fofCentralId', 'stellarMass', 'mvir', 'num_sat']

    write_table(nums, FILES['nums'](sn))

    sats = join(sats, nums['fofCentralId', 'num_sat'], 'fofCentralId')

    write_table(sats, FILES['sats'](sn))

    return nums, sats

def satkin_sim(sn):
    centrals, sats = satkin_combine(sn)
    nums, sats = satkin_get_nums(centrals, sats)

    del centrals, nums

    swp = SWProcedure(sn, False)
    swp.sats = sats
    swp.go()


class Procedure:
    def __init__(self, sn, observe, stdev_file, reg_file):
        self.sn = sn
        self.observe = observe

        self.sats = None
        self.stdev = None
        self.regression = None
        self.predictions = None

        self.stdev_file = stdev_file
        self.reg_file = reg_file


    def load_sats(self):
        if self.sats is None:
            fname = FILES['sats'](self.sn)
            print('Loading table', fname)
            t = Table.read(fname)
            t = t[np.where(t['stellarMass'] > 1e-4)]
            t['mvir'] = 10 + np.log10(t['mvir'])
            t['stellarMass'] = 10 + np.log10(t['stellarMass'])
            self.sats = t

    def load_stdev(self):
        if self.stdev is None:
            print('Loading table', self.stdev_file)
            self.stdev = Table.read(self.stdev_file)
    def write_stdev(self):
        print('Writing to', self.stdev_file)
        self.stdev.write(self.stdev_file, overwrite=True)

    def get_subset(self):
        raise NotImplementedError

    def regress(self):
        raise NotImplementedError
    def go(self):
        raise NotImplementedError

    def predict(self):
        p = self.get_regression()
        self.load_stdev()

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
            with open(self.reg_file) as f:
                self.regression = json.load(f)[str(self.sn)]
        return self.regression


class SWProcedure(Procedure):
    def __init__(self, sn, observe):
        Procedure.__init__(self, sn, observe,
                           reg_file=FILES['sw-reg'],
                           stdev_file=FILES['sw-sigmas'](sn, 'stellarMass', observe))

        self.mvir = None
        self.mvir_file = FILES['sw-sigmas'](sn, 'mvir', observe)

    def load_mvir(self):
        if self.mvir is None:
            print('Loading table', self.mvir_file)
            self.mvir = Table.read(self.mvir_file)
    def write_mvir(self):
        fname = FILES['sw-sigmas'](self.sn, 'mvir', self.observe)
        print('Writing to', fname)
        self.mvir.write(fname, overwrite=True)

    def get_subset(self):
        t = Table.read(FILES['nums'](self.sn))
        return t[np.where(t['num_sat'] > 2)]

    def store_stdev(self, stdev):
        self.stdev = stdev
    def store_mvir(self, mvir):
        self.mvir = mvir

    def bin(self, left, width, mname, store, write=None):
        print('Binning by', mname)
        self.sats['bin'] = np.floor_divide(self.sats[mname] - left, width)
        binned = self.sats.group_by('bin')

        print('Aggregating')
        res = binned['bin', 'vpec', 'vpecx', 'vpecy', 'vpecz'].groups.aggregate(np.std)

        res.rename_column('vpec', 'sigma')
        res.rename_column('vpecx', 'sigmax')
        res.rename_column('vpecy', 'sigmay')
        res.rename_column('vpecz', 'sigmaz')

        res['mean'] = binned['vpec'].groups.aggregate(np.mean)
        res['median'] = binned['vpec'].groups.aggregate(np.median)
        res[mname] = left + (res['bin'] + 0.5) * width
        res['N'] = binned['fofCentralId'].groups.aggregate(lambda x: len(np.unique(x)))

        store(res)

        if callable(write):
            write()

    def plot(self, mname, xlim, xlabel='$\log_{10}(M / M_\odot)$', sparse=333):
        print('Plotting...')
        fig, ax = plt.subplots(2, 1, sharex='col', squeeze=True,
                               gridspec_kw={'height_ratios': [3, 1],
                                            'hspace': 0})

        plt.sca(ax[0])
        r = np.random.randint(0, len(self.sats), int(len(self.sats) / sparse))
        plt.plot(self.sats[mname][r], self.sats['vpec'][r], ',', label='_nolegend_')

        X = self.stdev[mname]
        Y = self.stdev
        plt.plot(X, Y['sigma'], '-r', linewidth=3, label='$\sigma_{3D}$')
        plt.plot(X, Y['sigmax'], '-y', linewidth=2, label='$\sigma_{los}$')
        plt.plot(X, Y['sigmay'], '-y', linewidth=2, label='_nolegend_')
        plt.plot(X, Y['sigmaz'], '-y', linewidth=2, label='_nolegend_')
        plt.plot(X, Y['mean'], '-g', linewidth=3, label='$\langle |v| \\rangle$')
        plt.semilogy()

        plt.legend(loc='lower right')

        plt.ylabel('Speed, km/s')


        print('Plotting...')
        plt.sca(ax[1])
        plt.plot(X, Y['mean'] / Y['sigmax'])
        plt.plot(xlim, [np.sqrt(3)]*2, ':')
        plt.plot(xlim, [np.sqrt(8 / np.pi)]*2, ':')

        plt.xlabel(xlabel)
        plt.ylabel('$v_{3D}$ / $\sigma_{los}$')

        plt.ylim(1, 2)
        plt.xlim(xlim)

    def regress(self, write=True):
        self.load_mvir()

        print('Regressing...')
        indices, fit = fit_or(models.PowerLawModel(), self.mvir['sigmax'], x=np.power(10, self.mvir['mvir']))

        self.regression = {
            'exp': fit.params['exponent'].value,
            'err_exp': fit.params['exponent'].stderr,
            'amp': fit.params['amplitude'].value,
            'err_amp': fit.params['amplitude'].stderr,
            'indices': indices.tolist()
        }

        if write:
            self.write_regression()

    def go(self, write=True, plot=False):
        self.load_sats()

        if self.observe:
            self.sats = Observer.get().observe(self.sats)

        self.bin(PLOTDATA['bins']['stellarMass'][0], BINWIDTH, 'stellarMass', self.store_stdev, self.write_stdev)
        self.bin(PLOTDATA['bins']['mvir'][0], BINWIDTH, 'mvir', self.store_mvir, self.write_mvir)

        if self.observe:
            self.regress(write)

        if plot:
            for mname in ('stellarMass', 'mvir'):
                xlim = PLOTDATA['bins'][mname][0], PLOTDATA['bins'][mname][-1]

                self.plot(mname, xlim)
                fname = FILES['plot-sw-sigmas'](self.sn, mname, self.observe)
                plt.savefig(fname)
                print('Figure saved:', fname)
                plt.close()


class HWProcedure(Procedure):
    def __init__(self, sn, observe, stdev_file=None):
        Procedure.__init__(self, sn, observe,
                           reg_file=FILES['hw-reg'],
                           stdev_file=stdev_file)

        self.rms = None
        self.rms_file = FILES['hw-rms'](sn, observe)

    def bin(self, left, width, write=True):
        raise NotImplementedError

    def load_rms(self):
        if self.rms is None:
            print('Loading table', self.rms_file)
            self.rms = Table.read(self.rms_file)
    def write_rms(self):
        print('Writing to', self.rms_file)
        self.rms.write(self.rms_file, overwrite=True)

    def get_subset(self):
        self.load_rms()
        return self.rms

    def calculate(self, write=True):
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
        nums = Table.read(FILES['nums'](self.sn))
        rms = join(rms, nums, keys='fofCentralId')
        self.rms = rms

        if write:
            self.write_rms()

    def regress(self, write=True):
        self.load_rms()

        t = deal(self.rms)
        print('Regressing...')
        indices, fit = fit_or(models.PowerLawModel(), t['rmsx'], x=np.power(10, t['mv']))

        self.regression = {
            'exp': fit.params['exponent'].value,
            'err_exp': fit.params['exponent'].stderr,
            'amp': fit.params['amplitude'].value,
            'err_amp': fit.params['amplitude'].stderr
        }

        if write:
            self.write_regression()

    def go(self, write=True):
        self.load_sats()

        if self.observe:
            self.sats = Observer.get().observe(self.sats)

        self.calculate(write)
        try:
            self.bin(PLOTDATA['bins']['stellarMass'][0], BINWIDTH, write)
        except NotImplementedError:
            pass

        if not self.observe:
            self.regress(write)

class HWSigmaProcedure(HWProcedure):
    @staticmethod
    def collapse(binned, cols):
        raise NotImplementedError

    def bin(self, left, width, write=True):
        self.load_rms()

        print('Binning in stellarMass')
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
    def __init__(self, sn, observe):
        HWSigmaProcedure.__init__(self, sn, observe,
                                  stdev_file=FILES['hwp-sigmas'](sn, observe))

    @staticmethod
    def collapse(binned, cols):
        for col in cols:
            binned[col] = binned[col]**2

        res = binned[['bin'] + cols].groups.aggregate(np.mean)

        for col in cols:
            res[col] = np.sqrt(res[col])

        return res

class HWAProcedure(HWSigmaProcedure):
    def __init__(self, sn, observe):
        HWSigmaProcedure.__init__(self, sn, observe,
                                  stdev_file=FILES['hwa-sigmas'](sn, observe))

    @staticmethod
    def collapse(binned, cols):
        return binned[['bin'] + cols].groups.aggregate(np.mean)

class HWMProcedure(HWSigmaProcedure):
    def __init__(self, sn, observe):
        HWSigmaProcedure.__init__(self, sn, observe,
                                  stdev_file=FILES['hwm-sigmas'](sn, observe))

    @staticmethod
    def collapse(binned, cols):
        return binned[['bin'] + cols].groups.aggregate(np.median)

class HWNProcedure(HWProcedure):
    def bin(self, left, width, write=True):
        raise NotImplementedError

    def predict(self):
        p = self.get_regression()
        self.load_rms()

        t = deal(self.rms['stellarMass', 'mvir', 'rmsx'])
        t['mv'] = np.log10(t['rmsx'] / p['amp']) / p['exp']

        b = binX(t, 'ms', 'mv', bins=PLOTDATA['bins']['stellarMass'],
                 funcsx={'mean': np.mean},
                 funcsy={
                     'mean': np.mean, 'median': np.median, 'std': np.std,
                     'p16': percentile_wrapper(16),
                     'p84': percentile_wrapper(84)})

        self.predictions = {
            'median': Data(mstar=b[0]['mean'], mvir=b[1]['median'],
                           p16=b[1]['p16'], p84=b[1]['p84']),
            'mean': Data(mstar=b[0]['mean'], mvir=b[1]['mean'], err_mvir=b[1]['std'])
        }