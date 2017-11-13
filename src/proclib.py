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


class Procedure:
    def __init__(self, sn, observe):
        self.sn = sn
        self.observe = observe
        self.stdev = None
        self.regression = None
        self.reg_file = None
        self.predictions = None

    def restore(self):
        raise NotImplementedError
    def plot(self, *args, **kwargs):
        raise NotImplementedError
    def regress(self):
        raise NotImplementedError
    def predict(self):
        raise NotImplementedError
    def go(self):
        raise NotImplementedError

    def get_predictions(self):
        if self.predictions is None:
            self.predict()
        return self.predictions

    def write_regression(self):
        with open(self.reg_file, 'w') as f:
            f.write(json.dumps(self.regression))
    def get_regression(self):
        if self.regression is None:
            with open(self.reg_file) as f:
                self.regression = json.load(f)
        return self.regression


class SWProcedure(Procedure):
    def __init__(self, sn, observe):
        Procedure.__init__(self, sn, observe)
        self.sats = None

    def restore(self, mname='stellarMass'):
        if self.stdev is None:
            self.stdev = Table.read(FILES['sw-sigmas'](self.sn, mname, self.observe))

    def load_sats(self):
        if self.sats is None:
            fname = FILES['sats'](self.sn)
            print('Loading table', fname)
            t = Table.read(fname)
            t['mvir'] = 10 + np.log10(t['mvir'])
            t['stellarMass'] = 10 + np.log10(t['stellarMass'])
            t = t[np.where(t['stellarMass'] > 6)]
            self.sats = t

    def write_stdev(self, mname):
        fname = FILES['sw-sigmas'](self.sn, mname, self.observe)
        print('Writing to', fname)
        self.stdev.write(fname, overwrite=True)

    def bin(self, mname, left, width):
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

        self.stdev = res

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

    def go(self, mnames=('mvir', 'stellarMass'), write=True):
        self.load_sats()

        if self.observe:
            self.sats = Observer.get().observe(self.sats)

        for mname in mnames:
            xlim = PLOTDATA['bins'][mname][0], PLOTDATA['bins'][mname][-1]
            self.bin(mname, xlim[0], BINWIDTH)

            if write:
                self.write_stdev(mname)

            self.plot(mname, xlim)
            fname = FILES['plot-sw-sigmas'](self.sn, mname, self.observe)
            plt.savefig(fname)
            print('Figure saved:', fname)
            plt.close()

    def regress(self, write=True):
        self.restore('mvir')

        t = self.stdev
        indices, fit = fit_or(models.PowerLawModel(), t['sigmax'], x=np.power(10, t['mvir']))

        if write:
            t.write(FILES['sw-reg-mvir'], format='ascii.ecsv', overwrite=True)

        self.regression = {
            'exp': fit.params['exponent'].value,
            'err_exp': fit.params['exponent'].stderr,
            'amp': fit.params['amplitude'].value,
            'err_amp': fit.params['amplitude'].value,
            'indices': list(indices)
        }

        if write:
            self.write_regression()

    def predict(self):
        params = Regressor.get()[self.sn]
        sn, amp, err_amp, exp, err_exp = params

        data = Table.read(FILES['sw-sigmas'](self.sn, 'stellarMass', self.observe))

        data['mvir'] = np.log10(data['sigmax'] / amp) / exp
        # data['mvir'] = np.log10(np.power((data['sigmax'] / amp), 1 / exp))

        data['relerr_mvir'] = np.sqrt(
            ((1 / exp) * (err_amp / amp)) ** 2 + (
            np.log(data['mvir']) * (err_exp / exp)) ** 2
        )

        self.predictions = {
            'mean': Data(mvir=data['mvir'], mstar=data['stellarMass'],
                         err_mvir=data['relerr_mvir'] / np.log(10),
                         N=data['N'])
        }


class HWProcedure(Procedure):
    def __init__(self, sn, observe):
        Procedure.__init__(self, sn, observe)
        self.sats = None
        self.rms = None
        self.predicted = None

        self.reg_file = FILES['hw-reg-mvir'](sn)

    def restore(self):
        try:
            if self.stdev is None:
                self.stdev = Table.read(FILES['hw-sigmas'](self.sn, self.observe))
        except:
            pass

    def load_sats(self):
        if self.sats is None:
            fname = FILES['sats'](self.sn)
            print('Loading table', fname)
            t = Table.read(fname)
            t = t[np.where(t['stellarMass'] > 1e-4)]
            t['mvir'] = 10 + np.log10(t['mvir'])
            t['stellarMass'] = 10 + np.log10(t['stellarMass'])
            self.sats = t
    def load_rms(self):
        if self.rms is None:
            fname = FILES['hw-rms'](self.sn, self.observe)
            print('Loading table', fname)
            self.rms = Table.read(fname)

    def get_subset(self):
        self.load_rms()
        return self.rms

    def calculate(self):
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

    def bin(self, left, width):
        self.load_rms()

        print('Binning in stellarMass')
        self.rms['bin'] = np.floor_divide(self.rms['stellarMass'] - left, width)
        binned = self.rms.group_by('bin')

        cols = ['rms', 'rmsx', 'rmsy', 'rmsz']
        for col in cols:
            binned[col] = binned[col]**2

        print('Calculating mean variance...')
        res = binned[['bin'] + cols].groups.aggregate(np.mean)
        res['stellarMass'] = left + (res['bin'] + 0.5) * width
        res['N'] = binned.groups.indices[1:] - binned.groups.indices[:-1]

        self.stdev = res

    def regress(self, write=True):
        self.load_rms()

        t = deal(self.rms)
        print('Regressing...')
        indices, fit = fit_or(models.PowerLawModel(), t['rmsx'], x=np.power(10, t['mv']))

        self.regression = {
            'exp': fit.params['exponent'].value,
            'err_exp': fit.params['exponent'].stderr,
            'amp': fit.params['amplitude'].value,
            'err_amp': fit.params['amplitude'].value
        }

        if write:
            self.write_regression()

    def go(self, write=True):
        self.load_sats()

        if self.observe:
            self.sats = Observer.get().observe(self.sats)

        self.calculate()

        if write:
            fname = FILES['hw-rms'](self.sn, self.observe)
            print('Writing to', fname)
            self.rms.write(fname, overwrite=True)

        self.bin(PLOTDATA['bins']['stellarMass'][0], BINWIDTH)

        if write:
            fname = FILES['hwp-sigmas'](self.sn, self.observe)
            print('Writing to', fname)
            self.stdev.write(fname, overwrite=True)

        self.regress(write)

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
        self.predicted = b

        self.predictions = {
            'median': Data(mstar=b[0]['mean'], mvir=b[1]['median'],
                           p16=b[1]['p16'], p84=b[1]['p84']),
            'mean': Data(mstar=b[0]['mean'], mvir=b[1]['mean'], err_mvir=b[1]['std'])
        }