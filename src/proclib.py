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

    def restore(self):
        raise NotImplementedError
    def load(self):
        raise NotImplementedError
    def plot(self, *args, **kwargs):
        raise NotImplementedError
    def predict(self):
        raise NotImplementedError
    def go(self):
        raise NotImplementedError



class SWProcedure(Procedure):
    def __init__(self, sn, observe):
        Procedure.__init__(self, sn, observe)
        self.sats = None

    def restore(self, mname='stellarMass'):
        self.stdev = Table.read(FILES['sigmas'](sn, mname, self.observe))

    def load(self):
        print('Loading table...')
        t = Table.read(FILES['sats'](self.sn))
        t['mv'] = 10 + np.log10(t['mvir'])
        t['ms'] = 10 + np.log10(t['stellarMass'])
        t = t[np.where(t['stellarMass'] > 6)]
        self.sats = t

    def bin(self, mname, left, width):
        self.sats['bin'] = np.floor_divide(self.sats[mname] - left, width)
        binned = self.sats.group_by('bin')

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
        r = np.random.randint(int(len(sats) / sparse))
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
        plt.plot(X, Y / X)
        plt.plot(xlim, [np.sqrt(3)]*2, ':')
        plt.plot(xlim, [np.sqrt(8 / np.pi)]*2, ':')

        plt.xlabel(xlabel)
        plt.ylabel('$v_{3D}$ / $\sigma_{los}$')

        plt.ylim(1, 2)
        plt.xlim(xlim)

    def go(self, mnames=('mvir', 'stellarMass'), write=True):
        self.load()

        if self.observe:
            self.sats = Observer.get().observe(self.sats)

        for mname in mnames:
            xlim = PLOTDATA['bins'][mname][0], PLOTDATA['bins'][mname][-1]
            self.bin(mname, xlim[0], BINWIDTH)

            if write:
                self.stdev.write(FILES['sigmas'](self.sn, mname, self.observe), overwrite=True)

            self.plot(mname, xlim)
            fname = FILES['plotsigmas'](sn)
            plt.savefig(fname)
            print('Figure saved:', fname)
            plt.close()

    def predict(self):
        params = Regressor.get()[self.sn]
        sn, amp, err_amp, exp, err_exp = params

        data = Table.read(FILES['sigmas'](self.sn, 'stellarMass', self.observe))

        data['mvir'] = np.power((data['sigmax'] / amp), 1 / exp)

        data['relerr_mvir'] = np.sqrt(
            ((1 / exp) * (err_amp / amp)) ** 2 + (
            np.log(data['mvir']) * (err_exp / exp)) ** 2
        )

        return Data(mvir=data['mvir'], mstar=data['stellarMass'],
                    relerr_mvir=data['relerr_mvir'],
                    N=data['N'])