import os
from glob import glob
from multiprocessing import Pool, cpu_count
import numdifftools as nd
from scipy.optimize import minimize
from lib import *

TABLES = {'H15-cube': 'Henriques2015a..MRscPlanck1',
          'H15-cone': lambda n=1: 'Henriques2015a.cones.MRscPlanck1_BC03_{:0>3}'.format(n)}
cols_centrals = [
    'galaxyId',
    'x', 'y', 'z', 'velX', 'velY', 'velZ', 'stellarSpinX', 'stellarSpinY',
    'stellarSpinZ',
    'mvir', 'rvir', 'stellarMass', 'sfr',
]
cols_satellites = [
    'galaxyId', 'fofCentralId',
    'velX', 'velY', 'velZ',
    'distanceToCentralGalX', 'distanceToCentralGalY', 'distanceToCentralGalZ'
]

class Pipeline:
    def __init__(self, sn):
        self.sn = sn

    def _download(self, fname, typewhere, cols):
        print('Downloading', fname)
        where = 'snapnum={} AND type{}'.format(self.sn, typewhere)
        q = Query(cols, TABLES['H15-cube'], where)

        csvname = fname.replace('.fits', '.csv')

        d = Downloader(q, csvname,
                       isline=lambda line: not (line[0] == '#' or line[0] == 'g'))
        d.go()

        csv2fits(csvname, fname, cols)

        os.remove(csvname)

        # def download_centrals(self):
        #     self._download(FILES['cube'](self.sn, _sats=False), '=0', cols_centrals)
        # def download_satellites(self):
        #     self._download(FILES['cube'](self.sn, _sats=True), '>0', cols_satellites)

    def download(self, sats='both'):
        if (not sats) or (sats=='both'):
            self._download(FILES['cube'](self.sn, _sats=False), '=0', cols_centrals)
        if sats:
            self._download(FILES['cube'](self.sn, _sats=True), '>0', cols_satellites)

    def combine(self):
        print('Calculating satellite characteristics')
        centrals = load(FILES['cube'](self.sn, _sats=False))
        sats = load(FILES['cube'](self.sn, _sats=True))

        print('Joining')
        centrals.rename_column('galaxyId', 'fofCentralId')
        sats = join(sats, centrals, 'fofCentralId',
                    table_names=['s', 'c'])  # column names like name_s/name_c

        print('Calculating velocities and distances.')
        sats['vpecx'] = sats['velX_s'] - sats['velX_c']
        sats['vpecy'] = sats['velY_s'] - sats['velY_c']
        sats['vpecz'] = sats['velZ_s'] - sats['velZ_c']
        sats['vpec'] = np.sqrt(sats['vpecx'] ** 2
                               + sats['vpecy'] ** 2
                               + sats[
                                   'vpecz'] ** 2)  # magn of peculiar velocity

        sats = sats['fofCentralId', 'galaxyId',
                    'mvir', 'stellarMass',
                    'vpecx', 'vpecy', 'vpecz', 'vpec',
                    'distanceToCentralGalX', 'distanceToCentralGalY', 'distanceToCentralGalZ']
        sats.rename_column('distanceToCentralGalX', 'rx')
        sats.rename_column('distanceToCentralGalY', 'ry')
        sats.rename_column('distanceToCentralGalZ', 'rz')
        sats['rx'] = -sats['rx']
        sats['ry'] = -sats['ry']
        sats['rz'] = -sats['rz']
        sats['r'] = np.sqrt(sats['rx'] ** 2
                            + sats['ry'] ** 2
                            + sats['rz'] ** 2)

        write(sats, FILES['sats'](self.sn))

    def count(self):
        print('Counting satellites')
        centrals = load(FILES['cube'](self.sn, _sats=False))
        sats = load(FILES['sats'](self.sn))

        g = sats.group_by('fofCentralId').groups
        nums = Table(names=('fofCentralId', 'num_sat'),
                     data=(g.keys['fofCentralId'].data,
                           g.indices[1:] - g.indices[:-1]))

        nums = join(nums, centrals, 'fofCentralId')[
            'fofCentralId', 'stellarMass', 'mvir', 'num_sat']

        write(nums, FILES['nums'](self.sn))

        print('Joining')
        sats = join(sats, nums['fofCentralId', 'num_sat'], 'fofCentralId')

        write(sats, FILES['sats'](self.sn))

    def god_predict(self):
        print('Observing from God\'s viewpoint')
        nums = load(FILES['nums'](self.sn))

        t = deal(nums['mvir', 'stellarMass'])
        b = bin2d(t, 'mv', 'ms', )
        b.compute(PLOTDATA['bins']['mvir'], PLOTDATA['bins']['stellarMass'])
        b.save(FILES['3dbins'](self.sn))

    def manipulate(self):
        print('Manipulating the raw data')
        self.combine()
        self.count()

        unload(FILES['cube'](self.sn, _sats=False))

    def theplot(self):
        sats = load_sats(self.sn)

        for observe in (False, True):
            if observe:
                sats = Observer.get().observe(sats)

            for P, name in zip(procs, procnames):
                p = P(self.sn, observe)
                p.sats = sats
                d = p.get_predictions()
                if isinstance(d, Data):
                    d.save(FILES['data'](self.sn, name, observe))
                else:
                    for key in d:
                        d[key].save(FILES['data'](self.sn, name+'-'+key, observe))

    def __call__(self):
        print('Starting the pipeline for sn ', self.sn)
        self.download()
        self.manipulate()
        self.god_predict()
        self.theplot()
        print('Pipeline finished for sn ', self.sn)


class ConePipeline:
    @staticmethod
    def split(t):
        if 'fofCentralId' in t.colnames:
            iscen = t['galaxyId']
        centrals = t['galaxyId', 'stellarMass', 'z_app'][t['iscen']]
        centrals.rename_column('galaxyId', 'fofCentralId')

        sats = t['fofCentralId', 'z_app', 'd'][t['issat']]

        sats = join(sats, centrals, 'fofCentralId')

        sats['dz'] = sats['z_app_1'] - sats['z_app_2']
        sats.rename_column('z_app_2', 'z')

        return centrals, sats['fofCentralId', 'stellarMass', 'z', 'dz', 'd']

    @staticmethod
    def fit_cumgauss(dv, w=None, dodouble=True, interlopers=True, **kwargs):
        if w is None:
            w = np.full_like(dv, 1)
        dv = np.abs(dv)
        asort = np.argsort(dv)
        dv = dv[asort]
        w = w[asort]

        y = np.cumsum(w) / w.sum()
        x = dv / np.max(dv)

        mod = Cum2GaussConstModel(nan_policy='omit')
        params = mod.make_params(
            sigma1=0.1, A1=0,
            sigma2=dv[np.searchsorted(y, 0.68)], A2=0.9,
            B=0.1)
        params['sigma1'].vary = params['A1'].vary = False

        if not interlopers:
            params['B'].value = 0
            params['B'].vary = False

        try:
            res1 = mod.fit(y, x=dv, params=params)

            if dodouble:
                params2 = params.copy()
                params2['sigma1'].value = params2['sigma2'].value / 2
                params2['A1'].value = params2['A2'].value = params2[
                                                                'A2'].value / 2
                params2['sigma1'].vary = params2['A1'].vary = True

                try:
                    res2 = mod.fit(y, x=dv, params=params2, **kwargs)
                    if (res2.params['sigma'].stderr and
                            res2.params['sigma'].stderr < res2.params[
                                'sigma'].value):
                        return res2.params
                except:
                    pass
        except:
            return False

        return res1.params

    cols = ['galaxyId', 'stellarMass',
            'z_app', 'ra', 'dec', 'd_comoving']

    def __init__(self, cone,
                 nvircen=2, nvirsat=1, dvcen=3000, dvsat=3000, fM=1.0, d0=2.0,
                 binwidth=0.1, sampler='all', quiet=False):
        self.cone = cone

        self.nvircen = nvircen
        self.nvirsat = nvirsat
        self.dvcen = dvcen
        self.dvsat = dvsat
        self.fM = fM
        self.d0 = d0

        self.binwidth = binwidth

        if sampler in self.samplers:
            self.sampler = self.samplers[sampler]
        elif callable(sampler):
            self.sampler = sampler
        else:
            raise ValueError('sampler should be a callable or one of {}'.format(list(self.samplers.keys())))


        self.quiet = quiet

        self.sample = None
        self.cen = None
        self.sats = None
        self.res = None

    samplers = {
        'all': lambda cone: np.full(len(cone), True),
        'observe': lambda cone: cone['p'] > np.random.random(len(cone)),
        'bootstrap': lambda cone: (np.random.random(len(cone))*len(cone)).astype(int)
    }
    def create_sample(self):
        self.sample = self.cone[self.cols][self.sampler(self.cone)]

    def examine(self):
        res = galocate(self.sample,
                       nvircen=self.nvircen, nvirsat=self.nvirsat,
                       dvcen=self.dvcen, dvsat=self.dvsat,
                       fM=self.fM, d0=self.d0)

        self.sample['iscen'] = res['iscen']
        self.sample['fofCentralId'] = res['fofCentralId']
        self.sample['d']     = res['d']
        self.sample['issat'] = np.logical_and(~res['iscen'],
                                              res['fofCentralId'])

        self.cen, self.sats = self.split(self.sample)

        return self
    def set_examined(self, t):
        self.sample = t[self.cols + ['fofCentralId']]
        self.sample['iscen'] = t['galaxyId'] == t['fofCentralId']
        self.sample['issat'] = ~self.sample['iscen']
        self.sample['d'] = np.nan

        self.cen, self.sats = self.split(self.sample)
        return self

    def predict(self, **kwargs):
        g = Table()
        g['fofCentralId'] = self.sats['fofCentralId']
        g['dv'] = 3e5 * self.sats['dz'] / (1 + self.sats['z'])
        g['ms'] = 10 + np.log10(self.sats['stellarMass'])
        g['bin'] = g['ms'] // self.binwidth

        g = g.group_by('bin')

        ms = []
        sigma = []
        sigma_err = []
        N = []
        for i in range(len(g.groups.indices)-1):
            start = g.groups.indices[i]
            end = g.groups.indices[i+1]

            if end-start < 5:
                continue

            ms.append((g.groups.keys['bin'][i] + 0.5) * self.binwidth)
            if not self.quiet:
                print('ms={:.2f} [{}:{}] / {}'.format(ms[-1], start, end, len(g)))


            b = g['fofCentralId', 'dv'][start:end]

            norm = np.max(b['dv'])
            dv = b['dv'] / norm

            if not self.quiet:
                print('SW:')
            ressw = self.fit_cumgauss(dv, **kwargs)
            if ressw is False:
                ssw = shw = n = np.nan
                esw = np.inf
                ehw = np.inf
            else:
                ssw = ressw['sigma'].value
                esw = ressw['sigma'].stderr

                s1 = ressw['sigma1'].value
                a1 = ressw['A1'].value / (np.sqrt(2*np.pi) * s1)
                s2 = ressw['sigma2'].value
                a2 = ressw['A2'].value / (np.sqrt(2*np.pi) * s2)
                f = a1 * np.exp(-dv**2 / (2*s1**2)) + a2 * np.exp(-dv**2 / (2*s2**2))
                b['f'] = f / (0.5*ressw['B'].value + f)

                counts = b.group_by('fofCentralId').groups
                counts.keys['n'] = counts.indices[1:] - counts.indices[:-1]
                counts = counts.keys
                b = join(b, counts, 'fofCentralId')

                if not self.quiet:
                    print('HW:')
                reshw = self.fit_cumgauss(dv, b['f'] / b['n'], **kwargs)
                if reshw is False:
                    shw = np.nan
                    ehw = np.inf
                else:
                    shw = reshw['sigma'].value
                    ehw = reshw['sigma'].stderr


                counts = b['fofCentralId', 'f'].group_by('fofCentralId').groups.aggregate(np.nansum)
                counts.rename_column('f', 'n')
                n = np.nanmean(counts['n'])

            sigma.append((ssw*norm, shw*norm))
            sigma_err.append((esw*norm, ehw*norm))
            N.append(n)

        self.res = Table(dict(ms=ms, sigma=sigma, sigma_err=sigma_err, N=N))

        return self

    def plot(self):
        plt.errorbar(self.res['mv'], self.res['ms'], xerr=self.res['mv_err'], linestyle='', marker='o')
        plt.xlim(10, 15)
        plt.ylim(7, 12)

    def go(self):
        self.create_sample()
        self.examine()
        self.predict()
        return self


class ConeBootstrapper:
    cone = None
    oname = None
    pipelinekwargs = None

    @staticmethod
    def _init(fname, oname, pipelinekwargs):
        ConeBootstrapper.cone = load(fname)
        ConeBootstrapper.oname = oname
        ConeBootstrapper.pipelinekwargs = pipelinekwargs

    @staticmethod
    def _bootstrap_one(i, seed):
        np.random.seed(seed)
        print('Bootstrap', i)
        cp = ConePipeline(ConeBootstrapper.cone, quiet=True, **ConeBootstrapper.pipelinekwargs).go()
        write(cp.res, ConeBootstrapper.oname.format(i))

    @staticmethod
    def bootstrap(fname, oname, N, start=0, nproc=cpu_count(), **kwargs):
        with Pool(nproc,
                  initializer=ConeBootstrapper._init,
                  initargs=(fname, oname, kwargs)) as pool:
            pool.starmap(ConeBootstrapper._bootstrap_one,
                         zip(
                             range(start, start+N),
                             np.uint32(np.random.random(N)*2**32)
                         ))
            print('Bootstrap complete.')

    @staticmethod
    def pack(folder, fpattern='*', binwidth=0.1, rbins=(8, 12)):
        files = glob(os.path.join(folder, fpattern))

        ms = np.arange(rbins[0] + binwidth / 2, rbins[1], binwidth) // binwidth

        shape = (len(files), len(ms), 2)
        sigma = np.full(shape, np.nan)
        sigma_err = np.full(shape, np.inf)

        N = np.full(shape[:-1], np.nan)

        for i in range(len(files)):
            t = load(files[i])
            w = np.where(ms[:, np.newaxis] == (t['ms'].data // binwidth)[np.newaxis, :])

            sigma[i, w[0]] = t['sigma'][w[1]]
            sigma_err[i, w[0]] = t['sigma_err'][w[1]]

            N[i, w[0]] = t['N'][w[1]]

        return dict(sigma=sigma, sigma_err=sigma_err, N=N, ms=(ms+0.5)*binwidth)
