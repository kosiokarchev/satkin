import os
from multiprocessing import Process, cpu_count, current_process
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
        centrals = t['galaxyId', 'stellarMass', 'z_app'][t['iscen']]
        centrals.rename_column('galaxyId', 'fofCentralId')

        sats = t['cId', 'z_app', 'd'][t['issat']]
        sats.rename_column('cId', 'fofCentralId')

        sats = join(sats, centrals, 'fofCentralId')

        sats['dz'] = sats['z_app_1'] - sats['z_app_2']
        sats.rename_column('z_app_2', 'z')

        return centrals, sats['fofCentralId', 'stellarMass', 'z', 'dz', 'd']

    @staticmethod
    def fit_cumgauss(dv):
        dv = np.abs(dv)
        dv.sort()

        p0 = (300, 500 / len(dv), 200 / len(dv) / np.max(dv))

        y = (np.arange(len(dv)) + 1) / len(dv)
        func = lambda args: ((y - cumgauss(dv, *args)) ** 2).sum()
        res = minimize(func, p0, method='Nelder-Mead')

        if res.success:
            popt = res.x
            pcov = np.linalg.pinv(nd.Hessian(func)(popt))
            errs = np.sqrt(np.abs(np.diag(pcov)))
            return popt, errs
        else:
            return np.array(p0), np.array([np.inf]*3)

    s2m = HWAProcedure(34, False).get_regression()

    @staticmethod
    def sigma2mvir(sigma, sigma_err):
        mvir = np.log10(sigma / ConePipeline.s2m['amp']) / ConePipeline.s2m['exp']
        mvir_err = (1 / (np.log(10) * ConePipeline.s2m['exp'])) * sigma_err / sigma

        return mvir, mvir_err

    def __init__(self, cone, nvircen=2, nvirsat=1, dvcen=3000, dvsat=3000, binwidth=0.1, oname='res{}.fits', quiet=False, seed=None):
        self.cone = cone

        self.nvircen = nvircen
        self.nvirsat = nvirsat
        self.dvcen = dvcen
        self.dvsat = dvsat

        self.binwidth = binwidth

        self.oname = oname

        self.quiet = quiet
        self.seed = seed

        self.sample = None
        self.cen = None
        self.sats = None
        self.res = None

    def create_sample(self):
        np.random.seed(self.seed)
        s = self.cone['p'] > np.random.random(len(self.cone))
        self.sample = self.cone['galaxyId', 'fofCentralId',
                                'mvir', 'stellarMass', 'sfr',
                                'z_app', 'ra', 'dec', 'd_comoving'][s]

    def examine(self):
        res = galocate(self.sample,
                       nvircen=self.nvircen, nvirsat=self.nvirsat,
                       dvcen=self.dvcen, dvsat=self.dvsat)

        self.sample['iscen'] = res['iscen']
        self.sample['cId']   = res['fofCentralId']
        self.sample['d']     = res['d']
        self.sample['issat'] = np.logical_and(~res['iscen'],
                                              res['fofCentralId'])

        self.cen, self.sats = ConePipeline.split(self.sample)

    def predict(self):
        g = Table()
        g['dv'] = 3e5 * self.sats['dz'] / (1 + self.sats['z'])
        g['ms'] = 10 + np.log10(self.sats['stellarMass'])
        g['bin'] = g['ms'] // self.binwidth

        g = g.group_by('bin')

        ms = []
        sigma = []
        sigma_err = []
        for i in range(len(g.groups.indices)-1):
            start = g.groups.indices[i]
            end = g.groups.indices[i+1]

            if end-start < 5:
                continue

            ms.append((g.groups.keys['bin'][i] + 0.5) * self.binwidth)
            if not self.quiet:
                print('ms={:.2f} [{}:{}] / {}'.format(ms[-1], start, end, len(g)))
            dv = g['dv'][start:end]
            res = ConePipeline.fit_cumgauss(dv)
            sigma.append(res[0][0])
            sigma_err.append(res[1][0])

        self.res = Table(dict(ms=ms, sigma=sigma, sigma_err=sigma_err))

        self.res['mv'], self.res['mv_err'] = ConePipeline.sigma2mvir(self.res['sigma'], self.res['sigma_err'])

    def plot(self):
        plt.errorbar(self.res['mv'], self.res['ms'], xerr=self.res['mv_err'], linestyle='', marker='o')
        plt.xlim(10, 15)
        plt.ylim(7, 12)

    def go(self):
        self.create_sample()
        self.examine()
        self.predict()

    @staticmethod
    def bootstrap(fname, oname, N, start=0, nproc=cpu_count()):
        def bootstrap_one(i, n, seed):
            cone = load(fname)
            ConePipeline(cone, oname=oname, quiet=True, seed=seed)._bootstrap(i, n)
            print(current_process().name, 'has finished')

        n = np.ceil(N/nproc)
        procs = []
        for i in np.arange(start, start+N, n):
            proc = Process(target=bootstrap_one,
                           args=(int(i), int(n),
                                 np.uint32(np.random.random()*(2**31))))
            proc.start()
            procs.append(proc)

        for proc in procs:
            proc.join()

    def _bootstrap(self, start, n):
        for i in range(start, start+n):
            print('Starting bootstrap', i)
            self.go()
            write(self.res, self.oname.format(i))
