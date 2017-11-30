import os
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
                       isline=lambda line: not (
                       line[0] == '#' or line[0] == 'g'))
        d.go()

        csv2fits(csvname, fname, cols)

        os.remove(csvname)

        def download_centrals(self):
            self._download(FILES['cube'](self.sn, _sats=False), '=0', cols_centrals)
        def download_satellites(self):
            self._download(FILES['cube'](self.sn, _sats=True), '>0', cols_satellites)

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
