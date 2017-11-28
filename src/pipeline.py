import os
from lib import *

TABLES = {'H15-cube': 'Henriques2015a..MRscPlanck1',
          'H15-cone': lambda n=1: 'Henriques2015a.cones.MRscPlanck1_BC03_{:0>3}'.format(n)}
cols_central = [
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

    def download(self):
        for typewhere, cols, name in zip(('=0', '>0'),
                                         (cols_central, cols_satellites),
                                         ('centrals', 'satellites')):
            print('Downloading', name)
            where = 'snapnum={} AND type{}'.format(self.sn, typewhere)
            q = Query(cols, TABLES['H15-cube'], where)

            fname = FILES['cube'](self.sn).replace('.fits', '.csv')

            d = Downloader(q, fname,
                           isline=lambda line: not (line[0] == '#' or line[0] == 'g'))
            d.go()

            csv2fits(fname, FILES['cube'](self.sn), cols)

            os.remove(fname)