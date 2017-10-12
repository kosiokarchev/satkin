import subprocess


output_path = 'test.csv'

username = '...'
password = '...'
url = 'http://gavo.mpa-garching.mpg.de/MyMillennium/'

limit = 1000000
columns = ['galaxyId', 'type',
           'phKey', 'x', 'y', 'z', 'velX', 'velY', 'velZ',
           'redshift', 'lookBackTime',
           'rvir', 'bulgeSize', 'stellarDiskRadius',
           'mvir', 'coldGas', 'stellarMass', 'bulgeMass', 'diskMass', 'hotGas', 'ejectedMass', 'blackHoleMass', 'icmStellarMass',
           'sfr', 'sfrBulge',
           'centralMVir', 'centralRvir', 'distanceToCentralGalX', 'distanceToCentralGalY', 'distanceToCentralGalZ',
           'treeId', 'descendantId', 'mainLeafId', 'treeRootId', 'firstProgenitorId', 'nextProgenitorId', 'lastProgenitorId', 'haloId', 'subHaloId', 'fofCentralId', 'fofSubhaloId']
table = 'Henriques2015a..MRscPlanck1'
sql = 'SELECT TOP {limit} {cols} FROM {table}'

SQL = sql.format(limit=limit, cols=','.join(columns), table=table)
full_url = url + '?action=doQuery&SQL=' + SQL

wget = ('wget', '--http-user='+username, '--http-passwd='+password,'-O', output_path, full_url)
subprocess.run(wget, check=True)
