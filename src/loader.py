import sqlite3
import numpy as np

FNAMES = ['data/s34.csv', 'data/s38.csv', 'data/s45.csv']
# FNAMES = ['/users/kosio/s45s.csv']

db = sqlite3.connect('data/data.db')

with db:
    db.execute('''CREATE TABLE IF NOT EXISTS cube (
    galaxyId INT PRIMARY KEY,
    haloId INT, subHaloId INT, fofCentralId INT, fofSubhaloId INT,
    centralMVir REAL, centralRvir REAL, distanceToCentralGalX REAL, distanceToCentralGalY REAL, distanceToCentralGalZ REAL,
    type INT,
    x REAL, y REAL, z REAL, velX REAL, velY REAL, velZ REAL, stellarSpinX REAL, stellarSpinY REAL, stellarSpinZ REAL,
    mvir REAL, rvir REAL, stellarMass REAL, sfr REAL
    ) WITHOUT ROWID''')

for fname in FNAMES:
    print('Reading file '+fname)
    rows = np.loadtxt(fname, skiprows=1, delimiter=',')
    print('Found {} rows.'.format(len(rows)))
    with db:
        db.executemany('INSERT INTO cube VALUES ('+(','.join(['?']*24))+')', rows)
    print('File added to database.')