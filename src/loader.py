from __future__ import division, print_function

import sqlite3
import numpy as np
from astropy.table import Table


def reporter():
    print('.', end='', flush=True)


FNAME = 'data/cube.fits'

db = sqlite3.connect('data/data.db')
db.set_progress_handler(reporter, 10000)

with db:
    with open('src/sql/cube_create.sql') as f:
        query = f.read()
    db.execute(query)

print('Reading file '+FNAME)
t = Table.read(FNAME)
print('Found', len(t), 'rows.')

colnames = [col for col in t.columns]
ta = np.array([t[col].data for col in colnames]).transpose()
print('It is now an array with shape', ta.shape)

query = 'INSERT INTO cube ({}) VALUES ({})'.format(','.join(colnames), ','.join(['?']*len(colnames)))
# with db:
#     db.executemany(query, ta)
# print('File added to database.')