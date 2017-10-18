from __future__ import division, print_function

import sqlite3
from astropy.table import Table

db = sqlite3.connect('data/data.db')
db.execute('PRAGMA busy_timeout = 10000')

with db:
    with open('src/sql/cube_create.sql') as f:
        query = f.read()
    db.execute(query)

FNAME = 'data/cube.fits'

print('Reading file '+FNAME)
t = Table.read(FNAME)
print('Found {} rows.'.format(len(t)))

query = 'INSERT INTO cube ({}) VALUES ({})'.format(t.columns, ','.join(['?']*24))
with db:
    db.executemany(query, t)
print('File added to database.')