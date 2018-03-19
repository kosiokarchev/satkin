from funcs import s2h
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import cosmology as cosmo
from galocator_mod import galocator as _galocator


def galocate(t, nvircen=2, nvirsat=1, dvcen=3000, dvsat=3000, fM=1.0, d0=2):
    zmax = np.max(t['z_app'])
    dvcen = dvcen / 3e5
    dvsat = dvsat / 3e5
    sz = dvcen * (1 + zmax)

    t['rho'] = cosmo.Planck13.clone(H0=100).critical_density(t['z_app']).to('solMass/Mpc3').value
    t['rvir'] = (s2h(t['stellarMass']) / ((4*np.pi/3)*200*t['rho']))**(1/3)
    t['d_a'] = t['d_comoving'] / (1+t['z_app'])

    t['binz'] = t['z_app'] // sz
    t.sort('binz')
    zindex = np.append(np.unique(t['binz'], return_index=True)[1], len(t))

    dcen_max = nvircen * t['rvir'] / t['d_a']
    dzcen_max = dvcen * (1 + t['z_app'])

    dsat_max = nvirsat * t['rvir'] / t['d_a']
    dzsat_max = dvsat * (1 + t['z_app'])

    args = [t['ra'].data.astype(np.double), t['dec'].data.astype(np.double), t['z_app'].data.astype(np.double),
            t['stellarMass'].data.astype(np.double),
            dcen_max.astype(np.double), dzcen_max.astype(np.double),
            dsat_max.astype(np.double), dzsat_max.astype(np.double),
            t['galaxyId'].data.astype(np.int64),
            zindex,
            d0*np.max(dcen_max), fM]

    return _galocator(*args)

def satlocate(t, ccol, nvir=1, dv=3000):
    cid = ccol+'Id'

    c = t[t[ccol]]
    s = t[~t[ccol]]
    s[cid] = c['galaxyId'].dtype.type(0)

    ccoords = SkyCoord(c['ra'], c['dec'], unit='rad')
    scoords = SkyCoord(s['ra'], s['dec'], unit='rad')

    c['dzmax'] = (dv/3e5) * (1+c['z_app'])
    c['dmax'] = nvir * (1 + c['z_app']) * c['rvir'] / c['d_comoving']


    for z, dzmax, dmax, cc, gId, i in zip(c['z_app'], c['dzmax'], c['dmax'], ccoords, c['galaxyId'], range(len(c))):
        inzindices = np.where(np.abs(s['z_app'] - z) < dzmax)[0]
        satindices = inzindices[np.where(cc.separation(scoords[inzindices]).to('rad').value < dmax)]

        s[cid][satindices] = gId
        if i%1000 == 0: print(i)

    return s[s[cid] != 0]
