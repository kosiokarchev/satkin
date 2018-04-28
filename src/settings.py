import numpy as np

SNDATA = {
    30: {'z': 2.0, 'z_exact': 2.070046, 'gclr': '0.75', 'clr': 'firebrick', 'marker': 'o'},
    34: {'z': 1.5, 'z_exact': 1.481728, 'gclr': '0.66', 'clr': 'red', 'marker': 'D'},
    38: {'z': 1.0, 'z_exact': 1.041466, 'gclr': '0.5', 'clr': 'green', 'marker': 's'},
    45: {'z': 0.5, 'z_exact': 0.513287, 'gclr': '0.33', 'clr': 'blue', 'marker': '^'},
    58: {'z': 0.0, 'z_exact': 0.000264, 'gclr': '0.25','clr': 'purple', 'marker': 'v'},

    41: {'z': 0.8, 'z_exact': 0.783975, 'gclr': '0.4', 'clr': 'cyan'},
    50: {'z': 0.25, 'z_exact': 0.262623, 'gclr': '0.3', 'clr': 'magenta'}
}
THESNS = (30, 34, 38, 45, 58)

BINWIDTH = 0.1
PLOTDATA = {
    'bins': {
        'mvir': np.arange(10, 15 + BINWIDTH, BINWIDTH),
        'stellarMass': np.arange(7, 12 + BINWIDTH, BINWIDTH)
    }
}

FILES = {
    'cube': lambda _sn, _sats=False: 'data/cube{}{}0.fits'.format(_sn, '>' if _sats else '='),
    'masses': lambda _sn: 'data/masses{}.fits'.format(_sn),
    'nums': lambda _sn: 'data/nums{}.fits'.format(_sn),
    '3dbins': lambda _sn: 'data/3dbins{}.json'.format(_sn),
    'sats': lambda _sn: 'data/sats{}.fits'.format(_sn),

    'sw-reg': 'data/sw-reg.json',
    'sw-reg-mvir': lambda _sn: 'data/sw-reg-mvir{}.csv'.format(_sn),
    'sw-sigmas': lambda _sn, _observe=False, _mname='stellarMass': 'data/{}sw-sigmas-{}{}.csv'.format('observe/' if _observe else '', _mname, _sn),
    'plot-sw-sigmas': lambda _sn, _mname='stellarMass', _observe=False: 'plots/{}sigmas-{}{}.pdf'.format('observe/' if _observe else '', _mname, _sn),

    'hw-reg': 'data/hw-reg.json',
    'hw-rms': lambda _sn, _observe=False: 'data/{}hw-rms{}.fits'.format('observe/' if _observe else '', _sn),
    'hwp-sigmas': lambda _sn, _observe=False: 'data/{}hwp-sigmas{}.csv'.format('observe/' if _observe else '', _sn),
    'hwa-sigmas': lambda _sn, _observe=False: 'data/{}hwa-sigmas{}.csv'.format('observe/' if _observe else '', _sn),
    'hwm-sigmas': lambda _sn, _observe=False: 'data/{}hwm-sigmas{}.csv'.format('observe/' if _observe else '', _sn),

    'data': lambda _sn, _name, _observe=False: 'data/data/{}{}{}.json'.format(_name, '-obs' if _observe else '', _sn),

    'vdist': lambda _sn: 'data/vdist/satstats{}.fits'.format(_sn),

    'cone': lambda _i=1, _sats=None: 'cones/cone{:03}{}.fits'.format(_i, '' if _sats is None else '>0' if _sats else '=0'),
    'hst': '3d-hst/3d-hst.fits',
    'hst-sw': '3d-hst/theplot-sw.fits',

    'numbins': lambda _sn: 'data/nums/numbins{}.fits'.format(_sn),
    'numstats': lambda _sn: 'data/nums/numstats{}.fits'.format(_sn),

    'p_obs': 'cones/p_obs.npz',
    'theplot': lambda _sn: 'data/theplot{}.npy'.format(_sn),
    'theplot_data': 'data/theplot_data.pickle',
    'theplot_params': 'data/theplot_params.fits',
    'theplot_mvfitp': lambda _sn: 'data_aux/theplot_mvfit{}.fits'.format(_sn),
    'theplot_msfitp': lambda _sn: 'data_aux/theplot_msfit{}.fits'.format(_sn),
    'phase': lambda _sn: 'data/phase{}.npy'.format(_sn),
    'phase_obs': lambda _sn: 'cones/phase_obs{}.npy'.format(_sn),
    'N(mv)': 'data_aux/num_sat(mvir).sav',
    'No(mv)': lambda _sn: 'data/nums{}_obs_mv.fits'.format(_sn),
    'f(mv_sat)': 'data_aux/f(mvir_sat).sav',
    'f(mv)': lambda _sn: 'data_aux/f(mv){}.sav'.format(_sn)
}