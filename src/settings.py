import numpy as np

SNDATA = {
    34: {'z': 1.5, 'clr': 'red'},
    38: {'z': 1.0, 'clr': 'green'},
    45: {'z': 0.5, 'clr': 'blue'}
}

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
    'sw-sigmas': lambda _sn, _mname='stellarMass', _observe=False: 'data/{}sw-sigmas-{}{}.csv'.format('observe/' if _observe else '', _mname, _sn),
    'plot-sw-sigmas': lambda _sn, _mname='stellarMass', _observe=False: 'plots/{}sigmas-{}{}.pdf'.format('observe/' if _observe else '', _mname, _sn),

    'hw-reg': 'data/hw-reg.json',
    'hw-rms': lambda _sn, _observe=False: 'data/{}hw-rms{}.fits'.format('observe/' if _observe else '', _sn),
    'hwp-sigmas': lambda _sn, _observe=False: 'data/{}hwp-sigmas{}.csv'.format('observe/' if _observe else '', _sn),
    'hwa-sigmas': lambda _sn, _observe=False: 'data/{}hwa-sigmas{}.csv'.format('observe/' if _observe else '', _sn),
    'hwm-sigmas': lambda _sn, _observe=False: 'data/{}hwm-sigmas{}.csv'.format('observe/' if _observe else '', _sn),

    'data': lambda _sn, _name, _observe=False: 'data/data/{}{}{}.json'.format(_name, '-obs' if _observe else '', _sn),

    'vdist': lambda _sn: 'data/vdist/satstats{}.fits'.format(_sn)
}