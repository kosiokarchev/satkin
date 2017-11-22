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
    'cube': lambda sn, sats=False: 'data/cube{}{}0.fits'.format(sn, '>' if sats else '='),
    'masses': lambda sn: 'data/masses{}.fits'.format(sn),
    'nums': lambda sn: 'data/nums{}.fits'.format(sn),
    '3dbins': lambda sn: 'data/3dbins{}.json'.format(sn),
    'sats': lambda sn: 'data/sats{}.fits'.format(sn),

    'sw-reg': 'data/sw-reg.json',
    'sw-reg-mvir': lambda sn: 'data/sw-reg-mvir{}.csv'.format(sn),
    'sw-sigmas': lambda sn, mname='stellarMass', observe=False: 'data/{}sw-sigmas-{}{}.csv'.format('observe/' if observe else '', mname, sn),
    'plot-sw-sigmas': lambda sn, mname='stellarMass', observe=False: 'plots/{}sigmas-{}{}.pdf'.format('observe/' if observe else '', mname, sn),

    'hw-reg': 'data/hw-reg.json',
    'hw-rms': lambda sn, observe=False: 'data/{}hw-rms{}.fits'.format('observe/' if observe else '', sn),
    'hwp-sigmas': lambda sn, observe=False: 'data/{}hwp-sigmas{}.csv'.format('observe/' if observe else '', sn),
    'hwa-sigmas': lambda sn, observe=False: 'data/{}hwa-sigmas{}.csv'.format('observe/' if observe else '', sn),
    'hwm-sigmas': lambda sn, observe=False: 'data/{}hwm-sigmas{}.csv'.format('observe/' if observe else '', sn),

    'data': lambda sn, name, observe=False: 'data/data/{}{}{}.json'.format(name, '-obs' if observe else '', sn)
}