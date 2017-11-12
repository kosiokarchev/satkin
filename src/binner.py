from lib import *
import pickle


def binX(t, x, y, nbins=20, bins=None, funcsx=None, funcsy=None):
    if funcsx is None:
        funcsx = {'mean': np.mean, 'std': np.std}
    elif callable(funcsx):
        funcsx = {'stat': funcsx}

    if funcsy is None:
        funcsy = {'mean': np.mean, 'std': np.std}
    elif funcsy is True:
        funcsy = funcsx
    elif callable(funcsy):
        funcsy = {'stat': funcsy}


    if bins is not None:
        width = bins[1] - bins[0]
        left = bins[0]
        nbins = len(bins) - 1
    else:
        left = np.min(t[x])
        width = (np.max(t[x]) - left) / nbins

    t['bin'] = np.floor_divide(t[x] - left, width).astype(int)
    g = t.group_by('bin')

    statsx = []
    statsy = []
    for i in range(len(g.groups.keys)):
        if 0 <= g.groups.keys[i]['bin'] < nbins:
            statsx.append([f(g[x].groups[i]) for f in funcsx.values()])
            statsy.append([f(g[y].groups[i]) for f in funcsy.values()])

    return ({name: stat for name, stat in zip(funcsx.keys(), np.array(statsx).transpose())},
            {name: stat for name, stat in zip(funcsy.keys(), np.array(statsy).transpose())})

def bin2D(t, x, y, nbinsx=20, nbinsy=20, binsx=None, binsy=None, funcsx=None, funcsy=None):
    print('Binning in', x)
    binx = binX(t, x, y, nbinsx, binsx, funcsx, funcsy)

    print('Binning in', y)
    biny = binX(t, y, x, nbinsy, binsy, funcsx, funcsy)
    biny = biny[1], biny[0]

    return binx, biny


class bin2d:
    def __init__(self, t=None, x=None, y=None, binx=None, biny=None):
        self.binx = binx
        self.biny = biny

        self.t = t
        self.x = x
        self.y = y

    def load(self, fname):
        with open(fname, 'rb') as f:
            binx, biny = pickle.load(f, encoding='bytes')

        for b in binx, biny:
            for c in b:
                for key in tuple(c.keys()):
                    if isinstance(key, bytes):
                        c[key.decode()] = c[key]
                        del c[key]

        self.binx, self.biny = binx, biny
    def save(self, fname):
        with open(fname, 'wb') as f:
            pickle.dump((self.binx, self.biny), f)

    def compute(self, binsx, binsy):
        self.binx, self.biny = bin2D(
            self.t, self.x, self.y, binsx=binsx, binsy=binsy,
            funcsx={'mean': np.mean, 'std': np.std, 'N': len},
            funcsy={'mean': np.mean, 'std': np.std, 'median': np.median,
                    'p2': percentile_wrapper(2),
                    'p16': percentile_wrapper(16),
                    'p84': percentile_wrapper(84),
                    'p98': percentile_wrapper(98)})

    def plot_data(self, sparse=1, label='_nolegend_'):
        X, Y = self.t[self.x], self.t[self.y]
        if sparse > 1:
            r = (np.random.random(int(len(self.t) / sparse)) * len(self.t)).astype(int)
            X, Y = X[r], Y[r]
        plt.plot(X, Y, 'r,', zorder=-1, label=label)
    def plot_bins(self,
                  label1='x binning', label2='y binning',
                  clr1='blue', clr2='yellow',
                  narrow=True, wide=False, median=False, ebars=False, mean=True):
        def getkw(o, i=None):
            return {} if o is True else getkw(o[i]) if isinstance(o, (list, tuple)) else o


        defaults = [
            [ # wide
                dict(color=clr1, alpha=0.25, label='_nolegend_'),
                dict(color=clr2, alpha=0.25, label='_nolegend_')
            ],
            [ # narrow
                dict(color=clr1, alpha=0.5, label=label1),
                dict(color=clr2, alpha=0.5, label=label2),
            ],
            [ # median
                dict(marker='o', linestyle='None', color=clr1, markeredgecolor='red'),
                dict(marker='o', linestyle='None', color=clr2, markeredgecolor='red')
            ],
            [dict(color='black', linestyle='None')]*2, # std
            [dict(color='black')]*2  # mean
        ]

        ds = (wide, narrow, median, ebars, mean)
        def ev(i, j, func, *args):
            ukw = getkw(ds[j], i)
            if ukw is not False:
                kw = defaults[j][i]
                kw.update(ukw)
                func(*args, **kw)

        for i, (x, y), (xm, ym), fill in zip(
                (0, 1),
                (self.binx, self.biny),
                (('mean', 'median'), ('median', 'mean')),
                (plt.fill_between, plt.fill_betweenx)):
            fillx = x if i==0 else y
            filly = y if i==0 else x

            ev(i, 0, fill, fillx['mean'], filly['p2'], filly['p98'])
            ev(i, 1, fill, fillx['mean'], filly['p16'], filly['p84'])
            ev(i, 2, plt.plot, x[xm], y[ym])
            ev(i, 3, plt.errorbar, x['mean'], y['mean'], y['std'], x['std'])
            ev(i, 4, plt.plot, x['mean'], y['mean'])
    def plot_amatch(self, label='abundance matching'):
        X, Y = self.t[self.x].copy(), self.t[self.y].copy()
        X.sort(), Y.sort()
        plt.plot(X, Y, 'r-', label=label)
    def plot(self, label1, label2, sparse=1, bins=True, data=False, amatch=False):
        if data:
            self.plot_data(sparse)
        if bins:
            self.plot_bins(label1, label2)
        if amatch:
            self.plot_amatch()


def build3Dbins(sn):
    fname = FILES['masses'](sn)
    print('Loading', fname)
    t = Table.read(fname)
    t = deal(t)

    b = bin2d(t, 'mv', 'ms')
    b.compute(PLOTDATA['bins']['mvir'], PLOTDATA['bins']['stellarMass'])
    b.save(FILES['3dbins'](sn))