import json, pickle
import numpy as np
from scipy.interpolate import RegularGridInterpolator as interpnd
import lmfit

from matplotlib import pyplot as plt

from funcs import schechter_func, SchechterModel, two_lines_2, mv2s2, percentile, logpowerlaw

sn = 38
files = {
    'p_obs': 'cones/p_obs.npz',
    'f(mv_sat)': 'data_aux/f(mvir_sat).sav',
    'sigma2(mv)': 'data_aux/mvir2sigma2{}.sav'.format(sn),
    's2fits': 'data/theory/s2fits.pickle'
}
fits = {
    'f(mv_sat)': lmfit.model.load_modelresult(files['f(mv_sat)'],
                                              {'two_lines_2': two_lines_2}),
    'sigma2(mv)': lmfit.model.load_modelresult(files['sigma2(mv)'],
                                               {'logpowerlaw': logpowerlaw}),
}

s2fits = pickle.load(open(files['s2fits'], 'rb'))

def save_parameters(par, fname):
    with open(fname, 'w') as f:
        par.dump(f)
def load_parameters(fname):
    p = lmfit.Parameters()
    with open(fname, 'r') as f:
        p.load(f)
    return p
def load_last_params(fname):
    with open(fname) as f:
        p = json.load(f)
    return lmfit.Parameters().loads(json.dumps(p[-1]))
def extract_params(prefix, params, includeexpr=True):
    return {key[len(prefix):]: params[key].value
            for key in params
            if key.startswith(prefix)
            and (includeexpr or not params[key].expr)}

class ThePlot:
    @staticmethod
    def f_mv_v(mv, a, x0, b, cutoff):
        s = np.power(10, schechter_func(mv, a, x0, b, 0))
        return np.where(mv < cutoff, 0, s)

    _cache_fmv = None
    def f_mv_t(self, mv, a, x0, b, cutoff):
        return self._cache_fmv

    f_mv = f_mv_v

    def set_f_mv(self, mv, a, x0, b, cutoff, *args, **kwargs):
        self._cache_fmv = self.f_mv(mv, a, x0, b, cutoff).reshape((-1, 1))
        self.f_mv = self.f_mv_t

    @staticmethod
    def N(mv, a, N12, p, f):
        return np.power(10, a * (mv - 12) + N12)

    @staticmethod
    def sigma2(mv):
        return np.power(10, s2fits[38].eval(x=mv))
        # return mv2s2(mv)
        # return np.power(10, fits['sigma2(mv)'].eval(x=mv))
        # return np.power(10, s2_model.eval(s2_params, x=mv))

    @staticmethod
    def mean_ms(mv, mv0, ms0, a1, a2, scale):
        return two_lines_2(mv, mv0, ms0, a1, a2, scale)

    @staticmethod
    def spread(mv, a, s12):
        return a * (mv - 12) + s12

    def G(self, mv, ms, mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale, spread_a, spread_s12):
        s = self.spread(mv, spread_a, spread_s12)
        mms = self.mean_ms(mv, mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale)
        return np.exp(- (ms - mms) ** 2 / (2 * s ** 2)) / (
                    np.sqrt(2 * np.pi) * s)

    def theplot(self, mv, ms,
                mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale,
                spread_a, spread_s12,
                f_a, f_x0, f_b, f_cutoff):
        mv = mv.reshape((-1, 1))
        ms = ms.reshape((1, -1))
        ret = (self.f_mv(mv, f_a, f_x0, f_b, f_cutoff)
               * self.G(mv, ms,
                        mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale,
                        spread_a, spread_s12))
        return ret / ret.sum() / (mv[1, 0] - mv[0, 0]) / (ms[0, 1] - ms[0, 0])

    def get_f_p(self, mv, ms,
                mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale,
                spread_a, spread_s12,
                f_a, f_x0, f_b, f_cutoff):
        f = self.f_mv(mv, f_a, f_x0, f_b, f_cutoff)
        p = f * self.G(mv, ms,
                       mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale,
                       spread_a, spread_s12)
        return f, p

    @staticmethod
    def get_observables(p, N, s2, out):
        pN = p * N
        pNnorm = np.nansum(pN, axis=0)

        pP = p * (1 - np.exp(-N))
        pPnorm = np.nansum(pP, axis=0)

        out[0, :] = np.nansum(pN * s2, axis=0) / pNnorm
        out[1, :] = np.nansum(pP * s2, axis=0) / pPnorm
        out[2, :] = pNnorm / pPnorm

        out[np.isnan(out)] = 0
        return out

    def observables(self, out, mv, s2, ms,
                    mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale,
                    spread_a, spread_s12,
                    f_a, f_x0, f_b, f_cutoff,
                    N_a, N_N12):
        mv = mv.reshape((-1, 1))
        ms = ms.reshape((1, -1))
        s2 = s2.reshape((-1, 1))

        f, p = self.get_f_p(
            mv, ms,
            mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale,
            spread_a, spread_s12,
            f_a, f_x0, f_b, f_cutoff)
        N = self.N(mv, N_a, N_N12, p, f)

        return self.get_observables(p, N, s2, out)


class ObservedThePlot(ThePlot):
    f_mv_sat = fits['f(mv_sat)']

    po = np.load(files['p_obs'])

    po_ms_z = po['ho'].sum(axis=2) / po['ha'].sum(axis=2)
    po_ms_z = np.where(np.isnan(po_ms_z), 0, po_ms_z)
    po_ms_i = interpnd(((po['ms'][1:] + po['ms'][:-1]) / 2,
                        (po['z'][1:] + po['z'][:-1]) / 2),
                       po_ms_z, bounds_error=False, fill_value=0)

    def __init__(self, z, ms, mv):
        self.dms = ms[1] - ms[0]
        self.z = z
        self.po_ms = self.po_ms_i(np.array((ms, [self.z] * len(ms))).T)

        mv = mv.reshape((-1, 1))
        pars = self.f_mv_sat.best_values.copy()
        pars['x0'] = mv + pars['x0']
        self.f = np.power(10, self.f_mv_sat.model.func(mv, **pars))

    def p_mv0_mv(self, f_mv):
        ret = f_mv * self.f
        ret = ret / ret.sum()
        return np.where(np.isnan(ret), 0, ret)

    def fn(self, mv, p, f_mv):
        p_mv0_mv = self.p_mv0_mv(f_mv.ravel())
        fn = ((p_mv0_mv / p_mv0_mv.sum(axis=1, keepdims=True))[:, :, np.newaxis]
              * (p / p.sum(axis=1, keepdims=True))[np.newaxis, :, :]
              * self.po_ms[np.newaxis, np.newaxis, :])
        return np.nansum(fn, axis=(1, 2)).reshape(mv.shape)

    def N(self, mv, a, N12, p, f):
        return self.fn(mv, p, f) * super(ObservedThePlot, self).N(mv, a, N12, p, f)


class ThePlotModel:
    p_init = lmfit.Parameters()
    p_init.add('mms_mv0', 12, min=9, max=15, brute_step=0.1)    # 60
    p_init.add('mms_ms0', 10, min=7, max=12, brute_step=0.1)    # 50
    p_init.add('mms_a1', 2, min=0, max=5, brute_step=0.1)       # 50
    # p_init.add('adiff', 1.7, min=0, max=5)
    p_init.add('mms_a2', 0.3, min=0, max=5, brute_step=0.1)     # 50
    p_init.add('mms_scale', 1, min=0, max=5, vary=False)
    p_init.add('spread_a', 0, vary=False)
    p_init.add('spread_s12', 0.25, min=0, max=2, brute_step=0.1)  # 20


    p_init.add('f_a', -1.64, min=-5, max=-1, brute_step=0.1)     # 40
    p_init.add('f_x0', 11.69)
    p_init.add('f_b', 0.28, min=0, max=5)
    p_init.add('f_cutoff', 10.3)

    p_init.add('slope_1', p_init['mms_a1'] / p_init['mms_scale'], expr='mms_a1 / mms_scale', max=5)
    p_init.add('slope_2', p_init['mms_a1'] / p_init['mms_scale'], expr='mms_a2 / mms_scale', max=2)
    # p_init.add('spread_low', p_init['spread_a'] * (9-12) + p_init['spread_s12'],
    #            expr='spread_a * (9-12) + spread_s12', min=0)
    # p_init.add('spread_high', p_init['spread_a'] * (15-12) + p_init['spread_s12'],
    #            expr='spread_a * (15-12) + spread_s12', min=0)

    params_names = ['mms_mv0', 'mms_ms0', 'mms_a1', 'mms_a2', 'mms_scale', 'spread_a', 'spread_s12', 'f_a', 'f_x0', 'f_b', 'f_cutoff']
    independent_vars = ['mv', 'ms']

    def __init__(self, data, weights=None, mask=None, theplotter=ThePlot(), **kwargs):
        self.data = data
        self.indep = kwargs
        self.weights = weights

        self.funcargs = {key: kwargs[key] for key in self.independent_vars}
        for p in self.params_names:
            self.funcargs[p] = None

        self.tpc = theplotter
        self.mask = None if mask is None else mask.ravel()

        self.mini = None

    def fix_fmv(self, params):
        self.tpc.set_f_mv(mv=self.get_mv(),
                          a=params['f_a'], x0=params['f_x0'],
                          b=params['f_b'], cutoff=params['f_cutoff'])
        params['f_a'].vary = params['f_x0'].vary = params['f_b'].vary = params['f_cutoff'].vary = False

    def func(self, *args, **kwargs):
        return self.tpc.theplot(*args, **kwargs)

    def get_mv(self):
        return self.indep['mv']

    def eval(self, params):
        for key in self.params_names:
            self.funcargs[key] = params[key].value
        return self.func(**self.funcargs)

    def fit(self, params=None, method='leastsq', varyf=True, **kwargs):
        if params is None:
            params = self.p_init.copy()
        if not varyf:
            self.fix_fmv(params)
        self.mini = lmfit.Minimizer(self.residual, params, **kwargs)
        res = self.mini.minimize(method, params)
        self.tpc.f_mv = self.tpc.f_mv_v
        return res

    def chi2(self, params):
        r = self.residual(params) ** 2
        return np.sum(r[np.isfinite(r)])

    def residual(self, params):
        res = self.eval(params) - self.data
        if self.weights is not None:
            res *= self.weights
        res = res.ravel()
        return res if self.mask is None else res[self.mask]


class LogThePlotModel(ThePlotModel):
    def func(self, *args, **kwargs):
        return np.log10(super(LogThePlotModel, self).func(*args, **kwargs))


class ThePlotSigmasModel(ThePlotModel):
    p_init = ThePlotModel.p_init.copy()
    p_init.add('N_a', 0.9, min=0, max=5, brute_step=0.1)  #50
    p_init.add('N_N12', 0, min=-2, max=2, brute_step=0.1, vary=False)  #20
    p_init.add('fNo', 1.5, min=0)



    params_names = ThePlotModel.params_names + ['N_a', 'N_N12', 'fNo']
    independent_vars = ['ms']

    def __init__(self, dM=0.1, low=8, high=15, **kwargs):
        super(ThePlotSigmasModel, self).__init__(**kwargs)

        self.mv = (np.arange(low, high, dM) + dM/2)[:, np.newaxis]
        self.s2 = self.tpc.sigma2(self.mv)
        self.y = None

    def get_mv(self):
        return self.mv

    def func(self, ms,
             mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale,
             spread_a, spread_s12,
             f_a, f_x0, f_b, f_cutoff,
             N_a, N_N12, fNo):
        if self.y is None or self.y.shape[1] != len(ms):
            self.y = np.empty((3, len(ms)))
        self.tpc.observables(self.y, self.mv, self.s2,
                             ms[np.newaxis, :],
                             mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale,
                             spread_a, spread_s12,
                             f_a, f_x0, f_b, f_cutoff,
                             N_a, N_N12)
        self.y[2] *= fNo
        return self.y


class ThePlotLogSigmasModel(ThePlotSigmasModel):
    def func(self, *args, **kwargs):
        return np.log10(super(ThePlotLogSigmasModel, self).func(*args, **kwargs))


class FitIterCallback:
    def __init__(self, params, plotfig, plotter, doparams=True, param_keys=None, param_labels=None, doresid=True, saveparams=None):
        self.param_keys = tuple(params.keys()) if param_keys is None else param_keys
        self.param_labels = self.param_keys if param_labels is None else param_labels

        self.resids = [np.inf]
        self.params = [params]

        self.plotfig = plotfig
        self.fitlines = None
        self.plotter = plotter

        self.doparams = doparams
        if doparams:
            self.paramfig, self.paramax = plt.subplots(len(self.param_keys),
                                                       sharex='col', squeeze=True,
                                                       gridspec_kw=dict(hspace=0))
            for k, ax in zip(self.param_labels, self.paramax):
                ax.set_ylabel(k)
            self.plot_params(params)

        self.doresid = doresid
        if doresid:
            self.residfig = plt.figure()
            self.residax = self.residfig.gca()
            self.residax.set_yscale('log')
            self.residax.set_ylabel('$\chi^2$')

        self.saveparams = saveparams

    def plot_params(self, params, i=0):
        for k, ax in zip(self.param_keys, self.paramax):
            ax.plot(i, params[k].value, '.')
        self.paramfig.canvas.draw()

    def plot_resid(self, resid, i=0):
        self.residax.plot(i, resid, '.')
        self.residfig.canvas.draw()

    def plot_plot(self, r, params):
        if self.fitlines:
            for fl in self.fitlines:
                try:
                    fl.remove()
                except:
                    pass
        self.fitlines = self.plotter(r, params)
        self.plotfig.canvas.draw()

    def __call__(self, params, i, resid, *args, **kws):
        self.resids.append(np.sum(resid**2))
        self.params.append(params)

        print(i, ':', self.resids[-1])
        if self.saveparams:
            with open(self.saveparams, 'w') as f:
               f.write('[' + ', '.join([p.dumps() for p in self.params]) + ']')
        if self.doparams:
            self.plot_params(params, i)
        if self.doresid:
            self.plot_resid(self.resids[-1], i)

        self.plot_plot(resid, params)
        plt.pause(0.01)


def fit_theplot(mv, h, mms, sms, binwidth=0.01, smssub=(50, 350), plot=False):
    p = h / h.sum() / binwidth ** 2
    fmv = p.sum(axis=1) * binwidth

    params = ThePlotModel.p_init.copy()

    fmvparamnames = ('a', 'x0', 'b')
    fmvmod = SchechterModel(nan_policy='omit')
    fmvparams = fmvmod.make_params(c=0,
                                   **{pn: params['f_' + pn].value
                                      for pn in fmvparamnames})
    fmvfr = fmvmod.fit(np.log10(fmv), fmvparams, x=mv)

    mmsmod = lmfit.Model(ThePlot.mean_ms, nan_policy='omit')
    mmsparams = mmsmod.make_params(**extract_params('mms_', params))
    mmsparams.add('slope_1', mmsparams['a1'] / mmsparams['scale'],
                  expr='a1 / scale')
    mmsparams.add('slope_2', mmsparams['a2'] / mmsparams['scale'],
                  expr='a2 / scale')
    mmsfr = mmsmod.fit(mms, params=mmsparams, mv=mv)

    for prefix, pars in zip(('f_', 'mms_'), (fmvfr.params, mmsfr.params)):
        for pn in pars:
            if pars[pn].expr:
                params[pn].value = pars[pn].value
                params[pn].stderr = pars[pn].stderr
            else:
                params[prefix + pn] = pars[pn]

    smssample = sms[smssub[0]:smssub[1]]
    smsm = np.mean(smssample)
    smss = np.std(smssample)

    params['spread_s12'].value = smsm
    params['spread_s12'].stderr = smss

    if plot:
        plt.figure()
        plt.plot(mv, np.log10(fmv), '.')
        plt.plot(mv, fmvmod.eval(fmvfr.params, x=mv))

        plt.figure()
        plt.plot(mv, mmsmod.eval(mmsfr.params, mv=mv), '--')
        plt.errorbar(mv[5:], mms[5:], yerr=sms[5:], errorevery=10, capthick=1)

        plt.figure()
        plt.plot(mv, sms)
        plt.hlines((smsm, smsm + smss, smsm - smss), 11, 14,
                   linestyles=('--', ':', ':'))

    return params

def examine_theplot(p, mv, ms, binwidth):
    f_mv = p.sum(axis=1)
    f_ms = p.sum(axis=0)

    vals = (ms+0.5*binwidth)[np.newaxis, :]
    mean = np.sum(p * vals, axis=1) / f_mv
    mean2 = np.sum(p * vals**2, axis=1) / f_mv
    mms = mean
    sms = np.sqrt(mean2 - mean**2)
    pms = percentile(ms+0.5*binwidth, p / f_mv[:, np.newaxis], [0.16, 0.5, 0.84])

    vals = (mv+0.5*binwidth)[:, np.newaxis]
    mean = np.sum(p * vals, axis=0) / f_ms
    mean2 = np.sum(p * vals**2, axis=0) / f_ms
    mmv = mean
    smv = np.sqrt(mean2 - mean**2)
    pmv = percentile(mv+0.5*binwidth, (p / f_ms[np.newaxis, :]).T, [0.16, 0.5, 0.84])

    return dict(mms=mms, sms=sms, pms=pms,
                mmv=mmv, smv=smv, pmv=pmv)

