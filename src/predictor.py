import numpy as np
from scipy.integrate import quad
from scipy.interpolate import RegularGridInterpolator as interpnd
import lmfit

from matplotlib import pyplot as plt

from funcs import schechter_func, two_lines_2

sn = 38
files = {
    'p_obs': 'cones/p_obs.npz',
    'f(mv_sat)': 'data_aux/f(mvir_sat).sav',
    'sigma2(mv)': 'data_aux/mvir2sigma2{}.sav'.format(sn)
}
fits = {
    'f(mv_sat)': lmfit.model.load_modelresult(files['f(mv_sat)'],
                                              {'two_lines_2': two_lines_2}),
    'sigma2(mv)': lmfit.model.load_modelresult(files['sigma2(mv)'])
}


def save_parameters(par, fname):
    with open(fname, 'w') as f:
        par.dump(f)
def load_parameters(fname):
    p = lmfit.Parameters()
    with open(fname, 'r') as f:
        p.load(f)
    return p

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
    def N(mv, a, N_12, p, f):
        return np.power(10, a * (mv - 12) + N_12)

    @staticmethod
    def sigma2(mv):
        return np.power(10, fits['sigma2(mv)'].eval(x=mv))

    @staticmethod
    def mean_ms(mv, mv0, ms0, a1, a2, scale):
        return two_lines_2(mv, mv0, ms0, a1, a2, scale)

    @staticmethod
    def spread(mv, a, spread_12):
        return a * (mv - 12) + spread_12

    def G(self, mv, ms, mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale, spread_a,
          spread_12):
        s = self.spread(mv, spread_a, spread_12)
        mms = self.mean_ms(mv, mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale)
        return np.exp(- (ms - mms) ** 2 / (2 * s ** 2)) / (
                    np.sqrt(2 * np.pi) * s)

    def theplot(self, mv, ms,
                mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale,
                spread_a, spread_12,
                f_a, f_x0, f_b, f_cutoff):
        mv = mv.reshape((-1, 1))
        ms = ms.reshape((1, -1))
        ret = (self.f_mv(mv, f_a, f_x0, f_b, f_cutoff)
               * self.G(mv, ms,
                        mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale,
                        spread_a, spread_12))
        return ret / ret.sum() / (mv[1, 0] - mv[0, 0]) / (ms[0, 1] - ms[0, 0])

    def observables(self, out, mv, s2, ms,
                    mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale,
                    spread_a, spread_12,
                    f_a, f_x0, f_b, f_cutoff,
                    N_a, N_12):
        mv = mv.reshape((-1, 1))
        ms = ms.reshape((1, -1))

        f = self.f_mv(mv, f_a, f_x0, f_b, f_cutoff)
        p = f * self.G(mv, ms,
                       mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale,
                       spread_a, spread_12)
        N = self.N(mv, N_a, N_12, p, f)

        pN = p * N
        pNnorm = np.nansum(pN, axis=0)

        pP = p * (1 - np.exp(-N))
        pPnorm = np.nansum(pP, axis=0)

        out[0, :] = np.nansum(pN * s2, axis=0) / pNnorm
        out[1, :] = np.nansum(pP * s2, axis=0) / pPnorm
        out[2, :] = pNnorm / pPnorm

        out[np.isnan(out)] = 0
        return out


# class ThePlot:
#     @staticmethod
#     def f_mv_v(mv, a, x0, b, cutoff):
#         s = np.power(10, schechter_func(mv, a, x0, b, 0))
#         # s /= quad(lambda x: np.power(10, schechter_func(x, a, x0, b, 0)),
#         #           cutoff, np.inf)[0]
#         # s[~np.isfinite(s)] = 0
#         return np.where(mv < cutoff, 0, s)
#
#     _cache_fmv = None
#     @classmethod
#     def f_mv_t(cls, mv, a, x0, b, cutoff):
#         return cls._cache_fmv
#
#     f_mv = f_mv_v
#
#     @classmethod
#     def set_f_mv(cls, mv, a, x0, b, cutoff, *args, **kwargs):
#         cls._cache_fmv = cls.f_mv(mv, a, x0, b, cutoff)[:, np.newaxis]
#         cls.f_mv = cls.f_mv_t
#
#     @staticmethod
#     def N(mv, a, N_12, p):
#         return np.power(10, a * (mv-12) + N_12)
#
#     @staticmethod
#     def sigma2(mv):
#         return np.power(10, fits['sigma2(mv)'].eval(x=mv))
#
#     @staticmethod
#     def mean_ms(mv, mv0, ms0, a1, a2, scale):
#         return two_lines_2(mv, mv0, ms0, a1, a2, scale)
#
#     @staticmethod
#     def spread(mv, a, spread_12):
#         return a * (mv-12) + spread_12
#
#     @classmethod
#     def G(cls, mv, ms, mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale, spread_a, spread_12):
#         s = cls.spread(mv, spread_a, spread_12)
#         mms = cls.mean_ms(mv, mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale)
#         return np.exp(- (ms - mms) ** 2 / (2 * s ** 2)) / (np.sqrt(2 * np.pi) * s)
#
#     @classmethod
#     def theplot(cls, mv, ms,
#                 mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale,
#                 spread_a, spread_12,
#                 f_a, f_x0, f_b, f_cutoff):
#         if len(mv.shape) == 1:
#             mv = mv[:, np.newaxis]
#         if len(ms.shape) == 1:
#             ms = ms[np.newaxis, :]
#         ret = (cls.f_mv(mv, f_a, f_x0, f_b, f_cutoff)
#                * cls.G(mv, ms,
#                        mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale,
#                        spread_a, spread_12))
#         return ret / ret.sum() / (mv[1, 0] - mv[0, 0]) / (ms[0, 1] - ms[0, 0])
#
#     @classmethod
#     def sigmas2(cls, out, mv, s2, ms,
#                 mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale,
#                 spread_a, spread_12,
#                 f_a, f_x0, f_b, f_cutoff,
#                 N_a, N_12):
#         p = cls.theplot(mv, ms,
#                         mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale, spread_a, spread_12,
#                         f_a, f_x0, f_b, f_cutoff)
#
#         pN = p * cls.N(mv, N_a, N_12)
#         out[0, :] = (pN * s2).sum(axis=0) / pN.sum(axis=0)
#         out[1, :] = (p * s2).sum(axis=0) / p.sum(axis=0)
#
#         out[np.isnan(out)] = 0
#
#         return out
#
#     @classmethod
#     def observables(cls, out, mv, s2, ms,
#                     mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale,
#                     spread_a, spread_12,
#                     f_a, f_x0, f_b, f_cutoff,
#                     N_a, N_12):
#         p = cls.theplot(mv, ms,
#                         mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale, spread_a,
#                         spread_12,
#                         f_a, f_x0, f_b, f_cutoff)
#         pnorm = p.sum(axis=0)
#         pN = p * cls.N(mv, N_a, N_12)
#         pNnorm = pN.sum(axis=0)
#
#         out[0, :] = (pN * s2).sum(axis=0) / pNnorm
#         out[1, :] = (p * s2).sum(axis=0) / pnorm
#         out[2, :] = pNnorm / pnorm
#
#         out[np.isnan(out)] = 0
#         return out


class ObservedThePlot(ThePlot):
    f_mv_sat = fits['f(mv_sat)']

    po = np.load(files['p_obs'])

    po_ms_z = po['ho'].sum(axis=2) / po['ha'].sum(axis=2)
    po_ms_z = np.where(np.isnan(po_ms_z), 0, po_ms_z)
    po_ms_i = interpnd(((po['ms'][1:] + po['ms'][:-1]) / 2,
                        (po['z'][1:] + po['z'][:-1]) / 2),
                       po_ms_z, bounds_error=False, fill_value=0)

    def __init__(self, z, ms):
        self.z = z
        self.po_ms = self.po_ms_i(np.array((ms, [self.z] * len(ms))).T)

    @classmethod
    def p_mv0_mv(cls, mv, f_mv):
        pars = cls.f_mv_sat.best_values.copy()
        pars['x0'] = mv[:, np.newaxis] + pars['x0']
        logf = cls.f_mv_sat.model.func(mv[np.newaxis, :], **pars)
        ret = f_mv * np.power(10, logf)
        ret = ret / ret.sum()
        return np.where(np.isnan(ret), 0, ret)

    def fn(self, mv, p, f_mv):
        p_mv0_mv = self.p_mv0_mv(mv.ravel(), f_mv.ravel())
        fn = ((p_mv0_mv / p_mv0_mv.sum(axis=1, keepdims=True))[:, :, np.newaxis]
              * (p / p.sum(axis=1, keepdims=True))[np.newaxis, :, :]
              * self.po_ms[np.newaxis, np.newaxis, :])
        return np.nansum(fn, axis=(1, 2)).reshape(mv.shape)

    def N(self, mv, a, N_12, p, f):
        return self.fn(mv, p, f) * super(ObservedThePlot, self).N(mv, a, N_12, p, f)


class ThePlotModel:
    p_init = lmfit.Parameters()
    p_init.add('mms_mv0', 12, min=9, max=15, brute_step=0.1)    # 60
    p_init.add('mms_ms0', 10, min=7, max=12, brute_step=0.1)    # 50
    p_init.add('mms_a1', 2, min=0, max=5, brute_step=0.1)       # 50
    p_init.add('mms_a2', 0.3, min=0, max=5, brute_step=0.1)     # 50
    p_init.add('mms_scale', 1, min=0, max=5)
    p_init.add('spread_a', 0, vary=False)
    p_init.add('spread_12', 0.2, min=0, max=2, brute_step=0.1)  # 20

    p_init.add('f_a', -1.9, min=-5, max=-1, vary=False, brute_step=0.1)     # 40
    p_init.add('f_x0', 13.4, vary=False)
    p_init.add('f_b', 0.65, min=0, max=5, vary=False)
    p_init.add('f_cutoff', 10.5, vary=False)

    p_init.add('slope_1', p_init['mms_a1'] / p_init['mms_scale'], expr='mms_a1 / mms_scale', max=5)
    p_init.add('slope_2', p_init['mms_a1'] / p_init['mms_scale'], expr='mms_a2 / mms_scale', max=2)
    p_init.add('spread_low', p_init['spread_a'] * (9-12) + p_init['spread_12'],
               expr='spread_a * (9-12) + spread_12', min=0)
    p_init.add('spread_high', p_init['spread_a'] * (15-12) + p_init['spread_12'],
               expr='spread_a * (15-12) + spread_12', min=0)

    params_names = []
    for p in p_init.keys():
        if not p_init[p].expr:
            params_names.append(p)
    independent_vars = ['mv', 'ms']

    def __init__(self, data, weights=None, func=None, mask=None, theplotter=ThePlot(), **kwargs):
        self.data = data
        self.indep = kwargs
        self.weights = weights

        self.funcargs = {key: kwargs[key] for key in self.independent_vars}
        for p in self.params_names:
            self.funcargs[p] = None

        self.tpc = theplotter
        self.mask = None if mask is None else mask.ravel()
        self.func = self.tpc.theplot if func is None else func

    def get_mv(self):
        return self.indep['mv']

    def eval(self, params):
        for key in self.params_names:
            self.funcargs[key] = params[key].value
        return self.func(**self.funcargs)

    def fit(self, params=None, varyf=True, **kwargs):
        if params is None:
            params = self.p_init.copy()
        if not varyf:
            self.tpc.set_f_mv(mv=self.get_mv(),
                              a=params['f_a'], x0=params['f_x0'],
                              b=params['f_b'], cutoff=params['f_cutoff'])
            params['f_a'].vary = params['f_x0'].vary = params['f_b'].vary = params['f_cutoff'].vary = False
        res = lmfit.minimize(self.residual, params=params, **kwargs)
        self.tpc.f_mv = self.tpc.f_mv_v
        return res

    def chi2(self, params):
        return np.nansum(self.residual(params) ** 2)

    def residual(self, params):
        res = self.eval(params) - self.data
        if self.weights is not None:
            res *= self.weights
        res = res.ravel()
        return res if self.mask is None else res[self.mask]


class ThePlotSigmasModel(ThePlotModel):
    p_init = ThePlotModel.p_init.copy()
    p_init.add('N_a', 0.9, min=0, max=5, brute_step=0.1)  #50
    p_init.add('N_12', 0, min=-2, max=2, brute_step=0.1, vary=False)  #20
    p_init.add('fNo', 0.6, min=0, max=1, brute_step=0.1)  #10

    params_names = []
    for p in p_init.keys():
        if not p_init[p].expr:
            params_names.append(p)
    independent_vars = ['ms']

    def __init__(self, dM=0.1, low=8, high=15, **kwargs):
        super(ThePlotSigmasModel, self).__init__(func=self.func, **kwargs)

        self.mv = (np.arange(low, high, dM) + dM/2)[:, np.newaxis]
        self.s2 = self.tpc.sigma2(self.mv)
        self.y = None

    def get_mv(self):
        return self.mv

    def func(self, ms,
             mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale,
             spread_a, spread_12,
             f_a, f_x0, f_b, f_cutoff,
             N_a, N_12, fNo):
        if self.y is None or self.y.shape[1] != len(ms):
            self.y = np.empty((3, len(ms)))
        self.tpc.observables(self.y, self.mv, self.s2,
                             ms[np.newaxis, :],
                             mms_mv0, mms_ms0, mms_a1, mms_a2, mms_scale,
                             spread_a, spread_12,
                             f_a, f_x0, f_b, f_cutoff,
                             N_a, N_12)
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
        print(resid)
        self.residax.plot(i, resid, '.')
        self.residfig.canvas.draw()

    def plot_plot(self, r):
        if self.fitlines:
            for fl in self.fitlines:
                if fl in fl.axes.lines:
                    fl.axes.lines.remove(fl)
        self.fitlines = self.plotter(r)
        self.plotfig.canvas.draw()

    def __call__(self, params, iter, resid, *args, **kws):
        print('='*10, iter, '='*10)

        self.resids.append(np.sum(resid**2))
        self.params.append(params)

        if self.saveparams:
            with open(self.saveparams, 'w') as f:
               f.write('[' + ', '.join([p.dumps() for p in self.params]) + ']')

        if self.doparams:
            self.plot_params(params, iter)
        if self.doresid:
            self.plot_resid(self.resids[-1], iter)
        self.plot_plot(resid)
        plt.pause(0.01)