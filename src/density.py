import numpy as np
from scipy.optimize import minimize_scalar
from scipy.interpolate import interp2d

import numba


@numba.jit
def kernel(d):
    t2 = np.sum(d*d)
    return (4/np.pi)*(1-t2)**3 if t2<1 else 0


@numba.jit
def find_all_roots(x, y):
    roots = []
    for i in range(len(y)-1):
        if y[i] == 0:
            roots.append(x[i])
        elif (y[i] < 0 < y[i+1]) or (y[i] > 0 > y[i+1]):
            roots.append(x[i] + (x[i+1] - x[i]) * (-y[i]) / (y[i+1] - y[i]))
    return np.array(roots)


class D99_caustics:
    def __init__(self, data):
        self.data = data
        self.xs = None
        self.ys = None

        self.rho = None
        self.A = None

        self.dc = None
        self.cc = None

    def calculate_density(self, xs, ys, xtol=0.1, maxiter=None):
        print('Calculating the density...')
        self.dc = D99_density_calculator(self.data, xs, ys)
        
        options = {'xtol': xtol}
        if maxiter is not None:
            options['maxiter'] = maxiter
        m = minimize_scalar(self.dc.M,
                            bracket=(self.dc.h_c,
                                     0.9*self.dc.h_c),
                            options=options)
        if self.dc.h_c != m.x:
            self.dc.h_c = m.x
            self.rho = self.dc.get_density()
        else:
            self.rho = self.dc.rho

        self.xs = xs
        self.ys = ys

        print('Densisty calculated.')
        return self.rho

    def mock_calculate_density(self, xs, ys, *args, **kwargs):
        self.dc = D99_density_calculator(self.data, xs, ys)
        self.rho = self.dc.get_density()
        self.xs = xs
        self.ys = ys

        return self.rho

    def calculate_caustics(self):
        print('Determining caustics...')
        self.cc = D99_caustic_calculator(self.data, self.rho, self.xs, self.ys)

        m = minimize_scalar(self.cc.S, [0.5, 1])
        if self.cc.k != m.x:
            self.cc.k = m.x
            self.A = self.cc.get_A()
        else:
            print(self.cc.A)
            self.A = self.cc.A

        print('Caustics calculated.')
        return self.A


@numba.jitclass({
    'data': numba.typeof(np.array([[1.0]])),
    'N': numba.int_,
    'xs': numba.typeof(np.array([1.0])),
    'ys': numba.typeof(np.array([1.0])),
    'nx': numba.int_,
    'ny': numba.int_,
    'dx': numba.float64,
    'dy': numba.float64,
    'h_c': numba.float64,
    'h_opt': numba.float64,
    'r_h_opt': numba.float64,
    'l': numba.typeof(np.array([1.0])),
    'rho': numba.typeof(np.array([[1.0]]))
})
class D99_density_calculator:
    def __init__(self, data, xs, ys):
        self.data = data
        self.N = len(data)

        self.xs, self.ys = xs, ys
        self.nx = len(xs)
        self.ny = len(ys)
        self.dx = xs[1] - xs[0]
        self.dy = ys[1] - ys[0]

        self.h_c = 1

        # self.h_opt = (3.12 / self.N**(1/6)) * np.sqrt(self.data.var(axis=0).sum() / 2)
        self.h_opt = (3.12 / self.N**(1/6)) * np.sqrt((self.data[:, 0].var() + self.data[:, 1].var()) / 2)
        self.r_h_opt = 1 / self.h_opt
        f1 = np.array([self._density(self.data[i]) for i in range(len(self.data))])
        self.l = self.h_opt * np.sqrt(np.exp(np.sum(np.log(f1)) / self.N) / f1)

    def _density(self, x0):
        s = 0
        for i in range(len(self.data)):
            s += kernel((x0 - self.data[i]) * self.r_h_opt)
        return s * self.r_h_opt*self.r_h_opt / self.N
    def _adaptive_density(self, x0, skip=-1):
        s = 0
        for i in range(self.N):
            if i == skip: continue
            h = self.h_c * self.l[i]
            s += kernel((x0 - self.data[i]) / h) / h ** 2
        return s / (self.N if skip < 0 else self.N - 1)

    def get_density(self):
        s = np.empty((self.ny, self.nx))
        for i in range(self.nx):
            for j in range(self.ny):
                c = np.array((self.xs[i], self.ys[j]))
                s[j, i] = self._adaptive_density(c, -1)
        return s

    def M(self, h_c):
        # print('Trying h_c='+str(h_c))

        self.h_c = h_c
        self.rho = self.get_density()
        S = (self.rho * self.rho).sum() * self.dx * self.dy

        s = 0
        for i in range(self.N):
            s += self._adaptive_density(self.data[i], i)

        return S - (2 / self.N) * s


class D99_caustic_calculator:
    def __init__(self, data, rho, xs, ys):
        self.data = data

        self.dr = (xs[1] - xs[0]) / 10
        self.dv = (ys[1] - ys[0]) / 10
        self.rs = np.arange(xs[0], xs[-1], self.dr)
        self.vs = np.arange(ys[0], ys[-1], self.dv)

        self.rho = interp2d(xs, ys, rho, 'cubic')(self.rs, self.vs)
        self.phi = self.rho.sum(axis=0) * self.dv

        self.R = self.data[:, 0].mean()
        self.four_avg_v2 = 4 * (self.data[:, 1]**2).mean()
        self.maxi = int((self.R - self.rs[0]) / self.dr) + 1

        self.phis = self.phi[:self.maxi]
        self.sum_phis = self.phis.sum()

        self.k = None
        self.A = np.ndarray([len(self.rs)])
        self.caustics = np.ndarray((2, len(self.rs)))

    def get_A(self):
        rho_k = self.rho - self.k
        for i in range(len(self.A)):
            roots = find_all_roots(self.vs, rho_k[:, i])

            if len(roots):
                maxloc = self.vs[np.argmax(rho_k[:, i])]
                pos_root = roots[np.where((roots > 0) * (roots > maxloc))]
                pos_root = np.min(pos_root) if len(pos_root) else np.inf
                neg_root = roots[np.where((roots < 0) * (roots < maxloc))]
                neg_root = np.min(neg_root) if len(neg_root) else -np.inf
            else:
                pos_root = neg_root = 0

            self.caustics[0, i] = pos_root
            self.caustics[1, i] = neg_root

        self.A = np.min(np.abs(self.caustics), axis=0)

        return self.A


    def S(self, k):
        print('Trying k='+str(k))

        self.k = k

        self.get_A()
        avg_esc2 = (self.A[:self.maxi]**2 * self.phis).sum() / self.sum_phis

        return (avg_esc2 - self.four_avg_v2)**2


if __name__=="__main__":
    from matplotlib import pyplot as plt

    N = 500
    coords = np.ndarray((N, 2))
    coords[:, 0] = np.abs(np.random.normal(0, 0.5, size=N))
    coords[:, 1] = np.random.normal(0, 0.1, size=N)

    GRIDSIZE = 0.1
    BOUNDS = (0, 1.1, -0.5, 0.6)

    xs = np.arange(BOUNDS[0], BOUNDS[1], GRIDSIZE)
    ys = np.arange(BOUNDS[2], BOUNDS[3], GRIDSIZE)

    d = D99_caustics(coords)
    d.calculate_density(xs, ys)
    d.calculate_caustics()

    plt.xlim((BOUNDS[0], BOUNDS[1]-GRIDSIZE))
    plt.ylim((BOUNDS[2], BOUNDS[3]-GRIDSIZE))
    plt.plot(coords[:, 0], coords[:, 1], 'rx')
    plt.imshow(d.cc.rho, cmap='inferno', extent=(BOUNDS[0], BOUNDS[1]-GRIDSIZE, BOUNDS[3]-GRIDSIZE, BOUNDS[2]), aspect='auto')
    plt.colorbar()

    plt.plot(d.cc.rs, d.cc.caustics[0], 'w:')
    plt.plot(d.cc.rs, d.cc.caustics[1], 'w:')
    plt.plot(d.cc.rs, d.cc.A, 'w-')
    plt.plot(d.cc.rs, -d.cc.A, 'w-')
