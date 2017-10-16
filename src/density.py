import numpy as np
from scipy.optimize import minimize_scalar
from scipy.interpolate import interp2d


def kernel(d):
    t2 = np.sum(d*d)
    return (4/np.pi)*(1-t2)**3 if t2<1 else 0


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

        return self.rho

    def mock_calculate_density(self, xs, ys, *args, **kwargs):
        self.dc = D99_density_calculator(self.data, xs, ys)
        self.rho = self.dc.get_density()
        self.xs = xs
        self.ys = ys

        return self.rho

    def calculate_caustics(self):
        self.cc = D99_caustic_calculator(self.data, self.rho, self.xs, self.ys)

        m = minimize_scalar(self.cc.S, [0.5, 1])
        if self.cc.k != m.x:
            self.cc.k = m.x
            self.A = self.cc.get_A()
        else:
            print(self.cc.A)
            self.A = self.cc.A

        return self.A


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
        self.h_opt = (3.12 / self.N**(1/6)) * np.sqrt(np.sum(self.data.var(0)) / 2)
        self.r_h_opt = 1 / self.h_opt
        f1 = np.array([self._density(x0) for x0 in self.data])
        self.l = self.h_opt * np.sqrt(np.exp(np.sum(np.log(f1)) / self.N) / f1)

        self.rho = None

    def _density(self, x0):
        s = 0
        for p in self.data:
            s += kernel((x0 - p) * self.r_h_opt)
        return s * self.r_h_opt*self.r_h_opt / self.N
    def _adaptive_density(self, x0, skip=-1):
        s = 0
        for i in range(self.N):
            if i == skip: continue
            h = self.h_c * self.l[i]
            s += kernel((x0 - self.data[i]) / h) / h ** 2
        return s / (self.N if skip < 0 else self.N - 1)

    def get_density(self):
        s = np.ndarray((self.nx, self.ny))
        for i in range(self.nx):
            for j in range(self.ny):
                c = np.array((self.xs[i], self.ys[j]))
                s[i, j] = self._adaptive_density(c)
        return s.transpose()

    def M(self, h_c):
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
        self.k = k
        self.get_A()

        avg_esc2 = (self.A[:self.maxi]**2 * self.phis).sum() / self.sum_phis

        return (avg_esc2 - self.four_avg_v2)**2


if __name__=="__main__":
    from matplotlib import pyplot as plt

    N = 300
    coords = np.ndarray((N, 2))
    coords[:, 0] = np.abs(np.random.normal(0, 0.5, size=N))
    coords[:, 1] = np.random.normal(0, 0.1, size=N)

    GRIDSIZE = 0.1
    BOUNDS = (0, 1.1, -0.5, 0.6)

    xs = np.arange(BOUNDS[0], BOUNDS[1], GRIDSIZE)
    ys = np.arange(BOUNDS[2], BOUNDS[3], GRIDSIZE)

    d = D99_caustics(coords)
    d.mock_calculate_density(xs, ys)
    d.calculate_caustics()

    plt.xlim((BOUNDS[0], BOUNDS[1]-GRIDSIZE))
    plt.ylim((BOUNDS[2], BOUNDS[3]-GRIDSIZE))
    plt.plot(coords[:, 0], coords[:, 1], 'rx')
    plt.imshow(d.cc.rho, extent=(BOUNDS[0], BOUNDS[1]-GRIDSIZE, BOUNDS[3]-GRIDSIZE, BOUNDS[2]), aspect='auto')
    plt.colorbar()

    plt.plot(d.cc.rs[:len(d.cc.A)], d.cc.caustics[0], 'w-')
    plt.plot(d.cc.rs[:len(d.cc.A)], d.cc.caustics[1], 'w-')
