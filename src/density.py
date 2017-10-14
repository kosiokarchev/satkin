import numpy as np
from scipy.optimize import minimize_scalar


def kernel(d):
    t2 = np.sum(d*d)
    return (4/np.pi)*(1-t2)**3 if t2<1 else 0


class D99_caustics:
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
        print('Trying h_c=' + str(h_c))

        self.h_c = h_c
        self.rho = self.get_density()
        S = (self.rho * self.rho).sum() * self.dx * self.dy

        plt.imshow(self.rho)
        plt.pause(0.05)

        sum = 0
        for i in range(self.N):
            sum += self._adaptive_density(self.data[i], i)

        return S - (2 / self.N) * sum

    def calculate_density(self, xtol=0.1, maxiter=None):
        options = {'xtol': xtol}
        if maxiter is not None:
            options['maxiter'] = maxiter
        m = minimize_scalar(self.M,
                            bracket=(self.h_c, 0.9*self.h_c),
                            options=options)
        if self.h_c != m.x:
            self.h_c = m.x
            self.rho = self.get_density()
        return self.rho


if __name__=="__main__":
    from matplotlib import pyplot as plt

    N = 300
    R = np.random.rayleigh(0.5, size=N)
    PHI = np.random.rayleigh(2, size=N)
    coords = np.ndarray((N, 2))
    coords[:, 0] = R*np.cos(PHI)
    coords[:, 1] = R*np.sin(PHI)

    xs = np.arange(-1, 1.1, 0.05)
    ys = np.arange(-1, 1.1, 0.05)
    d = D99_caustics(coords, xs, ys)

    plt.plot(20 + 20 * coords[:, 0], 20 + 20 * coords[:, 1], 'rx')
    d.calculate_density()
    plt.imshow(d.rho)
