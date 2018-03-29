import numpy as np
from numpy import abs, sqrt, arccos as acos, arccosh as acosh, log as ln, pi, inf
from scipy.special import spence
from scipy.integrate import quad as integrate

__all__ = ['lokas']

pi2 = pi**2
Li = lambda x: spence(1-x)

sa = 4/3
sa2 = sa**2


g = lambda c: 1 / (ln(1+c) - c/(1+c))

C1 = lambda Rt, c: acos if Rt>1/c else acosh

beta_om = lambda s: s**2 / (s**2 + sa2)
def sigmar2_overvvir2_om(s, c):
    c2 = c**2
    cs = c*s
    ln1cs = ln(1+cs)

    f1 = g(c) * s * (1+cs)**2 / (2 * (s**2 + sa2))
    f2  = -c * sa2 / s
    f2 += -c2 * sa2 * ln(cs)
    f2 +=  c2 * sa2 * ln1cs * (1 + 1/cs**2 - 4/cs)
    f2 += -(1 + c2*sa2) * (1/(1+cs)**2 + 2*ln1cs/(1+cs))
    f2 +=  (1 + 3*c2*sa2) * (pi2/3 - 2/(1+cs) + ln1cs**2 + 2*Li(-cs))

    return f1*f2

beta0 = 0.25
beta_const = lambda s: beta0
def sigmar2_overvvir2_const(s, c):
    def integrand(x):
        cx1 = (1+c*x)
        return x**(2*beta0 - 3) * ln(cx1) / cx1**2 \
               - c*x**(2*beta0 - 2) / cx1**3

    return g(c) * (1 + c*s)**2 * s**(1 - 2*beta0) * integrate(integrand, s, inf)[0]

beta = beta_const
sigmar2_overvvir2 = sigmar2_overvvir2_const

def S2_overmvirvvir2gc(Rt, c):
    I1 = lambda s: sigmar2_overvvir2(s, c) \
                   * s / (1 + c*s)**2 \
                   * (1 - (2/3)*beta(s))
    I2 = lambda s: sigmar2_overvvir2(s, c) \
                   * sqrt(s**2 - Rt**2) / (1 + c*s)**2 \
                   * (beta(s) * (Rt**2 + 2*s**2) / (3*s**2) - 1)
    return c**2 * (integrate(I1, 0, inf)[0] + integrate(I2, Rt, inf)[0])

def Mp_overmvirgc(Rt, c):
    cRt = c*Rt
    return C1(Rt, c)(1/cRt) / sqrt(abs(cRt**2 - 1)) + ln(cRt/2)

def sigmaap2_overvvir2(Rt, c):
    return S2_overmvirvvir2gc(Rt, c) / Mp_overmvirgc(Rt, c)

lokas = np.vectorize(sigmaap2_overvvir2)
