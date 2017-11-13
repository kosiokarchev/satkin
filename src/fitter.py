from lib import *


def broken_powerlaw(x, a1=1, k1=1, a2=1, k2=1, x0=1):
    return np.where(x<x0, a1 * x**k1, a2 * x**k2)


if __name__ == '__main__':
    # sn = 34
    # t = Table.read(FILES['hw-rms'](sn))
    # t = deal(t)
    #
    # print('Fitting linearly')
    # indices_l, fit_l = fit_or(models.LinearModel(), np.log10(t['rmsx']), x=t['mv'])
    #
    # plt.figure()
    # plt.plot(t['mv'], t['rmsx'], ',')
    # plt.plot(t['mv'][indices_l], t['rmsx'][indices_l], ',')
    # plt.plot(t['mv'][indices_l], np.power(10, fit_l.best_fit))
    # plt.semilogy()
    #
    #
    # print('Fitting powerfully')
    # indices_p, fit_p = fit_or(models.PowerLawModel(), t['rmsx'], x=np.power(10, t['mv']))
    #
    # plt.figure()
    # plt.plot(t['mv'], t['rmsx'], ',')
    # plt.plot(t['mv'][indices_p], t['rmsx'][indices_p], ',')
    # plt.plot(t['mv'][indices_p], fit_p.best_fit)
    # plt.semilogy(t['mv'][indices_l], np.power(10, fit_l.best_fit))


    p = HWProcedure(34, False)



# def do(sn):
#     fname = 'data/sigmas-{}{}.fits'
#     t = Table.read(fname.format('mvir', sn))
#
#     t['mass'] = np.power(10, t['mass'])
#
#     indices, fit = fit_or(models.PowerLawModel(), t['sigmax'], x=t['mass'])
#
#     t.meta['exp'] = fit.params['exponent'].value
#     t.meta['err_exp'] = fit.params['exponent'].stderr
#     t.meta['amp'] = fit.params['amplitude'].value
#     t.meta['err_amp'] = fit.params['amplitude'].stderr
#     t.meta['indices'] = list(indices)
#
#     t.write('data/reg-mvir{}.csv'.format(sn), format='ascii.ecsv', overwrite=True)
#
#     return t


# if __name__ == '__main__':
#     for sn in (34, 38, 45):
#         print(sn)
#         t = do(sn)
#         print('exp: {0[exp]:.4f} ± {0[err_exp]}'.format(t.meta))
#         print('exp: {0[amp]:.4f} ± {0[err_amp]}'.format(t.meta))



# plt.plot(x, y, 'r.')
# plt.plot(x[i_p], y[i_p], 'gx')
# plt.plot(x[i_p], fit_p.best_fit, 'g--', label='power law')
# plt.plot(x[mask_l], y[mask_l], 'b+')
# plt.plot(x[i_l], np.exp(fit_l.best_fit), 'b--', label='linear in logarithms')
# plt.loglog()
#
# plt.legend()