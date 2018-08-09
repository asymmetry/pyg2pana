# Author: Chao Gu, 2018

import numpy as np
from scipy import constants, integrate, special

__all__ = ['radiate_inelastic_xs']

_alpha = constants.alpha
_alpha_pi = constants.alpha / np.pi
_ma = constants.value('atomic mass constant energy equivalent in MeV') / 1000
_me = constants.value('electron mass energy equivalent in MeV') / 1000

_za_mass = {
    (1, 1): 1.007940,
    (2, 4): 4.002602,
    (6, 12): 12.0107,
    (7, 14): 14.0067,
}


def radiate_inelastic_xs(func, z, a, e, ep, theta, tb, ta, *, args=()):
    """
    Return radiated inelastic cross section.

    Parameters
    ----------
    func : callable
        Non-radiated cross section function.
    z : int
        Atomic number.
    a : int
        Mass number.
    e : float or rank-1 array
        Energy of incident electron.
    ep : float or rank-1 array
        Energy of scattered electron.
    theta : float or rank-1 array
        Scattering angle.
    tb : float
        Radiation length before scattering.
    ta : float
        Radiation length after scattering.
    args : tuple, optional
        Extra arguments to pass to function, if any.

    References
    ----------
    S. Stein et al., Phys. Rev. D 12(1975)1884
    """

    es = e
    es, ep, theta = np.broadcast_arrays(es, ep, theta)

    # scalars
    de = 0.01  # (A83)
    m_t = _za_mass.get((z, a), a) * _ma
    logz13 = np.log(183 * np.power(z, -1 / 3))
    eta = np.log(1440 * np.power(z, -2 / 3)) / logz13  # (A46)
    b = 4 / 3 * (1 + 1 / 9 * ((z + 1) / (z + eta)) / logz13)  # (A45)
    t = tb + ta  # (A47)
    xi = _me / (2 * _alpha_pi) * t / ((z + eta) * logz13)  # (A52)

    # vectors
    sin2_theta_2 = np.sin(theta / 2)**2
    q2 = 4 * es * ep * sin2_theta_2
    r = (m_t + 2 * es * sin2_theta_2) / (m_t - 2 * ep * sin2_theta_2)

    # in scipy, spence is defined as \int_0^z log(t)/(1-t) dt
    # spence in the reference is spence(1 - z) here
    # so use 1 - cos2_theta_2 = sin2_theta_2
    spence = special.spence(sin2_theta_2)  # (A48)

    # (A57)
    def _tr(q2):
        logq2me = np.log(q2 / _me**2)
        return _alpha_pi * (logq2me - 1) / b

    # (A44)
    def _f(es, ep, q2, spence):
        logq2me = np.log(q2 / _me**2)
        ff = 1 + 0.5772 * b * t
        ff += 2 * _alpha_pi * (-14 / 9 + 13 / 12 * logq2me)
        ff -= _alpha_pi / 2 * (np.log(es / ep))**2
        ff += _alpha_pi * (np.pi**2 / 6 - spence)
        return ff

    # (A82), 1st term
    tr_q2 = _tr(q2)
    term1_1 = np.power(r * de / es, b * (tb + tr_q2))
    term1_2 = np.power(de / ep, b * (ta + tr_q2))
    term1_3 = 1 - xi / de / (1 - b * (t + 2 * tr_q2))
    xs = func(z, a, es, ep, theta, *args)
    term1 = term1_1 * term1_2 * term1_3 * _f(es, ep, q2, spence) * xs

    # (A54)
    def _phi(v):
        return 1 - v + 0.75 * v**2

    # (A82), 2nd term, integrand
    def term2_integrand(esp, es, ep, theta, q2, r, spence):
        tr_q2 = _tr(q2)
        term2_1 = np.power((es - esp) / (ep * r), b * (ta + tr_q2))
        term2_2 = np.power((es - esp) / es, b * (tb + tr_q2))
        term2_3 = b * (tb + tr_q2) / (es - esp) * _phi((es - esp) / es)
        term2_3 += xi / (2 * (es - esp)**2)
        xs = func(z, a, esp, ep, theta, *args)
        return term2_1 * term2_2 * term2_3 * _f(es, ep, q2, spence) * xs

    # (A82), 3rd term, integrand
    def term3_integrand(epp, es, ep, theta, q2, r, spence):
        tr_q2 = _tr(q2)
        term3_1 = np.power((epp - ep) / epp, b * (ta + tr_q2))
        term3_2 = np.power((epp - ep) * r / es, b * (tb + tr_q2))
        term3_3 = b * (ta + tr_q2) / (epp - ep) * _phi((epp - ep) / epp)
        term3_3 += xi / (2 * (epp - ep)**2)
        xs = func(z, a, es, epp, theta, *args)
        return term3_1 * term3_2 * term3_3 * _f(es, ep, q2, spence) * xs

    term2 = np.zeros_like(term1)
    term3 = np.zeros_like(term1)
    for i, arg_set in enumerate(zip(es, ep, theta, q2, r, spence)):
        term2[i], _ = integrate.quad(
            term2_integrand,
            ep[i] / (1 - 2 * ep[i] * sin2_theta_2[i] / m_t), # (A50)
            es[i] - r[i] * de,
            args=arg_set,
        )
        term3[i], _ = integrate.quad(
            term3_integrand,
            ep[i] + de,
            es[i] / (1 + 2 * es[i] * sin2_theta_2[i] / m_t), # (A51)
            args=arg_set,
        )

    return term1 + term2 + term3
