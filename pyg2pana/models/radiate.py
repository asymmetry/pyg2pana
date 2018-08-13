# Author: Chao Gu, 2018

import numpy as np
from scipy import constants, integrate, special

from ._tools import mass

__all__ = ['radiate_inelastic_xs']

_alpha = constants.alpha
_alpha_pi = constants.alpha / np.pi
_m_e = constants.value('electron mass energy equivalent in MeV') / 1000


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
    e : rank-1 array of float
        Energy of incident electron.
    ep : rank-1 array of float
        Energy of scattered electron.
    theta : rank-1 array of float
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
    m_t = mass(z, a)

    # scalars
    de = 0.005  # (A83)
    logz13 = np.log(183 * np.power(z, -1 / 3))
    eta = np.log(1440 * np.power(z, -2 / 3)) / logz13  # (A46)
    b = 4 / 3 * (1 + 1 / 9 * ((z + 1) / (z + eta)) / logz13)  # (A45)
    t = tb + ta  # (A47)
    xi = _m_e / (2 * _alpha_pi) * t / ((z + eta) * logz13)  # (A52)

    # vectors
    sin2_theta_2 = np.sin(theta / 2)**2
    q2 = 4 * es * ep * sin2_theta_2
    r = (m_t + 2 * es * sin2_theta_2) / (m_t - 2 * ep * sin2_theta_2)
    tr = _alpha_pi * (np.log(q2 / _m_e**2) - 1) / b  # (A57)

    # in scipy, spence is defined as \int_0^z log(t)/(1-t) dt
    # spence in the reference is spence(1 - z) here
    # so use 1 - cos2_theta_2 = sin2_theta_2
    spence = special.spence(sin2_theta_2)  # (A48)

    # (A44)
    def _f(es, ep, q2, spence):
        logq2me = np.log(q2 / _m_e**2)
        ff = 1 + 0.5772 * b * t
        ff += 2 * _alpha_pi * (-14 / 9 + 13 / 12 * logq2me)
        ff -= _alpha_pi / 2 * (np.log(es / ep))**2
        ff += _alpha_pi * (np.pi**2 / 6 - spence)
        return ff

    # (A82), 1st term
    term1_1 = np.power(r * de / es, b * (tb + tr))
    term1_2 = np.power(de / ep, b * (ta + tr))
    term1_3 = 1 - xi / de / (1 - b * (t + 2 * tr))
    xs = func(z, a, es, ep, theta, *args)
    term1 = term1_1 * term1_2 * term1_3 * _f(es, ep, q2, spence) * xs

    # (A54)
    def _phi(v):
        return 1 - v + 0.75 * v**2

    # (A82), 2nd term, integrand
    def term2_integrand(esp, es, ep, theta, q2, r, tr, spence):
        term2_1 = np.power((es - esp) / (ep * r), b * (ta + tr))
        term2_2 = np.power((es - esp) / es, b * (tb + tr))
        term2_3 = b * (tb + tr) / (es - esp) * _phi((es - esp) / es)
        term2_3 += xi / (2 * (es - esp)**2)
        xs = func(z, a, esp, ep, theta, *args)
        return term2_1 * term2_2 * term2_3 * _f(es, ep, q2, spence) * xs

    # (A82), 3rd term, integrand
    def term3_integrand(epp, es, ep, theta, q2, r, tr, spence):
        term3_1 = np.power((epp - ep) / epp, b * (ta + tr))
        term3_2 = np.power((epp - ep) * r / es, b * (tb + tr))
        term3_3 = b * (ta + tr) / (epp - ep) * _phi((epp - ep) / epp)
        term3_3 += xi / (2 * (epp - ep)**2)
        xs = func(z, a, es, epp, theta, *args)
        return term3_1 * term3_2 * term3_3 * _f(es, ep, q2, spence) * xs

    if np.isscalar(term1):
        term2, _ = integrate.quad(
            term2_integrand,
            ep / (1 - q2 / (2 * es * m_t)), # (A50)
            es - r * de,
            args=(es, ep, theta, q2, r, tr, spence),
            epsrel=1e-3,
        )
        term3, _ = integrate.quad(
            term3_integrand,
            ep + de,
            es / (1 + q2 / (2 * ep * m_t)), # (A51)
            args=(es, ep, theta, q2, r, tr, spence),
            epsrel=1e-3,
        )
    else:
        term2 = np.zeros_like(term1)
        term3 = np.zeros_like(term1)
        it = np.nditer(
            [es, ep, theta, q2, r, tr, spence, term2, term3],
            op_flags=[['readonly'], ['readonly'], ['readonly'], ['readonly'],
                      ['readonly'], ['readonly'], ['readonly'], ['writeonly'],
                      ['writeonly']],
        )
        for ies, iep, itheta, iq2, ir, itr, ispence, iterm2, iterm3 in it:
            iterm2[...], _ = integrate.quad(
                term2_integrand,
                iep / (1 - iq2 / (2 * ies * m_t)), # (A50)
                ies - ir * de,
                args=(ies, iep, itheta, iq2, ir, itr, ispence),
                epsrel=1e-3,
            )
            iterm3[...], _ = integrate.quad(
                term3_integrand,
                iep + de,
                ies / (1 + iq2 / (2 * iep * m_t)), # (A51)
                args=(ies, iep, itheta, iq2, ir, itr, ispence),
                epsrel=1e-3,
            )

    return term1 + term2 + term3
