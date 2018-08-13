# Author: Chao Gu, 2018

from functools import partial

import numpy as np
from scipy import constants, integrate

from .._tools import mass

__all__ = ['Elastic']

_alpha = constants.alpha
_me = constants.value('electron mass energy equivalent in MeV') / 1000
_m_p = constants.value('proton mass energy equivalent in MeV') * 1e-3
_inv_fm_to_gev = constants.hbar * constants.c / constants.e * 1e6
_gev_to_inv_fm = 1 / _inv_fm_to_gev
_inv_gev_to_fm = _inv_fm_to_gev
_inv_gev_to_mkb = _inv_gev_to_fm**2 * 1e4  # GeV^{-2} to microbarn

_res = 1e-3 / 3
_sqrt_2pi = np.sqrt(2 * np.pi)

_rho_limit = 20


class Elastic():
    """
    Calculate elastic cross section for a particular nucleus.

    Parameters
    ----------
    z : int
        Atomic number.
    a : int
        Mass number.
    """

    def __init__(self, z, a):
        self.z = z
        self.a = a
        self.m = mass(z, a)

        if (z, a) == (1, 1):
            from ._proton import ge as ge_proton, gm as gm_proton
            self.ff_func = partial(
                self._ff_sachs,
                ge_func=ge_proton,
                gm_func=gm_proton,
            )
        elif (z, a) in ((2, 4), (6, 12), (7, 14)):
            from ._density import get_density_func
            charge_density = get_density_func('charge', z, a)
            magnet_density = get_density_func('magnetization', z, a)
            if charge_density is not None:
                self.charge_density_0, _ = integrate.quad(
                    lambda x: charge_density(x) * x**2,
                    *(0, _rho_limit),
                )
            if magnet_density is not None:
                self.magnet_density_0, _ = integrate.quad(
                    lambda x: magnet_density(x) * x**2,
                    *(0, _rho_limit),
                )
            self.ff_func = partial(
                self._ff_density,
                charge_density_func=charge_density,
                magnet_density_func=magnet_density,
            )
        else:
            self.ff_func = self._ff

    def __call__(self, e, theta):
        """
        Calculate inelastic cross section for a particular nucleus.

        Parameters
        ----------
        e : float or rank-1 array of float
            Energy of incident electron.
        theta : float or rank-1 array of float
            Scattering angle.
        """

        if any(not np.isscalar(x) for x in (e, theta)):
            e, theta = np.broadcast_arrays(e, theta)

        result = self._xs(self.z, self.a, e, theta)
        return result

    def _xs(self, z, _, e, theta):
        sin2_theta_2 = np.sin(theta / 2)**2
        cos2_theta_2 = 1 - sin2_theta_2

        recoil = 1 / (1 + 2 * e / self.m * sin2_theta_2)
        el = e * recoil
        q2 = 4.0 * e * el * sin2_theta_2

        mott = (z * _alpha / (2 * e * sin2_theta_2))**2 * cos2_theta_2
        ff = self.ff_func(e, q2)

        return mott * recoil * ff * _inv_gev_to_mkb

    def _ff(self, _, q2):
        x_alpha = (self.z - 2) / 3
        q2_fm = q2 * _gev_to_inv_fm * _gev_to_inv_fm

        if self.z == 6:
            if q2_fm < 3.2:
                xa = 1.64
            elif q2_fm > 3.5:  # second diffraction minimum
                xa = 1.68
            else:
                xa = 0
        else:
            xa = 1.64

        result = 1 - x_alpha / (2 * (2 + 3 * x_alpha)) * q2_fm * xa**2
        result *= np.exp(-(q2_fm * xa**2) / 4)

        if result < 1.0e-6:
            return 1.0e-6

        return result

    def _ff_sachs(self, e, q2, ge_func, gm_func):
        el = e - q2 / (2 * self.m)
        sin2_theta_2 = q2 / (4 * e * el)
        cos2_theta_2 = 1 - sin2_theta_2
        tan2_theta_2 = sin2_theta_2 / cos2_theta_2
        ge = ge_func(q2)
        gm = gm_func(q2)
        tau = q2 / (4 * self.m**2)
        epsilon = 1 / (1 + 2 * (1 + tau) * tan2_theta_2)
        return (epsilon * ge**2 + tau * gm**2) / (epsilon * (1 + tau))

    def _ge_integrand(self, r, q, charge_density_func):
        result = charge_density_func(r)
        result *= r * np.sin(r * q) / (self.charge_density_0 * q)
        return result

    def _gm_integrand(self, r, q, magnet_density_func):
        result = magnet_density_func(r)
        result *= r * np.sin(r * q) / (self.magnet_density_0 * q)
        return result

    def _ff_density(self, e, q2, charge_density_func, magnet_density_func):
        q = np.sqrt(q2)
        q_fm = q * _gev_to_inv_fm

        if np.isscalar(q_fm):
            ge, gm = 0, 0
            if charge_density_func is not None:
                ge, _ = integrate.quad(
                    self._ge_integrand,
                    *(0, _rho_limit),
                    args=(q_fm, charge_density_func),
                )
            if magnet_density_func is not None:
                gm, _ = integrate.quad(
                    self._gm_integrand,
                    *(0, _rho_limit),
                    args=(q_fm, magnet_density_func),
                )
        else:
            ge = np.zeros_like(q_fm)
            gm = np.zeros_like(q_fm)
            if charge_density_func is not None:
                it = np.nditer(
                    [q_fm, ge],
                    op_flags=[['readonly'], ['writeonly']],
                )
                for iq_fm, ige in it:
                    ige[...], _ = integrate.quad(
                        self._ge_integrand,
                        *(0, _rho_limit),
                        args=(iq_fm, charge_density_func),
                    )
            if magnet_density_func is not None:
                it = np.nditer(
                    [q_fm, gm],
                    op_flags=[['readonly'], ['writeonly']],
                )
                for iq_fm, igm in it:
                    igm[...], _ = integrate.quad(
                        self._gm_integrand,
                        *(0, _rho_limit),
                        args=(iq_fm, magnet_density_func),
                    )

        el = e - q2 / (2 * self.m)
        sin2_theta_2 = q2 / (4 * e * el)
        cos2_theta_2 = 1 - sin2_theta_2
        tan2_theta_2 = sin2_theta_2 / cos2_theta_2
        tau = q2 / (4 * self.m**2)
        epsilon = 1 / (1 + 2 * (1 + tau) * tan2_theta_2)
        return (epsilon * ge**2 + tau * gm**2) / (epsilon * (1 + tau))
