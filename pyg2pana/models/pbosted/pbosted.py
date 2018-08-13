# Author: Chao Gu, 2018

import numpy as np

from . import _pbosted
from ..radiate import radiate_inelastic_xs as rad

__all__ = ['PBosted']


class PBosted():
    """
    Calculate inelastic cross section using Peter Bosted's fit result for a
    particular nucleus.

    Parameters
    ----------
    z : int
        Atomic number.
    a : int
        Mass number.
    radiate : bool
        Calculate non-radiated cross section if False (default), return
        radiated cross section if True.
    """

    def __init__(self, z, a, radiate=False):
        self.z = z
        self.a = a
        self.radiate = radiate

    def __call__(self, e, ep, theta, tb=0, ta=0):
        """
        Calculate inelastic cross section for a particular nucleus.

        Parameters
        ----------
        e : float or rank-1 array of float
            Energy of incident electron.
        ep : float or rank-1 array of float
            Energy of scattered electron.
        theta : float or rank-1 array of float
            Scattering angle.
        tb : float
            Radiation length before scattering.
        ta : float
            Radiation length after scattering.
        """

        if any(not np.isscalar(x) for x in (e, ep, theta)):
            e, ep, theta = np.broadcast_arrays(e, ep, theta)

        if self.radiate:
            result = rad(self._xs, self.z, self.a, e, ep, theta, tb, ta)
        else:
            result = self._xs(self.z, self.a, e, ep, theta)

        return result

    def _xs(self, z, a, e, ep, theta):
        if any(not np.isscalar(x) for x in (e, ep, theta)):
            return _pbosted.cal_xs_array(z, a, e, ep, theta)
        return _pbosted.cal_xs_scalar(z, a, e, ep, theta)
