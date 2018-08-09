# Author: Chao Gu, 2018

import numpy as np

from . import _pbosted
from ..radiate import radiate_inelastic_xs

__all__ = ['xs', 'xs_rad']


def xs(z, a, e, ep, theta):
    """
    Return non-radiated inelastic cross section.

    Parameters
    ----------
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
    """

    after_broadcast = np.broadcast_arrays(e, ep, theta)
    result = _pbosted.cal_xs(z, a, *after_broadcast)
    if result.shape == (1, ) and isinstance(e, (int, float)):
        return result[0]
    return result


def xs_rad(z, a, e, ep, theta, tb, ta):
    """
    Return radiated inelastic cross section.

    Parameters
    ----------
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
    """

    return radiate_inelastic_xs(xs, z, a, e, ep, theta, tb, ta)
