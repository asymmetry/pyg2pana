# Author: Chao Gu, 2018

# C. W. De Jager et al., At. Data Nucl. Data Tables, 14(1974)479

from functools import partial
import numpy as np

_pars_charge = {
    (2, 4): (
        ('3pF', (0.964, 0.322, 0.517), (0.59, 2.5)),
        ('3pF', (1.008, 0.327, 0.445), (0.7, 4.47)),
    ),
    (6, 12): (
        ('3pF', (2.355, 0.5224, -0.149), (0.25, 2.3)),
        ('MHO', (1.649, 1.247), (1.05, 4.01)),
        ('HO', (1.687, 1.067), (0.18, 0.70)),
        ('MHO', (1.672, 1.150), (1.04, 2.15)),
        ('HO', (1.692, 1.082), (0.29, 0.48)),
    ),
    (7, 14): (
        ('HO', (1.76, 1.234), (0.86, 1.62)),
        ('HO', (1.729, 1.291), (0.29, 0.46)),
    ),
}

_pars_magnet = {
    (2, 4): (),
    (6, 12): (),
    (7, 14): (('HO', 0.404, (1.61, 2.48, 3.36), (1.00, 1.80)))
}


def _3pf(r, c, z, w):
    return (1 + w * (r / c)**2) / (1 + np.exp((r - c) / z))


def _ho(r, c, z):
    return (1 + z * (r / c)**2) * np.exp(-(r / c)**2)


def get_density_func(par_set, z, a, par_id=0):
    if par_set == 'charge':
        pars = _pars_charge.get((z, a), ())
    elif par_set == 'magnetization':
        pars = _pars_magnet.get((z, a), ())

    if not pars:
        return None

    pars = pars[par_id]

    if pars[0] == '3pF':
        return partial(_3pf, c=pars[1][0], z=pars[1][1], w=pars[1][2])

    if pars[0] in ('HO', 'MHO'):
        return partial(_ho, c=pars[1][0], z=pars[1][1])
