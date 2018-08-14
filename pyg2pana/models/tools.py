# Author: Chao Gu, 2018

from scipy import constants

__all__ = ['mass']

_ma = constants.value('atomic mass constant energy equivalent in MeV') / 1000


def mass(z, a):
    return {
        (1, 1): 1.007940,
        (2, 4): 4.002602,
        (6, 12): 12.0107,
        (7, 14): 14.0067,
    }.get((z, a), a) * _ma
