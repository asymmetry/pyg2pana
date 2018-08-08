#!/usr/bin/env python3

import numpy

from . import _pbosted

__all__ = ['pbosted']


def pbosted(z, a, e, ep, theta):
    after_broadcast = numpy.broadcast_arrays(z, a, e, ep, theta)
    result = _pbosted.cal_xs(*after_broadcast)
    if result.shape == (1, ) and isinstance(z, (int, float)):
        return result[0]

    return result
