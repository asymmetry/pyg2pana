# Author: Chao Gu, 2018

# S. Venkat et al., Phys. Rev. C, 83(2011)015203

_m = 0.938783


def ge(q2):
    t = q2 / (4 * _m**2)
    result = 1 + 2.90966 * t - 1.11542229 * t**2 + 3.866171e-2 * t**3
    result /= (1 + 14.5187212 * t + 40.88333 * t**2 + 99.999998 * t**3 +
               4.579e-5 * t**4 + 10.3580447 * t**5)
    return result


def gm(q2):
    t = q2 / (4 * _m**2)
    result = 1 - 1.43573 * t + 1.19052066 * t**2 + 2.5455841e-1 * t**3
    result /= (1 + 9.70703681 * t + 3.7357e-4 * t**2 + 6.0e-8 * t**3 +
               9.9527277 * t**4 + 12.7977739 * t**5)
    return 2.792782 * result
