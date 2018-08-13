#!/usr/bin/env python3


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('pbosted', parent_package, top_path)
    config.add_extension(
        '_pbosted',
        sources=['pbosted.pyf', 'pbosted.f95', 'F1F209.f'],
    )
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup

    setup(**configuration(top_path='').todict())
