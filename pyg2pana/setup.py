#!/usr/bin/env python3


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('pyg2pana', parent_package, top_path)
    config.add_subpackage('models')
    config.add_data_files('g2p.db')
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup

    setup(**configuration(top_path='').todict())
