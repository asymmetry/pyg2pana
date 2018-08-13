#!/usr/bin/env python3


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(
        ignore_setup_xxx_py=True,
        assume_default_configuration=True,
        delegate_options_to_subpackages=True,
        quiet=True,
    )
    config.add_subpackage('pyg2pana')
    config.get_version('pyg2pana/version.py')
    return config


def setup_package():
    try:
        import numpy
    except ImportError:
        build_requires = ['numpy>=1.10.0']
    else:
        build_requires = []

    metadata = dict(
        name='pyg2pana',
        author='Chao Gu',
        author_email='guchao.pku@gmail.com',
        maintainer='Chao Gu',
        maintainer_email='guchao.pku@gmail.com',
        description=
        'Analysis software for Jefferson Lab Experiment E08-027 (g2p).',
        license='GPL-3.0',
        url='https://github.com/asymmetry/pyg2pana',
        classifiers=[
            'Development Status :: 1 - Planning',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Programming Language :: Fortran',
            'Programming Language :: Python',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Operating System :: MacOS',
            'Operating System :: POSIX',
            'Operating System :: Unix',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Physics',
            'Topic :: Utilities',
        ],
        platforms='Any',
        python_requires='>=3.4',
        setup_requires=build_requires,
        install_requires=build_requires,
    )

    metadata['configuration'] = configuration

    from setuptools import setup
    from numpy.distutils.core import setup
    setup(**metadata)


if __name__ == '__main__':
    setup_package()
