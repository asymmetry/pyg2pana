#!/usr/bin/env python3

import argparse as ap

import matplotlib.pyplot as plt
import numpy

from pyg2pana import Data

cuts = {
    'y': [-0.015, 0.025],
    't': [-0.01, 0.03],
    'p': [-0.015, 0.015],
    'sr': 0.5,
}

parser = ap.ArgumentParser(prog='subtract.py')
parser.add_argument('runs', nargs=2, help='production run and empty run')

args = vars(parser.parse_args())

run_p = int(args['runs'][0])
run_e = int(args['runs'][1])

p = Data('data/g2p_{}.npz'.format(run_p))
e = Data('data/g2p_{}.npz'.format(run_e), ref=run_p)

p.cuts = cuts
e.cuts = cuts

yield_p, _ = numpy.histogram(p.nu[p.cuts], bins=800, range=(-100, 1500))
yield_e, _ = numpy.histogram(e.nu[e.cuts], bins=800, range=(-100, 1500))

yield_p = yield_p * p.scale
yield_e = yield_e * e.scale

yield_sub = yield_p / p.charge - yield_e / e.charge

bin_centers = numpy.linspace(-99, 1499, 800)

plt.figure()
ax = plt.gca()
plt.xlabel(r'$\nu$')
plt.xlim(p.p0 * 0.94, p.p0 * 1.06)
plt.bar(bin_centers, yield_sub, width=2, fill=False)
plt.show()
