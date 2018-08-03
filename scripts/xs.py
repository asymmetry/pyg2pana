#!/usr/bin/env python3

from os.path import join

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy

from pyg2pana import Data, SimFile, configs

e0 = 2253.5
run_list = configs.l_22545000
correction = configs.corrections['l_22545000']

cuts = {
    'y': [-0.015, 0.025],
    't': [-0.01, 0.03],
    'p': [-0.015, 0.015],
    'sr': 0.5,
}

binning = {
    'bins': 1500,
    'range': (-100, 1400),
}

avogadro = 6.023e23  # 1/mol
charge_e = 1.602176487e-13  # uC
l_target = 2.83  # cm
density_target = 0.817  # g/cm^3
a_target = 17.031  # g/mol
factor = 1.0e33  # from cm^2/(sr*MeV) to ub/(sr*GeV)
dilution = 0.17


def get_weight(x, corr):
    return sum(x**i * par for i, par in enumerate(corr[::-1]))


p0_list = []
xs_list = []
exs_list = []
for key, value in run_list.items():
    p0 = float(key)
    p0_list.append(p0)
    print(p0)

    run = value['production'][0]

    data = Data(join('data', 'g2p_{}.npz'.format(run)))
    sim = SimFile(join('sim', 'sim_{}.npz'.format(run)))

    data.cuts = cuts
    sim.cuts = cuts

    acceptance = sim.get_acceptance('nu', **binning)

    w = get_weight(data.rec.d[data.cuts], correction)

    xs, _ = numpy.histogram(data.nu[data.cuts], **binning, weights=w)
    exs, _ = numpy.histogram(data.nu[data.cuts], **binning)
    exs[exs > 0] = 1 / numpy.sqrt(exs[exs > 0])

    xs = xs * data.scale / acceptance
    lumi = ((data.charge / charge_e) *
            (density_target / a_target * avogadro * l_target))
    xs = xs / lumi * factor * dilution

    exs = exs * xs

    xs_list.append(xs)
    exs_list.append(exs)

nu = numpy.linspace(
    binning['range'][0] +
    (binning['range'][1] - binning['range'][0]) / binning['bins'] / 2,
    binning['range'][1] -
    (binning['range'][1] - binning['range'][0]) / binning['bins'] / 2,
    binning['bins'],
)

with PdfPages('xs.pdf') as pdf:
    plt.figure(figsize=(8, 6))
    ax = plt.gca()
    plt.xlabel(r'$\nu$')
    plt.ylabel(r"$\frac{d\sigma}{d\Omega dE'}$")
    plt.xlim(-100, 1400)
    plt.ylim(0, 1.2e4)
    for p0, xs, exs in zip(p0_list, xs_list, exs_list):
        s = (nu < e0 - 0.95 * p0) & (nu > e0 - 1.05 * p0)
        plt.errorbar(nu[s], xs[s], exs[s], fmt='k.', markersize=3)
    pdf.savefig(bbox_inches='tight')

    plt.ylim(0, 1e3)
    pdf.savefig(bbox_inches='tight')

    plt.xlim(-20, 400)
    plt.ylim(0, 1.2e4)
    for p0, xs, exs in zip(p0_list, xs_list, exs_list):
        s = (nu < e0 - 0.95 * p0) & (nu > e0 - 1.05 * p0)
        plt.plot(nu[s], xs[s], 'r-')
    pdf.savefig(bbox_inches='tight')
