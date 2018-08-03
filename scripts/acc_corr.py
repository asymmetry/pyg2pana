#!/usr/bin/env python3

from os.path import join

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy
from scipy import interpolate

from pyg2pana import Data, SimFile, configs

e0 = 2253.5
run_list = configs.l_22545000

cuts = {
    'y': [-0.015, 0.025],
    't': [-0.01, 0.03],
    'p': [-0.015, 0.015],
    'sr': 0.5,
}

binning = {
    'bins': 750,
    'range': (-100, 1400),
}

p0_list = []
hist_list = []
for key, value in run_list.items():
    p0 = float(key)
    if p0 > 1800:
        continue

    p0_list.append(p0)
    print(p0)

    run = value['production'][0]

    data = Data(join('data', 'g2p_{}.npz'.format(run)))
    sim = SimFile(join('sim', 'sim_{}.npz'.format(run)))

    data.cuts = cuts
    sim.cuts = cuts

    acceptance = sim.get_acceptance('nu', **binning)

    hist, _ = numpy.histogram(data.nu[data.cuts], **binning)
    hist = hist * data.scale / data.charge / acceptance

    hist_list.append(hist)

nu = []
yield_ = []
for p0, hist in zip(p0_list, hist_list):
    nu.append(e0 - p0)
    center = int((e0 - p0 - (-100)) / 2)
    temp = 0
    for i in range(center - 5, center + 6):
        temp = temp + hist[i]
    yield_.append(temp / 11)

nu = numpy.array(nu)
yield_ = numpy.array(yield_)
yield_func = interpolate.interp1d(
    nu, yield_, kind='cubic', fill_value='extrapolate')

x = numpy.linspace(
    binning['range'][0] +
    (binning['range'][1] - binning['range'][0]) / binning['bins'] / 2,
    binning['range'][1] -
    (binning['range'][1] - binning['range'][0]) / binning['bins'] / 2,
    binning['bins'],
)

dp_min = 0.039  # 3.9%
dp_max = 0.039
dp_list = []
corr_list = []
for p0, hist in zip(p0_list[:-1], hist_list[:-1]):
    center = int((e0 - p0 - (-100)) / 2)
    min_ = int((e0 - p0 * (1 + dp_min) - (-100)) / 2)
    max_ = int((e0 - p0 * (1 - dp_max) - (-100)) / 2)

    x_select = x[min_:max_]
    dp = (e0 - p0 - x_select) / p0
    corr = yield_func(x_select) / hist[min_:max_]
    dp_list.append(dp)
    corr_list.append(corr)

dp = numpy.concatenate(dp_list)
corr = numpy.concatenate(corr_list)

dp = dp[(corr > 0.7) & (corr < 1.3)]
corr = corr[(corr > 0.7) & (corr < 1.3)]

pars = numpy.polyfit(dp, corr, 3)
print(pars)

with PdfPages('result.pdf') as pdf:
    plt.figure(figsize=(8, 6))
    ax = plt.gca()
    plt.xlabel(r'$\nu$')
    plt.xlim(-100, 1400)
    plt.ylim(0, 2e6)
    for hist in hist_list:
        plt.plot(x, hist, 'k.', markersize=2.5)
    x_fit = x[(x > 500) & (x < 1250)]
    plt.plot(x_fit, yield_func(x_fit), 'r-', linewidth=1.5)
    pdf.savefig(bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(8, 6))
    ax = plt.gca()
    plt.xlabel(r'$dp$')
    plt.xlim(-0.05, 0.05)
    plt.ylim(0.6, 1.2)
    color_list = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    for color, dp, corr in zip(color_list, dp_list, corr_list):
        plt.plot(dp, corr, color + '.', markersize=2.5)
    dp_fit = numpy.linspace(-0.036, 0.0385, 100)
    plt.plot(
        dp_fit,
        sum(dp_fit**i * par for i, par in enumerate(pars[::-1])),
        'r-',
        linewidth=1,
    )
    pdf.savefig(bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(8, 6))
    ax = plt.gca()
    plt.xlabel(r'$\nu$')
    plt.xlim(-100, 1400)
    plt.ylim(0, 2e6)
    for p0, hist in zip(p0_list, hist_list):
        dp = (e0 - p0 - x) / p0
        select = (dp < 0.039) & (dp > -0.039)
        hist[select] *= sum(
            dp[select]**i * par for i, par in enumerate(pars[::-1]))
        plt.plot(x, hist, 'k.', markersize=2.5)
    pdf.savefig(bbox_inches='tight')
    plt.close()
