#!/usr/bin/env python3

from os.path import join

from pyg2pana import Data, SimFile, configs

run_list = configs.l_22545000

for key, value in run_list.items():
    p_runs = value['production']
    e_runs = value['empty']

    for run in p_runs + e_runs:
        files = [join('data', 'g2p_{}.root'.format(run))]
        files += [
            join('data', 'g2p_{}_{}.root'.format(run, x)) for x in range(1, 3)
        ]

        data = Data(files)
        data.save(join('data', 'g2p_{}.npz'.format(run)))
        del data

    for run in p_runs:
        files = [
            join('sim', 'sim_{}_{}.root'.format(run, x)) for x in range(1, 11)
        ]

        sim = SimFile(files)
        sim.save(join('sim', 'sim_{}.npz'.format(run)))
        del data
