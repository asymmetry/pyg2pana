#!/usr/bin/env python3

import re
from os.path import exists, splitext

import numpy

from .run_db import RunDB

__all__ = ['SimFile']


class SimFile():

    def __init__(self, files, *args, **kwargs):
        self._db = None
        self._ref_db = None
        self._cuts = None

        self._var_list = ['bpm', 'rec', 'xs']

        if not isinstance(files, (list, tuple)):
            files = [files]

        if all(ext == '.root' for _, ext in map(splitext, files)):
            self.run = int(re.findall(r'sim_(\d+).*\.root', files[0])[0])
            self._load_root(files, *args, **kwargs)
        elif all(ext == '.npz' for _, ext in map(splitext, files)):
            self.run = int(re.findall(r'sim_(\d+).*\.npz', files[0])[0])
            self._load_numpy(files[0], *args, **kwargs)
        else:
            raise ValueError('bad filename')

        if self._db is None:
            self._db = RunDB(self.run)

        ref_run = kwargs.get('ref', None)
        if ref_run is not None:
            self._ref_db = RunDB(ref_run)
        else:
            self._ref_db = self._db

        self.e0 = self._db.beam_energy
        self.p0 = self._db.d1p * 1000

        self.range = {'d': 0.08, 't': 0.12, 'p': 0.08}

    def _load_root(self, files, *args, **kwargs):
        from root_numpy import tree2array
        import ROOT
        ROOT.gROOT.SetBatch(ROOT.kTRUE)

        tree = ROOT.TChain('T')
        for file_ in files:
            if exists(file_):
                tree.Add(file_)

        self.N = tree.GetEntries()

        cut = 'isgood>0.5'

        bpm_vars = ['bpm.l_x', 'bpm.l_y', 'bpm.l_t', 'bpm.l_p']
        rec_vars = ['rec.x', 'rec.t', 'rec.y', 'rec.p', 'rec.d']
        xs_vars = ['phys.react.xs']

        all_vars = tree2array(
            tree,
            bpm_vars + rec_vars + xs_vars,
            cut,
            stop=kwargs.get('stop', None),
        )

        self.bpm = numpy.rec.fromarrays(
            [all_vars[x] for x in bpm_vars],
            names='x,y,t,p',
        )
        self.rec = numpy.rec.fromarrays(
            [all_vars[x] for x in rec_vars],
            names='x,t,y,p,d',
        )
        self.xs = numpy.rec.fromarrays(
            [all_vars[x] for x in xs_vars],
            names='val',
        )

    def _load_numpy(self, file_):
        loaded = numpy.load(file_)
        for var in self._var_list:
            setattr(self, var, loaded[var].view(numpy.recarray))
        self.N = int(loaded['n'][0])

    def save(self, file_):
        if all(hasattr(self, x) for x in self._var_list):
            arrays = {x: getattr(self, x) for x in self._var_list}
            n = numpy.array([self.N])
            numpy.savez_compressed(file_, **arrays, n=n)
        else:
            raise ValueError('attributes do not exist')

    @property
    def cuts(self):
        if self._cuts is None:
            return numpy.ones_like(self.rec.d, dtype=bool)
        cut_t = ((self.rec.t > self._cuts['t'][0]) &
                 (self.rec.t < self._cuts['t'][1]))
        cut_p = ((self.rec.p > self._cuts['p'][0]) &
                 (self.rec.p < self._cuts['p'][1]))
        db = self._db if self._ref_db is None else self._ref_db
        cut_sr = ((self.bpm.x - db.sim_cut_x * 1e-3)**2 +
                  (self.bpm.y - db.sim_cut_y * 1e-3)**2 <
                  (db.sim_cut_r * self._cuts['sr'] * 1e-3)**2)
        return cut_t & cut_p & cut_sr

    @cuts.setter
    def cuts(self, value):
        if (isinstance(value, dict)
                and all(x in value.keys() for x in ['t', 'p'])):
            self._cuts = value
        else:
            raise ValueError('bad cuts')

    @property
    def nu(self):
        return self.e0 - self.p0 * (1 + self.rec.d)

    def get_acceptance(self, var, *args, **kwargs):
        vv = getattr(self, var)
        hist, _ = numpy.histogram(
            vv[self.cuts],
            bins=kwargs.get('bins'),
            range=kwargs.get('range'),
        )
        ave = numpy.average(hist[hist > 0])
        hist[hist < 0.8 * ave] = numpy.iinfo(numpy.int64).max
        norm = self.N / (
            self.range['t'] * self.range['p'] * self.range['d'] * self.p0)
        return hist / norm
