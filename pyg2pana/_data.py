# Author: Chao Gu, 2018

import re
from os.path import exists, splitext

import numpy as np

from ._run_db import RunDB

__all__ = ['Data']


class Data():
    """
    Data Reader
    -----------
    Load g2p data rootfiles, extract kinematics for each event and save these
    information into a numpy npz file.

    Parameters
    ----------
    files : sequence of str
        If all files are rootfiles, root_numpy module is used to extract the
        kinematics. If all files are npz files, they are directly loaded by
        numpy.
    """

    def __init__(self, files, *args, **kwargs):
        self._db = None
        self._ref_db = None
        self._cuts = None

        self._var_list = ['hel', 'bpm', 'sr', 'gold', 'rec']

        if not isinstance(files, (list, tuple)):
            files = [files]

        if all(ext == '.root' for _, ext in map(splitext, files)):
            self.run = int(re.findall(r'g2p_(\d+).*\.root', files[0])[0])
            self._load_root(files, *args, **kwargs)
        elif all(ext == '.npz' for _, ext in map(splitext, files)):
            self.run = int(re.findall(r'g2p_(\d+).*\.npz', files[0])[0])
            self._load_numpy(files[0])
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
        self.charge = self._db.charge
        self.charge_plus = self._db.charge_plus
        self.charge_minus = self._db.charge_minus

    def _load_root(self, files, **kwargs):
        from root_numpy import tree2array
        import ROOT
        ROOT.gROOT.SetBatch(ROOT.kTRUE)

        tree = ROOT.TChain('T')
        for file_ in files:
            if exists(file_):
                tree.Add(file_)

        self._db = RunDB(self.run)
        p0 = self._db.d1p * 1000

        if self.run < 20000:
            arm = 'L'
            pr1_cut = 'L.prl1.e>{}'.format(self._db.pr1_cut * p0)
            sum_cut = '(L.prl1.e+L.prl2.e)>{}'.format(self._db.sum_cut * p0)
        else:
            arm = 'R'
            pr1_cut = 'R.ps.e>{}'.format(self._db.pr1_cut * p0)
            sum_cut = '(R.ps.e+R.sh.e)>{}'.format(self._db.sum_cut * p0)
        cut = '&&'.join([
            'D{}.evtypebits!=0'.format(arm),
            '{}.tr.n==1'.format(arm),
            '{}.cer.asum_c>{}'.format(arm, self._db.cer_cut),
            pr1_cut,
            sum_cut,
            '{}rb.bpmavail>0.5'.format(arm),
        ])

        hel_vars = ['hel.{}.{}'.format(arm, x) for x in ['hel_act', 'error']]
        bpm_vars = [
            '{}rb.tgt_0_{}'.format(arm, x) for x in ['x', 'y', 'theta', 'phi']
        ]
        sr_vars = [
            '{}rb.Raster.rawcurSL.{}'.format(arm, x) for x in ['x', 'y']
        ]
        gold_vars = [
            '{}.gold.{}'.format(arm, x) for x in ['th', 'y', 'ph', 'dp']
        ]
        rec_vars = [
            '{}.rec.{}'.format(arm, x) for x in ['x', 'th', 'y', 'ph', 'dp']
        ]

        all_vars = tree2array(
            tree,
            hel_vars + bpm_vars + sr_vars + gold_vars + rec_vars,
            cut,
            stop=kwargs.get('stop', None),
        )

        self.hel = np.rec.fromarrays(
            [all_vars[x] for x in hel_vars],
            names='val,err',
        )
        self.bpm = np.rec.fromarrays(
            [all_vars[x] for x in bpm_vars],
            names='x,y,t,p',
        )
        self.sr = np.rec.fromarrays(
            [all_vars[x] for x in sr_vars],
            names='x,y',
        )
        self.gold = np.rec.fromarrays(
            [all_vars[x] for x in gold_vars],
            names='t,y,p',
        )
        self.rec = np.rec.fromarrays(
            [all_vars[x] for x in rec_vars],
            names='x,t,y,p,d',
        )

    def _load_numpy(self, file_):
        loaded = np.load(file_)

        for var in self._var_list:
            setattr(self, var, loaded[var].view(np.recarray))

    def save(self, file_):
        if all(hasattr(self, x) for x in self._var_list):
            arrays = {x: getattr(self, x) for x in self._var_list}
            np.savez_compressed(file_, **arrays)
        else:
            raise ValueError('attributes do not exist')

    @property
    def cuts(self):
        if self._cuts is None:
            return np.ones_like(self.rec.d, dtype=bool)
        cut_y = ((self.gold.y > self._cuts['y'][0]) &
                 (self.gold.y < self._cuts['y'][1]))
        cut_t = ((self.rec.t > self._cuts['t'][0]) &
                 (self.rec.t < self._cuts['t'][1]))
        cut_p = ((self.rec.p > self._cuts['p'][0]) &
                 (self.rec.p < self._cuts['p'][1]))
        db = self._db if self._ref_db is None else self._ref_db
        cut_sr = ((self.sr.x - db.slow_raster_cut_x)**2 +
                  (self.sr.y - db.slow_raster_cut_y)**2 <
                  (db.slow_raster_cut_r * self._cuts['sr'])**2)
        return cut_y & cut_t & cut_p & cut_sr

    @cuts.setter
    def cuts(self, value):
        if (isinstance(value, dict)
                and all(x in value.keys() for x in ['y', 't', 'p'])):
            self._cuts = value
        else:
            raise ValueError('bad cuts')

    @property
    def scale(self):
        return self._efficiency_prescale / (1 - self._db.deadtime)

    @property
    def scale_plus(self):
        return self._efficiency_prescale / (1 - self._db.deadtime_plus)

    @property
    def scale_minus(self):
        return self._efficiency_prescale / (1 - self._db.deadtime_minus)

    @property
    def _efficiency_prescale(self):
        efficiency = self._ref_db.one_track_eff / self._ref_db.all_track_eff
        efficiency *= self._ref_db.trigger_eff
        efficiency *= self._ref_db.cer_eff * self._ref_db.pr_eff
        if self.run < 20000:
            prescale = self._db.ps3
        else:
            prescale = self._db.ps1
        return prescale / efficiency

    @property
    def nu(self):
        return self.e0 - self.p0 * (1 + self.rec.d)
