#!/usr/bin/env python

from os.path import exists, splitext

import numpy


class Data(object):

    def __init__(self, files, *args, **kwargs):
        self.__var_list = ['hel', 'bpm', 'sr', 'gold', 'rec']

        if not isinstance(files, (list, tuple)):
            files = [files]

        if all(ext == '.root' for _, ext in map(splitext, files)):
            self._load_root(files, *args, **kwargs)
        elif len(files) == 1 and splitext(files[0])[1] == '.npz':
            self._load_numpy(files[0], *args, **kwargs)

    def _load_root(self, files, run_number, *args, **kwargs):
        from root_numpy import tree2array
        from ROOT import TChain

        from .run_db import RunDB

        t = TChain('T')
        for file_ in files:
            if (exists(file_)):
                t.Add(file_)

        db = RunDB(run_number)

        p0 = db.d1p * 1000
        cer_cut = db.cer_cut
        pr1_cut = db.pr1_cut * p0
        sum_cut = db.sum_cut * p0

        if run_number < 20000:
            arm = 'L'
            cut = 'DL.evtypebits!=0&&L.tr.n==1&&L.cer.asum_c>{0}&&L.prl1.e>{1}&&(L.prl1.e+L.prl2.e)>{2}'.format(cer_cut, pr1_cut, sum_cut)
        else:
            arm = 'R'
            cut = 'DR.evtypebits!=0&&R.tr.n==1&&R.cer.asum_c>{0}&&R.ps.e>{1}&&(R.ps.e+R.sh.e)>{2}'.format(cer_cut, pr1_cut, sum_cut)

        cut += '&&{0}rb.bpmavail>0.5'.format(arm)

        hel_vars = [x.format(arm) for x in ['hel.{}.hel_act', 'hel.{}.error']]
        bpm_vars = [x.format(arm) for x in ['{}rb.tgt_0_x', '{}rb.tgt_0_y', '{}rb.tgt_0_theta', '{}rb.tgt_0_phi']]
        sr_vars = [x.format(arm) for x in ['{}rb.Raster.rawcurSL.x', '{}rb.Raster.rawcurSL.y']]
        gold_vars = [x.format(arm) for x in ['{}.gold.th', '{}.gold.y', '{}.gold.ph']]
        rec_vars = [x.format(arm) for x in ['{}.rec.x', '{}.rec.th', '{}.rec.y', '{}.rec.ph', '{}.rec.dp']]

        stop = kwargs.get('stop', None)
        all_vars = tree2array(t, hel_vars + bpm_vars + sr_vars + gold_vars + rec_vars, cut, stop=stop)

        hel = numpy.rec.fromarrays([all_vars[x] for x in hel_vars], names='val,err')
        bpm = numpy.rec.fromarrays([all_vars[x] for x in bpm_vars], names='x,y,t,p')
        sr = numpy.rec.fromarrays([all_vars[x] for x in sr_vars], names='x,y')
        gold = numpy.rec.fromarrays([all_vars[x] for x in gold_vars], names='t,y,p')
        rec = numpy.rec.fromarrays([all_vars[x] for x in rec_vars], names='x,t,y,p,d')

        for var in self.__var_list:
            setattr(self, var, eval(var))

    def _load_numpy(self, file_):
        if exists(file_):
            loaded = numpy.load(file_)

            for var in self.__var_list:
                setattr(self, var, loaded[var].view(numpy.recarray))

    def save(self, file_):
        if all(hasattr(self, x) for x in self.__var_list):
            numpy.savez_compressed(file_, **{var: getattr(self, var) for var in self.__var_list})
        else:
            raise ValueError('attribute not exists')
