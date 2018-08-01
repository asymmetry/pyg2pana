#!/usr/bin/env python

from os.path import exists, splitext

import numpy


class SimResult(object):

    def __init__(self, files, *args, **kwargs):
        self.__var_list = ['bpm', 'rec', 'xs']

        if not isinstance(files, (list, tuple)):
            files = [files]

        if all(ext == '.root' for _, ext in map(splitext, files)):
            self._load_root(files, *args, **kwargs)
        elif len(files) == 1 and splitext(files[0])[1] == '.npz':
            self._load_numpy(files[0], *args, **kwargs)

    def _load_root(self, files, *args, **kwargs):
        from root_numpy import tree2array
        from ROOT import TChain

        t = TChain('T')
        for file_ in files:
            if (exists(file_)):
                t.Add(file_)

        cut = 'isgood>0.5'

        bpm_vars = ['bpm.l_x', 'bpm.l_y', 'bpm.l_t', 'bpm.l_p']
        rec_vars = ['rec.x', 'rec.t', 'rec.y', 'rec.p', 'rec.d']
        xs_vars = ['phys.react.xs']

        stop = kwargs.get('stop', None)
        all_vars = tree2array(t, bpm_vars + rec_vars + xs_vars, cut, stop=stop)

        bpm = numpy.rec.fromarrays([all_vars[x] for x in bpm_vars], names='x,y,t,p')
        rec = numpy.rec.fromarrays([all_vars[x] for x in rec_vars], names='x,t,y,p,d')
        xs = numpy.rec.fromarrays([all_vars[x] for x in xs_vars], names='val')

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
