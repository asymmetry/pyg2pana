#!/usr/bin/env python

import datetime
from os.path import join, dirname, realpath
import sqlite3
import time


class RunDB(object):

    def __init__(self, run):
        if not isinstance(run, int):
            raise TypeError("run must be type 'int'")

        self._con = sqlite3.connect(join(dirname(realpath(__file__)), 'g2p.db'))
        self._cur = self._con.cursor()

        variable_list = [
            ('RunQuality', 'int', 'run_quality'),
            ('RunStatus', 'int', 'run_status'),
            ('SeptaStatus', 'int', 'septa_status'),
            ('Ihwp', 'int', 'hwp_status'),
            ('IhwpSTD', 'float', 'hwp_status_error'),
            ('ps7', 'float', 'ps7'),
            ('ps8', 'float', 'ps8'),
            ('TargetEncoder', 'float', 'target_encoder'),
            ('TargetSTD', 'float', 'target_error'),
            ('Q1p', 'float', 'q1p'),
            ('Q1pSTD', 'float', 'q1p_error'),
            ('Q2p', 'float', 'q2p'),
            ('Q2pSTD', 'float', 'q2p_error'),
            ('Q3p', 'float', 'q3p'),
            ('Q3pSTD', 'float', 'q3p_error'),
            ('D1p', 'float', 'd1p'),
            ('D1pSTD', 'float', 'd1p_error'),
            ('SeptaI', 'float', 'septa_current'),
            ('SeptaSTD', 'float', 'septa_current_error'),
            ('Energy', 'float', 'beam_energy'),
            ('EnergySTD', 'float', 'beam_energy_error'),
            ('TEff', 'float', 'trigger_efficiency'),
            ('Deadtime', 'float', 'deadtime'),
            ('CerCut', 'float', 'cer_cut'),
            ('PR1Cut', 'float', 'pr1_cut'),
            ('SumCut', 'float', 'sum_cut'),
            ('BeamPol', 'float', 'beam_pol'),
            ('BeamPolStat', 'float', 'beam_pol_stat'),
            ('BeamPolSys', 'float', 'beam_pol_sys'),
            ('Bleedthrough', 'float', 'bleed_through'),
            ('OneTrackEff', 'float', 'one_track_eff'),
            ('AllTrackEff', 'float', 'all_track_eff'),
            ('AllTrackEffLow', 'float', 'all_track_eff_low'),
            ('AllTrackEffHigh', 'float', 'all_track_eff_high'),
            ('TargetPol', 'float', 'target_pol'),
            ('TargetPolError', 'float', 'target_pol_error'),
            ('TargetPolStat', 'float', 'target_pol_stat'),
            ('TargetPolSys', 'float', 'target_pol_sys'),
            ('DTPlus', 'float', 'dead_time_plus'),
            ('DTMinus', 'float', 'dead_time_minus'),
            ('LTAsym', 'float', 'live_time_asym'),
            ('ChargeAsym', 'float', 'charge_asym'),
            ('QPlus', 'float', 'charge_plus'),
            ('QMinus', 'float', 'charge_minus'),
            ('QTotal', 'float', 'charge_total'),
            ('TargetField', 'float', 'target_field'),
            ('TargetOrientation', 'int', 'target_orientation'),
            ('MaterialID', 'int', 'material_id'),
            ('ExpertC', 'string', 'expert_comment'),
            ('TargetCup', 'string', 'target_cup'),
            ('RunStartTime', 'time', 'run_start_time'),
            ('EntryTime', 'time', 'entry_time'),
            ('RunStopTime', 'time', 'run_stop_time'),
            ('CerDetEff', 'float', 'cer_eff'),
            ('PRDetEff', 'float', 'pr_eff'),
            ('Current', 'float', 'current'),
            ('thCutMin', 'float', 'theta_cut_min'),
            ('thCutMax', 'float', 'theta_cut_max'),
            ('phCutMin', 'float', 'phi_cut_min'),
            ('phCutMax', 'float', 'phi_cut_max'),
            ('xBeam', 'float', 'sim_cut_x'),
            ('yBeam', 'float', 'sim_cut_y'),
            ('rBeam', 'float', 'sim_cut_r'),
            ('xSR', 'float', 'slow_raster_cut_x'),
            ('ySR', 'float', 'slow_raster_cut_y'),
            ('rSR', 'float', 'slow_raster_cut_r'),
        ]

        if run > 20000:
            variable_list += [('ps1', 'float', 'ps1'), ('ps2', 'float', 'ps2')]
        else:
            variable_list += [('ps3', 'float', 'ps3'), ('ps4', 'float', 'ps4')]

        for variable in variable_list:
            if variable[1] == 'int':
                setattr(self, variable[2], self._int_search(variable[0], run))
            elif variable[1] == 'float':
                setattr(self, variable[2], self._float_search(variable[0], run))
            elif variable[1] == 'string':
                setattr(self, variable[2], self._string_search(variable[0], run))
            elif variable[1] == 'time':
                setattr(self, variable[2], self._time_search(variable[0], run))

    def _search(self, field, run):
        table = 'AnaInfoR' if run > 20000 else 'AnaInfoL'
        self._cur.execute('Select {} from {} where RunNumber = {}'.format(field, table, run))
        r = self._cur.fetchone()

        return r[0] if r is not None else None

    def _int_search(self, field, run):
        result = self._search(field, run)
        try:
            result = int(result)
        except (TypeError, ValueError):
            result = -999
        return result

    def _float_search(self, field, run):
        result = self._search(field, run)
        try:
            result = float(result)
        except (TypeError, ValueError):
            result = -999.0
        return result

    def _string_search(self, field, run):
        result = self._search(field, run)
        try:
            result = str(result)
        except (TypeError, ValueError):
            result = 'NULL'
        return result

    def _time_search(self, field, run):
        result = self._search(field, run)
        try:
            result = int(time.mktime(datetime.datetime.strptime(result, '%Y-%m-%d %H:%M:%S').timetuple()))
        except (TypeError, ValueError):
            result = -999
        return result
