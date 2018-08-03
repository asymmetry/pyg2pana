#!/usr/bin/env python3

import sqlite3
import time
from datetime import datetime
from os.path import dirname, join, realpath

__all__ = ['RunDB']


class RunDB():
    variables = [
        ('run_quality', 'RunQuality', 'int'),
        ('run_status', 'RunStatus', 'int'),
        ('septa_status', 'SeptaStatus', 'int'),
        ('hwp_status', 'Ihwp', 'int'),
        ('hwp_status_error', 'IhwpSTD', 'float'),
        ('ps1', 'ps1', 'float'),
        ('ps2', 'ps2', 'float'),
        ('ps3', 'ps3', 'float'),
        ('ps4', 'ps4', 'float'),
        ('ps7', 'ps7', 'float'),
        ('ps8', 'ps8', 'float'),
        ('target_encoder', 'TargetEncoder', 'float'),
        ('target_error', 'TargetSTD', 'float'),
        ('q1p', 'Q1p', 'float'),
        ('q1p_error', 'Q1pSTD', 'float'),
        ('q2p', 'Q2p', 'float'),
        ('q2p_error', 'Q2pSTD', 'float'),
        ('q3p', 'Q3p', 'float'),
        ('q3p_error', 'Q3pSTD', 'float'),
        ('d1p', 'D1p', 'float'),
        ('d1p_error', 'D1pSTD', 'float'),
        ('septa_current', 'SeptaI', 'float'),
        ('septa_current_error', 'SeptaSTD', 'float'),
        ('beam_energy', 'Energy', 'float'),
        ('beam_energy_error', 'EnergySTD', 'float'),
        ('trigger_eff', 'TEff', 'float'),
        ('deadtime', 'Deadtime', 'float'),
        ('cer_cut', 'CerCut', 'float'),
        ('pr1_cut', 'PR1Cut', 'float'),
        ('sum_cut', 'SumCut', 'float'),
        ('beam_pol', 'BeamPol', 'float'),
        ('beam_pol_stat_error', 'BeamPolStat', 'float'),
        ('beam_pol_sys_error', 'BeamPolSys', 'float'),
        ('bleed_through', 'Bleedthrough', 'float'),
        ('one_track_eff', 'OneTrackEff', 'float'),
        ('all_track_eff', 'AllTrackEff', 'float'),
        ('all_track_eff_low', 'AllTrackEffLow', 'float'),
        ('all_track_eff_high', 'AllTrackEffHigh', 'float'),
        ('target_pol', 'TargetPol', 'float'),
        ('target_pol_error', 'TargetPolError', 'float'),
        ('target_pol_stat_error', 'TargetPolStat', 'float'),
        ('target_pol_sys_error', 'TargetPolSys', 'float'),
        ('deadtime_plus', 'DTPlus', 'float'),
        ('deadtime_minus', 'DTMinus', 'float'),
        ('live_time_asym', 'LTAsym', 'float'),
        ('charge_asym', 'ChargeAsym', 'float'),
        ('charge_plus', 'QPlus', 'float'),
        ('charge_minus', 'QMinus', 'float'),
        ('charge', 'QTotal', 'float'),
        ('target_field', 'TargetField', 'float'),
        ('target_orientation', 'TargetOrientation', 'int'),
        ('material_id', 'MaterialID', 'int'),
        ('expert_comment', 'ExpertC', 'string'),
        ('target_cup', 'TargetCup', 'string'),
        ('run_start_time', 'RunStartTime', 'time'),
        ('entry_time', 'EntryTime', 'time'),
        ('run_stop_time', 'RunStopTime', 'time'),
        ('cer_eff', 'CerDetEff', 'float'),
        ('pr_eff', 'PRDetEff', 'float'),
        ('current', 'Current', 'float'),
        ('theta_cut_min', 'thCutMin', 'float'),
        ('theta_cut_max', 'thCutMax', 'float'),
        ('phi_cut_min', 'phCutMin', 'float'),
        ('phi_cut_max', 'phCutMax', 'float'),
        ('sim_cut_x', 'xBeam', 'float'),
        ('sim_cut_y', 'yBeam', 'float'),
        ('sim_cut_r', 'rBeam', 'float'),
        ('slow_raster_cut_x', 'xSR', 'float'),
        ('slow_raster_cut_y', 'ySR', 'float'),
        ('slow_raster_cut_r', 'rSR', 'float'),
    ]

    def __init__(self, run):
        if not isinstance(run, int):
            raise TypeError("run must be type 'int'")

        db_file = join(dirname(realpath(__file__)), 'g2p.db')
        self._con = sqlite3.connect(db_file)
        self._cur = self._con.cursor()

        for variable in self.__class__.variables:
            value = self._search(variable[1], run)
            if value is None:
                continue

            try:
                if variable[2] == 'int':
                    value = int(value)
                elif variable[2] == 'float':
                    value = float(value)
                elif variable[2] == 'string':
                    value = str(value)
                elif variable[2] == 'time':
                    value = datetime.strptime(value, '%Y-%m-%d %H:%M:%S')
                    value = int(time.mktime(value.timetuple()))
            except (TypeError, ValueError):
                value = None
            setattr(self, variable[0], value)

    def _search(self, field, run):
        table = 'AnaInfoR' if run > 20000 else 'AnaInfoL'
        command = 'Select {} from {}'.format(field, table)
        condition = 'where {} = {}'.format('RunNumber', run)

        try:
            self._cur.execute(command + ' ' + condition)
            r = self._cur.fetchone()
        except sqlite3.OperationalError:
            return None

        return r[0] if r is not None else None
