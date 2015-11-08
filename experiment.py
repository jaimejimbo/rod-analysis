"""
    Library for time evolution study.
"""

from system_state import *
import re
from methods import *
import multiprocessing as mp    #for using all cores
import os
from system_state import SystemState



class Experiment(object):
    """
        Has a list of system states, one for each t.
    """
    def __init__(self, system_states_list = None, diff_t = 1):
        """
            Creation of experiment object.
        """
        type_ = str(type(system_states_list))
        if re.match(r'.*NoneType.*', type_):
            self._states = []
        elif re.match(r'.*list.*', type_):
            self._states = system_states_list
        else:
            raise TypeError;
        self._diff_t = diff_t
        self._evolution_dictionaries = []

    def _reset(self):
        """
        Resets all variables, so they must be computed again.
        """
        self._evolution_dictionaries = []

    def _create_evolution_dict_keys(self):
        """
        Create evolucion dictionaries keys.
        Each key is a rod's id.
        """
        for dummy_index in range(len(self._states)-1):
            self._evolution_dictionaries.append({})
        for index in range(len(self._states)-1):
            state = self._states[index]
            evol_dict = self._evolution_dictionaries[index]
            for rod in state:
                evol_dict[rod.identifier] = set([])

    def _fill_evolution_dicts(self, max_speed, max_angle_diff):
        """
        Looks for rods that have only one possible predecessor.
        """
        for index in range(len(self._states)-1):
            initial_state = self._states[index]
            final_state = self._states[index+1]
            evol_dict = self._evolution_dictionaries[index]
            for initial_rod in initial_state:
                for final_rod in final_state:
                    distance = initial_rod.distance_to_rod(final_rod)
                    angle = abs(initial_rod.angle_between_rods(final_rod))
                    angle = min(angle, abs(180-angle))
                    speed = float(distance)/self._diff_t
                    if speed <= max_speed and angle <= max_angle_diff:
                        evol_dict[initial_rod.identifier] |= set([final_rod.identifier])

    def _clear_evolution_dicts(self):
        """
        Checks if there is only one possible evolution for the rods. If so, 
        delete that possible evolution from the rest of rods.
        """
        changed = False
        for evol_dict in self._evolution_dictionaries:
            for initial_rod_id in evol_dict.keys():
                final_rods = evol_dict[initial_rod_id]
                if len(final_rods) == 1:
                    changed = True
                    _remove_final_rod(evol_dict, initial_rod_id, final_rods)
        return changed

    def evolution_dictionaries(self, max_speed, max_angle_diff):
        """
        List of evolution dictionaries.
        Each dictionary has the form:
            {initial_rod_id1: set([final_rod_id11,final_rod_id12,...]),
             initial_rod_id2: set([final_rod_id21,final_rod_id22,...]),
             ...
             initial_rod_idN: set([final_rod_idN1,final_rod_idN2,...])}
        """
        if not len(self._evolution_dictionaries):
            self._create_evolution_dict_keys()
            self._fill_evolution_dicts(max_speed, max_angle_diff)
            self._clear_evolution_dicts()
        return self._evolution_dictionaries            


def _remove_final_rod(evol_dict, initial_rod_id, final_rod):
    """
    Remove final rod from the rest of rods' evolutions.
    """
    for rod_id in evol_dict.keys():
        if rod_id != initial_rod_id:
            evol_dict[rod_id] -= final_rod
