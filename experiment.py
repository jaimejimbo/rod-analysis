"""
    Library for time evolution study.
"""

from system_state import *
import re
from methods import *
#import multiprocessing as mp    #for using all cores
import os
from rod import Rod
from system_state import SystemState



class Experiment(object):
    """
        Has a list of system states, one for each t.
    """
    def __init__(self, system_states_list = None, t_diff = 1):
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
        self._t_diff = t_diff
        self._evolution_dictionaries = []

    def _reset(self):
        """
        Resets all variables, so they must be computed again.
        """
        self._evolution_dictionaries = []

    def _create_dict_keys(self):
        """
        Create evolucion dictionaries keys.
        Each key is a rod's id.
        """
        for state in self._states:
            self._evolution_dictionaries.append({})
        for index in range(len(self._states)):
            state = self._states[index]
            evol_dict = self._evolution_dictionaries[index]
            for rod in state:
                evol_dict[rod.identifier] = None

    @property
    def evolution_dictionaries(self):
        """
        List of evolution dictionaries.
        Each dictionary has the form:
            {initial_rod_id1: final_rod_id1,
             initial_rod_id2: final_rod_id2,
             ...
             initial_rod_idN: final_rod_idN,}
        """
        if not len(self._evolution_dictionaries):
            self._create_dict_keys()
        return self._evolution_dictionaries


