"""
    Library for time evolution study.
"""

import system_state as ss
import re

class Experiment(object):
    """
        Has a list of system states, one for each t.
    """
    def __init__(self, system_states_list = None, t_diff = 0):
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
