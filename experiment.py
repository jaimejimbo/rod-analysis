"""
    Library for time evolution study.
"""
from system_state import *
import re
from methods import *
import multiprocessing as mp    #for using all cores
import os
from system_state import SystemState
import Queue

class Experiment(object):
    """
        Has a list of system states, one for each t.
    """
    def __init__(self, system_states_name_list=None, system_states_list = None, diff_t = 1):
        """
            Creation of experiment object.
        """
        type_ = str(type(system_states_list))
        if re.match(r'.*NoneType.*', type_):
            self._states = []
        elif re.match(r'.*list.*', type_):
            self._states = system_states_list
        else:
            raise TypeError
        type_ = str(type(system_states_name_list))
        if re.match(r'.*NoneType.*', type_):
            self._state_numbers = []
        elif re.match(r'.*list.*', type_):
            self._state_numbers = [get_number_from_string(num)
                                   for num in system_states_name_list]
        else:
            raise TypeError
        self._states_dict = {}
        for index in range(len(self._state_numbers)):
            number = self._state_numbers[index]
            state = self._states[index]
            self._states_dict[number] = state
        self._diff_t = diff_t
        self._evolution_dictionaries = []
        self._conflictive_final_rods = []
        self._relative_dictionaries = []
        self._unjoined_initial_rods = []
        self._unjoined_final_rods = []
        self._final_rods = []
        self._initial_rods = []

    def __getitem__(self, state_num):
        """
            Get the state identified by state_num.
        """
        type_ = str(type(state_num))
        if re.match(r'.*str.*', type_):
            identifier = get_number_from_string(state_num)
        elif re.match(r'.*int.*', type_):
            identifier = state_num
        else:
            raise ValueError
        return self._states_dict[identifier]

    def _reset(self):
        """
            Resets all variables, so they must be computed again.
        """
        self._evolution_dictionaries = []
        self._conflictive_final_rods = []
        self._relative_dictionaries = []
        self._unjoined_initial_rods = []
        self._unjoined_final_rods = []
        self._final_rods = []
        self._initial_rods = []

    def _create_dict_keys(self):
        """
            Create evolucion dictionaries keys.
        Each key is a rod's id.
        """
        for dummy_index in range(len(self._states)-1):
            self._evolution_dictionaries.append({})
            self._relative_dictionaries.append({})
            self._conflictive_final_rods.append(set([]))
            self._final_rods.append(set([]))
            self._initial_rods.append(set([]))
        for index in range(len(self._states)-1):
            state = self._states[index]
            evol_dict = self._evolution_dictionaries[index]
            relative_dict = self._relative_dictionaries[index]
            initial_set = self._initial_rods[index]
            for rod in state:
                initial_set |= set([rod.identifier])
                rod_id = rod.identifier
                evol_dict[rod_id] = set([])
                relative_dict[rod_id] = {}
        for index in range(len(self._states)-2):
            final_set = self._final_rods[index]
            final_set |= self._initial_rods[index+1]

    def _fill_dicts(self, max_speed, max_angle_diff):
        """
            Looks for rods that have only one possible predecessor.
        """
        processes = []
        output_queue = mp.Queue()
        for index in range(len(self._states)-1):
            processes.append(mp.Process(target=self._fill_dicts_process_limited,
                                        args=(index, max_speed, max_angle_diff, output_queue)))
        run_processes(processes)
        try:
            while True:
                output_row = output_queue.get(False)
                self._evolution_dictionaries[output_row[0]] = output_row[1]
                self._relative_dictionaries[output_row[0]] = output_row[2]
        except Queue.Empty:
            pass

    def _fill_dicts_process(self, index, max_speed, max_angle_diff, output_queue, amount_of_rods=100):
        """
            Allows to create a process and use all cores.
        """
        initial_state = self._states[index]
        final_state = self._states[index+1]
        evol_dict = self._evolution_dictionaries[index]
        relative_dict = self._relative_dictionaries[index]
        for initial_rod in initial_state:
            initial_id = initial_rod.identifier
            start_id = initial_id-amount_of_rods/2
            end_id = initial_id+amount_of_rods/2
            available_final_rods = final_state.get_rods_range(start_id, end_id)
            for final_rod in final_state:
                final_id = final_rod.identifier
                distance = initial_rod.distance_to_rod(final_rod)
                angle = initial_rod.angle_between_rods(final_rod)
                speed = float(distance)/self._diff_t
                if speed <= max_speed and angle <= max_angle_diff:
                    evol_dict[initial_id] |= set([final_id])
                    relative_dict[initial_id][final_id] = (distance, angle, speed)
        output_queue.put([index, evol_dict, relative_dict])

    def _fill_dicts_process_limited(self, index, max_speed, max_angle_diff, output_queue, limit=5, amount_of_rods=100):
        """
            Allows to create a process and use all cores.
        It limits possible final rods amount.
        """
        initial_state = self._states[index]
        final_state = self._states[index+1]
        evol_dict = self._evolution_dictionaries[index]
        relative_dict = self._relative_dictionaries[index]
        for initial_rod in initial_state:
            initial_id = initial_rod.identifier
            start_id = initial_id-amount_of_rods/2
            end_id = initial_id+amount_of_rods/2
            available_final_rods = final_state.get_rods_range(start_id, end_id)
            speeds = []
            initial_kappa = initial_rod.kappa
            kappa_error = initial_state.kappa_error/2
            for final_rod in available_final_rods:
                final_kappa = final_rod.kappa
                #rods can't change kappa.
                if initial_kappa-kappa_error > final_kappa or initial_kappa+kappa_error < final_kappa:
                    continue
                final_id = final_rod.identifier
                distance = initial_rod.distance_to_rod(final_rod)
                angle = initial_rod.angle_between_rods(final_rod)
                speed = float(distance)/self._diff_t
                if speed <= max_speed and angle <= max_angle_diff:
                    speeds.append([speed, final_id])
                    speeds.sort()
                    if len(speeds) > limit:
                        highest_speed, rod_id = speeds.pop(-1)
                    else:
                        highest_speed, rod_id = speeds[-1]
                        evol_dict[initial_id] |= set([final_id])
                        relative_dict[initial_id][final_id] = (distance, angle, speed)
                        continue
                    if rod_id != final_id:
                        evol_dict[initial_id] |= set([final_id])
                        relative_dict[initial_id][final_id] = (distance, angle, speed)
                        evol_dict[initial_id] -= set([rod_id])
                        del relative_dict[initial_id][rod_id]
        output_queue.put([index, evol_dict, relative_dict])

    def _use_unique_evolutions(self, max_reps=50):
        """
            Checks if there is only one possible evolution for the rods. If so,
        delete that possible evolution from the rest of rods.
        """
        changed = True
        rep = 0
        while changed:
            rep += 1
            if rep > max_reps:
                break
            changes_queue = mp.Queue()
            output_queue = mp.Queue()
            processes = []
            for index in range(len(self._evolution_dictionaries)):
                processes.append(mp.Process(target=self._use_unique_evolutions_process,
                                            args=(index, changes_queue, output_queue)))
            run_processes(processes)
            changed = False
            try:
                while True:
                    system_changed = changes_queue.get(False)
                    if system_changed:
                        changed = True
            except Queue.Empty:
                pass
            try:
                while True:
                    output = output_queue.get(False)
                    self._evolution_dictionaries[output[0]] = output[1]
                    self._conflictive_final_rods[output[0]] |= output[2]
            except Queue.Empty:
                pass

    def _use_unique_evolutions_process(self, index, changes_queue, output_queue):
        """
            Process for method.
        """
        evol_dict = self._evolution_dictionaries[index]
        conflicts = set([])
        changed = False
        for initial_rod_id in evol_dict.keys():
            conflicts2 = set([])
            final_rods = evol_dict[initial_rod_id]
            if len(final_rods) == 1:
                changed, evol_dict, conflicts2 = self._remove_final_rod(index, initial_rod_id, final_rods)
            conflicts |= conflicts2
        changes_queue.put(changed)
        output_queue.put([index, evol_dict, conflicts])

    def _remove_final_rod(self, index, initial_rod_id, final_rod_id):
        """
            Remove final rod from the rest of rods' evolutions.
        """
        evol_dict = self._evolution_dictionaries[index]
        system_end = self._states[index+1]
        final_rod = system_end[list(final_rod_id)[0]]
        changed = False
        conflicts = set([])
        final_rod_set = set([final_rod])
        for rod_id in evol_dict.keys():
            last_changed = False
            final_rods = evol_dict[rod_id]
            counter = 0
            if rod_id != initial_rod_id and len(final_rods):
                final_rods -= final_rod_set
                last_changed = True
            if not len(final_rods) and last_changed:
                final_rods |= final_rod_set
                conflicts |= final_rod_set
                last_changed = False
            if last_changed:
                changed = True
        return changed, evol_dict, conflicts

    def _leave_only_closer(self):
        """
            Leaves only the closer rod of possible evolutions.
        It will repeat final rods!
        """
        output_queue = mp.Queue()
        selected_queue = mp.Queue()
        processes = []
        for index in range(len(self._evolution_dictionaries)):
            processes.append(mp.Process(target=self._leave_only_closer_process,
                                        args=(index, output_queue, selected_queue)))
        run_processes(processes)
        try:
            while True:
                output = output_queue.get(False)
                index = output[0]
                initial_rod_id = output[1]
                final_rod_id = output[2]
                evol_dict = self._evolution_dictionaries[index]
                evol_dict[initial_rod_id] = final_rod_id
        except Queue.Empty:
            pass
        try:
            while True:
                selected = selected_queue.get(False)
                index = selected[0]
                self._final_rods[index] -= selected[1]
        except Queue.Empty:
            pass

    def _leave_only_closer_process(self, index, output_queue, selected_queue):
        """
        Process.
        """
        selected = set([])
        evol_dict = self._evolution_dictionaries[index]
        for initial_rod_id in evol_dict.keys():
            final_rod_id = self._closer_rod(index, initial_rod_id, selected)
            output_queue.put([index, initial_rod_id, final_rod_id])
            selected |= set([final_rod_id])
        selected_queue.put([index, selected])

    def _closer_rod(self, index, initial_rod_id, selected):
        """
            If there are multiple choices,
        this erase all but the closest.
        """
        evol_dict = self._evolution_dictionaries[index]
        relative_dict = self._relative_dictionaries[index]
        min_distance = 1e100
        final_rod = None
        for final_rod_id in evol_dict[initial_rod_id]:
            relative_values = relative_dict[initial_rod_id][final_rod_id]
            distance = relative_values[0]
            if distance < min_distance and final_rod_id not in selected:
                final_rod = final_rod_id
                min_distance = distance
        return final_rod

    def evolution_dictionaries(self, max_speed=100, max_angle_diff=90):
        """
                List of evolution dictionaries.
        Each dictionary has the form:
            {initial_rod_id1: set([final_rod_id11,final_rod_id12,...]),
             initial_rod_id2: set([final_rod_id21,final_rod_id22,...]),
             ...
             initial_rod_idN: set([final_rod_idN1,final_rod_idN2,...])}
        """
        if not len(self._evolution_dictionaries):
            self._create_dict_keys()
            self._fill_dicts(max_speed, max_angle_diff)
            #self._use_unique_evolutions()
            #self._leave_only_closer()
            #self._join_rods_left()
        return self._evolution_dictionaries

    def _join_rods_left(self):
        """
        After using methods listed before, some rods are unjoined.
        This joins closest rods.
        """
        for index in range(len(self._evolution_dictionaries)-1):
            evol_dict = self._evolution_dictionaries[index]
            initial_rods = set([])
            for initial_rod in evol_dict.keys():
                if not evol_dict[initial_rod]:
                    initial_rods |= set([initial_rod])
            self._initial_rods[index] = initial_rods

    def average_quadratic_speed(self, max_speed=100, max_angle_diff=90):
        """
        """
        pass

    def average_quadratic_angular_speed(self, max_speed=100, max_angle_diff=90):
        """
        """
        pass


