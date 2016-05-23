"""
    Library for time evolution study.
"""
import re, methods, math, copy, gc, sys, os
import multiprocessing as mp
from matplotlib import animation
import matplotlib
import matplotlib.pyplot as plt
import numpy
import numpy as np
import scipy.optimize as optimization
import datetime
import settings, time
from sys import getsizeof
import inspect

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'
if settings.special_chars:
    WHITE_BLOCK = u'\u25A0'
else:
    WHITE_BLOCK = 'X'

DEBUG = True

class Experiment(object):
    """
        Has a list of system states, one for each t.
    """
    def __init__(self, system_states_name_list=None, kappas=None,
                system_states_list=None, diff_t=1, dates=None):
        """
            Creation of experiment object.
        """

        type_ = str(type(system_states_list))
        if re.match(r'.*NoneType.*', type_):
            self._states = []
            self._number_of_states = 0
        elif re.match(r'.*list.*', type_):
            type_2 = str(type(system_states_list[0]))
            if re.match(r'.*SystemState.*', type_2):
                self._states = [methods.compress(state, level=settings.medium_comp_level) for state in system_states_list]
            elif re.match(r'.*str.*', type_2):
                self._states = system_states_list
            else:
                raise TypeError
            self._number_of_states = len(system_states_list)
        else:
            raise TypeError

        type_ = str(type(system_states_name_list))
        if re.match(r'.*NoneType.*', type_):
            self._state_numbers = []
        elif re.match(r'.*list.*', type_):
            self._state_numbers = [methods.get_number_from_string(num)
                            for num in system_states_name_list]
        else:
            raise TypeError

        self._states_dict = {}
        for index in range(self._number_of_states):
            number = self._state_numbers[index]
            state = self._states[index]
            self._states_dict[number] = methods.compress(state,
                level=settings.medium_comp_level)
        if not dates:
            msg = "Time evolution needs dates file.\n"
            msg = "To create one run export_image_dates().\n"
            raise ValueError(msg)
        self._diff_t = diff_t
        self._dates = dates
        self._kappas = kappas
        self._reset()

    def _reset(self):
        """
            Resets all variables, so they must be computed again.
        """
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "reset"
        self._bursts_computed = False
        self._evolution_dictionaries = []
        self._conflictive_final_rods = []
        self._relative_dictionaries = []
        self._unjoined_initial_rods = []
        self._unjoined_final_rods = []
        self._final_rods = []
        self._initial_rods = []
        self._speeds_vectors = []
        self._speeds = []
        self._angular_speeds = []
        self._max_distance = None
        self._max_angle_diff = None
        self._limit = None
        self._amount_of_rods = None
        self._local_speeds = []
        self._local_average_quadratic_speeds = []
        self._local_average_quadratic_angular_speeds = []
        self._image_id_by_index = {}
        self._densities_array = []
        self._quad_speeds_array = []
        self._speeds_matrices = []
        self._cluster_areas = []
        self._bursts_groups = []
        self._divisions = None
        self._min_density = None
        self._max_density = None
        self._indices = None
        self._done = False
        self._total_cluster_areas = None
        self._popt = None
        self._std_dev = None
        self._covered_area_prop = None
        self._image_ids_done = False
        self._average_kappa = None
        self._kappa_dev = None
        self._average_rad = None
        self._average_rod_length = None
        self._average_rod_width = None
        self._number_of_rods = None
        self._number_of_rods_dev = None

    def __len__(self):
        """
        len magic method
        """
        return len(self._states)

    def __iter__(self):
        """
            Magic method for looping.
        """
        for state in self._states:
            yield methods.decompress(state, level=settings.medium_comp_level)

    def __getitem__(self, state_num):
        """
            Get the state identified by state_num.
        """
        type_ = str(type(state_num))
        if re.match(r'.*str.*', type_):
            identifier = methods.get_number_from_string(state_num)
        elif re.match(r'.*int.*', type_):
            identifier = state_num
        else:
            raise ValueError
        return methods.decompress(self._states_dict[identifier])

    def get(self, index):
        """
            Getter with index
        """
        state = self._states[index]
        state = methods.decompress(state, level=settings.medium_comp_level)
        return state

    def set_coef(self, value):
        """
            Changes coef for dividing in circles.
        """
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Setting coeficient"
        states = []
        output_queue = mp.Queue()
        processes = []
        for index in range(len(self)):
            process = mp.Process(target=self.set_coef_process,
                                 args=(index, output_queue, value))
            processes.append(process)
            states.append(None)
        num_processes = len(processes)
        running, processes_left = methods.run_processes(processes)
        finished_ = 0
        previous_time = datetime.datetime.now()
        counter = 0
        time_left = None
        times = []
        print " "
        while finished_ < num_processes:
            counter += 1
            finished_ += 1
            previous_time, counter, time_left = methods.print_progress(finished_, num_processes, counter,
                                                   times, time_left, previous_time)
            output_row = output_queue.get()
            index = output_row[0]
            compressed_state = output_row[1]
            self._states[index] = compressed_state
            if len(processes_left):
                new_process = processes_left.pop(0)
                time.sleep(settings.waiting_time)
                new_process.start()
        print(CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE)

    def set_coef_process(self, index, output_queue, value):
        """
            Process.
        """
        state = self.get(index)
        state.coef = value
        output = methods.compress(state,
                    level=methods.settings.medium_comp_level)
        output_queue.put([index, output])

    def set_coords(self, zone_coords):
        """
            Changes coef for dividing in circles.
        """
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Setting coords"
        states = []
        output_queue = mp.Queue()
        processes = []
        for index in range(len(self)):
            process = mp.Process(target=self.set_coords_process,
                                 args=(index, output_queue, zone_coords))
            processes.append(process)
            states.append(None)
        num_processes = len(processes)
        running, processes_left = methods.run_processes(processes)
        finished_ = 0
        previous_time = datetime.datetime.now()
        counter = 0
        time_left = None
        times = []
        print " "
        while finished_ < num_processes:
            counter += 1
            finished_ += 1
            previous_time, counter, time_left = methods.print_progress(finished_, num_processes, counter,
                                                   times, time_left, previous_time)
            output_row = output_queue.get()
            index = output_row[0]
            compressed_state = output_row[1]
            self._states[index] = compressed_state
            if len(processes_left):
                new_process = processes_left.pop(0)
                time.sleep(settings.waiting_time)
                new_process.start()
        print(CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE)

    def set_coords_process(self, index, output_queue, zone_coords):
        """
            Process.
        """
        state = self.get(index)
        state.set_coords(zone_coords)
        output = methods.compress(state,
                    level=methods.settings.medium_comp_level)
        output_queue.put([index, output])

    def homogenize_coords(self):
        """
            Homogenizes coords
        """
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Homogenizing coords for grid plot"
        state = self.get(0)
        zone_coords = state.zone_coords
        self.set_coords(zone_coords)
        

    def _create_dict_keys(self):
        """
            Create evolucion dictionaries keys.
        Each key is a rod's id.
        """
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Creating dict keys"
        for dummy_index in range(len(self)):
            self._evolution_dictionaries.append(methods.compress({}, level=settings.medium_comp_level))
            self._relative_dictionaries.append(methods.compress({}, level=settings.medium_comp_level))
            self._conflictive_final_rods.append(set([]))
            self._final_rods.append(set([]))
            self._initial_rods.append(set([]))
        processes = []
        output_queue = mp.Queue()
        for index in range(len(self)):
            process = mp.Process(target=self._create_dicts_keys_process,
                                args=(index, output_queue))
            processes.append(process)
        num_processes = len(processes)
        running, processes_left = methods.run_processes(processes)
        finished_ = 0
        previous_time = datetime.datetime.now()
        counter = 0
        time_left = None
        times = []
        print " "
        while finished_ < num_processes:
            counter += 1
            finished_ += 1
            previous_time, counter, time_left = methods.print_progress(finished_, num_processes,
                                    counter, times, time_left, previous_time)
            [index, evol_dict, rel_dict, initial_rods] = output_queue.get()
            self._evolution_dictionaries[index] = evol_dict
            self._relative_dictionaries[index] = rel_dict
            self._initial_rods[index] = initial_rods
            if len(processes_left):
                new_process = processes_left.pop(0)
                time.sleep(settings.waiting_time)
                new_process.start()
        for index in range(len(self)-1):
            self._final_rods[index] = self._initial_rods[index+1]
        print(CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE)


    def _create_dicts_keys_process(self, index, output_queue):
        """
        Process
        """
        state = self.get(index)
        if not state:
            msg = "State is not defined."
            raise TypeError(msg)
        evol_dict = methods.decompress(self._evolution_dictionaries[index], level=settings.medium_comp_level)
        relative_dict = methods.decompress(self._relative_dictionaries[index], level=settings.medium_comp_level)
        initial_rods = set([])
        for rod_ in state:
            initial_rods |= set([rod_.identifier])
            rod_id = rod_.identifier
            evol_dict[rod_id] = set([])
            relative_dict[rod_id] = {}
        evol_dict = methods.compress(evol_dict, level=settings.medium_comp_level)
        rel_dict = methods.compress(relative_dict, level=settings.medium_comp_level)
        output_queue.put([index, evol_dict, rel_dict, initial_rods])

    def _fill_dicts(self, max_distance, max_angle_diff,
                    limit=20, amount_of_rods=None):
        """
            Looks for rods that have only one possible predecessor.
        """
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Filling dictionaries (rod evolution)"
        (self._max_distance, self._max_angle_diff,
        self._limit, self._amount_of_rods) = (max_distance, max_angle_diff,
                                                limit, amount_of_rods)
        processes = []
        output_queue = mp.Queue()
        for index in range(len(self)-1):
            processes.append(mp.Process(target=self._fill_dicts_process_limited,
                                    args=(index, max_distance, max_angle_diff,
                                        output_queue, limit, amount_of_rods)))
        num_processes = len(processes)
        running, processes_left = methods.run_processes(processes)
        finished_ = 0
        previous_time = datetime.datetime.now()
        counter = 0
        time_left = None
        times = []
        print " "
        while finished_ < num_processes:
            counter += 1
            finished_ += 1
            previous_time, counter, time_left = methods.print_progress(finished_, num_processes,
                                    counter, times, time_left, previous_time)
            output_row = output_queue.get()
            index = output_row[0]
            self._evolution_dictionaries[index] = output_row[1]
            self._relative_dictionaries[index] = output_row[2]
            if len(processes_left):
                new_process = processes_left.pop(0)
                time.sleep(settings.waiting_time)
                new_process.start()
        print(CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE)

    def _fill_dicts_process_limited(self, index, max_distance,
                                max_angle_diff, output_queue,
                                limit=20, amount_of_rods=None):
        """
            Allows to create a process and use all cores.
        It limits possible final rods amount.
        """
        initial_state = self.get(index)
        final_state = self.get(index+1)
        assert type(initial_state) != type("string"), "initial state can't be a string."
        evol_dict = methods.decompress(self._evolution_dictionaries[index], level=settings.medium_comp_level)
        relative_dict = methods.decompress(self._relative_dictionaries[index], level=settings.medium_comp_level)
        error = True
        for initial_rod in initial_state:
            initial_id = initial_rod.identifier
            if not amount_of_rods:
                available_final_rods = list(final_state)
            else:
                start_id = initial_id-amount_of_rods/2
                end_id = initial_id+amount_of_rods/2
                available_final_rods = final_state._get_rods_range(start_id,
                                                                  end_id)
            speeds = []
            initial_kappa = initial_rod.kappa
            kappa_error = initial_state.kappa_error/2
            for final_rod in available_final_rods:
                final_id = final_rod.identifier
                distance = initial_rod.distance_to_rod(final_rod)
                angle = initial_rod.angle_between_rods(final_rod)
                angle = abs(min([angle, 180-angle]))
                speed = float(distance)/self._diff_t
                if distance <= max_distance and angle <= max_angle_diff:
                    speeds.append([speed, final_id])
                    speeds.sort()
                    if len(speeds) > limit:
                        dummy, rod_id = speeds.pop(-1)
                    else:
                        dummy, rod_id = speeds[-1]
                        evol_dict[initial_id] |= set([final_id])
                        relative_dict[initial_id][final_id] = (distance, angle)
                        continue
                    if rod_id != final_id:
                        evol_dict[initial_id] |= set([final_id])
                        relative_dict[initial_id][final_id] = (distance, angle)
                        evol_dict[initial_id] -= set([rod_id])
                        del relative_dict[initial_id][rod_id]
        for initial_id in relative_dict.keys():
            if len(relative_dict[initial_id]) == 1:
                relative_dict[initial_id] = relative_dict[initial_id].values()[0]
        output_queue.put([index, methods.compress(evol_dict, level=settings.medium_comp_level), methods.compress(relative_dict, level=settings.medium_comp_level)])

    def _remove_final_rod(self, index, initial_rod_id, final_rod_id):
        """
            Remove final rod from the rest of rods' evolutions.
        """
        evol_dict = methods.decompress(self._evolution_dictionaries[index], level=settings.medium_comp_level)
        system_end = self.get(index+1)
        final_rod = system_end[list(final_rod_id)[0]]
        changed = False
        conflicts = set([])
        final_rod_set = set([final_rod])
        for rod_id in list(evol_dict.keys()):
            last_changed = False
            final_rods = evol_dict[rod_id]
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

    def _leave_only_closer(self, max_distance=50):
        """
            Leaves only the closer rod of possible evolutions.
        It will repeat final rods!
        """
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Leaving closer rod"
        output_queue = mp.Queue()
        selected_queue = mp.Queue()
        processes = []
        max_distance = 200
        for index in range(len(self._evolution_dictionaries)):
            processes.append(mp.Process(target=self._leave_only_closer_process,
                                        args=(index, output_queue,
                                              selected_queue, max_distance)))
        num_processes = len(processes)
        running, processes_left = methods.run_processes(processes)
        finished = 0
        previous_time = datetime.datetime.now()
        counter = 0
        time_left = None
        times = []
        print " "
        while finished < num_processes:
            counter += 1
            finished += 1
            previous_time, counter, time_left = methods.print_progress(finished, num_processes,
                                    counter, times, time_left, previous_time)
            output = output_queue.get()
            selected = selected_queue.get()
            index = output[0]
            evol_dict = output[1]
            relative_dict = output[2]
            self._evolution_dictionaries[index] = evol_dict
            self._relative_dictionaries[index] = relative_dict
            index = selected[0]
            self._final_rods[index] -= selected[1]
            if len(processes_left):
                new_process = processes_left.pop(0)
                time.sleep(settings.waiting_time)
                new_process.start()
        print(CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE)

    def _leave_only_closer_process(self, index, output_queue,
                                selected_queue, max_distance):
        """
        Process.
        """
        selected = set([])
        evol_dict = methods.decompress(self._evolution_dictionaries[index], level=settings.medium_comp_level)
        relative_dict = methods.decompress(self._relative_dictionaries[index], level=settings.medium_comp_level)
        for initial_rod_id in list(evol_dict.keys()):
            final_rod_id, distance, angle_diff = self._closer_rod(index,
                                          initial_rod_id, selected, evol_dict,
                                          relative_dict, max_distance)
            evol_dict[initial_rod_id] = final_rod_id
            relative_dict[initial_rod_id] = (distance, angle_diff)
            selected |= set([final_rod_id])
        output_queue.put([index, methods.compress(evol_dict, level=settings.medium_comp_level), methods.compress(relative_dict, level=settings.medium_comp_level)])
        selected_queue.put([index, selected])

    def _closer_rod(self, index, initial_rod_id, selected, evol_dict, relative_dict_, max_distance=50):
        """
            If there are multiple choices,
        this erase all but the closest.
        """
        final_rods = evol_dict[initial_rod_id]
        relative_dict = relative_dict_[initial_rod_id]
        min_distance = max_distance
        final_rod = -1
        final_rods_list = list(final_rods)
        if len(final_rods_list) == 1:
            final_rod = final_rods_list[0]
            min_distance = relative_dict[0]
            angle_diff = relative_dict[1]
        elif len(final_rods_list) == 0:
            return None, None, None
        else:
            while True:
                try:
                    final_rod_id = final_rods_list.pop(0)
                except IndexError:
                    return None, None, None
                if final_rod_id not in selected:
                    final_rod = final_rod_id
                    rel = relative_dict[final_rod_id]
                    min_distance = rel[0]
                    angle_diff = rel[1]
                    break
        return final_rod, min_distance, angle_diff

    def compute_dictionaries(self, max_distance=100, max_angle_diff=90,
                            limit=10, amount_of_rods=None):
        """
                List of evolution dictionaries.
        Each dictionary has the form:
            {initial_rod_id1: set([final_rod_id11,final_rod_id12,...]),
             initial_rod_id2: set([final_rod_id21,final_rod_id22,...]),
             ...
             initial_rod_idN: set([final_rod_idN1,final_rod_idN2,...])}
        """
        tuple1 = (max_distance, max_angle_diff, limit, amount_of_rods)
        tuple2 = (self._max_distance, self._max_angle_diff,
                    self._limit, self._amount_of_rods)
        if tuple1 != tuple2:
            self._reset()
        if not len(self._evolution_dictionaries):
            self._create_dict_keys()
            self._fill_dicts(max_distance, max_angle_diff, limit=limit,
                                amount_of_rods=amount_of_rods)
            self._leave_only_closer(max_distance=max_distance)
            self._join_rods_left(max_distance=max_distance)
            #self._get_vectors()

    def vectors_dictionaries(self, max_distance=100, max_angle_diff=90,
                            limit=5, amount_of_rods=200):
        """
        Returns a dictionary {initial_rod: vector_to_final_rod}
        """
        self.compute_dictionaries(max_distance=max_distance,
                                  max_angle_diff=max_angle_diff,
                                  limit=5, amount_of_rods=200)
        return self._speeds_vectors

    @property
    def speeds_vectors(self):
        """
         Rods' speed vectors.
        """
        if not len(self._speeds_vectors):
            self._get_vectors()
        return self._speeds_vectors

    def _get_vectors(self):
        """
        Creates a dictionary {initial_rod: vector_to_final_rod}
        """
        output_queue = mp.Queue()
        processes = []
        for index in range(len(self._evolution_dictionaries)-1):
            processes.append(mp.Process(target=self._get_vectors_process,
                                        args=(index, output_queue)))
            self._speeds_vectors.append(0)
        num_processes = len(processes)
        running, processes_left = methods.run_processes(processes)
        finished = 0
        while finished < num_processes:
            finished += 1
            output = output_queue.get()
            index = output[0]
            speeds_vectors = output[1]
            self._speeds_vectors[index] = speeds_vectors
            if len(processes_left):
                new_process = processes_left.pop(0)
                time.sleep(settings.waiting_time)
                new_process.start()

    def _get_vectors_process(self, index, output_queue):
        """
        Process.
        """
        speeds_vectors = {}
        evol_dict = methods.decompress(self._evolution_dictionaries[index], level=settings.medium_comp_level)
        initial_state = self.get(index)
        final_state = self.get(index+1)
        keys = list(evol_dict.keys())
        if not keys[0]:
            return
        for initial_rod_id in keys:
            if initial_rod_id:
                final_rod_id = evol_dict[initial_rod_id]
                if final_rod_id:
                    initial_rod = initial_state[initial_rod_id]
                    final_rod = final_state[final_rod_id]
                    vector = initial_rod.vector_to_rod(final_rod)
                    speeds_vectors[initial_rod_id] = vector
        output_queue.put([index, methods.compress(speeds_vectors, level=settings.medium_comp_level)])

    def evolution_dictionaries(self, max_distance=100, max_angle_diff=90,
                                limit=5, amount_of_rods=200):
        """
            List of evolution dictionaries.
        Each dictionary has the form:
            {initial_rod_id1: set([final_rod_id11,final_rod_id12,...]),
             initial_rod_id2: set([final_rod_id21,final_rod_id22,...]),
             ...
             initial_rod_idN: set([final_rod_idN1,final_rod_idN2,...])}
        """
        self.compute_dictionaries(max_distance=max_distance,
                                  max_angle_diff=max_angle_diff,
                                  limit=limit,
                                  amount_of_rods=amount_of_rods)
        return self._evolution_dictionaries

    def _join_rods_left(self, max_distance=50):
        """
        After using methods listed before, some rods are unjoined.
        This joins closest rods.
        """
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Joining lonely rods"
        output_queue = mp.Queue()
        processes = []
        for index in range(len(self._evolution_dictionaries)-1):
            process = mp.Process(target=self._join_rods_left_process,
                                 args=(index, output_queue, max_distance))
            processes.append(process)
        num_processes = len(processes)
        running, processes_left = methods.run_processes(processes)
        finished = 0
        previous_time = datetime.datetime.now()
        counter = 0
        time_left = None
        times = []
        print " "
        while finished < num_processes:
            counter += 1
            finished += 1
            previous_time, counter, time_left = methods.print_progress(finished, num_processes,
                                    counter, times, time_left, previous_time)
            output = output_queue.get()
            index = output[0]
            evol_dict = output[1]
            relative_dict = output[2]
            self._evolution_dictionaries[index] = evol_dict
            self._relative_dictionaries[index] = relative_dict
            if len(processes_left):
                new_process = processes_left.pop(0)
                time.sleep(settings.waiting_time)
                new_process.start()
        print(CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE)


    def _join_rods_left_process(self, index, output_queue, max_distance=50):
        """
        Process for method.
        """
        evol_dict = methods.decompress(self._evolution_dictionaries[index], level=settings.medium_comp_level)
        relative_dict = methods.decompress(self._relative_dictionaries[index], level=settings.medium_comp_level)
        initial_rods = set([])
        for initial_rod_id in evol_dict.keys():
            if evol_dict[initial_rod_id] is None:
                initial_rods |= set([initial_rod_id])
        self._initial_rods[index] = initial_rods
        final_rods = self._final_rods[index]
        initial_state = self.get(index)
        final_state = self.get(index+1)
        for final_rod_id in final_rods:
            min_distance = max_distance
            selected_rod = None
            selected_rod_id = None
            final_rod = final_state[final_rod_id]
            for initial_rod_id in initial_rods:
                initial_rod = initial_state[initial_rod_id]
                distance = final_rod.distance_to_rod(initial_rod)
                if distance < min_distance:
                    min_distance = distance
                    selected_rod_id = initial_rod_id
                    selected_rod = initial_rod
            evol_dict[selected_rod_id] = final_rod_id
            angle_diff = None
            try:
                angle_diff = final_rod.angle_between_rods(selected_rod)
                angle_diff = min([angle_diff, 180-angle_diff])
            except AttributeError:
                min_distance = None
                angle_diff = None
            relative_dict[selected_rod_id] = (min_distance, angle_diff)
            initial_rods -= set([selected_rod_id])
        output_queue.put([index, methods.compress(evol_dict, level=settings.medium_comp_level), methods.compress(relative_dict, level=settings.medium_comp_level)])

    def _compute_speeds(self, max_distance, max_angle_diff, limit, amount_of_rods):
        """
        Get speeds and angular speeds.
        """
        tuple1 = (max_distance, max_angle_diff, limit, amount_of_rods)
        tuple2 = (self._max_distance, self._max_angle_diff,
                    self._limit, self._amount_of_rods)
        if tuple1 != tuple2:
            self._reset()
        self.compute_dictionaries(max_distance, max_angle_diff,
                                limit, amount_of_rods)
        if not len(self._speeds):
            print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Computing speeds"
            speeds_queue = mp.Queue()
            angular_speeds_queue = mp.Queue()
            processes = []
            for index in range(len(self._evolution_dictionaries)-1):
                process = mp.Process(target=self._compute_speeds_process,
                                     args=(index, speeds_queue,
                                            angular_speeds_queue))
                self._speeds.append(None)
                self._angular_speeds.append(None)
                processes.append(process)
            num_processes = len(processes)
            running, processes_left = methods.run_processes(processes)
            finished = 0
            previous_time = datetime.datetime.now()
            counter = 0
            time_left = None
            times = []
            print " "
            while finished < num_processes:
                counter += 1
                finished += 1
                previous_time, counter, time_left = methods.print_progress(finished, num_processes,
                                        counter, times, time_left, previous_time)
                index_speeds, speeds = speeds_queue.get()
                index_angular_speeds, angular_speeds = angular_speeds_queue.get()
                self._speeds[index_speeds] = speeds
                self._angular_speeds[index_angular_speeds] = angular_speeds
                if len(processes_left):
                    new_process = processes_left.pop(0)
                    time.sleep(settings.waiting_time)
                    new_process.start()
            print CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE

    def _compute_speeds_process(self, index, speeds_queue,
                                angular_speeds_queue):
        """
        Returns an array of speeds.
        """
        rel_dict = methods.decompress(self._relative_dictionaries[index], level=settings.medium_comp_level)
        speeds = {}
        angular_speeds = {}
        error = True
        for initial_rod_id in list(rel_dict.keys()):
            values = rel_dict[initial_rod_id]
            try:
                error &= (not bool(values[0]))
                speed = float(values[0])/self._diff_t
                angular_speed = float(values[1])/self._diff_t
                speeds[initial_rod_id] = speed
                angular_speeds[initial_rod_id] = angular_speed
            except TypeError:
                pass
        if error:
            speeds_queue.put([index, None])
            angular_speeds_queue.put([index, None])
            return
        speeds_queue.put([index, methods.compress(speeds, level=settings.medium_comp_level)])
        angular_speeds_queue.put([index, methods.compress(angular_speeds, level=settings.medium_comp_level)])

    def speeds(self, max_distance=100, max_angle_diff=90, limit=5,
                     amount_of_rods=200):
        """
        Returns [speeds, angular_speeds] where both outputs are array with one
        value for each rod.
        """
        self._compute_speeds(max_distance, max_angle_diff, limit, amount_of_rods)
        return [self._speeds, self._angular_speeds]

    def average_quadratic_speed(self, max_distance=100, max_angle_diff=90,
                                limit=5, amount_of_rods=200):
        """
        Returns average quadratic speeds
        """
        self._compute_speeds(max_distance, max_angle_diff, limit, amount_of_rods)
        output = []
        for index in range(len(self._speeds)):
            speeds = methods.decompress(self._speeds[index], level=settings.medium_comp_level)
            num_of_rods = len(speeds)
            output.append(0)
            for speed in list(speeds.values()):
                output[index] += speed**2/num_of_rods
        return output

    def average_quadratic_angular_speed(self, max_distance=100, max_angle_diff=90,
                                        limit=5, amount_of_rods=200):
        """
        Returns average quadratic angular speed
        """
        self._compute_speeds(max_distance, max_angle_diff, limit, amount_of_rods)
        output = []
        for index in range(len(self._angular_speeds)):
            num_of_rods = len(self._angular_speeds[index])
            output.append(0)
            for angular_speed in list(self._angular_speeds[index].values()):
                output[index] += angular_speed**2/num_of_rods
        return output

    def local_speeds(self, max_distance=100, max_angle_diff=90, limit=5,
                            amount_of_rods=200, divisions=5):
        """
        Returns local_speeds array.
        """
        self._compute_speeds(max_distance, max_angle_diff, limit, amount_of_rods)
        if len(self._local_speeds) == 0:
            self._compute_local_speeds(divisions)
        return self._local_speeds

    def _compute_local_speeds(self, divisions):
        """
        Creates an array of matrices. Each matrix's entry is a dictionariy such
        as {rod_id: (speed, angular_speed)}
        [index1_loc, index2_loc...]
            index1 = compress([subsys11_dic, subsys12_dic...,
                               subsys21_dic, subsys22_dic...])
                subsys1_dic = {rod_id: (speed, angular_speed)}
        """
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Computing local speeds"
        output_queue = mp.Queue()
        processes = []
        for index in range(len(self._evolution_dictionaries)-1):
            process = mp.Process(target=self._compute_local_speeds_process,
                                args=(index, output_queue, divisions))
            self._local_speeds.append(None)
            processes.append(process)
        num_processes = len(processes)
        running, processes_left = methods.run_processes(processes)
        finished = 0
        previous_time = datetime.datetime.now()
        counter = 0
        time_left = None
        times = []
        print " "
        while finished < num_processes:
            counter += 1
            finished += 1
            previous_time, counter, time_left = methods.print_progress(finished, num_processes,
                                    counter, times, time_left, previous_time)
            output = output_queue.get()
            index = output[0]
            speeds_matrix = output[1]
            self._local_speeds[index] = speeds_matrix
            if len(processes_left):
                new_process = processes_left.pop(0)
                time.sleep(settings.waiting_time)
                new_process.start()
        print CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE

    def _compute_local_speeds_process(self, index, output_queue, divisions):
        """
        Process
        """
        state = self.get(index)
        subgroups_matrix = state.subgroups_matrix(divisions)
        speeds_matrix = []
        speeds = methods.decompress(self._speeds[index], level=settings.medium_comp_level)
        angular_speeds = methods.decompress(self._angular_speeds[index], level=settings.medium_comp_level)
        for row in subgroups_matrix:
            speeds_row = []
            for subsystem in row:
                subsystem_dict = {}
                for rod in subsystem:
                    rod_id = rod.identifier
                    try:
                        speed = speeds[rod_id]
                        angular_speed = angular_speeds[rod_id]
                        subsystem_dict[rod_id] = (speed, angular_speed)
                    except KeyError:
                        pass
                speeds_row.append(subsystem_dict)
            speeds_matrix.append(speeds_row)
        output_queue.put([index, methods.compress(speeds_matrix, level=settings.medium_comp_level)])

    def _compute_local_average_speeds(self, max_distance=100, max_angle_diff=90,
                                     limit=5, amount_of_rods=200, divisions=5):
        """
        Compute local average speeds.
        """
        if not len(self._local_average_quadratic_speeds):
            print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Computing local average speeds"
            output_queue = mp.Queue()
            processes = []
            local_speeds = self.local_speeds(max_distance, max_angle_diff,
                                            limit, amount_of_rods, divisions)
            for index in range(len(self._evolution_dictionaries)-1):
                process = mp.Process(target=compute_local_average_speeds_process,
                                    args=(index, output_queue, local_speeds, divisions))
                self._local_average_quadratic_speeds.append(None)
                self._local_average_quadratic_angular_speeds.append(None)
                processes.append(process)
            num_processes = len(processes)
            running, processes_left = methods.run_processes(processes)
            finished = 0
            previous_time = datetime.datetime.now()
            counter = 0
            time_left = None
            times = []
            print " "
            while finished < num_processes:
                counter += 1
                finished += 1
                previous_time, counter, time_left = methods.print_progress(finished, num_processes,
                                        counter, times, time_left, previous_time)
                output = output_queue.get()
                index = output[0]
                self._local_average_quadratic_speeds[index] = output[1]
                self._local_average_quadratic_angular_speeds[index] = output[2]
                if len(processes_left):
                    new_process = processes_left.pop(0)
                    time.sleep(settings.waiting_time)
                    new_process.start()
            print CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE

    def local_average_quadratic_speed(self, max_distance=100, max_angle_diff=90,
                                        limit=5, amount_of_rods=200, divisions=5):
        """
        Returns a array of matrices. Each matrix represents
        local average quadratic speed values.
        """
        self._compute_local_average_speeds(max_distance, max_angle_diff, limit,
                                            amount_of_rods, divisions)
        return self._local_average_quadratic_speeds

    def local_average_quadratic_angular_speed(self, max_distance=100,
                                        max_angle_diff=90, limit=5,
                                        amount_of_rods=200, divisions=5):
        """
        Returns a array of matrices. Each matrix represents
        local average quadratic angular speed values.
        """
        self._compute_local_average_speeds(max_distance, max_angle_diff,
                                            limit, amount_of_rods, divisions)
        return self._local_average_quadratic_angular_speeds

    def density_and_quad_speed(self, max_distance=100, max_angle_diff=90,
                                limit=5, amount_of_rods=200, divisions=5):
        """
        Returns 2 arrays: density values and temperature
        """
        if (max_distance, max_angle_diff, limit, amount_of_rods) != (self._max_distance,
                    self._max_angle_diff, self._limit, self._amount_of_rods):
            self._reset()
        if not len(self._densities_array):
            laqs = self.local_average_quadratic_speed
            quad_speeds_array = laqs(max_distance=max_distance,
                                max_angle_diff=max_angle_diff,
                                limit=limit, amount_of_rods=amount_of_rods,
                                divisions=divisions)
            laqas = self.local_average_quadratic_angular_speed
            ang_speeds_array = laqas(max_distance=max_distance,
                                    max_angle_diff=max_angle_diff,
                                    limit=limit, amount_of_rods=amount_of_rods,
                                    divisions=divisions)
            output_queue = mp.Queue()
            processes = []
            for index in range(len(self)-1):
                process = mp.Process(target=self._density_and_quad_speed_process,
                            args=(index, output_queue, quad_speeds_array,
                                  ang_speeds_array, divisions))
                processes.append(process)
            num_processes = len(processes)
            running, processes_left = methods.run_processes(processes)
            finished = 0
            previous_time = datetime.datetime.now()
            counter = 0
            time_left = None
            times = []
            print " "
            while finished < num_processes:
                counter += 1
                finished += 1
                previous_time, counter, time_left = methods.print_progress(finished, num_processes,
                                        counter, times, time_left, previous_time)
                output = output_queue.get()
                densities.append(output[0])
                quad_speeds.append(output[1])
                self._speeds_matrices[output[3]](output[2])
                if len(processes_left):
                    new_process = processes_left.pop(0)
                    time.sleep(settings.waiting_time)
                    new_process.start()
            self._densities_array = densities
            self._quad_speeds_array = quad_speeds
            print(CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE)
        return [self._densities_array, self._quad_speeds_array]

    def _density_and_quad_speed_process(self, index, output_queue,
                                    quad_speeds_array, ang_speeds_array, divisions):
        """
        Process
        """
        quad_speeds = methods.decompress(quad_speeds_array[index], level=settings.medium_comp_level)
        ang_speeds = methods.decompress(ang_speeds_array[index], level=settings.medium_comp_level)
        state = self.get(index)
        subgroups = state.subgroups_matrix(divisions)
        densities = []
        speeds = []
        speeds_matrix = []
        for row in range(len(quad_speeds)):
            speeds_row = []
            for col in range(len(quad_speeds[row])):
                subgroup = subgroups[row][col]
                num_rods = subgroup.number_of_rods
                if num_rods > 0:
                    quad_speed = quad_speeds[row][col]
                    ang_speed = ang_speeds[row][col]
                    density = subgroup.density
                    total_quad_speed = float(ang_speed + quad_speed)/num_rods
                else:
                    density = 0
                    total_quad_speed = 0
                speeds_row.append(total_quad_speed)
                densities.append(density)
                speeds.append(total_quad_speed)
            speeds_matrix.append(speeds_row)
        output_queue.put([densities, speeds, speeds_matrix, index])

    def divide_systems_in_circles(self, divisions=5):
        """
        Divides all systems in circles.
        """
        if divisions != self._divisions and not self._done:
            print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Dividing systems in circles"
            self._done = True
            self._divisions = divisions
            processes = []
            output_queue = mp.Queue()
            states = []
            for index in range(len(self)):
                states.append(None)
                process = mp.Process(target=self.divide_system_in_circles_process,
                                     args=(divisions, index, output_queue))
                processes.append(process)
            num_processes = len(processes)
            running, processes_left = methods.run_processes(processes)
            finished = 0
            previous_time = datetime.datetime.now()
            counter = 0
            time_left = None
            times = []
            print " "
            while finished < num_processes:
                counter += 1
                finished += 1
                previous_time, counter, time_left = methods.print_progress(finished, num_processes,
                                        counter, times, time_left, previous_time)
                output = output_queue.get()
                index = output[0]
                state = output[1]
                self._states[index] = state
                if len(processes_left):
                    new_process = processes_left.pop(0)
                    time.sleep(settings.waiting_time)
                    new_process.start()
            print(CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE)

    def divide_system_in_circles_process(self, divisions, index, output_queue):
        """
            Process
        """
        state = self.get(index)
        state.divide_in_circles(divisions)
        state = methods.compress(state, level=methods.settings.medium_comp_level)
        output_queue.put([index, state])

    def create_density_video(self, divisions, folder, fps,
                                 number_of_bursts):
        """
        Creates a video of density's evolution.
        """
        
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Creating densities video"
        frames = len(self)
        function_name = 'plottable_density_matrix'
        kappas = self.kappas
        prop = self.average_covered_area_proportion[0]
        name = str(folder)+str(function_name)+"_K"+str(kappas)+".mp4"#+"prop"+str(round(100*prop,1))+'%.mp4'
        units = "Normalized occupied area [S.U.]"
        z_max, z_min = self._generic_scatter_animator(name, function_name, units,
                        divisions, fps=fps, number_of_bursts=number_of_bursts)
        self._min_density = z_min
        self._max_density = z_max

    @property
    def kappas(self):
        """
        Kappas of systems.
        """
        if not self._kappas:
            state = self.get(0)
            self._kappas = state.kappas
        return self._kappas

    def create_relative_g2_video(self, divisions, folder, fps,
                                 number_of_bursts):
        """
        Creates a video of correlation g2 evolution.
        """
        
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Creating g2 video"
        frames = len(self)
        function_name = 'correlation_g2_plot_matrix'
        kappas = self.kappas
        prop = self.average_covered_area_proportion[0]
        name = str(folder)+str(function_name)+"_K"+str(kappas)+".mp4"#+"prop"+str(round(100*prop,1))+'%.mp4'
        units = "g2 [S.U.]"
        self._generic_scatter_animator(name, function_name, units,
                        divisions, fps=fps, number_of_bursts=number_of_bursts)

    def create_relative_g4_video(self, divisions, folder, fps,
                                 number_of_bursts):
        """
        Creates a video of correlation g4 evolution.
        """
        
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Creating g4 video"
        frames = len(self)
        function_name = 'correlation_g4_plot_matrix'
        kappas = self.kappas
        prop = self.average_covered_area_proportion[0]
        name = str(folder)+str(function_name)+"_K"+str(kappas)+".mp4"#+"prop"+str(round(100*prop,1))+'%.mp4'
        units = "g4 [S.U.]"
        self._generic_scatter_animator(name, function_name, units,
                        divisions, fps=fps, number_of_bursts=number_of_bursts)

    def create_average_angle_video(self, divisions, folder, fps,
                                 number_of_bursts):
        """
        Creates a video of average angle evolution.
        """
        
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Creating average angle video"
        frames = len(self)
        function_name = 'plottable_average_angle_matrix'
        kappas = self.get(0).kappas
        prop = self.average_covered_area_proportion[0]
        name = str(folder)+str(function_name)+"_K"+str(kappas)+".mp4"#+"prop"+str(round(100*prop,1))+'%.mp4'
        units = "Average angle [grad]"
        self._generic_scatter_animator(name, function_name, units,
                        divisions, fps=fps, number_of_bursts=number_of_bursts)

    @property
    def radius(self):
        """
            Returns radius of systems (first)
        """
        return self.get(1).radius

    def _generic_scatter_animator(self, name, function_name, units,
                                    divisions, fps=15, number_of_bursts=10):
        """
        Generic animator
        """
        
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Plotting"
        groups = self.bursts_groups
        groups = copy.deepcopy(groups)
        bursts_ = len(groups)
        x_val, y_val, z_vals_avg, z_max, z_min = self.get_z_vals(groups, bursts_, function_name, divisions)
        frames = len(z_vals_avg)
        if settings.plot:
            methods.create_scatter_animation(x_val, y_val, z_vals_avg, divisions, z_max, z_min, units, name, self.radius)
        if settings.to_file:
            output_file_name = name + ".data"
            output_file = open(output_file_name, 'w')
            data = methods.compress([x_val, y_val, z_vals_avg, divisions, z_max, z_min, units, name, self.radius], level=9)
            output_file.write(data)
            output_file.close()
        return z_max, z_min

    def _generic_scatter_animator_process(self, divisions, index, output_queue, function_name):
        """
        Process
        """
        state = self.get(index)
        function = getattr(state, function_name)
        values = function(divisions)
        output_queue.put([index, values[0], values[1], values[2], methods.compress(state,
                level=settings.medium_comp_level)])

    def get_z_vals(self, groups, bursts_, function_name, divisions):
        """
            Get z values for function.
        """
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Getting function values: [" + str(function_name) +"]"
        z_vals = []
        z_vals_avg = []
        x_val = []
        y_val = []
        match1 = re.match(r'.*plottable_density.*', function_name)
        match2 = re.match(r'.*g[2|4].*', function_name)
        print ""
        if match2:
            z_max = 1
            z_min = -1
        elif match1:
            z_max = 1
            z_min = 0
        else:
            z_maxs = []
            z_mins = []
        finished_ = 0
        previous_time = datetime.datetime.now()
        times = []
        time_left = None
        counter = 0
        for group in groups:
            finished_ += 1
            counter += 1
            previous_time, counter, time_left = methods.print_progress(finished_, bursts_, counter, times, time_left, previous_time)
            output_queue = mp.Queue()
            processes = []
            for index in group:
                process = mp.Process(target=self._generic_scatter_animator_process,
                                     args=(divisions, index, output_queue, function_name))
                processes.append(process)
            num_processes = len(processes)
            running, processes_left = methods.run_processes(processes)
            finished__ = 0
            results = []
            while finished__ < num_processes:
                finished__ += 1
                output = output_queue.get()
                index = output[0]
                x_val = output[1]
                y_val = output[2]
                z_val = output[3]
                assert type(output[4]) == type("string"), "El estado tiene que ir comprimido"
                self._states[index] = output[4]
                if not (match2 or match1):
                    z_maxs.append(max(z_val))
                    z_mins.append(min(z_val))
                z_vals.append(z_val)
                if len(processes_left):
                    new_process = processes_left.pop(0)
                    time.sleep(settings.waiting_time)
                    new_process.start()
            z_vals_avg.append(methods.compress(methods.array_average(z_vals),
                                           level=settings.medium_comp_level))
        print CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE
        if not (match2 or match1):
            z_max = max(z_maxs)
            z_min = min(z_mins)
        return x_val, y_val, z_vals_avg, z_max, z_min

    def _get_image_ids(self):
        """
        Creates a dictionary with image ids referenced with indices.
        """
        if not len(self._image_id_by_index):
            print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Creating dictionary image index -> image id"
            processes = []
            self._image_ids_done = True
            output_queue = mp.Queue()
            self._image_id_by_index = {}
            for index in range(len(self)):
                process = mp.Process(target=self._get_image_ids_process, args=(index, output_queue))
                processes.append(process)
            num_processes = len(processes)
            running, processes_left = methods.run_processes(processes)
            finished = 0
            previous_time = datetime.datetime.now()
            counter = 0
            time_left = None
            times = []
            print ""
            while finished < num_processes:
                counter += 1
                finished += 1
                previous_time, counter, time_left = methods.print_progress(finished, num_processes, counter,
                                            times, time_left, previous_time)
                [index, output] = output_queue.get()
                self._image_id_by_index[index] = output
                if len(processes_left):
                    new_process = processes_left.pop(0)
                    time.sleep(settings.waiting_time)
                    new_process.start()
            print CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE

    def _get_image_ids_process(self, index, output_queue):
        """
        Process
        """
        state = self.get(index)
        output = methods.get_number_from_string(state.id_string)
        state = None
        output_queue.put([index, output])

    def _get_image_id(self, index):
        """
        Returns consecutive system images' ids.
        """
        if not self._image_ids_done:
            self._get_image_ids()
            self._image_ids_done = True
        return self._image_id_by_index[index]

    def create_videos(self, divisions=5, folder="./", fps=1,
                            max_distance=100, max_angle_diff=90, limit=5,
                            amount_of_rods=200, number_of_bursts=1, only_density=False):
        """
        Creates a video per property of the system that shows evolution.
        """
        self.divide_systems_in_circles(divisions)
        self.create_density_video(divisions, folder, fps, number_of_bursts)
        if not only_density:
            self.create_relative_g2_video(divisions, folder, fps, number_of_bursts)
            self.create_relative_g4_video(divisions, folder, fps, number_of_bursts)
        if not settings.ignore_temperature:
            self.create_temperature_video(divisions, folder, fps, max_distance,
                                   max_angle_diff, limit, amount_of_rods,
                                   number_of_bursts)

    def plottable_local_average_quadratic_speeds(self,
                                        max_distance=100,
                                        max_angle_diff=90, limit=5,
                                        amount_of_rods=200, divisions=5):
        """
        Returns plotable data.
        """
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Creating plottable matrix"
        quad_speeds = self.local_average_quadratic_speed(max_distance,
                                        max_angle_diff, limit,
                                        amount_of_rods, divisions)
        x_vals, y_vals, z_vals = [], [], []
        output_queue = mp.Queue()
        processes = []
        output_queue = mp.Queue()
        for index in range(len(self)-1):
            process = mp.Process(target=self.plottable_local_average_quadratic_speeds_process,
                                 args=(index, quad_speeds, max_distance, max_angle_diff, limit,
                                        amount_of_rods, divisions, output_queue))
            processes.append(process)
        num_processes = len(processes)
        running, processes_left = methods.run_processes(processes)
        finished = 0
        previous_time = datetime.datetime.now()
        counter = 0
        time_left = None
        times = []
        print " "
        while finished < num_processes:
            counter += 1
            finished += 1
            previous_time, counter, time_left = methods.print_progress(finished, num_processes,
                                    counter, times, time_left, previous_time)
            output = output_queue.get()
            x_vals.append(output[0])
            y_vals.append(output[1])
            z_vals.append(output[2])
            if len(processes_left):
                new_process = processes_left.pop(0)
                time.sleep(settings.waiting_time)
                new_process.start()
        print CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE
        return x_vals, y_vals, z_vals

    def plottable_local_average_quadratic_speeds_process(self, index, quad_speeds,
                                        max_distance, max_angle_diff, limit,
                                        amount_of_rods, divisions, output_queue):
        """
        Process
        """
        state = self.get(index)
        subgroups = state.subgroups_matrix(divisions)
        x_val, y_val, z_val = [], [], []
        for row_index in range(len(subgroups)):
            for col_index in range(len(subgroups[row_index])):
                subgroup = subgroups[row_index][col_index]
                quad_speed = methods.decompress(quad_speeds[index], level=settings.medium_comp_level)[row_index][col_index]
                center = subgroup.center
                center_x = center[0]
                center_y = center[1]
                x_val.append(center_x)
                y_val.append(center_y)
                z_val.append(quad_speed)
        output_queue.put([x_val, y_val, z_val])

    def create_cluster_histogram_video(self, max_distance=None,
                                    max_angle_diff=None, fps=15):
        """
            Creates a video of cluster length histogram.
        """
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Creating cluster histogram animation"
        kappa = self.kappas
        prop = self.average_covered_area_proportion[0]
        name = "cluster_hist_K"+str(int(kappa))+"prop"+str(round(100*prop,1))+"%.mp4"
        processes = []
        output_queue = mp.Queue()
        arrays = []
        fig = plt.figure()
        for index in range(len(self)):
            process = mp.Process(target=self.create_cluster_hist_video_process,
                                 args=(index, max_distance, max_angle_diff,
                                       output_queue))
            processes.append(process)
            arrays.append(None)
        num_processes = len(processes)
        running, processes_left = methods.run_processes(processes)
        finished = 0
        previous_time = datetime.datetime.now()
        counter = 0
        time_left = None
        times = []
        print " "
        while finished < num_processes:
            counter += 1
            finished += 1
            previous_time, counter, time_left = methods.print_progress(finished, num_processes,
                                    counter, times, time_left, previous_time)
            output = output_queue.get()
            index = output[0]
            array = output[1]
            arrays[index] = array
            if len(processes_left):
                new_process = processes_left.pop(0)
                time.sleep(settings.waiting_time)
                new_process.start()
        def animate(dummy_frame):
            """
            Animation function.
            """
            self._cluster_video_wrapper(arrays)
        frames = len(self)
        anim = animation.FuncAnimation(fig, animate, frames=frames)
        anim.save(name, writer=methods.WRITER, fps=fps)

    def _cluster_video_wrapper(self, arrays):
        """
        Wrapper.
        """
        try:
            plt.cla()
            plt.clf()
            array = arrays.pop(0)
            boundaries = [0+time*1 for time in range(400)]
            plt.xlim((0, 1000))
            plt.ylim((0, 1))
            plt.hist(array, bins=boundaries, normed=True)
            plt.suptitle("Cluster length histogram")
            plt.xlabel("Number of rods in cluster")
            plt.ylabel("Number of clusters (normalized)")
        except:
            pass

    def create_cluster_hist_video_process(self, index, max_distance, max_angle_diff,
                                        output_queue):
        """
        Process
        """
        state = self.get(index)
        array = state.cluster_lengths(max_distance=max_distance,
                                            max_angle_diff=max_angle_diff,
                                            min_size=1)
        for index2 in range(len(array)):
            array[index2] *= index2
        output_queue.put([index, array])
        return

    def create_speeds_vectors_video(self, divisions, folder, fps, max_distance,
                                 max_angle_diff, limit, amount_of_rods,
                                 number_of_bursts):
        """
            Creates a video of average speed vectors over subsystem.
        """
        vectors_matrices = self.average_speeds_vectors(
                                            divisions, max_distance,
                                            max_angle_diff)
        fig = plt.figure()
        kappa = self.kappas
        name = str(folder) + "speeds_vectors_K" + str(kappa)
        bursts_groups = copy.deepcopy(self.bursts_groups)
        end = False
        vectors_matrices_avg = []
        def animate(dummy_frame):
            """
            Animation function.
            """
            self._speeds_vectors_video_wrapper(vectors_matrices)    #XXX BUG
        frames = len(self)
        anim = animation.FuncAnimation(fig, animate, frames=frames)
        anim.save(name, writer=methods.WRITER, fps=fps)

    def _speeds_vectors_video_wrapper(self, vectors_matrix):
        """
            Wrapper
        """
        pass

    def average_speeds_vectors(self, divisions, max_distance, max_angle_diff):
        """
            Computes average speeds vectors
        """
        output_queue = mp.Queue()
        speeds = self.speeds_vectors
        processes = []
        vector_matrices = []
        speeds_num = len(self)-1
        for index in range(speeds_num):
            speeds_ = speeds[index]
            compressed_state = self._states[index]
            process = mp.Process(target=average_speeds_vectors_video_process,
                                 args=(divisions, index, compressed_state,
                                     speeds_, output_queue))
            processes.append(process)
            vector_matrices.append(None)
        num_processes = len(processes)
        running, processes_left = methods.run_processes(processes)
        finished = 0
        previous_time = datetime.datetime.now()
        counter = 0
        time_left = None
        times = []
        print " "
        while finished < num_processes:
            counter += 1
            finished += 1
            previous_time, counter, time_left = methods.print_progress(finished, num_processes,
                                    counter, times, time_left, previous_time)
            output = output_queue.get()
            index = output[0]
            vector_matrix = output[1]
            vector_matrices[index] = vector_matrix
            if len(processes_left):
                new_process = processes_left.pop(0)
                time.sleep(settings.waiting_time)
                new_process.start()
        vectors = []
        for index in range(speeds_num):
            state = self.get(index)
            vector_matrix = vector_matrices[index]

    def create_temperature_video(self, divisions, folder, fps,
                            max_distance, max_angle_diff,
                            limit, amount_of_rods, number_of_bursts):
        """
        Creates a video of temperature evolution.
        """
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Creating temperature video"
        x_vals, y_vals, z_vals = self.plottable_local_average_quadratic_speeds(
                                        max_distance, max_angle_diff, limit,
                                        amount_of_rods, divisions)
        x_vals = x_vals[0]
        y_vals = y_vals[0]
        bursts_groups = copy.deepcopy(self.bursts_groups)
        end = False
        z_vals_avg = []
        number_of_bursts *= 5
        print "--"*(len(inspect.stack())-1)+">"+"["+str(inspect.stack()[0][3])+"]: " + "Computing averages"
        while not end:
            groups = []
            average = None
            _z_vals = []
            for dummy_time in range(number_of_bursts):
                try:
                    group = bursts_groups.pop(0)
                    groups.append(group)
                except IndexError:
                    end = True
            if not len(groups):
                break
            for group in groups:
                for dummy_time in range(len(group)):
                    z_val = z_vals.pop(0)
                    new_vals = []
                    for val in z_val:
                        try:
                            new_vals.append(val*float(self.kappas))
                        except TypeError:
                            new_vals.append(None)
                    z_val = new_vals
                    _z_vals.append(z_val)
            try:
                average = methods.array_average(_z_vals)
            except IndexError:
                average = _z_vals
            z_vals_avg.append(methods.compress(average, level=settings.medium_comp_level))
        if not settings.to_file:
            fig = plt.figure()
        state = self.get(0)
        kappas = state.kappas
        state = None
        name = str(folder)+"Temperature"+str(kappas)+".mp4"
        z_maxs = []
        z_mins = []
        try:
            for z_val_ in z_vals_avg:
                z_val = methods.decompress(z_val_, level=settings.medium_comp_level)
                z_maxs.append(max(z_val))
                z_mins.append(min(z_val))
        except ValueError:
            print(z_vals_avg)
            raise ValueError
        z_max = max(z_maxs)
        z_min = min(z_mins)
        units = "Temperature [pixels^2/seg^2]"
        rad = self.radius
        if settings.plot:
            methods.create_scatter_animation(x_vals, y_vals, z_vals_avg, divisions, z_max, z_min, units, name, radius=rad)
        if settings.to_file:
            output_file_name = name + ".data"
            output_file = open(output_file_name, 'w')
            data = methods.compress([x_vals, y_vals, z_vals_avg, divisions, z_max, z_min, units, name, self.radius], level=9)
            output_file.write(data)
            output_file.close()

    def _temperature_video_wrapper(self, x_val, y_val, z_vals,
                                divisions, name, z_max, z_min):
        """
        Wrapper
        """
        try:
            z_val = methods.decompress(z_vals.pop(0), settings.medium_comp_level)
        except IndexError:
            return
        plt.cla()
        plt.clf()
        rad = float(self.radius*1.3)/divisions
        size = (rad/4)**2
        x_min = min(x_val)-rad*1.1
        x_max = max(x_val)+rad*1.1
        y_min = min(y_val)-rad*1.1
        y_max = max(y_val)+rad*1.1
        plt.xlim((x_min, x_max))
        plt.ylim((y_min, y_max))
        plt.scatter(x_val, y_val, s=size, c=z_val, marker='s',
                    vmax=z_max, vmin=z_min)
        plt.gca().invert_yaxis()
        cb = plt.colorbar()
        plt.xlabel("x [pixels]")
        plt.ylabel("y [pixels]")
        cb.set_label("Temperature [pixels^2/seg^2]")

    def _average_cluster_areas(self, z_vals, number_of_bursts=1,
                               min_size=3):
        """
        Wrapper
        """
        bursts_groups = copy.deepcopy(self.bursts_groups)
        end = False
        z_vals_avg = []
        index = 0
        indices = []
        while not end:
            initial_id = self._get_image_id(index)
            indices.append(index)
            groups = []
            average = None
            _z_vals = []
            for dummy_time in range(number_of_bursts):
                try:
                    group = bursts_groups.pop(0)
                    groups.append(group)
                except IndexError:
                    end = True
            if not len(groups):
                break
            for group in groups:
                for dummy_time in group:
                    index += 1
                    z_val = z_vals.pop(0)
                    _z_vals.append(z_val)
            try:
                average = sum(_z_vals)/float(len(_z_vals))
            except TypeError:
                not_none_vals = []
                for val in _z_vals:
                    if val is not None:
                        not_none_vals.append(val)
                try:
                    average = sum(not_none_vals)/float(len(not_none_vals))
                except ZeroDivisionError:
                    average = 0
            z_vals_avg.append(average)
        return z_vals_avg, indices

    def _get_cluster_areas(self, max_distance=None,
                    max_angle_diff=None, min_size=3):
        """
        Compute cluster areas for all states.
        """        
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Getting cluster areas"
        output_queue = mp.Queue()
        processes = []
        areas = []
        for index in range(len(self)):
            process = mp.Process(target=self._get_cluster_areas_process,
                                 args=(index, max_distance,
                                    max_angle_diff, output_queue, min_size))
            processes.append(process)
            areas.append(None)
        num_processes = len(processes)
        running, processes_left = methods.run_processes(processes)
        finished = 0
        previous_time = datetime.datetime.now()
        counter = 0
        time_left = None
        times = []
        print " "
        while finished < num_processes:
            counter += 1
            finished += 1
            previous_time, counter, time_left = methods.print_progress(finished, num_processes,
                                    counter, times, time_left, previous_time)
            output = output_queue.get()
            index = output[0]
            area = output[1]
            areas[index] = area
            if len(processes_left):
                new_process = processes_left.pop(0)
                time.sleep(settings.waiting_time)
                new_process.start()
        print CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE
        return areas

    def _get_cluster_areas_process(self, index, max_distance,
                    max_angle_diff, output_queue, min_size):
        """
        Process
        """
        state = self.get(index)
        area = state.total_area_of_clusters(max_distance=max_distance,
                            max_angle_diff=max_angle_diff, min_size=min_size)
        output_queue.put([index, area])
        return

    def cluster_areas(self, number_of_bursts=1, max_distance=None,
                    max_angle_diff=None, min_size=5):
        """
        Returns an array where each value is total cluster area in
        that moment.
        """
        if max_distance and max_angle_diff:
            self._reset()
        elif not max_distance and not max_angle_diff:
            pass
        else:
            msg = "Both or none to be defined: max_distance, max_angle_diff"
            raise ValueError(msg)
        if not len(self._cluster_areas):
            z_vals = self._get_cluster_areas(max_distance=max_distance,
                            max_angle_diff=max_angle_diff, min_size=min_size)
            areas, indices = self._average_cluster_areas(z_vals,
                                        number_of_bursts=number_of_bursts,
                                        min_size=min_size)
            self._cluster_areas = areas
            z_vals = self._get_cluster_areas(max_distance=max_distance,
                            max_angle_diff=max_angle_diff, min_size=0)
            norm_areas, indices = self._average_cluster_areas(z_vals,
                                        number_of_bursts=number_of_bursts,
                                        min_size=0)
            self._total_cluster_areas = norm_areas
            self._indices = indices
        return self._cluster_areas

    def get_order_evolution_coeficient(self, number_of_bursts=1, max_distance=None,
                    max_angle_diff=None, min_size=5):
        """
            Returns coeficient of order param evolution.
        """
        if self._popt is None or self._std_dev is None:
            print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Computing order evolution parameter"
            cluster_areas = self.cluster_areas(number_of_bursts=number_of_bursts,
                                       max_distance=max_distance,
                                       max_angle_diff=max_angle_diff,
                                       min_size=min_size)
            indices = self._indices
            self._compute_times(number_of_bursts)
            times = self._times
            times.pop(0)
            cluster_areas.pop(0)
            log_areas = []
            log_times = []
            for index in range(len(cluster_areas)):
                try:
                    log_areas.append(math.log(cluster_areas[index]))
                except ValueError:
                    continue
                try:
                    log_times.append(math.log(times[index]))
                except ValueError:
                    log_times.pop()
            top = 0
            bot = 0
            try:
                x_m = float(sum(log_times))/len(log_times)  
                y_m = float(sum(log_areas))/len(log_areas)
            except ZeroDivisionError:
                print "--"*(len(inspect.stack())-1)+">"+"["+str(inspect.stack()[0][3])+"]: Not enough data"
                return
            for index in range(len(log_areas)):
                top += (log_times[index]-x_m)*(log_areas[index]-y_m)
                bot += (log_times[index]-x_m)**2
            m = float(top)/bot
            b = y_m - m*x_m
        return m, b

    def plot_cluster_areas(self, number_of_bursts=1, max_distance=None,
                    max_angle_diff=None, min_size=10):
        """
            Plots cluster areas evolution.
        """
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Computing cluster areas"
        areas = self.cluster_areas(number_of_bursts=number_of_bursts,
                                   max_distance=max_distance,
                                   max_angle_diff=max_angle_diff,
                                   min_size=min_size)
        total_areas = self._total_cluster_areas
        assert total_areas, "Error"
        norm_areas = []
        self._compute_times(number_of_bursts=number_of_bursts)
        times = []
        for index in range(len(areas)):
            area = areas[index]
            total_area = total_areas[index]
            try:
                proportion = float(area)/total_area
            except ZeroDivisionError:
                proportion = 0
                continue
            times.append(self._times[index])
            norm_areas.append(proportion)
        if not settings.to_file:
            fig = plt.figure()
            plt.xlabel("time[seconds]")
            plt.ylabel("cluster area proportion")
        m, b = self.get_order_evolution_coeficient(number_of_bursts=number_of_bursts,
                                            max_distance=max_distance,
                                            max_angle_diff=max_angle_diff,
                                            min_size=min_size)
        line = []
        for time in times:
            if time != 0:
                line.append((math.e**b)*(time**m))
        name = "coef_K"+str(self.kappas)+".log"
        output = open(name, 'w')
        text = "Coeficient: "+str(m)+"\nIndep: "+str(b) +"\n"
        output.write(text)
        output.close()
        if not settings.to_file:
            plt.plot(times, line)
        else:
		    output_file_name = "linear_approx_K"+str(self.kappas)+".data"
		    output_file = open(output_file_name, 'w')
		    data = methods.compress([times, line], level=9)
		    output_file.write(data)
		    output_file.close()
        clust = open("cluster_areas.txt", "w")
        for index in range(len(times)):
            try:
                string = str(times[index])+"    "+str(norm_areas[index])+"\n"
                clust.write(string)
            except:
                pass
        clust.close()
        clust = open("cluster_areas_log.txt", "w")
        for index in range(len(times)):
            try:
                string = str(math.log(times[index]))+"    "+str(math.log(norm_areas[index]))+"\n"
                clust.write(string)
            except:
                pass
        clust.close()
        if not settings.to_file:
            try:
                plt.scatter(times, norm_areas)
            except ValueError:
                print(len(times), len(norm_areas))
                print(times, norm_areas)
                raise ValueError
            plt.grid(b=True, which='major', color='b', linestyle='-')
            plt.grid(b=True, which='minor', color='b', linestyle='--')
            plt.xscale('log')
            plt.yscale('log')
            file_name = "cluster_areas_K" + str(self.kappas) + ".png"
            plt.savefig(file_name)
        else:
	        output_file_name = "cluster_areas_K"+str(self.kappas)+".data"
	        output_file = open(output_file_name, 'w')
	        data = methods.compress([times, norm_areas], level=9)
	        output_file.write(data)
	        output_file.close()

    def _compute_times(self, number_of_bursts):
        """
            Computes times for experiment.
        """
        diff_t = 15*number_of_bursts
        id_0 = self._get_image_id(self._indices[0])
        times = []
        previous_time = 0
        for index_ in range(len(self._indices)):
            index = self._indices[index_]
            if not index:
                time = 0
            else:
                id_1 = self._get_image_id(index)
                index_0 = self._indices[index_-1]
                id_0 = self._get_image_id(index_0)
                time = methods.time_difference(self._dates, id_0, id_1)
            time += previous_time
            previous_time = time
            times.append(time)
        self._times = times

    @property
    def bursts_groups(self):
        """
        Returns a list of groups of indices that are in a row.
        """
        if not self._bursts_computed:
            self._bursts_computed = True
            print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "bursts_groups"
            groups = []
            group = []
            output_queue = mp.Queue()
            processes = []
            initial_id = self._get_image_id(0)
            for index in range(len(self._state_numbers)-1):
                final_id = self._get_image_id(index+1)
                date1 = self._dates[initial_id]
                date2 = self._dates[final_id]
                process = mp.Process(target=methods.are_in_burst_queue,
                                     args=(index, date1, date2, output_queue))
                processes.append(process)
                initial_id = final_id
            num_processes = len(processes)
            running, processes_left = methods.run_processes(processes)
            finished = 0
            previous_time = datetime.datetime.now()
            counter = 0
            time_left = None
            times = []
            results = []
            print ""
            while finished < num_processes:
                counter += 1
                finished += 1
                previous_time, counter, time_left = methods.print_progress(finished, num_processes,
                                    counter, times, time_left, previous_time)
                output = output_queue.get()
                index = output[0]
                burst = output[1]
                results.append([index, burst])
                if len(processes_left):
                    new_process = processes_left.pop(0)
                    new_process.start()
            for result in results:
                index = result[0]
                burst = result[1]
                if burst:
                    group.append(index)
                else:
                    group.append(index)
                    groups.append(group)
                    group = []
            self._bursts_groups = groups
            print(CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE)
        assert len(self._bursts_groups), "Fail"
        return self._bursts_groups

    def plot_average_temperature(self, max_distance, max_angle_diff, limit):
        """
            Average temperature over time.
        """
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Computing average temperature"

        self._compute_speeds(max_distance, max_angle_diff,
                        limit, None)
        indices = []
        average_speeds = []
        processes = []
        output_queue = mp.Queue()
        for index in range(len(self._speeds)):
            indices.append(index)
            process = mp.Process(target=self.average_speeds,
                                 args=(index, output_queue))
            processes.append(process)
            average_speeds.append(None)
        num_processes = len(processes)
        running, processes_left = methods.run_processes(processes)
        finished = 0
        previous_time = datetime.datetime.now()
        counter = 0
        time_left = None
        times = []
        results = []
        print ""
        while finished < num_processes:
            counter += 1
            finished += 1
            previous_time, counter, time_left = methods.print_progress(finished, num_processes,
                                counter, times, time_left, previous_time)
            output = output_queue.get()
            index = output[0]
            average_speed = output[1]
            average_speeds[index] = average_speed
            if len(processes_left):
                new_process = processes_left.pop(0)
                time.sleep(settings.waiting_time)
                new_process.start()
        print(CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE)
        if not settings.to_file:
            plt.figure()
            plt.plot(indices, average_speeds)
            name = "avg_temp_K"
            kappa = self.kappas
            name += str(kappa) + ".png"
            plt.xlabel("time [seconds]")
            plt.ylabel("average temperature [pixels^2 / s^2]")
            plt.savefig(name)
        else:
            name = "average_temperature_K"+str(self.kappas)+".data"
            output_file = open(name, 'w')
            data = methods.compress([indices, average_speeds], level=9)
            output_file.write(data)
            output_file.close()

    def average_speeds(self, index, output_queue):
        """
            Averages speeds of all rods.
        """
        speeds = methods.decompress(self._speeds[index], level=settings.medium_comp_level)
        number_of_rods = len(list(speeds.keys()))
        average_speed = 0
        for speed in list(speeds.values()):
            average_speed += speed*self.kappas
        average_speed /= number_of_rods
        output_queue.put([index, average_speed])

    @property
    def lost_rods_percentage(self):
        """
            Computes maximum percentage of rods lost in analysis.
        """
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Computing lost rods percentaje"
        number_of_rods = []
        for index in range(len(self)):
            state = self.get(index)
            number = state.number_of_rods
            number_of_rods.append(number)
        supposed_total_number = max(number_of_rods)
        loses = [float(supposed_total_number-number)/supposed_total_number
                 for number in number_of_rods]
        loses_2 = [lose**2 for lose in loses]
        loses_average = float(sum(loses))/len(loses)
        loses_2_average = float(sum(loses_2))/len(loses_2)
        std_dev = math.sqrt(loses_2_average-loses_average**2)
        return loses_average, std_dev

    @property
    def average_covered_area_proportion(self):
        """
            Returns average covered area proportion.
        """
        if not self._covered_area_prop:
            self._covered_area_prop = self._average_covered_area_proportion()
        return self._covered_area_prop

    def _average_covered_area_proportion(self):
        """
            Returns average covered area proportion.
        """
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Computing covered area proportion"
        props = []
        processes = []
        output_queue = mp.Queue()
        for index in range(len(self)):
            process = mp.Process(target=self._average_covered_area_prorportion_process,
                                 args=(index, output_queue))
            processes.append(process)
            props.append(None)
        num_processes = len(processes)
        running, processes_left = methods.run_processes(processes)
        finished = 0
        previous_time = datetime.datetime.now()
        counter = 0
        time_left = None
        times = []
        results = []
        print " "
        while finished < num_processes:
            counter += 1
            finished += 1
            previous_time, counter, time_left = methods.print_progress(finished, num_processes,
                                counter, times, time_left, previous_time)
            output = output_queue.get()
            index = output[0]
            prop = output[1]
            props[index] = prop
            if len(processes_left):
                new_process = processes_left.pop(0)
                time.sleep(settings.waiting_time)
                new_process.start()
        print(CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE)
        avg = float(sum(props))/len(props)
        std_dev_step1 = [(value-avg)**2 for value in props]
        std_dev_step2 = float(sum(std_dev_step1))/(len(std_dev_step1)-1)
        dev = math.sqrt(std_dev_step2)
        return avg, dev

    def _average_covered_area_prorportion_process(self, index, output_queue):
        """
        Process
        """
        state = self.get(index)
        prop = state.covered_area_proportion()
        output_queue.put([index, prop])

    def _compute_averages(self):
        """
        Computes experiment averages.
        """
        if not self._average_kappa:
            print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Computing experiment averages"
            kappas = []
            kappas_dev = []
            avg_rads = []
            avg_rod_lengths = []
            avg_rod_widths = []
            number_of_rods_ = []
            processes = []
            output_queue = mp.Queue()
            for index in range(len(self)):
                process = mp.Process(target=self._compute_averages_process,
                                     args=(index, output_queue))
                processes.append(process)
            num_processes = len(processes)
            running, processes_left = methods.run_processes(processes)
            finished = 0
            previous_time = datetime.datetime.now()
            counter = 0
            time_left = None
            times = []
            results = []
            print " "
            while finished < num_processes:
                counter += 1
                finished += 1
                previous_time, counter, time_left = methods.print_progress(finished, num_processes,
                                    counter, times, time_left, previous_time)
                [avg_kappa, kappa_dev, avg_rad,
                 avg_rod_length, avg_rod_width,
                 number_of_rods] = output_queue.get()
                kappas.append(avg_kappa)
                kappas_dev.append(kappa_dev)
                avg_rads.append(avg_rad)
                avg_rod_lengths.append(avg_rod_length)
                avg_rod_widths.append(avg_rod_width)
                number_of_rods_.append(number_of_rods)
                if len(processes_left):
                    new_process = processes_left.pop(0)
                    time.sleep(settings.waiting_time)
                    new_process.start()
            print(CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE)
            self._average_kappa = sum(kappas)*1.0/len(kappas)
            self._kappa_dev = sum(kappas_dev)*1.0/len(kappas_dev)
            self._average_rad = sum(avg_rads)*1.0/len(avg_rads)
            self._average_rod_length = sum(avg_rod_lengths)*1.0/len(avg_rod_lengths)
            self._average_rod_width = sum(avg_rod_widths)*1.0/len(avg_rod_widths)
            self._number_of_rods = sum(number_of_rods_)/float(len(number_of_rods_))
            std_dev_step1 = [(value-self._number_of_rods)**2 for value in number_of_rods_]
            std_dev_step2 = float(sum(std_dev_step1))/(len(std_dev_step1)-1)
            deviation = math.sqrt(std_dev_step2)
            self._number_of_rods_dev = deviation
        
    def _compute_averages_process(self, index, output_queue):
        """
        Process
        """
        state = self.get(index)
        avg_kappa = state.average_kappa
        kappa_dev = state.kappa_dev
        rad = state.radius
        rod_length = state.average_rod_length
        rod_width = state.average_rod_width
        number_of_rods = state.number_of_rods
        output_queue.put([avg_kappa, kappa_dev, rad, rod_length, rod_width, number_of_rods])

    @property
    def average_number_of_rods(self):
        """
            Returns average number of rods of all states.
        """
        self._compute_averages()
        return self._number_of_rods, self._number_of_rods_dev

    @property
    def average_kappa(self):
        """
            Returns average kappa of all states.
        """
        self._compute_averages()
        return self._average_kappa

    @property
    def average_kappa_dev(self):
        """
            Returns average kappa deviation of all states.
        """
        self._compute_averages()
        return self._kappa_dev


    @property
    def average_system_rad(self):
        """
            It returns average system rad
        """
        self._compute_averages()
        return self._average_rad

    @property
    def average_rod_length(self):
        """
            It returns average rod length.
        """
        self._compute_averages()
        return self._average_rod_length

    @property
    def average_rod_width(self):
        """
            It returns average rod width.
        """
        self._compute_averages()
        return self._average_rod_width

    def plot_rods(self, index):
        """
        Plot all rods of index-th system.
        """
        figure = plt.figure()
        fig = plt.gcf()
        state = self.get(index)
        x_vals, y_vals = state.rod_centers
        center_x, center_y = state.center
        rad = state.radius
        plt.scatter(center_x, center_y, color='b', s=10)
        plt.plot([center_x, center_x+rad],[center_y, center_y])
        plt.plot([center_x, center_x],[center_y, center_y+rad])
        plt.plot([center_x, center_x-rad],[center_y, center_y])
        plt.plot([center_x, center_x],[center_y, center_y-rad])
        plt.scatter(x_vals, y_vals, color='r')
        plt.show()


    def plot_number_of_rods_over_time(self):
        """
        Plot number of rods over time.
        """
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Creating graph... "
        images = [image for image in range(len(self))]
        number_of_rods = []
        output_queue = mp.Queue()
        processes = []
        for index in range(len(self)):
            process = mp.Process(target=self._plot_number_of_rods_process,
                                 args=(index, output_queue))
            processes.append(process)
            number_of_rods.append(None)
        num_processes = len(processes)
        running, processes_left = methods.run_processes(processes)
        finished = 0
        previous_time = datetime.datetime.now()
        counter = 0
        time_left = None
        times = []
        print " "
        while finished < num_processes:
            counter += 1
            finished += 1
            previous_time, counter, time_left = methods.print_progress(finished, num_processes,
                                counter, times, time_left, previous_time)
            [index, num_rods] = output_queue.get()
            number_of_rods[index] = num_rods
            if len(processes_left):
                new_process = processes_left.pop(0)
                #time.sleep(settings.waiting_time)
                new_process.start()
        print(CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE)
        figure = plt.figure()
        plt.plot(images, number_of_rods)
        name = "num_rods_time_K" + str(self.kappas) + ".png"
        plt.savefig(name)

    def _plot_number_of_rods_process(self, index, output_queue):
        """
        Process
        """
        state = self.get(index)
        num_rods = state.number_of_rods
        output_queue.put([index, num_rods])




def compute_local_average_speeds_process(index, output_queue, local_speeds, divisions):
    """
    Process
        [index1_loc, index2_loc...]
            index1_loc = compress([subsys1_dic, subsys2_dic...])
                subsys1_dic = {rod_id: (speed, angular_speed)}
    """
    speeds_matrix = []
    angular_speeds_matrix = []
    local_speeds_ = methods.decompress(local_speeds[index], level=settings.medium_comp_level)
    speeds_matrix = [[None for dummy_1 in range(divisions)] for dummy_2 in range(divisions)]
    angular_speeds_matrix = [[None for dummy_1 in range(divisions)] for dummy_2 in range(divisions)]
    for row in range(len(local_speeds_)):
        for col in range(len(local_speeds_[0])):
            subsys_dict = local_speeds_[row][col]
            subsys_quad_avg_speed = 0
            subsys_quad_avg_ang_speed = 0
            num_rods = len(subsys_dict)
            for speeds in subsys_dict.values():
                subsys_quad_avg_speed += float(speeds[0]**2)/num_rods
                subsys_quad_avg_ang_speed += float(speeds[1]**2)/num_rods
            speeds_matrix[row][col] = subsys_quad_avg_speed
            angular_speeds_matrix[row][col] = subsys_quad_avg_ang_speed
    output_queue.put([index, methods.compress(speeds_matrix, level=settings.medium_comp_level),
                    methods.compress(angular_speeds_matrix, level=settings.medium_comp_level)])

def average_speeds_vectors_video_process(divisions, index, compressed_state,
                                     speeds, output_queue):
    """
        Process
    """
    state = methods.decompress(compressed_state,
                        level=methods.settings.medium_comp_level)
    subgroups_matrix = state.subgroups_matrix(divisions)
    emp = [None, None]
    vectors_matrix = [[emp for dummy in range(len(subgroups_matrix[idx]))]
                      for idx in range(len(subgroups_matrix))]
    for row in range(len(subgroups_matrix)):
        for col in range(len(subgroups_matrix[row])):
            subgroup = subgroups_matrix[row][col]
            vector = [0, 0]
            number_of_rods = float(subgroup.number_of_rods)
            for rod_ in subgroup:
                rod_id = rod_.identifier
                rod_vector = list(speeds[rod_id])
                vector[0] += rod_vector[0]/number_of_rods
                vector[1] += rod_vector[1]/number_of_rods
            vectors_matrix[row][col] = vector
    output_queue.put([index, vectors_matrix])

