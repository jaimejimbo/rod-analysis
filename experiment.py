"""
    Library for time evolution study.
"""
import re, methods, math, copy, gc
import multiprocessing as mp
from matplotlib import animation
import matplotlib.pyplot as plt
import numpy
import numpy as np
import scipy.optimize as optimization



class Experiment(object):
    """
        Has a list of system states, one for each t.
    """
    def __init__(self, system_states_name_list=None,
                system_states_list=None, diff_t=1, dates=None):
        """
            Creation of experiment object.
        """

        type_ = str(type(system_states_list))
        if re.match(r'.*NoneType.*', type_):
            self._states = []
            self._number_of_states = 0
        elif re.match(r'.*list.*', type_):
            self._states = methods.compress(system_states_list)
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
            self._states_dict[number] = methods.compress(state)
        if not dates:
            msg = "Time evolution needs dates file.\n"
            msg = "To create one run export_image_dates().\n"
            raise ValueError(msg)

        self._diff_t = diff_t
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
        self._densities_array = []
        self._quad_speeds_array = []
        self._dates = dates
        self._speeds_matrices = []
        self._cluster_areas = []
        self._bursts_groups = []
        self._divisions = None
        self._min_density = None
        self._max_density = None
        self._indices = None
        self._total_cluster_areas = None
        #self._inversed = False
        _writer = animation.writers['ffmpeg']
        writer = _writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
        self._writer = writer

    def __iter__(self):
        """
            Magic method for looping.
        """
        for state in self._states:
            yield methods.decompress(state)


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


    def set_coef(self, value):
        """
            Changes coef for dividing in circles.
        """
        for state in self._states:
            state.coef = value

    def _reset(self):
        """
            Resets all variables, so they must be computed again.
        """
        self._evolution_dictionaries = []
        self._conflictive_final_rods = []
        self._relative_dictionaries = []
        self._unjoined_initial_rods = []
        self._unjoined_final_rods = []
        self._speeds_vectors = []
        self._final_rods = []
        self._initial_rods = []
        self._speeds = []
        self._angular_speeds = []
        self._max_distance = None
        self._max_angle_diff = None
        self._limit = None
        self._amount_of_rods = None
        self._local_speeds = []
        self._local_average_quadratic_speeds = []
        self._local_average_quadratic_angular_speeds = []
        self._densities_array = []
        self._quad_speeds_array = []
        self._speeds_matrices = []
        self._cluster_areas = []
        self._divisions = None
        self._min_density = None
        self._max_density = None
        self._indices = None
        self._total_cluster_areas = None


    def _create_dict_keys(self):
        """
            Create evolucion dictionaries keys.
        Each key is a rod's id.
        """
        for dummy_index in range(len(self._states)):
            self._evolution_dictionaries.append({})
            self._relative_dictionaries.append({})
            self._conflictive_final_rods.append(set([]))
            self._final_rods.append(set([]))
            self._initial_rods.append(set([]))
        for index in range(len(self._states)):
            state = self._states[index]
            if not state:
                msg = "State is not defined."
                raise TypeError(msg)
            evol_dict = self._evolution_dictionaries[index]
            relative_dict = self._relative_dictionaries[index]
            for rod in state:
                self._initial_rods[index] |= set([rod.identifier])
                rod_id = rod.identifier
                evol_dict[rod_id] = set([])
                relative_dict[rod_id] = {}
        for index in range(len(self._states)-1):
            self._final_rods[index] = self._initial_rods[index+1].copy()


    def _fill_dicts(self, max_distance, max_angle_diff,
                    limit=5, amount_of_rods=None):
        """
            Looks for rods that have only one possible predecessor.
        """
        (self._max_distance, self._max_angle_diff,
        self._limit, self._amount_of_rods) = (max_distance, max_angle_diff,
                                                limit, amount_of_rods)
        processes = []
        output_queue = mp.Queue()
        for index in range(len(self._states)-1):
            processes.append(mp.Process(target=self._fill_dicts_process_limited,
                                    args=(index, max_distance, max_angle_diff,
                                        output_queue, limit, amount_of_rods)))
        running, processes_left = methods.run_processes(processes)
        num_processes = len(running)
        finished = 0
        while finished < num_processes:
            finished += 1
            output_row = output_queue.get()
            index = output_row[0]
            self._evolution_dictionaries[index] = output_row[1]
            self._relative_dictionaries[index] = output_row[2]
            if len(processes_left):
                finished -= 1
                new_process = processes_left.pop(0)
                new_process.start()


    def _fill_dicts_process(self, index, max_distance, max_angle_diff,
                            output_queue, amount_of_rods=None):
        """
            Allows to create a process and use all cores.
        """
        initial_state = self._states[index]
        final_state = self._states[index+1]
        evol_dict = self._evolution_dictionaries[index]
        relative_dict = self._relative_dictionaries[index]
        for initial_rod in initial_state:
            initial_id = initial_rod.identifier
            if not amount_of_rods:
                available_final_rods = final_state._get_rods_range(0,
                                        final_state.number_of_rods)
            else:
                start_id = initial_id-amount_of_rods/2
                end_id = initial_id+amount_of_rods/2
                available_final_rods = final_state._get_rods_range(start_id,
                                                                    end_id)
            for final_rod in available_final_rods:
                final_id = final_rod.identifier
                distance = initial_rod.distance_to_rod(final_rod)
                angle = initial_rod.angle_between_rods(final_rod)
                angle = min([angle, 180-angle])
                speed = float(distance)/self._diff_t
                if distance <= max_distance and angle <= max_angle_diff:
                    evol_dict[initial_id] |= set([final_id])
                    relative_dict[initial_id][final_id] = (distance,
                                                        angle, speed)
        for initial_id in list(relative_dict.keys()):
            if len(relative_dict[initial_id]) == 1:
                relative_dict[initial_id] = list(relative_dict[initial_id].values())[0]
        output_queue.put([index, evol_dict, relative_dict])


    def _fill_dicts_process_limited(self, index, max_distance,
                                max_angle_diff, output_queue,
                                limit=5, amount_of_rods=None):
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
                if final_rod:
                    final_kappa = final_rod.kappa
                    cond1 = initial_kappa-kappa_error > final_kappa
                    cond2 = initial_kappa+kappa_error < final_kappa
                    if cond1 or cond2:
                        continue
                    final_id = final_rod.identifier
                    distance = initial_rod.distance_to_rod(final_rod)
                    angle = initial_rod.angle_between_rods(final_rod)
                    angle = min([angle, 180-angle])
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
        for initial_id in list(relative_dict.keys()):
            if len(relative_dict[initial_id]) == 1:
                relative_dict[initial_id] = list(relative_dict[initial_id].values())[0]
        output_queue.put([index, evol_dict, relative_dict])

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
        output_queue = mp.Queue()
        selected_queue = mp.Queue()
        processes = []
        for index in range(len(self._evolution_dictionaries)):
            processes.append(mp.Process(target=self._leave_only_closer_process,
                                        args=(index, output_queue,
                                              selected_queue, max_distance)))
        running, processes_left = methods.run_processes(processes)
        num_processes = len(running)
        finished = 0
        while finished < num_processes:
            finished += 1
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
                finished -= 1
                new_process = processes_left.pop(0)
                new_process.start()


    def _leave_only_closer_process(self, index, output_queue,
                                selected_queue, max_distance):
        """
        Process.
        """
        selected = set([])
        evol_dict = self._evolution_dictionaries[index]
        relative_dict = self._relative_dictionaries[index]
        for initial_rod_id in list(evol_dict.keys()):
            final_rod_id, distance, angle_diff = self._closer_rod(index,
                                          initial_rod_id, selected,
                                            max_distance)
            evol_dict[initial_rod_id] = final_rod_id
            relative_dict[initial_rod_id] = None
            if distance:
                relative_dict[initial_rod_id] = (distance, angle_diff)
            selected |= set([final_rod_id])
        output_queue.put([index, evol_dict, relative_dict])
        selected_queue.put([index, selected])


    def _closer_rod(self, index, initial_rod_id, selected, max_distance=50):
        """
            If there are multiple choices,
        this erase all but the closest.
        """
        evol_dict = self._evolution_dictionaries[index]
        final_rods = evol_dict[initial_rod_id]
        relative_dict = self._relative_dictionaries[index][initial_rod_id]
        min_distance = max_distance
        final_rod = None
        final_rod_list = list(final_rods)
        if len(final_rod_list) == 1:
            final_rod = final_rod_list[0]
            distance = relative_dict[0]
            angle_diff = relative_dict[1]
        elif len(final_rod_list) == 0:
            return None, None, None
        else:
            for final_rod_id in final_rods:
                relative_values = relative_dict[final_rod_id]
                distance = relative_values[0]
                cond1 = (distance < min_distance)
                cond2 = (final_rod_id not in selected)
                if cond1 and cond2:
                    final_rod = final_rod_id
                    min_distance = distance
            if final_rod:
                angle_diff = relative_dict[final_rod][1]
            else:
                final_rod = None
                min_distance = None
                angle_diff = None
        return final_rod, min_distance, angle_diff


    def compute_dictionaries(self, max_distance=100, max_angle_diff=90,
                            limit=5, amount_of_rods=200):
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
            self._get_vectors()

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
        running, processes_left = methods.run_processes(processes)
        num_processes = len(running)
        finished = 0
        while finished < num_processes:
            finished += 1
            output = output_queue.get()
            index = output[0]
            speeds_vectors = output[1]
            self._speeds_vectors[index] = speeds_vectors
            if len(processes_left):
                finished -= 1
                new_process = processes_left.pop(0)
                new_process.start()

    def _get_vectors_process(self, index, output_queue):
        """
        Process.
        """
        speeds_vectors = {}
        evol_dict = self._evolution_dictionaries[index]
        initial_state = self._states[index]
        final_state = self._states[index+1]
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
        output_queue.put([index, speeds_vectors])


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
        output_queue = mp.Queue()
        processes = []
        for index in range(len(self._evolution_dictionaries)-1):
            process = mp.Process(target=self._join_rods_left_process,
                                 args=(index, output_queue, max_distance))
            processes.append(process)
        running, processes_left = methods.run_processes(processes)
        num_processes = len(running)
        finished = 0
        while finished < num_processes:
            finished += 1
            output = output_queue.get()
            index = output[0]
            evol_dict = output[1]
            relative_dict = output[2]
            self._evolution_dictionaries[index] = evol_dict
            self._relative_dictionaries[index] = relative_dict
            if len(processes_left):
                finished -= 1
                new_process = processes_left.pop(0)
                new_process.start()


    def _join_rods_left_process(self, index, output_queue, max_distance=50):
        """
        Process for method.
        """
        evol_dict = self._evolution_dictionaries[index]
        initial_rods = set([])
        relative_dict = self._relative_dictionaries[index]
        self._initial_rods[index] = initial_rods
        final_rods = self._final_rods[index]
        initial_state = self._states[index]
        final_state = self._states[index+1]
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
                min_distance = -1
                angle_diff = -1
            relative_dict[selected_rod_id] = (min_distance, angle_diff)
            initial_rods -= set([selected_rod_id])
        output_queue.put([index, evol_dict, relative_dict])


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
            speeds_queue = mp.Queue()
            angular_speeds_queue = mp.Queue()
            processes = []
            for index in range(len(self._evolution_dictionaries)-1):
                process = mp.Process(target=self._compute_speeds_process,
                                     args=(index, speeds_queue,
                                            angular_speeds_queue))
                processes.append(process)
            running, processes_left = methods.run_processes(processes)
            num_processes = len(running)
            finished = 0
            while finished < num_processes:
                finished += 1
                speeds = speeds_queue.get()
                angular_speeds = angular_speeds_queue.get()
                self._speeds.append(speeds)
                self._angular_speeds.append(angular_speeds)
                if len(processes_left):
                    finished -= 1
                    new_process = processes_left.pop(0)
                    new_process.start()


    def _compute_speeds_process(self, index, speeds_queue,
                                angular_speeds_queue):
        """
        Returns an array of speeds.
        """
        rel_dict = self._relative_dictionaries[index]
        speeds = {}
        angular_speeds = {}
        for initial_rod_id in list(rel_dict.keys()):
            values = rel_dict[initial_rod_id]
            try:
                speed = float(values[0])/self._diff_t
                angular_speed = float(values[1])/self._diff_t
                speeds[initial_rod_id] = speed
                angular_speeds[initial_rod_id] = angular_speed
            except TypeError:
                pass
        speeds_queue.put(speeds)
        angular_speeds_queue.put(angular_speeds)


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
            num_of_rods = len(self._speeds[index])
            output.append(0)
            for speed in list(self._speeds[index].values()):
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
        """
        output_queue = mp.Queue()
        processes = []
        for index in range(len(self._evolution_dictionaries)-1):
            process = mp.Process(target=self._compute_local_speeds_process,
                                args=(index, output_queue, divisions))
            self._local_speeds.append({})
            processes.append(process)
        running, processes_left = methods.run_processes(processes)
        num_processes = len(running)
        finished = 0
        while finished < num_processes:
            finished += 1
            output = output_queue.get()
            index = output[0]
            speeds_matrix = output[1]
            self._local_speeds[index] = speeds_matrix
            if len(processes_left):
                finished -= 1
                new_process = processes_left.pop(0)
                new_process.start()


    def _compute_local_speeds_process(self, index, output_queue, divisions):
        """
        Process
        """
        state = self._states[index]
        subgroups_matrix = state.subgroups_matrix(divisions)
        speeds_matrix = []
        for row in subgroups_matrix:
            speeds_row = []
            for subsystem in row:
                subsystem_dict = {}
                for rod in subsystem:
                    rod_id = rod.identifier
                    try:
                        speed = self._speeds[index][rod_id]
                        angular_speed = self._angular_speeds[index][rod_id]
                        subsystem_dict[rod_id] = (speed, angular_speed)
                    except KeyError:
                        pass
                speeds_row.append(subsystem_dict)
            speeds_matrix.append(speeds_row)
        output_queue.put([index, speeds_matrix])


    def _compute_local_average_speeds(self, max_distance=100, max_angle_diff=90,
                                     limit=5, amount_of_rods=200, divisions=5):
        """
        Compute local average speeds.
        """
        self.local_speeds(max_distance, max_angle_diff, limit,
                            amount_of_rods, divisions)
        if not len(self._local_average_quadratic_speeds):
            output_queue = mp.Queue()
            processes = []
            local_speeds = self.local_speeds(max_distance, max_angle_diff,
                                            limit, amount_of_rods, divisions)
            for index in range(len(self._evolution_dictionaries)-1):
                local_speeds_ = local_speeds[index]
                process = mp.Process(target=compute_local_average_speeds_process,
                                    args=(index, output_queue, local_speeds_))
                self._local_average_quadratic_speeds.append(None)
                self._local_average_quadratic_angular_speeds.append(None)
                processes.append(process)
            running, processes_left = methods.run_processes(processes)
            num_processes = len(running)
            finished = 0
            while finished < num_processes:
                finished += 1
                output = output_queue.get()
                index = output[0]
                self._local_average_quadratic_speeds[index] = output[1]
                self._local_average_quadratic_angular_speeds[index] = output[2]
                if len(processes_left):
                    finished -= 1
                    new_process = processes_left.pop(0)
                    new_process.start()

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
            for index in range(len(self._states)-1):
                process = mp.Process(target=self._density_and_quad_speed_process,
                            args=(index, output_queue, quad_speeds_array,
                                  ang_speeds_array, divisions))
                processes.append(process)
            running, processes_left = methods.run_processes(processes)
            num_processes = len(running)
            finished = 0
            densities = []
            quad_speeds = []
            while finished < num_processes:
                finished += 1
                output = output_queue.get()
                densities.append(output[0])
                quad_speeds.append(output[1])
                self._speeds_matrices[output[3]](output[2])
                if len(processes_left):
                    finished -= 1
                    new_process = processes_left.pop(0)
                    new_process.start()
            self._densities_array = densities
            self._quad_speeds_array = quad_speeds
        return [self._densities_array, self._quad_speeds_array]


    def _density_and_quad_speed_process(self, index, output_queue,
                                    quad_speeds_array, ang_speeds_array, divisions):
        """
        Process
        """
        quad_speeds = quad_speeds_array[index]
        ang_speeds = ang_speeds_array[index]
        state = self._states[index]
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
        if divisions != self._divisions:
            print("Dividing systems in circles...")
            self._divisions = divisions
            processes = []
            output_queue = mp.Queue()
            for index in range(len(self._states)):
                process = mp.Process(target=self.divide_system_in_circles_process,
                                     args=(divisions, index, output_queue))
                processes.append(process)
            running, processes_left = methods.run_processes(processes, cpus=4)
            num_processes = len(processes)
            finished = 0
            while finished < num_processes:
                finished += 1
                if not len(running):
                    break
                output = output_queue.get()
                index = output[0]
                state = output[1]
                self._states[index] = state
                if len(processes_left):
                    new_process = processes_left.pop(0)
                    new_process.start()
            gc.collect()


    def divide_system_in_circles_process(self, divisions, index, output_queue):
        """
            Process
        """
        state = self._states[index]
        state.divide_in_circles(divisions)
        output_queue.put([index, state])

    def create_density_gif(self, divisions, folder, fps,
                                 number_of_bursts):
        """
        Creates a gif of density's evolution.
        """
        print("Creating densities animations...")
        #self.divide_systems_in_circles(divisions=divisions)
        frames = len(self._states)
        function_name = 'plottable_density_matrix'
        kappas = self._states[0].kappas
        name = str(folder)+str(function_name)+"_K"+str(kappas)+'.mp4'
        z_min, z_max = self._generic_scatter_animator(name, function_name,
                        divisions, fps=fps, number_of_bursts=number_of_bursts)
        self._min_density = z_min
        self._max_density = z_max

    def create_relative_g2_gif(self, divisions, folder, fps,
                                 number_of_bursts):
        """
        Creates a gif of correlation g2 evolution.
        """
        print("Creating g2 animations...")
        #self.divide_systems_in_circles(divisions=divisions)
        frames = len(self._states)
        function_name = 'correlation_g2_plot_matrix'
        #function_name = 'relative_g2_plot_matrix'
        kappas = self._states[0].kappas
        name = str(folder)+str(function_name)+"_K"+str(kappas)+'.mp4'
        self._generic_scatter_animator(name, function_name,
                        divisions, fps=fps, number_of_bursts=number_of_bursts)

    def create_relative_g4_gif(self, divisions, folder, fps,
                                 number_of_bursts):
        """
        Creates a gif of correlation g4 evolution.
        """
        print("Creating g4 animations...")
        #self.divide_systems_in_circles(divisions=divisions)
        frames = len(self._states)
#        function_name = 'relative_g4_plot_matrix'
        function_name = 'correlation_g4_plot_matrix'
        kappas = self._states[0].kappas
        name = str(folder)+str(function_name)+"_K"+str(kappas)+'.mp4'
        self._generic_scatter_animator(name, function_name,
                        divisions, fps=fps, number_of_bursts=number_of_bursts)

    def create_average_angle_gif(self, divisions, folder, fps,
                                 number_of_bursts):
        """
        Creates a gif of average angle evolution.
        """
        print("Creating average angle animations...")
        #self.divide_systems_in_circles(divisions=divisions)
        frames = len(self._states)
        function_name = 'plottable_average_angle_matrix'
        kappas = self._states[0].kappas
        name = str(folder)+str(function_name)+"_K"+str(kappas)+'.mp4'
        self._generic_scatter_animator(name, function_name,
                        divisions, fps=fps, number_of_bursts=number_of_bursts)


    def _generic_scatter_animator(self, name, function_name,
                                    divisions, fps=1, number_of_bursts=1):
        """
        Generic animator
        """
        fig = plt.figure()
        bursts_groups = copy.deepcopy(self.bursts_groups)
        z_vals = []
        z_vals_avg = []
        x_val = []
        y_val = []
        cont = True
        while cont:
            groups = []
            try:
                for burst in range(number_of_bursts):
                    group = bursts_groups.pop(0)
                    groups.append(group)
            except IndexError:
                cont = False
            for group in groups:
                for index in group:
                    state = self._states[index]
                    function = getattr(state, function_name)
                    x_val, y_val, z_val = function(divisions)
                    z_vals.append(z_val)
            z_vals_avg.append(methods.array_average(z_vals))
        frames = len(z_vals_avg)
        match = re.match(r'.*?g[2|4].*', function_name)
        if not match:
            z_maxs = []
            z_mins = []
            for z_val in z_vals_avg:
                z_maxs.append(max(z_val))
                z_mins.append(min(z_val))
            z_max = max(z_maxs)
            z_min = min(z_mins)
        else:
            z_max = 1
            z_min = 0
        def animate(dummy_frame):
            """
            Wrapper.
            """
            self._animate_scatter(x_val, y_val, z_vals_avg,
                                divisions, name, z_max, z_min)
        anim = animation.FuncAnimation(fig, animate, frames=frames)
        anim.save(name, writer=self._writer, fps=fps)
        return z_min, z_max

    def _animate_scatter(self, x_val, y_val, z_vals,
                                divisions, name, z_max, z_min):
        """
        Specific animator.
        """
        try:
            z_val = z_vals.pop(0)
        except IndexError:
            return
        plt.cla()
        plt.clf()
        rad = 2000.0/divisions
        size = (rad/8)**2
        x_min = min(x_val)*.9
        x_max = max(x_val)*1.1
        y_min = min(y_val)*.7
        y_max = max(y_val)*1.1
        plt.xlim((x_min, x_max))
        plt.ylim((y_min, y_max))
        plt.suptitle(name)
        plt.scatter(x_val, y_val, s=size, c=z_val, marker='s',
                    vmin=z_min, vmax=z_max)
        plt.gca().invert_yaxis()
        plt.colorbar()

    def _get_image_id(self, index):
        """
        Returns consecutive system images' ids.
        """
        image1_id_str = self._states[index].id_string
        image1_id = methods.get_number_from_string(image1_id_str)
        return image1_id

    def create_gifs(self, divisions=5, folder="./", fps=1,
                            max_distance=100, max_angle_diff=90, limit=5,
                            amount_of_rods=200, number_of_bursts=1):
        """
        Creates a gif per property of the system that shows evolution.
        """
        self.divide_systems_in_circles(divisions)
        processes = []
        """process = mp.Process(target=self.create_density_gif,
                             args=(divisions, folder, fps, number_of_bursts))
        processes.append(process)"""
        process = mp.Process(target=self.create_relative_g2_gif,
                             args=(divisions, folder, fps, number_of_bursts))
        processes.append(process)
        process = mp.Process(target=self.create_relative_g4_gif,
                             args=(divisions, folder, fps, number_of_bursts))
        processes.append(process)
        """process = mp.Process(target=self.create_temperature_gif,
                             args=(divisions, folder, fps, max_distance,
                               max_angle_diff, limit, amount_of_rods,
                               number_of_bursts))"""
        processes.append(process)
        running, processes_left = methods.run_processes(processes, cpus=1)
        while True:
            try:
                process = running.pop(0)
                process.join()
            except IndexError:
                pass
            try:
                new_process = processes_left.pop(0)
                new_process.start()
                running.append(new_process)
            except IndexError:
                for process in running:
                    process.join()
                break

    def plottable_local_average_quadratic_speeds(self,
                                        max_distance=100,
                                        max_angle_diff=90, limit=5,
                                        amount_of_rods=200, divisions=5):
        """
        Returns plotable data.
        """
        quad_speeds = self.local_average_quadratic_speed(max_distance,
                                        max_angle_diff, limit,
                                        amount_of_rods, divisions)
        #quad_ang_speeds = self.local_average_quadratic_angular_speed(max_distance,
        #                                max_angle_diff, limit,
        #                                amount_of_rods, divisions)
        x_vals, y_vals, z_vals = [], [], []
        for index in range(len(self._states)-1):
            state = self._states[index]
            subgroups = state.subgroups_matrix(divisions)
            x_val, y_val, z_val = [], [], []
            for row_index in range(len(subgroups)):
                for col_index in range(len(subgroups[row_index])):
                    subgroup = subgroups[row_index][col_index]
                    quad_speed = quad_speeds[index][row_index][col_index]
                    #ang_speed = quad_ang_speeds[index][row_index][col_index]
                    center = subgroup.center
                    center_x = center[0]
                    center_y = center[1]
                    #total_speed = quad_speed + ang_speed
                    x_val.append(center_x)
                    y_val.append(center_y)
                    #z_val.append(total_speed)
                    z_val.append(quad_speed)
            x_vals.append(x_val)
            y_vals.append(y_val)
            z_vals.append(z_val)
        return x_vals, y_vals, z_vals

    def create_cluster_histogram_gif(self, max_distance=None,
                                    max_angle_diff=None, fps=15):
        """
            Creates a gif of cluster length histogram.
        """
        kappa = self._states[0].average_kappa
        name = "cluster_hist_K"+str(int(kappa))+".mp4"
        processes = []
        output_queue = mp.Queue()
        arrays = []
        fig = plt.figure()
        for index in range(len(self._states)):
            process = mp.Process(target=self.create_cluster_hist_gif_process,
                                 args=(index, max_distance, max_angle_diff,
                                       output_queue))
            processes.append(process)
            arrays.append(None)
        running, processes_left = methods.run_processes(processes, cpus=4)
        num_processes = len(processes)
        finished = 0
        while finished < num_processes:
            finished += 1
            if not len(running):
                break
            output = output_queue.get()
            index = output[0]
            array = output[1]
            arrays[index] = array
            if len(processes_left):
                new_process = processes_left.pop(0)
                new_process.start()
        def animate(dummy_frame):
            """
            Animation function.
            """
            self._cluster_gif_wrapper(arrays)
        frames = len(self._states)
        anim = animation.FuncAnimation(fig, animate, frames=frames)
        anim.save(name, writer=self._writer, fps=fps)

    def _cluster_gif_wrapper(self, arrays):
        """
        Wrapper.
        """
        try:
            plt.cla()
            plt.clf()
            array = arrays.pop(0)
            boundaries = [0+time*10 for time in range(40)]
            plt.xlim((0, 200))
            plt.ylim((0,0.3))
            plt.hist(array, bins=boundaries, normed=True)
            plt.suptitle("Cluster length histogram")
        except:
            pass

    def create_cluster_hist_gif_process(self, index, max_distance, max_angle_diff,
                                        output_queue):
        """
        Process
        """
        state = self._states[index]
        array = state.cluster_lengths(max_distance=max_distance,
                                            max_angle_diff=max_angle_diff,
                                            min_size=1)
        output_queue.put([index, array])
        return

    def create_speeds_vectors_gif(self, divisions, folder, fps, max_distance,
                                 max_angle_diff, limit, amount_of_rods,
                                 number_of_bursts):
        """
            Creates a gif of average speed vectors over subsystem.
        """
        vectors_matrices = self.average_speeds_vectors(
                                            divisions, max_distance,
                                            max_angle_diff)
        fig = plt.figure()
        kappa = self._states[0].average_kappa
        name = str(folder) + "speeds_vectors_K" + str(kappa)
        bursts_groups = copy.deepcopy(self.bursts_groups)
        end = False
        vectors_matrices_avg = []
        def animate(dummy_frame):
            """
            Animation function.
            """
            self._speeds_vectors_gif_wrapper(vectors_matrix)
        frames = len(self._states)
        anim = animation.FuncAnimation(fig, animate, frames=frames)
        anim.save(name, writer=self._writer, fps=fps)

    def _speeds_vectors_gif_wrapper(self, vectors_matrix):
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
        speeds_num = len(self._states)-1
        for index in range(speeds_num):
            speeds_ = speeds[index]
            state = self._states[index]
            process = mp.Process(target=average_speeds_vectors_gif_process,
                                 args=(divisions, index, state,
                                     speeds_, output_queue))
            processes.append(process)
            vector_matrices.append(None)
        running, processes_left = methods.run_processes(processes, cpus=4)
        num_processes = len(processes)
        finished = 0
        while finished < num_processes:
            finished += 1
            if not len(running):
                break
            output = output_queue.get()
            index = output[0]
            vector_matrix = output[1]
            vector_matrices[index] = vector_matrix
            if len(processes_left):
                new_process = processes_left.pop(0)
                new_process.start()
        vectors = []
        for index in range(speeds_num):
            state = self._states[index]
            vector_matrix = vector_matrices[index]
            

    def create_temperature_gif(self, divisions, folder, fps,
                            max_distance, max_angle_diff,
                            limit, amount_of_rods, number_of_bursts):
        """
        Creates a gif of temperature evolution.
        """
        x_vals, y_vals, z_vals = self.plottable_local_average_quadratic_speeds(
                                        max_distance, max_angle_diff, limit,
                                        amount_of_rods, divisions)
        bursts_groups = copy.deepcopy(self.bursts_groups)
        end = False
        z_vals_avg = []
        number_of_bursts *= 5
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
                    _z_vals.append(z_val)
            try:
                average = methods.array_average(_z_vals)
            except IndexError:
                average = _z_vals
            z_vals_avg.append(average)
        fig = plt.figure()
        kappas = self._states[0].kappas
        name = str(folder)+"Temperature"+str(kappas)+".mp4"
        z_maxs = []
        z_mins = []
        try:
            for z_val in z_vals_avg:
                z_maxs.append(max(z_val))
                z_mins.append(min(z_val))
        except ValueError:
            print(z_vals_avg)
            raise ValueError
        z_max = max(z_maxs)
        z_min = min(z_mins)
        frames = len(z_vals_avg)
        x_val = x_vals[0]
        y_val = y_vals[0]
        def animate(dummy_frame):
            """
            Animation function.
            """
            self._temperature_gif_wrapper(x_val, y_val,
                                        z_vals_avg, divisions, name,
                                        z_max, z_min)
        anim = animation.FuncAnimation(fig, animate, frames=frames)
        anim.save(name, writer=self._writer, fps=fps)


    def _temperature_gif_wrapper(self, x_val, y_val, z_vals,
                                divisions, name, z_max, z_min):
        """
        Wrapper
        """
        try:
            z_val = z_vals.pop(0)
        except IndexError:
            return
        plt.cla()
        plt.clf()
        rad = 2000.0/divisions
        size = (rad/8)**2
        x_min = min(x_val)*.9
        x_max = max(x_val)*1.1
        y_min = min(y_val)*.7
        y_max = max(y_val)*1.1
        plt.xlim((x_min, x_max))
        plt.ylim((y_min, y_max))
        plt.suptitle(name)
        plt.scatter(x_val, y_val, s=size, c=z_val, marker='s',
                    vmax=z_max, vmin=z_min)
        plt.gca().invert_yaxis()
        plt.colorbar()

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
                average = sum(not_none_vals)/float(len(not_none_vals))
            z_vals_avg.append(average)
        return z_vals_avg, indices

    def _get_cluster_areas(self, max_distance=None,
                    max_angle_diff=None, min_size=3):
        """
        Compute cluster areas for all states.
        """
        output_queue = mp.Queue()
        processes = []
        areas = []
        for index in range(len(self._states)):
            state = self._states[index]
            process = mp.Process(target=self._get_cluster_areas_process,
                                 args=(index, state, max_distance,
                                    max_angle_diff, output_queue, min_size))
            processes.append(process)
            areas.append(None)
        running, processes_left = methods.run_processes(processes, cpus=4)
        num_processes = len(processes)
        finished = 0
        while finished < num_processes:
            finished += 1
            if not len(running):
                break
            output = output_queue.get()
            index = output[0]
            area = output[1]
            areas[index] = area
            if len(processes_left):
                new_process = processes_left.pop(0)
                new_process.start()
        return areas
        

    def _get_cluster_areas_process(self, index, state, max_distance,
                    max_angle_diff, output_queue, min_size):
        """
        Process
        """
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
        cluster_areas = self.cluster_areas(number_of_bursts=number_of_bursts,
                                   max_distance=max_distance,
                                   max_angle_diff=max_angle_diff,
                                   min_size=min_size)
        indices = self._indices
        self._compute_times(number_of_bursts)
        times = self._times
        times.pop(0)
        cluster_areas.pop(0)
        log_areas = numpy.array([math.log(area) for area in cluster_areas])
        log_times = numpy.array([math.log(time) for time in times])
        x_0 = numpy.array([0,0,0,0])
        function = lambda value, coef1, coef2: coef1 + coef2*value
        popt, pcov = optimization.curve_fit(function, log_times, log_areas)
        std_dev = numpy.sqrt(numpy.diag(pcov))
        return popt[0], popt[1], std_dev
        

    def plot_cluster_areas(self, number_of_bursts=1, max_distance=None,
                    max_angle_diff=None, min_size=10):
        """
            Plots cluster areas evolution.
        """
        areas = self.cluster_areas(number_of_bursts=number_of_bursts,
                                   max_distance=max_distance,
                                   max_angle_diff=max_angle_diff,
                                   min_size=min_size)
        total_areas = self._total_cluster_areas
        norm_areas = []
        for index in range(len(areas)):
            area = areas[index]
            total_area = total_areas[index]
            proportion = float(area)/total_area
            norm_areas.append(proportion)
        self._compute_times()
        times = self._times
        fig = plt.figure()
        plt.ylim((0,1))
        plt.xlabel("time[seconds]")
        plt.ylabel("cluster area proportion")
        plt.grid()
        try:
            plt.plot(times, norm_areas)
        except ValueError:
            print(len(times), len(areas))
            print(times, areas)
        plt.savefig("cluster_areas.png")

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
        if not len(self._bursts_groups):
            groups = []
            group = []
            for index in range(len(self._state_numbers)-1):
                initial_id  = self._get_image_id(index)
                final_id = self._get_image_id(index+1)
                burst = methods.are_in_burst(self._dates, initial_id,
                                             final_id)
                if burst:
                    group.append(index)
                else:
                    group.append(index)
                    groups.append(group)
                    group = []
            self._bursts_groups = groups
        return self._bursts_groups

    def plot_average_temperature(self, max_distance, max_angle_diff):
        """
            Average temperature over time.
        """
        self._compute_speeds(max_distance, max_angle_diff,
                        5, None)
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
        running, processes_left = methods.run_processes(processes, cpus=4)
        num_processes = len(processes)
        finished = 0
        while finished < num_processes:
            finished += 1
            if not len(running):
                break
            output = output_queue.get()
            index = output[0]
            average_speed = output[1]
            average_speeds[index] = average_speed
            if len(processes_left):
                new_process = processes_left.pop(0)
                new_process.start()
        plt.figure()
        plt.plot(indices, average_speeds)
        name = "avg_temp_K"
        kappa = int(self._states[0].average_kappa)
        name += str(kappa) + ".png"
        plt.savefig(name)

    def average_speeds(self, index, output_queue):
        """
            Averages speeds of all rods.
        """
        speeds = self._speeds[index]
        number_of_rods = len(list(speeds.keys()))
        average_speed = 0
        for speed in list(speeds.values()):
            average_speed += speed
        average_speed /= number_of_rods
        output_queue.put([index, average_speed])

    @property
    def lost_rods_percentage(self):
        """
            Computes maximum percentage of rods lost in analysis.
        """
        number_of_rods = []
        for state in self._states:
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
        

def compute_local_average_speeds_process(index, output_queue, local_speeds):
    """
    Process
    """
    speeds_matrix = []
    angular_speeds_matrix = []
    for row in local_speeds:
        speeds_row = []
        angular_speeds_row = []
        for dictionary in row:
            quadratic_speed = 0
            quadratic_angular_speed = 0
            num_rods = len(dictionary)
            for speeds in list(dictionary.values()):
                quadratic_speed += float(speeds[0]**2)/num_rods
                quadratic_angular_speed += float(speeds[1]**2)/num_rods
            speeds_row.append(quadratic_speed)
            angular_speeds_row.append(quadratic_angular_speed)
        speeds_matrix.append(speeds_row)
        angular_speeds_matrix.append(angular_speeds_row)
    output_queue.put([index, speeds_matrix, angular_speeds_matrix])

def average_speeds_vectors_gif_process(divisions, index, state,
                                     speeds, output_queue):
    """
        Process
    """
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

