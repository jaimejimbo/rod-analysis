"""
    Library for time evolution study.
"""
import re, methods
import multiprocessing as mp
from matplotlib import animation
import matplotlib.pyplot as plt


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
        elif re.match(r'.*list.*', type_):
            self._states = system_states_list
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
        for index in range(len(self._state_numbers)):
            number = self._state_numbers[index]
            state = self._states[index]
            self._states_dict[number] = state
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
        self._speeds = []
        self._angular_speeds = []
        self._max_speed = None
        self._max_angle_diff = None
        self._limit = None
        self._amount_of_rods = None
        self._local_speeds = []
        self._local_average_quadratic_speeds = []
        self._local_average_quadratic_angular_speeds = []
        self._densities_array = []
        self._quad_speeds_array = []
        self._dates = dates


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
        self._speeds = []
        self._angular_speeds = []
        self._max_speed = None
        self._max_angle_diff = None
        self._limit = None
        self._amount_of_rods = None
        self._local_speeds = []
        self._local_average_quadratic_speeds = []
        self._local_average_quadratic_angular_speeds = []
        self._densities_array = []
        self._quad_speeds_array = []


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
            evol_dict = self._evolution_dictionaries[index]
            relative_dict = self._relative_dictionaries[index]
            for rod in state:
                self._initial_rods[index] |= set([rod.identifier])
                rod_id = rod.identifier
                evol_dict[rod_id] = set([])
                relative_dict[rod_id] = {}
        for index in range(len(self._states)-1):
            self._final_rods[index] = self._initial_rods[index+1].copy()


    def _fill_dicts(self, max_speed, max_angle_diff,
                    limit=5, amount_of_rods=None):
        """
            Looks for rods that have only one possible predecessor.
        """
        (self._max_speed, self._max_angle_diff,
        self._limit, self._amount_of_rods) = (max_speed, max_angle_diff,
                                                limit, amount_of_rods)
        processes = []
        output_queue = mp.Queue()
        #for state_id in self._states_dict.keys():
        # Hay que cambiar la forma de acceder en el proceso.
        # o hacer un dict {index: identifier}
        for index in range(len(self._states)-1):
            processes.append(mp.Process(target=self._fill_dicts_process_limited,
                                    args=(index, max_speed, max_angle_diff,
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
                new_process = processes_left.pop()
                new_process.start()


    def _fill_dicts_process(self, index, max_speed, max_angle_diff,
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
                available_final_rods = final_state.get_rods_range(0,
                                        final_state.number_of_rods)
            else:
                start_id = initial_id-amount_of_rods/2
                end_id = initial_id+amount_of_rods/2
                available_final_rods = final_state.get_rods_range(start_id,
                                                                    end_id)
            for final_rod in available_final_rods:
                final_id = final_rod.identifier
                distance = initial_rod.distance_to_rod(final_rod)
                angle = initial_rod.angle_between_rods(final_rod)
                angle = min([angle, 180-angle])
                speed = float(distance)/self._diff_t
                if speed <= max_speed and angle <= max_angle_diff:
                    evol_dict[initial_id] |= set([final_id])
                    relative_dict[initial_id][final_id] = (distance,
                                                        angle, speed)
        for initial_id in relative_dict.keys():
            if len(relative_dict[initial_id]) == 1:
                value = relative_dict[initial_id].values()
                relative_dict[initial_id] = value
        output_queue.put([index, evol_dict, relative_dict])


    def _fill_dicts_process_limited(self, index, max_speed,
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
                available_final_rods = final_state.get_rods_range(start_id,
                                                                  end_id)
            speeds = []
            initial_kappa = initial_rod.kappa
            kappa_error = initial_state.kappa_error/2
            for final_rod in available_final_rods:
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
                if speed <= max_speed and angle <= max_angle_diff:
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
                value = relative_dict[initial_id].values()[0]
                relative_dict[initial_id] = value
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
        for rod_id in evol_dict.keys():
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
        for index in range(len(self._evolution_dictionaries)-1):
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
                new_process = processes_left.pop()
                new_process.start()


    def _leave_only_closer_process(self, index, output_queue,
                                selected_queue, max_distance):
        """
        Process.
        """
        selected = set([])
        evol_dict = self._evolution_dictionaries[index]
        relative_dict = self._relative_dictionaries[index]
        for initial_rod_id in evol_dict.keys():
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


    def compute_dictionaries(self, max_speed=100, max_angle_diff=90,
                            limit=5, amount_of_rods=200):
        """
                List of evolution dictionaries.
        Each dictionary has the form:
            {initial_rod_id1: set([final_rod_id11,final_rod_id12,...]),
             initial_rod_id2: set([final_rod_id21,final_rod_id22,...]),
             ...
             initial_rod_idN: set([final_rod_idN1,final_rod_idN2,...])}
        """
        tuple1 = (max_speed, max_angle_diff, limit, amount_of_rods)
        tuple2 = (self._max_speed, self._max_angle_diff,
                    self._limit, self._amount_of_rods)
        if tuple1 != tuple2:
            self._reset()
        if not len(self._evolution_dictionaries):
            self._create_dict_keys()
            self._fill_dicts(max_speed, max_angle_diff, limit=limit,
                                amount_of_rods=amount_of_rods)
            #self._use_unique_evolutions()
            self._leave_only_closer(max_distance=max_speed)
            self._join_left_rods(max_distance=max_speed)


    def evolution_dictionaries(self, max_speed=100, max_angle_diff=90,
                                limit=5, amount_of_rods=200):
        """
            List of evolution dictionaries.
        Each dictionary has the form:
            {initial_rod_id1: set([final_rod_id11,final_rod_id12,...]),
             initial_rod_id2: set([final_rod_id21,final_rod_id22,...]),
             ...
             initial_rod_idN: set([final_rod_idN1,final_rod_idN2,...])}
        """
        self.compute_dictionaries(max_speed=max_speed,
                                  max_angle_diff=max_angle_diff,
                                  limit=limit,
                                  amount_of_rods=amount_of_rods)
        return self._evolution_dictionaries


    def _join_left_rods(self, max_distance=50):
        """
        After using methods listed before, some rods are unjoined.
        This joins closest rods.
        """
        output_queue = mp.Queue()
        processes = []
        for index in range(len(self._evolution_dictionaries)-1):
            process = mp.Process(target=self._join_left_rods_process,
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
                new_process = processes_left.pop()
                new_process.start()


    def _join_left_rods_process(self, index, output_queue, max_distance=50):
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


    def _compute_speeds(self, max_speed, max_angle_diff, limit, amount_of_rods):
        """
        After using methods listed before, some rods are unjoined.
        This joins closest rods.
        """
        tuple1 = (max_speed, max_angle_diff, limit, amount_of_rods)
        tuple2 = (self._max_speed, self._max_angle_diff,
                    self._limit, self._amount_of_rods)
        if tuple1 != tuple2:
            self._reset()
        self.compute_dictionaries(max_speed, max_angle_diff,
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
                    new_process = processes_left.pop()
                    new_process.start()


    def _compute_speeds_process(self, index, speeds_queue,
                                angular_speeds_queue):
        """
        Returns an array of speeds.
        """
        rel_dict = self._relative_dictionaries[index]
        speeds = {}
        angular_speeds = {}
        for initial_rod_id in rel_dict.keys():
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


    def speeds(self, max_speed=100, max_angle_diff=90, limit=5,
                     amount_of_rods=200):
        """
        Returns [speeds, angular_speeds] where both outputs are array with one
        value for each rod.
        """
        self._compute_speeds(max_speed, max_angle_diff, limit, amount_of_rods)
        return [self._speeds, self._angular_speeds]


    def average_quadratic_speed(self, max_speed=100, max_angle_diff=90,
                                limit=5, amount_of_rods=200):
        """
        Returns average quadratic speeds
        """
        self._compute_speeds(max_speed, max_angle_diff, limit, amount_of_rods)
        output = []
        for index in range(len(self._speeds)):
            num_of_rods = len(self._speeds[index])
            output.append(0)
            for speed in self._speeds[index].values():
                output[index] += speed**2/num_of_rods
        return output


    def average_quadratic_angular_speed(self, max_speed=100, max_angle_diff=90,
                                        limit=5, amount_of_rods=200):
        """
        Returns average quadratic angular speed
        """
        self._compute_speeds(max_speed, max_angle_diff, limit, amount_of_rods)
        output = []
        for index in range(len(self._angular_speeds)):
            num_of_rods = len(self._angular_speeds[index])
            output.append(0)
            for angular_speed in self._angular_speeds[index].values():
                output[index] += angular_speed**2/num_of_rods
        return output


    def local_speeds(self, max_speed=100, max_angle_diff=90, limit=5,
                            amount_of_rods=200, rad=50):
        """
        Returns local_speeds array.
        """
        self._compute_speeds(max_speed, max_angle_diff, limit, amount_of_rods)
        if len(self._local_speeds) == 0:
            self._compute_local_speeds(rad)
        return self._local_speeds


    def _compute_local_speeds(self, rad):
        """
        Creates an array of matrices. Each matrix's entry is a dictionariy such
        as {rod_id: (speed, angular_speed)}
        """
        output_queue = mp.Queue()
        processes = []
        for index in range(len(self._evolution_dictionaries)-1):
            process = mp.Process(target=self._compute_local_speeds_process,
                                args=(index, output_queue, rad))
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
                new_process = processes_left.pop()
                new_process.start()


    def _compute_local_speeds_process(self, index, output_queue, rad):
        """
        Process
        """
        state = self._states[index]
        subgroups_matrix = state.subgroups_matrix(rad)
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


    def _compute_local_average_speeds(self, max_speed=100, max_angle_diff=90,
                                      limit=5, amount_of_rods=200, rad=50):
        """
        Compute local average speeds.
        """
        self.local_speeds(max_speed, max_angle_diff, limit,
                            amount_of_rods, rad)
        if not len(self._local_average_quadratic_speeds):
            output_queue = mp.Queue()
            processes = []
            local_speeds = self.local_speeds(max_speed, max_angle_diff,
                                            limit, amount_of_rods, rad)
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
                local_avg_quad_speeds_matrix = output[1]
                local_avg_quad_angular_speeds_matrix = output[2]
                self._local_average_quadratic_speeds[index] = local_avg_quad_speeds_matrix
                self._local_average_quadratic_angular_speeds[index] = local_avg_quad_angular_speeds_matrix
                if len(processes_left):
                    finished -= 1
                    new_process = processes_left.pop()
                    new_process.start()

    def local_average_quadratic_speed(self, max_speed=100, max_angle_diff=90,
                                        limit=5, amount_of_rods=200, rad=50):
        """
        Returns a array of matrices. Each matrix represents
        local average quadratic speed values.
        """
        self._compute_local_average_speeds(max_speed, max_angle_diff, limit,
                                            amount_of_rods, rad)
        return self._local_average_quadratic_speeds

    def local_average_quadratic_angular_speed(self, max_speed=100,
                                        max_angle_diff=90, limit=5,
                                        amount_of_rods=200, rad=50):
        """
        Returns a array of matrices. Each matrix represents
        local average quadratic angular speed values.
        """
        self._compute_local_average_speeds(max_speed, max_angle_diff,
                                            limit, amount_of_rods, rad)
        return self._local_average_quadratic_angular_speeds

    def density_and_quad_speed(self, max_speed=100, max_angle_diff=90,
                                limit=5, amount_of_rods=200, rad=50):
        """
        Returns 2 arrays: density values and temperature
        """
        if (max_speed, max_angle_diff, limit, amount_of_rods) != (self._max_speed,
                    self._max_angle_diff, self._limit, self._amount_of_rods):
            self._reset()
        if not len(self._densities_array):
            quad_speeds_array = self.local_average_quadratic_speed(max_speed=max_speed, max_angle_diff=max_angle_diff,
                                                          limit=limit, amount_of_rods=amount_of_rods, rad=rad)
            ang_speeds_array = self.local_average_quadratic_angular_speed(max_speed=max_speed, max_angle_diff=max_angle_diff,
                                                          limit=limit, amount_of_rods=amount_of_rods, rad=rad)
            output_queue = mp.Queue()
            processes = []
            for index in range(len(self._states)-1):
                process = mp.Process(target=self._density_and_quad_speed_process,
                            args=(index, output_queue, quad_speeds_array,
                                  ang_speeds_array, rad))
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
                if len(processes_left):
                    finished -= 1
                    new_process = processes_left.pop()
                    new_process.start()
            self._densities_array = densities
            self._quad_speeds_array = quad_speeds
        return [self._densities_array, self._quad_speeds_array]


    def _density_and_quad_speed_process(self, index, output_queue,
                                    quad_speeds_array, ang_speeds_array, rad):
        """
        Process
        """
        quad_speeds = quad_speeds_array[index]
        ang_speeds = ang_speeds_array[index]
        state = self._states[index]
        subgroups = state.subgroups_matrix(rad)
        densities = []
        speeds = []
        for row in range(len(quad_speeds)):
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
                densities.append(density)
                speeds.append(total_quad_speed)
        output_queue.put([densities, speeds])

    def create_density_gif(self, rad=50, folder="./", fps=1):
        """
        Creates a gif of density's evolution.
        """
        frames = len(self._states)
        function_name = 'plottable_density_matrix'
        kappas = self._states[0].kappas
        name = str(folder)+str(function_name)+"_K"+str(kappas)+'.gif'
        self._generic_scatter_animator(frames, name, function_name,
                            rad, fps=fps)

    def create_relative_g2_gif(self, rad=50, folder="./", fps=1):
        """
        Creates a gif of correlation g2 evolution.
        """
        frames = len(self._states)
        function_name = 'relative_g2_plot_matrix'
        kappas = self._states[0].kappas
        name = str(folder)+str(function_name)+"_K"+str(kappas)+'.gif'
        self._generic_scatter_animator(frames, name, function_name,
                            rad, fps=fps)

    def create_relative_g4_gif(self, rad=50, folder="./", fps=1):
        """
        Creates a gif of correlation g4 evolution.
        """
        frames = len(self._states)
        function_name = 'relative_g4_plot_matrix'
        kappas = self._states[0].kappas
        name = str(folder)+str(function_name)+"_K"+str(kappas)+'.gif'
        self._generic_scatter_animator(frames, name, function_name,
                            rad, fps=fps)

    def create_average_angle_gif(self, rad=50, folder="./", fps=1):
        """
        Creates a gif of average angle evolution.
        """
        frames = len(self._states)
        function_name = 'plottable_average_angle_matrix'
        kappas = self._states[0].kappas
        name = str(folder)+str(function_name)+"_K"+str(kappas)+'.gif'
        self._generic_scatter_animator(frames, name, function_name,
                            rad, fps=fps)

    def _generic_scatter_animator(self, frames, name, function_name,
                                    rad, fps=1):
        """
        Generic animator
        """
        fig = plt.figure()
        self._last_index = 0
        def animate(dummy_frame):
            """
            Wrapper.
            """
            if self._last_index != -1:
                self._last_index = self._animate_scatter(function_name,
                                                rad, self._last_index)
        anim = animation.FuncAnimation(fig, animate, frames=frames)
        anim.save(name, writer='imagemagick', fps=fps)

    def _animate_scatter(self, function_name, rad, last_index):
        """
        Specific animator.
        """
        dates = self._dates
        index = last_index
        number_of_states = len(self._states)
        if index == number_of_states-2:
            return -1
        image1_id, image2_id = self._get_image_ids(index)
        x_vals, y_vals, z_vals = [], [], []
        x_val, y_val, z_val = [], [], []
        while True:
            state = self._states[index]
            function = getattr(state, function_name)
            x_val, y_val, z_val = function(rad)
            z_vals.append(z_val)
            index += 1
            if index == number_of_states-2:
                break
            image1_id, image2_id = self._get_image_ids(index)
            if not methods.is_in_burst(dates, image1_id, image2_id):
                index += 1
                break
        if len(z_vals) > 1:
            z_val = methods.array_average(z_vals)
        if len(z_val):
            plt.cla()
            plt.clf()
            size = (0.8*rad/2.0)**2
            plt.scatter(x_val, y_val, c=z_val, s=size, marker='s')
            xmin = min(x_val)-rad
            xmax = max(x_val)+rad
            ymin = min(y_val)-rad
            ymax = max(y_val)+rad
            plt.xlim((xmin, xmax))
            plt.ylim((ymin, ymax))
            try:
                plt.colorbar()
            except TypeError as e:
                print e
                print "Colorbar exception: "
                print len(x_val), len(y_val), len(z_val)
                print "\n\n"
        return index

    def _get_image_ids(self, index):
        """
        Returns consecutive system images' ids.
        """
        image1_id_str = self._states[index].id_string
        image2_id_str = self._states[index+1].id_string
        image1_id = methods.get_number_from_string(image1_id_str)
        image2_id = methods.get_number_from_string(image2_id_str)
        return image1_id, image2_id

    def create_gifs(self, rad=50, folder="./", fps=1):
        """
        Creates a gif per property of the system that shows evolution.
        """
        self.create_density_gif(rad, folder, fps)
        self.create_relative_g2_gif(rad, folder, fps)
        self.create_relative_g4_gif(rad, folder, fps)



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
            for speeds in dictionary.values():
                quadratic_speed += float(speeds[0]**2)/num_rods
                quadratic_angular_speed += float(speeds[1]**2)/num_rods
            speeds_row.append(quadratic_speed)
            angular_speeds_row.append(quadratic_angular_speed)
        speeds_matrix.append(speeds_row)
        angular_speeds_matrix.append(angular_speeds_row)
    output_queue.put([index, speeds_matrix, angular_speeds_matrix])

