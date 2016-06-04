"""
        System State. Group of rods in a determined moment.
"""
import math, queue, matrix, copy
import multiprocessing as mp
import methods, rod, settings
import cPickle, zlib, datetime, time, inspect

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'
if settings.special_chars:
    WHITE_BLOCK = u'\u25A0'
else:
    WHITE_BLOCK = 'X'

class SystemState(object):
    """
        Group of rods in a moment.
    Each image has to be translated into a RodGroup (by this class?)
    """
    def __init__(self, kappas=10, real_kappas=10, allowed_kappa_error=.5,
            radius_correction_ratio=0,
            id_string="", zone_coords=None, rods=None):
        """
            Initialization
        """
        if not zone_coords:
            zone_coords = []
        self._rods = queue.Queue(rods)
        self._is_subsystem = False
        self._kappas = kappas
        self._real_kappas = real_kappas
        self._allowed_kappa_error = allowed_kappa_error
        self._radius_correction_ratio = radius_correction_ratio
        self._id_string = id_string
        self._gaussian_exp = {}
        self._real_radius = 70.0
        try:
            self._zone_coords = list(zone_coords)
            self._fixed_center_radius = True
        except IndexError:
            self._fixed_center_radius = False
        self._coef = None
        self._reset()

    def _reset(self):
        """
            Called when system is changed..
        reset all important values, so they must be
        computed again.
        """
        self._actual_subdivision = []
        self._subdivision_centers = []
        self._density_matrix = []
        self._correlation_Q2 = None
        self._correlation_Q4 = None
        self._correlation_Q2_subsystems = []
        self._correlation_Q4_subsystems = []
        self._relative_Q2_subsystems = []
        self._relative_Q4_subsystems = []
        #self._average_kappa = None
        #self._kappa_dev = None
        self._average_angle = None
        self._angle_matrix = []
        self._density = None
        self._relative_Q2 = None
        self._relative_Q4 = None
        self._closest_rod_matrix = []
        self._direction_matrix = []
        self._rods_dict = {}
        self._clusters = []
        self._cluster_checked_dict = {}
        self._clusters_max_distance = None
        self._clusters_max_angle_diff = None
        self._divisions = None
        self._fixed = True
        self._length, self._length_error = None, None
        if not self._fixed_center_radius:
            self._zone_cords = None

    @property
    def scale(self):
        """
        Gets system scale (real_rad(mm) / rad)
        """
        return self._real_radius*1.0/self._zone_coords[2]

    #def set_kappa(self, value):
    #    """
    #        Changes kappa of all rods in system.
    #    """
    #    for rod_ in self:
    #        rod_.kappa = value

    @property
    def divisions(self):
        """
            Returns number of divisions per axis.
        """
        return self._divisions

    @divisions.setter
    def divisions(self, value):
        """
            Changes divisions value.
        """
        self._divisions = value
        self._reset()

    @property
    def coef(self):
        """
            Returns coefficient for circle division.
        """
        return self._coef

    @coef.setter
    def coef(self, value):
        """
            Changes coef value.
        """
        self._reset()
        self._coef = value

    @property
    def rods_possitions(self):
        """
            Returns possitions of rods.
        """
        x_values = []
        y_values = []
        for rod_ in self:
            x_values.append(rod_.x_mid)
            y_values.append(rod_.y_mid)
        return x_values, y_values

    @property
    def id_string(self):
        """
        Returns the name of rod list.
        """
        return self._id_string

    @property
    def kappas(self):
        """
        Returns kappas of system.
        """
        return self._kappas

    @property
    def kappa_error(self):
        """
        Returns kappa error
        """
        return self._allowed_kappa_error

    def __getitem__(self, rod_id):
        """
        Magic method for [].
        """
        self.fill_dicts()
        try:
            return methods.decompress_rod(self._rods_dict[rod_id])
        except KeyError:
            msg = str(rod_id)
            msg += " " + str(self._rods_dict.keys())
            raise IndexError(msg)

    def get(self, index):
        """
        Returns rod with index.
        """
        return list(self._rods)[index]

    def _get_rods_range(self, initial_id, final_id):
        """
        Returns a list of rods between initial_id and final_id.
        """
        output = []
        for ident in range(initial_id, final_id+1):
            try:
                output.append(self[ident])
            except IndexError:
                pass
        return output

    def __iter__(self):
        """
        Magic method for in.
        """
        for rod_ in self._rods:
            rod__ = methods.decompress_rod(rod_)
            yield rod__

    def __list__(self):
        """
        Returns a list of rods
        """
        output = [methods.decompress_rod(rod_) for rod_ in self]
        return output

    def __len__(self):
        """
        Number of rods
        """
        return len(self._rods)

    @property
    def clone(self):
        """
        Returns a copy of this object.
        """
        _zone_coords = None
        if len(self._zone_coords):
            _zone_coords = self._zone_coords
        clone = SystemState(kappas=self._kappas, real_kappas = self._real_kappas,
                        allowed_kappa_error=self._allowed_kappa_error,
                        radius_correction_ratio=self._radius_correction_ratio,
                        id_string=self.id_string, zone_coords=_zone_coords,
                        rods=self._rods)
        return clone

    @property
    def rods_dictionary(self):
        """
            Returns a dictionary of rods ordered by id.
        """
        return self._rods_dict.copy()

    @property
    def number_of_rods(self):
        """
            Returns the number of rods.
        """
        return len(self)

    def put_rod(self, rod_):
        """
            Adds a rod to the group
        """
        self._rods.put(methods.compress_rod(rod_))

    def _get_rod(self):
        """
            Returns the first rod in the queue
        """
        rod_ = self._rods.get()
        self._rods.put(methods.decompress_rod(rod_))
        return rod_

    def _remove_rod(self, rod_):
        """
            Removes a rod from the group (queue object mod needed)
        """
        self._rods.delete(methods.compress_rod(rod_))

    def fill_dicts(self):
        """
            Fill dictionaries.
        """
        if not len(self._cluster_checked_dict):
            for rod_ in self:
                identifier = rod_.identifier
                self._rods_dict[identifier] = methods.compress_rod(rod_)
                self._cluster_checked_dict[identifier] = False

    def compute_center_and_radius(self):
        """
            Computes where the center of the system is and its
        radius.
        """
        if not self._fixed_center_radius:
            # There must be rods to make statistics.
            if not self.number_of_rods:
                msg = "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "center_and_radius can't be computed before adding rods"
                raise ValueError(msg)
            try:
                not_defined = len(self._zone_coords)
            except AttributeError:
                not_defined = True
            if not_defined:
                x_values = []
                y_values = []
                for rod_ in self:
                    x_values.append(rod_.x_mid)
                    y_values.append(rod_.y_mid)
                # center is the mean position of all particles.
                center_x = sum(x_values)*1.0/self.number_of_rods
                center_y = sum(y_values)*1.0/self.number_of_rods
                # radius is the average of maximum distances /2.
                radius = (max(x_values)-min(x_values)+max(y_values)-min(y_values))/4.0
                # radius *= (1-self._zone_coords[2]_correction_ratio)/4.0
                self._zone_coords[0] = center_x
                self._zone_coords[1] = center_y
                self._zone_coords[2] = radius
        self._zone_coords = (self._zone_coords[0], self._zone_coords[1], self._zone_coords[2])
        return self._zone_coords[0], self._zone_coords[1], self._zone_coords[2]

    @property
    def area(self):
        """
            Area of the subsystem.
        """
        return math.pi*self._real_radius**2

    @property
    def center(self):
        """
            Returns center of the system.
        """
        self.compute_center_and_radius()
        return self._zone_coords[0], self._zone_coords[1]

    @property
    def radius(self):
        """
            Returns radius of the system.
        """
        self.compute_center_and_radius()
        return self._zone_coords[2]

    @property
    def zone_coords(self):
        """
            Returns a tuple with center and radius.
        """
        self.compute_center_and_radius()
        return self._zone_coords

    def check_rods(self):
        """
            Check if rods are correct.
        """
        valid_rods = []
        initial_length = len(self._rods)
        for rod_ in self:
            valid = rod_.is_valid_rod(self._kappas,
                        self._allowed_kappa_error,
                        self.zone_coords)
            if valid:
                valid_rods.append(methods.compress_rod(rod_))
        self._rods = queue.Queue(valid_rods)
        final_length = len(self._rods)
        self._reset()

    def check_rods_with_length(self, length, length_error):
        """
            Check if rods has requiered length and if it is in the circle.
        """
        self._length, self._length_error = length, length_error
        valid_rods = []
        initial_length = len(self._rods)
        for rod_ in self:
            valid = rod_.is_valid_rod_length(length,
                        length_error, self._real_kappas,
                        self.zone_coords)
            if valid:
                valid_rods.append(methods.compress_rod(rod_))
        self._rods = queue.Queue(valid_rods)
        final_length = len(self._rods)
        self._reset()

    def divide_in_circles(self, divisions):
        """
            Creates a grid over system with
            "divisions" rows and columns.
        """
        if divisions != self._divisions:
            self._reset()
            self._divisions = divisions
            # Defining zone and distance between points.
            start_x = self.center[0]-self.radius
            end_x = self.center[0]+self.radius
            start_y = self.center[1]-self.radius
            diff = abs(start_x-end_x)/float(divisions)
            # Getting all possible x and y values.
            possible_x_values = [start_x + (times)*diff
                                 for times in range(divisions)]
            possible_y_values = [start_y + (times)*diff
                                 for times in range(divisions)]
            rad = diff*math.sqrt(2)*self._coef/2.0
            real_rad = self._real_radius*math.sqrt(2)*self._coef*1.0/divisions
            subsystems = self._subsystems(possible_x_values, possible_y_values,
                                          rad, diff, divisions, real_rad)
            self._actual_subdivision = subsystems

    def _separate_rods_by_coords(self, rods_list, possible_x_values,
                                    possible_y_values, rad, diff,
                                    divisions):
        """
            Separates rods by zones. It reduces the amount of steps that
        the programm makes.
        BUG: indices get higher values than allowed.
        """
        x_min = min(possible_x_values)
        y_min = min(possible_y_values)
        x_max = max(possible_x_values)
        y_max = max(possible_y_values)
        div_range = range(divisions)
        output = [[[] for dummy_1 in div_range] for dummy_2 in div_range]
        for rod__ in rods_list:
            rod_ = methods.decompress_rod(rod__)
            index_x = int((rod_.x_mid-x_min)/diff)
            index_y = int((rod_.y_mid-y_min)/diff)
            try:
                output[index_y][index_x].append(rod_)
            except IndexError:
                print index_y, index_x
                print len(output), len(output[index_y])
                raise IndexError
        return output

    def _subsystems(self, possible_x_values, possible_y_values, rad, diff,
                            divisions, real_rad):
        """
            Creates subsystems
        """
        subsystems = []
        y_vals = range(len(possible_y_values))
        x_vals = range(len(possible_x_values))
        centers = []
        for actual_y in possible_y_values:
            for actual_x in possible_x_values:
                centers.append((actual_x, actual_y))
        array2 = [self.center for dummy in range(len(centers))]
        distances = methods.compute_distances(centers, array2)
        array2 = None
        for index in range(len(centers)):
            distance = distances[index]
            center = centers[index]
            time.sleep(settings.waiting_time_process)
            if distance < self.radius:
                subsystem = SubsystemState(center, rad, self.zone_coords,
                                           self._rods, self._kappas, self._real_kappas,
                                           self._allowed_kappa_error, real_rad)
                subsystem.check_rods()
                subsystems.append(methods.compress_state(subsystem))
        return subsystems

    def _compute_density_matrix(self, divisions, normalized=False):
        """
            Computes density matrix of the system.
        """
        self.divide_in_circles(divisions)
        density = []
        subdivision = self._actual_subdivision
        for subsystem_ in subdivision:
            subsystem = methods.decompress_state(subsystem_)
            dens = subsystem.density
            if normalized:
                dens /= self.number_of_rods
            subdensity = [subsystem.center[0], subsystem.center[1]]
            subdensity.append(dens)
            density.append(subdensity)
        subdivision = None
        self._density_matrix = methods.compress_state(density)

    def plottable_density_matrix(self, divisions=50):
        """
            Returns 3 arrays: one for x, another for y and the last for values.
        """
        if len(self._density_matrix) == 0:
            self._compute_density_matrix(divisions)
        x_values = []
        y_values = []
        z_values = []
        density_matrix = methods.decompress_state(self._density_matrix)
        for row in density_matrix:
            x_values.append(row[0])
            y_values.append(row[1])
            z_values.append(row[2])
        return x_values, y_values, z_values

    def plottable_density_matrix_queue(self, divisions, index, output_queue):
        """
        Multiprocessing friendly function.
        """
        x_val, y_val, z_val = self.plottable_density_matrix(divisions=divisions)
        output_queue.put([index, x_val, y_val, z_val])

    def _create_subgroups_matrix(self, divisions):
        """
            Put subsystems in a matrix form.
        """
        if divisions != self._divisions:
            self._reset()
        if not len(self._subdivision_centers):
            self.divide_in_circles(divisions)
            act_sub = self._actual_subdivision
            actual_y = methods.decompress_state(act_sub[0]).center[1]
            row = []
            subgroups_matrix = []
            for index in range(len(act_sub)):
                element = methods.decompress_state(act_sub[index])
                element_y = element.center[1]
                if element_y != actual_y:
                    subgroups_matrix.append(row)
                    row = []
                    actual_y = element_y
                row.append(element)
            subgroups_matrix.append(row)
            act_sub = None
            self._subdivision_centers = methods.compress_state(subgroups_matrix)

    def subgroups_matrix(self, divisions):
        """
            Returns subgroups matrix
        """
        self._create_subgroups_matrix(divisions)
        return methods.decompress_state(self._subdivision_centers)

    def _compute_Q2_and_Q4(self):
        """
            Computes correlation_Q2 and correlation_Q4 values
        """
        num = self.number_of_rods
        if num in [0,1,2] or not self.area:
            self._correlation_Q2 = [-1000, -1000] #None
            self._correlation_Q4 = [-1000, -1000] #None
            return
        Q2_cos, Q2_sin = 0, 0
        Q4_cos, Q4_sin = 0, 0
        N = len(self)
        for rod_ in self:
            Q2_cos += float(math.cos(2*rod_.angle))/N
            Q2_sin += float(math.sin(2*rod_.angle))/N
            Q4_cos += float(math.cos(4*rod_.angle))/N
            Q4_sin += float(math.sin(4*rod_.angle))/N
        self._correlation_Q2 = [Q2_cos, Q2_sin]
        self._correlation_Q4 = [Q4_cos, Q4_sin]

    @property
    def correlation_Q2(self):
        """
            sqrt(<cos(2*angle)>^2+<sin(2*angle)>^2)
        """
        if not self._correlation_Q2:
            self._compute_Q2_and_Q4()
        return self._correlation_Q2

    @property
    def correlation_Q4(self):
        """
            sqrt(<cos(4*angle)>^2+<sin(4*angle)>^2)
        """
        if not self._correlation_Q4:
            self._compute_Q2_and_Q4()
        return self._correlation_Q4

    def _compute_Q2_Q4_matrices(self, divisions):
        """
            Computes correlation_Q2 and correlation_Q4 matrices for subgroups.
        """
        self.divide_in_circles(divisions)
        if self._correlation_Q2 is None or self._correlation_Q4 is None:
            subdivision = self._actual_subdivision
            for subsystem_ in subdivision:
                subsystem = methods.decompress_state(subsystem_)
                correlation_Q2 = [subsystem.center[0], subsystem.center[1]]
                correlation_Q4 = [subsystem.center[0], subsystem.center[1]]
                correlation_Q2.append(subsystem.correlation_Q2)
                correlation_Q4.append(subsystem.correlation_Q4)
                self._correlation_Q2_subsystems.append(methods.compress_state(correlation_Q2))
                self._correlation_Q4_subsystems.append(methods.compress_state(correlation_Q4))
            self._correlation_Q2_subsystems = self._correlation_Q2_subsystems
            self._correlation_Q4_subsystems = self._correlation_Q4_subsystems

    def correlation_Q2_plot_matrix(self, divisions):
        """
            Returns values for plotting correlation_Q2 matrix.
        """
        self._compute_Q2_Q4_matrices(divisions)
        x_values = []
        y_values = []
        z_values = []
        Q2_subsystems = self._correlation_Q2_subsystems
        for subsystem_ in Q2_subsystems:
            subsystem = methods.decompress_state(subsystem_)
            [x_val, y_val, z_val] = subsystem
            x_values.append(x_val)
            y_values.append(y_val)
            z_values.append(z_val)
        Q2_subsystems = None
        return x_values, y_values, z_values

    def correlation_Q2_plot_matrix_queue(self, divisions, index, output_queue):
        """
            Multiprocessing friendly function.
        """
        x_val, y_val, z_val = self.correlation_Q2_plot_matrix(divisions)
        output_queue.put([index, x_val, y_val, z_val])

    def correlation_Q4_plot_matrix(self, divisions):
        """
            Returns values for plotting correlation_Q2 matrix.
        """
        self._compute_Q2_Q4_matrices(divisions)
        x_values = []
        y_values = []
        z_values = []
        Q4_subsystems = self._correlation_Q4_subsystems
        for subsystem_ in Q4_subsystems:
            subsystem = methods.decompress_state(subsystem_)
            [x_val, y_val, z_val] = subsystem
            x_values.append(x_val)
            y_values.append(y_val)
            z_values.append(z_val)
        Q4_subsystems = None
        return x_values, y_values, z_values

    def correlation_Q4_plot_matrix_queue(self, divisions, index, output_queue):
        """
            Multiprocessing friendly function.
        """
        x_val, y_val, z_val = self.correlation_Q4_plot_matrix(divisions)
        output_queue.put([index, x_val, y_val, z_val])

    @property
    def average_kappa(self):
        """
            Returns kappa average of group.
        """
        #if not self._average_kappa:
        #    self._average_kappa = 0
        #    for rod_ in self:
        #        self._average_kappa += rod_.kappa
        #    self._average_kappa /= self.number_of_rods
        return self.kappas #self._average_kappa

    @property
    def kappa_dev(self):
        """
            Returns sqrt(<kappa^2> - <kappa>^2)
        """
        #if not self._kappa_dev:
        #    kappa2 = 0
        #    for rod_ in self:
        #        kappa2 += (rod_.kappa-self.average_kappa)**2
        #    kappa2 /= (self.number_of_rods-1)
        #    self._kappa_dev = math.sqrt(kappa2)
        return 0 #self._kappa_dev

    @property
    def average_angle(self):
        """
            Returns average angle of the system (if exists).
        """
        if not self._average_angle:
            if self.correlation_Q2 > 0.5 and self.correlation_Q4 < 0.3:
                angle = 0
                for rod_ in self:
                    angle2 = rod_.angle
                    if angle2 > math.pi:
                        angle2 -= math.pi
                    angle += angle2
                angle /= self.number_of_rods
                self._average_angle = angle
                return angle
            else:
                return None
        else:
            return self._average_angle

    def _compute_average_angle_matrix(self, divisions):
        """
            Computes average angle matrix
        """
        subdivision = self._actual_subdivision
        for subsystem_ in subdivision:
            subsystem = methods.decompress_state(subsystem_)
            row = [subsystem.center[0], subsystem.center[1]]
            row.append(subsystem.average_angle)
            self._angle_matrix.append(row)
        self._angle_matrix = methods.compress_state(self._angle_matrix)

    def plottable_average_angle_matrix(self, divisions):
        """
            Returns a plottable average angle matrix.
        """
        self._compute_average_angle_matrix(divisions)
        x_values = []
        y_values = []
        z_values = []
        angle_matrix = methods.decompress_state(self._angle_matrix)
        for subsystem in angle_matrix:
            x_values.append(subsystem[0])
            y_values.append(subsystem[1])
            z_values.append(subsystem[2])
        angle_matrix = None
        return x_values, y_values, z_values

    def plottable_average_angle_matrix_queue(self, divisions, index, output_queue):
        """
            Multiprocessing friendly function.
        """
        x_val, y_val, z_val = self.plottable_average_angle_matrix(divisions)
        output_queue.put([index, x_val, y_val, z_val])

    @property
    def angle_histogram(self):
        """
            Returns all angles in a list to make an histogram.
        """
        output = [rod_.angle for rod_ in self]
        return output

    def _get_closest_rod(self, rod_):
        """
            Returns closest rod in group to given rod.
        """
        distance = 1e100
        selected_rod = rod_
        for rod2 in self:
            if rod_ != rod2:
                new_distance = rod_.distance_to_rod(rod2)
                if new_distance < distance:
                    distance = new_distance
                    selected_rod = rod2
        return selected_rod

    def clusters(self, max_distance=None, max_angle_diff=None, min_size=3):
        """
            Returns clusters with min_size number of rods or more.
        """
        clusters = self._get_clusters(max_distance=max_distance,
                                  max_angle_diff=max_angle_diff)
        output = []
        for cluster in clusters:
            if len(cluster) >= min_size:
                output.append(cluster)
        return output

    def _get_clusters(self, max_distance=None, max_angle_diff=None):
        """
            Gets the cluster for rod.
        Recursive method.
        Angles in grad.
        """
        self.fill_dicts()
        if len(self._clusters) and not max_distance and not max_angle_diff:
            return self._clusters
        if not max_distance or not max_angle_diff:
            msg = "cluster method without args only valid "
            msg += "when previously computed."
            raise ValueError(msg)
        cond2 = (self._clusters_max_distance != max_distance)
        cond2 &= (self._clusters_max_angle_diff != max_angle_diff)
        if not len(self._clusters) or cond2:
            self._clusters_max_distance = max_distance
            self._clusters_max_angle_diff = max_angle_diff
            clusters = []
            list_of_rods = [rod_ for rod_ in self]
            rods_left = set(list_of_rods)
            for rod_ in self:
                if self._cluster_checked_dict[rod_.identifier]:
                    continue
                rods_left -= set([rod_])
                cluster = self._get_cluster_members(rod_,
                                max_distance, max_angle_diff)
                if len(cluster):
                    clusters.append(cluster)
            self._clusters = clusters
        return clusters

    def _get_cluster_members(self, reference_rod,
                                    max_distance, max_angle_diff):
        """
            Gets the closest neighbour to a rod that fulfill
        some conditions.
        Angles in grad.
        """
        rods = set([reference_rod.identifier])
        if self._cluster_checked_dict[reference_rod.identifier]:
            return set([])
        self._cluster_checked_dict[reference_rod.identifier] = True
        length = reference_rod.feret
        for rod_ in self:
            if self._cluster_checked_dict[rod_.identifier]:
                continue
            vector = reference_rod.vector_to_rod(rod_)
            distance = methods.vector_module(vector)
            angle_diff = rod_.angle_between_rods(reference_rod)
            vector_angle = methods.vector_angle(vector)
            distance_angle = vector_angle-reference_rod.angle
            diff = length-max_distance
            max_dist = max_distance+math.sin(distance_angle)*diff
            if angle_diff <= max_angle_diff and distance < max_dist:
                subrods = self._get_cluster_members(rod_,
                                               max_distance, max_angle_diff)
                rods |= subrods
        return rods

    def average_cluster_rod_num(self, max_distance=None,
                                max_angle_diff=None, min_size=3):
        """
            Gets the average number of rods in clusters.
        Angles in grad.
        """
        lengths = self.cluster_lengths(max_distance,
                                        max_angle_diff, min_size)
        try:
            return float(sum(lengths))/len(lengths)
        except ZeroDivisionError:
            print "No clusters detected."   

    def cluster_lengths(self, max_distance=None,
                            max_angle_diff=None, min_size=3):
        """
            Creates an array of lengths of clusters:
        [1,4,4,5,2] means that there are 1 cluster of lenght 1,
        2 of length 4...
        """
        lengths = []
        clusters = self.clusters(max_distance,
                                max_angle_diff, min_size)
        for cluster in clusters:
            lengths.append(len(cluster))
        return lengths

    def number_of_clusters(self, max_distance=None,
                            max_angle_diff=None, min_size=3):
        """
            Returns the number of clusters in the system.
        Angles in grad.
        """
        clusters = self.clusters(max_distance,
                                max_angle_diff, min_size)
        return len(clusters)

    def total_area_of_clusters(self, max_distance=None,
                                max_angle_diff=None, min_size=3):
        """
            Returns the area covered by clusters.
        """
        rod_area = self.rod_area
        rods_num = self.cluster_lengths(max_distance=max_distance,
                                               max_angle_diff=max_angle_diff,
                                               min_size=min_size)
        rods_num = sum(rods_num)
        if not rods_num or not rod_area:
            return 0
        total_area = rods_num*self.rod_area
        assert total_area, "Total area faillure: "+str(total_area)
        return total_area

    @property
    def rod_area(self):
        """
        Returns the area of a rod of this system.
        """
        first_rod = self.get(0)
        return first_rod.area

    def _compute_closest_rod_matrix(self):
        """
            Creates closer rod matrix:
        [[rod1, closest_rod_to_rod1],
        [rod2, closest_rod_to_rod2],
        ...
        [rodN, closest_rod_to_rodN]]
        """
        function = self._get_closest_rod
        closest_rod_matrix = []
        for rod_ in self:
            new_row = [rod_]
            closest_rod = function(rod_)
            if not closest_rod:
                continue
            new_row.append(closest_rod)
            closest_rod_matrix.append(new_row)
        self._closest_rod_matrix = methods.compress_state(closest_rod_matrix)
        closest_rod_matrix = None

    @property
    def closest_rod_matrix(self):
        """
            Returns closest rod matrix.
        See _compute_closest_rod_matrix for more info.
        """
        if len(self._closest_rod_matrix) == 0:
            self._compute_closest_rod_matrix()
        return methods.decompress_state(self._closest_rod_matrix)

    def closest_rod_dict(self):
        """
            Returns a dictionary of closest rods.
        """
        closest_rod_matrix = self.closest_rod_matrix
        dictionary = {}
        for pair in closest_rod_matrix:
            dictionary[pair[0].identifier] = pair[1]
        closest_rod_matrix = None
        return dictionary

    @property
    def average_angle_using_matrix(self):
        """
            Returns average angle of the system using direction matrices.
        Value in radians.
        """
        if len(self._direction_matrix) == 0:
            self._direction_matrix = matrix.zeros(2, 2)
            for rod_ in self:
                self._direction_matrix += rod_.direction_matrix
        eigen1, dummy_ = self._direction_matrix.diagonalize_2x2()
        return eigen1

    def covered_area_proportion(self):
        """
            Returns covered area proportion by rods.
        """
        total_area = self.area
        covered_area = 0
        for rod in self:
            covered_area += rod.area
        proportion = covered_area / total_area
        return proportion

    def covered_area_proportion_queue(self, index, output_queue):
        """
        Wrapper
        """
        output_queue.put([index, self.covered_area_proportion()])

    @property
    def average_rod_length(self):
        """
            It returns average rod length.
        """
        lengths = []
        for rod_ in self:
            lengths.append(rod_.feret)
        return float(sum(lengths))/len(lengths)

    @property
    def average_rod_width(self):
        """
            It returns average rod width.
        """
        widths = []
        for rod_ in self:
            widths.append(rod_.min_feret)
        return float(sum(widths))/len(widths)

    @property
    def rod_centers(self):
        """
            It returns centers of rods.
        """
        x_values = []
        y_values = []
        for rod_ in self:
            x_val, y_val = rod_.center
            x_values.append(x_val)
            y_values.append(y_val)
        return x_values, y_values

    def set_coords(self, zone_coords):
        """
        Changes zone_coords
        """
        self._zone_coords = zone_coords

    @property
    def plottable_rods(self):
        """
        Returns a list of plotable rods
        """
        x_0_list = []
        y_0_list = []
        x_f_list = []
        y_f_list = []
        scale = self.scale
        for rod_ in self:
            x_mid = rod_.x_mid
            y_mid = rod_.y_mid
            angle = rod_.angle
            length = float(rod_.feret)/scale
            x_0_list.append(x_mid - length*math.cos(angle)/2.0)
            x_f_list.append(x_mid + length*math.cos(angle)/2.0)
            y_0_list.append(y_mid - length*math.sin(angle)/2.0)
            y_f_list.append(y_mid + length*math.sin(angle)/2.0)
        return methods.compress([x_0_list, y_0_list, x_f_list, y_f_list])
        





class SubsystemState(SystemState):
    """
        Group of rods. Used to put all rods that are in a zone or
    have something in common.
    """

    def __init__(self, center, rad, zone_coords, rods, kappas, real_kappas, allowed_kappa_error, real_rad):
        """
            Initialization
        """
        self._subsystem_coords = (center[0], center[1], rad)
        SystemState.__init__(self, rods=rods, zone_coords=self._subsystem_coords,
                            kappas=kappas, real_kappas=real_kappas,
                            allowed_kappa_error=allowed_kappa_error)
        self._is_subsystem = True
        self._center = center
        self._zone_coords[2] = rad
        self._real_radius = real_rad# + real_kappas
        self._scale = real_rad*1.0/rad
        self._main_center = (zone_coords[0], zone_coords[1])
        self._main_rad = zone_coords[2]
        self._position_rad = methods.distance_between_points(self._main_center,
                                            self._center)
        self._area = methods.effective_area(self._zone_coords[2],
                                    self._position_rad, self._main_rad)
        self._real_area = methods.effective_area(self._real_radius, 
                                                 self._position_rad*self._scale,
                                                 self._main_rad*self._scale)
        self._gaussian_exp = {}

    def remove_all_rods(self):
        """
        Removes all rods of subsystem.
        """
        self._rods = []

    @property
    def center(self):
        """
            Center of the subsystem.
        """
        if not self._zone_coords[2]:
            self.compute_center_and_radius()
        return self._center

    @property
    def radius(self):
        """
            Radius of the subsystem
        """
        if not self._zone_coords[2]:
            self.compute_center_and_radius()
        return self._zone_coords[2]

    @property
    def real_area(self):
        """
        Area in mm of subsystem.
        """
        return self._real_area

    def _update_density(self):
        """
            Computes density of the group.
        """
        density = 0
        for rod_ in self:
            density += rod_.area*self._gaussian_exp[rod_.identifier]
        if not density or not self.area:
            self._density = 0
        else:
            self._density = float(density)/self.real_area

    @property
    def density(self):
        """
            Returns the density of the group.
        """
        if not self._density:
            self._update_density()
        return self._density

    def check_rods(self):
        """
            Check if rods are correct.
        """
        rods = []
        for rod_ in self:
            if rod_.is_in_circle(self.center, self.radius):
                distance = methods.distance_between_points(self.center, rod_.center)
                proportion = methods.norm_gaussian(distance, self.radius)
                self._gaussian_exp[rod_.identifier] = proportion
                rods.append(methods.compress_rod(rod_))
        self._reset()
        self._rods = queue.Queue(rods)

    def compute_center_and_radius(self):
        """
        Overrides
        """
        pass


def create_rods(folder="./", kappas=10, real_kappas=10, allowed_kappa_error=.3,
                radius_correction_ratio=0.1, file_range=[]):
    """
    Create one rod for each rod_data and for each file
    returns [RodGroup1, RodGroup2, ...]
    """
    print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Importing data"
    names = methods.get_file_names(folder=folder)
    if len(file_range) != 0:
        try:
            names = names[file_range[0]:file_range[1]]
        except:
            pass
    num_of_files = len(names)
    if not num_of_files:
        print "No files to import."
        raise ValueError
    states = [None for dummy_ in range(num_of_files)]
    processes = []
    states_queue = mp.Queue()
    for index in range(num_of_files):
        process = mp.Process(target=create_rods_process,
                            args=(kappas, real_kappas, allowed_kappa_error,
                            radius_correction_ratio, names,
                            index, states_queue))
        processes.append(process)
    num_processes = len(processes)
    running, processes_left = methods.run_processes(processes)
    finished_ = 0
    previous_time = datetime.datetime.now()
    counter = 0
    time_left = None
    times = []
    print " "
    empty_states = []
    while finished_ < num_processes:
        counter += 1
        finished_ += 1
        previous_time, counter, time_left, times = methods.print_progress(finished_, num_processes,
                                counter, times, time_left, previous_time)
        [index, state] = states_queue.get()
        if state is not None:
            states[index] = state
        else:
            empty_states.append(index)
        if len(processes_left):
            new_process = processes_left.pop(0)
            time.sleep(settings.waiting_time)
            new_process.start()
    for process in processes:
        if process.is_alive():
            process.terminate()
    for empty_state in empty_states:
        diff = 1
        while True:
            new_state = empty_state - diff
            if states[new_state] is not None:
                break
            diff += 1
        new_state = methods.decompress_state(states[new_state]).clone
        states[empty_state] = methods.compress_state(new_state)
    print(CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE)
    return names, states

def create_rods_process(kappas, real_kappas, allowed_kappa_error,
                        radius_correction_ratio, names,
                        index, states_queue):
    """
    Process of method.
    """
    name = names[index]
    file_ = open(name, 'r')

    state = SystemState(kappas=kappas, real_kappas=real_kappas,
                        allowed_kappa_error=allowed_kappa_error,
                        radius_correction_ratio=radius_correction_ratio,
                        id_string=name, zone_coords=settings.zone_coords)
    data = methods.import_data(file_)
    for dataline in data:
        try:
            parameters = tuple(dataline)
            new_rod = rod.Rod(parameters, kappa=real_kappas, real_length=real_kappas)
            state.put_rod(new_rod)
        except ValueError:
            pass
    file_.close()
    file_ = None
    if not state:
        states_queue.put([index, None])
        return
    state.compute_center_and_radius()
    state.check_rods()
    state = methods.compress_state(state)
    states_queue.put([index, state])

def create_rods_with_length(folder="./", length=10, length_error=.3, real_kappas=10,
                radius_correction_ratio=0.1, file_range=[]):
    """
    Create rods using rod length instead of kappa.
    """
    print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Importing data"
    names = methods.get_file_names(folder=folder)
    if len(file_range) != 0:
        try:
            names = names[file_range[0]:file_range[1]]
        except:
            pass
    num_of_files = len(names)
    if not num_of_files:
        print "No files to import."
        raise ValueError
    states = [None for dummy_ in range(num_of_files)]
    processes = []
    states_queue = mp.Queue()
    for index in range(num_of_files):
        process = mp.Process(target=create_rods_with_length_process,
                            args=(length, length_error, real_kappas,
                            radius_correction_ratio, names,
                            index, states_queue))
        processes.append(process)
    num_processes = len(processes)
    running, processes_left = methods.run_processes(processes)
    finished_ = 0
    previous_time = datetime.datetime.now()
    counter = 0
    time_left = None
    times = []
    print " "
    empty_states = []
    while finished_ < num_processes:
        counter += 1
        finished_ += 1
        previous_time, counter, time_left, times = methods.print_progress(finished_, num_processes,
                                counter, times, time_left, previous_time)
        [index, state] = states_queue.get()
        if state is not None:
            states[index] = state
        else:
            empty_states.append(index)
        if len(processes_left):
            new_process = processes_left.pop(0)
            time.sleep(settings.waiting_time)
            new_process.start()
    for process in processes:
        if process.is_alive():
            process.terminate()
    for empty_state in empty_states:
        diff = 1
        while True:
            new_state = empty_state - diff
            if states[new_state] is not None:
                break
            diff += 1
        new_state = methods.decompress_state(states[new_state]).clone
        states[empty_state] = methods.compress_state(new_state)
    print(CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE)
    return names, states

def create_rods_with_length_process(length, length_error, real_kappa,
                            radius_correction_ratio, names,
                            index, states_queue):
    """
    Process of method.
    """
    name = names[index]
    file_ = open(name, 'r')
    state = SystemState(kappas=real_kappa, real_kappas=real_kappa,
                        allowed_kappa_error=0,
                        radius_correction_ratio=radius_correction_ratio,
                        id_string=name, zone_coords=settings.zone_coords)
    data = methods.import_data(file_)
    if not len(data):
        states_queue.put([index, None])
        return
    for dataline in data:
        try:
            parameters = tuple(dataline)
            new_rod = rod.Rod(parameters, kappa=real_kappa, real_length=real_kappa)
            state.put_rod(new_rod)
        except ValueError:
            pass
    file_.close()
    file_ = None
    state.compute_center_and_radius()
    state.check_rods_with_length(length, length_error)
    state = methods.compress_state(state)
    states_queue.put([index, state])

