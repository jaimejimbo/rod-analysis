"""
        System State. Group of rods in a determined moment.
"""
import math, queue, matrix, copy
import multiprocessing as mp
import methods, rod


class SystemState(object):
    """
        Group of rods in a moment.
    Each image has to be translated into a RodGroup (by this class?)
    """
    def __init__(self, kappas=10, allowed_kappa_error=.5,
            radius_correction_ratio=0,
            id_string="", zone_coords=None, rods=None):
        """
            Initialization
        """
        if not zone_coords:
            zone_coords = []
        self._rods = queue.Queue(rods)
        self._is_subsystem = False
        self._number_of_rods = len(self._rods)
        self._kappas = kappas
        self._allowed_kappa_error = allowed_kappa_error
        self._radius_correction_ratio = radius_correction_ratio
        self._id_string = id_string
        self._actual_subdivision = []
        self._subdivision_centers = []
        self._density_matrix = []
        self._correlation_g2 = None
        self._correlation_g4 = None
        self._correlation_g2_subsystems = []
        self._correlation_g4_subsystems = []
        self._relative_g2_subsystems = []
        self._relative_g4_subsystems = []
        self._average_kappa = None
        self._kappa_dev = None
        self._average_angle = None
        self._angle_matrix = []
        self._density = None
        self._relative_g2 = None
        self._relative_g4 = None
        self._closest_rod_matrix = []
        self._direction_matrix = []
        self._rods_dict = {}
        self._clusters = []
        self._cluster_checked_dict = {}
        self._clusters_max_distance = None
        self._clusters_max_angle_diff = None
        self._divisions = None
        try:
            self._radius = zone_coords[2]
            self._center_x = zone_coords[0]
            self._center_y = zone_coords[1]
            self._zone_coords = zone_coords
            self._fixed_center_radius = True
        except IndexError:
            self._fixed_center_radius = False
            self._zone_cords = zone_coords
            self._center_x = None
            self._center_y = None
            self._radius = None

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
            return self._rods_dict[rod_id]
        except KeyError:
            msg = str(rod_id)
            msg += " " + str(self._rods_dict.keys())
            raise IndexError(msg)

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
            yield rod_

    def __list__(self):
        """
        Returns a list of rods
        """
        output = [rod_ for rod_ in self]
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
        if not len(self._zone_coords):
            _zone_coords = None
        clone = SystemState(kappas=self._kappas,
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
        self._number_of_rods = len(self._rods)
        return self._number_of_rods

    def _put_rod(self, rod_):
        """
            Adds a rod to the group
        """
        self._rods.put(rod_)

    def _get_rod(self):
        """
            Returns the first rod in the queue
        """
        rod_ = self._rods.get()
        self._rods.put(rod_)
        return rod_

    def _remove_rod(self, rod_):
        """
            Removes a rod from the group (queue object mod needed)
        """
        self._rods.delete(rod_)

    def _reset(self):
        """
            Called when system is changed..
        reset all important values, so they must be
        computed again.
        """
        self._actual_subdivision = []
        self._density_matrix = []
        self._correlation_g2 = None
        self._correlation_g4 = None
        self._correlation_g2_subsystems = []
        self._correlation_g4_subsystems = []
        self._average_kappa = None
        self._kappa_dev = None
        self._average_angle = None
        self._angle_matrix = []
        self._density = None
        self._relative_g2 = None
        self._subdivision_centers = []
        self._relative_g4 = None
        self._closest_rod_matrix = []
        self._relative_g2_subsystems = []
        self._relative_g4_subsystems = []
        self._direction_matrix = []
        self._clusters = []
        self._cluster_checked_dict = {}
        self._rods_dict = {}
        self._clusters_max_distance = None
        self._clusters_max_angle_diff = None
        self._divisions = None
        self._subdivision_centers = []
        if not self._fixed_center_radius:
            self._radius = None
            self._center_x = None
            self._center_y = None
            self._zone_coords = []

    def fill_dicts(self):
        """
            Fill dictionaries.
        """
        if not len(self._cluster_checked_dict):
            for rod_ in self._rods:
                identifier = rod_.identifier
                self._rods_dict[identifier] = rod_
                self._cluster_checked_dict[identifier] = False

    def _compute_center_and_radius(self):
        """
            Computes where the center of the system is and its
        radius.
        """
        #There must be rods to make statistics.
        if self.number_of_rods == 0:
            msg = "center_and_radius can't be computed before adding rods"
            raise ValueError(msg)
        try:
            not_defined = not len(self._zone_coords)
            not_defined = not_defined and not self._fixed_center_radius
        except AttributeError:
            not_defined = True
        if not_defined and not self._is_subsystem:
            x_values = []
            y_values = []
            for rod_ in list(self._rods):
                x_values.append(rod_.x_mid)
                y_values.append(rod_.y_mid)
            #center is the mean position of all particles.
            center_x = sum(x_values)*1.0/self.number_of_rods
            center_y = sum(y_values)*1.0/self.number_of_rods
            #radius is the average of maximum distances /2.
            radius = (max(x_values)-min(x_values)+max(y_values)-min(y_values))
            radius *= (1-self._radius_correction_ratio)/4.0
            self._center_x = center_x
            self._center_y = center_y
            self._radius = radius
            self._zone_coords = (center_x, center_y, radius)
        return self._center_x, self._center_y, self._radius

    @property
    def center(self):
        """
            Returns center of the system.
        """
        self._compute_center_and_radius()
        return self._center_x, self._center_y

    @property
    def radius(self):
        """
            Returns radius of the system.
        """
        self._compute_center_and_radius()
        return self._radius

    @property
    def zone_coords(self):
        """
            Returns a tuple with center and radius.
        """
        self._compute_center_and_radius()
        return self._zone_coords


    def check_rods(self):
        """
            Check if rods are correct.
        """
        zone_coords = [coord for coord in self.zone_coords]
        for rod_ in self._rods:
            valid = rod_.is_valid_rod(self._kappas,
                        self._allowed_kappa_error,
                        zone_coords)
            if not valid:
                self._remove_rod(rod_)
        self._reset()

    def divide_in_circles(self, divisions):
        """
        Wrapper
        """
        if self._divisions != divisions:
            self._reset()
        self._divide_in_circles(divisions)


    def _divide_in_circles(self, divisions):
        """
            Creates a grid over system with
            "divisions" rows and columns.
        """
        if divisions != self._divisions:
            self._compute_center_and_radius()
            self._divisions = divisions
            # Defining zone and distance between points.
            start_x = self.center[0]-self.radius
            end_x = self.center[0]+self.radius
            start_y = self.center[1]-self.radius
            diff = float(abs(start_x - end_x))/(divisions-1)
            # Getting all possible x and y values.
            possible_x_values = [start_x + (times)*diff
                                 for times in range(divisions)]
            possible_y_values = [start_y + (times)*diff
                                 for times in range(divisions)]
            rad = diff*math.sqrt(2)/2
            subsystems = self._subsystems(possible_x_values, possible_y_values,
                                          rad, diff)
            self._actual_subdivision = subsystems

    def _separate_rods_by_coords(self, rods_list, possible_x_values,
                                    possible_y_values, rad, diff):
        """
            Separates rods by zones. It reduces the amount of steps that
        the programm makes.
        """
        x_0 = min(possible_x_values)
        y_0 = min(possible_y_values)
        x_length = len(possible_x_values)
        y_length = len(possible_y_values)
        divisions = len(possible_x_values)
        output = [[[] for dummy_ in range(divisions)] for dummy_2 in range(divisions)]
        for rod in rods_list:
            index_x = int(float(rod.x_mid - x_0)/diff)
            index_y = int(float(rod.y_mid - y_0)/diff)
            output[index_y][index_x].append(rod)
            if index_x > 0:
                output[index_y][index_x-1].append(rod)
            if index_x < x_length-1:
                output[index_y][index_x+1].append(rod)
            if index_y > 0:
                output[index_y-1][index_x].append(rod)
            if index_y < y_length-1:
                output[index_y+1][index_x].append(rod)
        return output
            

    def _subsystems(self, possible_x_values, possible_y_values, rad, diff):
        """
            Creates subsystems
        """
        subsystems = []
        rods_list = copy.deepcopy(list(self._rods))
        rods_matrix = self._separate_rods_by_coords(rods_list, possible_x_values,
                                                    possible_y_values, rad, diff)
        y_vals = range(len(possible_y_values))
        x_vals = range(len(possible_x_values))
        for index_y in y_vals:
            actual_y = possible_y_values[index_y]
            for index_x in x_vals:
                actual_x = possible_x_values[index_x]
                center = (actual_x, actual_y)
                distance = methods.distance_between_points(center,
                                                     self.center)
                rods = rods_matrix[index_y][index_x]
                subsystem = SubsystemState(center, rad, self.zone_coords,
                                           rods, self._kappas,
                                           self._allowed_kappa_error)
                if distance < self._radius:
                    subsystem.check_rods()
                else:
                    subsystem.remove_all_rods()
                subsystems.append(subsystem)
        return subsystems

    def _compute_density_matrix(self, divisions, normalized=False,
                               divided_by_area=False):
        """
            Computes density matrix of the system.
        """
        self.divide_in_circles(divisions)
        density = []
        for subsystem in self._actual_subdivision:
            subdensity = [subsystem.center[0], subsystem.center[1]]
            dens = subsystem.density
            if divided_by_area:
                dens /= subsystem.area
            if normalized:
                dens /= self.number_of_rods
            subdensity.append(dens)
            density.append(subdensity)
        self._density_matrix = density

    def plottable_density_matrix(self, divisions=50):
        """
            Returns 3 arrays: one for x, another for y and the last for values.
        """
        if len(self._density_matrix) == 0:
            self._compute_density_matrix(divisions=divisions)
        x_values = []
        y_values = []
        z_values = []
        for row in self._density_matrix:
            x_values.append(row[0])
            y_values.append(row[1])
            z_values.append(row[2])
        return x_values, y_values, z_values

    def _create_subgroups_matrix(self, divisions):
        """
            Put subsystems in a matrix form.
        """
        if divisions != self._divisions:
            self._reset()
        if not len(self._subdivision_centers):
            self.divide_in_circles(divisions)
            act_sub = self._actual_subdivision
            actual_y = act_sub[0].center[1]
            row = []
            subgroups_matrix = []
            for index in range(len(act_sub)):
                element = act_sub[index]
                element_y = element.center[1]
                if element_y != actual_y:
                    subgroups_matrix.append(row)
                    row = []
                    actual_y = element_y
                row.append(element)
            subgroups_matrix.append(row)
            self._subdivision_centers = subgroups_matrix

    def subgroups_matrix(self, divisions):
        """
            Returns subgroups matrix
        """
        self._create_subgroups_matrix(divisions)
        return self._subdivision_centers

    def _compute_g2_and_g4(self):
        """
            Computes correlation_g2 and correlation_g4 values
        """
        cos2_av, sin2_av, cos4_av, sin4_av = 0, 0, 0, 0
        for rod_ in list(self._rods):
            angle = rod_.angle*math.pi/180.0
            cos2_av += math.cos(2*angle)/self.number_of_rods
            sin2_av += math.sin(2*angle)/self.number_of_rods
            cos4_av += math.cos(4*angle)/self.number_of_rods
            sin4_av += math.sin(4*angle)/self.number_of_rods
        self._correlation_g2 = math.sqrt(cos2_av**2+sin2_av**2)
        self._correlation_g4 = math.sqrt(cos4_av**2+sin4_av**2)

    @property
    def correlation_g2(self):
        """
            sqrt(<cos(2*angle)>^2+<sin(2*angle)>^2)
        """
        if not self._correlation_g2:
            self._compute_g2_and_g4()
        return self._correlation_g2

    @property
    def correlation_g4(self):
        """
            sqrt(<cos(4*angle)>^2+<sin(4*angle)>^2)
        """
        if not self._correlation_g4:
            self._compute_g2_and_g4()
        return self._correlation_g4

    def _compute_g2_g4_matrices(self, divisions):
        """
            Computes correlation_g2 and correlation_g4 matrices for subgroups.
        """
        self.divide_in_circles(divisions)
        if not self._correlation_g2 or not self._correlation_g2:
            for subsystem in self._actual_subdivision:
                correlation_g2 = [subsystem.center[0], subsystem.center[1]]
                correlation_g4 = [subsystem.center[0], subsystem.center[1]]
                correlation_g2.append(subsystem.correlation_g2)
                correlation_g4.append(subsystem.correlation_g4)
                self._correlation_g2_subsystems.append(correlation_g2)
                self._correlation_g4_subsystems.append(correlation_g4)

    def correlation_g2_plot_matrix(self, divisions):
        """
            Returns values for plotting correlation_g2 matrix.
        """
        self._compute_g2_g4_matrices(divisions)
        x_values = []
        y_values = []
        z_values = []
        for subsystem in self._correlation_g2_subsystems:
            x_values.append(subsystem[0])
            y_values.append(subsystem[1])
            z_values.append(subsystem[2])
        return x_values, y_values, z_values

    def correlation_g4_plot_matrix(self, divisions):
        """
            Returns values for plotting correlation_g2 matrix.
        """
        self._compute_g2_g4_matrices(divisions)
        x_values = []
        y_values = []
        z_values = []
        for subsystem in self._correlation_g4_subsystems:
            x_values.append(subsystem[0])
            y_values.append(subsystem[1])
            z_values.append(subsystem[2])
        #return self._transform_for_pcolor(z_values, rad)
        return x_values, y_values, z_values

    @property
    def average_kappa(self):
        """
            Returns kappa average of group.
        """
        if not self._average_kappa:
            self._average_kappa = 0
            for rod_ in list(self._rods):
                self._average_kappa += rod_.kappa
            self._average_kappa /= self.number_of_rods
        return self._average_kappa

    @property
    def kappa_dev(self):
        """
            Returns sqrt(<kappa^2> - <kappa>^2)
        """
        if not self._kappa_dev:
            kappa2 = 0
            for rod_ in list(self._rods):
                kappa2 += rod_.kappa**2
            kappa2 /= self.number_of_rods
            self._kappa_dev = math.sqrt(kappa2-self.average_kappa**2)
        return self._kappa_dev

    @property
    def average_angle(self):
        """
            Returns average angle of the system (if exists).
        """
        if not self._average_angle:
            if self.correlation_g2 > 0.5 and self.correlation_g4 < 0.3:
                angle = 0
                for rod_ in list(self._rods):
                    angle2 = rod_.angle
                    # angle = angle + pi (simetry)
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
        self.divide_in_circles(divisions)
        for subsystem in self._actual_subdivision:
            row = [subsystem.center[0], subsystem.center[1]]
            row.append(subsystem.average_angle)
            self._angle_matrix.append(row)

    def plottable_average_angle_matrix(self, divisions):
        """
            Returns a plottable average angle matrix.
        """
        self._compute_average_angle_matrix(divisions)
        x_values = []
        y_values = []
        z_values = []
        for subsystem in self._angle_matrix:
            x_values.append(subsystem[0])
            y_values.append(subsystem[1])
            z_values.append(subsystem[2])
        #return self._transform_for_pcolor(z_values, rad)
        return x_values, y_values, z_values

    @property
    def angle_histogram(self):
        """
            Returns all angles in a list to make an histogram.
        """
        output = [rod_.angle for rod_ in list(self._rods)]
        return output

    def _get_closest_rod(self, rod_):
        """
            Returns closest rod in group to given rod.
        """
        distance = 1e100
        selected_rod = rod_
        for rod2 in list(self._rods):
            if rod_ != rod2:
                new_distance = rod_.distance_to_rod(rod2)
                if new_distance < distance:
                    distance = new_distance
                    selected_rod = rod2
        return selected_rod

    def _get_cluster_members(self, reference_rod,
                                    max_distance, max_angle_diff):
        """
            Gets the closest neighbour to a rod that fulfill
        some conditions.
        Angles in grad.
        """
        rods = set([])
        self.fill_dicts()
        if self._cluster_checked_dict[reference_rod.identifier]:
            return rods
        self._cluster_checked_dict[reference_rod.identifier] = True
        for rod_ in list(self._rods):
            if self._cluster_checked_dict[rod_.identifier]:
                continue
            x_diff = rod_.x_mid-reference_rod.x_mid
            y_diff = rod_.y_mid-reference_rod.y_mid
            distance = math.sqrt(x_diff**2+y_diff**2)
            angle_diff = abs(rod_.angle-reference_rod.angle)
            angle_diff = min([angle_diff, 180-angle_diff])
            if angle_diff <= max_angle_diff and distance < max_distance:
                subrods = self._get_cluster_members(rod_, max_distance,
                                                    max_angle_diff)
                rods.add(rod_)
                rods |= subrods
                continue
            cond1 = distance <= 1.4*max_distance
            cond2 = angle_diff <= max_angle_diff/2.0
            if cond1 and cond2:
                subrods = self._get_cluster_members(rod_, max_distance,
                                                   max_angle_diff)
                rods.add(rod_)
                rods |= subrods
        return rods

    def clusters(self, max_distance=None, max_angle_diff=None):
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
        cond2 = cond2 and (self._clusters_max_angle_diff != max_angle_diff)
        if not len(self._clusters) or cond2:
            self._clusters_max_distance = max_distance
            self._clusters_max_angle_diff = max_angle_diff
            if len(self._cluster_checked_dict.keys()):
                self.fill_dicts()
            clusters = []
            for rod_ in self._rods:
                if self._cluster_checked_dict[rod_.identifier]:
                    continue
                cluster = self._get_cluster_members(rod_,
                                max_distance, max_angle_diff)
                if cluster:
                    clusters.append(cluster)
            self._clusters = methods.erase_length_one_elements(clusters)
        return self._clusters

    def average_cluster_rod_num(self, max_distance=None, max_angle_diff=None):
        """
            Gets the average number of rods in clusters.
        Angles in grad.
        """
        lengths = self.number_of_rods_in_cluster(max_distance, max_angle_diff)
        try:
            return float(sum(lengths))/len(lengths)
        except ZeroDivisionError:
            print "No clusters detected."

    def number_of_rods_in_cluster(self, max_distance=None, max_angle_diff=None):
        """
            Creates a list with the number of rods in each cluster.
        Angles in grad.
        """
        lengths = []
        clusters = self.clusters(max_distance, max_angle_diff)
        for cluster in clusters:
            lengths.append(len(cluster))
        return lengths

    def number_of_clusters(self, max_distance=None, max_angle_diff=None):
        """
            Returns the number of clusters in the system.
        Angles in grad.
        """
        clusters = self.clusters(max_distance, max_angle_diff)
        return len(clusters)

    def total_area_of_clusters(self, max_distance=None, max_angle_diff=None):
        """
            Returns the area covered by clusters.
        """
        clusters = self.clusters(max_distance, max_angle_diff)
        rod_area = self.rod_area
        rods_num = self.number_of_rods_in_cluster(max_distance=max_distance,
                                               max_angle_diff=max_angle_diff)
        rods_num = sum(rods_num)
        if not rods_num or not rod_area:
            return 0
        total_area = rods_num*self.rod_area
        return total_area

    @property
    def rod_area(self):
        """
        Returns the area of a rod of this system.
        """
        first_rod = list(self._rods)[0]
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
        for rod_ in list(self._rods):
            new_row = [rod_]
            closest_rod = function(rod_)
            if not closest_rod:
                continue
            new_row.append(closest_rod)
            closest_rod_matrix.append(new_row)
        self._closest_rod_matrix = closest_rod_matrix


    def closest_rod_matrix(self):
        """
            Returns closest rod matrix.
        See _compute_closest_rod_matrix for more info.
        """
        if len(self._closest_rod_matrix) == 0:
            self._compute_closest_rod_matrix()
        return self._closest_rod_matrix

    def closest_rod_dict(self):
        """
            Returns a dictionary of closest rods.
        """
        closest_rod_matrix = self.closest_rod_matrix
        dictionary = {}
        for pair in closest_rod_matrix:
            dictionary[pair[0].identifier] = pair[1]
        return dictionary

    @property
    def relative_g2(self):
        """
            sqrt(<cos(2*angle)>^2+<sin(2*angle)>^2)
        Now angle is the relative between rods.
        N^2
        """
        if self._relative_g2 == None:
            if not self.number_of_rods:
                return 0
            sin = 0
            cos = 0
            for row in self.closest_rod_matrix:
                angle = math.radians(row[0].angle_between_rods(row[1]))
                sin += math.sin(2*angle)
                cos += math.cos(2*angle)
            sin /= self.number_of_rods
            cos /= self.number_of_rods
            self._relative_g2 = math.sqrt(sin**2+cos**2)
        return self._relative_g2

    @property
    def relative_g4(self):
        """
            sqrt(<cos(4*angle)>^2+<sin(4*angle)>^2)
        Now angle is the relative between rods.
        N^2
        """
        if self._relative_g4 == None:
            if not self.number_of_rods:
                return 0
            sin = 0
            cos = 0
            for row in self.closest_rod_matrix:
                angle = math.radians(row[0].angle_between_rods(row[1]))
                sin += math.sin(4*angle)
                cos += math.cos(4*angle)
            sin /= self.number_of_rods
            cos /= self.number_of_rods
            self._relative_g4 = math.sqrt(sin**2+cos**2)
        return self._relative_g4

    def _compute_relative_g2_g4_mat(self, divisions):
        """
            Computes correlation_g2 and correlation_g4 matrices for subgroups.
        """
        self.divide_in_circles(divisions)
        len_correlation_g2 = len(self._relative_g2_subsystems)
        len_correlation_g4 = len(self._relative_g4_subsystems)
        if  len_correlation_g2 == 0 or len_correlation_g4 == 0:
            for subsystem in self._actual_subdivision:
                correlation_g2 = [subsystem.center[0], subsystem.center[1]]
                correlation_g4 = [subsystem.center[0], subsystem.center[1]]
                correlation_g2.append(subsystem.correlation_g2)
                correlation_g4.append(subsystem.correlation_g4)
                self._relative_g2_subsystems.append(correlation_g2)
                self._relative_g4_subsystems.append(correlation_g4)

    def relative_g2_plot_matrix(self, divisions):
        """
            Returns values for plotting correlation_g2 matrix.
        """
        self._compute_relative_g2_g4_mat(divisions)
        x_values = []
        y_values = []
        z_values = []
        for subsystem in self._relative_g2_subsystems:
            x_values.append(subsystem[0])
            y_values.append(subsystem[1])
            z_values.append(subsystem[2])
        #return self._transform_for_pcolor(z_values, rad)
        return x_values, y_values, z_values

    def relative_g4_plot_matrix(self, divisions):
        """
            Returns values for plotting correlation_g2 matrix.
        """
        self._compute_relative_g2_g4_mat(divisions)
        x_values = []
        y_values = []
        z_values = []
        for subsystem in self._relative_g4_subsystems:
            x_values.append(subsystem[0])
            y_values.append(subsystem[1])
            z_values.append(subsystem[2])
        #return self._transform_for_pcolor(z_values, rad)
        return x_values, y_values, z_values

    @property
    def average_angle_using_matrix(self):
        """
            Returns average angle of the system using direction matrices.
        Value in radians.
        """
        if len(self._direction_matrix) == 0:
            self._direction_matrix = matrix.zeros(2, 2)
            for rod_ in list(self._rods):
                self._direction_matrix += rod_.direction_matrix
        eigen1, dummy_ = self._direction_matrix.diagonalize_2x2()
        return eigen1
















class SubsystemState(SystemState):
    """
        Group of rods. Used to put all rods that are in a zone or
    have something in common.
    """

    def __init__(self, center, rad, zone_coords, rods, kappas, allowed_kappa_error):
        """
            Initialization
        """
        self._subsystem_coords = (center[0], center[1], rad)
        SystemState.__init__(self, rods=rods, zone_coords=self._subsystem_coords,
                            kappas=kappas, allowed_kappa_error=allowed_kappa_error)
        self._is_subsystem = True
        self._center = center
        self._rad = rad
        self._main_center = (zone_coords[0], zone_coords[1])
        self._main_rad = zone_coords[2]
        self._position_rad = methods.distance_between_points(self._main_center,
                                            self._center)
        self._area = methods.effective_area(self._rad,
                                    self._position_rad, self._main_rad)

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
        return self._center

    @property
    def radius(self):
        """
            Radius of the subsystem
        """
        return self._rad

    @property
    def area(self):
        """
            Area of the subsystem.
        """
        return self._area

    def _update_density(self):
        """
            Computes density of the group.
        """
        density = 0
        for rod_ in self:
            density += rod_.kappa
        if not self.area:
            self._density = 0
        else:
            self._density = float(density)/self.area

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
        zone_coords = [coord for coord in self._subsystem_coords]
        for rod_ in self._rods:
            valid = rod_.is_valid_rod(self._kappas,
                        self._allowed_kappa_error,
                        zone_coords)
            if not valid:
                self._remove_rod(rod_)
        self._reset()




def create_rods(folder="./", kappas=10, allowed_kappa_error=.3,
                radius_correction_ratio=0.1):
    """
    Create one rod for each rod_data and for each file
    returns [RodGroup1, RodGroup2, ...]
    """
    names, files = methods.import_files(folder=folder)
    if len(files) == 0:
        print "No files to import."
        raise ValueError
    states = [None for dummy_ in range(len(files))]
    processes = []
    states_queue = mp.Queue()
    for index in range(len(files)):
        process = mp.Process(target=create_rods_process,
                            args=(kappas, allowed_kappa_error,
                            radius_correction_ratio, names,
                            files, index, states_queue))
        processes.append(process)
    running, processes_left = methods.run_processes(processes)
    num_processes = len(running)
    finished = 0
    while finished < num_processes:
        finished += 1
        [index, state] = states_queue.get()
        states[index] = state
        if len(processes_left):
            finished -= 1
            new_process = processes_left.pop()
            new_process.start()
    return names, states




def create_rods_process(kappas, allowed_kappa_error,
                        radius_correction_ratio, names,
                        files, index, states_queue):
    """
    Process of method.
    """
    state = SystemState(kappas, allowed_kappa_error,
               radius_correction_ratio, names[index])
    data = methods.import_data(files[index])
    for dataline in data:
        try:
            parameters = tuple(dataline)
            new_rod = rod.Rod(parameters)
            state._put_rod(new_rod)
        except ValueError:
            print names[index]
    state.check_rods()
    states_queue.put([index, state])

