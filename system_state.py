"""
    System State. Group of rods in a determined moment.
"""

import math
from methods import is_in_circle, same_area_rad
from queue import Queue


class SystemState(object):
    """
    Group of rods in a moment.
    Each image has to be translated into a RodGroup (by this class?)
    """
    def __init__(self, kappas=10, allowed_kappa_error=.5,
                radius_correction_ratio=0,
                id_string="", zone_coords=None):
        """
        Initialization
        """
        self._rods = Queue()
        self._number_of_rods = 0
        self._kappas = kappas
        self._allowed_kappa_error = allowed_kappa_error
        self._radius_correction_ratio = radius_correction_ratio
        self.id_string = id_string
        self._rad_for_division = None
        self._actual_subdivision = []
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
        try:
            self._radius = zone_coords[2]
            self._center_x = zone_coords[0]
            self._center_y = zone_coords[1]
            self._zone_coords = zone_coords
            self._fixed_center_radius = True
        except TypeError:
            self._fixed_center_radius = False
            self._zone_cords = []

    @property
    def number_of_rods(self):
        """
        Returns the number of rods.
        """
        return self._number_of_rods

    def add_rod(self, rod):
        """
        Adds a rod to the group
        """
        self._rods.join(rod)
        self._number_of_rods += 1
        self.reset()

    def get_rod(self):
        """
        Returns the first rod in the queue
        The rod is removed of the group!
        """
        self._number_of_rods -= 1
        self.reset()
        return self._rods.get_next()

    def remove_rod(self, rod):
        """
        Removes a rod from the group (queue object mod needed)
        """
        self._rods.delete(rod)
        self._number_of_rods -= 1
        self.reset()

    def reset(self):
        """
        Called when system is changed..
        Reset all important values, so they must be
        computed again.
        """
        self._rad_for_division = None
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
        self._relative_g4 = None
        self._closest_rod_matrix = []
        self._relative_g2_subsystems = []
        self._relative_g4_subsystems = []
        if not self._fixed_center_radius:
            self._radius = None
            self._center_x = None
            self._center_y = None
            self._zone_coords = []

    def _compute_center_and_radius(self):
        """
        Computes where the center of the system is and its
        radius.
        """

            #There must be rods to make statistics.
        if self._number_of_rods == 0:
            msg = "center_and_radius can't be computed before adding rods"
            raise ValueError(msg)
        if not len(self._zone_coords) and not self._fixed_center_radius:
            x_values = []
            y_values = []
            for rod in list(self._rods):
                x_values.append(rod.x_mid)
                y_values.append(rod.y_mid)
            #center is the mean position of all particles.
            center_x = sum(x_values)*1.0/self._number_of_rods
            center_y = sum(y_values)*1.0/self._number_of_rods
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
        self._compute_center_and_radius()
        for rod in list(self._rods):
            valid = rod.is_valid_rod(self._kappas,
                        self._allowed_kappa_error,
                        self.zone_coords)
            if not valid:
                self.remove_rod(rod)

    def _divide_in_circles(self, rad):
        """
        Divides rods into groups contained in circles of given rad.
        """
        if self._rad_for_division == rad:
            return
        if (rad < 0) or (rad > self.radius):
            print "Use a correct radius (0<rad<main_rad)"
            raise ValueError
        self._rad_for_division = rad
        # Defining zone and distance between points.
        diff = int(rad/2)
        start_x = self.center[0]-self.radius
        end_x = self.center[0]+self.radius
        start_y = self.center[1]-self.radius
        end_y = self.center[1]+self.radius
        # Getting all possible x and y values.
        max_times = int((end_x-start_x)/diff+1)
    
        possible_x_values = [start_x + times*diff
                             for times in range(max_times)]
        max_times = int((end_y-start_y)/diff+1)
        possible_y_values = [start_y + times*diff
                             for times in range(max_times)]
        subsystems = self._subsystems(possible_x_values, possible_y_values,
                                      rad)
        self._actual_subdivision = subsystems

    def _subsystems(self, possible_x_values, possible_y_values, rad):
        """
        Creates subsystems
        """
        subsystems = []
        for actual_x in possible_x_values:
            for actual_y in possible_y_values:
                if is_in_circle(actual_x, actual_y,
                                self.center[0], self.center[1],
                                self.radius-rad):
                    x_diff = abs(actual_x-self.center[0])
                    y_diff = abs(actual_y-self.center[1])
                    pos_rad = math.sqrt(x_diff**2+y_diff**2)
                    same_area_radius = same_area_rad(rad, pos_rad, self.radius)
                    subsystem = SubsystemState((actual_x, actual_y),
                                               same_area_radius, math.pi*rad**2)
                    subsystem.add_rods(list(self._rods))
                    subsystems.append(subsystem)
        return subsystems

    def _compute_density_matrix(self, rad=100, normalized=False,
                               divided_by_area=False):
        """
        Computes density matrix of the system.
        """
        if self._rad_for_division == rad:
            return self._density_matrix
        self._divide_in_circles(rad)
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

    def plottable_density_matrix(self, rad):
        """
        Returns 3 arrays: one for x, another for y and the last for values.
        """
        if len(self._density_matrix) == 0:
            self._compute_density_matrix(rad=rad)
        x_values = []
        y_values = []
        z_values = []
        for row in self._density_matrix:
            x_values.append(row[0])
            y_values.append(row[1])
            z_values.append(row[2])
        return x_values, y_values, z_values

    def _compute_g2_and_g4(self):
        """
        Computes correlation_g2 and correlation_g4 values
        """
        cos2_av, sin2_av, cos4_av, sin4_av = 0, 0, 0, 0
        for rod in list(self._rods):
            angle = rod.angle*math.pi/180.0
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

    def _compute_g2_g4_matrices(self, rad):
        """
        Computes correlation_g2 and correlation_g4 matrices for subgroups.
        """
        if not self._correlation_g2 or not self._correlation_g2:
            self._divide_in_circles(rad)
            for subsystem in self._actual_subdivision:
                correlation_g2 = [subsystem.center[0], subsystem.center[1]]
                correlation_g4 = [subsystem.center[0], subsystem.center[1]]
                correlation_g2.append(subsystem.correlation_g2)
                correlation_g4.append(subsystem.correlation_g4)
                self._correlation_g2_subsystems.append(correlation_g2)
                self._correlation_g4_subsystems.append(correlation_g4)

    def correlation_g2_plot_matrix(self, rad):
        """
        Returns values for plotting correlation_g2 matrix.
        """
        self._compute_g2_g4_matrices(rad)
        xval = []
        yval = []
        zval = []
        for subsystem in self._correlation_g2_subsystems:
            xval.append(subsystem[0])
            yval.append(subsystem[1])
            zval.append(subsystem[2])
        return xval, yval, zval

    def correlation_g4_plot_matrix(self, rad):
        """
        Returns values for plotting correlation_g2 matrix.
        """
        self._compute_g2_g4_matrices(rad)
        xval = []
        yval = []
        zval = []
        for subsystem in self._correlation_g4_subsystems:
            xval.append(subsystem[0])
            yval.append(subsystem[1])
            zval.append(subsystem[2])
        return xval, yval, zval

    @property
    def average_kappa(self):
        """
        Returns kappa average of group.
        """
        if not self._average_kappa:
            self._average_kappa = 0
            for rod in list(self._rods):
                self._average_kappa += rod.kappa
            self._average_kappa /= self.number_of_rods
        return self._average_kappa

    @property
    def kappa_dev(self):
        """
        Returns sqrt(<kappa^2> - <kappa>^2)
        """
        if not self._kappa_dev:
            kappa2 = 0
            for rod in list(self._rods):
                kappa2 += rod.kappa**2
            kappa2 /= self.number_of_rods
            self._kappa_dev = math.sqrt(kappa2-self.average_kappa**2)
        return self._kappa_dev

    @property
    def average_angle(self):
        """
        Returns average angle of the system (if exists).
        """
        if not self._average_angle:
            if self.correlation_g2 > 0.3 and self.correlation_g4 < 0.3:
                angle = 0
                for rod in list(self._rods):
                    angle2 = rod.angle
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

    def _compute_average_angle_matrix(self, rad):
        """
        Computes average angle matrix
        """
        if self._rad_for_division == rad and len(self._angle_matrix):
            return
        self._divide_in_circles(rad)
        for subsystem in self._actual_subdivision:
            row = [subsystem.center[0], subsystem.center[1]]
            row.append(subsystem.average_angle)
            self._angle_matrix.append(row)

    def plottable_average_angle_matrix(self, rad):
        """
        Returns a plottable average angle matrix.
        """
        self._compute_average_angle_matrix(rad)
        xval = []
        yval = []
        zval = []
        for subsystem in self._angle_matrix:
            xval.append(subsystem[0])
            yval.append(subsystem[1])
            zval.append(subsystem[2])
        return xval, yval, zval

    @property
    def angle_histogram(self):
        """
        Returns all angles in a list to make an histogram.
        """
        output = [rod.angle for rod in list(self._rods)]
        return output

    def _get_closest_rod(self, rod):
        """
        Returns closest rod in group to given rod.
        """
        distance = 1e100
        selected_rod = rod
        for rod2 in list(self._rods):
            if rod != rod2:
                new_distance = rod.distance_to_rod(rod2)
                if new_distance < distance:
                    distance = new_distance
                    selected_rod = rod2
        return selected_rod

    def _compute_closest_rod_matrix(self):
        """
        Creates closer rod matrix:
        [[rod1, closest_rod_to_rod1],
        [rod2, closest_rod_to_rod2],
        ...
        [rodN, closest_rod_to_rodN]]
        """
        closest_rod_matrix = []
        for rod in list(self._rods):
            new_row = [rod]
            new_row.append(self._get_closest_rod(rod))
            closest_rod_matrix.append(new_row)
        self._closest_rod_matrix = closest_rod_matrix


    @property
    def closest_rod_matrix(self):
        """
        Returns closest rod matrix.
        See _compute_closest_rod_matrix for more info.
        """
        if len(self._closest_rod_matrix) == 0:
            self._compute_closest_rod_matrix()
        return self._closest_rod_matrix

    @property
    def relative_g2(self):
        """
        sqrt(<cos(2*angle)>^2+<sin(2*angle)>^2)
        Now angle is the relative between rods.
        N^2
        """
        if self._relative_g2 == None:
            sin = 0
            cos = 0
            for row in self.closest_rod_matrix:
                angle = row[0].angle_between_rods(row[1])
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
            sin = 0
            cos = 0
            for row in self.closest_rod_matrix:
                angle = row[0].angle_between_rods(row[1])
                sin += math.sin(4*angle)
                cos += math.cos(4*angle)
            sin /= self.number_of_rods
            cos /= self.number_of_rods
            self._relative_g4 = math.sqrt(sin**2+cos**2)
        return self._relative_g4

    def _compute_relative_g2_g4_mat(self, rad):
        """
        Computes correlation_g2 and correlation_g4 matrices for subgroups.
        """
        self._divide_in_circles(rad)
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

    def relative_g2_plot_matrix(self, rad):
        """
        Returns values for plotting correlation_g2 matrix.
        """
        self._compute_relative_g2_g4_mat(rad)
        xval = []
        yval = []
        zval = []
        for subsystem in self._relative_g2_subsystems:
            xval.append(subsystem[0])
            yval.append(subsystem[1])
            zval.append(subsystem[2])
        return xval, yval, zval

    def relative_g4_plot_matrix(self, rad):
        """
        Returns values for plotting correlation_g2 matrix.
        """
        self._compute_relative_g2_g4_mat(rad)
        xval = []
        yval = []
        zval = []
        for subsystem in self._relative_g4_subsystems:
            xval.append(subsystem[0])
            yval.append(subsystem[1])
            zval.append(subsystem[2])
        return xval, yval, zval

    def _cluster_finder(self, reference_rod, max_distance, max_angle_diff):
        """
        Gets the closest neighbour to a rod that fulfill some conditions.
        """
        selected_rod = None
        min_distance = max_distance
        for rod in self._rods:
            if rod == reference_rod:
                continue
            x_diff = rod.x_mid-reference_rod.x_mid
            y_diff = rod.y_mid-reference_rod.y_mid
            distance = math.sqrt(x_diff**2+y_diff**2)
            angle_diff = abs(rod.angle-reference_rod.angle)
            angle_diff = min([angle_diff, 180-angle_diff])
            if angle_diff <= max_angle_diff and distance < min_distance:
                selected_rod = rod
                min_distance = distance
            if distance <= 1.4*max_distance and angle_diff <= max_angle_diff/2.0:
                slope = y_diff / x_diff
                alpha = math.atan(abs(slope))
                theta = (math.pi*reference_rod.angle)/180
                if slope > 0:
                    beta = abs(alpha-theta)
                else:
                    beta = abs(math.pi+alpha-theta)
                perp_distance = distance*math.sin(beta)
                if perp_distance < min_distance:
                    selected_rod = rod
                    min_distance = perp_distance
        return selected_rod




class SubsystemState(SystemState):
    """
    Group of rods. Used to put all rods that are in a zone or
    have something in common.
    """

    def __init__(self, center, rad, area):
        """
        Initialization
        """
        SystemState.__init__(self)
        self._center = center
        self._rad = rad
        self._area = area

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
        self._density = self._number_of_rods

    @property
    def density(self):
        """
        Returns the density of the group.
        """
        if not self._density:
            self._update_density()
        return self._density

    def add_rods(self, rod_list):
        """
        Add all rods of the list that are inside the circle.
        """
        try:
            for rod in rod_list:
                if rod.is_in_circle(self.center, self.radius):
                    self.add_rod(rod)
        except TypeError:
            print "Use a rod list in add_rods method"

