"""
        System State. Group of rods in a determined moment.
"""
import math
import queue
import Queue
import matrix
import multiprocessing as mp
import re
import os

#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################

class Rod(object):
    """
    Rod object.
    """

    def __init__(self, (ID, area, xm, ym, major, minor,
                        angle, feret, feretx, ferety,
                        feretangle, minferet, xstart, ystart)):
        """
        Initialization of rod
        """                                     #Column
        self._id = int(ID)                      #0
        self._area = float(area)                #1
        self._x_mid = float(xm)                 #2
        self._y_mid = float(ym)                 #3
        self._major = float(major)              #4
        self._minor = float(minor)              #5
        self._angle = float(angle)              #6
        self._feret = float(feret)              #7
        self._feret_x = float(feretx)           #8
        self._feret_y = float(ferety)           #9
        self._feret_angle = float(feretangle)   #10
        self._min_feret = float(minferet)       #11
        self._x_start = float(xstart)           #12
        self._y_start = float(ystart)           #13
        self._hash = 0
        self._direction_matrix = matrix.zeros(2, 2)
        self._kappa = float(self.feret)/self.min_feret

    @property
    def feret(self):
        """
        Feret length.
        (wikipedia) The Feret diameter or Feret's diameter is a measure
        of an object size along a specified direction. In general, it can
        be defined as the distance between the two parallel planes restricting
        the object perpendicular to that direction. It is therefore also called
        the caliper diameter, referring to the measurement of the object size
        with a caliper. This measure is used in the analysis of particle sizes,
        for example in microscopy, where it is applied to projections of a
        three-dimensional (3D) object on a 2D plane. In such cases, the Feret
        diameter is defined as the distance between two parallel tangential
        lines rather than planes.[1][2]
        """
        return self._feret

    @property
    def center(self):
        """
        Returns position of the center of the rod.
        """
        return self._x_mid, self._y_mid

    def __eq__(self, rod2):
        """
        Check if a rod is the same as another rod.
        Rods must be of the same group.
        """
        return self.hash_ == rod2.hash_

    def __ne__(self, rod2):
        """
        != magic method
        """
        return not self == rod2

    def __repr__(self):
        """
        String transformation.
        """
        output = ""
        output += "id: "+str(self.identifier)+"\n"
        output += "center: "+str(self.center)+"\n"
        output += "angle: "+str(self.angle)+"\n"
        return output

    @property
    def identifier(self):
        """
        Returns an identification number.
        """
        return self._id

    @property
    def hash_(self):
        """
        Returns an unique number of this rod.
        Uses some of rod properties.
        """
        output = ""
        output += str(self.identifier)
        output += str(int(self.min_feret))
        output += str(int(self.x_mid))
        output += str(int(self.y_mid))
        output += str(int(self.kappa))
        return int(output)


    @property
    def min_feret(self):
        """
        Minimum Feret length.
        """
        return self._min_feret

    @property
    def x_mid(self):
        """
        Average x of rod.
        """
        return self._x_mid

    @property
    def y_mid(self):
        """
        Average y of rod.
        """
        return self._y_mid

    @property
    def kappa(self):
        """
        L/D of rod.
        """
        return self._kappa

    @property
    def angle(self):
        """
        Angle of rod.
        """
        return self._angle

    def is_in_circle(self, center, rad):
        """
        Checks if rod is in the circle defined by the given center and
        the given rad.
        """
        return is_in_circle(self.x_mid, self.y_mid,
                            center[0], center[1], rad)

    def has_valid_proportions(self, kappas, allowed_error):
        """
        Checks if rod has valid L/D (kappas are possibles values
        for L/D).
        """
        passed = []
        try:
            for kappa in kappas:
                condition = abs(self.kappa-kappa) < allowed_error
                passed.append(condition)
#kappa is not an array, so there is only 1 kappa.
        except TypeError:
            condition = abs(self.kappa-kappas) < allowed_error
            passed.append(condition)
        output = False
        for condition in passed:
            output = output or condition
        return output

    def is_valid_rod(self, kappas,
                    allowed_kappa_error,
                    zone_coords):
        """
        Check if rod is valid checking L/D and distance to center.
        TODO: If rods are near, kappa is not correct.
        """
        center = (zone_coords[0], zone_coords[1])
        is_in_main = self.is_in_circle(center, zone_coords[2])
        has_valid_proportions = self.has_valid_proportions(kappas,
                                                           allowed_kappa_error)
        output = is_in_main and has_valid_proportions
        return output

    def distance_to_rod(self, rod):
        """
        Returns the distance to another rod.
        """
        diff_x = abs(self.x_mid-rod.x_mid)
        diff_y = abs(self.y_mid-rod.y_mid)
        return math.sqrt(diff_x**2+diff_y**2)

    def angle_between_rods(self, rod):
        """
        Returns value of angle that formes this rod with another.
        """
        angle1 = abs(self.angle-rod.angle)
        angle2 = 180 - angle1
        return min(angle1, angle2)

    @property
    def direction_matrix(self):
        """
        Returns a matrix with the form:
        ex^2-1  ex*ey
        ex*ey   ey^2-1
        """
        if self._direction_matrix == matrix.zeros(2, 2):
            e_x = math.cos(self.angle)
            e_y = math.sin(self.angle)
            self._direction_matrix[0][0] = 2*e_x**2-1
            self._direction_matrix[1][1] = 2*e_y**2-1
            self._direction_matrix[1][0] = e_x*e_y
            self._direction_matrix[0][1] = e_x*e_y
        return self._direction_matrix

#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################

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
        self._rods = queue.Queue(rods)
        self._number_of_rods = len(self._rods)
        self._kappas = kappas
        self._allowed_kappa_error = allowed_kappa_error
        self._radius_correction_ratio = radius_correction_ratio
        self.id_string = id_string
        self._rad_of_division = None
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
        self._clusters = []
        self._rods_dict = {}
        self._cluster_checked_dict = {}
        self._clusters_max_distance = None
        self._clusters_max_angle_diff = None
        try:
            self._radius = zone_coords[2]
            self._center_x = zone_coords[0]
            self._center_y = zone_coords[1]
            self._zone_coords = zone_coords
            self._fixed_center_radius = True
        except TypeError:
            self._fixed_center_radius = False
            self._zone_cords = []
        except IndexError:
            print zone_coords
            raise IndexError

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
            msg += " " + str(self._rods_dict.keys())    #There are not keys! -> Dict not filled.
            raise IndexError(msg)

    def get_rods_range(self, initial_id, final_id):
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
        for rod in self._rods:
            yield rod

    def __list__(self):
        """
        Returns a list of rods
        """
        output = [rod for rod in self]
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

    def put_rod(self, rod):
        """
            Adds a rod to the group
        """
        self._rods.put(rod)
        self._reset()

    def get_rod(self):
        """
            Returns the first rod in the queue
        """
        rod = self._rods.get()
        self._rods.put(rod)
        return rod

    def remove_rod(self, rod):
        """
            Removes a rod from the group (queue object mod needed)
        """
        self._rods.delete(rod)
        self._reset()

    def _reset(self):
        """
            Called when system is changed..
        reset all important values, so they must be
        computed again.
        """
        self._rad_of_division = None
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
            for rod in self._rods:
                identifier = rod.identifier
                self._rods_dict[identifier] = rod
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
        if not len(self._zone_coords) and not self._fixed_center_radius:
            x_values = []
            y_values = []
            for rod in list(self._rods):
                x_values.append(rod.x_mid)
                y_values.append(rod.y_mid)
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


    def check_rods(self, all_cores=False):
        """
            Check if rods are correct.
        """
        zone_coords = [coord for coord in self.zone_coords]
        for rod in self._rods:
            valid = rod.is_valid_rod(self._kappas,
                        self._allowed_kappa_error,
                        zone_coords)
            if not valid:
                self.remove_rod(rod)
        

    def _divide_in_circles(self, rad):
        """
            Divides rods into groups contained in circles of given rad.
        """
        if self._rad_of_division == rad:
            return
        if (rad < 0) or (rad > self.radius):
            print "Use a correct radius (0<rad<main_rad)"
            raise ValueError
        self._rad_of_division = rad
        # Defining zone and distance between points.
        diff = 2*rad
        start_x = self.center[0]-self.radius
        end_x = self.center[0]+self.radius
        start_y = self.center[1]-self.radius
        end_y = self.center[1]+self.radius
        # Getting all possible x and y values.
        max_times = int(float(end_x-start_x)/diff+1)

        possible_x_values = [start_x + times*diff
                             for times in range(max_times)]
        max_times = int(float(end_y-start_y)/diff+1)
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
        for actual_y in possible_y_values:
            for actual_x in possible_x_values:
                if is_in_circle(actual_x, actual_y,
                                self.center[0], self.center[1],
                                self.radius-rad/2):
                    x_diff = abs(actual_x-self.center[0])
                    y_diff = abs(actual_y-self.center[1])
                    pos_rad = math.sqrt(x_diff**2+y_diff**2)
                    same_area_radius = same_area_rad(rad, pos_rad, self.radius)
                    subsystem = SubsystemState((actual_x, actual_y),
                                               same_area_radius, math.pi*rad**2)
                    subsystem.put_rods(list(self._rods))
                    subsystems.append(subsystem)
        return subsystems

    def _compute_density_matrix(self, rad, normalized=False,
                               divided_by_area=False):
        """
            Computes density matrix of the system.
        """
        if self._rad_of_division != rad:
            self._reset()
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
        #return self._transform_for_pcolor(z_values, rad)
        return x_values, y_values, z_values

    def _transform_for_pcolor(self, z_values, rad):
        """
            Transform arrays to a plotable set of arrays.
        """
        xmat = []
        ymat = []
        zmat = []
        matrix_ = self.subgroups_matrix(rad)
        index = 0
        for row in range(len(matrix_)):
            xrow = []
            yrow = []
            zrow = []
            for col in range(len(matrix[row])):
                element = matrix[row][col]
                xrow.append(element.center[0]-rad)
                yrow.append(element.center[1]-rad)
                zrow.append(z_values[index])
            xrow.append(element.center[0]+rad)
            yrow.append(element.center[1]+rad)
            xmat.append(xrow)
            ymat.append(yrow)
            zmat.append(zrow)
        return xmat, ymat, zmat

    def _create_subgroups_matrix(self, rad):
        """
            Put subsystems in a matrix form.
        """
        if self._rad_of_division != rad:
            self._reset()
        self._divide_in_circles(rad)
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
        self._subdivision_centers = subgroups_matrix

    def subgroups_matrix(self, rad):
        """
            Returns subgroups matrix
        """
        self._create_subgroups_matrix(rad)
        return self._subdivision_centers



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
        if self._rad_of_division != rad:
            self._reset()
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
        x_values = []
        y_values = []
        z_values = []
        for subsystem in self._correlation_g2_subsystems:
            x_values.append(subsystem[0])
            y_values.append(subsystem[1])
            z_values.append(subsystem[2])
        #return self._transform_for_pcolor(z_values, rad)
        return x_values, y_values, z_values

    def correlation_g4_plot_matrix(self, rad):
        """
            Returns values for plotting correlation_g2 matrix.
        """
        self._compute_g2_g4_matrices(rad)
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
            if self.correlation_g2 > 0.5 and self.correlation_g4 < 0.3:
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
        if self._rad_of_division != rad:
            self._reset()
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
        for rod in list(self._rods):
            if self._cluster_checked_dict[rod.identifier]:
                continue
            x_diff = rod.x_mid-reference_rod.x_mid
            y_diff = rod.y_mid-reference_rod.y_mid
            distance = math.sqrt(x_diff**2+y_diff**2)
            angle_diff = abs(rod.angle-reference_rod.angle)
            angle_diff = min([angle_diff, 180-angle_diff])
            if angle_diff <= max_angle_diff and distance < max_distance:
                subrods = self._get_cluster_members(rod, max_distance,
                                                    max_angle_diff)
                rods.add(rod)
                rods |= subrods
                continue
            if distance <= 1.4*max_distance and angle_diff <= max_angle_diff/2.0:
                subrods = self._get_cluster_members(rod, max_distance,
                                                   max_angle_diff)
                rods.add(rod)
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
            for rod in self._rods:
                if self._cluster_checked_dict[rod.identifier]:
                    continue
                cluster = self._get_cluster_members(rod,
                                max_distance, max_angle_diff)
                if cluster:
                    clusters.append(cluster)
            assert len(clusters) > 0, "no clusters detected"
            self._clusters = erase_length_one_elements(clusters)
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
        for cluster in self.clusters(max_distance, max_angle_diff):
            lengths.append(len(cluster))
        return lengths

    def number_of_clusters(self, max_distance=None, max_angle_diff=None):
        """
            Returns the number of clusters in the system.
        Angles in grad.
        """
        return len(self.clusters(max_distance, max_angle_diff))

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
        for rod in list(self._rods):
            new_row = [rod]
            closest_rod = function(rod)
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

    def _compute_relative_g2_g4_mat(self, rad):
        """
            Computes correlation_g2 and correlation_g4 matrices for subgroups.
        """
        if self._rad_of_division != rad:
            self._reset()
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
        x_values = []
        y_values = []
        z_values = []
        for subsystem in self._relative_g2_subsystems:
            x_values.append(subsystem[0])
            y_values.append(subsystem[1])
            z_values.append(subsystem[2])
        #return self._transform_for_pcolor(z_values, rad)
        return x_values, y_values, z_values

    def relative_g4_plot_matrix(self, rad):
        """
            Returns values for plotting correlation_g2 matrix.
        """
        self._compute_relative_g2_g4_mat(rad)
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
            for rod in list(self._rods):
                self._direction_matrix += rod.direction_matrix
        eigen1, dummy_ = self._direction_matrix.diagonalize_2x2()
        return eigen1

#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################

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
        self._density = self.number_of_rods

    @property
    def density(self):
        """
            Returns the density of the group.
        """
        if not self._density:
            self._update_density()
        return self._density

    def put_rods(self, rod_list):
        """
            Add all rods of the list that are inside the circle.
        """
        try:
            for rod in rod_list:
                if rod.is_in_circle(self.center, self.radius):
                    self.put_rod(rod)
        except TypeError:
            print "Use a rod list in put_rods method"


#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#####################################################################
#########################################################################
# METHODS
#########################################################################



def segment_area(rad, min_dist):
    """
    Computes the area of an intersection of a circle with a line
    (the part that doesn't have the center)
    rad: radius of the circle
    min_dist: minimum distance from small circle center to a line that joins
        both intersections of the circles.
    """
    min_dist = float(min_dist)
    if min_dist >= rad:
        return 0
    elif abs(min_dist) >= rad:
        return math.pi*rad**2
    phi = math.acos(abs(min_dist)/rad)
    #DEBUG#
    assert 0 <= phi <= math.pi/2, "Error in angle"
    #######
    if phi <= 1e-10:
        if min_dist > 0:
            return 0
        else:
            return math.pi*rad**2
    section = phi*rad**2
    #DEBUG#
    assert section >= 0, "segment_area: Section is negative"
    #######
    distance_between_intersections = 2*min_dist*math.tan(phi)
    #DEBUG#
    msg = "segment_area: distance between "
    msg += "intersections can't be greater than diameter"
    assert distance_between_intersections <= 2*rad, msg
    #######
    triangle_area = distance_between_intersections*min_dist/2.0
    #DEBUG#
    msg = "segment_area: Triangle area must be smaller than section area"
    msg += "\nRatio="+str(triangle_area*1.0/section)
    assert triangle_area < section, msg
    ######
    if min_dist >= 0:
        output = section - triangle_area
    else:
        output = math.pi*rad**2 - section + triangle_area
    #DEBUG#
    msg = "segment_area: Obtained area is negative. "
    msg += "Values: rad:"+str(rad)
    msg += " min_dist:"+str(min_dist)+" rat:"+str(min_dist/rad)
    msg += " phi:"+str(phi)+" area:"+str(output)
    assert output > 0, msg
    #######
    return output

#######################################################################
#######################################################################

def effective_area(small_rad, small_position_rad, main_rad):
    """
    Computes the area of the small circle intersected with main circle.
    There are errors in this method implementation.
    """
    # circle completely included in the bigger one
    if small_rad+small_position_rad <= main_rad:
        return math.pi*small_rad**2
    #DEBUG#
    assert small_position_rad <= main_rad, "Circle is outside the bigger one"
    #######
    min_dist = compute_min_dist(small_rad, small_position_rad, main_rad)
    if min_dist >= small_rad:
        return math.pi*small_rad**2
    elif abs(min_dist) >= small_rad:
        return 0
    min_dist_main = small_position_rad+min_dist
    correction = segment_area(main_rad, min_dist_main)
    #DEBUG#
    msg = "effective_area: Correction must be smaller than small circle's area"
    assert correction < math.pi*small_rad**2, msg
    #######
    section_area = segment_area(small_rad, min_dist)
    small_area = math.pi*small_rad**2
    #DEBUG#
    msg = "In the limit, h=-rad has to return total area"
    assert small_area == segment_area(small_rad, -small_rad), msg
    msg = "Correction too high: Ration: "+str(float(correction)/small_area)
    assert correction < small_area, msg
    #######
    output = math.pi*small_rad**2 - section_area + correction
    return output

#######################################################################
#######################################################################

def compute_min_dist(small_rad, small_position_rad, main_rad):
    """
    Computes the distance from small circle center to the line that joins both
    circles' intersections.
    """
    try:
        min_dist = (main_rad**2)-(small_position_rad**2)-(small_rad**2)
        min_dist = float(min_dist)
    except OverflowError:
        return small_rad*1.1
    min_dist /= (2*small_position_rad)
    if min_dist > small_rad:
        return small_rad
    if abs(min_dist) > small_rad:
        return -small_rad
    return min_dist

#######################################################################
#######################################################################

def same_area_rad(small_rad, small_position_rad,
                    main_rad, allowed_error_ratio=.2,
                    max_reps=10):
    """
    Computes a new radius. With that, effective area is the same small circle's.
    Better use binary search
    """
    # circle completely included in main
    if small_position_rad + small_rad <= main_rad:
        return small_rad
    # circle completely excluded of main
    if small_position_rad - small_rad > main_rad:
        print "External circle introduced"
        raise ValueError
    wanted_area = math.pi*small_rad**2
    allowed_error = wanted_area * allowed_error_ratio
    low_rad = small_rad
    high_rad = small_rad*10
    def area(rad, small_position_rad=small_position_rad, main_rad=main_rad):
        """
        Needed function for binary search, as only 1 arg is allowed.
        """
        try:
            return effective_area(rad, small_position_rad, main_rad)
        except OverflowError:
            return wanted_area*10
    actual_area = area(high_rad)
    while actual_area < wanted_area:
        high_rad *= 10
        actual_area = area(high_rad)
    return binary_search(low_rad, high_rad,
                        area, wanted_area,
                        allowed_error, max_reps)

#######################################################################
#######################################################################

def is_in_circle(point_x, point_y, center_x, center_y, rad):
    """
    Checks if a point is in a circle
    """
    diff_x = abs(point_x-center_x)
    diff_y = abs(point_y-center_y)
    distance = math.sqrt(diff_x**2 + diff_y**2)
    return distance <= rad

#######################################################################
#######################################################################

def binary_search(low, high, ordering_function, expected,
                  max_error_ratio=.3, max_reps=1e4):
    """
    Binary search algorithm.
    low is a point where ordering function is under expected
    high is a point where ordering function is over expected
    ordering function must be a function(raise TypeError)
    """
    max_error = expected*max_error_ratio
    error = max_error+1
    rep = 0
    while (error > max_error) or (rep < max_reps):
        rep += 1
        mid = float(low+high)/2
        try:
            actual_value = ordering_function(mid)
        except TypeError as exception:
            msg = "binary search: Ordering_function must be a function. "
            msg += "Introduced: "+str(ordering_function)
            print msg
            raise exception
        error = abs(actual_value-expected)
        if actual_value > expected:
            high = mid
        elif actual_value < expected:
            low = mid
        else:
            return mid
    return mid

#######################################################################
#######################################################################

def binary_order(array, ordering_id):
    """
    Orders an array using ordering_function.
    ordering_id must return an integer
    orders from low id to high id
    """
    #base part
    if len(array) <= 1:
        return array
    #recursive part (divide array in 2 parts and order those parts)
    length = len(array)
    array_1 = array[:length/2]
    array_2 = array[length/2:]
    array_1 = binary_order(array_1, ordering_id)
    array_2 = binary_order(array_2, ordering_id)
    ordered_array = []
    #merge part (take ordered arrays and ordered them into a bigger one)
    element_1 = array_1.pop(0)
    element_2 = array_2.pop(0)
    unemptied_array = None
    while True:
        id_1 = ordering_id(element_1)
        id_2 = ordering_id(element_2)
        if id_1 < id_2:
            ordered_array.append(element_1)
            try:
                element_1 = array_1.pop(0)
            except IndexError:
                ordered_array.append(element_2)
                unemptied_array = array_2
                break
        else:   #If id_1==id_2 order isn't important
            ordered_array.append(element_2)
            try:
                element_2 = array_2.pop(0)
            except IndexError:
                ordered_array.append(element_1)
                unemptied_array = array_1
                break
    for element in unemptied_array:
        ordered_array.append(element)
    return ordered_array

#######################################################################
#######################################################################

def erase_length_one_elements(group, minimum_length=2):
    """
    Erase elements of groups of lenght < minimum_length
    """
    new_group = []
    try:
        while True:
            element = group.pop()
            if len(element) >= minimum_length:
                new_group.append(element)
    except IndexError:
        pass
    return new_group

#######################################################################
#######################################################################

def import_files(folder="./", regular_expression=r'rods_[0-9]*'):
    """
    Import all files using glob and checking with reg exp.
    """
    if not re.match(r"\.?/?[a-zA-Z0-9\.]/*", folder):
        print "You must provide a folder like: \"./\" or \"/home/user/\""
        raise ValueError
    names = get_file_names(folder=folder, regular_expression=regular_expression)
    files = []
    for name in names:
        files.append(open(name, 'r'))
    return names, files


def get_file_names(folder="./", regular_expression=r'rods_[0-9]*'):
    """
    Get file names of the folder whose names pass regular expression
    """
    names = []
    reg1 = re.compile(regular_expression)
    extension = re.compile(r'.*\.png')
    files = os.listdir(folder)
    for _file in files:
        if reg1.match(_file) and not extension.match(_file):
            names.append(_file)
    return binary_order(names, get_number_from_string)

#######################################################################
#######################################################################

def get_number_from_string(name):
    """
    Gets the number in the name of the file.
    """
    reg = re.compile(r'\d+')
    found = reg.findall(name)
    output = ""
    for found_element in found:
        output += found_element
    return int(output)

#######################################################################
#######################################################################

def import_data(_file, split_char='\t', regular_expression=r'[0-9]\.?[0-9]*'):
    """
    Import data of a file
    Returns an array with data
    """
    if str(type(_file)) != "<type 'file'>":
        print "Passed file argument is not a file descriptor"
        raise ValueError
    reg_exp = re.compile(regular_expression)
    data = []
    try:
        for line in _file:
            #Removes new line char from the line
            dataline = re.split(split_char, line.rstrip('\n'))
            for element in dataline:
                if not reg_exp.match(element):
                    try:
                        dataline.remove(element)
                    except AttributeError:
                        print "Something went badly" + str(dataline)
            data.append(dataline)
    except ValueError:
        print "File must be passed opened to import_data."
    except TypeError:
        print "Perhaps the file passed is not in the right format."
        print _file
    except IndexError:
        print "Error importing files (empty list)"
        raise IndexError
    return data

#######################################################################
#######################################################################

def create_rods(folder="./", kappas=10, allowed_kappa_error=.3,
                radius_correction_ratio=0.1):
    """
    Create one rod for each rod_data and for each file
    returns [RodGroup1, RodGroup2, ...]
    """
    names, files = import_files(folder=folder)
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
    run_processes(processes)        #This part seem to take a lot of time.
    try:
        while True:
            [index, state] = states_queue.get(False)
            states[index] = state
    except Queue.Empty:
        pass    
    return names, states

#######################################################################
#######################################################################

def create_rods_process(kappas, allowed_kappa_error,
                        radius_correction_ratio, names,
                        files, index, states_queue):
    """
    Process of method.
    """
    state = SystemState(kappas, allowed_kappa_error,
               radius_correction_ratio, names[index])
    data = import_data(files[index])
    for dataline in data:
        parameters = tuple(dataline)
        new_rod = Rod(parameters)
        state.put_rod(new_rod)
    state.check_rods()
    states_queue.put([index, state])

#######################################################################
#######################################################################

def run_processes(processes, time_out=None):
    """
        Runs all processes using all cores.
    """
    running = []
    cpus = mp.cpu_count()
    try:
        #while True:
        for cpu in range(cpus):
            next_process = processes.pop()
            running.append(next_process)
            next_process.start()
    except IndexError:
        pass
    try:
        while True:
            process = running.pop()
            if process.is_alive():
                running.append(process)
            else:
                next_process = processes.pop()
                running.append(next_process)
                next_process.start()
    except IndexError:
        pass
    if not time_out:
        while True:
            for process in running:
                if not process.is_alive():
                    running.remove(process)
            if not len(running):
                break
    else:
        for process in running:
            process.join(time_out)

