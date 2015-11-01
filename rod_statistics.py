"""
Analyze rods data
"""
from queue import Queue
import re
#import multiprocessing as mp    #for using all cores
import math
import os
import matrix
from methods import *


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

    @property
    def id(self):
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
        output += str(self.id)
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
        angle1 = (self.angle-rod.angle)*math.pi/180
        angle2 = math.pi - (self.angle-rod.angle)*math.pi/180
        return min(angle1, angle2)












class SystemState(object):
    """
    Group of rods in a moment.
    Each image has to be translated into a RodGroup (by this class?)
    """
    def __init__(self, kappas=10, allowed_kappa_error=.5,
                radius_correction_ratio=0,
                id_string="", zone_coords=[]):
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
        self._g2 = None
        self._g4 = None
        self._g2_subsystems = []
        self._g4_subsystems = []
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
        except:
            self._fixed_center_radius = False

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
        self._g2 = None
        self._g4 = None
        self._g2_subsystems = []
        self._g4_subsystems = []
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

    def compute_center_and_radius(self):
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
        self.compute_center_and_radius()
        return self._center_x, self._center_y

    @property
    def center_x(self):
        """
        Returns center of the system. (x)
        """
        self.compute_center_and_radius()
        return self._center_x

    @property
    def center_y(self):
        """
        Returns center of the system. (y)
        """
        self.compute_center_and_radius()
        return self._center_y

    @property
    def radius(self):
        """
        Returns radius of the system.
        """
        self.compute_center_and_radius()
        return self._radius

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
        # Defining zone and distance between points.
        diff = int(rad/2)
        start_x = self.center_x-self.radius
        end_x = self.center_x+self.radius
        start_y = self.center_y-self.radius
        end_y = self.center_y+self.radius
        # Getting all possible x and y values.
        subsystems = []
        max_times = int((end_x-start_x)/diff+1)
        possible_x_values = [start_x + times*diff
                             for times in range(max_times)]
        max_times = int((end_y-start_y)/diff+1)
        possible_y_values = [start_y + times*diff
                             for times in range(max_times)]
        # Main loops.
        for actual_x in possible_x_values:
            for actual_y in possible_y_values:
                if is_in_circle(actual_x, actual_y,
                                self.center[0], self.center[1],
                                self.radius):
                    diff_x = abs(actual_x-self.center[0])
                    diff_y = abs(actual_y-self.center[1])
                    pos_rad = math.sqrt(diff_x**2+diff_y**2)
                    same_area_radius = same_area_rad(rad, pos_rad, self.radius)
                    subsystem = SubsystemState((actual_x, actual_y),
                                               same_area_radius)
                    subsystem.add_rods(list(self._rods))
                    subsystems.append(subsystem)
        self._actual_subdivision = subsystems

    def compute_density_matrix(self, rad=100, normalized=False,
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

    @property
    def plottable_density_matrix(self):
        """
        Returns 3 arrays: one for x, another for y and the last for values.
        """
        if len(self._density_matrix) == 0:
            self.compute_density_matrix()
        x_values = []
        y_values = []
        z_values = []
        for row in self._density_matrix:
            x_values.append(row[0])
            y_values.append(row[1])
            z_values.append(row[2])
        return x_values, y_values, z_values

    def compute_g2_and_g4(self):
        """
        Computes g2 and g4 values
        """
        cos2_av, sin2_av, cos4_av, sin4_av = 0, 0, 0, 0
        for rod in list(self._rods):
            angle = rod.angle*math.pi/180.0
            cos2_av += math.cos(2*angle)/self.number_of_rods
            sin2_av += math.sin(2*angle)/self.number_of_rods
            cos4_av += math.cos(4*angle)/self.number_of_rods
            sin4_av += math.sin(4*angle)/self.number_of_rods
        self._g2 = math.sqrt(cos2_av**2+sin2_av**2)
        self._g4 = math.sqrt(cos4_av**2+sin2_av**2)

    @property
    def g2(self):
        """
        sqrt(<cos(2*angle)>^2+<sin(2*angle)>^2)
        """
        if not self._g2:
            self.compute_g2_and_g4()
        return self._g2

    @property
    def g4(self):
        """
        sqrt(<cos(4*angle)>^2+<sin(4*angle)>^2)
        """
        if not self._g4:
            self.compute_g2_and_g4()
        return self._g4

    def compute_g2_g4_matrices(self, rad):
        """
        Computes g2 and g4 matrices for subgroups.
        """ 
        self._divide_in_circles(rad)
        for subsystem in self._actual_subdivision:
            g2 = [subsystem.center[0], subsystem.center[1]]
            g4 = [subsystem.center[0], subsystem.center[1]]
            g2.append(subsystem.g2)
            g4.append(subsystem.g4)
            self._g2_subsystems.append(g2)
            self._g4_subsystems.append(g4)

    @property
    def g2_plot_matrix(self):
        """
        Returns values for plotting g2 matrix.
        """
        xval = []
        yval = []
        zval = []
        for subsystem in self._g2_subsystems:
            xval.append(subsystem[0])
            yval.append(subsystem[1])
            zval.append(subsystem[2])
        return xval, yval, zval

    @property
    def g4_plot_matrix(self):
        """
        Returns values for plotting g2 matrix.
        """
        xval = []
        yval = []
        zval = []
        for subsystem in self._g4_subsystems:
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
            if self.g2 > 0.3 and self.g4 < 0.3:
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
    def compute_g2_g4_matrices(self, rad):
        """
        Computes g2 and g4 matrices for subgroups.
        """ 
        self._divide_in_circles(rad)
        len_g2 = len(self._g2_subsystems)
        len_g4 = len(self._g4_subsystems)
        if  len_g2 == 0 or len_g4 == 0:
            for subsystem in self._actual_subdivision:
                g2 = [subsystem.center[0], subsystem.center[1]]
                g4 = [subsystem.center[0], subsystem.center[1]]
                g2.append(subsystem.g2)
                g4.append(subsystem.g4)
                self._g2_subsystems.append(g2)
                self._g4_subsystems.append(g4)

    def compute_average_angle_matrix(self, rad):
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

    @property
    def plottable_average_angle_matrix(self):
        """
        Returns a plottable average angle matrix.
        """
        xval = []
        yval = []
        zval = []
        for subsystem in self._angle_matrix:
            xval.append(subsystem[0])
            yval.append(subsystem[1])
            zval.append(subsystem[2])
        return xval, yval, zval

    def compute_all_matrices(self, rad):
        """
        Computes all possible matrices.
        O(N^2)
        """
        self.compute_average_angle_matrix(rad)
        self.compute_g2_g4_matrices(rad)
        self.compute_density_matrix(rad=rad)
        self.compute_closest_rod_matrix()
        self.compute_relative_g2_g4_matrices(rad)

    @property
    def angle_histogram(self):
        """
        Returns all angles in a list to make an histogram.
        """
        output = [rod.angle for rod in list(self._rods)]
        return output

    def get_closest_rod(self, rod):
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

    def compute_closest_rod_matrix(self):
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
            new_row.append(self.get_closest_rod(rod))
            closest_rod_matrix.append(new_row)
        self._closest_rod_matrix = closest_rod_matrix


    @property
    def closest_rod_matrix(self):
        """
        Returns closest rod matrix.
        See compute_closest_rod_matrix for more info.
        """
        if len(self._closest_rod_matrix) == 0:
            self.compute_closest_rod_matrix()
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
        
    def compute_relative_g2_g4_matrices(self, rad):
        """
        Computes g2 and g4 matrices for subgroups.
        """ 
        self._divide_in_circles(rad)
        len_g2 = len(self._relative_g2_subsystems)
        len_g4 = len(self._relative_g4_subsystems)
        if  len_g2 == 0 or len_g4 == 0:
            for subsystem in self._actual_subdivision:
                g2 = [subsystem.center[0], subsystem.center[1]]
                g4 = [subsystem.center[0], subsystem.center[1]]
                g2.append(subsystem.g2)
                g4.append(subsystem.g4)
                self._relative_g2_subsystems.append(g2)
                self._relative_g4_subsystems.append(g4)

    @property
    def relative_g2_plot_matrix(self):
        """
        Returns values for plotting g2 matrix.
        """
        xval = []
        yval = []
        zval = []
        for subsystem in self._relative_g2_subsystems:
            xval.append(subsystem[0])
            yval.append(subsystem[1])
            zval.append(subsystem[2])
        return xval, yval, zval

    @property
    def relative_g4_plot_matrix(self):
        """
        Returns values for plotting g2 matrix.
        """
        xval = []
        yval = []
        zval = []
        for subsystem in self._relative_g4_subsystems:
            xval.append(subsystem[0])
            yval.append(subsystem[1])
            zval.append(subsystem[2])
        return xval, yval, zval








class SubsystemState(SystemState):
    """
    Group of rods. Used to put all rods that are in a zone or
    have something in common.
    """

    def __init__(self, center, rad):
        """
        Initialization
        """
        SystemState.__init__(self)
        self._center = center
        self._rad = rad
        self._area = float(math.pi*rad**2)

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










def import_files(folder="./", regular_expression='rods_[0-9]*'):
    """
    Import all files using glob and checking with reg exp.
    """
    if not re.match("\.?/?[a-zA-Z0-9\.]/*", folder):
        print "You must provide a folder like: \"./\" or \"/home/user/\""
        raise ValueError
    names = get_file_names(folder=folder, regular_expression=regular_expression)
    files = []
    for name in names:
        files.append(open(name, 'r'))
    return names, files












def get_file_names(folder="./", regular_expression='rods_[0-9]*'):
    """
    Get file names of the folder whose names pass regular expression
    """
    names = []
    reg1 = re.compile(regular_expression)
    extension = re.compile('.*\.png')
    files = os.listdir(folder)
    for _file in files:
        if reg1.match(_file) and not extension.match(_file):
            names.append(_file)
    return binary_order(names, file_name_ordering_function)












def file_name_ordering_function(name):
    """
    Gets the number in the name of the file.
    """
    reg = re.compile('\d+')
    found = reg.findall(name)
    output = ""
    for found_element in found:
        output += found_element
    return int(output)











def import_data(_file, split_char='\t', regular_expression='[0-9]\.?[0-9]*'):
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
    states = []
    for index in range(len(files)):
        _file = files[index]
        name = names[index]
        state = SystemState(kappas=kappas,
                   allowed_kappa_error=allowed_kappa_error,
                   radius_correction_ratio=radius_correction_ratio,
                   id_string=name)
        data = import_data(_file)
        for dataline in data:
            parameters = tuple(dataline)
            try:
                new_rod = Rod(parameters)
            except ValueError as exception:
                print parameters
                print exception.message
                raise ValueError
            state.add_rod(new_rod)
        state.compute_center_and_radius()
        state.check_rods()
        states.append(state)
    return names, states

