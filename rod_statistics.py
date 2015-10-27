"""
Analyze rods data
"""
from queue import Queue
import re
import multiprocessing as mp    #for using all cores
import math
import os
import matrix



#consts
RADIUS = 1704*0.5
CENTER_X = 519+RADIUS
CENTER_Y = 96+RADIUS
#ROD_LENGTH,ROD_DIAMETER #Now there are two types of rod
# Time consts
FRAME_INTERVAL = 1.0/3



class Rod(object):
    """
    Rod object.
    """

    def __init__(self, (ID, area, xm, ym, major, minor,
                        angle, feret, feretx, ferety,
                        feretangle, minferet, xstart, ystart)):
        """
        Initialization of rod
        """                                    #Column
        self._id = int(ID)                     #0
        self.area = float(area)                #1
        self.x_mid = float(xm)                 #2
        self.y_mid = float(ym)                 #3
        self.major = float(major)              #4
        self.minor = float(minor)              #5
        self.angle = float(angle)              #6
        self.feret = float(feret)              #7
        self.feret_x = float(feretx)           #8
        self.feret_y = float(ferety)           #9
        self.feret_angle = float(feretangle)   #10
        self.min_feret = float(minferet)       #11
        self.x_start = float(xstart)           #12
        self.y_start = float(ystart)           #13
        dif_x = abs(self.x_mid - CENTER_X)
        dif_y = abs(self.y_mid - CENTER_Y)
        self.distance_to_center = math.sqrt(dif_x**2+dif_y**2)
        self.kappa = float(self.feret)/self.min_feret

    def is_in_circle(self, center, rad):
        """
        Checks if rod is in circle.
        """
        return is_in_circle(self.x_mid, self.y_mid,
                            center[0], center[1], rad)

    def is_in_main(self, allowed_distance_from_border=0):
        """
        Checks if rod is in main.
        It deletes all rods that are near the border (delta_rad).
        """
        return self.is_in_circle((CENTER_X, CENTER_Y),
                            RADIUS-allowed_distance_from_border)

    def has_valid_proportions(self, kappas, allowed_error):
        """
        Checks if rod has valid proportions.
        """
        passed = []
        try:
            for kappa in kappas:
                condition = abs(self.kappa-kappa) < allowed_error
                passed.append(condition)
        except TypeError:
            condition = abs(self.kappa-kappas) < allowed_error
            passed.append(condition)
        output = False
        for condition in passed:
            output = output or condition
        return output

    def is_valid_rod(self, kappas,
                    allowed_kappa_error,
                    allowed_distance_from_border):
        """
        Checks if this is a rod looking at different factors
        If it is a group of two rods that are very near.
        Remove rods that are near the border.
        """
        is_in_main = self.is_in_main(allowed_distance_from_border)
        has_valid_proportions = self.has_valid_proportions(kappas,
                                                           allowed_kappa_error)
        output = is_in_main and has_valid_proportions
        return output






class SystemState(object):
    """
    Group of rods in a moment.
    Each image has to be translated into a RodGroup (by this class?)
    """
    def __init__(self, kappas=10, allowed_kappa_error=.5,
                allowed_distance_from_border=0,
                id_string=""):
        """
        Initialization
        """
        self._rods = Queue()
        self._number_of_particles = 0
        self._kappas = kappas
        self._allowed_kappa_error = allowed_kappa_error
        self._allowed_distance_from_border = allowed_distance_from_border
        self._rad_for_division = -1
        self._actual_subdivision = []
        self._changed = False
        self._density_matrix = []
        self.id_string = id_string

    def add_rod(self, rod):
        """
        Adds a rod to the group
        """
        if rod.is_valid_rod(self._kappas, self._allowed_kappa_error,
                            self._allowed_distance_from_border):
            self._rods.join(rod)
            self._number_of_particles += 1
            self._changed = True

    def get_rod(self):
        """
        Returns the first rod in the queue
        The rod is removed of the group!
        """
        self._number_of_particles -= 1
        self._changed = True
        return self._rods.get_next()

    def remove_rod(self, rod):
        """
        Removes a rod from the group (queue object mod needed)
        """
        self._rods.delete(rod)
        self._number_of_particles -= 1
        self._changed = True

    def _divide_in_circles(self, rad):
        """
        Divides rods into groups contained in circles.
        """
        if self._rad_for_division == rad and not self._changed:
            return
        if (rad < 0) or (rad > RADIUS):
            print "Use a correct radius (0<rad<main_rad)"
            raise ValueError
        diff = int(rad/2)
        start_x = CENTER_X-RADIUS
        end_x = CENTER_X+RADIUS
        start_y = CENTER_Y-RADIUS
        end_y = CENTER_Y+RADIUS
        subsystems = []
        max_times = int((end_x-start_x)/diff+1)
        possible_x_values = [start_x + times*diff
                             for times in range(max_times)]
        max_times = int((end_y-start_y)/diff+1)
        possible_y_values = [start_y + times*diff
                             for times in range(max_times)]
        for actual_x in possible_x_values:
            for actual_y in possible_y_values:
                if is_in_circle(actual_x, actual_y,
                                CENTER_X, CENTER_Y,
                                RADIUS):
                    diff_x = abs(actual_x-CENTER_X)
                    diff_y = abs(actual_y-CENTER_Y)
                    pos_rad = math.sqrt(diff_x**2+diff_y**2)
                    same_area_radius = same_area_rad(rad, pos_rad, RADIUS)
                    subsystem = SubsystemState((actual_x, actual_y),
                                               same_area_radius)
                    subsystem.add_rods(list(self._rods))
                    subsystems.append(subsystem)
        self._actual_subdivision = subsystems

    def compute_density_matrix(self, rad=100, normalized=False,
                               divided_by_area=False):
        """
        Computes density of the system.
        Returns a plotable matrix
        """
        if self._rad_for_division == rad and not self._changed:
            return self._density_matrix
        self._divide_in_circles(rad)
        density = []
        for subsystem in self._actual_subdivision:
            subdensity = [subsystem.center[0], subsystem.center[1]]
            dens = subsystem.get_density()
            if divided_by_area:
                dens /= subsystem.area
            if normalized:
                dens /= self._number_of_particles
            subdensity.append(dens)
            density.append(subdensity)
        self._density_matrix = density
        return self._density_matrix

    def density_matrix_for_plot(self):
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

    #def 



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
        self._density = -1
        self.center = center
        self.rad = rad
        self._area = float(math.pi*rad**2)

    @property
    def area(self):
        """
        getter for area
        """
        return self._area

    def _update_density(self):
        """
        Overrides.
        Computes density of the group.
        """
        self._density = self._number_of_particles


    def get_density(self):
        """
        Returns the density of the group
        Perhaps sometimes update_density method can be ignored
        """
        self._update_density()
        return self._density

    def add_rods(self, rod_list):
        """
        Add all rods of the list that are inside the circle.
        """
        try:
            for rod in rod_list:
                if rod.is_in_circle(self.center, self.rad):
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
                allowed_distance_from_border=0):
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
                   allowed_distance_from_border=allowed_distance_from_border,
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
        states.append(state)
    return names, states





def segment_area(rad, min_dist):
    """
    Computes the area of an intersection of a circle with a line
    (the part that doesn't have the center)
    rad: radius of the circle
    min_dist: minimum distance from small circle center to a line that joins
        both intersections of the circles.
    """
    min_dist = float(min_dist)
    if rad < abs(min_dist):
        message = "In segment_area:\n\th can't be greater than "
        message += "rad\nvalues:\trad="+str(rad)+"\n\th="+str(min_dist)
        raise ValueError(message)
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
    #DEBUG#
    msg = "Error in min_dist computation: [small_rad,min_dist] "
    msg += str([small_rad, abs(min_dist)])
    assert small_rad > abs(min_dist), msg
    #######
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





def compute_min_dist(small_rad, small_position_rad, main_rad):
    """
    Computes the distance from small circle center to the line that joins both
    circles' intersections.
    """
    min_dist = float((main_rad**2)-(small_position_rad**2)-(small_rad**2))
    min_dist /= (2*small_position_rad)
    return min_dist






def same_area_rad(small_rad, small_position_rad,
                    main_rad, allowed_error_ratio=.05,
                    max_reps=1e2):
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
        return effective_area(rad, small_position_rad, main_rad)
    actual_area = area(high_rad)
    while actual_area < wanted_area:
        high_rad *= 10
        actual_area = area(high_rad)
    return binary_search(low_rad, high_rad,
                        area, wanted_area,
                        allowed_error, max_reps)





def is_in_circle(point_x, point_y, center_x, center_y, rad):
    """
    Checks if a point is in a circle
    """
    diff_x = abs(point_x-center_x)
    diff_y = abs(point_y-center_y)
    distance = math.sqrt(diff_x**2 + diff_y**2)
    return distance <= rad





def binary_search(low, high, ordering_function, expected, max_error, max_reps):
    """
    Binary search algorithm.
    low is a point where ordering function is under expected
    high is a point where ordering function is over expected
    ordering function must be a function(raise TypeError)
    """
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





def binary_order(array, ordering_id):
    """
    Orders an array using ordering_function.
    ordering_id must return an integer
    orders from low id to high id
    """
    #base part
    if len(array) <= 1:
        return array
    #recursive part
    length = len(array)
    array_1 = array[:length/2]
    array_2 = array[length/2:]    
    array_1 = binary_order(array_1, ordering_id)
    array_2 = binary_order(array_2, ordering_id)
    ordered_array = []
    #merge part
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

