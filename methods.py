"""
    Needed methods for study rod systems.
"""

import multiprocessing as mp    #for using all cores
import math
import re
import os

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














def is_in_circle(point_x, point_y, center_x, center_y, rad):
    """
    Checks if a point is in a circle
    """
    diff_x = abs(point_x-center_x)
    diff_y = abs(point_y-center_y)
    distance = math.sqrt(diff_x**2 + diff_y**2)
    return distance <= rad











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
                            radius_correction_ratio, names, index,
                            states_queue))
        processes.append(process)
    run_processes(processes)
    try:
        [index, state] = states_queue.get(False)
        states[index] = state
    except Queue.Empty:
        pass    
    return names, states

def create_rods_process(kappas, allowed_kappa_error, radius_correction_ratio, names, index):
    """
    Process of method.
    """
    state = SystemState(kappas, allowed_kappa_error,
               radius_correction_ratio, names[index])
    data = import_data(files[index])
    for dataline in data:
        parameters = tuple(dataline)
        new_rod = Rod(parameters)
        state.add_rod(new_rod)
    state.check_rods()
    states_queue.put([index, state])

def run_processes(processes):
    """
        Runs all processes using all cores.
    """
    running = []
    try:
        while True:
            next_process = processes.pop()
            running.append(next_process)
            next_process.start()
    except IndexError:
        pass
    for process in running:
        process.join(30)

