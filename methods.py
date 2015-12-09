"""
Methods library.
"""
import re, math, os
import multiprocessing as mp
from PIL import Image

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
    assert 0 <= phi <= math.pi/2, "Error in angle"
    if phi <= 1e-10:
        if min_dist > 0:
            return 0
        else:
            return math.pi*rad**2
    section = phi*rad**2
    assert section >= 0, "segment_area: Section is negative"
    distance_between_intersections = 2*min_dist*math.tan(phi)
    msg = "segment_area: distance between "
    msg += "intersections can't be greater than diameter"
    assert distance_between_intersections <= 2*rad, msg
    triangle_area = distance_between_intersections*min_dist/2.0
    msg = "segment_area: Triangle area must be smaller than section area"
    msg += "\nRatio="+str(triangle_area*1.0/section)
    assert triangle_area < section, msg
    if min_dist >= 0:
        output = section - triangle_area
    else:
        output = math.pi*rad**2 - section + triangle_area
    msg = "segment_area: Obtained area is negative. "
    msg += "Values: rad:"+str(rad)
    msg += " min_dist:"+str(min_dist)+" rat:"+str(min_dist/rad)
    msg += " phi:"+str(phi)+" area:"+str(output)
    assert output > 0, msg
    return output




def effective_area(small_rad, small_position_rad, main_rad):
    """
    Computes the area of the small circle intersected with main circle.
    There are errors in this method implementation.
    """
    # circle completely included in the bigger one
    if small_rad+small_position_rad <= main_rad:
        return math.pi*small_rad**2
    # assert small_position_rad <= main_rad, "Circle is outside the bigger one"
    min_dist = compute_min_dist(small_rad, small_position_rad, main_rad)
    if min_dist >= small_rad:
        return math.pi*small_rad**2
    elif abs(min_dist) >= small_rad:
        return 0
    min_dist_main = small_position_rad+min_dist
    correction = segment_area(main_rad, min_dist_main)
    msg = "effective_area: Correction must be smaller than small circle's area"
    assert correction < math.pi*small_rad**2, msg
    section_area = segment_area(small_rad, min_dist)
    small_area = math.pi*small_rad**2
    msg = "In the limit, h=-rad has to return total area"
    assert small_area == segment_area(small_rad, -small_rad), msg
    msg = "Correction too high: Ration: "+str(float(correction)/small_area)
    assert correction < small_area, msg
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
    try:
        min_dist /= (2*small_position_rad)
    except:
        min_dist = main_rad*2
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



def distance_between_points(point1, point2):
    """
    Distance between points.
    """
    diff_x = abs(point1[0]-point2[0])
    diff_y = abs(point1[1]-point2[1])
    distance = math.sqrt(diff_x**2 + diff_y**2)
    return distance



def is_in_circle(point_x, point_y, center_x, center_y, radius):
    """
    Checks if a point is in a circle
    """
    point = [point_x, point_y]
    center = [center_x, center_y]
    return distance_between_points(point, center) <= radius




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
        else:
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


def get_file_names(folder="./", regular_expression=r'^rods_[0-9]{4}$'):
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


def imagej():
    """
    Creates and run imagej script to get data of images.
    """
    image_names = get_file_names(regular_expression=r'IMG_[0-9]{4}.JPG')
    numbers = [get_number_from_string(image) for image in image_names]
    try:
        start = min(numbers)
        end = max(numbers)
        imagej_template = open("./acquisition_v4.ijm", 'r')
        imagej_script = open("./imagej_script.ijm", 'w')
        current_folder = os.getcwd()
        for line in imagej_template:
            output_line = line
            if re.match(r'.*for.img.*', line):
                output_line = "for(img_num="
                output_line += str(start) + "; img_num<="
                output_line += str(end)+ "; img_num++){\n"
            if re.match(r'.*\"\.\".*', line):
                output_line = "\tfolder = \""+str(current_folder)+"\";"
            imagej_script.write(output_line)

        imagej_template.close()
        imagej_script.close()

        run_imagej = raw_input("Run imagej script?(Yn)")
        if run_imagej == "y" or not run_imagej:
            os.system("imagej ./imagej_script.ijm")
    except ValueError:
        pass



def get_image_dates(folder="./"):
    """
    Returns a dictionary:
    {image_number: date}
    """
    image_names = get_file_names(regular_expression=r'IMG_[0-9]{4}.JPG')
    numbers = [get_number_from_string(image) for image in image_names]
    dates = {}
    for index in range(len(image_names)):
        name = image_names[index]
        number = numbers[index]
        date = get_date_taken(name, folder=folder)
        dates[number] = date
    return dates

def export_image_dates(file_name="dates.txt", folder="./"):
    """
    Saves all dates in a file.
    Creates a file with dates of images:
    image_name      date
    """
    file_path = str(folder) + str(file_name)
    output_file = open(file_path, 'w')
    dates = get_image_dates(folder)
    for image_number in dates.keys():
        line = str(image_number) + "\t"
        date = dates[image_number]
        line += str(date) + "\n"
        output_file.write(line)
    output_file.close()

def import_image_dates(file_name="dates.txt", folder="./"):
    """
    Get dates from a file.
    """
    file_path = str(folder) + str(file_name)
    output_file = open(file_path, 'r')
    dates = {}
    for line in output_file:
        data = line.strip('\n')
        data = data.split('\t')
        dates[int(data[0])] = data[1]
    return dates

def get_date_taken(img_name, folder="./"):
    """
    Returns date info about a image.
    """
    image = str(folder) + str(img_name)
    return Image.open(image)._getexif()[36867]

def are_in_burst(dates, image_num_1, image_num_2):
    """
    Return True if image 1 and image 2 are in a burst.
    Image 2 must be the next to image 1
    """
    image_num_1 = int(image_num_1)
    image_num_2 = int(image_num_2)
    if image_num_2 != image_num_1+1:
        msg = "Image 2 must be the next to image 1."
        raise ValueError(msg)
    date1 = dates[image_num_1].split(' ')
    date2 = dates[image_num_2].split(' ')
    if date1[0] != date2[0]:
        return False
    time1 = date1[1].split(':')
    time2 = date2[1].split(':')
    time1 = [int(time_part) for time_part in time1]
    time2 = [int(time_part) for time_part in time2]
    for index in [2, 1]:
        part = time1[index]
        if int(part) >= 59:
            part = part-60
            time1[index-1] += 1
    if time2[2] != time1[2]+1 and time2[2] != time1[2]:
        return False
    if time2[1] != time1[1]:
        return False
    return True

def run_processes(processes, cpus=None):
    """
        Runs all processes using all cores.
    """
    running = []
    if not cpus:
        cpus = 4#int(mp.cpu_count()/2)
    try:
        for dummy in range(cpus):
            next_process = processes.pop()
            running.append(next_process)
            next_process.start()
    except IndexError:
        pass
    return running, processes

def array_average(array_of_arrays):
    """
    Gets average array value over a list of arrays.
    """
    number_of_arrays = len(array_of_arrays)
    array_length = len(array_of_arrays[0])
    output = [0 for dummy in range(array_length)]
    for array in array_of_arrays:
        if len(array) != array_length:
            print len(array), array_length
            msg = "Arrays are not of same length."
            raise ValueError(msg)
        for index in range(len(array)):
            output[index] += float(array[index])/number_of_arrays
    return output

def vector_module(vector):
    """
       Returns array's module. 
    """
    return math.sqrt(vector[0]**2+vector[1]**2)

def vector_angle(vector):
    """
        Returns array's vector.
    """
    try:
        angle = math.atan(vector[1]/float(vector[0]))
    except ZeroDivisionError:
        angle = math.pi
    return angle

