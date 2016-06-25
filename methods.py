"""
Methods library.
"""
import re, math, os
import multiprocessing as mp
from PIL import Image
import cPickle
import settings
import numpy
from matplotlib import animation
from scipy.interpolate import griddata
#from matplotlib.mlab import griddata
import datetime
import inspect
import lzo
import lz4
import zlib
import gzip

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'
CLEAR_LAST = CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE
if settings.special_chars:
    WHITE_BLOCK = u'\u25A0'
else:
    WHITE_BLOCK = 'X'
GZIP = 4
LZ4_FAST = 3
LZ4 = 2
LZO = 1
ZLIB = 0

using_cl = False
_writer = animation.writers['ffmpeg']
writer = _writer(fps=15, metadata=dict(artist='Jaime Perez Aparicio'), bitrate=1800)
WRITER = writer

try:
    import pyopencl as cl
    platform = cl.get_platforms()
    my_gpu_devices = platform[0].get_devices(device_type=cl.device_type.CPU)
    ctx = cl.Context(devices=my_gpu_devices)
    using_cl = True
except:
    using_cl = False

using_cl = False

def change_compression_coef(new_coef):
    """
    Changes compression algorithm coeficient.
    """
    os.system("del settings")
    settings = open("settings.py", "w")
    settings.write("default_comp_level = {0}".format(9))
    settings.close()
    import settings


def effective_area(small_rad, centers_distance, main_rad):
    """
    Computes the area of the small circle intersected with main circle.
    There are errors in this method implementation.
    """
    # circle completely included in the bigger one
    if small_rad+centers_distance <= main_rad:
        return math.pi*small_rad**2
    R, r, d = main_rad, small_rad, centers_distance
    if R<r:
        rr = R
        R = r
        r = rr
    try:
        q_1 = float(d**2+r**2-R**2)/(2.0*d*r)
        A_1 = float(r**2)*math.acos(q_1)
        q_2 = float(d**2-r**2+R**2)/(2.0*d*R)
        A_2 = float(R**2)*math.acos(q_2)
        corr = -.5*math.sqrt((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R))
        A = A_1+A_2+corr
        return A
    except ValueError:
        #print small_rad, centers_distance, main_rad
        return float('inf')




def compute_min_dist(small_rad, centers_distance, main_rad):
    """
    Computes the distance from small circle center to the line that joins both
    circles' intersections.
    """
    try:
        min_dist = (main_rad**2)-(centers_distance**2)-(small_rad**2)
        min_dist = float(min_dist)
    except OverflowError:
        return small_rad*1.1
    try:
        min_dist /= (2*centers_distance)
    except:
        min_dist = main_rad*2
    if min_dist > small_rad:
        return small_rad
    if abs(min_dist) > small_rad:
        return -small_rad
    return min_dist




def same_area_rad(small_rad, centers_distance,
                    main_rad, allowed_error_ratio=.2,
                    max_reps=10):
    """
    Computes a new radius. With that, effective area is the same small circle's.
    Better use binary search
    """
    if centers_distance + small_rad <= main_rad:
        return small_rad
    wanted_area = math.pi*small_rad**2
    allowed_error = wanted_area * allowed_error_ratio
    low_rad = small_rad
    high_rad = small_rad*10
    def area(rad, centers_distance=centers_distance, main_rad=main_rad):
        """
        Needed function for binary search, as only 1 arg is allowed.
        """
        try:
            return effective_area(rad, centers_distance, main_rad)
        except OverflowError:
            return float('inf')
    actual_area = area(high_rad)
    reps = 0
    while actual_area < wanted_area:
        high_rad *= 10
        actual_area = area(high_rad)
        if reps >= max_reps:
            break
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
    names = binary_order(names, get_number_from_string)
    return names




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
        print "Any image in folder"



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
    image_name\tdate
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

def are_in_burst(date1, date2):
    """
    Return True if image 1 and image 2 are in a burst.
    Image 2 must be the next to image 1
    """
    date1 = date1.split(' ')
    date2 = date2.split(' ')
    days1 = date1[0].split(':')
    days2 = date1[0].split(':')
    days1 = [int(date_part) for date_part in days1]
    days2 = [int(date_part) for date_part in days2]
    time1 = date1[1].split(':')
    time2 = date2[1].split(':')
    time1 = [int(time_part) for time_part in time1]
    time2 = [int(time_part) for time_part in time2]
    time1 = days1[0]*365*30*24*3600+days1[1]*30*24*3600+days1[2]*24*3600+time1[0]*3600+time1[1]*60+time1[2]
    time2 = days2[0]*365*30*24*3600+days2[1]*30*24*3600+days2[2]*24*3600+time2[0]*3600+time2[1]*60+time2[2]
    if abs(time1-time2)<=2:
        return True
    return False

def are_in_burst_queue(index, date1, date2, output_queue):
    """
    Multiprocessing friendly method.
    """
    output = are_in_burst(date1, date2)
    output_queue.put([index, output])

def time_difference(dates, image_num_1, image_num_2):
    """
        Returns the difference of time between 2 images.
    """
    image_num_1 = int(image_num_1)
    image_num_2 = int(image_num_2)
    date1 = dates[image_num_1].split(' ')
    date2 = dates[image_num_2].split(' ')
    date_1 = date1[0].split(':')
    date_2 = date2[0].split(':')
    for index in range(len(date_1)):
        if date_1[index] != date_2[index]:
            return 15*60
    time1 = date1[1].split(':')
    time2 = date2[1].split(':')
    time1 = [int(time_part) for time_part in time1]
    time2 = [int(time_part) for time_part in time2]
    time1 = time1[0]*3600+time1[1]*60+time1[2]
    time2 = time2[0]*3600+time2[1]*60+time2[2]
    return abs(time2-time1)

def run_processes(processes, cpus=None):
    """
        Runs all processes using all cores.
    """
    running = []
    if not cpus or cpus is None:
        cpus = settings.cpus
        if not cpus:
            cpus = int(mp.cpu_count())
    try:
        for dummy in range(cpus):
            next_process = processes.pop(0)
            try:
                next_process.start()
            except AssertionError:
                pass
            else:
                running.append(next_process)
    except IndexError:
        pass
    return running, processes

def array_average(array_of_arrays):
    """
    Gets average array value over a list of arrays.
    """
    if not using_cl:
        number_of_arrays = len(array_of_arrays)
        if not number_of_arrays:
            raise ValueError("any array to compute average")
        array_length = len(array_of_arrays[0])
        if not array_length:
            return sum(array_of_arrays)/number_of_arrays
        output = [0 for dummy in range(array_length)]
        array_of_arrays_2 = []
        for array in array_of_arrays:
            valid_array = True
            for index in range(len(array)):
                if array[index] is None:
                    number_of_arrays -= 1
                    valid_array = False
            if valid_array:
                array_of_arrays_2.append(array)
        number_of_arrays = len(array_of_arrays_2)
        for array in array_of_arrays_2:
            if len(array) != array_length:
                print len(array), array_length
                msg = "Arrays are not of same length."
                raise ValueError(msg)
            for index in range(len(array)):
                if array[index] is None:
                    output[index] = None
                    break
                else:
                    output[index] += float(array[index])/number_of_arrays
        return output
    else:
        number_of_arrays = len(array_of_arrays)
        if not number_of_arrays:
            return None
        array_length = len(array_of_arrays[0])
        if not array_length:
            return sum(array_of_arrays)/number_of_arrays
        output = array_of_arrays[0]
        for index in range(number_of_arrays):
            output = sum_arrays_with_cl(output, array_of_arrays[index])
        return normalize_opencl(output, number_of_arrays)

def array_average_sqrt(array_of_arrays):
    """
    Computes sqrt after averagin.
    """
    array = array_average(array_of_arrays)
    array_ = []
    for value in array:
        try:
            sqrtvalue = math.sqrt(value)
        except ValueError:
            sqrtvalue = -1000
        array_.append(sqrtvalue)
    return array

def array_average_N(array_of_arrays):
    """
    Average for data of type: [[[val000, val001], [val010, val011], ...],
                               [[val100, val101], [val110, val111], ...],
                               ...
    """
    N = len(array_of_arrays[0][0])
    M = len(array_of_arrays)
    L = len(array_of_arrays[0])
    array_of_arrays_ = [[[] for dummy_ in range(M)] for dummy in range(N)]
    output = [0 for dummy in range(L)]
    for index in range(len(array_of_arrays)):
        array_ = array_of_arrays[index]
        for col in range(N):
            for index_2 in range(len(array_)):
                array_of_arrays_[col][index].append(array_[index_2][col])
    output_queue = mp.Queue()
    processes = []
    for col in range(N):
        array_of_arrays__ = array_of_arrays_[col]
        process = mp.Process(target=array_average_N_process,
                             args=(col, output_queue, array_of_arrays__))
        processes.append(process)
    num_processes = len(processes)
    running, processes_left = run_processes(processes)
    finished = 0
    while finished < num_processes:
        finished += 1
        [col, array_of_arrays___] = output_queue.get()
        array_of_arrays_[col] = array_of_arrays___
        if len(processes_left):
            new_process = processes_left.pop(0)
            #time.sleep(settings.waiting_time)
            new_process.start()
    return array_of_arrays_

def array_average_N_process(col, output_queue, array_of_arrays):
    """
    Process
    """
    output_queue.put([col, array_average(array_of_arrays)])

def sum_arrays_with_cl(array1, array2):
    """
        Sums 2 arrays with GPU.
    """
    mf = cl.mem_flags
    a_array = numpy.array(array1).astype(numpy.float32)
    b_array = numpy.array(array2).astype(numpy.float32)
    a_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=a_array)
    b_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=b_array)
    dest_buf = cl.Buffer(ctx, mf.WRITE_ONLY, b_array.nbytes)
    prg = cl.Program(ctx, """
    __kernel void sum(__global const float *a,
    __global const float *b, __global float *c)
    {
      int gid = get_global_id(0);
      c[gid] = a[gid] + b[gid];
    }
    """).build()
    prg.sum(queue, a_array.shape, None, a_buf, b_buf, dest_buf)
    a_plus_b = numpy.empty_like(a_array)
    cl.enqueue_copy(queue_cl, a_plus_b, dest_buf)
    return list(a_plus_b)

def array_x_scalar_cl(array, value):
    """
        Multiplies an array by a scalar.
    """
    mf = cl.mem_flags
    np_array = numpy.array(array).astype(numpy.float32)
    scalar = numpy.float32(value)
    a_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np_array)
    dest_buf = cl.Buffer(ctx, mf.WRITE_ONLY, np_array.nbytes)
    prg = cl.Program(ctx, """
    __kernel void norm(__global const float *a,
                      const float b, __global float *c)
    {
      int gid = get_global_id(0);
      c[gid] = a[gid]*b;
    }
    """).build()
    prg.norm(queue, np_array.shape, None, a_buf, scalar, dest_buf)
    output = numpy.empty_like(np_array)
    cl.enqueue_copy(queue_cl, output, dest_buf)
    return list(output)

def vector_distance(vector1, vector2):
    """
        Returns distance between vectors
    """
    diff = []
    assert len(vector1)==len(vector2), "Use vectors of same dimension\n\n"
    for index in range(len(vector1)):
        diffi = vector1[index]-vector2[index]
        diff.append(diffi)
    return vector_module(diff)

def vector_module(vector):
    """
       Returns array's module.
    """
    output = 0
    for value in vector:
        output += value**2
    return math.sqrt(output)

def vector_angle(vector):
    """
        Returns array's vector.
    """
    try:
        angle = math.atan(vector[1]/float(vector[0]))
    except ZeroDivisionError:
        angle = math.pi
    return angle

def compress(obj, level=settings.default_comp_level, method=settings.default_comp):
    """
    Compress data of an object.
    """
    if level is None:
        return obj
    if method==ZLIB:
        dumps = cPickle.dumps(obj, -1)
        #print "Compressing with ZLIB. Level "+str(level)
        compressed = zlib.compress(dumps, level)
    elif method==LZO:
        #level = max([level,1])
        dumps = cPickle.dumps(obj, -1)
        compressed = lzo.compress(dumps, level)
    elif method==LZ4_FAST:
        dumps = cPickle.dumps(obj, -1)
        compressed = lz4.compress_fast(dumps, level)
    elif method==LZ4:
        dumps = cPickle.dumps(obj, -1)
        compressed = lz4.compress(dumps)#, level)
    else:
        compressed = obj
    return compressed

def decompress(obj, level=settings.default_comp_level, method=settings.default_comp):
    """
    Decompress an object.
    """
    if level is None:
        return obj
    if method==ZLIB:
        dumps = zlib.decompress(obj)#, level)
    elif method==LZO:
        #level = max([level,1])
        dumps = lzo.decompress(obj, level)
    elif method==LZ4_FAST:
        dumps = lz4.decompress(obj)
    elif method==LZ4:
        dumps = lz4.decompress(obj)
    else:
        return obj
    data = cPickle.loads(dumps)
    return data

def compress_state(obj):
    """
    wrapper
    """
    return compress(obj, level=settings.internal_level)

def decompress_state(obj):
    """
    wrapper
    """
    return decompress(obj, level=settings.internal_level)

def compress_rod(obj):
    """
    wrapper
    """
    return compress(obj, level=settings.rod_level)

def decompress_rod(obj):
    """
    wrapper
    """
    return decompress(obj, level=settings.rod_level)

def gaussian(distance, sigma=settings.sigma):
    """
    Returns gaussian probability for distance.
    """
    if not sigma:
        return 1
    value = math.exp(-(distance**2)/(2.0*sigma**2))
    return value

def norm_gaussian(distance, rad=float('inf'), sigma=settings.sigma):
    """
    Returns gaussian prob for distance (normalized to circle).
    """
    if sigma==0:
        return 1
    norm = 1.0/(math.sqrt(2.0*math.pi)*sigma*math.erf(rad*1.0/(math.sqrt(2)*sigma)))
    return gaussian(distance, sigma=sigma)*norm

def compute_distances(array1, array2):
    """
    Wrapper
    """
    if using_cl:
        return compute_distances_cl(array1, array2)
    else:
        distances = []
        for index in range(len(array1)):
            point1 = array1[index]
            point2 = array2[index]
            distances.append(distance_between_points(point1, point2))
        return distances

def compute_distances_cl(array1, array2):
    """
    Computes distances with gpu
    """
    mf = cl.mem_flags
    x_array_1 = []
    y_array_1 = []
    x_array_2 = []
    y_array_2 = []
    for value in array1:
        x_array_1.append(value[0])
        y_array_1.append(value[1])
    for value in array2:
        x_array_2.append(value[0])
        y_array_2.append(value[1])
    x_array_1 = numpy.array(x_array_1).astype(numpy.float32)
    x_array_2 = numpy.array(x_array_2).astype(numpy.float32)
    y_array_1 = numpy.array(y_array_1).astype(numpy.float32)
    y_array_2 = numpy.array(y_array_2).astype(numpy.float32)
    a_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=x_array_1)
    b_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=x_array_2)
    c_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=y_array_1)
    d_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=y_array_2)
    dest_buf = cl.Buffer(ctx, mf.WRITE_ONLY, x_array_1.nbytes)
    prg = cl.Program(ctx, """
    __kernel void sum(__global const float *x1,
    __global const float *x2, __global float *y1,
    __global const float *y2, __global float *res)
    {
      int gid = get_global_id(0);
      res[gid] = sqrt(pown((x1[gid]-x2[gid]),2)+pown((y1[gid]-y2[gid]),2));
    }
    """).build()
    prg.sum(queue, x_array_1.shape, None, a_buf, b_buf, c_buf, d_buf, dest_buf)
    output = numpy.empty_like(x_array_1)
    cl.enqueue_copy(queue_cl, output, dest_buf)
    return list(output)

import matplotlib.pyplot as plt

def animate_scatter(x_val, y_val, z_vals,
                        divisions, name, z_max, z_min, units, radius, title):
    """
    Specific animator.
    """
    try:
        z_val = decompress(z_vals.pop(0), level=settings.default_comp_level)
        plt.cla()
        plt.clf()
        rad = float(radius*1.3)/divisions
        size = (rad/4)**2
        x_min = min(x_val)-rad*1.1
        x_max = max(x_val)+rad*1.1
        y_min = min(y_val)-rad*1.1
        y_max = max(y_val)+rad*1.1
        plt.xlim((x_min, x_max))
        plt.ylim((y_min, y_max))
        scat = plt.scatter(x_val, y_val, s=size, c=z_val, marker='s',
                                            vmin=z_min, vmax=z_max)
        scat.cmap.set_under('w')
        scat.cmap.set_over('w')
        plt.gca().invert_yaxis()
        #step = int((z_max-z_min)*10.0)/100.0
        #ticks = [int((z_min + index*step)*10.0)/10.0 for index in range(10+1)]
        try:
            cb = plt.colorbar()#ticks=ticks)
        except TypeError:
            print z_val
        plt.xlabel("x [pixels]")
        plt.ylabel("y [pixels]")
        if not title is None:
            plt.suptitle(title)
        cb.set_label(units)
    except IndexError:
        pass

def create_scatter_animation(x_val, y_val, z_vals_avg, divisions, z_max, z_min, units, name, radius=800, fps=15, title=None):
    """
    Creates animation from data.
    """
    print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Creating animation"
    fig = plt.figure()
    num_frames = len(z_vals_avg)
    def animate(dummy_frame):
        """
        Wrapper.
        """
        animate_scatter(x_val, y_val, z_vals_avg,
                            divisions, name, z_max, z_min, units, radius, title)
    anim = animation.FuncAnimation(fig, animate, frames=num_frames)
    anim.save(name, writer=WRITER, fps=fps)

def animate_vector_map(x_val, y_val, u_vals, v_vals, units, name, radius, scale):
    """
    Specific animator.
    """
    try:
        u_val = decompress(u_vals.pop(0), level=settings.default_comp_level)
        v_val = decompress(v_vals.pop(0), level=settings.default_comp_level)
    except IndexError:
        return
    plt.cla()
    plt.clf()
    rad = radius/3.0
    size = (rad/4)**2
    x_min = min(x_val)-rad*1.1
    x_max = max(x_val)+rad*1.1
    y_min = min(y_val)-rad*1.1
    y_max = max(y_val)+rad*1.1
    plt.xlim((x_min, x_max))
    plt.ylim((y_min, y_max))
    try:
        plt.quiver(x_val, y_val, u_val, v_val, label=units, pivot="mid")#, units='x')#, scale=rad*scale*100)
    except TypeError:
        print (x_val, y_val, u_val, v_val)
    plt.gca().invert_yaxis()
    plt.xlabel("x [pixels]")
    plt.ylabel("y [pixels]")
    #plt.legend()

def create_vector_map(x_vals, y_vals, u_vals, v_vals, units, name, scale, radius=800, fps=15):
    """
    Creates animation from data.
    """
    print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Creating vector animation"
    fig = plt.figure()
    frames = len(u_vals)
    def animate(dummy_frame):
        """
        Wrapper.
        """
        animate_vector_map(x_vals, y_vals, u_vals, v_vals, units, name, radius, scale)
    anim = animation.FuncAnimation(fig, animate, frames=frames)
    anim.save(name, writer=WRITER, fps=fps)


def import_and_plot(source, radius=None, level=9):
    """
    Imports data from a compressed file and plots it.
    """
    print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Plotting " + str(source)
    src = open(source, 'r')
    name = source[:-5]
    values = decompress(src.read(), level=level)
    x_val, y_val, z_vals_avg = values[0], values[1], values[2]
    divisions = values[3]
    z_max, z_min = values[4], values[5]
    units = values[6]
    src.close()
    create_scatter_animation(x_val, y_val, z_vals_avg, divisions, z_max, z_min, units, name)

def needed_rods(wanted_long_prop, wanted_area_prop, area, long_area, short_area, longs, shorts):
    """
    Computes needed long rods / short rods to be added/removed.
    """
    long_rod_area = float(long_area*prop_long)/longs
    short_rod_area = float(short_area*prop_short)/shorts
    needed_longs = int(wanted_long_prop*wanted_area_prop*area/long_rod_area)
    needed_shorts = int(needed_longs*long_rod_area*(1-wanted_long_prop)/(wanted_long_prop*short_rod_area))
    return needed_longs, needed_shorts


def plot_all_data_files(radius=None, level=9, folder="./"):
    """
    Import all .data files and plot its data.
    """
    print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Plotting data from files"
    regular_expression = re.compile(r'.*\.data')
    dens_re = re.compile(r'.*dens.*\.data')
    q2_q4_re = re.compile(r'.*q[2-4].*\.data')
    temp_re = re.compile(r'.*temp.*\.data')
    cluster_re = re.compile(r'.*cluster.*\.data')
    files = get_file_names(folder=folder, regular_expression=regular_expression)
    for file_ in files:
        if dens_re.match(file_):
            print "Ploting "+file_
            import_and_plot(file_)
        elif q2_q4_re.match(file_):
            print "Ploting "+file_
            import_and_plot(file_)
        else:
            print "Ignoring "+file_+". Method not implemented yet."



def reset_dates_ids(folder="./", start=0):
    """
    When imgs goes over id limit, they restart in start(default 0), so
    """
    print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Resetting dates ids"
    src = folder + "dates.txt"
    input_ = open(src, 'r')
    output_file = folder + "dates_reordered.txt"
    output_ = open(output_file, 'w')
    id_ = start
    for line in input_:
        date = line.split("\t")[1]
        line = ""
        if id_<10:
            line += "000"
        elif id_<100:
            line += "00"
        elif id_<1000:
            line += "0"
        line += str(id_) + "\t" + str(date)
        output_.write(line)
        id_ += 1
    output_.close()
    input_.close()

def print_progress(done, total, counter, times, time_left, previous_time, counter_refresh=settings.counter_refresh):
    """
    Print progress of stack of tasks.
    """
    msg = "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]: "
    left = total-done
    now = datetime.datetime.now()
    seconds_passed = (now-previous_time).total_seconds()
    times.append(seconds_passed)
    progress = int(done*100/total)
    previous_time = now
    string = msg + " %d%%  " % (progress)
    perten = progress/10.0
    string += "["
    prog = int(perten*4)
    string += WHITE_BLOCK*prog
    string += " "*(40-prog)
    string += "]"
    if counter >= counter_refresh:
        counter = 0
    try:
        subtimes = []
        if counter == 0:
            subtimes = times
        else:
            subtimes = times[:-counter]
        avg_time = sum(subtimes)*1.0/(len(times)-counter)
        time_left = int(left*avg_time/60)
    except ZeroDivisionError:
        avg_time = None
        time_left = None
    if not time_left is None and not avg_time is None:
        if time_left:
            string += "\t" + str(time_left) + " minutes"
        else:
            string += "\t" + str(int(left*avg_time)) + " seconds"
    print CLEAR_LAST
    print string
    return previous_time, counter, time_left, times

def rods_differences(shorts, longs, wlprop):
    """
    Computes how many rods must be included.
    """
    diff_l = int(((wlprop-1)*longs+wlprop*shorts)*1.0/(1+wlprop))
    diff_s = -2*diff_l
    Nl = diff_l + longs
    Ns = diff_s + shorts
    return diff_l, diff_s, Nl, Ns

def rods_animation(rods, colours, x_lim, y_lim, zone_coords, name="rods.mp4", fps=1):
    """
    Plot rods with colours
    """
    print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Creating rods animation"
    fig = plt.figure()
    frames = len(rods[0])
    rods_12, rods_6 = rods[0], rods[1]
    #_export_rods(rods_12, 12)
    #_export_rods(rods_6, 6)
    decompressed_rods_12 = []
    decompressed_rods_6 = []
    output_queue = mp.Queue()
    processes = []
    print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Decompressing"
    for index in range(frames):
        rods_12_ = rods_12[index]
        rods_6_ = rods_6[index]
        process = mp.Process(target=_decompress_rods_process,
                            args=(index, output_queue, rods_12_, rods_6_))
        processes.append(process)
        decompressed_rods_12.append(None)
        decompressed_rods_6.append(None)
    num_processes = len(processes)
    running, processes_left = run_processes(processes)
    finished = 0
    previous_time = datetime.datetime.now()
    counter = 0
    time_left = None
    times = []
    print " "
    while finished < num_processes:
        counter += 1
        finished += 1
        previous_time, counter, time_left, times = print_progress(finished, num_processes,
                            counter, times, time_left, previous_time)
        [index, decompressed_rods_12_, decompressed_rods_6_] = output_queue.get()
        decompressed_rods_12[index] = decompressed_rods_12_
        decompressed_rods_6[index] = decompressed_rods_6_
        if len(processes_left):
            new_process = processes_left.pop(0)
            #time.sleep(settings.waiting_time)
            new_process.start()
    print CLEAR_LAST
    print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Animating"
    def animate(dummy_frame):
        """
        Wrapper.
        """
        try:
            rods_12_ = decompressed_rods_12.pop(0)
            rods_6_  = decompressed_rods_6.pop(0)
            rods_ = [rods_12_, rods_6_]
        except IndexError:
            return
        animate_rods(rods_, colours, x_lim, y_lim, zone_coords)
    anim = animation.FuncAnimation(fig, animate, frames=frames)
    anim.save(name, writer=WRITER, fps=fps)

def _decompress_rods_process(index, output_queue, rods_12, rods_6):
    """
    Process
    """
    decomp_12 = decompress(rods_12)
    decomp_6 = decompress(rods_6)
    output_queue.put([index, decomp_12, decomp_6])

def animate_rods(rods, colours, x_lim, y_lim, zone_coords):
    """
    Specific animator.
    """
    rods_12, rods_6 = rods[0], rods[1]
    plt.cla()
    plt.clf()
    plt.xlim(x_lim)
    plt.ylim(y_lim)
    #circle1 = plt.Circle((zone_coords[0], zone_coords[1]), zone_coords[2]+4, color='black')
    #circle2 = plt.Circle((zone_coords[0], zone_coords[1]), zone_coords[2], color='white')
    #plt.gca().add_artist(circle1)
    #plt.gca().add_artist(circle2)
    # rods_12 = [x_0_list, y_0_list, x_f_list, y_f_list]
    plt.plot([rods_12[0], rods_12[2]], [rods_12[1], rods_12[3]], c=colours[0], linewidth=1.5)
    plt.plot([rods_6[0], rods_6[2]], [rods_6[1], rods_6[3]], c=colours[1], linewidth=1.5)
    plt.gca().invert_yaxis()
    #plt.gca().invert_xaxis()
    plt.xlabel("x [pixels]")
    plt.ylabel("y [pixels]")

def _export_rods(rods, kappa):
    """
    Exports rods to file.
    """
    output_queue = mp.Queue()
    processes = []
    print "\n\n"
    print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Exporting data"
    for index in range(len(rods)):
        rods_ = rods[index]
        process = mp.Process(target=_export_rods_process,
                            args=(index, output_queue, rods_, kappa))
        processes.append(process)
    num_processes = len(processes)
    running, processes_left = run_processes(processes)
    finished = 0
    previous_time = datetime.datetime.now()
    counter = 0
    time_left = None
    times = []
    print " "
    while finished < num_processes:
        counter += 1
        finished += 1
        previous_time, counter, time_left, times = print_progress(finished, num_processes,
                            counter, times, time_left, previous_time)
        output_queue.get()
        if len(processes_left):
            new_process = processes_left.pop(0)
            #time.sleep(settings.waiting_time)
            new_process.start()
    print CLEAR_LAST

def _export_rods_process(index, output, rods, kappa):
    """
    process
    """
    name = str(index) + "_"+str(kappa)+".data"
    file_ = open(name, 'w')
    rods_ = decompress(rods)
    for index in range(len(rods_[0])):
        data = [rods_[0][index], rods[1][index], rods_[2][index], rods_[3][index]]
        line = str(data[0]) + "\t" + str(data[1]) + "\n"
        line += str(data[2]) + "\t" + str(data[3]) + "\n"
        line += "\n"
        file_.write(line)
    file_.close()
    output.put(None)

import copy

def order_param_animation(matrices_12, matrices_6, divisions, bursts_groups, bursts_times):
    """
    Computes order param.
    """
    print "\n\n"
    print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Computing order param and exporting"
    output_queue = mp.Queue()
    processes = []
    groups = copy.deepcopy(bursts_groups)
    bursts_ = len(groups)
    param_order_matrices = []
    for index in range(len(matrices_12)):
        matrix_12 = matrices_12[index]
        matrix_6 = matrices_6[index]
        name = str(index) + "_order_param.data"
        process = mp.Process(target=_order_param_process,
                            args=(index, output_queue, matrix_12, matrix_6, name))
        processes.append(process)
        param_order_matrices.append(None)
    num_processes = len(processes)
    running, processes_left = run_processes(processes)
    finished = 0
    previous_time = datetime.datetime.now()
    counter = 0
    time_left = None
    times = []
    print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Getting values"
    print " "
    x_val, y_val, z_vals = [], [], []
    while finished < num_processes:
        counter += 1
        finished += 1
        previous_time, counter, time_left, times = print_progress(finished, num_processes,
                            counter, times, time_left, previous_time)
        index, x_val, y_val, z_val = output_queue.get()
        z_vals.append(z_val)
        if len(processes_left):
            new_process = processes_left.pop(0)
            #time.sleep(settings.waiting_time)
            new_process.start()
    print CLEAR_LAST
    if settings.plot:
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Computing averages"
        z_vals_avg = []
        output_queue = mp.Queue()
        for index in range(len(groups)):
            group = groups[index]
            z_vals_avg.append(None)
            process = mp.Process(target=_compute_z_averages_process,
                                args=(index, output_queue, group, z_vals))
            processes.append(process)
            param_order_matrices.append(None)
        num_processes = len(processes)
        running, processes_left = run_processes(processes)
        finished = 0
        previous_time = datetime.datetime.now()
        counter = 0
        time_left = None
        times = []
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Getting values"
        print " "
        while finished < num_processes:
            counter += 1
            finished += 1
            previous_time, counter, time_left, times = print_progress(finished, num_processes,
                                counter, times, time_left, previous_time)
            [index__, avg_z_val] = output_queue.get()
            z_vals_avg[index__] = avg_z_val
            if len(processes_left):
                new_process = processes_left.pop(0)
                #time.sleep(settings.waiting_time)
                new_process.start()
        print CLEAR_LAST
        frames = len(z_vals_avg)
        assert frames, "Not values! \n\n"
        z_max = 1
        z_min = -1
        units = "normalized [S.U.]"
        name = "order_param.mp4"
        radius = 800
        title = "Order parameter"
        z_vals_copy = copy.deepcopy(z_vals_avg)
        print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Plotting"
        create_scatter_animation(x_val, y_val, z_vals_copy, divisions, z_max, z_min, units, name, radius, title=title)
        if settings.order_param_lengths:
            print "--"*(len(inspect.stack())-2)+">"+"["+str(inspect.stack()[1][3])+"]->["+str(inspect.stack()[0][3])+"]: " + "Computing contour lengths"
            lengths = []
            for index in range(len(z_vals_avg)):
                #fig = plt.figure()
                z_val = decompress(z_vals_avg[index])
                x_grid = numpy.linspace(min(x_val), max(x_val), settings.grid_length)
                y_grid = numpy.linspace(min(y_val), max(y_val), settings.grid_length)
                x_val = numpy.array(x_val)
                y_val = numpy.array(y_val)
                z_val = numpy.array(z_val)
                len(x_grid)
                z_grid = griddata((x_val, y_val), z_val, (x_grid[None,:], y_grid[:,None]), method='cubic')
                levels = [0]
                cs = plt.contour(x_grid, y_grid, z_grid, linewidths=1.25, colors='k', levels=levels)
                length = 0  
                x0,y0 =  cs.allsegs[0][0][0]
                startx = x0
                starty = y0
                length = 0
                for coords in cs.allsegs[0][0][1:]:
                    x1,y1 =  coords[0], coords[1]
                    length += numpy.sqrt((x1-x0)**2 + (y1-y0)**2)
                    x0,y0 = x1,y1
                length += numpy.sqrt((startx-x0)**2 + (starty-y0)**2)
                lengths.append(length)
            fig = plt.figure()
            plt.scatter(bursts_times, lengths)
            plt.savefig("cluster_boundaries_length.png")


def _compute_z_averages_process(index, output_queue, group, z_vals):
    """
    Process
    """
    average = []
    for index_ in group:
        average.append(z_vals[index_])
    average = compress(array_average(average), level=settings.default_comp_level)
    output_queue.put([index, average])


def _order_param_process(index, output_queue, matrix_12, matrix_6, name):
    """
    process
    """
    file_ = open(name, 'w')
    z_vals = []
    matrix_12 = decompress(matrix_12)
    matrix_6 = decompress(matrix_6)
    x_val_12 = matrix_12[0]
    y_val_12 = matrix_12[1]
    for row_index in range(len(matrix_12[0])):
        n_12 = matrix_12[2][row_index]
        n_6 = matrix_6[2][row_index]
        x_val = x_val_12[row_index]
        y_val = y_val_12[row_index]
        if n_12+n_6 == 0:
            order_param = -1000
        else:
    	    order_param = float(n_12-n_6)/(n_12+n_6) #float(n_12)#float(n_12-n_6+1)/2
        line = str(x_val)+"\t"+str(y_val)+"\t"+str(order_param)+"\n"
        file_.write(line)
        z_vals.append(order_param)
    file_.close()
    output_queue.put([index, x_val_12, y_val_12, z_vals])

    
