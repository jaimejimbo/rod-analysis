"""
Analyze rods data
"""
from queue import Queue
from glob import glob
import re
import multiprocessing as mp    #for using all cores
import math



#consts
RADIUS = 1704*0.5
CENTER_X = 519+RADIUS
CENTER_Y = 96+RADIUS
KAPPA_1 = 4
KAPPA_2 = 12
#ROD_LENGTH,ROD_DIAMETER #Now there are two types of rod
# Time consts
FRAME_INTERVAL = 1.0/3



class Rod(object):
    """
    Rod object.
    """
    
    def __init__(self, (ID, area, xm, ym, major, minor, angle, feret, feretx, ferety, feretangle, minferet, xstart, ystart)):
        """
        Initialization of rod
        """                                    #Column
        self.id = int(ID)                      #0
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

    def is_valid_rod(self):
        """
        Checks if this is a rod looking at different factors
        If it is a group of two rods that are very near.
        Remove rods that are near the boundrie.
        """
        return True
    





class SystemState(object):
    """
    Group of rods in a moment.
    Each image has to be translated into a RodGroup (by this class?)
    """
    def __init__(self):
        """
        Initialization
        """
        self._rods = Queue()    #I use a Queue, but perhaps is better a tree or a recursive list
        self._number_of_particles = 0

    def add_rod(self, rod):
        """
        Adds a rod to the group
        """
        self._rods.join(rod)
        self._number_of_particles += 1

    def get_rod(self):
        """
        Returns the first rod in the queue
        The rod is removed of the group!
        """
        self._number_of_particles -= 1
        return self._rods.get_next()

    def remove_rod(self, rod):
        """
        Removes a rod from the group (queue object mod needed)
        """
        self._rods.delete(rod)
        self._number_of_particles -= 1
    
    def _compute_density(self):
        """
        Computes density of the system
        """
        pass



class RodGroup(SystemState):
    """
    Group of rods. Used to put all rods that are in a zone or
    have something in common.
    """
    
    def __init__(self, area):
        """
        Initialization
        """
        self._rods = Queue()    #I use a Queue, but perhaps is better a tree or a recursive list
        self._density = -1
        self._area = area


    def _update_density(self):
        """
        Overrides.
        Computes density of the group.
        """
        self._density = self._number_of_particles * 1.0/self._area

    def get_density(self):
        """
        Returns the density of the group
        Perhaps sometimes update_density method can be ignored
        """
        self._update_density()
        return self._density







def import_files(folder="./", _glob='rods_*', regular_expression='rods_[0-9]*'):
    """
    Import all files using glob and checking with reg exp.
    """
    reg1 = re.compile(regular_expression)
    check_if_data = re.compile(reg1)
    names = []
    for data_file in glob(_glob):
        if reg1.match(data_file):
            names.append(data_file)
    files = []
    for name in names:
        files.append(open(name))
    return files





def import_data(_file, split_char='\t', regular_expression='[0-9]\.?[0-9]*'):
    """
    Import data of a file
    Returns an array with data
    """
    reg_exp = re.compile(regular_expression)
    data = []
    try:
        for line in _file:            
            dataline = re.split(split_char, line.rstrip('\n')) # Removes new line char from the line
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
    return data





def create_rods(folder="./"):
    """
    Create one rod for each rod_data and for each file
    returns [RodGroup1, RodGroup2, ...]
    """
    files = import_files(folder=folder)
    states = []
    for _file in files:
        state = State()
        data = import_data(_file)
        for dataline in data:
            parameters = tuple(dataline)
            try:
                new_rod = Rod(parameters)
            except ValueError as e:
                print parameters
                print e.message
                raise ValueError
            state.add_rod(new_rod)
        states.append(rod_group)
    return rod_groups




#
#   IN DEVELOPMENT
#  Allows the creation of circles all with same area.
#  Computes the area of circles' intesection. 
#  Computes a radius of a cricle that, when intesecting with the main circle, has the same area than the others
#

def segment_area(r,h):
    """
    Computes the area of a small circle intersecting with the bigger one.
    r is the radius of the small circle.
    h is the distance from the center of the small circle to the intersection.
    is the circle much smaller? (R of the main circle -> infinity?)
    """
    return r**2 * math.acos(h/r) - h*sqrt(r**2-h**2)

def effective_area(r,r_pos, R):
    """
    r: radius of small circle
    R: radius of main circle
    r_pos: position of the small circle
    """
    h = (r_pos**2-r**2+R**2)/(2*r_pos)
    if h>=r_pos: 
        return math.pi*r**2 - segment_area(r,h-r_pos)+segment_area(R,h)
    else:
        return segment_area(r,r_pos-h)+segment_area(R,h)
