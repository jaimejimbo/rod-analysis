"""
Analyze rods data
"""
from queue import Queue
from glob import glob
import re
import multiprocessing as mp    #for using all cores

class Rod(object):
    """
    Rod object.
    """
    
    def __init__(self, (ID, area, xm, ym, major, minor, angle, feret, feretx, ferety, feretangle, minferet, xstart, ystart)):
        """
        Initialization of rod
        """
        self.id = ID
        self.area = area
        self.xm = xm
        self.ym = ym
        self.major = major
        self.minor = minor
        self.angle = angle
        self.feret = feretx
        self.ferety = ferety
        self.feretangle = feretangle
        self.minferet = minferet
        self.xstart = xstart
        self.ystart = ystart

    def check_if_rod(self):
        """
        Checks if this is a rod looking at different factors
        """
        return True
    


class RodGroup(object):
    """
    Group of rods.
    Each image has to be translated into a RodGroup (by this class?)
    """
    def __init__(self):
        """
        Initialization
        """
        self._rods = Queue()    #I use a Queue, but perhaps is better a tree or a recursive list

    def add_rod(self, rod):
        """
        Adds a rod to the group
        """
        self._rods.join(rod)

    def remove_rod(self, rod):
        """
        Removes a rod from the group (queue object mod needed)
        """
        pass
    

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

def import_data(_file, split_char='\t', regular_expression='[0-9]\.[0-9]*'):
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
    files = import_files()
    rod_groups = []
    for _file in files:
        rod_group = RodGroup()
        data = import_data(_file)
        for dataline in data:
            parameters = tuple(dataline)
            new_rod = Rod(parameters)
            rod_group.add_rod(new_rod)
        rod_groups.append(rod_group)
    return rod_groups
