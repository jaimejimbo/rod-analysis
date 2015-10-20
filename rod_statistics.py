"""
Analyze rods data
"""
from queue import Queue
from glob import glob
import re

class Rod(object):
    """
    Rod object.
    """
    
    def __init__(self, ARGS):
        """
        Initialization of rod
        """
        pass

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
    

def import_files(folder="./", _glob='rods_*', regular_expression='rods_[0-9]{1,5}'):
    """
    Import all files using glob and checking with reg exp.
    """
    reg1 = re.compile(regular_expression)
    check_if_data = re.compile(reg1)
    names = []
    for data_file in glob(_glob):
        if reg1.match(data_file):             #check if file is a data file
            names.append(data_file)
    files = []
    for name in names:
        files.append(open(name))
    return files

def import_data(file_, regular_expression='*'):
    """
    Import data of a file (only works with rod data)
    """
    reg_exp = re.compile(regular_expression)
    rod = None                     #Rod(prop1, prop2)
    return rod
