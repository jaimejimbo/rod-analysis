"""
Analyze rods data
"""
from queue import Queue

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
	

def import_files(folder, ext):
	"""
	Import all files with given ext
	"""
	files = []
	return files

def import_data(file_):
	"""
	Import data of a file (only works with rod data)
	"""
	reg_exp = None			#Regular expression to parse data file
	rod = None 				#Rod(prop1, prop2)
	return rod
