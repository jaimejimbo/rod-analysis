"""
Analyze rods data
"""
import re
#import multiprocessing as mp    #for using all cores
import os
from rod import Rod
from methods import binary_order
from system_state import SystemState



def import_files(folder="./", regular_expression='rods_[0-9]*'):
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


def get_file_names(folder="./", regular_expression='rods_[0-9]*'):
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
    return binary_order(names, file_name_ordering_function)



def file_name_ordering_function(name):
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
    states = []
    for index in range(len(files)):
        state = SystemState(kappas, allowed_kappa_error,
                   radius_correction_ratio, names[index])
        data = import_data(files[index])
        for dataline in data:
            parameters = tuple(dataline)
            new_rod = Rod(parameters)
            state.add_rod(new_rod)
        state.check_rods()
        states.append(state)
    return names, states

