import re, os, sys
import numpy as np


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
    return list(np.sort(names))

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
    #if str(type(_file)) != "<type 'file'>":
    #    print("Passed file argument is not a file descriptor")
    #    raise ValueError
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
                        print(str("Something went badly" + str(dataline)))
            data.append(dataline)
    except ValueError:
        print("File must be passed opened to import_data.")
    except TypeError:
        print("Perhaps the file passed is not in the right format.")
        print(_file)
    except IndexError:
        print("Error importing files (empty list)")
        raise IndexError
    return data


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
                    
cwd = os.getcwd()

data = []
str_ = "getting file names..."
print(str_)
names = get_file_names()

str_ = "importing data..."
print(str_)
for name in names:
    file = open(name, 'r')
    data.append(import_data(file))
    file.close()

#(experiment_id, file_id, ID, area, xmid, ymid, major, minor, angle, feret, feretx, ferety, feretangle, minferet, xstart, ystart)

sqldb = ''#input("sqlite3 db src-> ")
experiment_id = get_number_from_string(cwd)
if sqldb == '':
    sqldb = "../rods.db"
if experiment_id == '':
    raise ValueError

import sqlite3 as sql

str_ = "inserting data..."
print(str_)
conn = sql.connect(sqldb)
cursor = conn.cursor()
for index in range(len(data)):
    row = data[index]
    file_id = get_number_from_string(names[index])
    for rod in row:
        data_ = [experiment_id, file_id]
        for value in rod:
            data_.append(value)
        #print(len(data_))
        data_ = tuple(data_)
        try:
            cursor.execute("insert into datos values(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",data_)
        except sql.ProgrammingError:
            print(data_)

conn.commit()
#conn.execute("VACUUM")
conn.close()    
