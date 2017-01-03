import re, os, sys
import numpy as np


def time_lapse(date1, date2):
    """
    Time lapse in seconds between dates.
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
    return abs(time2-time1)

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

sqldb = ''#input("sqlite3 db src-> ")
experiment_id = get_number_from_string(cwd)
if sqldb == '':
    sqldb = "../rods.db"
if experiment_id == '':
    raise ValueError

dates = import_image_dates()
keys = list(dates.keys())
orig = keys[0]

import sqlite3 as sql

str_ = "inserting data..."
print(str_)
conn = sql.connect(sqldb)
cursor = conn.cursor()

for key in keys:
    file_id = key
    tim_lap = time_lapse(dates[key], dates[orig])
    cursor.execute("insert into times values(?,?,?)", (experiment_id, file_id, tim_lap))

conn.commit()
conn.close()    
