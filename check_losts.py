#!/usr/bin/python
"""

"""

"""
Settings
"""
#0-9 or None to disable
low_comp_level = None
medium_comp_level = 0
strong_comp_level = None
#Reduces cpu usage (0 for disable)
waiting_time = 0
if not waiting_time:
    waiting_time = 0    
waiting_time_proccess = 0
if not waiting_time:
    waiting_time_process = 0  
#None uses all cpus      
cpus = None
# 1 or True, 0 or False
create_videos = 0
clusters = 1
avg_temp = 1
order_param_exp = 1
lost_percentage = 1
run_imagej = 0
run_props = 1
run_graphs = 0
run_check_dim = 0
get_image_dates = 0
to_file = 1
plot = 0
only_density = 0
# variables
coef = 5
divisions = 50
#sigma = sigma_coef * subsystem_rad
sigma_coef = None
changing_props = 0
discard_exceptions = 0
special_chars = 0
import re


"""
Script
"""
def get_coords_from_imagej():
    """
        Gets center and radius from imagej script.
    """
    coords = []
    imagej_script = open("acquisition_v4.ijm", "r")
    reg_exp = re.compile(".*makeOval.*")
    for line in imagej_script:
        if reg_exp.match(line):
            data_line = line
            break
    zone_coords_re = re.compile(".*makeOval\((.*),(.*),(.*),(.*)\);.*")
    results_ = zone_coords_re.search(data_line)
    results = [results_.group(index+1) for index in range(4)]
    imagej_script.close()
    rad = float(results[2])/2
    center_x = float(results[1])+rad 
    center_y = float(results[0])+rad
    zone_coords_ = [center_x, center_y, rad]
    return zone_coords_

def set_coords_in_imagej(zone_coords_):
    """
        Sets values in imagej script.
    """
    imagej_script = open("acquisition_v4.ijm", "w")
    reg_exp = re.compile(".*makeOval.*")
    center_x = zone_coords[0]
    center_y = zone_coords[1]
    rad = zone_coords[2]
    values = "("+str(center_y-rad)+","+str(center_x-rad)+","+str(rad*2)+","+str(rad*2)+");"
    for line in imagej_script:
        if reg_exp.match(line):
            line = "\tmakeOval"+values
    imagej_script.close()
    
dates = None
        

#zone_coords = []
#set_coords_in_imagej()
zone_coords = get_coords_from_imagej()
zone_coords = [zone_coords[1], zone_coords[0], zone_coords[2]]
#print zone_coords
#raw_input("")
settings = open("settings.py", "w")
settings.write("low_comp_level = {0}".format(low_comp_level))
settings.write("\n")
settings.write("medium_comp_level = {0}".format(medium_comp_level))
settings.write("\n")
settings.write("strong_comp_level = {0}".format(strong_comp_level))
settings.write("\n")
settings.write("cpus = {0}".format(cpus))
settings.write("\n")
settings.write("zone_coords = {0}".format(zone_coords))
settings.write("\n")
settings.write("to_file = {0}".format(to_file))
settings.write("\n")
settings.write("plot = {0}".format(plot))
settings.write("\n")
settings.write("waiting_time = {0}".format(waiting_time))
settings.write("\n")
settings.write("special_chars = {0}".format(special_chars))
settings.write("\n")
settings.write("waiting_time_process = {0}".format(waiting_time_process))
settings.write("\n")
rad = zone_coords[2]
if sigma_coef:
    gaussian_sigma = sigma_coef*float(rad)/divisions
else:
    gaussian_sigma = None
settings.write("sigma = {0}".format(gaussian_sigma))
settings.write("\n")
settings.close()

import experiment, system_state, methods, settings
import os, gc, math, time#, guppy


if get_image_dates:
    os.system("rm dates.txt")
    dates = methods.get_image_dates()
    methods.export_image_dates()
else:
    dates = methods.import_image_dates()
    
    
def plot_nums(kappa, real_kappa, error):
    msg = "\t\t\tK"+str(real_kappa)+"\t\t\t"
    print msg
    names, rod_groups = system_state.create_rods(kappas=kappa, real_kappas=real_kappa,
                                            allowed_kappa_error=error)
    experiment_ = experiment.Experiment(system_states_name_list=names, kappas=real_kappa,
                                        system_states_list=rod_groups,
                                        dates=dates, diff_t=5/3.0)
    experiment_.plot_number_of_rods_over_time()


def plot_nums_length(length, length_error, real_kappa):
    """
        Gets proportions of long/shorts rods.
    """
    msg = "\t\t\tK"+str(real_kappa)+"\t\t\t"
    names, rod_groups = system_state.create_rods_with_length(length=length, length_error=length_error,
                                            real_kappas=real_kappa)
    experiment_ = experiment.Experiment(system_states_name_list=names, kappas=real_kappa,
                                        system_states_list=rod_groups,
                                        dates=dates, diff_t=5/3.0)
    experiment_.plot_number_of_rods_over_time()

#plot_nums(50, 50, 50)
plot_nums_length(160, 20, 12)
plot_nums_length(80, 20, 6)
#plot_nums(8, 6, 5)
