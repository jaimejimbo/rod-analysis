#!/usr/bin/python -i

"""
Main script.
"""

"""
Settings
"""
#0-9 or None to disable
low_comp_level = None
medium_comp_level = None
strong_comp_level = 0 
#None uses all cpus      
cpus = 12
# 1 or True, 0 or False
run_17 = 0
run_5 = 1
run_all = 1
create_videos = 1
clusters = 0
avg_temp = 0
order_param_exp = 0
lost_percentage = 0
run_imagej = 0
# variables
coef = 3
divisions = 30










"""
Script
"""


settings = open("settings.py", "w")
settings.write("low_comp_level = {0}".format(low_comp_level))
settings.write("\n")
settings.write("medium_comp_level = {0}".format(medium_comp_level))
settings.write("\n")
settings.write("strong_comp_level = {0}".format(strong_comp_level))
settings.write("\n")
settings.write("cpus = {0}".format(cpus))
settings.write("\n")
settings.close()

import experiment, system_state, methods, settings
import os, gc

os.system("clear && clear")

#try:
if True:

    dates = None
    if run_imagej:
        print "Running imagej..."
        methods.imagej()
        dates = methods.get_image_dates()
        methods.export_image_dates()
    else:
        dates = methods.import_image_dates()

    if run_17:
        print "\t\t\tK17\t\t\t"
        names, rod_groups_17 = system_state.create_rods(kappas=17.5,
                                                allowed_kappa_error=1.5)
        experiment_17 = experiment.Experiment(system_states_name_list=names,
                                            system_states_list=rod_groups_17,
                                            dates=dates, diff_t=5/3.0)
        experiment_17.set_coef(coef)
        
        if create_videos:
            experiment_17.divide_systems_in_circles(divisions=divisions)
            experiment_17.create_videos(divisions=divisions, fps=10, max_distance=10, max_angle_diff=5,
                                     number_of_bursts=1)
        if clusters:
            experiment_17.plot_cluster_areas(number_of_bursts=3, max_distance=100,
                    max_angle_diff=10, min_size=10)
        if avg_temp:
            experiment_17.plot_average_temperature(100, 10)
        if order_param_exp:
            indep, param, std_dev = experiment_17.get_order_evolution_coeficient(number_of_bursts=3, max_distance=50,
                    max_angle_diff=10, min_size=10)
            print "Order parameter: "+str(param)+" Standard deviation: "+str(std_dev)
        if lost_percentage:
            percentage, std_dev = experiment_17.lost_rods_percentage
            print "Rods lost: "+str(percentage)+"%"+" Standard deviation: "+str(std_dev)
    try:
        names = None
        experiment_17 = None #methods.compress(experiment_17, level=settings.strong_comp_level)
        rod_groups_17 = None
        gc.collect()
    except NameError:
        pass
    print "\n"
    if run_5:
        print "\t\t\tK5\t\t\t"
        names, rod_groups_5 = system_state.create_rods(kappas=5.5,
                                                allowed_kappa_error=.5)
        experiment_5 = experiment.Experiment(system_states_name_list=names,
                                            system_states_list=rod_groups_5,
                                            dates=dates, diff_t=5/3.0)
        experiment_5.set_coef(coef)
        if create_videos:
            experiment_5.divide_systems_in_circles(divisions=divisions)
            experiment_5.create_videos(divisions=divisions, fps=10, max_distance=10, max_angle_diff=5,
                                     number_of_bursts=1)
        if clusters:
            experiment_5.plot_cluster_areas(number_of_bursts=5, max_distance=100,
                    max_angle_diff=10, min_size=20)
        if avg_temp:
            experiment_5.plot_average_temperature(100, 10)
        if lost_percentage:
            percentage, std_devar = experiment_5.lost_rods_percentage
            print "Rods lost: "+str(percentage)+"%"+" Standard deviation: "+str(std_dev)
    try:
        names = None
        experiment_5 = None #methods.compress(experiment_5, level=settings.strong_comp_level)
        rod_groups_5 = None
        gc.collect()
    except NameError:
        pass
    print "\n"
    if run_all:
        print "\t\t\tK all\t\t\t"
        names, rod_groups_all = system_state.create_rods(kappas=10,
                                                allowed_kappa_error=10)
        experiment_all = experiment.Experiment(system_states_name_list=names,
                                            system_states_list=rod_groups_all,
                                            dates=dates, diff_t=5/3.0)
        experiment_all.set_coef(coef)
        if create_videos:
            experiment_all.divide_systems_in_circles(divisions=divisions)
            experiment_all.create_videos(divisions=divisions, fps=10, max_distance=10, max_angle_diff=5,
                                     number_of_bursts=1)
        if clusters:
            experiment_all.plot_cluster_areas(number_of_bursts=5, max_distance=100,
                    max_angle_diff=10, min_size=20)
        if avg_temp:
            experiment_all.plot_average_temperature(100, 10)
        if lost_percentage:
            percentage, std_devar = experiment_all.lost_rods_percentage
            print "Rods lost: "+str(percentage)+"%"+" Standard deviation: "+str(std_dev)
    try:
        names = None
        experiment_all = None #methods.compress(experiment_all, level=settings.strong_comp_level)
        rod_groups_all = None
        gc.collect()
    except NameError:
        pass

    #os.system("bash tomp4script.sh")
"""except KeyboardInterrupt:
    os.system("sudo pkill -f main.py")
except AssertionError as error:
    print "Assertion error({0})".format(error)
    os.system("sudo pkill -f main.py")"""
print "Type exit() to exit."
def exit():
    print "Give root password to kill all remaining zombie processes"
    os.system("sudo pkill -f main.py")
