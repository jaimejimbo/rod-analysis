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
run_12 = 1
run_4 = 1
run_all = 1
create_videos = 1
clusters = 1
avg_temp = 1
order_param_exp = 1
lost_percentage = 1
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

    if run_12:
        print "\t\t\tK12\t\t\t"
        names, rod_groups_12 = system_state.create_rods(kappas=17.5, real_kappas=12,
                                                allowed_kappa_error=1.5)
        experiment_12 = experiment.Experiment(system_states_name_list=names, kappas=12,
                                            system_states_list=rod_groups_17,
                                            dates=dates, diff_t=5/3.0)
        experiment_12.set_coef(coef)
        
        if create_videos:
            experiment_12.divide_systems_in_circles(divisions=divisions)
            experiment_12.create_videos(divisions=divisions, fps=10, max_distance=10, max_angle_diff=5,
                                     number_of_bursts=1)
        if clusters:
            experiment_12.plot_cluster_areas(number_of_bursts=3, max_distance=100,
                    max_angle_diff=10, min_size=10)
            experiment_12.create_cluster_histogram_video(max_distance=100, max_angle_diff=10, fps=15)
        if avg_temp:
            experiment_12.plot_average_temperature(100, 10)
        if order_param_exp:
            indep, param, std_dev = experiment_12.get_order_evolution_coeficient(number_of_bursts=3, max_distance=50,
                    max_angle_diff=10, min_size=10)
            print "Order parameter: "+str(param)+" Standard deviation: "+str(std_dev)
        if lost_percentage:
            percentage, std_dev = experiment_12.lost_rods_percentage
            print "Rods lost: "+str(percentage)+"%"+" Standard deviation: "+str(std_dev)
    try:
        names = None
        experiment_12 = None #methods.compress(experiment_12, level=settings.strong_comp_level)
        rod_groups_12 = None
        gc.collect()
    except NameError:
        pass
    print "\n"
    if run_4:
        print "\t\t\tK4\t\t\t"
        names, rod_groups_4 = system_state.create_rods(kappas=5.5, real_kappas=4,
                                                allowed_kappa_error=.5)
        experiment_4 = experiment.Experiment(system_states_name_list=names, kappas=4,
                                            system_states_list=rod_groups_5,
                                            dates=dates, diff_t=5/3.0)
        experiment_4.set_coef(coef)
        if create_videos:
            experiment_4.divide_systems_in_circles(divisions=divisions)
            experiment_4.create_videos(divisions=divisions, fps=10, max_distance=10, max_angle_diff=5,
                                     number_of_bursts=1)
        if clusters:
            #experiment_4.plot_cluster_areas(number_of_bursts=5, max_distance=100,
            #        max_angle_diff=10, min_size=20)
            pass
        if avg_temp:
            experiment_4.plot_average_temperature(100, 10)
        if lost_percentage:
            percentage, std_devar = experiment_4.lost_rods_percentage
            print "Rods lost: "+str(percentage)+"%"+" Standard deviation: "+str(std_dev)
    try:
        names = None
        experiment_4 = None #methods.compress(experiment_4, level=settings.strong_comp_level)
        rod_groups_4 = None
        gc.collect()
    except NameError:
        pass
    print "\n"
    if run_all:
        print "\t\t\tK all\t\t\t"
        names, rod_groups_all = system_state.create_rods(kappas=10, real_kappas=1,
                                                allowed_kappa_error=10)
        experiment_all = experiment.Experiment(system_states_name_list=names, kappas=1,
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
