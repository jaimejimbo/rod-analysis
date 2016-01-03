#!/usr/bin/python -i

"""
Main script.
"""

default_comp_level = 0
strong_comp_level = 1       #0-9

settings = open("settings.py", "w")
settings.write("default_comp_level = {0}".format(default_comp_level))
settings.write("\n")
settings.write("strong_comp_level = {0}".format(strong_comp_level))
settings.write("\n")
settings.close()

import experiment, system_state, methods
import os, gc

os.system("clear && clear")

#try:
if True:
    run_imagej = False

    dates = None
    if run_imagej:
        print "Running imagej..."
        methods.imagej()
        dates = methods.get_image_dates()
        methods.export_image_dates()
    else:
        dates = methods.import_image_dates()

    run_17 = True
    run_5 = True
    run_all = True
    create_videos = True
    clusters = 1
    avg_temp = 1
    order_param_exp = 1
    lost_percentage = 1

    coef = 5
    divisions = 30

    if run_17:
        print "K17"
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
        experiment_17 = methods.compress(experiment_17, level=9)
        rod_groups_17 = None
        gc.collect()
    except NameError:
        pass

    if run_5:
        print "K5"
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
        experiment_5 = methods.compress(experiment_5, level=9)
        rod_groups_5 = None
        gc.collect()
    except NameError:
        pass

    if run_all:
        print "K all"
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
        experiment_all = methods.compress(experiment_all, level=9)
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
