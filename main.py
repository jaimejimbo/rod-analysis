#!/usr/bin/python -i

"""
Main script.
"""

import experiment, system_state, methods
import os, gc

try:
    run_imagej = False

    dates = None
    if run_imagej:
        methods.imagej()
        dates = methods.get_image_dates()
        methods.export_image_dates()
    else:
        dates = methods.import_image_dates()

    run_17 = False
    run_5 = True
    run_all = False
    create_gifs = True
    clusters = False
    avg_temp = False
    order_param_exp = False
    lost_percentage = False

    coef = 5
    divisions = 30

    if run_17:
        names, rod_groups_17 = system_state.create_rods(kappas=17.5,
                                                allowed_kappa_error=1.5)
        experiment_17 = experiment.Experiment(system_states_name_list=names,
                                            system_states_list=rod_groups_17,
                                            dates=dates, diff_t=5/3.0)
        experiment_17.set_coef(coef)
        
        if create_gifs:
            experiment_17.divide_systems_in_circles(divisions=divisions)
            experiment_17.create_gifs(divisions=divisions, fps=10, max_distance=10, max_angle_diff=5,
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
        del names, experiment_17, rod_groups_17
        gc.collect()
    except NameError:
        pass

    if run_5:
        names, rod_groups_5 = system_state.create_rods(kappas=5.5,
                                                allowed_kappa_error=.5)
        experiment_5 = experiment.Experiment(system_states_name_list=names,
                                            system_states_list=rod_groups_5,
                                            dates=dates, diff_t=5/3.0)
        experiment_5.set_coef(coef)
        if create_gifs:
            experiment_5.divide_systems_in_circles(divisions=divisions)
            experiment_5.create_gifs(divisions=divisions, fps=10, max_distance=10, max_angle_diff=5,
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
        del names, experiment_5, rod_groups_5
        gc.collect()
    except NameError:
        pass

    if run_all:
        names, rod_groups_all = system_state.create_rods(kappas=10,
                                                allowed_kappa_error=10)
        experiment_all = experiment.Experiment(system_states_name_list=names,
                                            system_states_list=rod_groups_all,
                                            dates=dates, diff_t=5/3.0)
        experiment_all.set_coef(coef)
        if create_gifs:
            experiment_all.divide_systems_in_circles(divisions=divisions)
            experiment_all.create_gifs(divisions=divisions, fps=10, max_distance=10, max_angle_diff=5,
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
        del names, experiment_all, rod_groups_all
        gc.collect()
    except NameError:
        pass

    #os.system("bash tomp4script.sh")
except KeyboardInterrupt:
    os.system("sudo pkill -f main.py")
except AssertionError:
    os.system("sudo pkill -f main.py")
print "Type exit() to exit."
def exit():
    print "Give root password to kill all remaining zombie processes"
    os.system("sudo pkill -f main.py")
