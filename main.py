#!/usr/bin/python -i

"""
Main script.
"""

import experiment, system_state, methods
import os

try:
    run_imagej = False

    dates = None
    if run_imagej:
        methods.imagej()
        dates = methods.get_image_dates()
        methods.export_image_dates()
    else:
        dates = methods.import_image_dates()

    #names, rod_groups_5 = system_state.create_rods(kappas=5.5,
    #                                        allowed_kappa_error=.5)
    names, rod_groups_17 = system_state.create_rods(kappas=17.5,
                                            allowed_kappa_error=1.5)
    #experiment_5 = experiment.Experiment(system_states_name_list=names,
    #                                    system_states_list=rod_groups_5,
    #                                    dates=dates, diff_t=5/3.0)
    experiment_17 = experiment.Experiment(system_states_name_list=names,
                                        system_states_list=rod_groups_17,
                                        dates=dates, diff_t=5/3.0)

    #experiment_5.set_coef(5)
    experiment_17.set_coef(5)
    #experiment_5.divide_systems_in_circles(divisions=30)
    #experiment_17.divide_systems_in_circles(divisions=30)
    #experiment_5.create_gifs(divisions=30, fps=10, max_distance=10, max_angle_diff=5,
    #                         number_of_bursts=1)
    experiment_17.create_gifs(divisions=30, fps=10, max_distance=10, max_angle_diff=5,
                             number_of_bursts=1)
    #cluster_areas5 = experiment_5.cluster_areas(max_distance=50,
    #                    max_angle_diff=30, min_size=3)
    #cluster_areas17 = experiment_17.cluster_areas(max_distance=50,
    #                    max_angle_diff=30, min_size=3)
    #print cluster_areas5
    #print cluster_areas17
    """x_val, y_val = experiment_5._states[1].rods_possitions
    import matplotlib.pyplot as plt
    plt.figure()
    plt.scatter(x_val, y_val)
    plt.gca().invert_yaxis()
    plt.show()"""

    os.system("bash tomp4script.sh")
except KeyboardInterrupt:
    os.system("sudo pkill -f main.py")
    raise KeyboardInterrupt
print "Give root password to kill all remaining zombie processes"
os.system("sudo pkill -f main.py")
