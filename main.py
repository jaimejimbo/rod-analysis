#!/usr/bin/python -i

"""
Main script.
"""

import experiment, system_state, methods
import os

#try:
run_imagej = False

dates = None
if run_imagej:
    methods.imagej()
    dates = methods.get_image_dates()
    methods.export_image_dates()
else:
    dates = methods.import_image_dates()

run_17 = True
run_5 = False
run_all = False
create_gifs = False
clusters = True

if run_17:
    names, rod_groups_17 = system_state.create_rods(kappas=17.5,
                                            allowed_kappa_error=1.5)
    experiment_17 = experiment.Experiment(system_states_name_list=names,
                                        system_states_list=rod_groups_17,
                                        dates=dates, diff_t=5/3.0)
    experiment_17.set_coef(5)
    
    if create_gifs:
        experiment_17.divide_systems_in_circles(divisions=30)
        experiment_17.create_gifs(divisions=30, fps=10, max_distance=10, max_angle_diff=5,
                                 number_of_bursts=1)
    if clusters:
        #try:
        #cluster_areas17 = experiment_17.cluster_areas(max_distance=50,
        #                    max_angle_diff=30, min_size=3)
        for state in experiment_17._states:
            print state.total_area_of_clusters(max_distance=100, max_angle_diff=30)
        """except:
            print "Exception in user code:"
            print '-'*60
            traceback.print_exc(file=sys.stdout)
            print '-'*60
            print e.errno, e.strerror"""

if run_5:
    names, rod_groups_5 = system_state.create_rods(kappas=5.5,
                                            allowed_kappa_error=.5)
    experiment_5 = experiment.Experiment(system_states_name_list=names,
                                        system_states_list=rod_groups_5,
                                        dates=dates, diff_t=5/3.0)
    experiment_5.set_coef(5)
    if create_gifs:
        experiment_5.divide_systems_in_circles(divisions=30)
        experiment_5.create_gifs(divisions=30, fps=10, max_distance=10, max_angle_diff=5,
                                 number_of_bursts=1)
    if clusters:
        cluster_areas5 = experiment_5.cluster_areas(max_distance=50,
                            max_angle_diff=30, min_size=3)

if run_all:
    names, rod_groups_all = system_state.create_rods(kappas=10,
                                            allowed_kappa_error=10)
    experiment_all = experiment.Experiment(system_states_name_list=names,
                                        system_states_list=rod_groups_all,
                                        dates=dates, diff_t=5/3.0)
    experiment_all.set_coef(5)
    if create_gifs:
        experiment_all.divide_systems_in_circles(divisions=30)
        experiment_all.create_gifs(divisions=30, fps=10, max_distance=10, max_angle_diff=5,
                                 number_of_bursts=1)
    if clusters:
        cluster_areas_all = experiment_all.cluster_areas(max_distance=50,
                            max_angle_diff=30, min_size=3)

os.system("bash tomp4script.sh")
#except:
#    os.system("sudo pkill -f main.py")
print "Give root password to kill all remaining zombie processes"
os.system("sudo pkill -f main.py")
