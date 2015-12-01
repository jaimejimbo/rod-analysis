#!/usr/bin/python -i

"""
Main script.
"""

import experiment
import system_state
import methods

run_imagej = False

dates = None
if run_imagej:
    methods.imagej()
    dates = methods.get_image_dates()
    methods.export_image_dates()
else:
    dates = methods.import_image_dates()

names, rod_groups_5 = system_state.create_rods(kappas=5.5,
                                        allowed_kappa_error=.5)
names, rod_groups_17 = system_state.create_rods(kappas=17,
                                        allowed_kappa_error=.5)
experiment_5 = experiment.Experiment(system_states_name_list=names,
                                    system_states_list=rod_groups_5,
                                    dates=dates, diff_t=5/3.0)
experiment_17 = experiment.Experiment(system_states_name_list=names,
                                    system_states_list=rod_groups_17,
                                    dates=dates, diff_t=5/3.0)

#print experiment_5._states[0].plottable_density_matrix(divisions=10)
experiment_5.create_gifs(divisions=10, fps=1, max_distance=10, max_angle_diff=5)
experiment_17.create_gifs(divisions=10, fps=1, max_distance=10, max_angle_diff=5)
#print experiment_5.cluster_area(max_distance=10, max_angle_diff=2)
#print experiment_17.cluster_area(max_distance=10, max_angle_diff=2)
"""x_val, y_val = experiment_5._states[1].rods_possitions
import matplotlib.pyplot as plt
plt.figure()
plt.scatter(x_val, y_val)
plt.gca().invert_yaxis()
plt.show()"""
