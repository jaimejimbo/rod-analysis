#!/usr/bin/python -i

"""
Main script.
"""

import experiment
import system_state
import methods

def main():
    """
    Main method.
    """
    run_imagej = False

    dates = None
    if run_imagej:
        methods.imagej()
        dates = methods.get_image_dates()
        methods.export_image_dates()
    else:
        dates = methods.import_image_dates()

    names, rod_groups_5 = system_state.create_rods(kappas=5.5,
                                            allowed_kappa_error=2)
    names, rod_groups_17 = system_state.create_rods(kappas=17,
                                            allowed_kappa_error=2)
    experiment_5 = experiment.Experiment(system_states_name_list=names,
                                        system_states_list=rod_groups_5,
                                        dates=dates, diff_t=5/3.0)
    experiment_17 = experiment.Experiment(system_states_name_list=names,
                                        system_states_list=rod_groups_17,
                                        dates=dates, diff_t=5/3.0)

    experiment_5.create_gifs(divisions=30)
    experiment_17.create_gifs(divisions=30)

main()
