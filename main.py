#!/usr/bin/python -i

from experiment import *
from base_classes import *

imagej()
kappas = [5.5, 17]
kappa_error = 1
names, rod_groupsK5 = create_rods(kappas=5.5, allowed_kappa_error=2)
names, rod_groupsK17 = create_rods(kappas=17, allowed_kappa_error=2)
experimentK5 = Experiment(system_states_name_list=names, system_states_list=rod_groupsK5)
experimentK17 = Experiment(system_states_name_list=names, system_states_list=rod_groupsK17)

experimentK5.create_gifs(rad=50)

