#!/usr/bin/python
from methods import *
import re
import os
import math
from methods import *
from experiment import *
from mpl_toolkits.mplot3d import Axes3D
import pylab
from scipy import interpolate
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt



#Gets rods_####.txt archives.
image_names = get_file_names(regular_expression=r'IMG_[0-9]{4}.JPG')
numbers = [get_number_from_string(image) for image in image_names]
start = min(numbers)
end = max(numbers)

imagej_template = open("./acquisition_v4.ijm", 'r')
imagej_script = open("./imagej_script.ijm", 'w')

current_folder = os.getcwd()

for line in imagej_template:
    output_line = line
    if re.match(r'.*for.img.*', line):
        output_line = "for(img_num="
        output_line += str(start) + "; img_num<="
        output_line += str(end)+ "; img_num++){\n"
    if re.match(r'.*\"\.\".*', line):
        output_line = "\tfolder = \""+str(current_folder)+"\";"
    imagej_script.write(output_line)

imagej_template.close()
imagej_script.close()

os.system("imagej ./imagej_script.ijm --headless")

#Imports data
kappas = [6, 17]
kappa_error = 2
names, rod_groups1 = create_rods(kappas=6, allowed_kappa_error=2)
names, rod_groups2 = create_rods(kappas=17, allowed_kappa_error=2)
experiment1 = Experiment(system_states_name_list=names, system_states_list=rod_groups1)
experiment2 = Experiment(system_states_name_list=names, system_states_list=rod_groups2)
print experiment1.evolution_dictionaries
