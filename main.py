#!/usr/bin/python -i

while True:
    try:
        import math
        from experiment import *
        from base_classes import *
        from mpl_toolkits.mplot3d import Axes3D
        import pylab
        from scipy import interpolate
        from matplotlib import cm
        from matplotlib.ticker import LinearLocator, FormatStrFormatter
        import matplotlib.pyplot as plt
        import os
        break
    except:
        import os
        print "Installing dependancies (Ubuntu 14.04 or simmilar)"
        os.system("sudo apt-get install python-dev python-pip libfreetype6-dev python-matplotlib python-numpy python-scipy")


#Gets rods_####.txt archives.
image_names = get_file_names(regular_expression=r'IMG_[0-9]{4}.JPG')
numbers = [get_number_from_string(image) for image in image_names]
try:
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

    run_imagej = raw_input("Run imagej script?(Yn)")
    if run_imagej == "y" or not run_imagej:
        os.system("imagej ./imagej_script.ijm")
except ValueError:
    pass

#Imports data
kappas = [5.5, 17]
kappa_error = 1
names, rod_groupsK5 = create_rods(kappas=5.5, allowed_kappa_error=2)
names, rod_groupsK17 = create_rods(kappas=17, allowed_kappa_error=2)
experimentK5 = Experiment(system_states_name_list=names, system_states_list=rod_groupsK5)
experimentK17 = Experiment(system_states_name_list=names, system_states_list=rod_groupsK17)
experimentK5.compute_dictionaries()
experimentK17.compute_dictionaries()






