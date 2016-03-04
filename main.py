#!/usr/bin/python -i

"""
Main script.
"""

"""
Settings
"""
#0-9 or None to disable
low_comp_level = None
medium_comp_level = 1
strong_comp_level = None
#None uses all cpus      
cpus = None
# 1 or True, 0 or False
create_videos = 1
clusters = 1
avg_temp = 1
order_param_exp = 1
lost_percentage = 1
run_imagej = 1
run_props = 0
run_graphs = 1
run_check_dim = 0
get_image_dates = 0
# variables
coef = 3
divisions = 30


import re







"""
Script
"""
def get_coords_from_imagej():
    """
        Gets center and radius from imagej script.
    """
    coords = []
    imagej_script = open("acquisition_v4.ijm", "r")
    reg_exp = re.compile(".*makeOval.*")
    for line in imagej_script:
        if reg_exp.match(line):
            data_line = line
            break
    zone_coords_re = re.compile(".*makeOval\((.*),(.*),(.*),(.*)\);.*")
    results_ = zone_coords_re.search(data_line)
    results = [results_.group(index+1) for index in range(4)]
    imagej_script.close()
    rad = float(results[2])/2
    center_x = float(results[1])+rad 
    center_y = float(results[0])+rad
    zone_coords_ = [center_x, center_y, rad]
    return zone_coords_

def set_coords_in_imagej(zone_coords_):
    """
        Sets values in imagej script.
    """
    imagej_script = open("acquisition_v4.ijm", "w")
    reg_exp = re.compile(".*makeOval.*")
    center_x = zone_coords[0]
    center_y = zone_coords[1]
    rad = zone_coords[2]
    values = "("+str(center_y-rad)+","+str(center_x-rad)+","+str(rad*2)+","+str(rad*2)+");"
    for line in imagej_script:
        if reg_exp.match(line):
            line = "\tmakeOval"+values
    imagej_script.close()

#zone_coords = []
#set_coords_in_imagej()
zone_coords = get_coords_from_imagej()
zone_coords = [zone_coords[1], zone_coords[0], zone_coords[2]]
#print zone_coords
#raw_input("")
settings = open("settings.py", "w")
settings.write("low_comp_level = {0}".format(low_comp_level))
settings.write("\n")
settings.write("medium_comp_level = {0}".format(medium_comp_level))
settings.write("\n")
settings.write("strong_comp_level = {0}".format(strong_comp_level))
settings.write("\n")
settings.write("cpus = {0}".format(cpus))
settings.write("\n")
settings.write("zone_coords = {0}".format(zone_coords))
settings.write("\n")
settings.close()

import experiment, system_state, methods, settings
import os, gc, math#, guppy

#hp = guppy.hpy()
#hp.setrelheap()

os.system("clear && clear")

#methods.export_image_dates()

if True:
    try:
        dates = None
            
        if get_image_dates:
            os.system("rm dates.txt")
            dates = methods.get_image_dates()
            methods.export_image_dates()
        else:
            dates = methods.import_image_dates()

        if run_imagej:
            print "Running imagej..."
            methods.imagej()


        def run_default(kappa, real_kappa, error):
            msg = "\t\t\tK"+str(real_kappa)+"\t\t\t"
            print msg
            names, rod_groups = system_state.create_rods(kappas=kappa, real_kappas=real_kappa,
                                                    allowed_kappa_error=error)
            experiment_ = experiment.Experiment(system_states_name_list=names, kappas=real_kappa,
                                                system_states_list=rod_groups,
                                                dates=dates, diff_t=5/3.0)
            experiment_.set_coef(coef)
            if create_videos:
                experiment_.create_videos(divisions=divisions, fps=10, max_distance=10, max_angle_diff=5,
                                         number_of_bursts=1)
            if clusters:
                experiment_.plot_cluster_areas(number_of_bursts=5, max_distance=100,
                        max_angle_diff=10, min_size=20)
            if avg_temp:
                experiment_.plot_average_temperature(100, 10)
            if lost_percentage:
                percentage, std_dev = experiment_.lost_rods_percentage
	        print "Rods lost: "+str(percentage)+"%"+" Standard deviation: "+str(std_dev)
            try:
                names = None
                experiment_ = None
                rod_groups = None
                gc.collect()
            except NameError:
                pass
            print "\n"

        def run_prop(kappa, real_kappa, error):
            msg = "\t\t\tK"+str(real_kappa)+"\t\t\t"
            names, rod_groups = system_state.create_rods(kappas=kappa, real_kappas=real_kappa,
                                                    allowed_kappa_error=error)
            experiment_ = experiment.Experiment(system_states_name_list=names, kappas=real_kappa,
                                                system_states_list=rod_groups,
                                                dates=dates, diff_t=5/3.0)
            experiment_.set_coef(coef)
            msg = str(experiment_.average_kappa) + "\t|\t" + str(experiment_.average_kappa_dev)
            rad = experiment_.average_system_rad
            print msg
            length = experiment_.average_rod_length
            width = experiment_.average_rod_width
            cov_area, cov_area_dev = experiment_.average_covered_area_proportion
            num_rods, num_rods_dev = experiment_.average_number_of_rods
            return cov_area, cov_area_dev, int(num_rods), int(num_rods_dev), rad
        

        if run_props:
            print "Computing area proportions..."
            print "kappa\t\t\tdeviation"
            prop_long, prop_long_dev, longs, longs_dev, rad1 = run_prop(15, 12, 3)
            prop_short, prop_short_dev, shorts, shorts_dev, rad2 = run_prop(7.8, 6, 2)
            areal = rad1**2*math.pi
            areas = rad2**2*math.pi
            rod_areal = float(areal*prop_long)/longs
            rod_areas = float(areas*prop_short)/shorts
            rad = float(rad1+rad2)/2
            area = math.pi*rad**2
            total_prop = prop_long + prop_short
            total_prop_err = prop_short_dev + prop_long_dev
            #print area
            print "Total area proportion: ", round(100*total_prop, 2), "% ", round(100*total_prop_err, 2), "%"
            print "Long proportion (area): ", round(100*float(prop_long)/total_prop, 2), "% ",round(100*prop_long_dev/total_prop, 2), "%"
            print "Number of long rods: ", longs, longs_dev
            print "Number of short rods: ", shorts, shorts_dev
            print "Total number of rods: ", longs+shorts
            waprop = input("Wanted area proportion: ")
            wlprop = input("Wanted longs proportion: ")
            Nl = int(wlprop*waprop*area/rod_areal)
            Ns = int(Nl*rod_areal*(1-wlprop)/(wlprop*rod_areas))
            print "Needed long rods: ", Nl
            print "Difference: ", Nl-longs
            print "Needed short rods: ", Ns
            print "Difference: ", Ns-shorts
            
        
        if run_graphs:
            print "Creating graphs..."
            run_default(15, 12, 3)
            run_default(7.8, 6, 2)

        if run_check_dim:
            names, rod_groups = system_state.create_rods(kappas=10, real_kappas=12,
                                                    allowed_kappa_error=10)
            experiment_ = experiment.Experiment(system_states_name_list=names, kappas=12,
                                                system_states_list=rod_groups,
                                                dates=dates, diff_t=5/3.0)
            experiment_.set_coef(coef)
            experiment_.plot_rods(0)


        #os.system("bash tomp4script.sh")
    except KeyboardInterrupt:
        os.system("sudo pkill -f main.py")
    except AssertionError as error:
        print "Assertion error({0})".format(error)
        os.system("sudo pkill -f main.py")
#print hp.heap()
print "Type exit() to exit."
def exit():
    print "Give root password to kill all remaining zombie processes"
    os.system("sudo pkill -f main.py")
