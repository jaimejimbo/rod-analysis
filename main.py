#!/usr/bin/python -i

"""
Main script.
"""

"""
Settings
"""
#0-9 or None to disable
low_comp_level = 0
medium_comp_level = 0
strong_comp_level = 0 
#None uses all cpus      
cpus = None
# 1 or True, 0 or False
run_12 = 0
run_4 = 0
run_all = 0
create_videos = 0
clusters = 0
avg_temp = 0
order_param_exp = 0
lost_percentage = 0
run_imagej = 1
run_props = 1
# variables
coef = 3
divisions = 30










"""
Script
"""


settings = open("settings.py", "w")
settings.write("low_comp_level = {0}".format(low_comp_level))
settings.write("\n")
settings.write("medium_comp_level = {0}".format(medium_comp_level))
settings.write("\n")
settings.write("strong_comp_level = {0}".format(strong_comp_level))
settings.write("\n")
settings.write("cpus = {0}".format(cpus))
settings.write("\n")
settings.close()

import experiment, system_state, methods, settings
import os, gc, math

os.system("clear && clear")

try:
#if True:

    dates = None
    if run_imagej:
        print "Running imagej..."
        methods.imagej()
        dates = methods.get_image_dates()
        methods.export_image_dates()
    else:
        dates = methods.import_image_dates()
        

    def run_default(kappa, real_kappa, error):
        msg = "\t\t\tK"+str(real_kappa)+"\t\t\t"
        names, rod_groups = system_state.create_rods(kappas=kappa, real_kappas=real_kappa,
                                                allowed_kappa_error=error)
        experiment_ = experiment.Experiment(system_states_name_list=names, kappas=real_kappa,
                                            system_states_list=rod_groups,
                                            dates=dates, diff_t=5/3.0)
        experiment_.set_coef(coef)
        if create_videos:
            experiment_.divide_systems_in_circles(divisions=divisions)
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
        return cov_area, cov_area_dev, int(num_rods), int(num_rods_dev), rad, length, width
        
    if run_props:
        print "Computing area proportions..."
        print "kappa\t\t\tdeviation"
        prop_long, prop_long_dev, longs, longs_dev, rad1, lengthl, widthl = run_prop(15, 12, 3)
        prop_short, prop_short_dev, shorts, shorts_dev, rad2, lengths, widths = run_prop(7.8, 6, 2)
        rod_areal = lengthl*widthl
        rod_areas = lengths*widths
        rad = float(rad1+rad2)/2
        area = math.pi*rad**2
        total_prop = prop_long + prop_short
        #print area
        waprop = input("Wanted area proportion: ")
        wlprop = input("Wanted longs proportion: ")
        Nl = int(wlprop*waprop*area/rod_areal)
        Ns = int(Nl*rod_areal*(1-wlprop)/(wlprop*rod_areas))
        print "Total area proportion: ", total_prop
        print "Long proportion (area): ", float(prop_long)/total_prop, prop_long_dev/total_prop
        print "Number of long rods: ", longs, longs_dev
        print "Needed long rods: ", Nl
        print "Difference: ", Nl-longs
        print "Number of short rods: ", shorts, shorts_dev
        print "Needed short rods: ", Ns
        print "Difference: ", Ns-shorts
        print "Total number of rods: ", longs+shorts
        
    

    #os.system("bash tomp4script.sh")
except KeyboardInterrupt:
    os.system("sudo pkill -f main.py")
except AssertionError as error:
    print "Assertion error({0})".format(error)
    os.system("sudo pkill -f main.py")
print "Type exit() to exit."
def exit():
    print "Give root password to kill all remaining zombie processes"
    os.system("sudo pkill -f main.py")
