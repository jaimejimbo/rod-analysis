import unittest
import rod_statistics
import os
import math
from methods import *
from experiment import *
from mpl_toolkits.mplot3d import Axes3D
import pylab
import numpy
from scipy import interpolate
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

class TestRod(unittest.TestCase):
    """
    Group of tests for rod_statistics library
    """
    
    def SetUp(self):
        """
        Long message when debugging.
        """
        self.longMessage = True

    def test_import_data(self):
        """
        Checks if data is parsed as expected.
        """
        try:
            os.remove("prueba1.txt")
            os.remove("prueba2.txt")
            os.remove("prueba3.txt")
        except:
            pass
        input3 = open('prueba3.txt','w+')
        input3.write("Hasta pronto")
        input3.close()
        input3 = open('prueba3.txt','r')
        data = rod_statistics.import_data(input3, " ", '[a-zA-Z]*?')
        self.assertEqual(data, [["Hasta","pronto"]],"Import data failed with strings and \" \" separator. Obtained: "+str(data))
        input3.close()
        input2 = open('prueba2.txt','w+')
        input2.write("hola-hasta-pronto")
        input2.close()
        input2 = open('prueba2.txt','r')
        data = rod_statistics.import_data(input2, "-", '[a-zA-Z]*?')
        self.assertEqual(data, [["hola","hasta","pronto"]],"Import data failed with strings and \"-\" separator. Obtained: "+str(data))        
        input2.close()
        input2 = open('prueba2.txt','w+')
        input2.write("hola\thasta\tpronto\nnueva\tlinea")
        input2.close()
        input2 = open('prueba2.txt','r')
        data = rod_statistics.import_data(input2, "\t", '[a-zA-Z]*?')
        self.assertEqual(data, [["hola","hasta","pronto"],["nueva","linea"]],"Import data failed with strings and \"\\t\" separator. Obtained: "+str(data))        
        input2.close()
        input1 = open('prueba1.txt','w+')
        input1.write("0.1123\t0.12312\t123123.123\n5.4\t323.123\t213.555\n5\t4\t3")
        input1.close()
        input1 = open('prueba1.txt','r')
        data = rod_statistics.import_data(input1, "\t", '[0-9]*?\.?[0-9]*')
        self.assertEqual(data, [["0.1123","0.12312","123123.123"],["5.4","323.123","213.555"],["5","4","3"]],"Import data failed with strings and \"\\t\" separator. Obtained: "+str(data))        
        input1.close()
        os.remove("prueba1.txt")
        os.remove("prueba2.txt")
        os.remove("prueba3.txt")

    def test_CreateRods(self):
        """
        Checks rod group creation
        """
        pass

    def test_segment_area(self):
        """
        Tests for segment area method
        """
        total_area = math.pi
        half_circle_area = total_area/2
        computed_area = segment_area(1,0)
        self.assertEqual(computed_area, half_circle_area, msg="Error in segment area computing #1. Obtained: "+str(computed_area)+" Expected: "+str(half_circle_area))
        computed_area = segment_area(1,0.5)
        self.assertAlmostEqual(computed_area, half_circle_area/2, msg="Error in segment area computing #2. Obtained: "+str(computed_area)+" Expected: "+str(half_circle_area/2), delta=half_circle_area*.2)
        computed_area = segment_area(1,0.999999999)
        self.assertAlmostEqual(computed_area, 0, msg="Error in segment area computing #3. Obtained: "+str(computed_area)+" Expected: "+str(0), delta=0.0005)
        computed_area = segment_area(1e6,1e6-.1)
        self.assertAlmostEqual(computed_area, 0, msg="Error in segment area computing #4. Obtained: "+str(computed_area)+" Expected: "+str(0), delta=((1e-4)*math.pi*(1e6)**2))
        computed_area = segment_area(1e6,1e6-1e-10)
        self.assertAlmostEqual(computed_area, 0, msg="Error in segment area computing #5. Obtained: "+str(computed_area)+" Expected: "+str(0), delta=((1e-8)*math.pi*(1e6)**2))
        computed_area = segment_area(1,-0.1)
        self.assertAlmostEqual(computed_area, half_circle_area, msg="Error in segment area computing #6. Obtained: "+str(computed_area)+" Expected: "+str(half_circle_area), delta=half_circle_area*.2)
        computed_area = segment_area(1,-0.5)
        expected = total_area*3.0/4
        self.assertAlmostEqual(computed_area, expected, msg="Error in segment area computing #7. Obtained: "+str(computed_area)+" Expected: "+str(expected), delta=(expected)*.1)
        computed_area = segment_area(1,-0.99999)
        self.assertAlmostEqual(computed_area, total_area, msg="Error in segment area computing #8. Obtained: "+str(computed_area)+" Expected: "+str(0), delta=0.01)

    def test_compute_min_dist(self):
        """
        Tests for h computation method
        """
        computed = compute_min_dist(1,1e10-.5, 1e10)
        expected = .5
        self.assertAlmostEqual(computed,expected,delta=.1,msg="Error in compute_min_dist #1: Obtained: "+str(computed)+" Expected: "+str(expected))
        computed = compute_min_dist(1,1e10, 1e10)
        expected = 0
        self.assertAlmostEqual(computed,expected,delta=.1,msg="Error in compute_min_dist #2: Obtained: "+str(computed)+" Expected: "+str(expected))
        computed = compute_min_dist(1,1e10+.5, 1e10)
        expected = -.5
        self.assertAlmostEqual(computed,expected,delta=.1,msg="Error in compute_min_dist #3: Obtained: "+str(computed)+" Expected: "+str(expected))
        

    def test_effective_area(self):
        """
        Tests for effective area method
        """
        total_area = math.pi
        half_circle_area = total_area/2
        computed_area = effective_area(1,1e10,1e10)
        self.assertAlmostEqual(computed_area, half_circle_area, msg="Error in effective_area computing #1. Obtained: "+str(computed_area)+" Expected: "+str(half_circle_area), delta=half_circle_area*3e-2)
        computed_area = effective_area(1,10,10)
        self.assertAlmostEqual(computed_area, half_circle_area, msg="Error in effective_area computing #2. Obtained: "+str(computed_area)+" Expected: "+str(half_circle_area), delta=half_circle_area*3e-2)

    def test_same_area_rad(self):
        """
        Tests for same area radius method 
        """
        small_rad, small_position_rad, main_rad = 1,9,10
        rad = same_area_rad(small_rad, small_position_rad, main_rad)
        self.assertAlmostEqual(rad, small_rad, delta=small_rad*0.05, msg="Error in same_area_rad #1. Initial rad should be returned. Obtained: "+str(rad))
        small_rad, small_position_rad, main_rad = 1,9.5,10
        rad = same_area_rad(small_rad, small_position_rad, main_rad)
        expected = small_rad*1.5
        self.assertAlmostEqual(rad, expected, delta=expected*.3, msg="Error in same_area_rad #2. Obtained: "+str(rad)+" Expected: "+str(expected))
        obtained = effective_area(rad,9.5,10)
        expected = math.pi
        self.assertAlmostEqual(obtained,expected,delta=math.pi*0.1,msg="Areas should be equivalent: Obtained: "+str(obtained)+" Expected: "+str(expected))
        small_rad, small_position_rad, main_rad = 1,1e6-.3, 1e6
        obtained = same_area_rad(small_rad, small_position_rad, main_rad)  
        expected = 1.5
        self.assertAlmostEqual(obtained, expected, delta=expected*.3, msg="Error in same_area_rad #3. Obtained: "+str(obtained)+" Expected: "+str(expected))

    def test_binary_search(self):
        """
        Checks binary search method.
        """
        a = [0,1,2,3,4,5,6,7,8,9,10]
        def extract(index):
            return a[int(index)]
        expected = 5
        obtained = int(binary_search(0,10,extract,expected,.1,10))
        self.assertEqual(obtained, expected, "Error in binary search #1 Obtained: "+str(obtained)+" Expected:"+str(expected))
        expected = 7
        obtained = int(binary_search(0,10,extract,expected,.1,10))
        self.assertEqual(obtained, expected, "Error in binary search #2 Obtained: "+str(obtained)+" Expected:"+str(expected))
        expected = 8
        obtained = int(binary_search(0,10,extract,expected,.1,10))
        self.assertEqual(obtained, expected, "Error in binary search #3 Obtained: "+str(obtained)+" Expected:"+str(expected))

    def test_binary_order(self):
        """
        Tests binary order method.
        """
        def ordering_id(element):
            return element
        a = [9,8,7,6,5,4,3,2,1,0,-1,-2]
        ordered_a = rod_statistics.binary_order(a, ordering_id)
        self.assertEqual(str(ordered_a),"[-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]","Error in ordering function #1 Obtained: "+str(ordered_a)+" Expected: [-2,-1,0,1,2,3,4,5,6,7,8,9]")

    def test_SystemState(self):
        """
        Checks rod groups and rod class.
        """
        names, rod_groups = rod_statistics.create_rods(folder="../rod-analysis", kappas=20, allowed_kappa_error=20, radius_correction_ratio=.1)
        original_rod_num = rod_groups[0].number_of_rods
        rod_group = rod_groups[0].clone
        clone_rod_num = rod_group.number_of_rods
        num_of_rods = rod_groups[1].number_of_rods
        rg = rod_groups[1]
        rod_group = rod_groups[1]
        average_length = rod_group.average_cluster_rod_num(30,10)
        expected = 2
        msg = "length must be greater than 2, obtained: "+str(average_length)
        msg += "\n"
        msg += str(rod_group.clusters()) + "\n" + str(average_length)
        self.assertGreaterEqual(average_length, expected, msg)
        names, rod_groups = rod_statistics.create_rods(folder="../rod-analysis", kappas=20, allowed_kappa_error=20, radius_correction_ratio=.1)
        rod_group = rod_groups[0].clone
        original_rod_num = rod_groups[0].number_of_rods
        clone_rod_num = rod_group.number_of_rods
        msg = "Original must not change when changing clone! Computed:" + str(clone_rod_num)+ " Expected:" + str(original_rod_num-1)
        self.assertEqual(clone_rod_num, original_rod_num,msg)
        rod = rod_group.get_rod()
        clone_rod_num = rod_group.number_of_rods
        msg = "Seems that rod are not correctly poped. Computed:" + str(clone_rod_num)+ " Expected:" + str(original_rod_num-1)
        self.assertEqual(clone_rod_num, original_rod_num-1,msg)
        """rad = 50
        x, y, z = rod_group.plottable_density_matrix(rad)
        plt.figure(1)
        plt.scatter(x,y,c=z,marker='s',s=(rad/2.0)**2)
        plt.colorbar()
        name = names[1]
        name += "all.png"
        plt.savefig(name)
        x, y, z = rod_group.relative_g2_plot_matrix(rad)
        plt.figure(4)
        plt.scatter(x,y,c=z,marker='s',s=(rad/2.0)**2)
        plt.colorbar()
        name = names[1]
        name += "all_g2.png"
        plt.savefig(name)
        names, rod_groups = rod_statistics.create_rods(folder="../rod-analysis", kappas=5.5, allowed_kappa_error=0.3, radius_correction_ratio=.1)
        rod_group = rod_groups[1]
        x, y, z = rod_group.plottable_density_matrix(rad=rad)
        plt.figure(2)
        plt.scatter(x,y,c=z,marker='s',s=(rad/2.0)**2)
        plt.colorbar()
        name = names[1]
        name += "K5.png"
        plt.savefig(name)
        names, rod_groups = rod_statistics.create_rods(folder="../rod-analysis", kappas=17, allowed_kappa_error=2, radius_correction_ratio=.1)
        rod_group = rod_groups[1] 
        x, y, z = rod_group.plottable_density_matrix(rad=rad)
        plt.figure(3)
        plt.scatter(x,y,c=z,marker='s',s=(rad/2.0)**2)
        plt.colorbar()
        name = names[1]
        name += "K17.png"
        plt.savefig(name)"""

    def test_Experiment(self):
        """
        Test experiment library.
        """
        names, rod_groups = rod_statistics.create_rods(folder="../rod-analysis", kappas=20, allowed_kappa_error=20, radius_correction_ratio=.1)
        experiment = Experiment(rod_groups)
        evo_dicts = experiment.evolution_dictionaries(20,5)
        print evo_dicts

