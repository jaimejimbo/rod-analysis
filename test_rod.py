import unittest
import rod_statistics
import os
import math

class TestRod(unittest.TestCase):
    """
    Group of tests for rod_statistics library
    """
    
    def SetUp(self):
        """
        Long message when debugging.
        """
        self.longMessage = True

    def test_import_files(self):
        """
        Checks if files are imported as expected.
        """
        try:
            os.remove("prueba1.txt")
            os.remove("prueba2.txt")
            os.remove("prueba3.txt")
        except:
            pass
        input1 = open('prueba1.txt','w+')
        input1.write("Hola")
        input1.close()
        files = rod_statistics.import_files(_glob="*.txt", regular_expression='^prueba1\.txt')
        line1 = files[0].readline()
        self.assertEqual(line1.split(), ["Hola"], "Error importing singular file. Obtained: " + str(line1.split()))
        input2 = open('prueba2.txt','w+')
        input2.write("Adios\n")
        input2.close()
        files = rod_statistics.import_files(_glob="*.txt", regular_expression='^prueba*')
        line1 = files[0].readline()
        line2 = files[1].readline()
        self.assertEqual(line1.split(), ["Adios"], "Error importing 2 files. Obtained: " + str(line1.split()))
        self.assertEqual(line2.split(), ["Hola"], "Error importing 2 files. Obtained: " + str(line2.split()))
        input3 = open('prueba3.txt','w+')
        input3.write("Hasta pronto\n")
        input3.close()
        files = rod_statistics.import_files(_glob="*.txt", regular_expression='^prueba*')
        line1 = files[0].readline()
        line2 = files[1].readline()
        line3 = files[2].readline()
        self.assertEqual(line1.split(), ["Adios"], "Error importing 3 files. Obtained: " + str(line1.split()))
        self.assertEqual(line3.split(), ["Hola"], "Error importing 3 files. Obtained: " + str(line3.split()))
        self.assertEqual(line2.split(), ["Hasta","pronto"], "Error importing 3 files. Obtained: " + str(line2.split()))
        os.remove("prueba1.txt")
        os.remove("prueba2.txt")
        os.remove("prueba3.txt")


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
        rod_grps = rod_statistics.create_rods()
        rod = rod_grps[0].get_rod()
        self.assertEqual(str(rod.id), "1", "Errors when importing data. ID obtained: "+str(rod.id))
        self.assertEqual(str(rod.area), "578.0", "Errors when importing data. Area obtained: "+str(rod.area))
        rod = rod_grps[1].get_rod()
        self.assertEqual(str(rod.id), "1", "Errors when importing data. ID obtained: "+str(rod.id))
        self.assertEqual(str(rod.area), "578.0", "Errors when importing data. Area obtained: "+str(rod.area))

    def test_Rod(self):
        """
        Checks rod object and methods.
        """
        pass

    def test_SubsystemState(self):
        """
        Checks rod groups.
        """
        pass

    def test_SystemState(self):
        """
        Checks rod groups.
        """
        pass

    def test_segment_area(self):
        """
        Tests for segment area method
        """
        half_circle_area = math.pi/2
        computed_area = rod_statistics.segment_area(1,0)
        self.assertEqual(computed_area, half_circle_area, msg="Error in sement area computing. Obtained: "+str(computed_area)+" Expected: "+str(half_circle_area))
        computed_area = rod_statistics.segment_area(1,0.5)
        self.assertAlmostEqual(computed_area, half_circle_area/2, msg="Error in sement area computing. Obtained: "+str(computed_area)+" Expected: "+str(half_circle_area/2), delta=half_circle_area*.05)
        computed_area = rod_statistics.segment_area(1e10,1e10-1)
        self.assertAlmostEqual(computed_area, 0, msg="Error in sement area computing. Obtained: "+str(computed_area)+" Expected: "+str(0), delta=10)
        

    def test_effective_area(self):
        """
        Tests for effective area method
        """
        computed_area = rod_statistics.effective_area(1,99.5,100)
        half_circle_area = math.pi
        self.assertAlmostEqual(computed_area, half_circle_area/2, msg="Error in sement area computing. Obtained: "+str(computed_area)+" Expected: "+str(half_circle_area/2), delta=half_circle_area*.05)
        

