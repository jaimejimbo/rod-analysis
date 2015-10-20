import unittest
import rod_statistics
import os

class TestRod(unittest.TestCase):
    """
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

    def test_Rod(self):
        """
        Checks rod object and methods.
        """
        pass

    def test_RodGroup(self):
        """
        Checks rod groups.
        """
        pass

    def test_CreateRods(self):
        """
        Checks rod groups.
        """
        rod_grps = rod_statistics.create_rods()
        print rod_grps
