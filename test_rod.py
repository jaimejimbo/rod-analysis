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
            os.remove("prueba.txt")
            os.remove("prueba2.txt")
            os.remove("prueba3.txt")
        except:
            pass
        input1 = open('prueba.txt','w+')
        input1.write("Hola")
        input1.close()
        files = rod_statistics.import_files(_glob="*.txt", regular_expression='^prueba\.txt')
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
        os.remove("prueba.txt")
        os.remove("prueba2.txt")
        os.remove("prueba3.txt")


    def test_import_data(self):
        """
        Checks if data is parsed as expected.
        """
        pass

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
        pass

