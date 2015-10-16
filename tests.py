import unittest
import matrix
import rod_statistics

class TestMatrix(unittest.TestCase):
    """
    Test Matrix class    
    """

    def setUp(self):
        pass
    
    def test_getset(self):
        """
        Getter and Setter tests. 
        """
        a = matrix.zeros(1,1)
        a[0][0] = 1
        self.assertEqual(a[0][0], 1, "Setter not working.")

    def test_str(self):
        """
        Getter and Setter tests. 
        """
        a = matrix.zeros(3,3)
        self.assertEqual(a.__str__(), "\n0.0\t0.0\t0.0\n0.0\t0.0\t0.0\n0.0\t0.0\t0.0\n", "str test.")
    
    def test_sumprod(self):
        """
        Sum and prod tests.
        """
        b = matrix.zeros(3,3)
        b[1][1] = 1.0
        a = matrix.zeros(3,3)
        self.assertEqual(str(a+b), str(b))
        self.assertEqual(str(a*b), str(a))
      
    def test_putout(self):  
        """
        ...
        """
        b = matrix.Matrix([[1,2,3],[4,5,6],[7,8,9]])
        c = b.delete(1,1)
        d = matrix.Matrix([[1,3],[7,9]])
        self.assertEqual(str(c), str(d))
    
    def test_det(self):
        """
        Determinant tests
        """
        a = matrix.Matrix([[1,1,1],[1,1,1],[1,1,1]])
        self.assertEqual(str(a.det(a)), str(0), "Error in det")
    
    def test_inverse(self):
        """
        Inversion tests
        """
        a = matrix.Matrix([[1,0,1],[0,1,0],[1,0,-1]])
        ident = a.identity()
        inv = a.inverse()
        self.assertEqual(str(a*inv),str(ident), "Error in inverse")
    
    def test_strmatrix(self):
        """
        Generalization of matrix values
        """
        a = matrix.Matrix([["Hola"]])




