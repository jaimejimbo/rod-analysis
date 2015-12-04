import math
import random

class Matrix:
    """
    Create a matrix
    """
    
    ###########
    # Methods #
    ###########
    
    def __init__(self, data = [[]]):
        """
        Initialization method Matrix()
        """
        if data == [[]]:
            self._rows = 0
            self._cols = 0
            self._data = data
        else:
            self._rows = len(data)
            self._cols = len(data[0])
            self._data = data
         
    def __str__(self):
        """
        Conversion to string str()
        """
        #output = "["
        output = "\n"
        for row in range(self._rows):
            #output += "["
            for col in range(self._cols):
                output += str(self._data[row][col])
                if col < self._cols-1:
                    output += "\t"
            #output += "]"
            if row < self._rows-1:
                output += "\n"
        #output += "]"
        output += "\n"
        return output
        
        
    def __getitem__(self, pos):
        """
        called when element[][]
        """
        return self._data[pos]
        
    def __len__(self):
        """
        Method for len()
        """
        return len(self._data)
        
    def __eq__(self, matrix):
        """
        Method for ==
        """
        for row in range(self._rows):
            for col in range(self._cols):
                if self._data[row][col] != matrix[row][col]:
                    return False
        return True
        
    def __ne__(self, matrix):
        """
        Method for !=
        """
        return not(self == matrix)
    
    
    def __matrix__(self, matrix):
        print "In da hood"
        
    def __trunc__(self, matrix_vector):
        """
        Matrix(vector) -> matrix
        """
        rows = len(matrix_vector)+1
        cols = len(matrix_vector[0])+1
        output = Matrix(matrix_vector)
        return output
        
    def __pow__(self, exponent):
        """
        A^n
        """
        A = self.copy()
        B = self.copy()
        for dummy_index in range(exponent-1):
            B = B * A
        return B
        
        
    def __add__(self,b):
        """
        Sum of two matrixes (same dimensions) A+B
        """
        adim = self.get_dim()
        bdim = b.get_dim()
        
        try:
            if adim != bdim:
                raise DimError
            else:
                output = zeros(adim[0], adim[1])
                for row in range(adim[0]):
                    for col in range(adim[1]):
                        avalue = self._data[row][col]
                        bvalue = b[row][col]
                        output[row][col] = avalue + bvalue
            return output
        except:
            print "NotAMatrix"
        
    def __radd__(self,b):
        """
        Right addition
        """
        return self + b
        
    def __sub__(self,b):
        """
        Substraction        
        """
        return self + (-1)*b
        
    def __rsub__(self,b):
        """
        Right substraction
        """
        return (-1)*(self-b)
    
    def __mul__(self,input_):
        """
        Product of 2 matrixes A*B
        """
        self_rows = self.get_rows()
        self_cols = self.get_cols()
        try:
            input_rows = input_.get_rows()
            input_cols = input_.get_cols()
            
            if self_cols != input_rows:
                raise DimError
            else:
                output = zeros(self_rows, input_cols)
                for self_row in range(self_rows):
                    for input_col in range(input_cols):
                        for index in range(self_cols):
                            output[self_row][input_col] += self._data[self_row][index] * input_[index][input_col]
            return output
        except:
            #Not a matrix
            output = zeros(self_rows, self_cols)
            for row in range(self_rows):
                for col in range(self_cols):
                    output[row][col] = input_ * 1.0 * self._data[row][col]
            return output
            
    def __rmul__(self,input_):
        """
        Right multiplication B*A
        """
        return self*input_ 
        
    def __div__(self,input_):
        """
        Division A/B
        """
        try:
            return self*input_.inverse()
        except AttributeError:
            new_num = 1.0/input_
            return self*new_num
            
    def __rdiv__(self,input_):
        """
        Right division B/A
        """
        return input_*self.inverse()
        
    ### TESTING    
    def __iter__(self):
        """
        for element in matrix -- returns row,col,element
        it ables someone to do:
            for element in matrix:
                matrix[element[0]][element[1]] = function(element[2])
        """
        for row in range(self.get_rows()):
            for col in range(self.get_cols()):
                yield row, col, self[row][col]
                
    def __call__(self,input1,input2):
        """
        called when A(row,col) (A(1,1) -> A[0][0])
        """
        return self._data[input1-1][input2-1]
        
    def to_float(self):
        """
        To float conversion
        """
        output = self.copy()
        for row in range(self.get_rows()):
            for col in range(self.get_cols()):
                output[row][col] = float(self[row][col])
        return output
    
    def get_rows(self):
        """
        Getter for rows
        """
        return self._rows
        
    def get_cols(self):
        """
        Getter for cols
        """
        return self._cols
        
    def get_dim(self):
        """
        Returns (rows,cols)
        """
        return self._rows, self._cols 
            
    def det(self, matrix):
        """
        Returns the determinant
        """
        matrix_dim = matrix.get_dim()
        if matrix_dim[0] != matrix_dim[1]:
            raise NotSquaredError
        _det = 0
        if matrix_dim[0] == 1:
            return matrix[0][0]
        for index in range(matrix_dim[0]):
            _det += (-1)**index * matrix[index][0] * self.det(matrix.delete(index,0))
        return _det
        
    def delete(self,row,col):
        """
        Delete an row and a col and returns the result
        """
        rows = self.get_rows()-1
        cols = self.get_cols()-1
        output = zeros(rows, cols)
        for row_ in range(rows):
            for col_ in range(cols):
                adding_row = int(row_ >= row)
                adding_col = int(col_ >= col)
                row__ = row_ + adding_row
                col__ = col_ + adding_col
                output[row_][col_] = self._data[row__][col__]
        return output
                    
    def inverse(self):
        """
        Inverse of the matrix
        """
        try:
            return self.adj().traspose()/self.det(self)
        except ZeroDivisionError:
            print "Not invertible matrix (null determinant)!"
            print self
            raise NotInvertibleError
        except NotSquaredError:
            print "Not squared matrix (not squared)!"
            print self
            raise NotInvertibleError
            
    def identity(self):
        """
        Returns identity element in actual space
        """
        cols = self.get_cols()
        output = zeros(cols, cols)
        for pos in range(cols):
            output[pos][pos] = 1.0
        return output

        
    def copy(self):
        """
        Returns a copy of the matrix
        """
        rows = self.get_rows()
        cols = self.get_cols()
        output=zeros(rows, cols)
        for row in range(rows):
            for col in range(cols):
                output[row][col] = self._data[row][col]
        return output
        
    def traspose(self):
        """
        returns the trasposed of the matrix
        """
        output = zeros(self._cols, self._rows)
        for row in range(self._rows):
            for col in range(self._cols):
                output[col][row] = self._data[row][col]
        return output
        
    def adj(self):
        """
        returns de adjunt of the matrix
        """
        output = self.copy()
        rows = self._rows
        cols = self._cols
        for row in range(rows):
            for col in range(cols):
                output[row][col] = (-1)**(row+col) * self.det(self.delete(row,col))
        return output

    def diagonalize_2x2(self):
        """
        Returns 2 eigenvalues.
        |a b|
        |c d|
        eigen eq : lambda^2 - (a+d)*lambda + a*d-c*b
        """
        coef1 = -(self[0][0]+self[1][1])
        coef2 = self.det(self)
        sqrt = math.sqrt(coef1**2-4*coef2)
        lambda1 = float(coef+sqrt)/2
        lambda2 = float(coef-sqrt)/2
        return lambda1, lambda2
        

    @property
    def data(self):
        """
        Allows to transform a matrix to python vectors representation.
        """
        return self._data
            
def zeros(rows = 1, cols = 1):
    """
    Create a zeros matrix
    """
    A = Matrix([[0.0 for dummy_index in range(cols)] for dummy_index2 in range(rows)])
    return A
    
def random(rows = 1, cols = 1, min_ = 0, max_ = 10):
    """
    Create a random matrix
    """
    A = Matrix([[random.ranint(min_,max_) for dummy_index in range(cols)] for dummy_index2 in range(rows)])
    return A
            
class DimError(Exception):
    pass

class NotSquaredError(Exception):
    pass

class NotInvertibleError(Exception):
    pass
