import matrix
import poc_simpletest

def test1():
    test_obj = poc_simpletest.TestSuite()
    
    a = matrix.Matrix(1,1)
    a[0][0] = 1
    test_obj.run_test(a[0][0], 1, "Error in getter / setter")
    
    a = matrix.Matrix(3,3)
    test_obj.run_test(a.__str__(), "0.0\t0.0\t0.0\n0.0\t0.0\t0.0\n0.0\t0.0\t0.0", "Error in getter / setter")
    
    b = matrix.Matrix(3,3)
    b[1][1] = 1.0
    
    test_obj.run_test(str(a+b), str(b), "Error in sum")
    test_obj.run_test(str(a*b), str(a), "Error in prod")
        
    b = matrix.Matrix(3,3,[[1,2,3],[4,5,6],[7,8,9]])
    c = b.delete(1,1)
    d = matrix.Matrix(2,2,[[1,3],[7,9]])
    test_obj.run_test(str(c), str(d), "Error in put_out")
    
    a = matrix.Matrix(3,3,[[1,1,1],[1,1,1],[1,1,1]])
    test_obj.run_test(str(a.det(a)), str(0), "Error in det")
    
    a = matrix.Matrix(3,3,[[1,0,1],[0,1,0],[1,0,-1]])
    ident = a.identity()
    inv = a.inverse()
    test_obj.run_test(str(a*inv),str(ident), "Error in inverse")
    
    print a^3
    
    a = matrix.Matrix(1,1,[["Hola"]])
    
        
    test_obj.report_results()

test1()