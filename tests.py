import unittest
from test_matrix import TestMatrix
from test_rod import TestRod
from StringIO import StringIO
from pprint import pprint

stream = StringIO()
runner = unittest.TextTestRunner(stream = stream)

def suite():
    suites = []
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestMatrix))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestRod))
    return unittest.TestSuite(suites)

result = runner.run(suite())
print 'Tests run ', result.testsRun
print 'Errors ', result.errors
pprint(result.failures)
stream.seek(0)
print 'Test output\n', stream.read()

