#!/usr/bin/python
import unittest
from test_matrix import TestMatrix
from test_rod import TestRod
from test_queue import TestQueue
from StringIO import StringIO
from pprint import pprint

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

stream = StringIO()
runner = unittest.TextTestRunner(stream = stream)

def suite():
    suites = []
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestMatrix))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestRod))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestQueue))
    return unittest.TestSuite(suites)

result = runner.run(suite())
print 'Tests run ', result.testsRun
print bcolors.WARNING + 'ERRORS ' + bcolors.FAIL + bcolors.BOLD + str(result.errors)
pprint(result.failures)
stream.seek(0)
if result.errors==[] and len(result.failures)==0: print bcolors.OKGREEN + bcolors.BOLD + str("NO ERRORS") + bcolors.ENDC
print bcolors.ENDC+'Test output\n', stream.read() + bcolors.ENDC
raw_input("Press Enter to end.")

