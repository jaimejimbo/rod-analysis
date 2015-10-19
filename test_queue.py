import unittest
import queue

class TestQueue(unittest.TestCase):
    """
    """

    
    def SetUp(self):
        self.longMessage = True


    def test_Queue(self):
        """
        Queue tests
        """
        a = queue.Queue()
        a.join(1)
        a.join(2)
        a.join(3)
        a.join(4)
        a.join(5)
        self.assertEqual(a.get_next(), 1, "Check get_next function. Obtained:"+str(a.get_next()))
        self.assertEqual(a.get_pos(3), 1, "Check get_pos function Obtained:"+str(a.get_pos(3)))

