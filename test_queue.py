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
        self.assertEqual(len(a), 5, "Length method doesn't work. Obtained:"+str(len(a)))
        next = a.get_next()
        self.assertEqual(next, 1, "Check get_next function. Obtained:"+str(next))
        self.assertEqual(len(a), 4, "Length method doesn't work. Obtained:"+str(len(a)))
        self.assertEqual(a.get_pos(2), 1, "Check get_pos function Obtained:"+str(a.get_pos(3)))
        a.delete(3)
        self.assertEqual(len(a), 3, "Length method doesn't work. Obtained:"+str(len(a)))

