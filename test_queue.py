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
        a.put(1)
        a.put(2)
        a.put(3)
        a.put(4)
        a.put(5)
        self.assertEqual(len(a), 5, "Length method doesn't work. Obtained:"+str(len(a)))
        next = a.get()
        self.assertEqual(next, 1, "Check get function. Obtained:"+str(next))
        self.assertEqual(len(a), 4, "Length method doesn't work. Obtained:"+str(len(a)))
        self.assertEqual(a.get_possition(2), 1, "Check get_possition function Obtained:"+str(a.get_possition(3)))
        a.delete(3)
        self.assertEqual(len(a), 3, "Length method doesn't work. Obtained:"+str(len(a)))

