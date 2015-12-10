# -*- coding: utf-8 -*-
"""
Class:    Queue class
Author:   Jaime Perez Aparicio
Mail:     jaime.91@hotmail.es
License:  GPL
"""
import re

class Queue(object):
    """
    Queue class: objects get in through one side and left through the other
    """
    
    def __init__(self, queue=None):
        """
        Creator
        """
        input_type = str(type(queue))
        if re.match(r".*'NoneType'.*",input_type):
            self._data = []
        elif re.match(r'.*list.*',input_type):
            self._data = [element for element in queue]
        elif re.match(r'.*Queue.*',input_type):
            self._data = [element for element in queue]
        else:
            msg = "Queue input data must be list or another"
            msg += " queue. "+str(type(queue))
            raise ValueError(msg)

    @property     
    def clone(self):
        """
        Returns a clone of the queue.
        """
        clone = Queue(self._data)
        return clone

    def __iter__(self):
        """
        To list converter.
        """
        for element in self._data:
            yield element
        
    def put(self, element):
        """
        Adds something to the Queue
        """
        self._data.append(element)
        
    def get(self):
        """
        Returns next element
        """
        return self._data.pop(0)
        
    def get_possition(self, element):
        """
        Returns position of the element
        """
        return self._data.index(element)

    def __set__(self):
        """
            Magic method for set conversion.
        """
        return set(self._data)

    def __len__(self):
        """
        Returns length of the queue
        """
        return len(self._data)

    def delete(self, element):
        """
        Deletes an element from the queue.
        """
        pos = self.get_possition(element)
        del self._data[pos]

