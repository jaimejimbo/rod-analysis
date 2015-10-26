# -*- coding: utf-8 -*-
"""
Class:    Queue class
Author:   Jaime Perez Aparicio
Mail:     jaime.91@hotmail.es
License:  GPL
"""

class Queue(object):
    """
    Queue class: objects get in through one side and left through the other
    """
    
    def __init__(self):
        """
        Creator
        """
        self._data = []

    def __list__(self):
        """
        To list converter.
        """
        return self._data
        
    def join(self, element):
        """
        Adds something to the Queue
        """
        self._data.append(element)
        
    def get_next(self):
        """
        Returns next element
        """
        return self._data.pop(0)
        
    def get_pos(self, element):
        """
        Returns position of the element
        """
        return self._data.index(element) + 1

    def __len__(self):
        """
        Returns length of the queue
        """
        return len(self._data)

    def delete(self, element):
        """
        Deletes an element from the queue.
        """
        pos = self.get_pos(element)
        del self._data[pos-1]
