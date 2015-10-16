# -*- coding: utf-8 -*-
"""
Class:    Recursive list
Author:   Jaime Pérez Aparicio
Mail:     jaime.91@hotmail.es
License:  GPL
"""

import tree

def RecursiveList(Tree):
    """
    Recursive list
    """
    
    def __init__(self, value):
        """
        Creator
        """
        super(value)

    def join(self, other_list):
        """
        Join a list
        """
        last_element = other_list.last_element()
        self._parent = last_element
        other_list.last_element().add_son(self)
    
    def last_element(self):
        """
        Returns the last element of the list
        """
        if len(self.sons) == 0:
            return self
        if len(self.sons) == 1:
            return self.sons[0].last_element()
        else:
            raise IndexError
    
    def first_element(self):
        """
        Returns the first element of the list
        """
        if self.parent == None:
            return self
        else:
            return self.parent.first_element()
    
    def leave(self):
        """
        Leaves list
        """
        parent = self.parent
        son = self.sons[0]
        parent.remove_son(self)
        parent.add_son(son)
        son.parent = parent
    
    def __str__(self):
        """
        Prints the whole list
        """
        first_element = self.first_element()
        output = first_element.get_output()
        return output

    def get_output(self):
        """
        Function for string
        """
        output = ""
        output += str(self.data)
        if len(self.sons) == 1:
            output += son.get_output()
        return output

