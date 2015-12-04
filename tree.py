# -*- coding: utf-8 -*-
"""
Class:    Tree class
Author:   Jaime Pérez Aparicio
Mail:     jaime.91@hotmail.es
License:  GPL
"""

class Tree(object):
    """
    Tree class. It allows to create Tree objects.
    """
    
    
    
    
    def __init__(self, value = None, input_parent = None, input_sons = [], repeated_allowed = False):
        """
        Init function
        """
        global ID
        self._parent = input_parent        
        self._repeated_allowed = False
        self._sons = []
        for son in input_sons:
            self.add_son(son)
        #Stablish the parent for sons
        for child in input_sons:
            child.parent = self
        self._value = value
        #Add this to parents' sons
        if input_parent != None:
            input_parent.add_son(self)
        self._repeated_allowed = repeated_allowed

        
        
        
        
    @property
    def parent(self):
        """
        Gets the parent
        """
        return self._parent
    
    @parent.setter
    def parent(self, value):
        """
        Changes the parent
        """
        self._parent = value




    @property
    def sons(self):
        """
        Gets the sons
        """
        return self._sons
        
    def add_son(self, son):
        """
        Add a son
        """
        if self._repeated_allowed:
            self._sons.append(son)
            son.parent = self
        elif not son in self:
            if self.find_branch_up(son) == [-1,-1]:
                if id(self) != id(son):
                    self._sons.append(son)
                    son.parent = self

    def clear_sons(self):
        """
        Clean sons list
        """
        self._sons = []

    def remove_son(self, input_son):
        """
        Removes son of the sons list
        """
        self.sons.remove(input_son)




    @property
    def value(self):
        """
        Gets the value
        """
        return self._value
        
    @value.setter
    def value(self,_value):
        """
        Changes the value
        """
        self._value = _value
        
        
        
        
        
        
    def __str__(self):
        """
        Allows to print the Tree, beginning in this branch
        """
        output = "["
        output += str(self._value)
        if len(self.sons) == 1:
            output += str(self.sons[0])
        elif len(self.sons) > 1:
            output += "["
            for child in self.sons:
                next_str = str(child)
                next_str = next_str[1:-1]
                output += next_str
                output += ", "
            output = output[:-2]
            output += "]"
        output += "]"
        return output

    def __eq__(self, input_tree):
        """
        Operator ==
        """
        equal = False
        if id(self) == id(input_tree):
            return True
        if id(self) == id(None) or id(input_tree) == id(None):
            return False
        if self.value == input_tree.value:
            if len(self.sons) != len(input_tree.sons):
                return False
            if len(self.sons) == 0:
                return True
            for child1 in self.sons:
                equal = False
                for child2 in input_tree.sons:
                    if child1 == child2:
                        equal = True
                if not equal:
                    return False
            return True
        else:
            return False
            
    def __neq__(self, input_tree):
        """
        Operator != 
        """
        return not self == input_tree

    def __contains__(self, input_tree):
        """
        Overrides in operator
        """
        if len(self.find_branch_down(input_tree)) == 0:
            return False
        return True
        
    def __getitem__(self, index):
        """
        Returns branch
        """
        return self.sons[index]






    def copy(self):
        """
        Returns a copy of this tree with head here BUGS
        """
        #main and copy should not be connected
        output_tree = Tree(value = self.value)
        for child in self.sons:
            new_son = child.copy()
            new_son.parent = output_tree
            #I use copy in next line for fixing pointing issues
            output_tree.add_son(new_son.copy())
        return output_tree
        
    def remove_repeated(self):
        """
        Removes repeated branches in the tree
        If everything is alright, this function should not be needed.
        """
        try:
            for index1 in range(len(self.sons)-1):
                for index2 in range(index1 + 1 , len(self.sons)):
                    son1, son2 = self.sons[index1], self.sons[index2]
                    if id(son1) == id(son2):
                        self.remove_son(son2)
                        index2 -= 1
                son1.remove_repeated()
        except:
            pass
        
    def remove_loops(self):
        """
        Removes loops in tree with ID
        """
        self.remove_repeated()
        self.find_and_delete(self)
        for child1 in self.sons:
            child1.find_and_delete(self)
            for child2 in self.sons:
                if id(child1) != id(child2):
                    child1.find_and_delete(child2)
            #child1.remove_loops()
        
            
    def find_and_delete(self, input_tree):
        """
        Delete son with same id that input_tree
        """
        for child in self.sons:
            if id(child) == id(self):
                self.remove_son(child)
        
        
    def find_branch_down(self, input_tree):
        """
        Finds a branch inside a branch with ID
        """
        output = []
        input_id = id(input_tree)
        next_id = 0
        for child in self.sons:
            index = self.sons.index(child)
            if id(child) == id(input_tree):
                output.append(index)
            else:
                output2 = child.find_branch_down(input_tree)
                if len(output2) > 0:
                    output.append(index)
                    output.append(output2)
        return output

    def find_branch_up(self, input_tree):
        """
        Finds a branch upside a branch with ID [number_of_jumps, son_number]
            number_of_jumps: how far is the parent from self (1-> self's parent)
            son_number: index of the son
        """
        output = [-1, -1]
        if self.parent == None:
            return output
        if id(self.parent) == id(input_tree):
            output = [0, -1]
        for index in range(len(self.parent.sons)):
            sibiling = self.parent.sons[index]
            if id(sibiling) != id(self):
                if id(sibiling) == id(input_tree):
                    output = [0, index]
        if output[0] == -1:
            output2 = self.parent.find_branch_up(input_tree)
            if output2[0] != -1:
                output[0] = output2[0] + 1
                output[1] = output2[1]
        return output
