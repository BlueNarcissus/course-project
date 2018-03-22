#!/usr/bin/env python
""" class module to create all the legal steps of chess """

""" base class with abtract methods"""
class piece(object):
    def __init__(self, x0, y0):
        (self._x0, self._y0) = (x0, y0)
    
    # a legal move should not be the current position
    def is_legal(self, x, y):
        pass
    
    # use self.is_legal() method to check all possible moves
    def all_legal(self):
        pass

""" four subclasses """

# A Rook can move any number of steps forwards, backwards, left, or right
class Rook(piece):
    def is_legal(self, x, y):
        if (self._x0, self._y0) != (x, y) and (x- self._x0) * (y - self._y0) == 0:
            return True
        return False
    
    def all_legal(self):
        list = []
        for i in range(1, 9):
            if self._x0 != i :
                list.append((i, self._y0))
        for j in range(1, 9):
            if self._y0 != j :
                list.append((self._x0, j))
        return list

# A Bishop can move any number of steps diagonally in a straight line
class Bishop(piece):
    def is_legal(self, x, y):
        if (self._x0, self._y0) != (x, y) and abs(x- self._x0) == abs(y - self._y0):
            return True
        return False
    
    def all_legal(self):
        list = []
        for i in range(1, 9):
            for j in range(1, 9):
                if (self._x0, self._y0) != (i, j) and abs(i- self._x0) == abs(j - self._y0):
                    list.append((i,j))
        return list

# A King can move one step diagonally, forwards, backwards, left, or right
class King(piece):
    # Check a legal move
    def is_legal(self, x, y):
        if (self._x0, self._y0) != (x, y) and abs(x- self._x0) <= 1 and abs(y - self._y0)<= 1:
            return True
        return False
    
    def all_legal(self):
        list = []
        for i in range(1, 9):
            for j in range(1, 9):
                if (self._x0, self._y0) != (i, j) and abs(i- self._x0) <= 1 and abs(j - self._y0)<= 1:
                    list.append((i,j))
        return list

# A Queen can move any number of steps diagonally, forwards, backwards, left, or right
class Queen(piece):
    def is_legal(self, x, y):
        if (self._x0, self._y0) != (x, y):
            if abs(x- self._x0) == abs(y - self._y0) or (x- self._x0) * (y - self._y0) == 0:
                return True
        return False
    
    def all_legal(self):
        list = []
        for i in range(1, 9):
            for j in range(1, 9):
                if (self._x0, self._y0) != (i, j):
                    if abs(i- self._x0) == abs(j - self._y0) or (i- self._x0) * (j - self._y0) == 0:
                        list.append((i,j))
        return list
