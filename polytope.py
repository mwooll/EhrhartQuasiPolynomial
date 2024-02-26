from copy import deepcopy
import numpy as np

from boundingbox import BoundingBox

from unittest import TestCase, main

class Polytope:
    def __init__(self, vertices):
        """
        
        @vertices : nested list of same-length lists containing numbers
            inner lists define the vertices
            
        Note: the resulting polytope will be the convex hull of the vertices
        
        """
        self.vertices = np.array(vertices)
        
        # number of rows = number of vertices
        # number of columns = dimension of vertices
        self.num, self.dim = self.vertices.shape
        
        #self.origin = np.zeros(self.dim)
        
        
    def shift_to_origin(self):
        """
        shifts the polytope such that one vertex lies on the origin
        and that the first coordinate of no vertex is negative
        """        
        # subtracting the point with the smallest (signed) value in the 
        # first coordinate gives us the desired shift
        lowest = self.vertices[self.vertices[:, 0].argsort()[0]]
        return self - lowest

    def needs_to_be_shifted(self):
        if (0 in np.linalg.norm(self.vertices, axis=1) and 
                all(row[0] >= 0 for row in self.vertices)):
            return False
        return True

    def get_bounding_box(self):
        # might be helpful for the rational scalar case
        # min_value = np.array(np.floor(np.min(self.vertices, axis=0)), dtype=int)
        # max_value = np.array(np.ceil (np.max(self.vertices, axis=0)), dtype=int)
        
        min_value = np.min(self.vertices, axis=0)
        max_value = np.max(self.vertices, axis=0)
        return BoundingBox(min_value, max_value)
            
    """
    standard dunder functions
    """
    def __str__(self):
        return f"Polytope with vertices: \n{self.vertices}"

    def __repr__(self):
        return f"Polytope(\n{self.vertices})"


    """
    Math support
        makes extensive use of numpy.ndarray's shape casting
    """
    def __mul__(self, number):
        """
        @number : must be convertible to float
        """
        return Polytope(number*self.vertices)

    __rmul__ = __mul__
    
    def __sub__(self, point):
        """
        @point : must be a np.array of the same width as self.vertices
        """
        return Polytope(self.vertices - point)
        
    def __add__(self, point):
        """ 
        @point : must be a np.array of the same width as self.vertices
        """
        return Polytope(self.vertices + point)
    
    # do not want:
    # __radd__ = __add__
    # breaks numpy's shape casting


class TestPolytope(TestCase):
    def test_bounding_box(self):
        poly = Polytope([[0, 0], [1, 0], [0, 1]])
        self.assertEqual(poly.get_bounding_box(), 
                         BoundingBox(np.array([0, 0]), np.array([1, 1])))
    
if __name__ == "__main__":
    main()
    
    
    vertices = [[1, 0, 2], [3, 0, 3], [2, 1, 2], [4, 2, 1]]
    poly = Polytope(np.array(vertices))
    print(vertices)
    for point in poly.bounding_box():
        print(point)
        
    
