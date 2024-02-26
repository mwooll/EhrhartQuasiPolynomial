import numpy as np

from polytope import Polytope

class ConvexPolytope(Polytope):
    def __init__(self, vertices):
        """
        
        @vertices : 2-dimensional np.array or nested list of same-length lists containing numbers
            rows/inner lists define the vertices
            
        Note: 
            the resulting polytope will be the convex hull of the vertices
            and therefore the order of the points does not matter
        
        """
        super().__init__(vertices)
    
    def ehrhart_polynomial(self):
        shifted = self.shift_to_origion()
        box = shifted.bounding_box()
        
        for point in box:
            pass
        return
    
    def __contains__(self, point):
        if self.dim != len(point):
            return False
        
        
        
        
if __name__ == "__main__":
    vertices = [[1, 2, 3], [3, 4, 5], [1, 2, 3]]
    poly = ConvexPolytope(vertices)
    print(poly + np.array([1, 1, 1]))
