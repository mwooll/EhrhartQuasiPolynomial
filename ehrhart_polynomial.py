import numpy as np
from itertools import product

import sage.all
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.geometry.polyhedron.constructor import Polyhedron
# from sage.all import QQ, PolynomialRing, Polyhedron


from boundingbox import BoundingBox, get_bounding_box

from unittest import TestCase, main


R = PolynomialRing(QQ, "x")
x = R.gen()


def ehrhart_polynomial(vertices):
    y_values = points_contained(vertices)

    interpolation_points = [(k+1, y) for k, y in enumerate(y_values)]
    polynomial = R.lagrange_polynomial(interpolation_points)

    return polynomial

def points_contained(vertices):
    base_poly = Polyhedron(vertices)
    base_box = get_bounding_box(vertices)
    dimension = base_box.dim

    points_contained = [0]*(dimension+1)
    for k in range(1, dimension+2):
        poly = k*base_poly
        box = k*base_box

        contained = 0
        for point in box:
            if point in poly:
                contained += 1 

        points_contained[k-1] = contained

    return points_contained

def get_itertool_bounding_box(vertices):
    dimension = len(vertices[0])
    columns = [[vertex[d] for vertex in vertices]
               for d in range(dimension)]
    mins = [min(col) for col in columns]
    maxs = [max(col) for col in columns]
    return product(*[range(mins[d], maxs[d]+1) for d in range(dimension)])

class TestEhrhartPolynomial(TestCase):
    def test_ehrhart_polynomial(self):
        point_poly = ehrhart_polynomial([(0, 0, 0, 0, 0)])
        self.assertEqual(point_poly, 1)

        axis = ehrhart_polynomial([[-3], [3]])
        self.assertEqual(axis, 6*x + 1)

        triangle = ehrhart_polynomial([[0, 0], [1, 0], [0, 1]])
        self.assertEqual(triangle, x**2/2 + 1.5*x + 1)

        flat_tri = ehrhart_polynomial([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        self.assertEqual(flat_tri, x**2/2 + 1.5*x + 1)

        triangle_pyramid = ehrhart_polynomial([[0, 0, 0], [1, 0, 0],
                                                [0, 1, 0], [0, 0, 1]])
        self.assertEqual(triangle_pyramid, x**3/6 + x**2 + 11/6*x + 1)

        square_poly = ehrhart_polynomial([[0, 0], [1, 0], [1, 1], [0, 1]])
        self.assertEqual(square_poly, x**2 + 2*x + 1) # (x + 1)**2

        vertices = BoundingBox(np.zeros(3, dtype=int), np.ones(3, dtype=int)).vertices
        cube_poly = ehrhart_polynomial(vertices)
        self.assertEqual(cube_poly, x**3 + 3*x**2 + 3*x + 1) # (x + 1)**3


    def test_points_contained(self):
        point = points_contained([[0, 0, 0, 0, 0]])
        self.assertEqual(point, [1, 1, 1, 1, 1, 1])

        square = points_contained([[0, 0], [1, 0], [1, 1], [0, 1]])
        self.assertEqual(square, [4, 9, 16])

        vertices = BoundingBox(np.zeros(3, dtype=int), np.ones(3, dtype=int)).vertices
        cube = points_contained(vertices)
        self.assertEqual(cube, [8, 27, 64, 125])

        axis = points_contained([[-3], [3]])
        self.assertEqual(axis, [7, 13])

        triangle = points_contained([[0, 0], [1, 0], [0, 1]])
        self.assertEqual(triangle, [3, 6, 10])

        flat_tri = points_contained([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        self.assertEqual(flat_tri, [3, 6, 10, 15])

def itertool_bounding_box():
    vertices = [[0, 0, 0], [2, 0, 0], [0, -1, 0], [0, 0, 1]]
    box = get_itertool_bounding_box(vertices)
    for k in box:
        print(k)
if __name__ == "__main__":
   main()
   # itertool_bounding_box()