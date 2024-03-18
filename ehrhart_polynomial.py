from itertools import product

import sage.all
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.geometry.polyhedron.constructor import Polyhedron

from unittest import TestCase, main


R = PolynomialRing(QQ, "x")
x = R.gen()


def ehrhart_polynomial(vertices):
    y_values = points_contained(vertices)

    interpolation_points = [(k+1, y) for k, y in enumerate(y_values)]
    polynomial = R.lagrange_polynomial(interpolation_points)

    return polynomial

def points_contained(vertices):
    dimension = len(vertices[0])

    base_poly = Polyhedron(vertices)
    base_min, base_max = get_bounding_extrema(vertices, dimension)

    points_contained = [0]*(dimension+1)
    for k in range(1, dimension+2):
        poly = k*base_poly
        box = get_bounding_box(base_min, base_max, k)

        contained = 0
        for point in box:
            if point in poly:
                contained += 1 

        points_contained[k-1] = contained

    return points_contained

def get_bounding_extrema(vertices, dimension):
    columns = [[vertex[d] for vertex in vertices]
               for d in range(dimension)]
    mins = [min(col) for col in columns]
    maxs = [max(col) for col in columns]
    return mins, maxs

def get_bounding_box(mins, maxs, factor):
    return product(*[range(factor*mini, factor*maxi + 1)
                     for mini, maxi in zip(mins, maxs)])


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

        vertices = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1),
                    (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)]
        cube_poly = ehrhart_polynomial(vertices)
        self.assertEqual(cube_poly, x**3 + 3*x**2 + 3*x + 1) # (x + 1)**3


    def test_points_contained(self):
        point = points_contained([[0, 0, 0, 0, 0]])
        self.assertEqual(point, [1, 1, 1, 1, 1, 1])

        square = points_contained([[0, 0], [1, 0], [1, 1], [0, 1]])
        self.assertEqual(square, [4, 9, 16])

        vertices = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1),
                    (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)]
        cube = points_contained(vertices)
        self.assertEqual(cube, [8, 27, 64, 125])

        axis = points_contained([[-3], [3]])
        self.assertEqual(axis, [7, 13])

        triangle = points_contained([[0, 0], [1, 0], [0, 1]])
        self.assertEqual(triangle, [3, 6, 10])

        flat_tri = points_contained([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        self.assertEqual(flat_tri, [3, 6, 10, 15])

    def test_get_bounding_extreme(self):
        vertices = [[0, 1, 2], [4, -1, 2]]
        expected = ([0, -1, 2], [4, 1, 2])
        self.assertEqual(get_bounding_extrema(vertices, 3), expected)

        vertex = [[0, 0, 0, 0, 0]]
        expected = ([0, 0, 0, 0, 0], [0, 0, 0, 0, 0])
        self.assertEqual(get_bounding_extrema(vertex, 5), expected)

        vertices = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1),
                    (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)]
        expected = ([0, 0, 0], [1, 1, 1])
        self.assertEqual(get_bounding_extrema(vertices, 3), expected)

    def test_get_bounding_box(self):
        mins, maxs = [0, -1, -2], [1, 2, 3]
        expected = [(0, 0, 0)]
        self.assertEqual(list(get_bounding_box(mins, maxs, 0)), expected)

        mins, maxs = [-1, 2], [0, 3]
        expected = [(-1, 2), (-1, 3), (0, 2), (0, 3)]
        self.assertEqual(list(get_bounding_box(mins, maxs, 1)), expected)

        mins, maxs = [0, 0, 0], [1, 1, 1]
        expected = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1),
                    (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)]
        self.assertEqual(list(get_bounding_box(mins, maxs, 1)),
                         expected)

        mins, maxs = [-5], [5]
        expected = [(k,) for k in range(-10, 11)]
        self.assertEqual(list(get_bounding_box(mins, maxs, 2)), expected)

if __name__ == "__main__":
   main()
