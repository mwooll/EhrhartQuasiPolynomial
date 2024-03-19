from ehrhart_polynomial import ehrhart_polynomial, points_contained, \
    get_bounding_box, get_bounding_extrema

import sage.all
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from unittest import TestCase, main


R = PolynomialRing(QQ, "x")
x = R.gen()

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