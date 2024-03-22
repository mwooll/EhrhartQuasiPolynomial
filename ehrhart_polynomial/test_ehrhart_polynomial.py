from ehrhart_polynomial import ehrhart_polynomial , \
    points_contained, get_bounding_box, get_bounding_extrema, \
    simplify_vertices, drop_dimensions, drop_constant_dimensions, \
    scale_down_vertices

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
        point = points_contained([[0, 0, 0, 0, 0]], False)
        self.assertEqual(point, [1, 1, 1, 1, 1, 1])

        square = points_contained([[0, 0], [1, 0], [1, 1], [0, 1]], False)
        self.assertEqual(square, [4, 9, 16])

        vertices = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1),
                    (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)]
        cube = points_contained(vertices, False)
        self.assertEqual(cube, [8, 27, 64, 125])

        axis = points_contained([[-3], [3]], False)
        self.assertEqual(axis, [7, 13])

        triangle = points_contained([[0, 0], [1, 0], [0, 1]], False)
        self.assertEqual(triangle, [3, 6, 10])

        flat_tri = points_contained([[0, 0, 0], [1, 0, 0], [0, 1, 0]], False)
        self.assertEqual(flat_tri, [3, 6, 10, 15])

    def test_points_contained_simplified(self):
        flat_tri = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        triangle = [[0, 0], [1, 0], [0, 1]]

        self.assertEqual(points_contained(triangle, True),
                         points_contained(triangle, False))
        self.assertEqual(points_contained(flat_tri, True),
                         points_contained(triangle, False))

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

    def test_simplify_vertices(self):
        triangle = [[0, 0], [1, 0], [1, 1]]
        mins, maxs = get_bounding_extrema(triangle, 2)
        self.assertEqual(simplify_vertices(triangle, 2),
                         (triangle, mins, maxs, 2, 1))
                         
        scaled_extra_dim = [[0, 0, 1], [3, 0, 1], [3, 3, 1]]
        self.assertEqual(simplify_vertices(scaled_extra_dim, 3),
                         (triangle, mins, maxs, 2, 3))

    def test_drop_dimensions(self):
        to_reduce = [[0, 0], [2, 0], [1, 1]]
        all_True = [True, True]
        self.assertEqual(drop_dimensions(to_reduce, all_True), to_reduce)
        
        mixed_bools = [True, False]
        self.assertEqual(drop_dimensions(to_reduce, mixed_bools), [[0], [2], [1]])

        all_False = [False, False]
        self.assertEqual(drop_dimensions(to_reduce, all_False), [[], [], []])

    def test_drop_constant_dimensions(self):
        simple = [[0, 0], [1, 0], [0, 1]]
        expected = (simple, *get_bounding_extrema(simple, 2), 2)
        self.assertEqual(drop_constant_dimensions(simple, 2), expected)

        extra_dim = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        self.assertEqual(drop_constant_dimensions(extra_dim, 3), expected)

    def test_scale_down_vertices(self):
        simple = [[0, 0], [1, 0], [0, 1]]
        self.assertEqual(scale_down_vertices(simple), (simple, 1))

        scaled = [[0, 0], [2, 0], [0, 2]]
        self.assertEqual(scale_down_vertices(scaled), (simple, 2))


if __name__ == "__main__":
   main()