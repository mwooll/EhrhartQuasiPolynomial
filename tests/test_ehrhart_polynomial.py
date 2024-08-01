from ehrhart_polynomial import (ehrhart_polynomial,
                                points_contained_sequence, points_contained,
                                get_period, get_bounding_extrema,
                                get_bounding_box, get_bounding_box_rational,
                                simplify_vertices, drop_constant_dimensions,
                                drop_dimensions, scale_down_vertices,
                                IntegerPeriodicFunctionRing,
                                QuasiPolynomialRing)

import sage.all
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.geometry.polyhedron.constructor import Polyhedron

from unittest import TestCase, main


R = PolynomialRing(QQ, "x")
x = R.gen()

QPR = QuasiPolynomialRing(QQ)

class TestEhrhartPolynomial(TestCase):
    def test_ehrhart_polynomial(self):
        # integral polytopes
        point_poly = ehrhart_polynomial([(0, 0, 0, 0, 0)])
        self.assertEqual(point_poly, QPR([1]))

        axis = ehrhart_polynomial([[-3], [3]])
        self.assertEqual(axis, QPR([1, 6]))

        triangle = ehrhart_polynomial([[0, 0], [1, 0], [0, 1]])
        self.assertEqual(triangle, QPR([1, 1.5, 0.5]))

        flat_tri = ehrhart_polynomial([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        self.assertEqual(flat_tri, QPR([1, 1.5, 0.5]))

        triangle_pyramid = ehrhart_polynomial([[0, 0, 0], [1, 0, 0],
                                                [0, 1, 0], [0, 0, 1]])
        self.assertEqual(triangle_pyramid, QPR([1, 11/6, 1, 1/6]))

        square_poly = ehrhart_polynomial([[0, 0], [1, 0], [1, 1], [0, 1]])
        self.assertEqual(square_poly, QPR([1, 2, 1])) # (x + 1)**2

        vertices = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1),
                    (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)]
        cube_poly = ehrhart_polynomial(vertices)
        self.assertEqual(cube_poly, QPR([1, 3, 3, 1]))# (x + 1)**3

        # rational polytopes
        rational_point = [[1/7, 1/2]]
        values = [1] + [0]*13
        self.assertEqual(ehrhart_polynomial(rational_point), QPR([values]))

        rational_triangle = [(0, 0), (3/2, 0), (0, 1/3)]
        self.assertEqual(ehrhart_polynomial(rational_triangle),
                         QPR([[1, 3/4], 1, 1/4]))


    def test_ehrhart_polynomial_simplified(self):
        triangle = [[0, 0], [1, 0], [0, 1]]
        self.assertEqual(ehrhart_polynomial(triangle, True),
                         ehrhart_polynomial(triangle, False))

        bloated = [[0, 0, 2], [3, 0, 2], [0, 3, 2]]
        self.assertEqual(ehrhart_polynomial(bloated, True),
                         ehrhart_polynomial(bloated, False))

        giant_rectangle = [[0, 0], [100, 0], [100, 50], [0, 50]]
        self.assertEqual(ehrhart_polynomial(giant_rectangle, True),
                         QPR([1, 150, 5000]))


    # points contained
    def test_points_contained_sequence(self):
        # integral
        point = points_contained_sequence([[0, 0, 0, 0, 0]], False)
        self.assertEqual(point,
                         ([1, 1, 1, 1, 1, 1], 1, 1))

        square = points_contained_sequence([[0, 0], [1, 0], [1, 1], [0, 1]], False)
        self.assertEqual(square, ([4, 9, 16], 1, 1))

        vertices = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1),
                    (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)]
        cube = points_contained_sequence(vertices, False)
        self.assertEqual(cube, ([8, 27, 64, 125], 1, 1))

        axis = points_contained_sequence([[-3], [3]], False)
        self.assertEqual(axis, ([7, 13], 1, 1))

        triangle = points_contained_sequence([[0, 0], [1, 0], [0, 1]], False)
        self.assertEqual(triangle, ([3, 6, 10], 1, 1))

        flat_tri = points_contained_sequence([[0, 0, 0], [1, 0, 0], [0, 1, 0]], False)
        self.assertEqual(flat_tri, ([3, 6, 10, 15], 1, 1))

        # rational
        half_unit_square = [(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]
        self.assertEqual(points_contained_sequence(half_unit_square, False),
                         ([1, 4, 4, 9, 9, 16], 1, 2))

        rational_triangle = [(0, 0), (3/2, 0), (0, 1/3)]
        self.assertEqual(points_contained_sequence(rational_triangle, True),
                         ([2, 4, 6, 9, 12, 16, 20, 25, 30,
                           36, 42, 49, 56, 64, 72, 81, 90, 100], 1, 6))

    def test_points_contained_sequence_simplified(self):
        triangle = [[0, 0], [1, 0], [0, 1]]
        flat_tri = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]

        self.assertEqual(points_contained_sequence(triangle, True),
                         points_contained_sequence(triangle, False))
        self.assertEqual(points_contained_sequence(flat_tri, True),
                         points_contained_sequence(triangle, False))

        bloated_square = [[0, 0], [3, 0], [3, 3], [3, 3]]
        self.assertEqual(points_contained_sequence(bloated_square, True),
                         (points_contained_sequence(triangle, False)[0], 3, 1))

    def test_points_contained(self):
        vertices = [[0, 0], [1, 0], [0, 1]]
        triangle = Polyhedron(vertices)

        triangle_box = get_bounding_box(*get_bounding_extrema(vertices, 2), 1)
        self.assertEqual(points_contained(triangle, triangle_box), 3)

        triangle_box_2 = get_bounding_box(*get_bounding_extrema(vertices, 2), 2)
        self.assertEqual(points_contained(triangle*2, triangle_box_2), 6)

        empty_poly = [[1/2, 0], [0, 1/2]]
        square_box = get_bounding_box([0, 0], [1, 1], 1)
        self.assertEqual(points_contained(empty_poly, square_box), 0)


    # period
    def test_get_period(self):
        cube = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1),
                (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)]
        self.assertEqual(get_period(cube), 1)

        simple = [(0, 0), (1/2, 0), (3/2, 1), (1, 1)]
        self.assertEqual(get_period(simple), 2)

        fractions = [(1/2, 0), (0, 1/3), (1/4, 1), (1, 1/5)]
        self.assertEqual(get_period(fractions), 60)


    # bounding box
    def test_get_bounding_extreme(self):
        vertices = [[0, 1, 2], [4, -1, 2]]
        expected = ([0, -1, 2], [4, 1, 2])
        self.assertEqual(get_bounding_extrema(vertices, 3), expected)

        vertex = [[0, 0, 0, 0, 0]]
        expected = ([0, 0, 0, 0, 0], [0, 0, 0, 0, 0])
        self.assertEqual(get_bounding_extrema(vertex, 5), expected)

        cube = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1),
                (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)]
        expected = ([0, 0, 0], [1, 1, 1])
        self.assertEqual(get_bounding_extrema(cube, 3), expected)

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

    def test_get_bounding_box_rational(self):
        mins, maxs = [0, 1/2, 1/3], [1, 1, 1/3]
        self.assertEqual(list(get_bounding_box_rational(mins, maxs, 1)), [])
        self.assertEqual(list(get_bounding_box_rational(mins, maxs, 2)), [])

        self.assertEqual(list(get_bounding_box_rational(mins, maxs, 3)),
                         list(get_bounding_box((0, 2, 1), (3, 3, 1), 1)))

        self.assertEqual(list(get_bounding_box_rational(mins, maxs, 6)),
                         list(get_bounding_box((0, 3, 2), (6, 6, 2), 1)))

        self.assertEqual(list(get_bounding_box_rational([-1/2], [1/3], 7)),
                         [(-3,), (-2,), (-1,), (0,), (1,), (2,)])

    # simplify polytope
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