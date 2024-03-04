import numpy as np

from util import scalar_check, point_check
from boundingbox import BoundingBox

from unittest import TestCase, main

class ConvexPolytope:
    def __init__(self, vertices):
        if not isinstance(vertices, np.ndarray):
            try:
                vertices = np.array(vertices)
            except:
                raise TypeError("vertices must be a numpy.ndarray or convertible to "
                                + f"a numpy.ndarray but is of type {type(vertices)}")
        if vertices.ndim != 2:
            raise ValueError("vertices need to be 2-dimensional but is of " +
                             f"dimension {vertices.ndim}")
        self.vertices = vertices

        self.vertex_set = self._get_vertex_set(vertices)

        self.box = self.get_bounding_box()

        # number of rows = number of vertices
        # number of columns = dimension of vertices
        self.num, self.dim = self.vertices.shape

    """
    methods called in __init__
    """
    def _get_vertex_set(self, vertices):
        if isinstance(vertices, (tuple, list)):
            try:
                vertex_set = set(vertices)
                return vertex_set
            except TypeError:
                pass

        vertex_set = set()
        for row in self.vertices:
            vertex_set.add(tuple(row))
        return vertex_set

    """
    main functionality
    """
    def shift_to_origin(self):
        """
        Returns a new ConvexPolytope such that one vertex lies on the origin
        and that the first coordinate of no vertex is negative
        """
        # subtracting the point with the smallest (signed) value in the
        # first coordinate gives us the desired shift
        lowest = self.vertices[self.vertices[:, 0].argsort()[0]]
        return self.shift(-lowest)

    def needs_to_be_shifted(self):
        if (0 in np.linalg.norm(self.vertices, axis=1) and 
                all(row[0] >= 0 for row in self.vertices)):
            return False
        return True

    def get_bounding_box(self):
        min_value = np.array(np.floor(np.min(self.vertices, axis=0)), dtype=int)
        max_value = np.array(np.ceil (np.max(self.vertices, axis=0)), dtype=int)
        return BoundingBox(min_value, max_value)

    def ehrhart_polynomial(self):
        shifted = self.shift_to_origion()
        box = shifted.bounding_box()

        for point in box:
            pass
        return

    def __contains__(self, point):
        point_check(point, self.dim)

        if point not in self.box:
            return False


    """
    standard dunder methods
    """
    def __str__(self):
        return f"ConvexPolytope with vertices: \n{self.vertices}"

    def __repr__(self):
        return f"ConvexPolytope(\n{self.vertices})"

    """
    comparison operators
    """
    def comparable(self, other):
        return isinstance(other, ConvexPolytope) and self.dim == other.dim

    def compare(self, other):
        if not isinstance(other, ConvexPolytope):
            raise TypeError("other needs to be of type ConvexPolytope " +
                            f" but is of type {type(other)}")
        if self.dim != other.dim:
            raise ValueError("self and other need to have the same dimension " +
                             f"but have {self.dim} and {other.dim} respectively")

    def __eq__(self, other):
        return self.comparable(other) and self.vertex_set == other.vertex_set

    # def __le__(self, other): return False

    # def __lt__(self, other): return False

    """
    math support
    """
    def __mul__(self, number):
        """
        Returns a scaled ConvexPolytope
        """
        scalar_check(number)
        return ConvexPolytope(number*self.vertices)

    __rmul__ = __mul__

    def __neg__(self):
        return (-1)*self

    def __truediv__(self, number):
        """
        Returns a scaled ConvexPolytope
        """
        scalar_check(number)
        return ConvexPolytope(self.vertices/number)

    def shift(self, point):
        """
        Returns a shifted ConvexPolytope
        """
        point_check(point, self.dim)
        return ConvexPolytope(self.vertices + point)

    # def union(self, other): return


class TestConvexPolytope(TestCase):
    """
    __init__
    """
    def test_init_errors(self):
        poly = ConvexPolytope([(0, 0), (1, 0), (0, 1)])
        self.assertRaises(TypeError, ConvexPolytope.__init__, poly)

        self.assertRaises(ValueError, ConvexPolytope, 1)
        self.assertRaises(ValueError, ConvexPolytope, [1])
        self.assertRaises(ValueError, ConvexPolytope, np.array(1))

        self.assertRaises(ValueError, ConvexPolytope, [(-1, -1, -1)])

    """
    main functionality
    """
    def test_needs_to_be_shifted(self):
        poly = ConvexPolytope([(0, 0), (1, 0), (0, 1)])
        self.assertFalse(poly.needs_to_be_shifted())

        square = ConvexPolytope([(-1, -1), (-1, 1), (1, -1), (1, 1)])
        self.assertTrue(square.needs_to_be_shifted())

        triangle = ConvexPolytope([(-1, -1), (-1, 1), (0, 0)])
        self.assertTrue(triangle.needs_to_be_shifted())

    def test_shift_to_origin(self):
        poly = ConvexPolytope([(0, 0), (1, 0), (0, 1)])
        self.assertEqual(poly.shift_to_origin(), poly)

        square = ConvexPolytope([(1, 1), (1, 2), (2, 1), (2, 2)])
        shifted_square = ConvexPolytope([(0, 0), (0, 1), (1, 0), (1, 1)])
        self.assertEqual(square.shift_to_origin(), shifted_square)

        triangle = ConvexPolytope([(-1, -1), (-1, 1), (0, 0)])
        shifted_triangle = ConvexPolytope([(0, 0), (0, 2), (1, 1)])
        self.assertEqual(triangle.shift_to_origin(), shifted_triangle)

    def test_bounding_box(self):
        poly = ConvexPolytope([(0, 0, 0), (1, 0, -1), (0, 1, 2)])
        self.assertEqual(poly.get_bounding_box(), 
                         BoundingBox(np.array([0, 0, -1]), np.array([1, 1, 2])))

        rational = ConvexPolytope([(0, 1.68), (3.1, -0.1)])
        self.assertEqual(rational.get_bounding_box(),
                         BoundingBox(np.array([0, -1]), np.array([4, 2])))

    """
    comparison operators
    """
    def test_comparable_false(self):
        poly = ConvexPolytope([(0, 0), (1, 0), (0, 1)])
        self.assertFalse(poly.comparable(4))
        self.assertFalse(poly.comparable((0, 0)))

        wrong_dimension = ConvexPolytope([(0, 0, 0)])
        self.assertFalse(poly.comparable(wrong_dimension))
        self.assertFalse(wrong_dimension.comparable(poly))

    def test_comparable_true(self):
        poly = ConvexPolytope([(0, 0), (1, 0), (0, 1)])
        larger = ConvexPolytope([(0, 0), (3, 0), (0, 3)])

        self.assertTrue(poly.comparable(larger))
        self.assertTrue(poly.comparable(poly))

    def test_compare_error(self):
        poly = ConvexPolytope([(0, 0), (1, 0), (0, 1)])
        wrong_dimension = ConvexPolytope([(0, 0, 0)])

        self.assertRaises(TypeError, poly.compare, 3)
        self.assertRaises(ValueError, poly.compare, wrong_dimension)

    def test_eq(self):
        poly = ConvexPolytope([(0, 0), (1, 0), (0, 1)])
        other = ConvexPolytope([(0, 0), (1, 0), (0, 1)])
        larger = ConvexPolytope([(0, 0), (3, 0), (0, 3)])

        self.assertFalse(poly == larger)
        self.assertTrue(poly == other)

    """
    math support
    """
    def test_mul(self):
        poly = ConvexPolytope([(0, 0), (1, 0), (0, 1)])
        larger = ConvexPolytope([(0, 0), (3, 0), (0, 3)])
        smaller = ConvexPolytope([(0, 0), (0.5, 0), (0, 0.5)])

        self.assertEqual(poly*3, larger)
        self.assertEqual(poly*0.5, smaller)
        self.assertEqual(poly*1, poly)

    def test_neg(self):
        poly = ConvexPolytope([(0, 0), (1, 0), (0, 1)])
        negative = ConvexPolytope([(0, 0), (-1, 0), (0, -1)])

        self.assertEqual(-poly, negative)
        self.assertEqual(-(-poly), poly)

    def test_truediv(self):
        poly = ConvexPolytope([(0, 0), (2, 0), (1, 1)])
        smaller = ConvexPolytope([(0, 0), (1, 0), (0.5, 0.5)])

        self.assertEqual(poly/2, smaller)
        self.assertEqual(poly/1, poly)

    def test_shift(self):
        poly = ConvexPolytope([(0, 0), (1, 0), (0, 1)])
        shifted = ConvexPolytope([(1, 1), (2, 1), (1, 2)])

        self.assertEqual(poly.shift(np.array([1, 1])), shifted)
        self.assertEqual(poly.shift(np.array([0, 0])), poly)


if __name__ == "__main__":
    main()
