import numpy as np

from unittest import TestCase, main

class BoundingBox:
    def __init__(self, min_values, max_values):
        if not isinstance(min_values, np.ndarray):
            raise TypeError("min_values and max_values need to be of type numpy.ndarray. "+
                            "but are of type {type(min_values)} and {type(max_values)}")
        if np.shape(min_values) != np.shape(max_values):
            raise ValueError("min_values and max_values need to be of the same shape")
        if (len(np.shape(min_values)) != 1):
            raise ValueError("min_values and max_values need to be 1-dimensional. " +
                             "but have shapes {np.shape(min_values)} and {np.shape(max_values)}")
        if not isinstance(min_values[0], np.int64) or not isinstance(max_values[0], np.int64):
            raise ValueError("dtype of min_values and max_values needs to be a subclass of " +
                             f"numpy.int64, but are {type(min_values[0])} and {type(max_values[0])}")
        self.dim = len(min_values)

        self.min, self.max = self._check_min_max(min_values, max_values)

        self.points = self._generate_points()
        self.index = -1

        self.vertices = self._generate_vertices()

    """
    functions called in __init__
    """
    def _check_min_max(self, mins, maxs):
        if all(mins <= maxs):
            return mins, maxs

        combined = np.array([mins, maxs])
        return np.min(combined, axis=0), np.max(combined, axis=0)

    def _generate_points(self):
        prev_layer = [(x,) for x in range(self.min[0], self.max[0]+1)]
        for d in range(1, self.dim):
            next_layer = []
            for point in prev_layer:
                for y in range(self.min[d], self.max[d]+1):
                    next_layer.append(point + (y,))
            prev_layer = next_layer
        return next_layer

    def _generate_vertices(self):
        prev_vertices = [(self.min[0],), (self.max[0],)]
        for d in range(1, self.dim):
            next_vertices = []
            for corner in prev_vertices:
                next_vertices.extend([corner + (self.min[d],),
                                     corner + (self.max[d],)])
            prev_vertices = next_vertices
        return next_vertices

    """
    basically list methods to get the points
    """
    def __contains__(self, other):
        if not hasattr(other, "__len__"):
            return False
        if len(other) != self.dim:
            return False

        return all(self.min[k] <= other[k] <= self.max[k] 
                   for k in range(self.dim))

    def __len__(self):
        return len(self.points)

    def __next__(self):
        self.index += 1
        return self.points[self.index]

    """
    standard dunder methods
    """
    def __str__(self):
        return f"{self.dim}-dimensional rectangular cuboid bounded by:\n{self.vertices}"

    def __repr__(self):
        return f"BoundingBox({self.min}, {self.max})"

    """
    comparison operators
    """
    def comparable(self, other):
        return isinstance(other, BoundingBox) and self.dim == other.dim
    
    def compare(self, other):
        if not isinstance(other, type(self)):
            raise TypeError("other needs to be of type BoundingBox " +
                            f" but is of type {type(other)}")
        if self.dim != other.dim:
            raise ValueError("self and other need to have the same dimension " +
                             f"but have {self.dim} and {other.dim} respectively")

    def __eq__(self, other):
        return self.comparable(other) and \
            all(other.min == self.min) and all(other.max == self.max)

    def __le__(self, other):
        """
        Returns whether self is a subset of other.
        """
        return self.comparable(other) and \
            all(other.min <= self.min) and all(self.max <= other.max)

    def __lt__(self, other):
        """
        Returns whether self is a true subset of other
        """
        return self.comparable(other) and \
            all(other.min < self.min) and all(self.max < other.max)

    """
    math support
    """
    def __mul__(self, number):
        return BoundingBox(number*self.min, number*self.max)

    __rmul__ = __mul__

    def __neg__(self):
        return self * (-1)

    def __truediv__(self, number):
        return BoundingBox(self.min/number, self.max/number)

    def __add__(self, other):
        """
        Returns the smallest BoundingBox which contains both self and other
        """
        self.compare(other)
        new_min = np.min(np.array([self.min, other.min]), axis=0)
        new_max = np.max(np.array([self.max, other.max]), axis=0)
        return BoundingBox(new_min, new_max)

    def intersection(self, other):
        """
        Returns the largest BoundingBox contained in both self and other
        """
        self.compare(other)
        new_min = np.max(np.array([self.min, other.min]), axis=0)
        new_max = np.min(np.array([self.max, other.max]), axis=0)
        return BoundingBox(new_min, new_max)

    def union(self, other):
        """
        Returns a set which contains all points of self and other
        """
        if not hasattr(other, "points"):
            raise AttributeError("other needs to have the attribute 'points'")
        return set(self.points + other.points)

    def shift(self, point):
        """
        type safe ersion of _shift
        """
        if not isinstance(point, np.ndarray) or np.shape(point) != (self.dim,):
            return ValueError("point needs to be a 1-dimensional " +
                              f"numpy.ndarray with a width of {self.dim}")
        return BoundingBox(self.min + point, self.max + point)

class TestBoundingBox(TestCase):
    """
    points and vertices
    """
    def test_unit_square(self):
        mins = np.array([0, 0])
        maxs = np.array([1, 1])
        unit_square = BoundingBox(mins, maxs)

        self.assertEqual(unit_square.points, [(0, 0), (0, 1), (1, 0), (1, 1)])

    def test_x_axis(self):
        mins = np.array([0, 0])
        maxs = np.array([3, 0])
        x_axis = BoundingBox(mins, maxs)

        expected_set = set([point for point in x_axis.points])
        self.assertEqual(expected_set, set([(0, 0), (1, 0), (2, 0), (3, 0)]))

    def test_unit_cube_inverted(self):
        maxs = np.zeros(3, dtype=int)
        mins = np.ones(3, dtype=int)
        unit_cube = BoundingBox(mins, maxs)

        expected = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1),
                    (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)]
        self.assertEqual(unit_cube.points, expected)

    def test_vertices_unit_cube(self):
        mins = np.zeros(3, dtype=int)
        maxs = np.ones(3, dtype=int)
        box = BoundingBox(mins, maxs)

        self.assertEqual(box.points, box.vertices)

    """
    list methods
    """
    def test_contains(self):
        mins = np.zeros(3, dtype=int)
        maxs = np.ones(3, dtype=int)
        box = BoundingBox(mins, maxs)

        self.assertTrue((0, 0, 0) in box)
        self.assertFalse(0 in box)
        self.assertFalse((2, 2, 2) in box)
        self.assertTrue((1, 0, 0.5) in box)

    def test_len(self):
        mins = np.zeros(3, dtype=int)
        maxs = np.ones(3, dtype=int)
        box = BoundingBox(mins, maxs)
        
        self.assertEqual(len(box), 2**3)

    def test_next(self):
        mins = np.array([-2, 2, 0, 5])
        maxs = np.array([-1, 3, 1, 5])
        box = BoundingBox(mins, maxs)

        actual = []
        for k in range(5):
            point = next(box)
            actual.append(point)
        expected = [(-2, 2, 0, 5), (-2, 2, 1, 5),
                    (-2, 3, 0, 5), (-2, 3, 1, 5),
                    (-1, 2, 0, 5)]
        self.assertEqual(actual, expected)

    """
    comparison operators
    """
    def test_comparable_false(self):
        mins = np.zeros(3, dtype=int)
        maxs = np.ones(3, dtype=int)
        wrong_min = np.zeros(4, dtype=int)
        wrong_max = np.ones(4, dtype=int)

        box1 = BoundingBox(mins, maxs)
        box2 = BoundingBox(wrong_min, wrong_max)

        self.assertFalse(box1.comparable("box2"))
        self.assertFalse(box1.comparable(box2))

    def test_comparable_true(self):
        mins = np.zeros(3, dtype=int)
        maxs = np.ones(3, dtype=int)
        max2 = np.ones(3, dtype=int)*2

        box1 = BoundingBox(mins, maxs)
        box2 = BoundingBox(mins, mins)
        box3 = BoundingBox(mins, max2)

        self.assertTrue(box1.comparable(box2))
        self.assertTrue(box3.comparable(box1))

    def test_compare_error(self):
        mins = np.zeros(3, dtype=int)
        maxs = np.ones(3, dtype=int)
        cube = BoundingBox(mins, maxs)
        
        self.assertRaises(TypeError, cube.compare, "string")
        
        wrong_min = np.zeros(4, dtype=int)
        wrong_max = np.ones(4, dtype=int)
        wrong_cube = BoundingBox(wrong_min, wrong_max)
        
        self.assertRaises(ValueError, cube.compare, wrong_cube)

    def test_comparisons_true(self):
        mins = np.zeros(3, dtype=int)
        maxs = np.ones(3, dtype=int)
        max3 = np.ones(3, dtype=int)*3

        box1 = BoundingBox(mins, maxs)
        box2 = BoundingBox(mins, maxs)
        box3 = BoundingBox(-maxs, max3)

        self.assertTrue(box1 == box2)
        self.assertTrue(box1 <= box2)
        self.assertTrue(box1 >= box2)
        self.assertTrue(box1 < box3)
        self.assertTrue(box3 > box2)

    def test_comparisons_false(self):
        mins = np.zeros(3, dtype=int)
        maxs = np.ones(3, dtype=int)

        other_mins = -np.ones(3, dtype=int)
        other_maxs = np.zeros(3, dtype=int)
        box1 = BoundingBox(mins, maxs)
        box2 = BoundingBox(other_mins, other_maxs)

        self.assertFalse(box1 == box2)
        self.assertFalse(box1 <= box2)
        self.assertFalse(box1 >= box2)
        self.assertFalse(box1 < box1)
        self.assertFalse(box1 > box2)

    """
    math operations
    """
    def test_mul(self):
        mins = -np.ones(3, dtype=int)
        maxs = np.ones(3, dtype=int)
        box = BoundingBox(mins, maxs)

        mul_mins = -np.ones(3, dtype=int)*2
        mul_maxs = np.ones(3, dtype=int)*2
        mul_box = BoundingBox(mul_mins, mul_maxs)

        self.assertEqual(2*box, mul_box)
        self.assertEqual((-1)*box, box)

    def test_neg(self):
        mins = -np.ones(3, dtype=int)
        maxs = np.ones(3, dtype=int)
        box = BoundingBox(mins, maxs)

        self.assertEqual(-box, box)

        origin = np.zeros(3, dtype=int)
        positive = BoundingBox(origin, maxs)
        negative = BoundingBox(mins, origin)

        self.assertEqual(-negative, positive)

    def test_add(self):
        mins = np.zeros(3, dtype=int)
        maxs = np.ones(3, dtype=int)*2
        box = BoundingBox(mins, maxs)

        other_mins = -np.ones(3, dtype=int)*2
        other_maxs = np.ones(3, dtype=int)
        other_box = BoundingBox(other_mins, other_maxs)

        result = BoundingBox(other_mins, maxs)
        self.assertEqual(box + other_box, result)
        self.assertEqual(other_box.__add__(box), result)

    def test_union(self):
        mins = np.zeros(2, dtype=int)
        maxs = np.ones(2, dtype=int)*2
        box = BoundingBox(mins, maxs)

        other_mins = -np.ones(2, dtype=int)*2
        other_maxs = np.ones(2, dtype=int)
        other_box = BoundingBox(other_mins, other_maxs)

        union = set(box.points + other_box.points)
        self.assertEqual(box.union(other_box), union)
        self.assertEqual(other_box.union(box), union)

    def test_intersection(self):
        mins = np.array([0, -1], dtype=int)
        maxs = np.ones(2, dtype=int)
        box = BoundingBox(mins, maxs)

        other_mins = np.array([-1, 0], dtype=int)
        other_maxs = np.array([0, 2], dtype=int)
        other_box = BoundingBox(other_mins, other_maxs)

        intersection = BoundingBox(np.zeros(2, dtype=int), np.array([0, 1], dtype=int))
        self.assertEqual(box.intersection(other_box), intersection)
        self.assertEqual(other_box.intersection(box), intersection)

    def test_shift(self):
        mins = np.zeros(3, dtype=int)
        maxs = np.ones(3, dtype=int)
        larger_maxs = np.ones(3, dtype=int)*2
        
        poly = BoundingBox(mins, maxs)
        shifted = BoundingBox(maxs, larger_maxs)
        
        self.assertEqual(poly.shift(maxs), shifted)
        self.assertEqual(shifted.shift(-maxs), poly)

if __name__ == "__main__":
   main()
