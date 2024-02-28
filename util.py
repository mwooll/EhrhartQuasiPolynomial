import numpy as np

from unittest import TestCase, main

def scalar_check(number):
    try:
        float(number)

    except TypeError:
        raise TypeError("number needs to be convertible to float")

    except ValueError:
        raise ValueError("number needs to be convertible to float")

def point_check(point, dimension):
    if dimension == 1:
        scalar_check(point)
        return

    if not isinstance(point, np.ndarray):
        raise TypeError("point needs to be a subclass numpy.ndarray " +
                        f"but is of type {type(point)}")
    if np.shape(point) != (dimension,):
        raise ValueError("point needs to have shape {(dimension,)} " +
                         f"but has shape {np.shape(point)}")

class TestUtil(TestCase):
    def test_scalar_check(self):
        self.assertRaises(TypeError, scalar_check, np.zeros(3))
        self.assertRaises(ValueError, scalar_check, "word")

        self.assertIsNone(scalar_check(3.0))
        self.assertIsNone(scalar_check(3))
        self.assertIsNone(scalar_check("3"))
        self.assertIsNone(scalar_check(np.ones(1)))

    def test_point_check(self):
        self.assertRaises(TypeError, point_check, 3, 2)
        self.assertRaises(ValueError, point_check, np.zeros(1), 2)

        self.assertIsNone(point_check(3, 1))
        self.assertIsNone(point_check(np.zeros(3), 3))
        
if __name__ == "__main__":
    main()