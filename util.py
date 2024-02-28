import numpy as np

from unittest import TestCase, main

def scalar_check(number):
    try:
        float(number)

    except TypeError:
        raise TypeError("number needs to be convertible to float")

    except ValueError:
        raise ValueError("number needs to be convertible to float")

class TestUtil(TestCase):
    def test_scalar_check(self):
        self.assertRaises(TypeError, scalar_check, np.zeros(3))
        self.assertRaises(ValueError, scalar_check, "word")

        self.assertIsNone(scalar_check(3.0))
        self.assertIsNone(scalar_check(3))
        self.assertIsNone(scalar_check("3"))
        self.assertIsNone(scalar_check(np.ones(1)))
        
if __name__ == "__main__":
    main()