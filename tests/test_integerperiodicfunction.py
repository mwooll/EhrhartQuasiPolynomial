from ehrhart_polynomial.integerperiodicfunction import IntegerPeriodicFunctionRing

from sage.all import *
from sage.rings.rational_field import QQ

from unittest import TestCase, main

class TestIntegerPeriodicFunction(TestCase):
    def setUp(self):
        self.ipfr = IntegerPeriodicFunctionRing(QQ)
        self.ipf = self.ipfr([1, 2, 3])
        self.zero = self.ipfr.zero()

    def test_zero(self):
        self.assertEqual(self.ipfr(), self.ipfr([0]))
        self.assertEqual(self.ipfr(), self.zero)

    def test_calculate_period(self):
        self.assertEqual(self.ipf.period(), 3)
        self.assertEqual(self.ipfr([4, 4, 4, 4]).period(), 1)
        self.assertEqual(self.ipfr([1, 2, 1, 2, 1, 2]).period(), 2)

    def test_call(self):
        self.assertEqual(self.zero(1), 0)
        self.assertEqual(self.ipf(0), 1)
        self.assertEqual(self.ipf(3), 1)
        self.assertEqual(self.ipf(7), 2)

    def test_repr(self):
        self.assertEqual(self.ipf.__repr__(),
                         'IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [1, 2, 3])')

    def test_str(self):
        self.assertEqual(str(self.ipf),
                         "IntegerPeriodicFunction over Rational Field given by"
                         + "\n\t1 if k%3 == 0"
                         + "\n\t2 if k%3 == 1"
                         + "\n\t3 if k%3 == 2")


    def test_neg(self):
        self.assertEqual(-self.ipf, self.ipfr([-1, -2, -3]))

    def test_add(self):
        self.assertEqual(1 + self.ipf, self.ipfr([2, 3, 4]))
        self.assertEqual(self.ipf + QQ(1/2), self.ipfr([3/2, 5/2, 7/2]))

        self.assertEqual(self.ipf + self.zero, self.ipf)
        self.assertEqual(self.ipf + self.ipf, self.ipfr([2, 4, 6]))

        self.assertEqual(self.ipf + self.ipfr([0, 1]),
                         self.ipfr([1, 3, 3, 2, 2, 4]))

    def test_sub(self):
        self.assertEqual(self.ipf - 1, self.ipfr([0, 1, 2]))
        self.assertEqual(QQ(1/2) - self.ipf, self.ipfr([-1/2, -3/2, -5/2]))

        self.assertEqual(self.ipf - self.zero, self.ipf)
        self.assertEqual(self.ipf - self.ipf, self.zero)

        self.assertEqual(self.ipf - self.ipfr([0, 1]),
                         self.ipfr([1, 1, 3, 0, 2, 2]))

    def test_mul(self):
        self.assertEqual(self.ipf*1, self.ipf)
        self.assertEqual(0*self.ipf, self.zero)
        self.assertEqual(self.ipf*2, self.ipfr([2, 4, 6]))
        self.assertEqual(QQ(-4/3)*self.ipf, self.ipfr([-4/3, -8/3, -4]))

        self.assertEqual(self.ipf*self.ipf, self.ipfr([1, 4, 9]))
        self.assertEqual(self.ipf*self.ipfr([1, -1]),
                         self.ipfr([1, -2, 3, -1, 2, -3]))


if __name__ == "__main__":
    main(failfast=True)   