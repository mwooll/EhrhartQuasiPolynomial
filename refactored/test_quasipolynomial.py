from quasipolynomial import QuasiPolynomialRing
from integerperiodicfunction import IntegerPeriodicFunctionRing

from sage.all import *
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from unittest import TestCase, main


class TestQuasiPolynomial(TestCase):
    def setUp(self):
        self.QPR = QuasiPolynomialRing(QQ)
        self.IPFR = IntegerPeriodicFunctionRing(QQ)

        self.poly = self.QPR([1, 2, 3])
        self.zero = self.QPR([0])

        self.qp = self.QPR([1, [2, 3], [1, 2, 3]])
        self.ipf = self.IPFR([1, 0])

    # __init__
    def test_zero(self):
        self.assertEqual(self.QPR(), self.QPR([0]))
        self.assertEqual(self.QPR.zero(), self.QPR([0]))

    def test_reduce_coefs(self):
        self.assertEqual(self.zero, self.QPR([0, 0, 0, 0]))
        self.assertEqual(self.poly, self.QPR([1, 2, 3, 0, 0]))
        self.assertEqual(self.QPR([1, 0]).degree(), 0)

    def test_calculate_period(self):
        self.assertEqual(self.zero.period(), 1)
        self.assertEqual(self.poly.period(), 1)
        self.assertEqual(self.qp.period(), 6)

    # call
    def test_call(self):
        self.assertEqual(self.zero(1), 0)
        self.assertEqual(self.poly(2), 17)

        self.assertEqual(self.qp(1), 6)
        self.assertEqual(self.qp(3), 19)

    # dunder methods
    def test_repr(self):
        pass

    def test_str(self):
        self.assertEqual(str(self.poly),
                          "QuasiPolynomial given by \n" +
                          "[1] + [2]*k + [3]*k^2")
        self.assertEqual(str(self.qp),
                         "QuasiPolynomial given by \n"
                         + "[1] + [2, 3]*k + [1, 2, 3]*k^2")

    # math support
    def test_neg(self):
        self.assertEqual(-self.poly, self.QPR([-1, -2, -3]))
        self.assertEqual(-self.zero, self.zero)

    def test_add(self):
        self.assertEqual(self.poly + self.zero, self.poly)
        self.assertEqual(self.zero + self.poly, self.poly)
        self.assertEqual(self.zero + self.zero, self.zero)

        self.assertEqual(self.poly + 1, self.QPR([2, 2, 3]))
        self.assertEqual(2 + self.poly, self.QPR([3, 2, 3]))

        self.assertEqual(self.qp + 4, self.QPR([5, [2, 3], [1, 2, 3]]))

        result = self.QPR([2, [4, 5], [4, 5, 6]])
        self.assertEqual(self.qp + self.poly, result)
        self.assertEqual(self.poly + self.qp, result)

        summation = self.QPR([[2, 1], 2, 3])
        self.assertEqual(self.poly + self.ipf, summation)
        self.assertEqual(self.ipf + self.poly, summation)

    def test_sub(self):
        self.assertEqual(self.poly - self.zero, self.poly)
        self.assertEqual(self.zero - self.poly, -self.poly)
        self.assertEqual(self.poly - self.poly, self.zero)

        self.assertEqual(self.poly - 1, self.QPR([0, 2, 3]))
        self.assertEqual(2 - self.poly, self.QPR([1, -2, -3]))

    def test_mul(self):
        self.assertEqual(self.poly*0, self.zero)
        self.assertEqual(self.poly*1, self.poly)
        self.assertEqual(self.poly*(-2), self.QPR([-2, -4, -6]))
        self.assertEqual(self.poly*QQ(1/2), self.QPR([1/2, 1, 3/2]))

        self.assertEqual(self.poly * self.zero, self.zero)
        self.assertEqual(self.poly*self.poly, self.QPR([1, 4, 10, 12, 9]))

        product = self.QPR([1, [4, 5], [8, 11, 10, 10, 9, 12],
                            [8, 13, 12, 11, 10, 15], [3, 6, 9]])
        self.assertEqual(self.poly*self.qp, product)
        self.assertEqual(self.qp*self.poly, product)

        result = self.QPR([[1, 0], [2, 0], [3, 0]])
        self.assertEqual(self.poly*self.ipf, result)
        self.assertEqual(self.ipf*self.poly, result)



if __name__ == "__main__":
    main(failfast=True)   