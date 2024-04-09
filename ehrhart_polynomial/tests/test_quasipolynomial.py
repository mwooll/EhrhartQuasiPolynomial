from quasipolynomial import QuasiPolynomial
from integerperiodicfunction import IntegerPeriodicFunction

from unittest import TestCase, main


class TestQuasiPolynomial(TestCase):
    def setUp(self):
        self.poly = QuasiPolynomial([1, 2, 3])
        self.zero = QuasiPolynomial([0])

        coefs = [1, IntegerPeriodicFunction([2, 3]), IntegerPeriodicFunction([1, 2, 3])]
        self.ipf = QuasiPolynomial(coefs)

    # __init__
    def test_default_constructor_parameters(self):
        self.assertEqual(QuasiPolynomial(), QuasiPolynomial([0]))

    def test_reduce_coefs(self):
        self.assertEqual(self.zero, QuasiPolynomial([0, 0, 0, 0]))
        self.assertEqual(self.poly, QuasiPolynomial([1, 2, 3, 0, 0]))
        self.assertEqual(QuasiPolynomial([1, 0]).degree, 0)

    def test_calculate_period(self):
        self.assertEqual(self.zero.period, 1)
        self.assertEqual(self.poly.period, 1)
        self.assertEqual(self.ipf.period, 6)

    # call
    def test_call(self):
        self.assertEqual(self.zero(1), 0)
        self.assertEqual(self.poly(2), 17)

        self.assertEqual(self.ipf(1), 6)
        self.assertEqual(self.ipf(3), 19)

    # dunder methods
    def test_repr(self):
        self.assertEqual(eval(repr(self.poly)), self.poly)
        self.assertEqual(eval(repr(self.ipf)), self.ipf)

    def test_str(self):
        self.assertEqual(str(self.poly),
                          "QuasiPolynomial given by \n" +
                          "1 + 2*k + 3*k^2")
        self.assertEqual(str(self.ipf),
                         "QuasiPolynomial given by \n"
                         + "1 + IntegerPeriodicFunction([2, 3])(k)*k + "
                         + "IntegerPeriodicFunction([1, 2, 3])(k)*k^2")

    # math support
    def test_neg(self):
        self.assertEqual(-self.poly, QuasiPolynomial([-1, -2, -3]))
        self.assertEqual(-self.zero, self.zero)

    def test_add(self):
        self.assertEqual(self.poly + self.zero, self.poly)
        self.assertEqual(self.zero + self.poly, self.poly)
        self.assertEqual(self.zero + self.zero, self.zero)

        self.assertEqual(self.poly + 1, QuasiPolynomial([2, 2, 3]))
        self.assertEqual(2 + self.poly, QuasiPolynomial([3, 2, 3]))

        self.assertEqual(self.ipf + 4,
                         QuasiPolynomial([5,
                                          IntegerPeriodicFunction([2, 3]),
                                          IntegerPeriodicFunction([1, 2, 3])
                                         ]))

        result = QuasiPolynomial([2,
                                  IntegerPeriodicFunction([4, 5]),
                                  IntegerPeriodicFunction([4, 5, 6]),
                                 ])
        self.assertEqual(self.ipf + self.poly, result)
        self.assertEqual(self.poly + self.ipf, result)

        summation = QuasiPolynomial([IntegerPeriodicFunction([1, 2]), 2, 3])
        self.assertEqual(self.poly + IntegerPeriodicFunction([0, 1]),
                         summation)
        self.assertEqual(IntegerPeriodicFunction([0, 1]) + self.poly,
                         summation)

    def test_sub(self):
        self.assertEqual(self.poly - self.zero, self.poly)
        self.assertEqual(self.zero - self.poly, -self.poly)
        self.assertEqual(self.poly - self.poly, self.zero)

        self.assertEqual(self.poly - 1, QuasiPolynomial([0, 2, 3]))
        self.assertEqual(2 - self.poly, QuasiPolynomial([1, -2, -3]))

    def test_mul(self):
        self.assertEqual(self.poly*0, self.zero)
        self.assertEqual(self.poly*1, self.poly)
        self.assertEqual(self.poly*(-2), QuasiPolynomial([-2, -4, -6]))
        self.assertEqual(self.poly*(1/2), QuasiPolynomial([1/2, 1, 3/2]))

        self.assertEqual(self.poly * self.zero, self.zero)
        self.assertEqual(self.poly*self.poly, QuasiPolynomial([1, 4, 10, 12, 9]))

        product = QuasiPolynomial([1,
                                   IntegerPeriodicFunction([4, 5]),
                                   IntegerPeriodicFunction([8, 11, 10, 10, 9, 12]),
                                   IntegerPeriodicFunction([8, 13, 12, 11, 10, 15]),
                                   IntegerPeriodicFunction([3, 6, 9])
                                  ])
        self.assertEqual(self.poly*self.ipf, product)
        self.assertEqual(self.ipf*self.poly, product)

        result = QuasiPolynomial([IntegerPeriodicFunction([0, 1]),
                                  IntegerPeriodicFunction([0, 2]),
                                  IntegerPeriodicFunction([0, 3])])
        self.assertEqual(self.poly*IntegerPeriodicFunction([0, 1]), result)
        self.assertEqual(IntegerPeriodicFunction([0, 1])*self.poly, result)

    def test_truediv(self):
        self.assertEqual(self.poly/2, QuasiPolynomial([1/2, 1, 3/2]))
        self.assertEqual(self.poly/(-3/2), QuasiPolynomial([-2/3, -4/3, -6/3]))
        self.assertRaises(NotImplementedError, self.poly.__rtruediv__, 2)
        self.assertRaises(NotImplementedError, self.poly.__rtruediv__, self.poly)


if __name__ == "__main__":
    main()   