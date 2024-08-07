from ehrhart_quasi_polynomial import ehrhart_quasi_polynomial, QuasiPolynomialRing

import sage.all
from sage.arith.misc import GCD
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.matrix.constructor import Matrix

from unittest import TestCase, main

R = PolynomialRing(QQ, "k")
k = R.gen()

def triangle_formula(m, n):
    m = abs(m)
    n = abs(n)
    return 1 + (m + n + GCD(m, n))/2*k + m*n/2*k**2


class TestTriangleFormula(TestCase):
    def create_triangle(self, m, n):
        return Polyhedron([(0, 0), (m, 0), (0, n)])

    def compute_ehrhart_polynomial(self, m, n):
        triangle = self.create_triangle(m, n)

        if triangle.is_empty():
            return 0

        vertices = triangle.Vrepresentation()
        ehr_pol = ehrhart_quasi_polynomial(vertices)
        return ehr_pol

    def compare(self, m, n):
        self.assertEqual(self.compute_ehrhart_polynomial(m, n),
                         triangle_formula(m, n))

    def test_negatives(self):
        self.compare(-5, 2)
        self.compare(3, -2)
        self.compare(-4, -1)

    def test_m_larger(self):
        self.compare(5, 0)
        self.compare(3, 2)

    def test_m_smaller(self):
        self.compare(0, 7)
        self.compare(2, 3)

    def test_equal_mn(self):
        self.compare(0, 0)
        self.compare(1, 1)
        self.compare(10, 10)

if __name__ == "__main__":
    main()
