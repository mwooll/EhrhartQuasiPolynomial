from ehrhart_quasi_polynomial import ehrhart_quasi_polynomial, QuasiPolynomialRing

import sage.all
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.geometry.polyhedron.constructor import Polyhedron

from unittest import TestCase, main

R = PolynomialRing(QQ, "k")
k = R.gen()

def quadrilateral_formula(m, n):
    """
    This function is an explicit general formula to compute the Ehrhart polynomial
    of the polytope defined by the inequalities:
        0 <= x,
        0 <= y < n,
        x + y <= m
    parametrized by m and n.
    """
    if m < 0 or n < 0:
        return 0
    if m == 0:
        return 1
    if n == 0:
        return 1 + m*k
    if m <= n:
        return 1 + 3*m/2*k + m**2*k**2/2
    if m > n:
        return 1 + (m + n/2)*k + (m*n - n**2/2)*k**2


class TestQuadrilateral(TestCase):
    def create_polytope(self, m, n):
        """
        Generates the wanted quadrilateral.
        """
        return Polyhedron(ieqs = [[0, 1, 0], [0, 0, 1], [m, -1, -1], [n, 0, -1]])

    def compute_ehrhart_polynomial(self, m, n):
        polytope = self.create_polytope(m, n)
        
        if polytope.is_empty():
            return 0
        
        vertices = polytope.Vrepresentation()
        ehr_pol = ehrhart_quasi_polynomial(vertices)
        return ehr_pol

    def compare(self, m, n):
        self.assertEqual(self.compute_ehrhart_polynomial(m, n),
                         quadrilateral_formula(m, n))
        
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
    # qpr = QuasiPolynomialRing(QQ)
    # one = qpr([1])
    # print(one == R([1]))