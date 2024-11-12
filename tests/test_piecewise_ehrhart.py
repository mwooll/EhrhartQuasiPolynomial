from itertools import combinations_with_replacement

from ehrhart_quasi_polynomial.ehrhart_piecewise import (
    PiecewiseEhrhartQuasiPolynomial as PEQP,
    secondary_fan,
    create_polytope_from_matrix,
    _compute_change_of_basis_matrices,
    _compute_periods,
    _compute_off_sets,
    _generate_cone_points,
    _nudge_off_set
    )

from sage.geometry.cone import Cone
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import free_module_element
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ

from unittest import TestCase, main


class TestEhrhartPiecewise(TestCase):
    @classmethod
    def setUpClass(cls):
        
        class PEQP_Proxy():
            _A = Matrix([[-1, 0], [0, -1], [1, 1], [0, 1]])
            _create_polytope_from_matrix = PEQP._create_polytope_from_matrix
            _R = PolynomialRing(QQ, 2, "x")

        cls.peqp = PEQP_Proxy()

    def tearDown(self):
        pass

    def test_compute_change_of_basis_matrices(self):
        orth = tuple(free_module_element(orth) for orth
                      in [(1, 0, 1, -1), (0, 1, 0, 1)])
        lin = tuple(free_module_element(lin) for lin
                    in [(1, 0, -1, 0), (0, 1, -1, -1)])
        K, R = _compute_change_of_basis_matrices(orth, lin)

        self.assertEqual(K*orth[0], free_module_element([1, 0]))
        self.assertEqual(K*orth[1], free_module_element([0, 1]))
        self.assertEqual(K*lin[0], free_module_element([0, 0]))
        self.assertEqual(K*lin[1], free_module_element([0, 0]))

        self.assertEqual(R*free_module_element([1, 0]), orth[0])
        self.assertEqual(R*free_module_element([0, 1]), orth[1])

    def test_process_secondary_fan(self):
        PEQP._process_secondary_fan(self.peqp)

        expected_rays = tuple(free_module_element(ray) for ray
                              in [(-1, 2, -1, 3), (1, 3, 1, 2), (2, 1, 2, -1)])
        self.assertEqual(self.peqp._rays, expected_rays)

        expected_orth = tuple(free_module_element(orth) for orth
                              in [(1, 0, 1, -1), (0, 1, 0, 1)])
        self.assertEqual(self.peqp._orth_vectors, expected_orth)

        expected_lin = tuple(free_module_element(lin) for lin
                             in [(1, 0, -1, 0), (0, 1, -1, -1)])
        self.assertEqual(self.peqp._lin_vectors, expected_lin)

    def test_generate_cone_dicts(self):
        PEQP._process_secondary_fan(self.peqp)
        cone_dicts = PEQP._generate_cone_dicts(self.peqp)

        self.assertEqual(len(cone_dicts), 2)
        for key in ["rays", "cone", "scaled_rays", "off_sets"]:
            self.assertTrue(key in cone_dicts[0])

    def test_create_polytope_from_matrix(self):
        A = Matrix([[-1, 0], [0, -1], [1, 1]])
        b = (0, 0, 1)
        actual = create_polytope_from_matrix(A, b)
        expected = Polyhedron([[0, 0], [0, 1], [1, 0]])
        self.assertEqual(expected, actual)

        A = Matrix([[-1, 0], [0, -1], [1, 1], [0, 1]])
        b = (0, 0, 5, 3)
        actual = create_polytope_from_matrix(A, b)
        expected = Polyhedron([[0, 0], [5, 0], [2, 3], [0, 3]])
        self.assertEqual(expected, actual)

        actual_unsafe = PEQP._create_polytope_from_matrix(self.peqp, b)
        self.assertEqual(expected, actual_unsafe)

    def test_compute_periods(self):
        A = Matrix([[-1, 0], [0, -1], [2, 3]])
        points = [(0, 0, k) for k in range(7)]
        actual = _compute_periods(A, points)
        expected = [1, 6, 3, 2, 3, 6, 1]
        self.assertEqual(expected, actual)

    def test_generate_cone_points(self):
        rays = [free_module_element(b)
                for b in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
        actual = _generate_cone_points(rays, 10)
        expected = [free_module_element(point) for point in
                    [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
                    (2, 0, 0), (1, 1, 0), (1, 0, 1),
                    (0, 2, 0), (0, 1, 1), (0, 0, 2)]]
        self.assertEqual(expected, actual)
        cone = Cone(rays)
        for point in actual:
            self.assertTrue(point in cone)

    def test_interpolate_off_set_poly(self):
        x1, x2 = self.peqp._R.gens()
        rays = tuple(free_module_element(orth) for orth
                     in [(-1, 2), (1, 3)])
        actual = PEQP._interpolate_off_set_poly(self.peqp,
                                                _generate_cone_points(rays, 6),
                                                2, 0*rays[0])
        expected = x1 + x2
        self.assertEqual(expected, actual)

    def test_compute_off_sets(self):
        orth = tuple(free_module_element(orth) for orth
                      in [(1, 0, 1, -1), (0, 1, 0, 1)])
        lin = tuple(free_module_element(lin) for lin
                    in [(1, 0, -1, 0), (0, 1, -1, -1)])
        M, I, K, R = _compute_change_of_basis_matrices(orth, lin)

        self.fail()

if __name__ == "__main__":
    # main()

    A = Matrix([[-1, 0], [0, -1], [1, 1]])
    p = PEQP(A)
    print(p._cone_dicts[0])

    num_int_points = lambda A, b: len(create_polytope_from_matrix(A, b).integral_points())

    def test_combinations(p, nrows):
        for b in combinations_with_replacement(range(-10, 11), nrows):
            expected = num_int_points(A, b)
            actual = p(b)
            if actual != expected:
                print(b, expected, actual)

    dim = A.nrows()
    # test_combinations(p, dim)

    def test_Q_off_sets():
        V = IntegralLattice(Matrix.identity(dim))
        W = V.sublattice(Matrix(p._orth_vectors + p._lin_vectors))
        Q = V.quotient(W)
        for a in Q:
            lift = a.lift() + sum(p._rays)
            print(a, lift, num_int_points(A, lift))

            print(p._investigate(lift))
            proj = p._projection_matrix*lift
            inv = p._projection_inverse*proj
            print(lift - inv, num_int_points(A, inv))
            print()

    test_Q_off_sets()
    