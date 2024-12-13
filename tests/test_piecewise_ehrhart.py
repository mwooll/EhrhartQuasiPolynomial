from ehrhart_quasi_polynomial.ehrhart_piecewise import (
    PiecewiseEhrhartQuasiPolynomial as PEQP,
    create_polytope_from_matrix,
    secondary_fan,
    _process_fan_vectors,
    _compute_change_of_basis_matrices,
    _hat_denominator,
    _generate_cone_points)

from sage.geometry.polyhedron.constructor import Polyhedron
from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import free_module_element
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


    def test_compute_change_of_basis_matrices(self):
        rays = tuple(free_module_element(ray) for ray
                      in [(1, 0, 1, -1), (0, 1, 0, 1)])
        lin = tuple(free_module_element(lin) for lin
                    in [(1, 0, -1, 0), (0, 1, -1, -1)])
        K, M = _compute_change_of_basis_matrices(rays, lin)

        self.assertEqual(K*free_module_element([1, 0, 0, 0]), rays[0])
        self.assertEqual(K*free_module_element([0, 0, 1, 0]), lin[0])

        self.assertEqual(M*rays[0], free_module_element([1, 0, 0, 0]))
        self.assertEqual(M*(lin[0] + lin[1]), free_module_element([0, 0, 1, 1]))
    

    def test_process_secondary_fan(self):
        PEQP._process_secondary_fan(self.peqp)

        expected_rays = tuple(free_module_element(ray) for ray
                              in [(-1, 2, -1, 3), (1, 3, 1, 2), (2, 1, 2, -1)])
        self.assertEqual(self.peqp._rays, expected_rays)

        expected_lin = tuple(free_module_element(lin) for lin
                             in [(1, 0, -1, 0), (0, 1, -1, -1)])
        self.assertEqual(self.peqp._lin_vectors, expected_lin)

    def test_process_fan_vectors(self):
        sec_fan = secondary_fan(self.peqp._A)
        expected = tuple(free_module_element(lin) for lin
                    in [(1, 0, -1, 0), (0, 1, -1, -1)])
        actual = _process_fan_vectors(sec_fan.fan_dict["LINEALITY_SPACE"])
        self.assertEqual(expected, actual)

    def test_generate_cone_dicts(self):
        PEQP._process_secondary_fan(self.peqp)
        cone_dicts = PEQP._generate_cone_dicts(self.peqp)

        self.assertEqual(len(cone_dicts), 2)
        for key in ["cone", "scaled_rays", "quotient"]:
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

    def test_simplex_points(self):
        actual = _generate_cone_points(3, 2)
        expected = tuple(free_module_element(point) for point in
                    [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
                    (2, 0, 0), (1, 1, 0), (1, 0, 1),
                    (0, 2, 0), (0, 1, 1), (0, 0, 2)])
        self.assertEqual(expected, actual)


if __name__ == "__main__":
    main()

    # A = Matrix([[-1, 0], [0, -1], [1, 1]])
    # PEQP(A)