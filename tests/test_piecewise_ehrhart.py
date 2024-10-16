from itertools import combinations_with_replacement

from ehrhart_quasi_polynomial.ehrhart_piecewise import (
    PiecewiseEhrhartQuasiPolynomial as PEQP,
    secondary_fan,
    create_polytope_from_matrix,
    _compute_change_of_basis_matrix,
    _compute_periods,
    _generate_cone_points,
    _move_projected_point_inside_cone
    )

from sage.geometry.cone import Cone
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import free_module_element

from unittest import TestCase, main


class TestEhrhartPiecewise(TestCase):
    @classmethod
    def setUpClass(cls):
        
        class PEQP_Proxy():
            def __init__(self):
                self._A = Matrix([[-1, 0], [0, -1], [1, 1], [0, 1]])

        cls.peqp = PEQP_Proxy()

    def tearDown(self):
        pass

# =============================================================================
#     def test_compute_change_of_basis_matrix(self):
#         orth = tuple(free_module_element(orth) for orth
#                       in [(1, 0, 1, -1), (0, 1, 0, 1)])
#         lin = tuple(free_module_element(lin) for lin
#                     in [(1, 0, -1, 0), (0, 1, -1, -1)])
#         matrices = _compute_change_of_basis_matrix(orth, lin)
#         R4_to_ol, R4_to_o, _ = matrices
#         # expected = Matrix([])
#         # self.assertEqual(expected, actual)
# =============================================================================


    def test_process_secondary_fan(self):
        PEQP._process_secondary_fan(self.peqp)

        self.assertEqual(self.peqp._amb_dim, 4)

        expected_rays = [free_module_element(ray) for ray
                         in [(-1, 2, -1, 3), (1, 3, 1, 2), (2, 1, 2, -1)]]
        self.assertEqual(self.peqp._rays, expected_rays)

        expected_max_cones = {4: [[0, 1], [1, 2]]}
        self.assertEqual(self.peqp._maximal_cones, expected_max_cones)

        expected_orth = tuple(free_module_element(orth) for orth
                              in [(1, 0, 1, -1), (0, 1, 0, 1)])
        self.assertEqual(self.peqp._orth_vectors, expected_orth)

        expected_lin = tuple(free_module_element(lin) for lin
                             in [(1, 0, -1, 0), (0, 1, -1, -1)])
        self.assertEqual(self.peqp._lin_vectors, expected_lin)

        # self.assertEqual(self.change_of_basis_matrix, [])

# =============================================================================
#     def test_generate_cone_dicts(self):
#         PEQP._process_secondary_fan(self.peqp)
#         PEQP._generate_cone_dicts(self.peqp)
#         cone_dicts = self.peqp._cone_dicts
# 
#         self.assertEqual(first, second)
# =============================================================================
        

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

# =============================================================================
#     def test_move_projected_point_inside_cone(self):
#         pass
# =============================================================================

if __name__ == "__main__":
    main()
