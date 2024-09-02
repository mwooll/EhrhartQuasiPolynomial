from itertools import combinations, combinations_with_replacement

from .ehrhart_quasi_polynomial import (ehrhart_quasi_polynomial,
                                       get_period, get_gcd)

from sage.arith.functions import lcm
from sage.functions.other import factorial, floor
from sage.geometry.cone import Cone
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.interfaces.gfan import gfan
from sage.matrix.constructor import Matrix as create_matrix
from sage.modules.free_module_element import free_module_element
from sage.rings.polynomial.groebner_fan import PolyhedralFan
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.structure.element import Matrix


class PiecewiseEhrhartQuasiPolynomial():
    def __init__(self, A):
        self._A = A
        self._sec_fan = self._secondary_fan()
        self._amb_dim = self._sec_fan.ambient_dim()
        self._origin = free_module_element([0]*self._amb_dim)

        self._compute_change_of_basis_matrices()
        self._max_cones = self._get_max_cones()
        self._compute_offsets()
        # self._qps, self._R, self._S = self._compute_piecewise()
        # self._str_var = [str(x) for x in self._R.gens()]        

    def _secondary_fan(self):
        """
        Returns the secondary fan associated with ``self._A``
        For details see the documentation of the gfan software package:
            https://users-math.au.dk/~jensen/software/gfan/gfan.html
        """
        gfan_input = "{" + ", ".join(str(row) for row in self._A.rows()) + "}"
        return PolyhedralFan(gfan(gfan_input, "secondaryfan"))

    def _get_max_cones(self):
        fan_rays = self._sec_fan.rays()
        max_cones = []
        for ray_lists in self._sec_fan.maximal_cones().values():
            for ray_list in ray_lists:
                cone = Cone(fan_rays[idx] for idx in ray_list)
                max_cones.append(cone)
        return max_cones

    def _compute_piecewise(self):
        num_variables = self._A.nrows()
        max_degree = self._sec_fan.lineality_dim()

        R = PolynomialRing(QQ, "x", num_variables)
        x_vars = R.gens()
        
        # TODO: remove t
        S = PolynomialRing(R, "t")
        t = S.gen()

        quasi_polynomials = []
        for idx, cone in enumerate(self._max_cones):
            cone_rays = cone.rays()

            red_cone, dimensions_removed = self._remove_zero_dimensions(cone_rays)
            actual_num_variables = num_variables - len(dimensions_removed)

            needed_points = factorial(actual_num_variables + max_degree)//(
                factorial(actual_num_variables)*factorial(max_degree) )
            cone_points = self._generate_cone_points(cone_rays, needed_points)

            expected_degree = [max_degree if dim not in dimensions_removed else 0
                               for dim in range(self._amb_dim)]

            cone_dict = {}
            for off_set in self._off_sets[idx]:
                off_cone_points = [p + off_set for p in cone_points]
                term_dict = self._get_term_dict(self._A, off_cone_points)
                cone_poly = 0
                for degree, terms in term_dict.items():
                    interpolated = R.interpolation(expected_degree, cone_points, terms)
                    cone_poly += interpolated(*self._projection(x_vars))*t**degree
                cone_dict[tuple(off_set)] = cone_poly

            quasi_polynomials.append(cone_dict)

        return quasi_polynomials, R, S

    def _estimate_period(self, cone):
        # TODO: find better estimate for each cone
        # since if the polytope is empty in a cone we could skip all calculations
        return lcm(a for a in self._A.list() if a)

    def _compute_offsets(self):
        self._scaled_rays = []
        self._off_sets = []
        for cone in self._max_cones:
            scaled_rays = []
            for ray in cone.rays():
                den = get_period(self._create_polytope_from_matrix(ray).Vrepresentation())
                scaled_rays.append(den*ray)
            self._scaled_rays.append(scaled_rays)

            vertices = scaled_rays[:]
            for dim in range(2, len(scaled_rays)+1):
                combis = combinations(scaled_rays, dim)
                vertices += [sum(combi, self._origin) for combi in combis]

            polytope = Polyhedron([self._origin] + vertices)

            integral_points = polytope.integral_points()
            without_vertices = [p for p in integral_points if p not in vertices]
            off_set_points = [point for point in without_vertices
                              if not any(point - ray in integral_points
                                         for ray in scaled_rays)]
            self._off_sets.append(off_set_points)

    def _rat_period(self, ray):
        original_polytope = self._create_polytope_from_matrix(ray)
        period = get_period(original_polytope.Vrepresentation())
        scaled_polytope = self._create_polytope_from_matrix([period*r for r in ray])
        gcd = get_gcd(scaled_polytope.Vrepresentation())
        return period/gcd

    def _remove_zero_dimensions(self, rays):
        zero_dimensions = []
        for dim, val in enumerate(rays[0]):
            if val == 0:
                if all(ray[dim] == 0 for ray in rays[1:]):
                    zero_dimensions.append(dim)
        reduced_cone = Cone([[val for dim, val in enumerate(ray)
                              if dim not in zero_dimensions] for ray in rays])
        return reduced_cone, zero_dimensions

    def _generate_cone_points(self, cone_rays, number):
        points = []
        points_len = 0
        index = 1
        while points_len <= number:
            new_points = [sum(combi) for combi
                          in combinations_with_replacement(cone_rays, index)]
            points += new_points
            index += 1
            points_len += len(new_points)
        return points[:number]

    def _get_term_dict(self, A_matrix, points):
        polynomials = []
        for point in points:
            polytope = self._create_polytope_from_matrix(point)
            ehr_poly = ehrhart_quasi_polynomial(polytope.Vrepresentation())
            polynomials.append(ehr_poly)

        max_degree = max(poly.degree() for poly in polynomials)
        terms = {d: self._get_terms_of_order(polynomials, d)
                 for d in range(max_degree + 1)}
        return terms

    def _create_polytope_from_matrix(self, b):
        """
        Unsafe version of ``create_polytope_from_matrix``
        to increase performance during lengthy computations.
        """
        inequalities = [[b[k]] + list(-self._A.rows()[k]) for k in range(self._A.nrows())]
        return Polyhedron(ieqs = inequalities)

    def _get_terms_of_order(self, polynomials, order):
        return [poly.coefficients()[order].constants()[0]
                if poly.degree() >= order else 0 for poly in polynomials]

    def _compute_change_of_basis_matrices(self):
        orth_vectors = []
        for vec in self._sec_fan.fan_dict["ORTH_LINEALITY_SPACE"]:
            orth_vectors.append(free_module_element([int(val) for val in vec.split(" ")]))
        self._orth_vectors = orth_vectors
        self._orth_dim = len(orth_vectors)

        lin_vectors = []
        for vec in self._sec_fan.fan_dict["LINEALITY_SPACE"]:
            lin_vectors.append(free_module_element([int(val) for val in vec.split(" ")]))

        self.change_of_basis_matrix = create_matrix(orth_vectors + lin_vectors)
        self.change_of_basis_inverse = self.change_of_basis_matrix.inverse()

    def _projected_coordinates(self, point):
        return free_module_element(point)*self.change_of_basis_inverse

    def _projection(self, point):
        new_representation = self._projected_coordinates(point)
        eval_point = sum(new_representation[k]*o_vec for k, o_vec in enumerate(self._orth_vectors))
        return eval_point

    def _find_offset(self, point, cone_idx):
        projected = self._projected_coordinates(point)
        off_set = sum( (projected[k] - int(projected[k]))*o_vec
                      for k, o_vec in enumerate(self._scaled_rays[cone_idx]))
        return tuple(int(val) for val in off_set)

    def __call__(self, point):
        return self.evaluate(point)

    def evaluate(self, point):
        if len(point) != self._amb_dim:
            raise ValueError("Dimension of ``point`` needs to be equal to the ambient"
                             f" dimension of ``self`` which is {self._amb_dim}")

        proj = self._projection(point)
        for idx, cone in enumerate(self._max_cones):
            if proj in cone:
                keywords = {x: point[k] for k, x in enumerate(self._str_var)}

                off_set = self._find_offset(point, idx)
                return self._qps[idx][off_set](**keywords)
        return 0

    def __repr__(self):
        return f"PiecewiseEhrhartQuasiPolynomial({self._A})"

    def matrix(self):
        return self._A

    def secondary_fan(self):
        return self._sec_fan

    def domains(self):
        return self._max_cones

    def functions(self):
        return self._qps


def create_polytope_from_matrix(A, b):
    """
    Return the polytope whose faces are defined by
        Ax <= b
    where the inequality is understood componentwise

    EXAMPLES::

        sage: from ehrhart_quasi_polynomial import create_polytope_from_matrix
        sage: A = Matrix([[-1, 0], [0, -1], [1, 1]]); b = [0, 0, 1]
        sage: poly = create_polytope_from_matrix(A, b)
        sage: poly
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices
        sage: poly.Vrepresentation()
        (A vertex at (1, 0), A vertex at (0, 1), A vertex at (0, 0))
    """
    if not isinstance(A, Matrix):
        raise TypeError("A needs to be an instance of sage.structure.element"
                        f".Matrix, but has type {type(A)}")
    if not hasattr(b, "__getitem__"):
        raise AttributeError("'b' needs to implement '__getitem__', " +
                             "needed for indexing: b[0]")
    if A.nrows() != len(b):
        raise ValueError("Dimensions of 'A' and 'b' need to be compatible, need "
                         f"A.nrows == len(b), but have {A.nrows()} != {len(b)}")
    inequalities = [[b[k]] + list(-A.rows()[k]) for k in range(A.nrows())]
    return Polyhedron(ieqs = inequalities)
