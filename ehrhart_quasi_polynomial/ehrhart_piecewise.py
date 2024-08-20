from itertools import combinations_with_replacement

from .ehrhart_quasi_polynomial import ehrhart_quasi_polynomial

from sage.arith.functions import lcm
from sage.functions.other import factorial
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
    """
    Compute the piecewise Ehrhart quasi-polynomials of the polytopes defined by
        Ax <= b
    where the inequality is understood componentwise
    on the chambers of the secondary fan associated with the matrix ``A``.

    """
    def __init__(self, A):
        self._A = A
        self._sec_fan = self._secondary_fan()

        self._projection_generator()
        self._max_cones, self._qps, self._R, self._S = self._compute_piecewise()
        self._str_var = [str(x) for x in self._R.gens()]
        self._amb_dim = self._sec_fan.ambient_dim()

    def _secondary_fan(self):
        """
        Returns the secondary fan associated with ``self._A``
        For details see the documentation of the gfan software package:
            https://users-math.au.dk/~jensen/software/gfan/gfan.html
        """
        gfan_input = "{" + ", ".join(str(row) for row in self._A.rows()) + "}"
        return PolyhedralFan(gfan(gfan_input, "secondaryfan"))

    def _compute_piecewise(self):
        num_variables = self._A.nrows()
        max_degree = self._sec_fan.lineality_dim()

        R = PolynomialRing(QQ, "x", num_variables)
        x_vars = R.gens()
        S = PolynomialRing(R, "t")
        t = S.gen()

        fan_rays = self._sec_fan.rays()

        max_cones = []
        quasi_polynomials = []
        for ray_lists in self._sec_fan.maximal_cones().values():
            for ray_list in ray_lists:
                cone = Cone([fan_rays[idx] for idx in ray_list])
                period = self._estimate_period(cone)
                max_cones.append(cone)

                red_cone, dimensions_removed = self._remove_zero_dimensions(cone.rays())
                actual_num_variables = num_variables - len(dimensions_removed)

                needed_points = factorial(actual_num_variables + max_degree)//(
                    factorial(actual_num_variables)*factorial(max_degree) )*period
                cone_points = self._generate_cone_points(cone, needed_points)
                term_dict = self._get_term_dict(self._A, cone_points)

                interpolation_variables = [x for dim, x in enumerate(x_vars)
                                           if dim not in dimensions_removed]
                T = PolynomialRing(QQ, names=interpolation_variables)
                cone_points = [[p for dim, p in enumerate(point)
                                if dim not in dimensions_removed]
                               for point in cone_points]

                cone_poly = 0
                for degree, terms in term_dict.items():
                    cone_poly += T.interpolation(degree, cone_points, terms)*t**degree

                quasi_polynomials.append(cone_poly)
        return max_cones, quasi_polynomials, R, S

    def _estimate_period(self, chamber):
        # TODO: find better estimate for each chamber
        # since if the polytope is empty in a chamber we could skip all calculations
        return lcm(a for a in self._A.list() if a)

    def _remove_zero_dimensions(self, rays):
        zero_dimensions = []
        for dim, val in enumerate(rays[0]):
            if val == 0:
                if all(ray[dim] == 0 for ray in rays[1:]):
                    zero_dimensions.append(dim)
        reduced_cone = Cone([[val for dim, val in enumerate(ray)
                              if dim not in zero_dimensions ]for ray in rays])
        return reduced_cone, zero_dimensions

    def _generate_cone_points(self, cone, number):
        vectors = [free_module_element(list(ray)) for ray in cone.rays()]

        points = []
        points_len = 0
        index = 1
        while points_len <= number:
            new_points = [sum(combi) for combi
                          in combinations_with_replacement(vectors, index)]
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

    def _projection_generator(self):
        orth_vectors = []
        for vec in self._sec_fan.fan_dict["ORTH_LINEALITY_SPACE"]:
            orth_vectors.append(free_module_element([int(val) for val in vec.split(" ")]))

        lin_vectors = []
        for vec in self._sec_fan.fan_dict["LINEALITY_SPACE"]:
            lin_vectors.append(free_module_element([int(val) for val in vec.split(" ")]))

        change_of_basis_matrix = create_matrix(orth_vectors + lin_vectors)
        inverse = change_of_basis_matrix.inverse()

        def projection(point):
            new_representation = free_module_element(point)*inverse
            eval_point = sum(new_representation[k]*o_vec for k, o_vec in enumerate(orth_vectors))
            return eval_point

        self.projection = projection


    def __call__(self, point):
        return self.evaluate(point)

    def evaluate(self, point):
        proj = self.projection(point)
        for idx, cone in enumerate(self._max_cones):
            if proj in cone:
                keywords = {x: proj[k] for k, x in enumerate(self._str_var)}
                return self._qps[idx](**keywords)
        return 0

    def matrix(self):
        return self._A

    def secondary_fan(self):
        return self._sec_fan

    def domains(self):
        return self._max_cones

    def quasi_polynomials(self):
        return self._qps


def create_polytope_from_matrix(A, b):
    """
    Returns the polytope whose faces are defined by
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
        raise ValueError("Dimensions of 'A' and 'b' need to be compatible, need "+
                         f"A.nrows == len(b), but have {A.nrows()} != {len(b)}")
    inequalities = [[b[k]] + list(-A.rows()[k]) for k in range(A.nrows())]
    return Polyhedron(ieqs = inequalities)
