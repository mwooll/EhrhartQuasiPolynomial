from itertools import combinations_with_replacement

from .ehrhart_quasi_polynomial import get_period

from sage.functions.other import ceil, factorial
from sage.geometry.cone import Cone
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.interfaces.gfan import gfan
from sage.matrix.constructor import Matrix as create_matrix
from sage.modules.free_module import span
from sage.modules.free_module_element import free_module_element
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.groebner_fan import PolyhedralFan
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.structure.element import Matrix


class PiecewiseEhrhartQuasiPolynomial():
    def __init__(self, A):
        self._A = A
        self._process_secondary_fan()

        self._cone_dicts = self._generate_cone_dicts()

        self._compute_piecewise()

    def _process_secondary_fan(self):
        sec_fan = secondary_fan(self._A)
        
        self._fan_dim = sec_fan.dim()
        self._amb_dim = sec_fan.ambient_dim()
        self._rays = tuple(free_module_element(ray) for ray in sec_fan.rays())
        self._scalars = _compute_periods(self._A, self._rays)
        self._maximal_cones = sec_fan.maximal_cones()

        self._orth_vectors = tuple(free_module_element([int(val) for val in vec.split(" ")])
                              for vec in sec_fan.fan_dict["ORTH_LINEALITY_SPACE"])
        self._lin_vectors = tuple(free_module_element([int(val) for val in vec.split(" ")])
                              for vec in sec_fan.fan_dict["LINEALITY_SPACE"])
        M, I, K, R = _compute_change_of_basis_matrices(self._orth_vectors, self._lin_vectors)
        self._M = M
        self._I = I
        self._projection_matrix = K
        self._projection_inverse = R

        self._projected_rays = tuple(K*ray for ray in self._rays)

        self._num_variables = len(self._orth_vectors)
        self._R = PolynomialRing(QQ, "x", self._num_variables)

    def _generate_cone_dicts(self):
        dictionaries = []
        for ray_lists in  self._maximal_cones.values():
            for ray_list in ray_lists:
                scalars = [self._scalars[idx] for idx in ray_list]

                cone_dic = {}
                cone_dic["rays"] = [self._projected_rays[idx] for idx in ray_list]
                cone_dic["cone"] = Cone(cone_dic["rays"] + self._lin_vectors) # changed
                cone_dic["scaled_rays"] = [scale*ray
                                           for scale, ray in zip(scalars, cone_dic["rays"])]

                # cone_dic["quotient"] = _get_off_sets(self)
                cone_dic["quotient"] = _compute_off_sets(self._projection_matrix,
                                                         cone_dic["scaled_rays"])

                dictionaries.append(cone_dic)
        return dictionaries

    def _compute_piecewise(self):
        max_degree = self._A.ncols()
        needed_points = factorial(self._num_variables + max_degree)//(
            factorial(self._num_variables)*factorial(max_degree) )

        for idx, cone_dict in enumerate(self._cone_dicts):
            cone_points = _generate_cone_points(cone_dict["scaled_rays"], needed_points)

            cone = cone_dict["cone"]
            ray_sum = free_module_element(sum(cone_dict["scaled_rays"]))

            polynomials = {}
            for off_set in cone_dict["quotient"]:
                nudge = _nudge_off_set(off_set.lift(), cone, ray_sum)
                points = [p + nudge for p in cone_points]
                # print(points)
                poly = self._interpolate_off_set_poly(points, max_degree)
                polynomials[off_set] = poly

            self._cone_dicts[idx]["polynomials"] = polynomials

    def _interpolate_off_set_poly(self, points, degree):
        num_integral_points = []
        for point in points:
            polytope = self._create_polytope_from_matrix(self._projection_inverse*point)
            num_integral_points.append(len(polytope.integral_points()))
        # print(num_integral_points)
        off_set_poly = self._R.interpolation(degree, points, num_integral_points)

        return off_set_poly

    def _create_polytope_from_matrix(self, b):
        """
        Unsafe version of ``create_polytope_from_matrix``
        to increase performance during lengthy computations.
        """
        inequalities = [[b[k]] + list(-self._A.rows()[k]) for k in range(self._A.nrows())]
        return Polyhedron(ieqs = inequalities)

    def __call__(self, point):
        return self.evaluate(point)

    def evaluate(self, point):
        if len(point) != self._amb_dim:
            raise ValueError("Dimension of ``point`` needs to be equal to the ambient"
                             f" dimension of ``self`` which is {self._amb_dim}.")

        proj = self._projection_matrix*free_module_element(point)
        for cone_dict in self._cone_dicts:
            if proj in cone_dict["cone"]:
                off_set = cone_dict["quotient"](proj)
                return cone_dict["polynomials"][off_set](*proj)
        return 0

    def _investigate(self, point):
        if len(point) != self._amb_dim:
            raise ValueError("Dimension of ``point`` needs to be equal to the ambient"
                             f" dimension of ``self`` which is {self._amb_dim}.")

        proj = self._projection_matrix*free_module_element(point)
        for cone_dict in self._cone_dicts:
            if proj in cone_dict["cone"]:
                off_set = cone_dict["quotient"](proj)
                return (proj, off_set, cone_dict["polynomials"][off_set](*proj))
        return 0

    def __repr__(self):
        return f"PiecewiseEhrhartQuasiPolynomial({self._A})"


def create_polytope_from_matrix(A, b):
    """
    Return the polytope whose faces are defined by
        Ax <= b
    where the inequality is understood componentwise

    EXAMPLES::

        sage: from ehrhart_quasi_polynomial.ehrhart_piecewise import create_polytope_from_matrix
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

def secondary_fan(A):
    """
    Return the secondary fan associated to ``A``
    For details see the documentation of the gfan software package:
        https://users-math.au.dk/~jensen/software/gfan/gfan.html

    EXAMPLES:

        sage: from ehrhart_quasi_polynomial.ehrhart_piecewise import secondary_fan
        sage: A = Matrix([[-1, 0], [0, -1], [1, 1]])
        sage: fan = secondary_fan(A); fan
        Polyhedral fan in 3 dimensions of dimension 3
        sage: fan.ambient_dim(), fan.dim(),
        (3, 3)
        sage: fan.rays()
        [[1, 1, 1]]
        sage: fan.cones()
        {3: [[0]]}
    """
    gfan_input = "{" + ", ".join(str(row) for row in A.rows()) + "}"
    return PolyhedralFan(gfan(gfan_input, "secondaryfan"))

def _compute_change_of_basis_matrices(orth_vectors, lin_vectors):
    """
    Return the change of basis matrix from the standard basis to the basis
    (``orth_vectors``, ``lineality_vectors``) along with its inverse and the projection
    into the space spanned by just ``orth_vectors`` and its "inverse".
    """
    M = create_matrix(list(orth_vectors) + list(lin_vectors))
    I = M.inverse()

    r = len(orth_vectors)
    m = len(lin_vectors)
    O = create_matrix.block([[create_matrix.identity(r),
                              create_matrix.zero(ZZ, nrows=r, ncols=m)]],
                            subdivide=False)
    K = O*I.transpose()*get_period(I) #not enough
    R = create_matrix(orth_vectors).transpose()*get_period(I)
    return M, I, K, R

def _compute_periods(A, points):
    """
    Return a list of the periods of `P_A(b)` for all `b` in ``points``
    An error will be raised if some ``b`` is not compatible with ``A``,
    see ``create_polytope_from_matrix`` for details.

    TESTS::
        
        sage: from ehrhart_quasi_polynomial.ehrhart_piecewise import _compute_periods
        sage: A = Matrix([[-1, 0], [0, -1], [2, 1]])
        sage: points = [(0, 0, 0), (0, 0, 1), (0, 0, 2)]
        sage: _compute_periods(A, points)
        [1, 2, 1]
    """
    return [get_period(create_polytope_from_matrix(A, b).Vrepresentation())
            for b in points]

def _compute_off_sets(matrix, projected_rays):
    standard_basis = create_matrix.identity(matrix.ncols()).rows()
    projected_basis = [matrix*vec for vec in standard_basis] + list(projected_rays)
    V = span(ZZ, projected_basis)
    W = span(ZZ, projected_rays)
    return V.quotient(W)

def _generate_cone_points(scaled_cone_rays, number):
    """
    Return a list of length ``number`` containing integer combinations of the
    ``scaled_cone_rays``.

    TESTS::

        sage: from ehrhart_quasi_polynomial.ehrhart_piecewise import _generate_cone_points
        sage: rays = [free_module_element([1, 0, 0]), free_module_element([0, 1, 1])]
        sage: _generate_cone_points(rays, 5)
        [(0, 0, 0), (1, 0, 0), (0, 1, 1), (2, 0, 0), (1, 1, 1)]
    """
    origin = scaled_cone_rays[0]*0
    points = [origin]
    points_len = 0
    index = 1
    while points_len <= number:
        new_points = [sum(combi) for combi
                      in combinations_with_replacement(scaled_cone_rays, index)]
        points += new_points
        index += 1
        points_len += len(new_points)
    return points[:number]

def _nudge_off_set(off_set, cone, ray_sum):
    nudge = off_set
    while nudge not in cone:
        nudge += ray_sum
    return nudge