from itertools import product, combinations_with_replacement

from .ehrhart_quasi_polynomial import get_period, get_gcd

from sage.geometry.cone import Cone
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.interfaces.gfan import gfan
from sage.matrix.constructor import Matrix as create_matrix
from sage.modules.free_module_element import free_module_element
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice
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

        self._amb_dim = sec_fan.ambient_dim()
        self._ZZ_lattice = IntegralLattice(create_matrix.identity(self._amb_dim))

        self._rays = tuple(free_module_element(ray) for ray in sec_fan.rays())
        self._scalars = _hat_denominator(self._A, self._rays)

        self._maximal_cones = sec_fan.maximal_cones()

        self._lin_vectors = _process_fan_vectors(sec_fan.fan_dict["LINEALITY_SPACE"])
        self._minus_lin = tuple(-lin for lin in self._lin_vectors)

        self._max_degree = self._A.ncols()
        self._num_variables = self._amb_dim - len(self._lin_vectors)
        self._R = PolynomialRing(QQ, "x", self._num_variables)


    def _generate_cone_dicts(self):
        cone_dicts = []
        for ray_lists in  self._maximal_cones.values():
            for ray_list in ray_lists:
                cone_dict = {}
                cone_dict["scaled_rays"] = tuple(self._scalars[idx]*self._rays[idx]
                                                for idx in ray_list)
                cone_basis = cone_dict["scaled_rays"] + self._lin_vectors
                cone_dict["cone"] = Cone(cone_basis + self._minus_lin)

                cone_lattice = self._ZZ_lattice.sublattice(cone_basis)
                cone_dict["quotient"] = self._ZZ_lattice.quotient(cone_lattice)
                cone_dict["lifts"] = _determine_lifts(cone_dict["quotient"])
                # print(cone_dict["lifts"])

                cone_dict["basis"] = _get_basis_vectors(self._lin_vectors,
                                                       cone_dict["scaled_rays"],
                                                       self._amb_dim)

                K, M = _compute_change_of_basis_matrices(cone_dict["basis"],
                                                         self._lin_vectors)
                cone_dict["change_of_basis_matrix"] = K
                cone_dict["change_of_basis_inverse"] = M
                cone_dict["lift_matrix"] = create_matrix(cone_dict["basis"]).T

                cone_dicts.append(cone_dict)
        return cone_dicts

    def _compute_piecewise(self):
        points = _generate_cone_points(self._num_variables, self._max_degree)
        cone_points = {0: points}
        for cone_dict in self._cone_dicts:
            ray_sum = sum(cone_dict["scaled_rays"])

            polynomials = {}
            for rep, lift in cone_dict["lifts"].items():
                mult = self._nudge_off_set(lift, cone_points,
                                           cone_dict["cone"], ray_sum)
                off_set_poly = self._get_off_set_poly(lift, cone_points,
                                                      cone_dict["lift_matrix"], mult)
                polynomials[rep] = off_set_poly

            cone_dict["polynomials"] = polynomials

    def _nudge_off_set(self, lift, cone_points, cone, ray_sum):
        mult = 0
        
        while lift not in cone:
            mult += 1
            lift += ray_sum

        if mult not in cone_points:
            mults = free_module_element(mult for k in range(self._num_variables))
            higher_points = [b + mults for b in cone_points[0]]
            cone_points[mult] = higher_points

        return mult

    def _get_off_set_poly(self, lift, cone_points, lift_matrix, mult):
        num_integral_points = []
        for point in cone_points[mult]:
            polytope_point = lift_matrix*point + lift
            polytope = self._create_polytope_from_matrix(polytope_point)
            num_integral_points.append(len(polytope.integral_points()))

        off_set_poly = self._R.interpolation(self._max_degree,
                                             cone_points[mult],
                                             num_integral_points)
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
    
        for k, cone_dict in enumerate(self._cone_dicts):
            if point in cone_dict["cone"]:
                rep = tuple(cone_dict["quotient"](point))
                lift = cone_dict["lifts"][rep]
                eval_p = ( cone_dict["change_of_basis_inverse"]*(
                    free_module_element(point) - lift) )[:self._num_variables]

                result = cone_dict["polynomials"][rep](*eval_p)
                return result
        return 0

    def __repr__(self):
        return f"PiecewiseEhrhartQuasiPolynomial({self._A})"


def create_polytope_from_matrix(A, b):
    """
    Return the polytope whose faces are defined by
        Ax <= b
    <=>
         0 <= b - Ax
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
    inequalities = tuple([b[k]] + list(-A.rows()[k]) for k in range(A.nrows()))
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
        sage: fan.maximal_cones()
        {3: [[0]]}
    """
    gfan_input = "{" + ", ".join(str(row) for row in A.rows()) + "}"
    return PolyhedralFan(gfan(gfan_input, "secondaryfan"))

def _process_fan_vectors(vectors):
    return tuple(free_module_element([int(val) for val in vec.split(" ")])
                 for vec in vectors)

def _hat_denominator(A, points):
    """
    Return a list of the den_hat of `P_A(b)` for all `b` in ``points``
    An error will be raised if some ``b`` is not compatible with ``A``,
    see ``create_polytope_from_matrix`` for details.

    TESTS::
        
        sage: from ehrhart_quasi_polynomial.ehrhart_piecewise import _hat_denominator
        sage: A = Matrix([[-1, 0], [0, -1], [2, 1]])
        sage: points = [(0, 0, 0), (0, 0, 1), (0, 0, 2)]
        sage: _hat_denominator(A, points)
        (1, 4, 1)
    """
    polytopes = [create_polytope_from_matrix(A, b).Vrepresentation()
                 for b in points]
    return tuple(get_period(poly) / get_gcd(poly)
                 for poly in polytopes)

def _determine_lifts(quotient):
    representatives = product(*(range(k) for k in quotient.invariants()))
    lift_matrix = create_matrix([gen.lift() for gen in quotient.gens()]).T
    lifts = {tuple(rep): lift_matrix*free_module_element(rep)
             for rep in representatives}
    return lifts

def _get_basis_vectors(lin_vectors, scaled_rays, amb_dim):
    if len(lin_vectors) + len(scaled_rays) == amb_dim:
        return scaled_rays

    raise NotImplementedError("Need to construct basis of R^n out of n+k vectors.")

def _compute_change_of_basis_matrices(cone_basis, lin_vectors):
    """
    Return the change of basis matrix from the standard basis to the basis
    (``lineality_vectors``, ``cone_basis``) along with its inverse.
    """
    K = create_matrix(cone_basis + lin_vectors).T
    M = K.inverse()
    return K, M

def _generate_cone_points(dimension, scaler):
    """
    Return all integral points of a simplex of dimension ``dim`` scaled
    by ``scaler``. Always returns exactly `dim + scaler choose dim` points.

    TESTS::

        sage: from ehrhart_quasi_polynomial.ehrhart_piecewise import _generate_cone_points
        sage: _generate_cone_points(2, 2)
        ((0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2))
    """
    origin = free_module_element(0 for k in range(dimension))
    points = (origin, )

    basis = [free_module_element(int(k == d) for k in range(dimension))
             for d in range(dimension)]

    points += tuple(sum(combi) for index in range(1, scaler+1)
                    for combi in combinations_with_replacement(basis, index))
    return points
