from itertools import combinations_with_replacement

from .ehrhart_quasi_polynomial import get_period, get_gcd

from sage.functions.other import factorial
from sage.geometry.cone import Cone
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.interfaces.gfan import gfan
from sage.matrix.constructor import Matrix as create_matrix
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

        self._amb_dim = sec_fan.ambient_dim()
        self._ZZ_lattice = IntegralLattice(create_matrix.identity(self._amb_dim))

        self._rays = tuple(free_module_element(ray) for ray in sec_fan.rays())
        self._scalars = _hat_denominator(self._A, self._rays)

        self._maximal_cones = sec_fan.maximal_cones()

        self._lin_vectors = _process_fan_vectors(sec_fan.fan_dict["LINEALITY_SPACE"])
        self._minus_lin = tuple(-lin for lin in self._lin_vectors)

        self._num_variables = self._amb_dim - len(self._lin_vectors)
        self._R = PolynomialRing(QQ, "x", self._num_variables)
        self._S = PolynomialRing(QQ, "y", self._amb_dim)


    def _generate_cone_dicts(self):
        dictionaries = []
        for ray_lists in  self._maximal_cones.values():
            for ray_list in ray_lists:
                cone_dic = {}
                cone_dic["scaled_rays"] = tuple(self._scalars[idx]*self._rays[idx]
                                                for idx in ray_list)
                cone_basis = cone_dic["scaled_rays"] + self._lin_vectors
                cone_dic["cone"] = Cone(cone_basis + self._minus_lin)

                cone_lattice = self._ZZ_lattice.sublattice(cone_basis)
                cone_dic["quotient"] = self._ZZ_lattice.quotient(cone_lattice)
                cone_dic["basis"] = _get_basis_vectors(self._lin_vectors,
                                                       cone_dic["scaled_rays"],
                                                       self._amb_dim)

                K, M, = _compute_change_of_basis_matrices(cone_dic["basis"],
                                                          self._lin_vectors)
                cone_dic["change_of_basis_matrix"] = K
                cone_dic["change_of_basis_inverse"] = M
                cone_dic["lift_matrix"] = create_matrix(cone_dic["basis"]).T

                dictionaries.append(cone_dic)
        return dictionaries

    def _compute_piecewise(self):
        max_degree = self._A.ncols()
        needed_points = factorial(self._num_variables + max_degree + 1)//(
            factorial(self._num_variables)*factorial(max_degree) )

        points = _generate_cone_points(self._num_variables, needed_points)
        for cone_dict in self._cone_dicts:
            cone_points = {0: (points, [tuple(b) for b in points])}
            ray_sum = sum(cone_dict["scaled_rays"])

            polynomials = {}
            for off_set in cone_dict["quotient"]:
                nudge = off_set.lift()
                mult = 0
                while nudge not in cone_dict["cone"]:
                    mult += 1
                    nudge += ray_sum

                if mult not in cone_points:
                    mults = free_module_element([mult for k in range(self._num_variables)])
                    higher_points = [b + mults for b in cone_points[0][0]]
                    cone_points[mult] = (higher_points, [tuple(b) for b in higher_points])

                num_integral_points = []
                for point in cone_points[mult][0]:
                    polytope_point = cone_dict["lift_matrix"]*point + off_set.lift()
                    polytope = self._create_polytope_from_matrix(polytope_point)
                    num_integral_points.append(len(polytope.integral_points()))

                # print(max_degree, cone_points[mult][1], num_integral_points)
                off_set_poly = self._R.interpolation(max_degree,
                                                     cone_points[mult][1],
                                                     num_integral_points)
                polynomials[off_set] = off_set_poly

            cone_dict["polynomials"] = polynomials


    def _create_polytope_from_matrix(self, b):
        """
        Unsafe version of ``create_polytope_from_matrix``
        to increase performance during lengthy computations.
        """
        inequalities = [[b[k]] + list(-self._A.rows()[k]) for k in range(self._A.nrows())]
        return Polyhedron(ieqs = inequalities)

    # def _transform_polynomial(self, cone_dict, cone_polynomial):
    #     coefficients = [cone_dict["change_of_basis_inverse"]*vector
    #                     for vector in cone_dict["basis"]]
    #     new_variables = [free_module_element(self._S.gens())*coeff
    #                      for coeff in coefficients]
    #     transformed_polynomial = cone_polynomial(*new_variables)
    #     return transformed_polynomial

    def __call__(self, point):
        return self.evaluate(point)

    def evaluate(self, point):
        if len(point) != self._amb_dim:
            raise ValueError("Dimension of ``point`` needs to be equal to the ambient"
                              f" dimension of ``self`` which is {self._amb_dim}.")
    
        for cone_dict in self._cone_dicts:
            if point in cone_dict["cone"]:
                # print(True)
                off_set = cone_dict["quotient"][cone_dict["quotient"](point)[0]]
                eval_p = cone_dict["change_of_basis_inverse"]*(
                    free_module_element(point) - off_set.lift())
                return cone_dict["polynomials"][off_set](*eval_p[:self._num_variables])
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
        sage: fan.cones()
        {3: [[0]]}
    """
    gfan_input = "{" + ", ".join(str(row) for row in A.rows()) + "}"
    return PolyhedralFan(gfan(gfan_input, "secondaryfan"))

def _process_fan_vectors(vectors):
    return tuple(free_module_element([int(val) for val in vec.split(" ")])
                 for vec in vectors)

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
    return tuple(get_period(create_polytope_from_matrix(A, b).Vrepresentation())
                 for b in points)

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
        [1, 2, 1]

        sage: points = [()]
        sage: _hat_denominator(A, points)
        []
    """
    polytopes = [create_polytope_from_matrix(A, b).Vrepresentation()
                 for b in points]
    return tuple(get_period(poly) / get_gcd(poly)
                 for poly in polytopes)

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

def _generate_cone_points(dimension, number):
    """
    Return a tuple of length ``number`` containing integral points of
    dimension ``dimension``.

    TESTS::

        sage: from ehrhart_quasi_polynomial.ehrhart_piecewise import _generate_cone_points
        sage: _generate_cone_points(3, 5)
        [(0, 0, 0), (1, 0, 0), (0, 1, 1), (2, 0, 0), (1, 1, 1)]
    """
    origin = free_module_element(0 for k in range(dimension))
    points = (origin, )
    points_len = 0

    basis = [free_module_element(int(k == d) for k in range(dimension))
             for d in range(dimension)]

    index = 1
    while points_len <= number:
        new_points = tuple(sum(combi) for combi in
                           combinations_with_replacement(basis, index))
        points += new_points
        index += 1
        points_len += len(new_points)
    return points[:number]
