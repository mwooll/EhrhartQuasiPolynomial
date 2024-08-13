from itertools import combinations_with_replacement

from .ehrhart_quasi_polynomial import ehrhart_quasi_polynomial

import sage.all
from sage.functions.other import factorial
from sage.geometry.cone import Cone
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.interfaces.gfan import gfan
from sage.modules.free_module_element import free_module_element
from sage.rings.polynomial.groebner_fan import PolyhedralFan
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.structure.element import Matrix


def piecewise_ehrhart_quasi_polynomial(A, alter=None):
    """
    Computes the piecewise Ehrhart Quasi-Polynomial associate with ``A``.
    Return a piecewise quasi-polynomial where the domains are the maximal cones
    of the secondary fan of A and the functions on the domains are the 
    Ehrhart Quasi-Polynomials of A on those maximal cones.

    EXAMPLES::
        sage: from ehrhart_quasi_polynomial.ehrhart_piecewise import piecewise_ehrhart_quasi_polynomial
        sage: A = Matrix([[-1, 0], [0, -1], [1, 1], [0, 1]]) 
        sage: # eqp = piecewise_ehrhart_quasi_polynomial(A); epq

        sage: # epq.domains()
    """
    # check if inputs are compatible
    num_rows = A.nrows()
    if isinstance(alter, dict):
        for key, val in alter.items():
            if key >= num_rows:
                raise ValueError("'alter' cannot have keys which are larger or equal " +
                                     f"than the number of rows of `A`, which is {num_rows}")
            if val not in QQ:
                raise ValueError(f"values of 'alter' need to be rationals, got {val} of type {type(val)}")

    sec_fan = secondary_fan(A)
    estimated_period = _estimate_period(A, sec_fan)

    max_cones, piecewies_qps = _compute_piecewise(A, sec_fan, estimated_period, alter)
    return max_cones, piecewies_qps

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
    return _create_polytope_from_matrix(A, b) 


def secondary_fan(A):
    """
    Returns the secondary fan of the matrix ``A``.
    For computational details see the documentation of the gfan software package:
        https://users-math.au.dk/~jensen/software/gfan/gfan.html
    """
    gfan_input = "{" + ", ".join(str(row) for row in A.rows()) + "}"
    return PolyhedralFan(gfan(gfan_input, "secondaryfan"))

def _estimate_period(A, sec_fan): return 1

def _compute_piecewise(A_matrix, sec_fan, period=1, alter=None):
    num_variables = A_matrix.nrows()
    max_degree = sec_fan.lineality_dim()

    R = PolynomialRing(QQ, "x", num_variables)
    x_vars = R.gens()
    S = PolynomialRing(R, "k")
    k = S.gen()

    fan_rays = sec_fan.rays()

    if isinstance(alter, dict):
        fan_rays = _alter_rays(fan_rays, alter, num_variables)

    max_cones = []
    quasi_polynomials = []
    for ray_lists in sec_fan.maximal_cones().values():
        for ray_list in ray_lists:
            cone = Cone([fan_rays[idx] for idx in ray_list])
            max_cones.append(cone)

            red_cone, dimensions_removed = _remove_zero_dimensions(cone.rays())
            actual_num_variables = num_variables - len(dimensions_removed)
    
            needed_points = factorial(actual_num_variables + max_degree)//(
                factorial(actual_num_variables)*factorial(max_degree) )
            cone_points = _generate_cone_points(cone, needed_points)
            term_dict = _get_term_dict(A_matrix, cone_points)

            interpolation_variables = [x for dim, x in enumerate(x_vars)
                                       if dim not in dimensions_removed]
            T = PolynomialRing(QQ, names=interpolation_variables)
            cone_points = [[p for dim, p in enumerate(point)
                            if dim not in dimensions_removed]
                           for point in cone_points]

            cone_poly = 0
            for degree, terms in term_dict.items():
                cone_poly += T.interpolation(max_degree, cone_points, terms)*k**degree

            quasi_polynomials.append(cone_poly)
    return max_cones, quasi_polynomials

def _alter_rays(rays, alter):
    altered_rays = []
    for ray in rays:
        altered_ray = ray
        for key, val in alter.items():
            altered_ray[key] = val
        altered_rays.append(altered_ray)
    return altered_rays

def _remove_zero_dimensions(rays):
    zero_dimensions = []
    for dim, val in enumerate(rays[0]):
        if val == 0:
            if all(ray[dim] == 0 for ray in rays[1:]):
                zero_dimensions.append(dim)
    reduced_cone = Cone([[val for dim, val in enumerate(ray) if dim not in zero_dimensions ]for ray in rays])
    return reduced_cone, zero_dimensions

def _generate_cone_points(cone, number):
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

def _get_term_dict(A_matrix, points):
    polynomials = []
    for point in points:
        polytope = _create_polytope_from_matrix(A_matrix, point)
        ehr_poly = ehrhart_quasi_polynomial(polytope.Vrepresentation())
        polynomials.append(ehr_poly)

    max_degree = max(poly.degree() for poly in polynomials)
    terms = {d: _get_terms_of_order(polynomials, d)for d in range(max_degree + 1)}
    return terms

def _create_polytope_from_matrix(A, b):
    inequalities = [[b[k]] + list(-A.rows()[k]) for k in range(A.nrows())]
    return Polyhedron(ieqs = inequalities)

def _get_terms_of_order(polynomials, order):
    return [poly.coefficients()[order].constants()[0]
            if poly.degree() >= order else 0 for poly in polynomials]
