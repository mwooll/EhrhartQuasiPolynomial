from itertools import product
from collections import deque

from .quasipolynomial import QuasiPolynomialRing

import sage.all
from sage.arith.misc import gcd
from sage.arith.functions import lcm
from sage.functions.other import ceil, floor
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.geometry.polyhedron.constructor import Polyhedron

R = PolynomialRing(QQ, "x", sparse=False)
x = R.gen()

QPR = QuasiPolynomialRing(QQ)

# calculate ehrhart polynomial
def ehrhart_polynomial(vertices, simplify=True):
    y_values, period, scale_factor = _points_contained_sequence(vertices, simplify)

    interpolation_points = [(k+1, y) for k, y in enumerate(y_values)]
    polynomial = _interpolate_polynomial(interpolation_points, period, scale_factor)

    return polynomial


# interpolate polynomial
def _interpolate_polynomial(points_sequence, period, scale_factor):
    r"""
    Computes the Quasi-Polynomial which evaluates to 'points_sequence' and has
    the specified period.
    Return a 'ehrhart_polynomial.quasipolynomial.QuasiPolynomialElement'

    TESTS:

        sage: from ehrhart_polynomial.ehrhart_polynomial import (
        ....:    _points_contained_sequence, _interpolate_polynomial)
        sage: vertices = [[0, 0], [2, 0], [2, 2], [0, 2]]
        sage: sequence, period, factor = _points_contained_sequence(vertices, True)
        sage: points_sequence = [(k+1, y) for k, y in enumerate(sequence)]
        sage: _interpolate_polynomial(points_sequence, period, factor)
        QuasiPolynomialElement(Ring of Quasi-Polynomials over Rational Field, [[1], [4], [4]])
        sage: vertices = [[0, 0], [1/2, 1/3]]
        sage: sequence, period, factor = _points_contained_sequence(vertices, False)
        sage: points_sequence = [(k+1, y) for k, y in enumerate(sequence)]
        sage: _interpolate_polynomial(points_sequence, period, factor)
        QuasiPolynomialElement(Ring of Quasi-Polynomials over Rational Field, [[1, 5/6, 2/3, 1/2, 1/3, 1/6], [1/6]])
    """
    if period == 1: # integral polytopes
        polynomial = R.lagrange_polynomial(points_sequence)
        polynomial = polynomial(scale_factor*x)
        coefs = [float(c) for c in polynomial.coefficients(sparse=False)]
        return QPR(coefs)

    polynomials = [0]*period # rational polytopes
    for k in range(period):
        period_points = points_sequence[k::period]
        polynomials[k] = R.lagrange_polynomial(period_points)

    return _construct_quasipolynomial(polynomials, period)

def _construct_quasipolynomial(polynomials, period):
    r"""
    Construct a Quasi-Polynomial out of 'polynomials' with the specified period.
    Return a 'ehrhart_polynomial.quasipolynomial.QuasiPolynomialElement'

    TESTS::

        sage: from ehrhart_polynomial.ehrhart_polynomial import _construct_quasipolynomial
        sage: x = PolynomialRing(QQ, "x").gen()
        sage: polys = [x+1, 2*x+2, x+3, 2*x+4]
        sage: _construct_quasipolynomial(polys, 4)
        QuasiPolynomialElement(Ring of Quasi-Polynomials over Rational Field, [[4, 1, 2, 3], [2, 1]])
    """
    polynomials = deque(polynomials)
    polynomials.rotate(1)

    degrees = [poly.degree() for poly in polynomials]
    max_degree = max(degrees)
    periodic_coefficients = [0]*(max_degree + 1)

    for degree in range(max_degree + 1):
        periodic_values = [0]*period
        for index, poly in enumerate(polynomials):
            if degree <= degrees[index]:
                periodic_values[index] = poly.coefficients()[degree]
        periodic_coefficients[degree] = periodic_values

    return QPR(periodic_coefficients)

# points contained
def _points_contained_sequence(vertices, simplify):
    r"""
    Compute sequence [(k, a_k)] where a_k is the number of integral points of 'vertices'*k,
    and k goes from 0 to (period of vertices)*(dimension of vertices), if 'simplify'
    is True, then k goes to (period of vertices)*(dimension of simplified vertices)
    Return the computed sequence, 'scale_factor' from '_scale_down_vertices()'
    and the period of 'vertices'

    Note: if 'vertices' are truly rational, then 'simplify' has no effect.

    TESTS::

        sage: from ehrhart_polynomial.ehrhart_polynomial import _points_contained_sequence
        sage: _points_contained_sequence([[0, 0], [2, 0], [2, 2], [0, 2]], False)
        ([9, 25, 49], 1, 1)
        sage: _points_contained_sequence([[0, 0], [2, 0], [2, 2], [0, 2]], True)
        ([4, 9, 16], 1, 2)
        sage: _points_contained_sequence([[0, 0], [1/2, 0], [0, 1/2]], False)
        ([1, 3, 3, 6, 6, 10], 2, 1)
        sage: _points_contained_sequence([[0, 0], [1/2, 0], [0, 1/2]], True)
        ([1, 3, 3, 6, 6, 10], 2, 1)
    """
    dimension = len(vertices[0])
    polytope_period = _get_period(vertices)

    if simplify and polytope_period == 1:
        result = _simplify_vertices(vertices, dimension)
        vertices, base_min, base_max, dimension, scale_factor = result
    else:
        base_min, base_max = _get_bounding_extrema(vertices, dimension)
        scale_factor = 1

    base_poly = Polyhedron(vertices)

    if polytope_period == 1:
        iterations = dimension + 1
        bounding_box = _get_bounding_box
    else:
        iterations = polytope_period*(dimension + 1)
        bounding_box = _get_bounding_box_rational

    counting_sequence = [0]*(iterations)

    poly = Polyhedron(vertices)*0
    for k in range(iterations):
        poly += base_poly
        box = bounding_box(base_min, base_max, k+1)

        counting_sequence[k] = _points_contained(poly, box)

    return counting_sequence, polytope_period, scale_factor

def _points_contained(poly, box):
    r"""
    Return the number of points in 'box' which are inside 'poly'.

    TESTS::

        sage: from ehrhart_polynomial.ehrhart_polynomial import (
        ....:       _points_contained, _get_bounding_extrema,
        ....:       _get_bounding_box, _get_bounding_box_rational)
        sage: int_vertices = [[0, 0], [0, 2], [2, 0]]
        sage: int_box = _get_bounding_box(*_get_bounding_extrema(int_vertices, 2), 1)
        sage: _points_contained(Polyhedron(int_vertices), int_box)
        6
        sage: rat_vertices = [[0, 0], [0, 7/2], [5/4, 0]]
        sage: rat_box = _get_bounding_box_rational(*_get_bounding_extrema(rat_vertices, 2), 2)
        sage: _points_contained(Polyhedron(rat_vertices)*2, rat_box)
        15
    """
    return sum(int(point in poly) for point in box)

# period
def _get_period(vertices):
    r"""
    Return the period of 'vertices'. The period of a set of vertices
    is defined as the lcm of the denominators of all vertices.

    TESTS::

        sage: from ehrhart_polynomial.ehrhart_polynomial import _get_period
        sage: _get_period([[0, 0], [1, 0], [1, 1], [0, 1]])
        1
        sage: _get_period([[-1/2, 0], [2/3, 0], [0, 1/5]])
        30
        sage: _get_period([[1, 1/2], [2, 1/4], [3, 1/8]])
        8
    """
    denominators = [QQ(coordinate).denominator() for vertex in vertices
                    for coordinate in vertex]
    period = lcm(denominators)
    return period

# bounding box
def _get_bounding_extrema(vertices, dimension):
    r"""
    Return the bounding box which encompasses 'vertices'.

    TESTS::

        sage: from ehrhart_polynomial.ehrhart_polynomial import _get_bounding_extrema
        sage: _get_bounding_extrema([[0, 0, 0], [-1, 0, 1]], 3)
        ([-1, 0, 0], [0, 0, 1])
        sage: _get_bounding_extrema([[-1, 1], [1, -1]], 2)
        ([-1, -1], [1, 1])
    """
    columns = [[vertex[d] for vertex in vertices]
               for d in range(dimension)]
    mins = [min(col) for col in columns]
    maxs = [max(col) for col in columns]
    return mins, maxs

def _get_bounding_box(mins, maxs, factor):
    r"""
    Return an iterator of the integral points between 'mins'*'factor' and
    'maxs'*'factor'.

    TESTS::

        sage: from ehrhart_polynomial.ehrhart_polynomial import (_get_bounding_box,
        ....:                                                    _get_bounding_extrema)
        sage: mins, maxs = _get_bounding_extrema([[0, 0], [1, 1]], 2)
        sage: list(_get_bounding_box(mins, maxs, 2)) # doctest: +NORMALIZE_WHITESPACE
        [(0, 0), (0, 1), (0, 2),
         (1, 0), (1, 1), (1, 2),
         (2, 0), (2, 1), (2, 2)]
    """
    return product(*[range(factor*mini, factor*maxi + 1)
                     for mini, maxi in zip(mins, maxs)])

def _get_bounding_box_rational(mins, maxs, factor):
    r"""
    Return an iterator of the integral points between ceil('mins'*'factor') and
    floor('maxs'*'factor').

    TESTS::

        sage: from ehrhart_polynomial.ehrhart_polynomial import (_get_bounding_box_rational,
        ....:                                                    _get_bounding_extrema)
        sage: mins, maxs = _get_bounding_extrema([[0, 0], [1/2, 3/2]], 2)
        sage: list(_get_bounding_box_rational(mins, maxs, 2))
        [(0, 0), (0, 1), (0, 2), (0, 3),
         (1, 0), (1, 1), (1, 2), (1, 3)]
    """
    return product(*[range(ceil(factor*mini), floor(factor*maxi) + 1)
                     for mini, maxi in zip(mins, maxs)])

# simplify polytope
def _simplify_vertices(vertices, dimension):
    r"""
    Simplify the vertices as much as possible:
        - drop constant dimensions
        - scale down remaining dimensions
    
    Return simplified vertices, minimums and maximums over all dimensions, the new
    dimension of the vertices and the factor by which the vertices were scaled down.

    Note: this function and the ones this function calls are somewhat designed for
    speed and therefore do not test any input for validity. Hence, it is crucial
    that all elements of 'vertices' have the same dimension and to get the wanted
    bahavior 'dimension' should really be the dimension of the vertices.

    TESTS::

        sage: from ehrhart_polynomial.ehrhart_polynomial import _simplify_vertices
        sage: _simplify_vertices([[0, 0, 1], [0, 0, 2]], 3)
        ([[1], [2]], [1], [2], 1, 1)
        sage: _simplify_vertices([[5, 25], [15, 10]], 2)
        ([[1, 5], [3, 2]], [1, 2], [3, 5], 2, 5)
        sage: _simplify_vertices([[7, 5], [7, 0], [7, 15]], 2)
        ([[1], [0], [3]], [0], [3], 1, 5)
    """
    vertices, mins, maxs, new_dim = _drop_constant_dimensions(vertices, dimension)
    new_vertices, scale_factor = _scale_down_vertices(vertices)

    new_mins = [mini//scale_factor for mini in mins]
    new_maxs = [maxi//scale_factor for maxi in maxs]
    return new_vertices, new_mins, new_maxs, new_dim, scale_factor

def _drop_constant_dimensions(vertices, dimension):
    r"""
    Drop all dimensions <= dimension of 'vertices' where every vertex has the same value.

    Return cleaned vertices, minimum and maximum over all vertices and dimension,
    and the new dimension of the vertices.

    TESTS::

        sage: from ehrhart_polynomial.ehrhart_polynomial import _drop_constant_dimensions
        sage: _drop_constant_dimensions([[0, 1, 2, 5], [0, 1, 4, 5]], 4)
        ([[2], [4]], [2], [4], 1)
        sage: _drop_constant_dimensions([[0, -1, 0], [0, 1, 3], [0, 2, 2]], 3)
        ([[-1, 0], [1, 3], [2, 2]], [-1, 0], [2, 3], 2)

    """
    columns = [[vertex[d] for vertex in vertices]
               for d in range(dimension)]
    mins = [min(col) for col in columns]
    maxs = [max(col) for col in columns]

    not_equal = [mins[d] != maxs[d] for d in range(dimension)]

    vertices = _drop_dimensions(vertices, not_equal)
    mins, maxs = _drop_dimensions([mins, maxs], not_equal)
    new_dimension = sum(not_equal)
    return vertices, mins, maxs, new_dimension

def _drop_dimensions(to_reduce, keep_filter):
    r"""
    Return 'to_reduce' where the dimensions are filtered according to 'keep_filter'

    TESTS::

        sage: from ehrhart_polynomial.ehrhart_polynomial import _drop_dimensions
        sage: _drop_dimensions([[0, 1, 2], [3, 4, 5]], [True, False, True])
        [[0, 2], [3, 5]]
        sage: _drop_dimensions([[0, 1, 2], [3, 4, 5]], [1, 1])
        [[0, 1], [3, 4]]
    """
    return [[val for val, keep in zip(vertex, keep_filter) if keep]
            for vertex in to_reduce]

def _scale_down_vertices(vertices):
    r"""
    Computes 'scale_factor' as gcd of all coordinates of the vertices.
    Return vertices divided by scale_factor and 'scale_factor'.

    TESTS::

        sage: from ehrhart_polynomial.ehrhart_polynomial import _scale_down_vertices
        sage: _scale_down_vertices([[0, 0], [5, 10], [15, 20]])
        ([[0, 0], [1, 2], [3, 4]], 5)
        sage: _scale_down_vertices([[0, 1], [3, 5]])
        ([[0, 1], [3, 5]], 1)
    """
    scale_factor = gcd(num for vertex in vertices for num in vertex)
    vertices = [[num//scale_factor for num in vertex] for vertex in vertices]
    return vertices, scale_factor
