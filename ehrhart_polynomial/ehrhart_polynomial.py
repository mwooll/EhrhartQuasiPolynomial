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

R = PolynomialRing(QQ, "x")
x = R.gen()

QPR = QuasiPolynomialRing(QQ)

# calculate ehrhart polynomial
def ehrhart_polynomial(vertices, simplify=True):
    y_values, scale_factor, period = _points_contained_sequence(vertices, simplify)

    interpolation_points = [(k+1, y) for k, y in enumerate(y_values)]
    polynomial = _interpolate_polynomial(interpolation_points, period, scale_factor)

    return polynomial


# interpolate polynomial
def _interpolate_polynomial(points, period, scale_factor):
    if period == 1: # integral polytopes
        polynomial = R.lagrange_polynomial(points)
        polynomial = polynomial(scale_factor*x)
        coefs = [float(c) for c in polynomial.coefficients(sparse=False)]
        return QPR(coefs)

    polynomials = [0]*period # rational polytopes
    for k in range(period):
        period_points = points[k::period]
        polynomials[k] = R.lagrange_polynomial(period_points)

    return _construct_quasipolynomial(polynomials, period, scale_factor)

def _construct_quasipolynomial(polynomials, period, scale_factor):
    polynomials = deque(polynomials)
    polynomials.rotate(1)
    polynomials = [poly(scale_factor*x) for poly in polynomials]

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

    return counting_sequence, scale_factor, polytope_period

def _points_contained(poly, box):
    contained = 0
    for point in box:
        if point in poly:
            contained += 1
    return contained


# period
def _get_period(vertices):
    denominators = [QQ(coordinate).denominator() for vertex in vertices
                    for coordinate in vertex]
    period = lcm(denominators)
    return period


# bounding box
def _get_bounding_extrema(vertices, dimension):
    r"""
    Return the bounding box which encompasses 'vertices'
    """
    columns = [[vertex[d] for vertex in vertices]
               for d in range(dimension)]
    mins = [min(col) for col in columns]
    maxs = [max(col) for col in columns]
    return mins, maxs

def _get_bounding_box(mins, maxs, factor):
    return product(*[range(factor*mini, factor*maxi + 1)
                     for mini, maxi in zip(mins, maxs)])

def _get_bounding_box_rational(mins, maxs, factor):
    return product(*[range(ceil(factor*mini), floor(factor*maxi) + 1)
                     for mini, maxi in zip(mins, maxs)])


# simplify polytope
def _simplify_vertices(vertices, dimension):
    r"""
    Simplifies the vertices as much as possible.
    Drop constant dimensions and scales down the remaining ones if possible.
    
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
