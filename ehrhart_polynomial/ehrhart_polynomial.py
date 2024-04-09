from itertools import product
from collections import deque

from .integerperiodicfunction import IntegerPeriodicFunction
from .quasipolynomial import QuasiPolynomial

import sage.all
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.arith.misc import gcd
from sage.arith.functions import lcm
from sage.functions.other import ceil, floor

R = PolynomialRing(QQ, "x")
x = R.gen() 


# calculate ehrhart polynomial
def ehrhart_polynomial(vertices, simplify=False):
    y_values, scale_factor, period = points_contained_sequence(vertices, simplify)

    interpolation_points = [(k+1, y) for k, y in enumerate(y_values)]
    polynomial = interpolate_polynomial(interpolation_points, period, scale_factor)

    return polynomial

def interpolate_polynomial(points, period, scale_factor):
    if period == 1: # integral polytopes
        polynomial =  R.lagrange_polynomial(points)
        return polynomial(scale_factor*x)

    polynomials = [0]*period
    for k in range(period):
        period_points = points[k::period]
        polynomials[k] = R.lagrange_polynomial(period_points)

    return construct_quasipolynomial(polynomials, period, scale_factor)

def construct_quasipolynomial(polynomials, period, scale_factor):
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
        periodic_coefficients[degree] = IntegerPeriodicFunction(periodic_values)

    return QuasiPolynomial(periodic_coefficients)

# points contained
def points_contained_sequence(vertices, simplify):
    dimension = len(vertices[0])
    polytope_period = get_period(vertices)

    if simplify and polytope_period == 1:
        result = simplify_vertices(vertices, dimension)
        vertices, base_min, base_max, dimension, scale_factor = result
    else:
        base_min, base_max = get_bounding_extrema(vertices, dimension)
        scale_factor = 1

    base_poly = Polyhedron(vertices)

    if polytope_period == 1:
        iterations = dimension + 1
        bounding_box = get_bounding_box
    else:
        iterations = polytope_period*(dimension + 1)
        bounding_box = get_bounding_box_rational

    counting_sequence = [0]*(iterations)

    poly = Polyhedron(vertices)*0
    for k in range(iterations):
        poly += base_poly
        box = bounding_box(base_min, base_max, k+1)

        counting_sequence[k] = points_contained(poly, box)

    return counting_sequence, scale_factor, polytope_period

def points_contained(poly, box):
    contained = 0
    for point in box:
        if point in poly:
            contained += 1
    return contained


# period
def get_period(vertices):
    denominators = [QQ(coordinate).denominator() for vertex in vertices
                    for coordinate in vertex]
    period = lcm(denominators)
    return period


# bounding box
def get_bounding_extrema(vertices, dimension):
    columns = [[vertex[d] for vertex in vertices]
               for d in range(dimension)]
    mins = [min(col) for col in columns]
    maxs = [max(col) for col in columns]
    return mins, maxs

def get_bounding_box(mins, maxs, factor):
    return product(*[range(factor*mini, factor*maxi + 1)
                     for mini, maxi in zip(mins, maxs)])

def get_bounding_box_rational(mins, maxs, factor):
    return product(*[range(ceil(factor*mini), floor(factor*maxi) + 1)
                     for mini, maxi in zip(mins, maxs)])


# simplify polytope
def simplify_vertices(vertices, dimension):
    vertices, mins, maxs, new_dim = drop_constant_dimensions(vertices, dimension)
    new_vertices, scale_factor = scale_down_vertices(vertices)

    new_mins = [mini//scale_factor for mini in mins]
    new_maxs = [maxi//scale_factor for maxi in maxs]
    return new_vertices, new_mins, new_maxs, new_dim, scale_factor

def drop_dimensions(to_reduce, keep_filter):
    return [[val for val, keep in zip(vertex, keep_filter) if keep]
            for vertex in to_reduce]

def drop_constant_dimensions(vertices, dimension):
    columns = [[vertex[d] for vertex in vertices]
               for d in range(dimension)]
    mins = [min(col) for col in columns]
    maxs = [max(col) for col in columns]

    not_equal = [mins[d] != maxs[d] for d in range(dimension)]

    vertices = drop_dimensions(vertices, not_equal)
    mins, maxs = drop_dimensions([mins, maxs], not_equal)
    new_dimension = sum(not_equal)
    return vertices, mins, maxs, new_dimension

def scale_down_vertices(vertices):
    scale_factor = gcd(num for vertex in vertices for num in vertex)
    vertices = [[num//scale_factor for num in vertex] for vertex in vertices]
    return vertices, scale_factor

if __name__ == "__main__":
    rat = [(0, 0), (3/2, 0), (0, 1/3)]
    print(ehrhart_polynomial(rat))