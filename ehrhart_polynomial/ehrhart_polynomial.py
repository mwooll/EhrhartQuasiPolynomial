from itertools import product

#import sage.all
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.geometry.polyhedron.constructor import Polyhedron

R = PolynomialRing(QQ, "x")
x = R.gen()


def ehrhart_polynomial(vertices):
    y_values = points_contained(vertices)

    interpolation_points = [(k+1, y) for k, y in enumerate(y_values)]
    polynomial = R.lagrange_polynomial(interpolation_points)

    return polynomial

def points_contained(vertices):
    dimension = len(vertices[0])

    base_poly = Polyhedron(vertices)
    base_min, base_max = get_bounding_extrema(vertices, dimension)

    points_contained = [0]*(dimension+1)
    for k in range(1, dimension+2):
        poly = k*base_poly
        box = get_bounding_box(base_min, base_max, k)

        contained = 0
        for point in box:
            if point in poly:
                contained += 1 

        points_contained[k-1] = contained

    return points_contained

def get_bounding_extrema(vertices, dimension):
    columns = [[vertex[d] for vertex in vertices]
               for d in range(dimension)]
    mins = [min(col) for col in columns]
    maxs = [max(col) for col in columns]
    return mins, maxs

def get_bounding_box(mins, maxs, factor):
    return product(*[range(factor*mini, factor*maxi + 1)
                     for mini, maxi in zip(mins, maxs)])

