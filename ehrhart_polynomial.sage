import numpy as np
from math import factorial
from scipy.interpolate import lagrange

from boundingbox import BoundingBox, get_bounding_box

from unittest import TestCase, main

x = var("x")

def ehrhart_polynomial(vertices):
    y_values = points_contained(vertices)

    # length = len(y_values)
    # x_values = [x for x in range(1, len(y_values))]
    polynomial = interpolate_polynomial(y_values)

    return polynomial

def points_contained(vertices):
    base_poly = Polyhedron(vertices)
    base_box = get_bounding_box(vertices)
    dimension = base_box.dim

    points_contained = [0]*(dimension+1)
    for k in range(1, dimension+2):
        poly = k*base_poly
        box = k*base_box

        contained = 0
        for point in box:
            if point in poly:
                contained += 1 

        points_contained[k-1] = contained

    return points_contained

# interpolate polynomial from sequence
def interpolate_polynomial(sequence):
    degree = 0
    differences = sequence[:]

    while len(set(differences)) > 1:
        differences = [differences[k] - differences[k-1] for k in range(1, len(differences))]
        degree += 1

    if degree == 0:
        return sequence[0]

    factor = differences[0]/factorial(degree)

    new_sequence = [sequence[k] - factor*(k+1)^degree for k in range(len(sequence))]

    return factor*x^degree + interpolate_polynomial(new_sequence)
    
# mkae it interpolate faster
def interpolate_polynomial_faster(sequence):
    degree = 0
    differences = sequence[:]

    while len(set(differences)) > 1:
    # while differences[0] != 0: # makes a RunTimeError occur... no idea why
        differences = [differences[k] - differences[k-1] for k in range(1, len(differences))]
        degree += 1

    if degree == 0:
        return sequence[0]

    factor = differences[0]/factorial(degree)

    new_sequence = [sequence[k] - factor*(k+1)^degree for k in range(degree)]

    return factor*x^degree + interpolate_polynomial(new_sequence)

interpolate_polynomial = interpolate_polynomial_faster


class TestEhrhartPolynomial(TestCase):
    def test_ehrhart_polynomial(self):
        point_poly = ehrhart_polynomial([[0, 0, 0, 0, 0]])
        self.assertEqual(point_poly, 1)

        axis = ehrhart_polynomial([[-3], [3]])
        self.assertEqual(axis, 6*x + 1)
    
        triangle = ehrhart_polynomial([[0, 0], [1, 0], [0, 1]])
        self.assertEqual(triangle, x^2/2 + 1.5*x + 1)
    
        flat_tri = ehrhart_polynomial([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        self.assertEqual(flat_tri, x^2/2 + 1.5*x + 1)

        triangle_pyramid = ehrhart_polynomial([[0, 0, 0], [1, 0, 0],
                                                [0, 1, 0], [0, 0, 1]])
        self.assertEqual(triangle_pyramid, x^3/6 + x^2 + 11/6*x + 1)
    
        square_poly = ehrhart_polynomial([[0, 0], [1, 0], [1, 1], [0, 1]])
        self.assertEqual(square_poly, x^2 + 2*x + 1) # (x + 1)^2
    
        vertices = BoundingBox(np.zeros(3, dtype=int), np.ones(3, dtype=int)).vertices
        cube_poly = ehrhart_polynomial(vertices)
        self.assertEqual(cube_poly, x^3 + 3*x^2 + 3*x + 1) # (x + 1)^3
    

    def test_points_contained(self):
        point = points_contained([[0, 0, 0, 0, 0]])
        self.assertEqual(point, [1, 1, 1, 1, 1, 1])

        square = points_contained([[0, 0], [1, 0], [1, 1], [0, 1]])
        self.assertEqual(square, [4, 9, 16])

        vertices = BoundingBox(np.zeros(3, dtype=int), np.ones(3, dtype=int)).vertices
        cube = points_contained(vertices)
        self.assertEqual(cube, [8, 27, 64, 125])

        axis = points_contained([[-3], [3]])
        self.assertEqual(axis, [7, 13])

        triangle = points_contained([[0, 0], [1, 0], [0, 1]])
        self.assertEqual(triangle, [3, 6, 10])

        flat_tri = points_contained([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        self.assertEqual(flat_tri, [3, 6, 10, 15])

    def test_predict_poly_rec(self):
        actual = interpolate_polynomial( (1, 1, 1, 1) )
        self.assertEqual(1, actual)
        
        actual = interpolate_polynomial( (1, 2, 3, 4) )
        self.assertEqual(x, actual)

        points = [x^3 + 2*x - 5 for x in range(1, 5)]
        actual = interpolate_polynomial(points)
        self.assertEqual(x^3 + 2*x - 5, actual)

        points = [x^17 - x^16 + x^15 + x^3 - 2*x^2 + 102*x for x in range(1, 19)]
        actual = interpolate_polynomial(points)
        self.assertEqual(x^17 - x^16 + x^15 + x^3 - 2*x^2 + 102*x, actual)

if __name__ == "__main__":
    main()
    