from ehrhart_quasi_polynomial.ehrhart_piecewise import PiecewiseEhrhartQuasiPolynomial

from sage.matrix.constructor import Matrix

A = Matrix([[-1, 0], [0, -1], [1, 3]])
peqp = PiecewiseEhrhartQuasiPolynomial(A)
print(peqp)

from ehrhart_quasi_polynomial.ehrhart_piecewise import create_polytope_from_matrix
from itertools import combinations_with_replacement
ehr = lambda b: len(create_polytope_from_matrix(A, b).integral_points())

for b in combinations_with_replacement(range(-10, 10), peqp._amb_dim):
    e = ehr(b)
    p = peqp(b)
    if e != p:
        print(b, e, p)
print(b)