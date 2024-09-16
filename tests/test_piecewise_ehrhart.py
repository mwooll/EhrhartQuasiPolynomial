from sage.matrix.constructor import Matrix
from itertools import combinations_with_replacement

from ehrhart_quasi_polynomial import (ehrhart_quasi_polynomial,
                                      create_polytope_from_matrix,
                                      PiecewiseEhrhartQuasiPolynomial)

if __name__ == "__main__":
    A = Matrix([[-1, 0], [0, -1], [1, 1], [0, 1]])
    ehr = lambda b: ehrhart_quasi_polynomial(create_polytope_from_matrix(A, b).Vrepresentation())(1)

    peqp = PiecewiseEhrhartQuasiPolynomial(A)
    
    # print(peqp._sec_fan.fan_dict)
    # print()
    # print(peqp._cone_dicts)
    
    

    for p in combinations_with_replacement(range(-5, 4), 4):
        expected = ehr(p)
        actual = peqp(p)
        if actual != expected:
            print(p, expected, actual)
            print(peqp.investigate(p))
