from sage.matrix.constructor import Matrix

from ehrhart_quasi_polynomial import (ehrhart_quasi_polynomial,
                                      create_polytope_from_matrix,
                                      PiecewiseEhrhartQuasiPolynomial)

if __name__ == "__main__":
    A = Matrix([[-1, 0], [0, -1], [1, 2], [0, 1]])
    b = [0, 0, 5, 5]
    # ehr = ehrhart_quasi_polynomial(create_polytope_from_matrix(A, b).Vrepresentation())
    # print(ehr)

    peqp = PiecewiseEhrhartQuasiPolynomial(A)
    print(peqp._qps)
    # print(peqp(b))
