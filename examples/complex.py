from ehrhart_quasi_polynomial.ehrhart_piecewise import PiecewiseEhrhartQuasiPolynomial as PEQP

from sage.matrix.constructor import Matrix


if __name__ == "__main__":
    A = Matrix([[-1, 0], [0, -1], [1, 1], [0, 1]])
    p = PEQP(A)
    print(p._cone_dicts)