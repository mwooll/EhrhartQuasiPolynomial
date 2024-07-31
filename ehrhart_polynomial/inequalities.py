import sage.all
from sage.calculus.var import var
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.symbolic.relation import solve

from sage.matrix.constructor import Matrix as create_matrix
from sage.modules.free_module_element import free_module_element # creates a vector
from sage.structure.element import Vector, Matrix # baseclasses


def create_polyhedron(A, b):
    """
    creates a polyhedron whose faces are defined by
        Ax <= b
    where the inequality is understood componentwise
    """
    if not isinstance(A, Matrix):
        raise TypeError("A needs to be an instance of sage.structure.element"
                        f".Matrix, but has type {type(A)}")
    if not hasattr(b, "__getitem__"):
        raise AttributeError("b needs to implement '__getitem__', " +
                             "needed for indexing: b[0]")
    if A.nrows() != len(b):
        raise ValueError("Dimensions of A and b need to be compatible, need "+
                         f"A.nrows = len(b), but have {A.nrows()} != {len(b)}")

    inequalities = [[b[k]] + list(-A.rows()[k]) for k in range(A.nrows())]
    print(inequalities)
    return Polyhedron(ieqs = inequalities)

