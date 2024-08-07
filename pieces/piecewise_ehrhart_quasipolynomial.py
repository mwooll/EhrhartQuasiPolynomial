from sage.all import *
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.structure.element import Matrix

from ehrhart_quasi_polynomial import secondary_fan, ehrhart_quasi_polynomial



def piecewise_ehrhart_quasipolynomial(A_matrix, sec_fan=None):
    if sec_fan is None:
        sec_fan = secondary_fan(A_matrix, False, "gfan_files")

    components = _compute_piecewiese(A_matrix, sec_fan)
    return components
    
def _compute_piecewise(A_matrix, sec_fan):
    pass

def create_polytope_from_matrix(A, b):
    """
    Creates and returns the polytope whose faces are defined by
        Ax <= b
    where the inequality is understood componentwise

    EXAMPLES::

        sage: from ehrhart_quasi_polynomial import create_polytope_from_matrix
        sage: A = Matrix([[-1, 0], [0, -1], [1, 1]]); b = [0, 0, 1]
        sage: poly = create_polytope_from_matrix(A, b)
        sage: poly
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices
        sage: poly.Vrepresentation()
        (A vertex at (1, 0), A vertex at (0, 1), A vertex at (0, 0))
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
    return Polyhedron(ieqs = inequalities)

if __name__ == "__main__":
    A_matrix = [[1], [-1]]
    sec_fan = secondary_fan(A_matrix, False, "gfan_files")