from itertools import combinations_with_replacement

from .ehrhart_quasi_polynomial import get_period, get_gcd

# from sage.arith.functions import lcm
from sage.functions.other import factorial
from sage.geometry.cone import Cone
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.interfaces.gfan import gfan
from sage.matrix.constructor import Matrix as create_matrix
from sage.matrix.special import identity_matrix
from sage.modules.free_module_element import free_module_element
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.groebner_fan import PolyhedralFan
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.structure.element import Matrix


class PiecewiseEhrhartQuasiPolynomial():
    def __init__(self, A):
        self._A = A
        self._sec_fan = self._secondary_fan()

        self._amb_dim = self._sec_fan.ambient_dim()
        self._orth_dim = self._amb_dim - self._sec_fan.lineality_dim()
        self._origin = free_module_element([0]*self._amb_dim)
        self._amb_lattice = IntegralLattice(create_matrix(ZZ, self._amb_dim,
                                                          lambda x, y: int(x == y)))

        self._compute_change_of_basis_matrices()
        self._cone_dicts = self._get_cone_dicts(self._sec_fan)
        self._compute_piecewise()

    def _secondary_fan(self):
        """
        Returns the secondary fan associated with ``self._A``
        For details see the documentation of the gfan software package:
            https://users-math.au.dk/~jensen/software/gfan/gfan.html
        """
        gfan_input = "{" + ", ".join(str(row) for row in self._A.rows()) + "}"
        return PolyhedralFan(gfan(gfan_input, "secondaryfan"))

    def _compute_change_of_basis_matrices(self):
        orth_vectors = []
        for vec in self._sec_fan.fan_dict["ORTH_LINEALITY_SPACE"]:
            orth_vectors.append(free_module_element([int(val) for val in vec.split(" ")]))
        self._orth_vectors = orth_vectors
        self._orth_dim = len(orth_vectors)

        lin_vectors = []
        for vec in self._sec_fan.fan_dict["LINEALITY_SPACE"]:
            lin_vectors.append(free_module_element([int(val) for val in vec.split(" ")]))
        self._lin_vectors = lin_vectors
        self._lin_lattice = self._amb_lattice.sublattice(self._lin_vectors)

        self.change_of_basis_matrix = create_matrix(orth_vectors + lin_vectors)
        self.change_of_basis_inverse = self.change_of_basis_matrix.inverse()

    def _get_cone_dicts(self, sec_fan):
        fan_rays = [free_module_element(ray) for ray in self._sec_fan.rays()]
        ray_scalars = self._compute_ray_scalars(fan_rays)

        dictionaries = []
        for ray_lists in self._sec_fan.maximal_cones().values():
            for ray_list in ray_lists:
                cone_dic = {}
                cone_dic["rays"] = [fan_rays[idx] for idx in ray_list]
                cone_dic["scaled_rays"] = [ray_scalars[idx]*fan_rays[idx]
                                           for idx in ray_list]

                lattice_basis = cone_dic["scaled_rays"] + self._lin_vectors
                cone_dic["lattice"] = self._amb_lattice.sublattice(lattice_basis)
                cone_dic["quotient"] = self._amb_lattice.quotient(cone_dic["lattice"])

                cone_dic["cone"] = Cone(cone_dic["rays"])
                dictionaries.append(cone_dic)
        return dictionaries

    def _compute_ray_scalars(self, rays):
        ray_scalars = []
        for ray in rays:
            den = get_period(self._create_polytope_from_matrix(ray).Vrepresentation())
            ray_scalars.append(den)
        return ray_scalars

    def _compute_piecewise(self):
        num_variables = self._amb_dim - self._sec_fan.lineality_dim()
        max_degree = num_variables

        R = PolynomialRing(QQ, "x", num_variables)
        x_vars = R.gens()

        for idx, cone_dict in enumerate(self._cone_dicts):
            scaled_rays = cone_dict["scaled_rays"]

            needed_points = factorial(num_variables + max_degree)//(
                factorial(num_variables)*factorial(max_degree) )

            lattice_basis = [free_module_element(int(k==dim)*invar for k in range(num_variables))
                             for dim, invar in enumerate(cone_dict["quotient"].invariants())]
            lattice_points = self._generate_lattice_points(lattice_basis, needed_points)

            polynomials = {}
            for off_set in cone_dict["quotient"]:
                off_set = tuple(off_set)
                off_cone_points = [p + free_module_element(off_set) for p in lattice_points]

                num_integral_points = []
                for point in off_cone_points:
                    lifted = sum(p*scaled_rays[idx] for idx, p in enumerate(point))
                    polytope = self._create_polytope_from_matrix(lifted)
                    num_integral_points.append(len(polytope.integral_points()))
                    
                interpolated = R.interpolation(max_degree, off_cone_points,
                                               num_integral_points)
                polynomials[off_set] = interpolated#(*self._projection(x_vars))
            self._cone_dicts[idx]["polynomials"] = polynomials

    def _rat_period(self, ray):
        original_polytope = self._create_polytope_from_matrix(ray)
        period = get_period(original_polytope.Vrepresentation())
        scaled_polytope = self._create_polytope_from_matrix([period*r for r in ray])
        gcd = get_gcd(scaled_polytope.Vrepresentation())
        return period/gcd

    def _generate_cone_points(self, scaled_cone_rays, number):
        points = [self._origin]
        points_len = 0
        index = 1
        while points_len <= number:
            new_points = [sum(combi) for combi
                          in combinations_with_replacement(scaled_cone_rays, index)]
            points += new_points
            index += 1
            points_len += len(new_points)
        return points[:number]

    def _generate_lattice_points(self, basis_vecs, number):
        sub_points = []
        points_len = 0
        index = 1
        while points_len <= number:
            new_points = [sum(combi) for combi
                          in combinations_with_replacement(basis_vecs, index)]
            sub_points += new_points
            index += 1
            points_len += len(new_points)
        return sub_points[:number]


    def _create_polytope_from_matrix(self, b):
        """
        Unsafe version of ``create_polytope_from_matrix``
        to increase performance during lengthy computations.
        """
        inequalities = [[b[k]] + list(-self._A.rows()[k]) for k in range(self._A.nrows())]
        return Polyhedron(ieqs = inequalities)

    def _projected_coordinates(self, point):
        return (free_module_element(point)*self.change_of_basis_inverse)[:self._orth_dim]

    def _projection(self, point):
        new_representation = self._projected_coordinates(point)
        eval_point = sum(new_representation[k]*o_vec
                          for k, o_vec in enumerate(self._orth_vectors))
        return eval_point

    def __call__(self, point):
        return self.evaluate(point)

    def evaluate(self, point):
        if len(point) != self._amb_dim:
            raise ValueError("Dimension of ``point`` needs to be equal to the ambient"
                             f" dimension of ``self`` which is {self._amb_dim}.")

        proj = self._projection(point)
        for cone_dict in self._cone_dicts:
            if proj in cone_dict["cone"]:
                off_set = tuple(cone_dict["quotient"](point))
                return cone_dict["polynomials"][off_set](tuple(self._projected_coordinates(point)))
        return 0

    def investigate(self, point):
        if len(point) != self._amb_dim:
            raise ValueError("Dimension of ``point`` needs to be equal to the ambient"
                             f" dimension of ``self`` which is {self._amb_dim}.")

        proj = self._projection(point)
        for k, cone_dict in enumerate(self._cone_dicts):
            if proj in cone_dict["cone"]:
                off_set = tuple(cone_dict["quotient"](point))
                return (k, proj, off_set, cone_dict["polynomials"][off_set],
                        cone_dict["polynomials"][off_set](tuple(self._projected_coordinates(point))))
        return 0


    def __repr__(self):
        return f"PiecewiseEhrhartQuasiPolynomial({self._A})"

    def matrix(self):
        return self._A

    def secondary_fan(self):
        return self._sec_fan

    def domains(self):
        return self._max_cones

    def functions(self):
        return self._qps


def create_polytope_from_matrix(A, b):
    """
    Return the polytope whose faces are defined by
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
        raise AttributeError("'b' needs to implement '__getitem__', " +
                             "needed for indexing: b[0]")
    if A.nrows() != len(b):
        raise ValueError("Dimensions of 'A' and 'b' need to be compatible, need "
                         f"A.nrows == len(b), but have {A.nrows()} != {len(b)}")
    inequalities = [[b[k]] + list(-A.rows()[k]) for k in range(A.nrows())]
    return Polyhedron(ieqs = inequalities)