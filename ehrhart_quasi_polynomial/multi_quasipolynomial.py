from .integerperiodicfunction import IntegerPeriodicFunctionRing
from .quasipolynomial import QuasiPolynomialRing

from sage.arith.functions import lcm
from sage.categories.pushout import ConstructionFunctor
from sage.categories.commutative_rings import CommutativeRings
from sage.matrix.constructor import Matrix as create_matrix
from sage.misc.cachefunc import cached_method
from sage.rings.ring import CommutativeRing
from sage.structure.element import RingElement, Matrix
from sage.structure.unique_representation import UniqueRepresentation


class MultiQuasiPolynomialElement(RingElement):
    """
    An element of the Ring of Quasi-Polynomial.

    This class should not be used to construct elements, rather construct an
    instance of the parent class ``MultiQuasiPolynomialRing`` and let that
    construct the elements, as in the examples below.


    EXAMPLES::
        
        sage: from ehrhart_quasi_polynomial.multi_quasipolynomial import MultiQuasiPolynomialRing
        sage: mqpr = MultiQuasiPolynomialRing(QQ, 2)
        sage: mqpr({(0, 0): })
    """
    def __init__(self, parent, coefs):
        """
        INPUT:

        - ``parent`` -- instance of ``MultiQuasiPolynomialRing``
        - ``coef_matrix`` -- (sparse) dictionary {index: value}
        """
        self.base = parent.base()

        if not isinstance(coefs, dict):
            raise TypeError("`coefs`` must be an instance of `dict` "
                            f"but is of type {type(coefs)}")
        self._coef_matrix = []
        
        self._coefficients, self._degree = self._reduce_coefficients(coefs)
        self._period = self._calculate_peroid()

        RingElement.__init__(self, parent)

    def _reduce_coefficients(self, coefficients):
        degree = len(coefficients)-1
        while degree > 0 and coefficients[degree] == 0:
            degree -= 1

        coefficients = coefficients[:degree+1]
        return coefficients, degree

    def _calculate_peroid(self):
        periods = [coef.period() for coef in self._coefficients]
        period = lcm(periods)
        return period

    def __call__(self, value):
        """
        Return ``self`` evaluated at ``value``.

        EXAMPLES::

            sage: from ehrhart_quasi_polynomial.multi_quasipolynomial import MultiQuasiPolynomialRing
            sage: mqpr = MultiQuasiPolynomialRing(QQ, 2)

        """
        result = 0
        for power, coef in enumerate(self._coefficients):
            result += coef(value) * value**power
        return result

    def degree(self):
        """
        Return the degree of ``self``.
        """
        return self._degree

    def period(self):
        """
        Return the period of ``self``.
        """
        return self._period

    def __str__(self):
        function_str = "MultiQuasiPolynomial given by \n"
        function_str += f"{self._coefficients[0].constants()}"
        for power, coef in enumerate(self._coefficients[1:]):
            function_str += f" + {coef.constants()}*t" + f"^{power+1}"*(power>0)
        return function_str

    def __repr__(self):
        return f"MultiQuasiPolynomialElement({self.parent()}, {self._coef_matrix})"

    def __eq__(self, other):
        """
        Return whether ``self`` and ``other`` are considered equal in ``self.parent()``.
        """
        return NotImplemented

    def __bool__(self):
        """
        Return whether ``self`` is a non-zero element of the ring ``self.parent()``.
        """
        return NotImplemented

    def _neg_(self):
        """
        Return the additive inverse of ``self``.
        """
        return NotImplemented

    def _add_(self, other):
        """
        Ring addition
        """
        return NotImplemented

    def _sub_(self, other):
        """
        """
        return self.__add__(-other)

    def _mul_(self, other):
        """
        Ring multiplication and scalar multiplication
        """
        return NotImplemented


class MultiQuasiPolynomialRing(UniqueRepresentation, CommutativeRing):
    Element = MultiQuasiPolynomialElement
    def __init__(self, base_ring, num_var=1):
        """
        INPUT:

        - ``base_ring`` -- the base ring for the underlying ``IntegerPeriodicFunctionRing``,
            needs to be a commutative ring
        - ``num_var`` -- the number of variables (default : ``1``)
            if ``1`` is given, an instance of the (univariate) Quasi-Polynomial Ring will be returned.
        """
        if num_var == 1:
            return QuasiPolynomialRing(base_ring)

        self._ngens = num_var
        base = IntegerPeriodicFunctionRing(base_ring)
        if base not in CommutativeRings():
            raise ValueError(f"{base} is not a commutative ring.")
        CommutativeRing.__init__(self, base)

    def _repr_(self):
        """
        TESTS::

            sage: from ehrhart_quasi_polynomial.multi_quasipolynomial import MultiQuasiPolynomialRing
            sage: MultiQuasiPolynomialRing(QQ, 2)
            Ring of Multi Quasi-Polynomials in x0, x1 over Rational Field
            sage: MultiQuasiPolynomialRing(SR, 3)
            Ring of Multi Quasi-Polynomials in x0, x1, x2 over Symbolic Ring
        """
        return f"Ring of Multi Quasi-Polynomials in {self.gens()} over {self.base().base_ring()}"

    @cached_method
    def gen(self, n=0):
        """
        Return the 'n-th' generator of ``self``.
        """
        if n >= self._ngens:
            raise IndexError(f"'{n}-th' generator is not defined")
        return NotImplemented

    def ngens(self):
        """
        Return the number of generators of ``self``.
        """
        return self._ngens

    def gens(self):
        """
        Return the generators of ``self``.
        """
        return tuple(self.gen(i) for i in range(self._ngens))

    def base_ring(self):
        """
        Return the base ring of the IntegerPeriodicFunctionRing of ``self``.

        EXAMPLES::

            sage: from ehrhart_quasi_polynomial.multi_quasipolynomial import MultiQuasiPolynomialRing
            sage: MultiQuasiPolynomialRing(QQ, 2).base_ring()
            Rational Field
            sage: MultiQuasiPolynomialRing(SR, 3).base_ring()
            Symbolic Ring
        """
        return self.base().base_ring()

    def characteristic(self):
        """
        Return the characteristic of the IntegerPeriodicFunctionRing of ``self``.

        EXAMPLES::

            sage: from ehrhart_quasi_polynomial.multi_quasipolynomial import MultiQuasiPolynomialRing
            sage: MultiQuasiPolynomialRing(QQ, 2).characteristic()
            0
            sage: MultiQuasiPolynomialRing(IntegerModRing(19), 3).characteristic()
            19
        """
        return self.base().characteristic()

    def _element_constructor_(self, *args, **kwds):
        if len(args) != 1:
            return self.element_class(self, *args, **kwds)
        x = args[0]
        try:
            if hasattr(x, "constants"):
                x = [x.constants()]
                P = x[0].parent()
            elif hasattr(x, "__iter__"):
                P = x[0].parent()
            else:
                P = x.parent()
            if P.is_subring(self.base().base_ring()):
                P = self
        except AttributeError:
            P = self
        return self.element_class(P, x, **kwds)

    def _coerce_map_from_(self, S):
        """
        Return whether there is a coercion map from ``S`` to ``self``.
        If so ``self(s)`` should work for all ``s in S``.
        """
        if self.base().base_ring().has_coerce_map_from(S):
            return True
        if isinstance(S, IntegerPeriodicFunctionRing):
            return True

    def construction(self):
        """
        Return the construction functor corresponding to ``self``.
        """
        return MultiQuasiPolynomialFunctor(self._reduction[1][1:], self._reduction[2]), self.base().base_ring()

    def is_integral_domain(self):
        """
        Return ``False``, since quasi-polynomial rings are never integral domains.

        TESTS::

            sage: from ehrhart_quasi_polynomial.multi_quasipolynomial import MultiQuasiPolynomialRing
            sage: MultiQuasiPolynomialRing(QQ, 2).is_integral_domain()
            False
        """
        return False

    def is_unique_factorization_domain(self):
        """
        Return ``False``, since quasi-polynomial rings are never unique factorization domains (UFD).

        TESTS::

            sage: from ehrhart_quasi_polynomial.multi_quasipolynomial import MultiQuasiPolynomialRing
            sage: MultiQuasiPolynomialRing(QQ, 3).is_unique_factorization_domain()
            False
        """
        return False


class MultiQuasiPolynomialFunctor(ConstructionFunctor):
    rank = 10
    def __init__(self, args=None, kwds=None):
        self.args = args or ()
        self.kwds = kwds or {}
        ConstructionFunctor.__init__(self, CommutativeRings(), CommutativeRings())

    def _apply_functor(self, R):
        return MultiQuasiPolynomialRing(R, *self.args, **self.kwds)

    def merge(self, other):
        if isinstance(other, type(self)):
            return self
