from .integerperiodicfunction import IntegerPeriodicFunctionRing

from sage.arith.functions import lcm
from sage.categories.pushout import ConstructionFunctor
from sage.categories.commutative_rings import CommutativeRings
from sage.misc.cachefunc import cached_method
from sage.rings.ring import CommutativeRing
from sage.structure.element import RingElement
from sage.structure.unique_representation import UniqueRepresentation


class QuasiPolynomialElement(RingElement):
    """
    An element of the Ring of Quasi-Polynomial.

    This class should not be used to construct elements, rather construct an
    instance of the parent class ``QuasiPolynomialRing`` and let that
    construct the elements, as in the examples below.


    EXAMPLES::
        
        sage: from ehrhart_quasi_polynomial.quasipolynomial import QuasiPolynomialRing
        sage: qpr = QuasiPolynomialRing(QQ)
        sage: qpr()
        QuasiPolynomialElement(Ring of Quasi-Polynomials over Rational Field, [[0]])
        sage: print(qpr([[0, 1], 2, 3])) # doctest: +NORMALIZE_WHITESPACE
        QuasiPolynomial given by
        [0, 1] + [2]*t + [3]*t^2
    """
    def __init__(self, parent, coefficients=None):
        """
        INPUT:

        - ``parent`` -- instance of ``QuasiPolynomialRing``
        - ``coefficients`` -- iterable of elements of ``parent.base()`` (default : ``None``)
        """
        base = parent.base()
        if coefficients is None:
            coefficients = [base(0)]
        elif not hasattr(coefficients, "__iter__"):
            coefficients = [base(coefficients)]
        else:
            coefficients = [base(c) for c in coefficients]
        self._coefficients, self._degree = self._reduce_coefficients(coefficients)
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

            sage: from ehrhart_quasi_polynomial.quasipolynomial import QuasiPolynomialRing
            sage: qpr = QuasiPolynomialRing(QQ)
            sage: q = qpr([2, [0, 1]])
            sage: q(0)
            2
            sage: q(1)
            3
            sage: q(2)
            2
        """
        result = 0
        for power, coef in enumerate(self._coefficients):
            result += coef(value) * value**power
        return result

    def coefficients(self):
        """
        Return the coefficients of ``self``.

        EXAMPLES::

            sage: from ehrhart_quasi_polynomial.quasipolynomial import QuasiPolynomialRing
            sage: qpr = QuasiPolynomialRing(QQ)
            sage: qpr([2, [0, 1]]).coefficients()
            [IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [2]), IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [0, 1])]
        """
        return self._coefficients

    def degree(self):
        """
        Return the degree of ``self``.

        EXAMPLES::

            sage: from ehrhart_quasi_polynomial.quasipolynomial import QuasiPolynomialRing
            sage: qpr = QuasiPolynomialRing(QQ)
            sage: qpr([2, [0, 1]]).degree()
            1
            sage: qpr().degree()
            0
        """
        return self._degree

    def period(self):
        """
        Return the period of ``self``.

        EXAMPLES::

            sage: from ehrhart_quasi_polynomial.quasipolynomial import QuasiPolynomialRing
            sage: qpr = QuasiPolynomialRing(QQ)
            sage: qpr().period()
            1
            sage: qpr([[1, 2, 3], [1, 2]]).period()
            6
        """
        return self._period

    def __str__(self):
        function_str = "QuasiPolynomial given by \n"
        function_str += f"{self._coefficients[0].constants()}"
        for power, coef in enumerate(self._coefficients[1:]):
            function_str += f" + {coef.constants()}*t" + f"^{power+1}"*(power>0)
        return function_str

    def __repr__(self):
        coefficients = [coef.constants() for coef in self._coefficients]
        return f"QuasiPolynomialElement({self.parent()}, {coefficients})"

    def __eq__(self, other):
        """
        Return whether ``self`` and ``other`` are considered equal.

        TESTS::

            sage: from ehrhart_quasi_polynomial.quasipolynomial import QuasiPolynomialRing
            sage: qpr = QuasiPolynomialRing(QQ)
            sage: qpr.zero() == qpr() == qpr([0])
            True
            sage: qpr([[0, 1]]) == qpr([0, 1])
            False
            sage: qpr([1]) == 1
            True
        """
        if isinstance(other, self.__class__):
            return self._coefficients == other.coefficients()
        elif hasattr(other, "coefficients"):
            return self._coefficients == other.coefficients()
        else:
            return self._degree == 0 and self._coefficients[0] == other

    def __bool__(self):
        """
        Return whether ``self`` is a non-zero element of the ring.

        TESTS::

            sage: from ehrhart_quasi_polynomial.quasipolynomial import QuasiPolynomialRing
            sage: qpr = QuasiPolynomialRing(QQ)
            sage: bool(qpr.zero())
            False
            sage: bool(qpr([[0, 1]]))
            True
        """
        return (self._degree != 0 or self._period != 1
                or bool(self._coefficients[0]))

    def _neg_(self):
        """
        Return the additive inverse of ``self``.

        TESTS::

            sage: from ehrhart_quasi_polynomial.quasipolynomial import QuasiPolynomialRing
            sage: qpr = QuasiPolynomialRing(QQ)
            sage: -qpr([[0, 1], 2])
            QuasiPolynomialElement(Ring of Quasi-Polynomials over Rational Field, [[0, -1], [-2]])
        """
        return self.__class__(self.parent(), [-c for c in self._coefficients])

    def _add_(self, other):
        """
        Ring addition

        TESTS::

            sage: from ehrhart_quasi_polynomial.quasipolynomial import QuasiPolynomialRing
            sage: qpr = QuasiPolynomialRing(QQ)
            sage: q = qpr([[1, 2], [1, 2]])
            sage: q + 1
            QuasiPolynomialElement(Ring of Quasi-Polynomials over Rational Field, [[2, 3], [1, 2]])
            sage: r = qpr([[1, 2, 3], [2, 1]])
            sage: q + r
            QuasiPolynomialElement(Ring of Quasi-Polynomials over Rational Field, [[2, 4, 4, 3, 3, 5], [3]])
        """
        other_degree = other.degree()
        other_coefs = other.coefficients()

        len_coefs = max(self._degree, other_degree)+1
        add_coefficients = [self.parent().base().zero()]*len_coefs
        for idx, coef in enumerate(self._coefficients):
            add_coefficients[idx] += coef
        for idx, coef in enumerate(other_coefs):
            add_coefficients[idx] += coef
        return self.__class__(self.parent(), add_coefficients)

    def _sub_(self, other):
        """
        TESTS::

            sage: from ehrhart_quasi_polynomial.quasipolynomial import QuasiPolynomialRing
            sage: qpr = QuasiPolynomialRing(QQ)
            sage: q = qpr([[1, 2], [1, 2]])
            sage: q - 1
            QuasiPolynomialElement(Ring of Quasi-Polynomials over Rational Field, [[0, 1], [1, 2]])
            sage: r = qpr([[1, 2, 3], [2, 1]])
            sage: q - r
            QuasiPolynomialElement(Ring of Quasi-Polynomials over Rational Field, [[0, 0, -2, 1, -1, -1], [-1, 1]])
        """
        return self.__add__(-other)

    def _mul_(self, other):
        """
        Ring multiplication and scalar multiplication

        TESTS::

            sage: from ehrhart_quasi_polynomial.quasipolynomial import QuasiPolynomialRing
            sage: qpr = QuasiPolynomialRing(QQ)
            sage: q = qpr([[1, 2], [1, 2]])
            sage: q * 2
            QuasiPolynomialElement(Ring of Quasi-Polynomials over Rational Field, [[2, 4], [2, 4]])
            sage: r = qpr([[1, 2, 3], [2, 1]])
            sage: q * r # doctest: +NORMALIZE_WHITESPACE
            QuasiPolynomialElement(Ring of Quasi-Polynomials over Rational Field, 
                                   [[1, 4, 3, 2, 2, 6], [3, 6, 5, 4, 4, 8], [2]])
        """
        mul_coefficients = [self.parent().base().zero()]*(self._degree + other.degree() + 1)
        for self_power, self_coef in enumerate(self._coefficients):
            for  other_power, other_coef in enumerate(other.coefficients()):
                mul_coefficients[self_power + other_power] += self_coef*other_coef
        return self.__class__(self.parent(), mul_coefficients)


class QuasiPolynomialRing(UniqueRepresentation, CommutativeRing):
    """
    
    """
    Element = QuasiPolynomialElement
    def __init__(self, base_ring):
        """
        INPUT:

        - ``base_ring`` -- the base ring for the underlying ``IntegerPeriodicFunctionRing``,
            needs to be a commutative ring
        """
        base = IntegerPeriodicFunctionRing(base_ring)
        if base not in CommutativeRings():
            raise ValueError(f"{base} is not a commutative ring.")
        CommutativeRing.__init__(self, base)

    def _repr_(self):
        """
        TESTS::

            sage: from ehrhart_quasi_polynomial.quasipolynomial import QuasiPolynomialRing
            sage: QuasiPolynomialRing(QQ)
            Ring of Quasi-Polynomials over Rational Field
            sage: QuasiPolynomialRing(SR)
            Ring of Quasi-Polynomials over Symbolic Ring
        """
        return f"Ring of Quasi-Polynomials over {self.base().base_ring()}"

    @cached_method
    def gen(self, n=0):
        """
        Return the indeterminate generator of ``self``.

        EXAMPLES::

            sage: from ehrhart_quasi_polynomial.quasipolynomial import QuasiPolynomialRing
            sage: qpr = QuasiPolynomialRing(QQ)
            sage: x = qpr.gen(); x
            QuasiPolynomialElement(Ring of Quasi-Polynomials over Rational Field, [[0], [1]])

        An identical generator is always returned.

        ::

            sage: x is qpr.gen()
            True
        """
        if n != 0:
            raise IndexError("generator 'n' is not defined")
        return self.element_class(self, [0, 1])

    def base_ring(self):
        """
        Return the base ring of the IntegerPeriodicFunctionRing of ``self``.

        EXAMPLES::

            sage: from ehrhart_quasi_polynomial.quasipolynomial import QuasiPolynomialRing
            sage: QuasiPolynomialRing(QQ).base_ring()
            Rational Field
            sage: QuasiPolynomialRing(SR).base_ring()
            Symbolic Ring
        """
        return self.base().base_ring()

    def characteristic(self):
        """
        Return the characteristic of the IntegerPeriodicFunctionRing of ``self``.

        EXAMPLES::

            sage: from ehrhart_quasi_polynomial.quasipolynomial import QuasiPolynomialRing
            sage: QuasiPolynomialRing(QQ).characteristic()
            0
            sage: QuasiPolynomialRing(IntegerModRing(19)).characteristic()
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
        return QuasiPolynomialFunctor(self._reduction[1][1:], self._reduction[2]), self.base().base_ring()

    def is_integral_domain(self):
        """
        Return ``False``, since quasi-polynomial rings are never integral domains.

        TESTS::

            sage: from ehrhart_quasi_polynomial.quasipolynomial import QuasiPolynomialRing
            sage: qpr = QuasiPolynomialRing(QQ)
            sage: qpr.is_integral_domain()
            False
        """
        return False

    def is_unique_factorization_domain(self):
        """
        Return ``False``, since quasi-polynomial rings are never unique factorization domains (UFD).

        TESTS::

            sage: from ehrhart_quasi_polynomial.quasipolynomial import QuasiPolynomialRing
            sage: qpr = QuasiPolynomialRing(QQ)
            sage: qpr.is_unique_factorization_domain()
            False
        """
        return False


class QuasiPolynomialFunctor(ConstructionFunctor):
    rank = 10
    def __init__(self, args=None, kwds=None):
        self.args = args or ()
        self.kwds = kwds or {}
        ConstructionFunctor.__init__(self, CommutativeRings(), CommutativeRings())

    def _apply_functor(self, R):
        return QuasiPolynomialRing(R, *self.args, **self.kwds)

    def merge(self, other):
        if isinstance(other, type(self)):
            return self
