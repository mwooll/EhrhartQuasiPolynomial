from integerperiodicfunction import IntegerPeriodicFunctionRing

from sage.arith.functions import lcm

from sage.categories.pushout import ConstructionFunctor
from sage.categories.commutative_rings import CommutativeRings

from sage.rings.ring import CommutativeRing

from sage.structure.coerce_maps import CallableConvertMap
from sage.structure.element import RingElement
from sage.structure.unique_representation import UniqueRepresentation



class QuasiPolynomialElement(RingElement):
    r"""
    An element of the Ring of Quasi-Polynomial.

    This class should not be used to construct elements, rather construct an
    instance of the parent class 'QuasiPolynomialRing' and let that
    construct the elements, as in the examples below.


    EXAMPLES::
        
        >>> from quasipolynomial import *
        >>> from sage.rings.rational_field import QQ
        >>> qpr = QuasiPolynomialRing(QQ)
        >>> qpr()
        QuasiPolynomialElement(Ring of Quasi Polynomials over Rational Field, [[0]])
        >>> print(qpr([[0, 1], 2, 3])) # doctest: +NORMALIZE_WHITESPACE
        QuasiPolynomial given by
        [0, 1] + [2]*k + [3]*k^2
    """
    def __init__(self, parent, coefficients=None):
        r"""
        INPUT:
            - parent : instance of 'QuasiPolynomialRing'
            - coefficients : iterable of elements of 'parent.base()' (default=None)
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
        r"""
        Return self evaluated at 'value'

        EXAMPLES::

            >>> from quasipolynomial import *
            >>> from sage.rings.rational_field import QQ
            >>> qpr = QuasiPolynomialRing(QQ)
            >>> q = qpr([2, [0, 1]])
            >>> q(0)
            2
            >>> q(1)
            3
            >>> q(2)
            2
        """
        result = 0
        for power, coef in enumerate(self._coefficients):
            result += coef(value) * value**power
        return result

    def coefficients(self):
        r"""
        Return the coefficients of self

        EXAMPLES::

            >>> from quasipolynomial import *
            >>> from sage.rings.rational_field import QQ
            >>> qpr = QuasiPolynomialRing(QQ)
            >>> qpr([2, [0, 1]]).coefficients()
            [IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [2]), IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [0, 1])]
        """
        return self._coefficients

    def degree(self):
        r"""
        Return the degree of self

        EXAMPLES::

            >>> from quasipolynomial import *
            >>> from sage.rings.rational_field import QQ
            >>> qpr = QuasiPolynomialRing(QQ)
            >>> qpr([2, [0, 1]]).degree()
            1
            >>> qpr().degree()
            0
        """
        return self._degree

    def period(self):
        r"""
        Return the period of self

        EXAMPLES::

            >>> from quasipolynomial import *
            >>> from sage.rings.rational_field import QQ
            >>> qpr = QuasiPolynomialRing(QQ)
            >>> qpr().period()
            1
            >>> qpr([[1, 2, 3], [1, 2]]).period()
            6
        """
        return self._period

    def __str__(self):
        function_str = "QuasiPolynomial given by \n"
        function_str += f"{self._coefficients[0].constants()}"
        for power, coef in enumerate(self._coefficients[1:]):
            function_str += f" + {coef.constants()}*k" + f"^{power+1}"*(power>0)
        return function_str

    def __repr__(self):
        coefficients = [coef.constants() for coef in self._coefficients]
        return f"QuasiPolynomialElement({self.parent()}, {coefficients})"

    def __eq__(self, other):
        r"""
        Return whether self and other are considered equal in the ring

        TESTS::

            >>> from quasipolynomial import *
            >>> from sage.rings.rational_field import QQ
            >>> qpr = QuasiPolynomialRing(QQ)
            >>> qpr.zero() == qpr() == qpr([0])
            True
            >>> qpr([[0, 1]]) == qpr([0, 1])
            False
        """
        return self._coefficients == other.coefficients()

    def __bool__(self):
        r"""
        Return whether self is a non-zero element of the ring.

        TESTS::

            >>> from quasipolynomial import *
            >>> from sage.rings.rational_field import QQ
            >>> qpr = QuasiPolynomialRing(QQ)
            >>> bool(qpr.zero())
            False
            >>> bool(qpr([[0, 1]]))
            True
        """
        return (self._degree != 0 or self._period != 1
                or bool(self._coefficients[0]))

    def _neg_(self):
        r"""
        Return the additive inverse of self

        TESTS::

            >>> from quasipolynomial import *
            >>> from sage.rings.rational_field import QQ
            >>> qpr = QuasiPolynomialRing(QQ)
            >>> -qpr([[0, 1], 2])
            QuasiPolynomialElement(Ring of Quasi Polynomials over Rational Field, [[0, -1], [-2]])
        """
        return self.__class__(self.parent(), [-c for c in self._coefficients])

    def _add_(self, other):
        r"""
        Ring addition

        TESTS::

            >>> from quasipolynomial import *
            >>> from sage.rings.rational_field import QQ
            >>> qpr = QuasiPolynomialRing(QQ)
            >>> q = qpr([[1, 2], [1, 2]])
            >>> q + 1
            QuasiPolynomialElement(Ring of Quasi Polynomials over Rational Field, [[2, 3], [1, 2]])
            >>> r = qpr([[1, 2, 3], [2, 1]])
            >>> q + r
            QuasiPolynomialElement(Ring of Quasi Polynomials over Rational Field, [[2, 4, 4, 3, 3, 5], [3]])
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
        r"""
        TESTS::

            >>> from quasipolynomial import *
            >>> from sage.rings.rational_field import QQ
            >>> qpr = QuasiPolynomialRing(QQ)
            >>> q = qpr([[1, 2], [1, 2]])
            >>> q - 1
            QuasiPolynomialElement(Ring of Quasi Polynomials over Rational Field, [[0, 1], [1, 2]])
            >>> r = qpr([[1, 2, 3], [2, 1]])
            >>> q - r
            QuasiPolynomialElement(Ring of Quasi Polynomials over Rational Field, [[0, 0, -2, 1, -1, -1], [-1, 1]])
        """
        return self.__add__(-other)

    def _mul_(self, other):
        r"""
        Ring multiplication and scalar multiplication

        TESTS::

            >>> from quasipolynomial import *
            >>> from sage.rings.rational_field import QQ
            >>> qpr = QuasiPolynomialRing(QQ)
            >>> q = qpr([[1, 2], [1, 2]])
            >>> q * 2
            QuasiPolynomialElement(Ring of Quasi Polynomials over Rational Field, [[2, 4], [2, 4]])
            >>> r = qpr([[1, 2, 3], [2, 1]])
            >>> q * r # doctest: +NORMALIZE_WHITESPACE
            QuasiPolynomialElement(Ring of Quasi Polynomials over Rational Field, 
                                   [[1, 4, 3, 2, 2, 6], [3, 6, 5, 4, 4, 8], [2]])
        """
        mul_coefficients = [self.parent().base().zero()]*(self._degree + other.degree() + 1)
        for self_power, self_coef in enumerate(self._coefficients):
            for  other_power, other_coef in enumerate(other.coefficients()):
                mul_coefficients[self_power + other_power] += self_coef*other_coef
        return self.__class__(self.parent(), mul_coefficients)


class QuasiPolynomialRing(UniqueRepresentation, CommutativeRing):
    Element = QuasiPolynomialElement
    def __init__(self, base_ring, num_variables=1):
        base = IntegerPeriodicFunctionRing(base_ring)
        if base not in CommutativeRings():
            raise ValueError(f"{base} is not a commutative ring.")
        CommutativeRing.__init__(self, base)

    def _repr_(self):
        return f"Ring of Quasi Polynomials over {self.base().base_ring()}"

    def base_ring(self):
        r"""
        Return the base ring of the underlying IntegerPeriodicFunctionRing.

        EXAMPLES::

            >>> from quasipolynomial import *
            >>> from sage.rings.rational_field import QQ
            >>> QuasiPolynomialRing(QQ).base_ring()
            Rational Field
            >>> from sage.symbolic.ring import SR
            >>> QuasiPolynomialRing(SR).base_ring()
            Symbolic Ring
        """
        return self.base().base_ring()

    def characteristic(self):
        r"""
        Return the characteristic of the underlying IntegerPeriodicFunctionRing.

        EXAMPLES::

            >>> from quasipolynomial import *
            >>> from sage.rings.rational_field import QQ
            >>> QuasiPolynomialRing(QQ).characteristic()
            0
            >>> from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
            >>> QuasiPolynomialRing(IntegerModRing(19)).characteristic()
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
        r"""
        Return whether there is a coercion map from S to self.
        If so "self(s)" should work for all s in S
        """
        if self.base().base_ring().has_coerce_map_from(S):
            return True
        if isinstance(S, IntegerPeriodicFunctionRing):
            return True

    def construction(self):
        r"""
        Return a ConstructionFunctor
        """
        return QuasiPolynomialFunctor(self._reduction[1][1:], self._reduction[2]), self.base()

    def is_integral_domain(self):
        r"""
        Return whether self is an integral domain, which is False

        TESTS::

            >>> from integerperiodicfunction import *
            >>> from sage.rings.rational_field import QQ
            >>> ipfr = IntegerPeriodicFunctionRing(QQ)
            >>> ipfr.is_integral_domain()
            False
        """
        return False

    def is_unique_factorization_domain(self):
        r"""
        Return whether self is a unique factorization domain (UFD), which is False

        TESTS::

            >>> from integerperiodicfunction import *
            >>> from sage.rings.rational_field import QQ
            >>> ipfr = IntegerPeriodicFunctionRing(QQ)
            >>> ipfr.is_unique_factorization_domain()
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


def run_tests():
    run_doctests()

    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
    from sage.rings.integer_ring import ZZ
    from sage.rings.rational_field import QQ
    from sage.symbolic.ring import SR
    
    run_TestSuites([IntegerModRing(19), ZZ, QQ, SR])

def run_doctests():
    print("Doctests:\n")
    import doctest

    doctest.testmod()

def run_TestSuites(base_rings):
    print("\n\nTestSuites\n")
    from sage.misc.sage_unittest import TestSuite

    for ring in base_rings:
        if ring in CommutativeRings():
            print(f"\nTesting IPFR with {ring}")
            ipfr = QuasiPolynomialRing(ring)
            TestSuite(ipfr)

if __name__ == "__main__":
    run_tests()