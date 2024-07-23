from sage.arith.functions import lcm
from sage.arith.misc import factor

from sage.categories.pushout import ConstructionFunctor
from sage.categories.rings import Rings

import doctest
from sage.misc.sage_unittest import TestSuite

from sage.rings.ring import Ring

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import RingElement


###     Element     ###
class IntegerPeriodicFunctionElement(RingElement):
    r"""
    An element of the Ring of Integer Periodic Functions.
    Throughout
    
    This class should not be used to construct elements, rather construct an
    instance of the parent class 'IntegerPeriodicFunctionRing' and let that
    construct the elements, as in the examples below.


    EXAMPLES::
        
        >>> from integerperiodicfunction import *
        >>> from sage.rings.rational_field import QQ
        >>> ipfr = IntegerPeriodicFunctionRing(QQ)
        >>> ipfr.zero()
        IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [0])
        >>> print(ipfr([1, 2, 3]))
        IntegerPeriodicFunction over Rational Field given by
            1 if k%3 == 0
            2 if k%3 == 1
            3 if k%3 == 2
    """
    def __init__(self, parent, constants=None):
        r"""
        INPUT:
            - parent : instance of "IntegerPeriodicFunctionRing"
            - constants : iterable of (default=None)
        """
        base = parent.base()
        if constants is None:
            constants = [base(0)]
        elif not hasattr(constants, "__iter__"):
            constants = [base(constants)]
        else:
            constants = [base(c) for c in constants]
        self.period = self._calculate_period(constants)
        self.constants = constants[:self.period]
        RingElement.__init__(self, parent)


    # methods called in __init__
    def _calculate_period(self, constants):
        r"""
        Computes the period of the Integer Periodic Function.

        
        """
        length = len(constants)

        if length == 1:
            return 1

        period = length
        factors = factor(length)
        for prime, power in factors._Factorization__x:
            div = self._get_period_divisor(constants, length, prime, power)
            period //= div
        
        return period

    def _get_period_divisor(self, constants, length, prime, power):
        last_div = 1
        fact = prime
        for k in range(power):
            if constants[:length//fact]*fact != constants:
                return last_div
            
            last_div = fact
            fact *= prime
        return last_div


    # use class elements as functions
    def __call__(self, k):
        r"""
        Returns the value at k%period.
        
        INPUT:
            - k : integer


        EXAMPLES::
            
            >>> from integerperiodicfunction import *
            >>> from sage.rings.rational_field import QQ
            >>> ipfr = IntegerPeriodicFunctionRing(QQ)
            >>> p = ipfr([1, 2, 3])
            >>> p(1)
            2
            >>> p(11)
            3
            >>> p(-6)
            1
        """
        try:
            return self.constants[k%self.period]
        except TypeError:
            raise TypeError("Integer periodic functions can only be evaluated at integer values.")


    def __getitem__(self, k):
        r"""
        Returns the value at k%period or for all elements in the slice.

        INPUT:
            - k : integer or slice

        EXAMPLES::
            
            >>> from integerperiodicfunction import *
            >>> from sage.rings.rational_field import QQ
            >>> ipfr = IntegerPeriodicFunctionRing(QQ)
            >>> p = ipfr([1, 2, 3])
            >>> p[1]
            2
            >>> p[99]
            1
            >>> p[0:6]
            [1, 2, 3, 1, 2, 3]
            >>> p[3:0:-1]
            [1, 3, 2]
        """
        if isinstance(k, slice):
            step = k.step if k.step else 1
            return [self.constants[i%self.period] for i in range(k.start, k.stop, step)]

        try:
            return self.constants[k%self.period]
        except TypeError:
            raise TypeError("Integer periodic functions can only be evaluated at integer values.")


    # representation
    def _repr_(self):
        return f"IntegerPeriodicFunctionElement({self.parent()}, {self.constants})"

    def __str__(self):
        function_str = f"IntegerPeriodicFunction over {self.parent().base_ring()} given by"
        for ind, val in enumerate(self.constants):
            function_str += f"\n\t{val} if k%{self.period} == {ind}"
        return function_str

    """
    def coefficient_repr(self, variable):
        if self.period == 1:
            if self.constants[0] == 0:
                return ""
            elif self.constants[0] == 1:
                return variable
            elif self.constants[0] == -1:
                return f"-{variable}"
            return f"{self.constants[0]}{variable}"
        else:
            return f"IPF{self.constants}({variable})"
    """

    # comparison
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.period == other.period and self.constants == other.constants
        else:
            return self.period == 1 and self.constants[0] == other

        return False

    # math support
    def __bool__(self):
        r"""
        Returns whether self is the zero element of the ring.
        
        Examples::

            >>> from integerperiodicfunction import *
            >>> from sage.rings.rational_field import QQ
            >>> ipfr = IntegerPeriodicFunctionRing(QQ)
            >>> bool(ipfr.zero())
            False
            >>> bool(ipfr([1, 2, 3]))
            True
        """
        return not (self.period == 1 and self.constants[0] == 0)

    def _neg_(self):
        r"""
        Returns the additive inverse of self

        Tests::

            >>> from integerperiodicfunction import *
            >>> from sage.rings.rational_field import QQ
            >>> ipfr = IntegerPeriodicFunctionRing(QQ)
            >>> -ipfr([1, 2, 3])
            IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [-1, -2, -3])
        """
        return self.__class__(self.parent(), [-c for c in self.constants])

    def _add_(self, other):
        r"""
        Ring addition
        
        Tests::

            >>> from integerperiodicfunction import *
            >>> from sage.rings.rational_field import QQ
            >>> ipfr = IntegerPeriodicFunctionRing(QQ)
            >>> ipfr([1, 2, 3]) + ipfr([0, 1])
            IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [1, 3, 3, 2, 2, 4])
            >>> ipfr([1, 2, 3]) + 2
            IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [3, 4, 5])
            >>> ipfr([1, 2, 3]) + ipfr([3, 2, 1])
            IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [4])
        """
        add_period = lcm(self.period, other.period)
        add_constants = [self[k] + other[k] for k in range(add_period)]
        return self.__class__(self.parent(), add_constants)

    def _radd_(self, other):
        r"""
        Tests::

            >>> from integerperiodicfunction import *
            >>> from sage.rings.rational_field import QQ
            >>> ipfr = IntegerPeriodicFunctionRing(QQ)
            >>> 2 + ipfr([1, 2, 3]) 
            IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [3, 4, 5])
        """
        return self.__add__(other)

    def _sub_(self, other):
        r"""
        Tests::

            >>> from integerperiodicfunction import *
            >>> from sage.rings.rational_field import QQ
            >>> ipfr = IntegerPeriodicFunctionRing(QQ)
            >>> ipfr([1, 2, 3]) - ipfr([0, 1])
            IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [1, 1, 3, 0, 2, 2])
            >>> ipfr([1, 2, 3]) - 2
            IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [-1, 0, 1])
        """
        return self.__add__(-other)

    def _rsub_(self, other):
        r"""
        Tests::

            >>> from integerperiodicfunction import *
            >>> from sage.rings.rational_field import QQ
            >>> ipfr = IntegerPeriodicFunctionRing(QQ)
            >>> 2 - ipfr([1, 2, 3]) 
            IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [1, 0, -1])
        """
        return (-self).__add__(other)

    def _mul_(self, other):
        r"""
        Ring multiplication and scalar multiplication
        
        Tests::

            >>> from integerperiodicfunction import *
            >>> from sage.rings.rational_field import QQ
            >>> ipfr = IntegerPeriodicFunctionRing(QQ)
            >>> ipfr([1, 2, 3]) * ipfr([0, 1])
            IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [0, 2, 0, 1, 0, 3])
            >>> ipfr([1, 2, 4]) * ipfr([4, 2, 1])
            IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [4])
            >>> ipfr([1, 2, 3]) * 2
            IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [2, 4, 6])
        """
        mul_period = lcm(self.period, other.period)
        mul_constants = [self[k]*other[k] for k in range(mul_period)]
        return self.__class__(self.parent(), mul_constants)

    def _rmul_(self, other):
        r"""
        Tests::

            >>> from integerperiodicfunction import *
            >>> from sage.rings.rational_field import QQ
            >>> ipfr = IntegerPeriodicFunctionRing(QQ)
            >>> 2 * ipfr([1, 0, 1])
            IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [2, 0, 2])
        """
        return self.__mul__(other)


###     Parent     ###
class IntegerPeriodicFunctionRing(UniqueRepresentation, Ring):
    # needed for automatic coercion
    Element = IntegerPeriodicFunctionElement
    def __init__(self, base):
        if base not in Rings():
            raise ValueError(f"{base} is not a ring.")
        Ring.__init__(self, base)                         
            
    def _repr_(self):
        return f"Ring of Integer Periodic Functions over {self.base()}"

    def base_ring(self):
        r"""
        Returns the base_ring of the base of self

        Examples::

            >>> from integerperiodicfunction import *
            >>> from sage.rings.rational_field import QQ
            >>> IntegerPeriodicFunctionRing(QQ).base_ring()
            Rational Field
            >>> from sage.symbolic.ring import SR
            >>> IntegerPeriodicFunctionRing(SR).base_ring()
            Symbolic Ring
        """
        return self.base().base_ring()

    def characteristic(self):
        r"""
        Returns the characteristic of the base of self

        Examples::

            >>> from integerperiodicfunction import *
            >>> from sage.rings.rational_field import QQ
            >>> IntegerPeriodicFunctionRing(QQ).characteristic()
            0
            >>> from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
            >>> IntegerPeriodicFunctionRing(IntegerModRing(19)).characteristic()
            19
        """
        return self.base().characteristic()

    def _element_constructor_(self, *args, **kwds):
        r"""
        Handles the automatic coercion of objects to IntegerPeriodicFunctionElement

        Examples::

            >>> from integerperiodicfunction import *
            >>> from sage.rings.rational_field import QQ
            >>> ipfr = IntegerPeriodicFunctionRing(QQ)
            >>> ipfr(1)
            IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [1])
            >>> ipfr([1, 2, 3])
            IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [1, 2, 3])
            >>> ipfr([4, 4, 4, 4])
            IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [4])
            >>> ipfr()
            IntegerPeriodicFunctionElement(Ring of Integer Periodic Functions over Rational Field, [0])
        """
        if len(args) != 1:
            return self.element_class(self, *args, **kwds)
        x = args[0]
        try:
            if hasattr(x, "__iter__"):
                P = x[0].parent()
            else:
                P = x.parent()
            if P.is_subring(self):
                P = self
        except AttributeError:
            P = self
        return self.element_class(P, x, **kwds)

    def _coerce_map_from_(self, S):
        r"""
        Returns whether there is a coercion map from S to self.
        If so "self(s)" should work for all s in S

        Tests::

            >>> from integerperiodicfunction import *
            >>> from sage.rings.rational_field import QQ
            >>> ipfr = IntegerPeriodicFunctionRing(QQ)
            >>> ipfr._coerce_map_from_(QQ)
            True
            >>> from sage.rings.integer_ring import ZZ
            >>> ipfr._coerce_map_from_(ZZ)
            True
            >>> from sage.symbolic.ring import SR
            >>> print(ipfr._coerce_map_from_(SR))
            None
        """
        if self.base().has_coerce_map_from(S):
            return True

    def construction(self):
        r"""
        Returns a ConstructionFunctor
        """
        return IntegerPeriodicFunctionFunctor(self._reduction[1][1:], self._reduction[2]), self.base()


###     Functor     ###
class IntegerPeriodicFunctionFunctor(ConstructionFunctor):
    rank = 10
    def __init__(self, args=None, kwds=None):
        self.args = args or ()
        self.kwds = kwds or {}
        ConstructionFunctor.__init__(self, Rings(), Rings())

    def _apply_functor(self, R):
        return IntegerPeriodicFunctionRing(R, *self.args, **self.kwds)

    def merge(self, other):
        if isinstance(other, type(self)):
            return self


if __name__ == "__main__":
    from sage.rings.rational_field import QQ
    
    ipfr = IntegerPeriodicFunctionRing(QQ)
    TestSuite(ipfr).run()
    
    doctest.testmod()