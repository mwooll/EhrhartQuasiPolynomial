from sage.arith.functions import lcm
from sage.arith.misc import factor

from sage.categories.pushout import ConstructionFunctor
from sage.categories.rings import Rings

from sage.misc.sage_unittest import TestSuite

from sage.rings.ring import Ring
from sage.rings.rational_field import QQ

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import RingElement


###     Parent     ###
class IntegerPeriodicFunctionRing(UniqueRepresentation, Ring):
    def __init__(self, base):
        if base not in Rings():
            raise ValueError(f"{base} is not a ring.")
        Ring.__init__(self, base)                         
            
    def _repr_(self):
        return f"Ring of Integer Periodic Functions over {self.base()}"
    def base_ring(self):
        return self.base().base_ring()
    def characteristic(self):
        return self.base().characteristic()

    def _element_constructor_(self, *args, **kwds):
        if len(args) != 1:
            return self.element_class(self, *args, **kwds)
        x = args[0]
        try:
            P = x.parent()
            if P.is_subring(QQ):
                P = self
        except AttributeError:
            P = self
        return self.element_class(P, x, **kwds)

    def _coerce_map_from_(self, S):
        if self.base().has_coerce_map_from(S):
            return True

    def construction(self):
        return IntegerPeriodicFunctionFunctor(self._reduction[1][1:], self._reduction[2]), self.base()


###     Element     ###
class IntegerPeriodicFunctionElement(RingElement):
    def __init__(self, parent, constants=None):
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

    """
    methods called in __init__
    """
    def _calculate_period(self, constants):
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

    """
    use class object as function
    """
    def __call__(self, k):
        return self.constants[k%self.period]

    """
    standard dunder methods
    """
    def _repr_(self):
        return f"IntegerPeriodicFunctionElement({self.parent()}, {self.constants})"

    def _str_(self):
        function_str = "IntegerPeriodicFunction over {self.parent().base_ring()} given by"
        for ind, val in enumerate(self.constants):
            function_str += f"\n\t{val} if k%{self.period} == {ind}"
        return function_str

    def coefficient_repr(self, variable):
        if self.period == 1:
            return str(self.constants[0])
        else:
            return f"{repr(self)}({variable})"

    """
    comparison operators
    """
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.period == other.period and self.constants == other.constants
        else:
            return self.period == 1 and self.constants[0] == other

        return False

    """
    math support
    """
    def __bool__(self):
        return not (self.period == 1 and self.constants[0] == 0)

    def _neg_(self):
        return self.__class__(self.parent(), [-c for c in self.constants])

    def _add_(self, other):
        add_period = lcm(self.period, other.period)
        add_constants = [self(k) + other(k) for k in range(add_period)]
        return self.__class__(self.parent(), add_constants)

    def _radd_(self, other):
        value = self.__add__(other)
        if value is not NotImplemented:
            return value

        raise NotImplementedError(self.NotImplementedError_message(other))

    def _sub_(self, other):
        return self.__add__(-other)

    def _rsub_(self, other):
        return (-self).__add__(other)

    def _mul_(self, other):
        mul_period = lcm(self.period, other.period)
        mul_constants = [self(k)*other(k) for k in range(mul_period)]
        return self.__class__(self.parent(), mul_constants)

    def _rmul_(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if not isinstance(other, self.__class__):
            try:
                other = self.__class__(self.parent(), other)
            except:
                raise TypeError(f"unsupported operand parent(s) for +: {other.__class__} and {self.__class__}")
        div_period = lcm(self.period, other.period)
        div_constants = [self(k)/other(k) for k in range(div_period)]
        return self.__class__(self.parent(), div_constants)

    def __rtruediv__(self, other):
        rdiv_period = lcm(self.period, other.period)
        rdiv_constants = [other(k)/self(k) for k in range(rdiv_period)]
        return self.__class__(self.parent(), rdiv_constants)


# needed for automatic coercion
IntegerPeriodicFunctionRing.Element = IntegerPeriodicFunctionElement

###     Functor     ###
class IntegerPeriodicFunctionFunctor(ConstructionFunctor):
    rank = 10
    def __init__(self, args=None, kwds=None):
        self.args = args or ()
        self.kwds = kwds or {}
        ConstructionFunctor.__init__(self, Rings(), Rings())
    def _apply_functor(self, R):
        return IntegerPeriodicFunctionRing(R,*self.args,**self.kwds)
    def merge(self, other):
        if isinstance(other, type(self)):
            return self


if __name__ == "__main__":
    ipf = IntegerPeriodicFunctionRing(QQ)
    TestSuite(ipf).run(verbose=True, catch = False)