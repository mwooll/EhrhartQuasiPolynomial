import numbers

from sage.calculus.var import var
from sage.categories.pushout import ConstructionFunctor
from sage.categories.commutative_rings import CommutativeRings
from sage.rings.ring import CommutativeRing
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import RingElement

import doctest
from sage.misc.sage_unittest import TestSuite


# from .gfan import secondary_fan
# from .integerperiodicfunction import IntegerPeriodicFunctionRing

class PiecewisePolynomialElement(RingElement):
    def __init__(self, parent, constituents=None, extend_value=0):
        r"""
        INPUT:
            - parent : instance of "IntegerPeriodicFunctionRing"
            - constituents : iterable of pairs of domains and polynomials (default=None)
        """
        base = parent.base()
        amb_dim = parent.ambient_dimension()

        if constituents is None:
            constituents = [base, base(0)]
        elif not hasattr(constituents, "__iter__"):
            constituents = [base(constituents)]
        else:
            constituents = [base(c) for c in constituents]
        
        RingElement.__init__(self, parent)

    def __call__(self, value):
        if hasattr(value, "__iter__"):
            return [self(val) for val in value]

        
        return 0

class PiecewisePolynomialRing(UniqueRepresentation, CommutativeRing):
    def __init__(self, base, num_var=1, extend=True):
        if base not in CommutativeRings():
            raise ValueError(f"{base} is not a commutative ring.")
        CommutativeRing.__init__(self, base)

        self._variables = self.create_variables(num_var)
        self._amb_dim = num_var
        self._extend = extend

    def _create_variables(self, num_var):
        if not isinstance(num_var, numbers.Integral):
            raise TypeError("'num_var' must be an instance of 'numbers.Integral', "
                            + "but is of type {type(num_var)}")
        if num_var <= 0:
            raise ValueError(f"Cannot construct PiecewisePolynomialRing with {num_var} variables")
        
        return var(", ".join([f"x{k}" for k in range(num_var)]))

    def _repr_(self):
        return f"Ring of Piecewise Polynomials {self.variables} over {self.base()}"

    def base_ring(self):
        r"""
        Returns the base_ring of the base of self
        """
        return self.base().base_ring()

    def characteristic(self):
        r"""
        Returns the characteristic of the base of self
        """
        return self.base().characteristic()

    def gens(self):
        r"""
        Returns the generators of self
        """
        return self.variables

    def ambient_dimension(self):
        r"""
        Returns the dimension of the ambient space
        """
        return self.amb_dim


class PiecewisePolynomialFunctor(ConstructionFunctor):
    rank = 10
    def __init__(self, args=None, kwds=None):
        self.args = args or ()
        self.kwds = kwds or {}
        ConstructionFunctor.__init__(self, CommutativeRings(), CommutativeRings())

    def _apply_functor(self, R):
        return PiecewisePolynomialRing(R, *self.args, **self.kwds)

    def merge(self, other):
        if isinstance(other, type(self)):
            return self
    
if __name__ == "__main__":
    from sage.rings.rational_field import QQ
    
    pwfr = PiecewisePolynomialRing(QQ)
    # TestSuite().run()
    # doctest.testmod()
    print(pwfr.base_ring())
    print(pwfr.characteristic())