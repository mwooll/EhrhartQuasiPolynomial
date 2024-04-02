from math import lcm
from sympy.ntheory import factorint

class IntegerPeriodicFunction:
    def __init__(self, constants=None):
        if constants is None:
            constants = [0]
        self.period = self._calculate_period(constants)
        self.constants = constants[:self.period]

    """
    methods called in __init__
    """
    def _calculate_period(self, constants):
        length = len(constants)

        if length == 1:
            return 1

        period = length
        factors = factorint(length)
        for prime, power in factors.items():
            div = self._get_period_divisor(constants, length, prime, power)
            period //= div
        
        return period

    def _get_period_divisor(self, constants, length, prime, power):
        last_div = 1
        factor = prime
        for k in range(power):
            if constants[:length//factor]*factor != constants:
                return last_div
            
            last_div = factor
            factor *= prime
        return last_div

    """
    use class object as function
    """
    def __call__(self, k):
        return self.constants[k%self.period]

    """
    standard dunder methods
    """
    def __repr__(self):
        return f"IntegerPeriodicFunction({self.constants})"

    def __str__(self):
        function_str = "IntegerPeriodicFunction given by"
        for ind, val in enumerate(self.constants):
            function_str += f"\n\t{val} if k%{self.period} == {ind}"
        return function_str

    """
    comparison operators
    """
    def __eq__(self, other):
        if isinstance(other, IntegerPeriodicFunction):
            return self.constants == other.constants
        if isinstance(other, (int, float)):
            return self.period == 1 and self.constants[0] == other
        
        return False

    """
    math support
    """
    def __neg__(self):
        return IntegerPeriodicFunction([-c for c in self.constants])

    def __add__(self, other):
        if isinstance(other, (int, float)):
            return IntegerPeriodicFunction([c + other for c in self.constants])
        elif isinstance(other, IntegerPeriodicFunction):
            add_period = lcm(self.period, other.period)
            add_constants = [self(k) + other(k) for k in range(add_period)]
            return IntegerPeriodicFunction(add_constants)
        else:
            raise TypeError("other must be instance of IntegerPeriodicFunction, "
                            + "int or float, but is has type {type(other)}")

    __radd__ = __add__

    def __sub__(self, other):
        return self.__add__(-other)

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return IntegerPeriodicFunction([c*other for c in self.constants])
        elif isinstance(other, IntegerPeriodicFunction):
            mul_period = lcm(self.period, other.period)
            mul_constants = [self(k)*other(k) for k in range(mul_period)]
            return IntegerPeriodicFunction(mul_constants)
        else:
            raise TypeError("other must be instance of IntegerPeriodicFunction, "
                            + "int or float, but is has type {type(other)}")

    __rmul__ = __mul__

    def __truediv__(self, other):
        return self.__mul__(1/other)

    def __rtruediv__(self, other):
        raise NotImplementedError
