from math import lcm

from integerperiodicfunction import IntegerPeriodicFunction

class QuasiPolynomial():
    def __init__(self, coefs=None):
        if coefs is None:
            coefs = [0]
        else:
            if not isinstance(coefs, (tuple, list)):
                return TypeError("coefs must be instance of list or tuple, or be None"
                                 + f", but has type {type(coefs)}")

        for index, coef in enumerate(coefs):
            if isinstance(coef, (int, float)):
                coefs[index] = IntegerPeriodicFunction([coef])
            elif isinstance(coef, IntegerPeriodicFunction):
                pass
            else:
                raise ValueError("elements of coefs must be instances of int, "
                                 + " float or IntegerPeriodicFunction but got "
                                 + f"{coef} of type {type(coef)}")

        self.coefs, self.degree = self._reduce_coefs(coefs)
        self.period = self._calculate_peroid()

    """
    methods called in __init__
    """
    def _reduce_coefs(self, coefs):
        degree = len(coefs)
        while coefs[degree-1] == 0 and degree > 1:
            degree -= 1

        if degree == 0:
            coefs = [0]
            degree = 1
        else:
            coefs = coefs[:degree]

        return coefs, degree

    def _calculate_peroid(self):
        periods = (coef.period for coef in self.coefs)
        period = lcm(*periods)
        return period

    """
    use class obejct as function
    """
    def __call__(self, value):
        result = 0
        for power, coef in enumerate(self.coefs):
            result += coef(value) * value**power
        return result

    """
    standard dunder methods
    """
    def __str__(self):
        function_str = "QuasiPolynomial given by \n"
        function_str += f"{repr(self.coefs[0])}(k)"
        for power, coef in enumerate(self.coefs[1:]):
            function_str += f" + {repr(coef)}(k)*k^{power+1}"
        return function_str

    def __repr__(self):
        return f"QuasiPolynomial({self.coefs})"

    """
    comparison operators
    """
    def __eq__(self, other):
        if isinstance(other, QuasiPolynomial):
            return self.coefs == other.coefs
        return False

    """
    math support
    """
    def __neg__(self):
        return QuasiPolynomial([-c for c in self.coefs])

    def __add__(self, other):
        if isinstance(other, QuasiPolynomial):
            add_coefs = []
            if self.degree > other.degree:
                sum_coefs = self.coefs.copy()
                add_coefs = other.coefs
            else:
                sum_coefs = other.coefs.copy()
                add_coefs = self.coefs
    
            for index, coef in enumerate(add_coefs):
                sum_coefs[index] += coef
    
            return QuasiPolynomial(sum_coefs)

        if isinstance(other, (int, float)):
            new_coefs = self.coefs.copy()
            new_coefs[0] += other
            return QuasiPolynomial(new_coefs)

        else:
            return TypeError("other must be instance of QuasiPolynomial, int"
                             f" or float, but has type {type(other)}")

    __radd__ = __add__

    def __sub__(self, other):
        return self.__add__(-other)

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return QuasiPolynomial([other*c for c in self.coefs])
            
        if isinstance(other, QuasiPolynomial):
            mul_coefs = [0]*(self.degree + other.degree + 1)
            for self_power, self_coef in enumerate(self.coefs):
                for  other_power, other_coef in enumerate(other.coefs):
                    mul_coefs[self_power + other_power] += self_coef*other_coef
            return QuasiPolynomial(mul_coefs)
        else:
            return TypeError("other must be instance of int, float or" +
                             f" QuasiPolynomial but has type {type(other)}")
    __rmul__ = __mul__

    def __truediv__(self, other):
        return self.__mul__(1/other)
    
    def __rtruediv__(self, other):
        raise NotImplementedError
