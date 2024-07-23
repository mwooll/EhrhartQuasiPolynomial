from math import lcm

from .integerperiodicfunction import IntegerPeriodicFunction

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

# commutative algebra fromework from adm-cycles?

class QuasiPolynomial():
    def __init__(self, coefficients=None):
        if coefficients is None:
            coefficients = [0]
        else:
            if not isinstance(coefficients, (tuple, list)):
                return TypeError("coefficients must be instance of list or tuple, or be None"
                                 + f", but has type {type(coefficients)}")

        for index, coef in enumerate(coefficients):
            if isinstance(coef, (int, float)):
                coefficients[index] = IntegerPeriodicFunction([coef])

            elif isinstance(coef, IntegerPeriodicFunction):
                pass

            elif hasattr(coef, "__int__"):
                coefficients[index] = IntegerPeriodicFunction([int(coef)])

            elif hasattr(coef, "__float__"):
                coefficients[index] = IntegerPeriodicFunction([float(coef)])

            else:
                raise TypeError("elements of coefficients must be instances of "
                                + "IntegerPeriodicFunction, int or float,"
                                + " or implement __int__ or __float__"
                                + f"but is has type {type(coef)}")

        self.coefficients, self.degree = self._reduce_coefficients(coefficients)
        self.period = self._calculate_peroid()

    """
    methods called in __init__
    """
    def _reduce_coefficients(self, coefficients):
        degree = len(coefficients)-1
        while degree > 0 and coefficients[degree] == 0:
            degree -= 1

        coefficients = coefficients[:degree+1]
        return coefficients, degree

    def _calculate_peroid(self):
        periods = (coef.period for coef in self.coefficients)
        period = lcm(*periods)
        return period

    """
    use class obejct as function
    """
    def __call__(self, value):
        result = 0
        for power, coef in enumerate(self.coefficients):
            result += coef(value) * value**power
        return result

    def coefficients(self):
        return self.coefficients

    """
    standard dunder methods
    """
    def __str__(self):
        function_str = "QuasiPolynomial given by \n"
        function_str += f"{self.coefficients[0].coefficient_repr('k')}"
        for power, coef in enumerate(self.coefficients[1:]):
            function_str += f" + {coef.coefficient_repr('k')}*k" + f"^{power+1}"*(power>0)
        return function_str

    def __repr__(self):
        return f"QuasiPolynomial({self.coefficients})"

    """
    comparison operators
    """
    def __eq__(self, other):
        if isinstance(other, QuasiPolynomial):
            return self.coefficients == other.coefficients

        if isinstance(other, int):
            return self.degree == 0 and self.coefficients[0] == other
        return False


    """
    math support
    """
    def NotImplementedError_message(self, other):
        return ("other must be instance of QuasiPolynomial, "
                + "IntegerPeriodicFunction, int or float or implement",
                + f"__int__ or __float__ but has type {type(other)}")

    def __neg__(self):
        return QuasiPolynomial([-c for c in self.coefficients])

    def __add__(self, other):
        if isinstance(other, QuasiPolynomial):
            add_coefficients = []
            if self.degree > other.degree:
                sum_coefficients = self.coefficients.copy()
                add_coefficients = other.coefficients
            else:
                sum_coefficients = other.coefficients.copy()
                add_coefficients = self.coefficients

            for index, coef in enumerate(add_coefficients):
                sum_coefficients[index] += coef

            return QuasiPolynomial(sum_coefficients)

        elif isinstance(other, (int, float, IntegerPeriodicFunction)):
            return QuasiPolynomial([c + other*(index==0)
                                    for index, c in enumerate(self.coefficients)])

        elif hasattr(other, "__int__"):
            return QuasiPolynomial([c + int(other)*(index==0)
                                    for index, c in enumerate(self.coefficients)])

        elif hasattr(other, "__float__"):
            return QuasiPolynomial([c + float(other)*(index==0)
                                    for index, c in enumerate(self.coefficients)])

        else:
            return NotImplemented

    def __radd__(self, other):
        value = self.__add__(other)
        if value is not NotImplemented:
            return value

        raise NotImplementedError(self.NotImplementedError_message(other))

    def __sub__(self, other):
        return self.__add__(-other)

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        if isinstance(other, (int, float, IntegerPeriodicFunction)):
            return QuasiPolynomial([other*c for c in self.coefficients])

        elif isinstance(other, QuasiPolynomial):
            mul_coefficients = [0]*(self.degree + other.degree + 1)
            for self_power, self_coef in enumerate(self.coefficients):
                for  other_power, other_coef in enumerate(other.coefficients):
                    mul_coefficients[self_power + other_power] += self_coef*other_coef
            return QuasiPolynomial(mul_coefficients)

        elif hasattr(other, "__int__"):
            return QuasiPolynomial([c*int(other) for c in self.coefficients])

        elif hasattr(other, "__float__"):
            return QuasiPolynomial([c*float(other) for c in self.coefficients])

        else:
            return NotImplemented

    def __rmul__(self, other):
        value = self.__mul__(other)
        if value is not NotImplemented:
            return value

        raise NotImplementedError(self.NotImplementedError_message(other))

    def __truediv__(self, other):
        return self.__mul__(1/other)

    def __rtruediv__(self, other):
        raise NotImplementedError("dividing by a QuasiPolynomial seems illegal")
