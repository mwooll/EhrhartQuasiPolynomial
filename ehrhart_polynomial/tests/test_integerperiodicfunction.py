from integerperiodicfunction import IntegerPeriodicFunction

from unittest import TestCase, main

class TestIntegerPeriodicFunction(TestCase):
    def setUp(self):
        self.ipf = IntegerPeriodicFunction([1, 2, 3])
        self.zero = IntegerPeriodicFunction([0])

    # __init__
    def test_default_constructor_parameters(self):
        self.assertEqual(IntegerPeriodicFunction(), IntegerPeriodicFunction([0]))

    def test_calculate_period(self):
        self.assertEqual(self.ipf.period, 3)
        self.assertEqual(IntegerPeriodicFunction([4, 4, 4, 4]).period, 1)
        self.assertEqual(IntegerPeriodicFunction([1, 2, 1, 2, 1, 2]).period, 2)

    # call
    def test_call(self):
        self.assertEqual(self.zero(1), 0)
        self.assertEqual(self.ipf(0), 1)
        self.assertEqual(self.ipf(3), 1)
        self.assertEqual(self.ipf(7), 2)

    # standard dunder methods
    def test_repr(self):
        self.assertEqual(eval(repr(self.ipf)), self.ipf)

    def test_str(self):
        self.assertEqual(str(self.ipf),
                         "IntegerPeriodicFunction given by"
                         + "\n\t1 if k%3 == 0"
                         + "\n\t2 if k%3 == 1"
                         + "\n\t3 if k%3 == 2")

    # math support
    def test_neg(self):
        self.assertEqual(-self.ipf, IntegerPeriodicFunction([-1, -2, -3]))

    def test_add(self):
        self.assertEqual(self.ipf + 1, IntegerPeriodicFunction([2, 3, 4]))
        self.assertEqual(1/2 + self.ipf, IntegerPeriodicFunction([1.5, 2.5, 3.5]))

        self.assertEqual(self.ipf + self.zero, self.ipf)
        self.assertEqual(self.ipf + self.ipf, IntegerPeriodicFunction([2, 4, 6]))

        self.assertEqual(self.ipf + IntegerPeriodicFunction([0, 1]),
                         IntegerPeriodicFunction([1, 3, 3, 2, 2, 4]))

    def test_sub(self):
        self.assertEqual(self.ipf - 1, IntegerPeriodicFunction([0, 1, 2]))
        self.assertEqual(1/2 - self.ipf, IntegerPeriodicFunction([-0.5, -1.5, -2.5]))

        self.assertEqual(self.ipf - self.zero, self.ipf)
        # self.assertEqual(self.ipf - self.ipf, self.zero)

        self.assertEqual(self.ipf - IntegerPeriodicFunction([0, 1]),
                         IntegerPeriodicFunction([1, 1, 3, 0, 2, 2]))

    def test_mul(self):
        self.assertEqual(self.ipf*1, self.ipf)
        # self.assertEqual(0*self.ipf, self.zero)
        self.assertEqual(self.ipf*2, IntegerPeriodicFunction([2, 4, 6]))
        self.assertEqual(-4/3*self.ipf, IntegerPeriodicFunction([-4/3, -8/3, -4]))

        self.assertEqual(self.ipf*self.ipf, IntegerPeriodicFunction([1, 4, 9]))
        self.assertEqual(self.ipf*IntegerPeriodicFunction([1, -1]),
                         IntegerPeriodicFunction([1, -2, 3, -1, 2, -3]))

    def test_truediv(self):
        self.assertEqual(self.ipf/2, IntegerPeriodicFunction([0.5, 1, 1.5]))
        self.assertEqual(self.ipf/(-4/5), IntegerPeriodicFunction([-5/4, -5/2, -15/4]))

if __name__ == "__main__":
    main()   