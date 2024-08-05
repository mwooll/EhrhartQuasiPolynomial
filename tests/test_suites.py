from ehrhart_polynomial import IntegerPeriodicFunctionRing, QuasiPolynomialRing

from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.symbolic.ring import SR

from sage.misc.sage_unittest import TestSuite

def IPFR_TestSuites():
    print("TestSuites for IPFR:")
    for ring in [IntegerModRing(19), ZZ, QQ, SR]:
        print(f"\n\tTesting IntegerPeriodicFunctionRing with {ring}")
        ipfr = IntegerPeriodicFunctionRing(ring)
        TestSuite(ipfr).run()
    print("\n")

def QPR_TestSuites():
    print("TestSuites for QPR:")
    for ring in [IntegerModRing(19), ZZ, QQ, SR]:
        print(f"\n\tTesting QuasiPolynomialRing with {ring}")
        qpr = QuasiPolynomialRing(ring)
        TestSuite(qpr).run()
    print("\n")

if __name__ == "__main__":
    IPFR_TestSuites()
    QPR_TestSuites()
