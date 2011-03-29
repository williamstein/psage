from sage.all import polygen, NumberField, EllipticCurve

from psage.modform.hilbert.sqrt5.sqrt5 import primes_of_bounded_norm, F

a = F.gen()

import psage.ellcurve.lseries.sqrt5 as sqrt5
import psage.modform.hilbert.sqrt5.sqrt5_fast as sqrt5_fast

def test_ap_via_enumeration(B=1000):
    _,_,_,a4,a6 = E.short_weierstrass_model().a_invariants()
    a4 = 6912*a - 5643
    a6 = -131328*a + 298566
    E = EllipticCurve([a4, a6])
    D = E.discriminant()
    for P in primes_of_bounded_norm(F, B):
        if D.valuation(P) == 0:
            R = sqrt5_fast.ResidueRing(P, 1)
            ap0 = sqrt5.ap_via_enumeration(R(a4), R(a6))
            k = P.residue_field(); E0 = E.change_ring(k); ap1 = k.cardinality() + 1 - E0.cardinality()
            assert ap0 == ap1
        

