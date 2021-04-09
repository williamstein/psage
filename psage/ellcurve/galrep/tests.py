from __future__ import absolute_import
from builtins import range
def test_nonsurj(v=list(range(1,50))):
    """
    For each non CM curve of conductor in the list v, compute the
    primes where the representation isn't surjective using both galrep
    and Sage, and make sure the answers agree.
    """
    from sage.all import cremona_curves, Integer
    from .wrapper import GalRep
    G = GalRep()
    for E in cremona_curves(v):
        if E.has_cm(): continue
        a = E.galois_representation().non_surjective()
        F = E.short_weierstrass_model()
        b = G.non_surjective_primes(Integer(F.a4()), Integer(F.a6()))
        if a != b:
            raise RuntimeError("Test failed for {0}!".format(E.cremona_label()))
    
