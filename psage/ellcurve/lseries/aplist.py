#################################################################################
#
# (c) Copyright 2011 William Stein
#
#  This file is part of PSAGE.
#
#  PSAGE is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
# 
#  PSAGE is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################

"""
Some general and probably very, very slow code for computing L-series
of elliptic curves over general number fields.

NOTE: This code could probably get moved into Sage at some point,
though it's not at all clear the API below is the best.  For example,
instead of storing keys as ideals, I store them as pairs (Norm,
reduced gens), since (1) I was running into bugs where equal ideals
have different hashes and (2) this format is much easier to look at.
"""


from sage.all import ZZ

def ap(E, p):
    """
    INPUT:
        - `E` -- an elliptic curve over a number field, assumed in
          minimal Weierstrass form
        - `p` -- a prime of the ring of integers of base field of the
          curve
    OUTPUT:
        - the number `a_p(E)`.

    NOTE: This should be really slow.
    """
    k = p.residue_field()
    t = E.local_data(p)
    if t.has_good_reduction():
        Ebar = E.change_ring(k)
        q = k.cardinality()
        return ZZ(q + 1 - Ebar.cardinality())
    elif t.has_split_multiplicative_reduction():
        return ZZ(1)
    elif t.has_nonsplit_multiplicative_reduction():
        return ZZ(-1)
    else:
        return ZZ(0)

def primes_of_bounded_norm(F, B):
    r"""
    Returns an iterator over all prime ideals of norm `\leq B` in the ring
    of integers of the number field F, with the primes ordered by
    their norm.

    INPUT:
        - `F` -- a number field
        - `B` -- a positive integer
    OUTPUT:
        iterator

    NOTE: This could be really slow, since it iterates over all
    ideals, and takes only those that are prime.
    """
    v = F.ideals_of_bdd_norm(B)
    for nm in sorted(v.keys()):
        X = v[nm]
        if len(X)>0 and ZZ(nm).is_prime_power():
            if X[0].is_prime():
                for p in X: yield p

def ap_list(E, B, primes=False):
    r"""
    The Dirichlet coefficients `a_p` of the L-series attached to this
    elliptic curve, for all primes `p` with `N(p) \leq B`, where
    `N(p)` is the norm of `p`.

    INPUT:
        - `E` -- an elliptic curve over a number field, assumed in
          minimal Weierstrass form
        - `B` -- integer
        - ``primes`` -- bool (default: False); if True, also return
          corresponding list of primes up to norm n.

    OUTPUT: list of integers

    NOTE: This should be really slow.
    """
    P = list(primes_of_bounded_norm(E.base_field(), B))
    v = [ap(E,p) for p in P]
    if primes:
        return v, P
    return v

def ap_dict(E, B):
    r"""
    The Dirichlet coefficients `a_p` of the L-series attached to this
    elliptic curve, for all primes `p` with `N(p) \leq n`, where
    `N(p)` is the norm of `p`.

    INPUT:
        - `E` -- an elliptic curve over a number field, assumed in
          minimal Weierstrass form
        - `n` -- integer

    OUTPUT: dictionary mapping reduced rep of primes (N(p),p.reduced_gens()) to integers

    NOTE: This should be really slow.
    """
    P = list(primes_of_bounded_norm(E.base_field(), B))
    return dict([(reduced_rep(p),ap(E,p)) for p in P])

def reduced_rep(I):
    return (I.norm(), I.gens_reduced())

def an_dict_from_ap(ap, N, B):
    r"""
    Give a dict ``ap`` of the `a_p`, for primes with norm up to `B`,
    return the dictionary giving all `a_I` for all ideals `I` up to
    norm `B`.

    NOTE: This code is specific to Dirichlet series of elliptic
    curves.

    INPUT:
        - ``ap`` -- dictionary of ap, as output, e.g., by the ap_dict function
        - `N` -- ideal; conductor of the elliptic curve
        - `B` -- positive integer
    
    OUTPUT: dictionary mapping reduced rep of primes (N(p),p.reduced_gens()) of ideals to integers

    NOTE: This should be really, really slow.  It's really just a toy
    reference implementation.
    """
    from sage.all import prod   # used below
 
    F = N.number_field()
    A = F.ideals_of_bdd_norm(B)

    an = dict(ap)
    
    for n in sorted(A.keys()):
        X = A[n]
        for I in X:
            if reduced_rep(I) in an:
                # prime case, already done
                pass
            else:
                # composite case
                fac = I.factor()
                if len(fac) == 0:
                    # unit ideal
                    an[reduced_rep(I)] = ZZ(1)
                elif len(fac) > 1:
                    # not a prime power, so just multiply together
                    # already known Dirichlet coefficients, for
                    # prime power divisors (which are all known).
                    an[reduced_rep(I)] = prod(an[reduced_rep(p**e)] for p, e in fac)
                else:
                    p, e = fac[0]
                    # a prime power
                    if p.divides(N):
                        # prime divides level
                        an[reduced_rep(I)] = an[reduced_rep(p)]**e
                    else:
                        # prime doesn't divide conductor: a_{p^e} = a_p*a_{p^(e-1)} - Norm(p)*a_{p^(e-2)}
                        assert e >= 2
                        an[reduced_rep(I)] = (an[reduced_rep(p)] * an[reduced_rep(p**(e-1))]
                                                - p.norm()*an[reduced_rep(p**(e-2))])
                        
    return an
                    
            


def an_dict(E, B):
    r"""
    Give an elliptic curve `E` over a number field, return dictionary
    giving the Dirichlet coefficient `a_I` for ideals of norm up to `B`.

    INPUT:
        - ``ap`` -- dictionary of ap, as output, e.g., by the ap_dict function
        - `N` -- ideal; conductor of the elliptic curve
        - `B` -- positive integer
    
    OUTPUT: dictionary mapping reduced rep of ideals (N(p),p.reduced_gens()) to integers

    NOTE: This should be really, really slow.  It's really just a toy
    reference implementation.
    """
    return an_dict_from_ap(ap_dict(E, B), E.conductor(), B)
    

def test1(B=50):
    """
    Tests that the functions all run without crashing over a specific number field.
    Does not test that the output is correct.  That should be in
    another test.
    """
    from sage.all import polygen, QQ, NumberField, EllipticCurve
    x = polygen(QQ,'x')
    F = NumberField(x**2 - x - 1,'a'); a = F.gen()
    E = EllipticCurve([1,a+1,a,a,0])
    ap(E,F.ideal(3))
    primes_of_bounded_norm(F,B)
    ap_list(E,B)
    assert len(ap_list(E,B,primes=True)) == 2
    apd = ap_dict(E,B)
    reduced_rep(F.ideal(3))
    assert an_dict(E,B) == an_dict_from_ap(apd, E.conductor(), B)

def _test_an_dict_over_Q(ainvs, B=100):
    """
    Test that the an_dict function works and gives the correct answer
    for an elliptic curve defined over QQ, by computing using the
    generic code in this file, and comparing with the output of Sage's
    anlist function for rational elliptic curves.
    """
    from sage.all import polygen, QQ, NumberField, EllipticCurve
    x = polygen(QQ,'x')
    F = NumberField(x - 1,'a'); a = F.gen()
    E = EllipticCurve(F, ainvs)
    EQ = EllipticCurve(QQ, ainvs)
    v = EQ.anlist(B)
    an = an_dict(E, B)
    for i, j in an.items():
        assert j == v[i[0]]

def test_an_dict_over_Q():
    """
    Fully test correctness of an_dict for a few curves over QQ.
    """
    _test_an_dict_over_Q([1,2,3,4,5], 50)
    _test_an_dict_over_Q([0,1], 100)  # j = 0
    _test_an_dict_over_Q([4,0], 100)  # j = 1728
    
