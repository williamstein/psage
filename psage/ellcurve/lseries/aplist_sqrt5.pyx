#################################################################################
#
# (c) Copyright 2011 William Stein
#
#  This file is part of PSAGE
#
#  PSAGE is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
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
Frobenius Traces over Q(sqrt(5))

This module implements functionality for computing traces `a_P` of
Frobenius for an elliptic curve over Q(sqrt(5)) efficiently.

EXAMPLES::

    sage: from psage.ellcurve.lseries.aplist_sqrt5 import aplist
    sage: from psage.number_fields.sqrt5.misc import F, a
    sage: E = EllipticCurve([1,-a,a,a+2,a-5])
    sage: aplist(E, 60)
    [-3, -1, 5, -4, 5, 8, 1, 9, -10, 10, 10, -4, 6, -10, 4, 3]

The `a_P` in the above list exactly correspond to those output by the primes_of_bounded_norm function::

    sage: from psage.number_fields.sqrt5 import primes_of_bounded_norm
    sage: primes_of_bounded_norm(60)
    [2a, 5a, 3a, 11a, 11b, 19a, 19b, 29a, 29b, 31a, 31b, 41a, 41b, 7a, 59a, 59b]
"""

from sage.libs.gmp.mpz cimport (mpz_t, mpz_set, mpz_set_si, mpz_init, mpz_clear, mpz_fdiv_ui)

from sage.rings.integer cimport Integer

from psage.libs.smalljac.wrapper1 import elliptic_curve_ap

from psage.number_fields.sqrt5.prime import primes_of_bounded_norm

from psage.modform.hilbert.sqrt5.sqrt5_fast cimport ResidueRing_abstract, residue_element
from psage.modform.hilbert.sqrt5.sqrt5_fast import ResidueRing

def short_weierstrass_invariants(E):
    """
    Compute the invariants of a short Weierstrass form of E.

    The main motivation for this function is that it doesn't require
    constructing an elliptic curve like the short_weierstrass_model
    method on E does.

    INPUT:
        - E -- an elliptic curve

    OUTPUT:
        - two elements A, B of the base field of the curve such that
          E is isomorphic to `y^2 = x^3 + Ax + B`.

    EXAMPLES::

        sage: from psage.ellcurve.lseries.aplist_sqrt5 import short_weierstrass_invariants
        sage: from psage.number_fields.sqrt5.misc import F, a
        sage: E = EllipticCurve([1,-a,a,a+2,a-5])
        sage: short_weierstrass_invariants(E)
        (1728*a + 2133, 101952*a - 206874)
        sage: E.short_weierstrass_model()
        Elliptic Curve defined by y^2 = x^3 + (1728*a+2133)*x + (101952*a-206874) over Number Field in a with defining polynomial x^2 - x - 1
    """
    a1, a2, a3, a4, a6 = E.a_invariants()
    # Compute short Weierstrass form directly (for speed purposes)
    if a1 or a2 or a3:
        b2 = a1*a1 + 4*a2; b4 = a1*a3 + 2*a4; b6 = a3**2 + 4*a6
        if b2:
            c4 = b2**2 - 24*b4; c6 = -b2**3 + 36*b2*b4 - 216*b6
            A = -27*c4; B = -54*c6
        else:
            A = 8*b4; B = 16*b6
    else:
        A = a4; B = a6

    return A, B

def aplist(E, bound, fast_only=False):
    """
    Compute the traces of Frobenius `a_P` of the elliptic curve E over Q(sqrt(5)) for
    all primes `P` with norm less than the given bound.

    INPUT:
        - `E` -- an elliptic curve given by an integral (but not
          necessarily minimal) model over the number field with
          defining polynomial `x^2-x-1`.
        - ``bound`` -- a nonnegative integer
        - ``fast_only`` -- (bool, default: False) -- if True, only the
          `a_P` that can be computed efficiently are computed, and the
          rest are left as None

    EXAMPLES::

        sage: from psage.ellcurve.lseries.aplist_sqrt5 import aplist
        sage: K.<a> = NumberField(x^2-x-1); E = EllipticCurve([1,-a,a,a-1,a+3])
        sage: aplist(E,60)
        [1, -2, 2, -4, -2, 0, -5, 0, -5, 0, 2, 11, 12, 10, -11, -6]
        sage: from psage.number_fields.sqrt5 import primes_of_bounded_norm
        sage: primes_of_bounded_norm(60)
        [2a, 5a, 3a, 11a, 11b, 19a, 19b, 29a, 29b, 31a, 31b, 41a, 41b, 7a, 59a, 59b]

    If you give the fast_only option, then only the `a_P` for which
    the current implemenation can compute quickly are computed::
    
        sage: aplist(E,60,fast_only=True)
        [None, None, None, -4, -2, 0, -5, 0, -5, 0, 2, 11, 12, 10, -11, -6]

    This is fast enough that up to a million should only take a few
    seconds::
    
        sage: t = cputime(); v = aplist(E,10^6,fast_only=True)
        sage: assert cputime(t) < 15, "too slow!"

    TESTS::

    We test some edge cases::

        sage: aplist(E, 0)
        []
        sage: L.<b> = NumberField(x^2-5); F = EllipticCurve([1,-b,b,0,0])
        sage: aplist(F, 50)
        Traceback (most recent call last):
        ...
        ValueError: E must have base field with defining polynomial x^2-x-1
        sage: E = EllipticCurve([1,-a,a,a+2,a/7])
        sage: aplist(E, 50)
        Traceback (most recent call last):
        ...
        ValueError: coefficients of the input curve must be algebraic integers

        sage: K.<a> = NumberField(x^2 - x - 1)
        sage: E = EllipticCurve(K, [1, 1, 1, 0, 0])
        sage: aplist(E, 50)
        [-3, 1, 1, -4, -4, 4, 4, -2, -2, 0, 0, 10, 10, -14]
    """
    if list(E.base_field().defining_polynomial()) != [-1,-1,1]:
        raise ValueError, "E must have base field with defining polynomial x^2-x-1"
    A, B = short_weierstrass_invariants(E)
    primes = primes_of_bounded_norm(bound)
    cdef list v = aplist_short(A, B, primes)
    if not fast_only:
        aplist_remaining_slow(E, v, primes)
    return v

def aplist_remaining_slow(E, v, primes):
    """
    Compute -- using possibly very slow methods -- the `a_P` in the
    list v that are currently set to None.

    Here v and primes should be two lists, where primes is a list of
    Cython Prime objects, as output by the primes_of_bounded_norm
    function.  This function currently does not assume anything about
    the model of E being minimal or even integral.

    INPUT:
        - `E` -- elliptic curve over Q(sqrt(5))
        - `v` -- list of ints or None
        - `primes` -- list of Prime objects
        
    OUTPUT:
        - modified v in place by changing all of the None's to ints

    EXAMPLES::

        sage: from psage.ellcurve.lseries.aplist_sqrt5 import aplist
        sage: K.<a> = NumberField(x^2-x-1); E = EllipticCurve([1,-a,a,a-1,a+3])
        sage: v = aplist(E, 60, fast_only=True); v
        [None, None, None, -4, -2, 0, -5, 0, -5, 0, 2, 11, 12, 10, -11, -6]
        sage: from psage.number_fields.sqrt5 import primes_of_bounded_norm
        sage: primes = primes_of_bounded_norm(60)
        sage: from psage.ellcurve.lseries.aplist_sqrt5 import aplist_remaining_slow
        sage: aplist_remaining_slow(E, v, primes); v
        [1, -2, 2, -4, -2, 0, -5, 0, -5, 0, 2, 11, 12, 10, -11, -6]

    We can also use this (slow) function to compute all entries of v.  This must give
    the same answer as above, of course::
    
        sage: v = [None]*len(primes)
        sage: aplist_remaining_slow(E, v, primes)
        sage: v
        [1, -2, 2, -4, -2, 0, -5, 0, -5, 0, 2, 11, 12, 10, -11, -6]
    """
    if len(v) != len(primes):
        raise ValueError, "input lists v and primes must have the same length"
    cdef Py_ssize_t i
    F = E.global_minimal_model()
    N = F.conductor()
    for i in range(len(v)):
        if v[i] is None:
            P = primes[i].sage_ideal()
            m = N.valuation(P)
            if m >= 2:
                v[i] = 0
            elif m == 1:
                if F.has_split_multiplicative_reduction(P):
                    v[i] = 1
                else:
                    v[i] = -1
            else:
                k = P.residue_field()
                v[i] = k.cardinality() + 1 - F.change_ring(k).cardinality()

def aplist_short(A, B, primes):
    """
    Return list of `a_P` values that can be quickly computed for the
    elliptic curve `y^2=x^3+AX+B` and the given list of primes.  See
    the important warnings in the INPUT section below.
    
    INPUT:
        - A, B -- algebraic integers in Q(sqrt(5)), which we assume
          (without checking!) is defined by the polynomial `x^2-x-1`.
          If you violate this assumption, you will get nonsense.
        - ``primes`` -- list of Prime objects
    OUTPUT:
        - list of ints or None

    EXAMPLES::

        sage: from psage.ellcurve.lseries.aplist_sqrt5 import aplist_short
        sage: from psage.number_fields.sqrt5 import primes_of_bounded_norm
        sage: K.<a> = NumberField(x^2-x-1)
        sage: aplist_short(2*a+5, 3*a-7, primes_of_bounded_norm(100))
        [None, None, None, 2, None, 2, 0, -6, -6, -4, -8, 3, 1, 10, -10, 0, 0, -4, -3, -16, 14, -8, 8, 15]
    """
    # Store the short Weierstrass form coefficients as a 4-tuple of GMP ints, 
    # where (Ax,Ay) <--> Ax + Ay * (1+sqrt(5))/2.
    cdef mpz_t Ax, Ay, Bx, By
    mpz_init(Ax); mpz_init(Ay); mpz_init(Bx); mpz_init(By)
    v = A._coefficients()  # slow
    cdef Integer z
    cdef int N
    try:
        N = len(v)
        if N == 0:
            mpz_set_si(Ax, 0)
            mpz_set_si(Ay, 0)
        elif N == 1:
            z = Integer(v[0]); mpz_set(Ax, z.value)
            mpz_set_si(Ay, 0)
        else:
            z = Integer(v[0]); mpz_set(Ax, z.value)
            z = Integer(v[1]); mpz_set(Ay, z.value)
        v = B._coefficients()
        N = len(v)
        if N == 0:
            mpz_set_si(Bx, 0)
            mpz_set_si(By, 0)
        elif N == 1:
            z = Integer(v[0]); mpz_set(Bx, z.value)
            mpz_set_si(By, 0)
        else:
            z = Integer(v[0]); mpz_set(Bx, z.value)
            z = Integer(v[1]); mpz_set(By, z.value)
    except TypeError:
        raise ValueError, "coefficients of the input curve must be algebraic integers"

    v = compute_ap(Ax, Ay, Bx, By, primes)

    mpz_clear(Ax); mpz_clear(Ay); mpz_clear(Bx); mpz_clear(By)
    return v

from psage.number_fields.sqrt5.prime cimport Prime

cdef object compute_split_trace(mpz_t Ax, mpz_t Ay, mpz_t Bx, mpz_t By, long p, long r):
    """
    Return the trace of Frobenius for the elliptic curve `E` given by
    `y^2 = x^3 + Ax + B` at the split prime `P=(p, a-r)`.  Here `a` is
    a root of `x^2-x-1`, and `A=Ax+a*Ay (mod P)`, and `B=Bx+a*By (mod
    P)`.  If the model for E has bad reduction, instead return None.

    This is used internally by the aplist and aplist_short function. It
    has no external Python API.

    INPUT:
        - Ax, Ay, Bx, By -- MPIR integers
        - p, r -- C long

    OUTPUT:
        - Python object -- either an int or None
    """
    cdef long a, b
    
    # 1. Reduce A, B modulo the prime to get a curve over GF(p)
    a = mpz_fdiv_ui(Ax, p) + mpz_fdiv_ui(Ay, p)*r
    b = mpz_fdiv_ui(Bx, p) + mpz_fdiv_ui(By, p)*r

    # 2. Check that the resulting curve is nonsingular
    if (-64*((a*a) % p) *a - 432*b*b)%p == 0:
        # curve is singular at p -- leave this to some other algorithm
        return None

    # 3. Compute the trace using smalljac
    return elliptic_curve_ap(a, b, p)

cdef object compute_inert_trace(mpz_t Ax, mpz_t Ay, mpz_t Bx, mpz_t By, long p):
    """
    Return the trace of Frobenius for the elliptic curve `E` given by
    `y^2 = x^3 + Ax + B` at the inert prime `P=(p)`.  Here `a` is
    a root of `x^2-x-1`, and `A=Ax+a*Ay (mod P)`, and `B=Bx+a*By (mod
    P)`.  If the model for E has bad reduction, instead return None.

    This is used internally by the aplist and aplist_short function. It
    has no external Python API.

    INPUT:
        - Ax, Ay, Bx, By -- MPIR integers
        - p -- C long

    OUTPUT:
        - Python object -- either an int or None
    """
    if p <= 3 or p >= 170:
        # no point using this naive enumeration in these cases
        return None
    
    from psage.number_fields.sqrt5.misc import F
    
    cdef ResidueRing_abstract R = ResidueRing(F.ideal([p]), 1)
    cdef residue_element a4, a6, x, z, w
    cdef Py_ssize_t i
    
    a4[0] = mpz_fdiv_ui(Ax, p); a4[1] = mpz_fdiv_ui(Ay, p)
    a6[0] = mpz_fdiv_ui(Bx, p); a6[1] = mpz_fdiv_ui(By, p)

    # Check if curve is singular, i.e., that this is nonzero: 64*a4^3 + 432*a6^2
    R.pow(z, a4, 3)
    R.mul(w, a6, a6)
    R.coerce_from_nf(x, F(64))
    R.mul(z, z, x)
    R.coerce_from_nf(x, F(432))
    R.mul(w, w, x)
    R.add(z, z, w)
    if R.element_is_0(z):
        # singular curve
        return None

    cdef long cnt = 1  # point at infinity
    for i in range(R.cardinality()):
        R.unsafe_ith_element(x, i)
        R.mul(z, x, x)   # z = x*x
        R.mul(z, z, x)   # z = x^3
        R.mul(w, a4, x)  # w = a4*x
        R.add(z, z, w)   # z = z + w = x^3 + a4*x
        R.add(z, z, a6)  # z = x^3 + a4*x + a6
        if R.element_is_0(z):
            cnt += 1
        elif R.is_square(z):
            cnt += 2
    ap = R.cardinality() + 1 - cnt
    return ap

    
cdef compute_ap(mpz_t Ax, mpz_t Ay, mpz_t Bx, mpz_t By, primes):
    """
    Return list of `a_P` values that can be quickly computed for the
    elliptic curve `y^2=x^3+AX+B` and the given list of primes.  Here
    A and B are defined by Ax,Ay,Bx,By exactly as in the function
    compute_split_traces in this file (so see its documentation).

    INPUT:
        - Ax, Ay, Bx, By -- MPIR integers
        - primes -- list of Prime objects

    OUTPUT:
        - list whose entries are either an int or None
    """
    cdef Prime P
    cdef list v = []
    for P in primes:
        if P.is_split():
            ap = compute_split_trace(Ax, Ay, Bx, By, P.p, P.r)
        elif P.is_ramified():
            ap = None
        else: # inert
            ap = compute_inert_trace(Ax, Ay, Bx, By, P.p)
        v.append(ap)
    return v
            
