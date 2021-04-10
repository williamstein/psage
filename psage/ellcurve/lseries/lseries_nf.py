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

"""
Computing L-series of elliptic curves over general number fields.

EXAMPLES::

    sage: from psage.ellcurve.lseries.lseries_nf import lseries_dokchitser
    sage: K.<a> = NumberField(x^2-x-1); E = EllipticCurve([0,-a,a,0,0])
    sage: L = lseries_dokchitser(E,32); L
    Dokchitser L-function of Elliptic Curve defined by y^2 + a*y = x^3 + (-a)*x^2 over Number Field in a with defining polynomial x^2 - x - 1
    sage: L(1)
    0.422214159
    sage: L.taylor_series(1, 5)
    0.422214159 + 0.575883865*z - 0.102163427*z^2 - 0.158119743*z^3 + 0.120350688*z^4 + O(z^5)

The sign of the functional equation is numerically determined when constructing the L-series::

    sage: L.eps
    1
    

AUTHORS:
    - William Stein
    - Adam Sorkin (very early version: http://trac.sagemath.org/sage_trac/ticket/9402)

"""
from __future__ import absolute_import
from __future__ import division

from builtins import range
from past.utils import old_div
import math

from sage.all import (PowerSeriesRing, Integer, factor, QQ, ZZ,
                      RDF, RealField, Dokchitser, prime_range, prod, pari)
from sage.rings.all import is_NumberField

from psage.ellcurve.lseries.helper import extend_multiplicatively_generic

##################################################################################
# Optimized Special Case: Q(sqrt(5))
##################################################################################
def anlist_over_sqrt5(E, bound):
    """
    Compute the Dirichlet L-series coefficients, up to and including
    a_bound.  The i-th element of the return list is a[i].

    INPUT:
        - E -- elliptic curve over Q(sqrt(5)), which must have
          defining polynomial `x^2-x-1`.
        - ``bound`` -- integer

    OUTPUT:
        - list of integers of length bound + 1

    EXAMPLES::

        sage: from psage.ellcurve.lseries.lseries_nf import anlist_over_sqrt5
        sage: K.<a> = NumberField(x^2-x-1); E = EllipticCurve([0,-a,a,0,0])
        sage: v = anlist_over_sqrt5(E, 50); v
        [0, 1, 0, 0, -2, -1, 0, 0, 0, -4, 0, 3, 0, 0, 0, 0, 0, 0, 0, 5, 2, 0, 0, 0, 0, -4, 0, 0, 0, 11, 0, -6, 0, 0, 0, 0, 8, 0, 0, 0, 0, -1, 0, 0, -6, 4, 0, 0, 0, -6, 0]
        sage: len(v)
        51

    This function isn't super fast, but at least it will work in a few
    seconds up to `10^4`::

        sage: t = cputime()
        sage: v = anlist_over_sqrt5(E, 10^4)
        sage: assert cputime(t) < 5
    """
    from . import aplist_sqrt5
    from psage.number_fields.sqrt5.prime import primes_of_bounded_norm, Prime

    # Compute all of the prime ideals of the ring of integers up to the given bound
    primes = primes_of_bounded_norm(bound+1)

    # Compute the traces of Frobenius: this is supposed to be the hard part
    v      = aplist_sqrt5.aplist(E, bound+1)

    # Compute information about the primes of bad reduction, in
    # particular the integers i such that primes[i] is a prime of bad
    # reduction.
    bad_primes = set([Prime(a.prime()) for a in E.local_data()])


    # We compute the local factors of the L-series as power series in ZZ[T].
    P = PowerSeriesRing(ZZ, 'T')
    T = P.gen()
    # Table of powers of T, so we don't have to compute T^4 (say) thousands of times.
    Tp = [T**i for i in range(5)]

    # For each prime, we write down the local factor.
    L_P = []
    for i, P in enumerate(primes):
        inertial_deg = 2 if P.is_inert() else 1
        a_p = v[i]
        if P in bad_primes:
            # bad reduction
            f = 1 - a_p*Tp[inertial_deg]
        else:
            # good reduction
            q = P.norm()
            f = 1 - a_p*Tp[inertial_deg] + q*Tp[2*inertial_deg]
        L_P.append(f)

    # Use the local factors of the L-series to compute the Dirichlet
    # series coefficients of prime-power index.
    coefficients = [0,1] + [0]*(bound-1)
    i = 0
    while i < len(primes):
        P = primes[i]
        if P.is_split():
            s = L_P[i] * L_P[i+1]
            i += 2
        else:
            s = L_P[i]
            i += 1
        p = P.p
        # We need enough terms t so that p^t > bound
        accuracy_p = int(math.floor(old_div(math.log(bound),math.log(p)))) + 1
        series_p = s.add_bigoh(accuracy_p)**(-1)
        for j in range(1, accuracy_p):
            coefficients[p**j] = series_p[j]

    # Using multiplicativity, fill in the non-prime power Dirichlet
    # series coefficients.
    extend_multiplicatively_generic(coefficients)
    return coefficients


##################################################################################
# General case over number fields.  Largely untouched from trac
# ticket.  Once Q(sqrt(5)) case above is more optimized, revisit the
# three *_over_nf functions below and refactor and improve.
##################################################################################

def get_factor_over_nf(curve, prime_ideal, prime_number, conductor, accuracy):
    """                                                                         
    Returns the inverse of the factor corresponding to the given prime
    ideal in the Euler product expansion of the L-function at
    prime_ideal. Unless the accuracy doesn't need this expansion, and
    then returns 1 in power series ring.
    """
    P = PowerSeriesRing(ZZ, 'T')
    T = P.gen()
    q = prime_ideal.norm()
    inertial_deg = Integer(q).ord(prime_number)
    if inertial_deg > accuracy: 
        return P(1)
    if prime_ideal.divides(conductor):
        a = curve.local_data(prime_ideal).bad_reduction_type()
        L = 1 - a*(T**inertial_deg)
    else:
        discriminant = curve.discriminant()
        if prime_ideal.divides(discriminant):
            a = q + 1 - curve.local_minimal_model(prime_ideal).reduction(prime_ideal).count_points()
        else:
            a = q + 1 - curve.reduction(prime_ideal).count_points()
        L = 1 - a*(T**inertial_deg) + q*(T**(2*inertial_deg))
    return L

def get_coeffs_p_over_nf(curve, prime_number, accuracy=20 , conductor=None):
    """                                                                         
    Computes the inverse of product of L_prime on all primes above prime_number,
    then returns power series of L_p up to desired accuracy.                    
    But will not return power series if need not do so (depends on accuracy).   
    """
    if conductor is None:
        conductor = curve.conductor()
    primes = curve.base_field().prime_factors(prime_number)
    series_p = [get_factor_over_nf(curve, prime_id, prime_number, conductor, accuracy) for prime_id in primes]
    return (prod(series_p).O(accuracy))**(-1)

def anlist_over_nf(E, bound):
    """
    Caution: This method is slow, especially for curves of high
    conductor, or defined over number fields of high degree. 
    The method is to take the Euler product form, and retrieve 
    the coefficients by expanding each product factor as a power 
    series. The bottleneck is counting points over good reductions.

    TODO: Cache this method: it is computed when initializing the 
    class dokchitser, if cached would have .num_coeffs() of a_i stored.

    EXAMPLE::

        sage: K.<i> = NumberField(x^2+1) 
        sage: E = EllipticCurve(K,[0,-1,1,0,0])
        sage: from psage.ellcurve.lseries.lseries_nf import anlist_over_nf
        sage: anlist_over_nf(E, 20)
        [0, 1, -2, 0, 2, 2, 0, 0, 0, -5, -4, 0, 0, 8, 0, 0, -4, -4, 10, 0, 4]
    """    
    conductor = E.conductor()
    coefficients = [0,1] + [0]*(bound-1)
    for p in prime_range(bound+1):
        accuracy_p = int(math.floor(old_div(math.log(bound),math.log(p)))) + 1
        series_p = get_coeffs_p_over_nf(E, p, accuracy_p, conductor)
        for i in range(1, accuracy_p):
            coefficients[p**i] = series_p[i]
    extend_multiplicatively_generic(coefficients)
    return coefficients


########################################################################

def anlist(E, bound):
    r"""                                                                         
    Compute the Dirichlet L-series coefficients, up to and including
    a_bound.  The i-th element of the return list is a[i].

    INPUT:
        - E -- elliptic curve over any number field or the rational numbers
        - ``bound`` -- integer

    OUTPUT:
        - list of integers of length bound + 1

    EXAMPLES::

        sage: from psage.ellcurve.lseries.lseries_nf import anlist
        sage: anlist(EllipticCurve([1,2,3,4,5]),40)
        [0, 1, 1, 0, -1, -3, 0, -1, -3, -3, -3, -1, 0, 1, -1, 0, -1, 5, -3, 4, 3, 0, -1, -6, 0, 4, 1, 0, 1, -2, 0, 2, 5, 0, 5, 3, 3, 7, 4, 0, 9]
        sage: K.<a> = NumberField(x^2-x-1)
        sage: anlist(EllipticCurve(K,[1,2,3,4,5]),40)
        [0, 1, 0, 0, -3, -3, 0, 0, 0, -6, 0, -2, 0, 0, 0, 0, 5, 0, 0, 8, 9, 0, 0, 0, 0, 4, 0, 0, 0, -4, 0, 4, 0, 0, 0, 0, 18, 0, 0, 0, 0]
        sage: K.<i> = NumberField(x^2+1)
        sage: anlist(EllipticCurve(K,[1,2,3,4,5]),40)
        [0, 1, 1, 0, -1, -6, 0, 0, -3, -6, -6, 0, 0, 2, 0, 0, -1, 10, -6, 0, 6, 0, 0, 0, 0, 17, 2, 0, 0, -4, 0, 0, 5, 0, 10, 0, 6, 14, 0, 0, 18]

    Note that the semantics of anlist agree with the anlist method on
    elliptic curves over QQ in Sage::

        sage: from psage.ellcurve.lseries.lseries_nf import anlist
        sage: E = EllipticCurve([1..5])
        sage: v = E.anlist(10); v
        [0, 1, 1, 0, -1, -3, 0, -1, -3, -3, -3]
        sage: len(v)
        11
        sage: anlist(E, 10)
        [0, 1, 1, 0, -1, -3, 0, -1, -3, -3, -3]
    """
    if E.base_field() == QQ:
        # Rational numbers -- use code built into Sage
        v = E.anlist(bound)  
    elif list(E.base_field().defining_polynomial()) == [-1,-1,1]:
        # An optimized special case -- Q(sqrt(5))
        v = anlist_over_sqrt5(E, bound)
    else:
        # General number field -- very slow in general
        v = anlist_over_nf(E, bound)
    return v


def lseries_dokchitser(E, prec=53):
    """
    Return the Dokchitser L-series object associated to the elliptic
    curve E, which may be defined over the rational numbers or a
    number field.  Also prec is the number of bits of precision to
    which evaluation of the L-series occurs.
    
    INPUT:
        - E -- elliptic curve over a number field (or QQ)
        - prec -- integer (default: 53) precision in *bits*

    OUTPUT:
        - Dokchitser L-function object
    
    EXAMPLES::

    A curve over Q(sqrt(5)), for which we have an optimized implementation::
    
        sage: from psage.ellcurve.lseries.lseries_nf import lseries_dokchitser
        sage: K.<a> = NumberField(x^2-x-1); E = EllipticCurve([0,-a,a,0,0])
        sage: L = lseries_dokchitser(E); L
        Dokchitser L-function of Elliptic Curve defined by y^2 + a*y = x^3 + (-a)*x^2 over Number Field in a with defining polynomial x^2 - x - 1
        sage: L(1)
        0.422214159001667
        sage: L.taylor_series(1,5)
        0.422214159001667 + 0.575883864741340*z - 0.102163426876721*z^2 - 0.158119743123727*z^3 + 0.120350687595265*z^4 + O(z^5)

    Higher precision::
    
        sage: L = lseries_dokchitser(E, 200)
        sage: L(1)
        0.42221415900166715092967967717023093014455137669360598558872

    A curve over Q(i)::

        sage: K.<i> = NumberField(x^2 + 1)
        sage: E = EllipticCurve(K, [1,0])
        sage: E.conductor().norm()
        256
        sage: L = lseries_dokchitser(E, 10)
        sage: L.taylor_series(1,5)
        0.86 + 0.58*z - 0.62*z^2 + 0.19*z^3 + 0.18*z^4 + O(z^5)

    More examples::

        sage: lseries_dokchitser(EllipticCurve([0,-1,1,0,0]), 10)(1)
        0.25
        sage: K.<i> = NumberField(x^2+1)
        sage: lseries_dokchitser(EllipticCurve(K, [0,-1,1,0,0]), 10)(1)
        0.37
        sage: K.<a> = NumberField(x^2-x-1)
        sage: lseries_dokchitser(EllipticCurve(K, [0,-1,1,0,0]), 10)(1)
        0.72

        sage: E = EllipticCurve([0,-1,1,0,0])
        sage: E.quadratic_twist(2).rank()
        1
        sage: K.<d> = NumberField(x^2-2)
        sage: L = lseries_dokchitser(EllipticCurve(K, [0,-1,1,0,0]), 10)
        sage: L(1)
        0
        sage: L.taylor_series(1, 5)
        0.58*z + 0.20*z^2 - 0.50*z^3 + 0.28*z^4 + O(z^5)

    You can use this function as an algorithm to compute the sign of the functional
    equation (global root number)::

        sage: E = EllipticCurve([1..5])
        sage: E.root_number()
        -1
        sage: L = lseries_dokchitser(E,32); L
        Dokchitser L-function of Elliptic Curve defined by y^2 + x*y = x^3 - x^2 + 4*x + 3 over Rational Field
        sage: L.eps
        -1

    Over QQ, this isn't so useful (since Sage has a root_number
    method), but over number fields it is very useful::

        sage: K.<a> = NumberField(x^2 - x - 1)
        sage: E1=EllipticCurve([0,-a-1,1,a,0]); E0 = EllipticCurve([0,-a,a,0,0])
        sage: lseries_dokchitser(E1, 16).eps
        -1
        sage: E1.rank()
        1
        sage: lseries_dokchitser(E0, 16).eps
        1
        sage: E0.rank()
        0
    """
    # The code assumes in various places that we have a global minimal model,
    # for example, in anlist_sqrt5 above.
    E = E.global_minimal_model() 

    # Check that we're over a number field.
    K = E.base_field()
    if not is_NumberField(K):
        raise TypeError("base field must be a number field")

    # Compute norm of the conductor -- awkward because QQ elements have no norm method (they should).
    N = E.conductor()
    if K != QQ:
        N = N.norm()

    # We guess the sign epsilon in the functional equation to be +1
    # first.  If our guess is wrong then we just choose the other
    # possibility.
    epsilon = 1

    # Define the Dokchitser L-function object with all parameters set:
    L = Dokchitser(conductor = N * K.discriminant()**2,
                   gammaV = [0]*K.degree() + [1]*K.degree(),
                   weight = 2, eps = epsilon, poles = [], prec = prec)

    # Find out how many coefficients of the Dirichlet series are needed
    # to compute to the requested precision.
    n = L.num_coeffs()
    # print "num coeffs = %s"%n


    # Compute the Dirichlet series coefficients
    coeffs = anlist(E, n)[1:]

    # Define a string that when evaluated in PARI defines a function
    # a(k), which returns the Dirichlet coefficient a_k.
    s = 'v=%s; a(k)=v[k];'%coeffs

    # Actually tell the L-series / PARI about the coefficients.
    L.init_coeffs('a(k)', pari_precode = s)      

    # Test that the functional equation is satisfied.  This will very,
    # very, very likely if we chose the sign of the functional
    # equation incorrectly, or made any mistake in computing the
    # Dirichlet series coefficients. 
    tiny = max(1e-8, old_div(1.0,2**(prec-1)))
    if abs(L.check_functional_equation()) > tiny:
        # The test failed, so we try the other choice of functional equation.
        epsilon *= -1
        L.eps = epsilon

        # It is not necessary to recreate L -- just tell PARI the different sign.
        L._gp_eval('sgn = %s'%epsilon)

        # Once again, verify that the functional equation is
        # satisfied.  If it is, then we've got it.  If it isn't, then
        # there is definitely some other subtle bug, probably in computed
        # the Dirichlet series coefficients.  
        if abs(L.check_functional_equation()) > tiny: 
            raise RuntimeError("Functional equation not numerically satisfied for either choice of sign")
        
    L.rename('Dokchitser L-function of %s'%E)
    return L
