#################################################################################
#
# (c) Copyright 2010 William Stein
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
TODO:

  [ ] make ModularSymbolMap work with Cremona modular symbols too.
  [ ] normalization
  
Can do both of these by changing to:
   - use E.modular_symbol(use_eclib=True) instead of dual rational
     space and multiply (but make both an option, since for abelian
     varieties the current code is more flexible.
   - just figure out the denominator once we know the values on P^1,
   - actually take into account denominator in this file, which we don't do!
"""

include 'stdsage.pxi'
include 'interrupt.pxi'

from sage.rings.all import ZZ, Zp, Qp, Integers, infinity, binomial

from sage.rings.integer cimport Integer
from modular_symbol_map cimport ModularSymbolMap


cdef class pAdicLseries:
    cdef object E
    cdef int p, prec
    cdef long normalization
    cdef long _alpha
    cdef long alpha_inv[64]
    cdef long p_pow[64]
    cdef long *teich
    cdef public ModularSymbolMap modsym

    def __cinit__(self):
        self.teich = NULL

    def __dealloc__(self):
        if self.teich:
            sage_free(self.teich)
    
    def __init__(self, E, p, algorithm='eclib'):
        """
        INPUT:
            - E -- an elliptic curve over QQ
            - p -- a prime of good ordinary reduction with E[p] surjective
            - algorithm -- str (default: 'eclib')
               - 'eclib' -- use elliptic curve modular symbol computed using eclib
               - 'sage' -- use elliptic curve modular symbol computed using sage
               - 'modsym' -- use sage modular symbols directly (!! not correctly normalized !!)

        EXAMPLES::

            sage: from psage.modform.rational.padic_elliptic_lseries_fast import pAdicLseries
            sage: E = EllipticCurve('389a'); L = pAdicLseries(E, 7)
            sage: L.series()
            O(7) + O(7)*T + (5 + O(7))*T^2 + (3 + O(7))*T^3 + (6 + O(7))*T^4 + O(T^5)
            sage: L.series(n=3)
            O(7^2) + O(7^2)*T + (5 + 4*7 + O(7^2))*T^2 + (3 + 6*7 + O(7^2))*T^3 + (6 + 3*7 + O(7^2))*T^4 + O(T^5)
            sage: L.series(n=4)
            O(7^3) + O(7^3)*T + (5 + 4*7 + 5*7^2 + O(7^3))*T^2 + (3 + 6*7 + 4*7^2 + O(7^3))*T^3 + (6 + 3*7 + 3*7^2 + O(7^3))*T^4 + O(T^5)

        We can also get the series just mod 7::

            sage: L.series_modp()
            5*T^2 + 3*T^3 + 6*T^4 + O(T^5)
            sage: L.series_modp().list()
            [0, 0, 5, 3, 6]        

        We use Sage modular symbols instead of eclib's::
        
            sage: L = pAdicLseries(E, 7, algorithm='sage')
            sage: L.series(n=3)
            O(7^2) + O(7^2)*T + (5 + 4*7 + O(7^2))*T^2 + (3 + 6*7 + O(7^2))*T^3 + (6 + 3*7 + O(7^2))*T^4 + O(T^5)

        When we use algorithm='modsym', we see that the normalization
        is not correct (as is documented above -- no attempt is made
        to normalize!)::
        
            sage: L = pAdicLseries(E, 7, algorithm='modsym')
            sage: L.series(n=3)
            O(7^2) + O(7^2)*T + (6 + 5*7 + O(7^2))*T^2 + (5 + 6*7 + O(7^2))*T^3 + (3 + 5*7 + O(7^2))*T^4 + O(T^5)
            sage: 2*L.series(n=3)
            O(7^2) + O(7^2)*T + (5 + 4*7 + O(7^2))*T^2 + (3 + 6*7 + O(7^2))*T^3 + (6 + 3*7 + O(7^2))*T^4 + O(T^5)

        These agree with the L-series computed directly using separate code in Sage::
        
            sage: L = E.padic_lseries(7)
            sage: L.series(3)
            O(7^5) + O(7^2)*T + (5 + 4*7 + O(7^2))*T^2 + (3 + 6*7 + O(7^2))*T^3 + (6 + 3*7 + O(7^2))*T^4 + O(T^5)            
        
        """
        self.E = E
        self.p = p

        # prec = biggest n such that p^n <= 2^63, so n = log_p(2^63)
        self.prec = ZZ(2**63).exact_log(p)
        
        assert Integer(p).is_pseudoprime(), "p (=%s) must be prime"%p
        assert E.is_ordinary(p), "p (=%s) must be ordinary for E"%p

        if algorithm == 'eclib':
            f = E.modular_symbol(sign=1, use_eclib=True)
            self.modsym = ModularSymbolMap(f)
        elif algorithm == 'sage':
            f = E.modular_symbol(sign=1, use_eclib=False)
            self.modsym = ModularSymbolMap(f)
        elif algorithm == "modsym":
            A = E.modular_symbol_space(sign=1)
            self.modsym = ModularSymbolMap(A)
        else:
            raise ValueError, "unknown algorithm '%s'"%algorithm
            
        # the denom must be a unit, given our hypotheses (and assumptions in code!)
        assert ZZ(self.modsym.denom)%self.p != 0, "internal modsym denominator must be a p(=%s)-unit, but it is %s"%(self.p, self.modsym.denom)
            
        # Compute alpha:
        f = ZZ['x']([p, -E.ap(p), 1])
        G = f.factor_padic(p, self.prec)
        Zp = None
        for pr, e in G:
            alpha = -pr[0]
            if alpha.valuation() < 1:
                Zp = alpha.parent()
                self._alpha = alpha.lift()
                break
        assert Zp is not None, "bug -- must have unit root (p=%s)"%p

        cdef int n
        K = Qp(self.p, self.prec//2)

        # Compute array of powers of inverse of alpha, modulo p^prec,
        # for each power up to prec/2.
        u = 1
        self.alpha_inv[0] = u
        for n in range(self.prec//2):
            u *= alpha
            self.alpha_inv[n+1] = K(1/u).lift()   # coerce to K to lower precision

        # Make a table of powers of p up to prec
        ppow = 1
        self.p_pow[0] = ppow
        for n in range(self.prec):
            ppow *= self.p
            self.p_pow[n+1] = ppow

        # Make a table of Teichmuller lifts
        self.teich = <long*> sage_malloc(sizeof(long)*self.p)
        self.teich[0] = 0
        for n in range(1,p):
            self.teich[n] = K.teichmuller(n).lift()

        # Compute normalization
        self.normalization = ZZ(self.E.real_components()).inverse_mod(self.p_pow[self.prec//2])

    def __repr__(self):
        return "%s-adic L-series of %s"%(self.p, self.E)

    def alpha(self):
        return self._alpha
        
    cpdef long modular_symbol(self, long a, long b):
        cdef long v[1]
        self.modsym.evaluate(v, a, b)
        return v[0]

    cpdef long measure(self, long a, int n):
        # TODO: case when p divides level is different
        return (self.alpha_inv[n] * self.modular_symbol(a, self.p_pow[n])
                   - self.alpha_inv[n+1] * self.modular_symbol(a, self.p_pow[n-1]))

    def _series(self, int n, prec, ser_prec=5, bint verb=0):
        """
        EXAMPLES::
        
            sage: import psage.modform.rational.padic_elliptic_lseries_fast as p; L = p.pAdicLseries(EllipticCurve('389a'),5)
            sage: f = L._series(2, 3, ser_prec=6); f.change_ring(Integers(5^3))
            73*T^4 + 42*T^3 + 89*T^2 + 120*T
            sage: f = L._series(3, 4, ser_prec=6); f.change_ring(Integers(5^3))
            36*T^5 + 53*T^4 + 47*T^3 + 99*T^2
            sage: f = L._series(4, 5, ser_prec=6); f.change_ring(Integers(5^3))
            61*T^5 + 53*T^4 + 22*T^3 + 49*T^2
            sage: f = L._series(5, 6, ser_prec=6); f.change_ring(Integers(5^3))
            111*T^5 + 53*T^4 + 22*T^3 + 49*T^2
        """
        cdef long a, b, j, s, gamma_pow, gamma, pp

        assert prec >= n, "prec (=%s) must be as large as approximation n (=%s)"%(prec, n)
        assert prec <= self.prec // 2, "requested precision (%s) too large (max: %s)"%(prec, self.prec//2)

        pp = self.p_pow[prec]

        gamma = 1 + self.p
        gamma_pow = 1

        R = Integers(pp)['T']
        T = R.gen()
        one_plus_T_factor = R(1)
        L = R(0)
        one_plus_T = 1+T

        _sig_on
        for j in range(self.p_pow[n-1]):
            s = 0
            for a in range(1, self.p):
                b = self.teich[a] * gamma_pow
                s += self.measure(b, n)
            L += (s * one_plus_T_factor).truncate(ser_prec)
            one_plus_T_factor = (one_plus_T*one_plus_T_factor).truncate(ser_prec)
            gamma_pow = (gamma_pow * gamma)%pp
            #if verb: print j, s, one_plus_T_factor, gamma_pow

        _sig_off
        return L * self.normalization

    def _prec_bounds(self, n, ser_prec):
        pn  = Integer(self.p_pow[n-1])
        enj = infinity
        res = [enj]
        for j in range(1, ser_prec):
            bino = binomial(pn,j).valuation(self.p)
            if bino < enj:
                enj = bino
            res.append(enj)
        return res
    
    def series(self, int n=2, prec=None, ser_prec=5, int check=True):
        """
        EXAMPLES::

            sage: from psage.modform.rational.padic_elliptic_lseries_fast import pAdicLseries
            sage: E = EllipticCurve('389a'); L = pAdicLseries(E,5)
            sage: L.series()
            O(5) + O(5)*T + (4 + O(5))*T^2 + (2 + O(5))*T^3 + (3 + O(5))*T^4 + O(T^5)
            sage: L.series(1)
            O(T^0)
            sage: L.series(3)
            O(5^2) + O(5^2)*T + (4 + 4*5 + O(5^2))*T^2 + (2 + 4*5 + O(5^2))*T^3 + (3 + O(5^2))*T^4 + O(T^5)
            sage: L.series(2, 8)
            O(5^6) + O(5)*T + (4 + O(5))*T^2 + (2 + O(5))*T^3 + (3 + O(5))*T^4 + O(T^5)
            sage: L.series(3, 8)
            O(5^6) + O(5^2)*T + (4 + 4*5 + O(5^2))*T^2 + (2 + 4*5 + O(5^2))*T^3 + (3 + O(5^2))*T^4 + O(T^5)
            sage: L.series(3, ser_prec=8)
            O(5^2) + O(5^2)*T + (4 + 4*5 + O(5^2))*T^2 + (2 + 4*5 + O(5^2))*T^3 + (3 + O(5^2))*T^4 + (1 + O(5))*T^5 + O(5)*T^6 + (4 + O(5))*T^7 + O(T^8)
        """
        if prec is None: prec = n+1
        p = self.p
        if check:
            assert self.E.galois_representation().is_surjective(p), "p (=%s) must be surjective for E"%p
        f = self._series(n, prec, ser_prec)
        aj = f.list()
        R = Zp(p, prec)
        if len(aj) > 0:
            bounds = self._prec_bounds(n, ser_prec)
            aj = [R(aj[0], prec-2)] + [R(aj[j], bounds[j]) for j in range(1,len(aj))]
        # make unknown coefficients show as 0 precision.
        aj.extend([R(0,0) for _ in range(ser_prec-len(aj))])
        ser_prec = min([ser_prec] + [i for i in range(len(aj)) if aj[i].precision_absolute() == 0])
        L = R[['T']](aj, ser_prec)
        return L / self.modsym.denom
        
    def series_modp(self, int n=2, ser_prec=5, int check=True):
        """
        EXAMPLES::

            sage: from psage.modform.rational.padic_elliptic_lseries_fast import pAdicLseries
            sage: E = EllipticCurve('389a')
            sage: L = pAdicLseries(E,5)
            sage: L.series_modp()
            4*T^2 + 2*T^3 + 3*T^4 + O(T^5)
            sage: L.series_modp(3, 8)
            4*T^2 + 2*T^3 + 3*T^4 + T^5 + 4*T^7 + O(T^8)
            sage: L.series_modp(2, 20)
            4*T^2 + 2*T^3 + 3*T^4 + O(T^5)
        """
        L = self.series(n=n, ser_prec=ser_prec, check=check)
        return L.change_ring(Integers(self.p))
