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

include 'stdsage.pxi'
include 'interrupt.pxi'

from sage.rings.all import ZZ, Qp, Integers

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
    
    def __init__(self, E, p):
        self.E = E
        self.p = p

        # prec = biggest n such that p^n <= 2^63, so n = log_p(2^63)
        self.prec = ZZ(2**63).exact_log(p)
        
        assert Integer(p).is_pseudoprime()
        assert E.is_ordinary(p)
        A = E.modular_symbol_space(sign=1)
        self.modsym = ModularSymbolMap(A)


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
        assert Zp is not None

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
    
    def series0(self, int n=2, prec=5, ser_prec=10):
        """
        EXAMPLES::
        
            sage: import psage.modform.rational.padic_elliptic_lseries_fast as p; L = p.pAdicLseries(EllipticCurve('389a'),5)
            sage: f = L.series0(2, ser_prec=6); f.change_ring(Integers(5^3))
            73*T^4 + 42*T^3 + 89*T^2 + 120*T
            sage: f = L.series0(3, ser_prec=6); f.change_ring(Integers(5^3))
            36*T^5 + 53*T^4 + 47*T^3 + 99*T^2
            sage: f = L.series0(4, ser_prec=6); f.change_ring(Integers(5^3))
            61*T^5 + 53*T^4 + 22*T^3 + 49*T^2
            sage: f = L.series0(5, ser_prec=6); f.change_ring(Integers(5^3))
            111*T^5 + 53*T^4 + 22*T^3 + 49*T^2
        """
        cdef long a, b, j, s, gamma_pow, gamma, pp

        assert prec < self.prec // 2

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
                b = (self.teich[a] * gamma_pow)%pp
                s += self.measure(b, n)
            L += (s * one_plus_T_factor).truncate(ser_prec)
            one_plus_T_factor = (one_plus_T*one_plus_T_factor).truncate(ser_prec)
            gamma_pow = (gamma_pow * gamma)%pp

        _sig_off
        return L * self.normalization
