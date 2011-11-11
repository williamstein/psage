"""

 (c) Copyright 2009-2010 Salman Baig and Chris Hall

 This file is part of ELLFF

 ELLFF is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 ELLFF is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

# Wrapper code for ellff library

include "stdsage.pxi"
include "cdefs.pxi"

from sage.libs.ntl.ntl_ZZ_decl cimport *
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.ntl_ZZ_p_decl cimport *
from sage.libs.ntl.ntl_ZZ_pX_decl cimport *
from sage.libs.ntl.ntl_ZZ_pE_decl cimport *
from sage.libs.ntl.ntl_ZZ_pEX_decl cimport *
from sage.libs.ntl.ntl_ZZ_pX cimport ntl_ZZ_pX
from sage.libs.ntl.ntl_ZZ_pE cimport ntl_ZZ_pE
from sage.libs.ntl.ntl_ZZ_pEX cimport ntl_ZZ_pEX
from sage.libs.ntl.ntl_lzz_pX_decl cimport *
from sage.libs.ntl.ntl_lzz_pX cimport *
from sage.libs.ntl.ntl_lzz_pContext cimport *
from sage.libs.ntl.ntl_lzz_pContext import ntl_zz_pContext
from sage.rings.all import PolynomialRing, GF, ZZ
from sage.rings.ring import is_Field
from sage.rings.integer import is_Integer
from sage.rings.arith import is_prime
from sage.rings.integer cimport Integer
#from sage.rings.fraction_field_element import FractionFieldElement, is_FractionFieldElement
from sage.rings.polynomial.polynomial_element import is_Polynomial
from sage.schemes.elliptic_curves.ell_generic import is_EllipticCurve
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.structure.sage_object import SageObject

####################
# PSAGE imports
from psage.function_fields import is_FunctionField, FunctionField
from psage.function_fields.function_field import RationalFunctionField
from psage.function_fields.function_field_element import FunctionFieldElement_rational
####################

import psage.ellff.euler_database as edb

####################

cdef extern from "ntl_wrap.h":
    void ZZ_to_mpz(mpz_t *, ZZ_c *)

####################

cdef extern from "NTL/lzz_pE.h":
    ctypedef struct zz_pEContext_c "struct zz_pEContext":
        void (*restore)()

    zz_pEContext_c* zz_pEContext_new "new zz_pEContext"(zz_pX_c p)
    void zz_pEContext_delete "delete "(zz_pEContext_c *mem)

    void zz_pEContext_restore(zz_pEContext_c *c)

    ctypedef struct zz_pE_c "struct zz_pE":
        pass

    zz_pE_c* zz_pE_new "new zz_pE"()
    void zz_pE_delete "delete "(zz_pE_c *mem)
    void zz_pE_conv "conv"(zz_pE_c x, zz_pX_c a)
    zz_pX_c zz_pE_rep "rep"(zz_pE_c z)

####################

cdef extern from "NTL/lzz_pEX.h":
    ctypedef struct zz_pEX_c "struct zz_pEX":
        pass

    zz_pEX_c* zz_pEX_new "new zz_pEX"()
    void zz_pEX_delete "delete "(zz_pEX_c *mem)

    void zz_pEX_SetCoeff "SetCoeff"(zz_pEX_c x, long i, zz_pE_c a)

####################

cdef extern from "lzz_pEratX.h":
    ctypedef struct zz_pEratX_c "zz_pEratX":
        pass

    zz_pEratX_c* zz_pEratX_new "new zz_pEratX"(zz_pEX_c num, zz_pEX_c den)
    void zz_pEratX_delete "delete "(zz_pEratX_c *mem)

####################

cdef extern from "euler.h":
    void __euler_table "euler_table"(long *table, long min_tau, long max_tau, int euler_deg) except+
    void __twist_table "twist_table"(zz_pEX_c f, long *untwisted_table, long *twisted_table, long min_tau, long max_tau)
    void __pullback_table "pullback_table"(zz_pEX_c finite_disc, zz_pEX_c infinite_disc, zz_pEratX_c f, long *base_table, long *pullback_table, long min_tau, long max_tau)
    void __sum_table "sum_table"(ZZ_c b_n, long *table, long min_tau, long max_tau)
    void __compute_b_n "compute_b_n"(ZZ_c b_n)
    void __compute_c_n "compute_c_n"(ZZ_c **b, ZZ_c **c, int n)

####################

cdef extern from "jacobi.h":
    void __jacobi_sum "jacobi_sum"(ZZ_c ***sum, long d)

####################

cdef extern from "ell_surface.h":
    ctypedef struct ell_surfaceInfoT_c "ell_surfaceInfoT":
        int     sign
        int     deg_L
        int     constant_f

    void ell_surface_init "ell_surface::init"(zz_pEratX_c a1, zz_pEratX_c a2, zz_pEratX_c a3, zz_pEratX_c a4, zz_pEratX_c a6)

    ctypedef struct ell_surface_c "ell_surface":
        int     sign
        int     deg_L
        int     constant_f

    void ell_surface_init "ell_surface::init"(zz_pEratX_c a1, zz_pEratX_c a2, zz_pEratX_c a3, zz_pEratX_c a4, zz_pEratX_c a6)
    ell_surfaceInfoT_c* ell_surface_getSurfaceInfo "ell_surface::getSurfaceInfo"()

    void ell_get_j_invariant "get_j_invariant"(ZZ_pEX_c j_num, ZZ_pEX_c j_den)
    void ell_get_an "get_an"(ZZ_pEX_c a4, ZZ_pEX_c a6)
    void ell_get_reduction "get_reduction"(ZZ_pEX_c **divisors, int finite)
    void ell_get_disc "get_disc"(ZZ_pEX_c disc, int finite_f)

####################

cdef extern from "helper.h":
    cdef void __get_modulus "get_modulus"(zz_pX_c pi_1, zz_pX_c pi_2, zz_pX_c a, int p, int d1, int d2)
    cdef void init_NTL_ff(int p, int d, int precompute_inverses, int precompute_square_roots, int precompute_legendre_char, int precompute_pth_frobenius_map)

####################

cdef class _ellff_EllipticCurve_c:
    cdef ZZ_c   **b, **b_star, **c, **c_star, **__b, **__c
    cdef int    bc_len

    cdef long   **trace_tables, *trace_table_sizes
    cdef int    num_trace_tables

    def __init__(self):
        pass

    def __cinit__(self):
        # allocate memory for tables of traces.  initially each
        # table is empty (and will need to be allocated).
        self.num_trace_tables = 15
        self.trace_tables = \
            <long**>sage_malloc(sizeof(long*)*self.num_trace_tables)
        self.trace_table_sizes  = \
            <long*>sage_malloc(sizeof(long)*self.num_trace_tables)
        for i in range(self.num_trace_tables):
            self.trace_tables[i] = NULL
            self.trace_table_sizes[i] = -1

        # allocate memory for b and c coefficients (for L-function)
        self.bc_len = 2*self.num_trace_tables
        self.b = <ZZ_c **>sage_malloc(sizeof(ZZ_c*)*(self.bc_len+1))
        self.c = <ZZ_c **>sage_malloc(sizeof(ZZ_c*)*(self.bc_len+1))
        self.b_star = \
            <ZZ_c **>sage_malloc(sizeof(ZZ_c*)*(self.bc_len+1))
        self.c_star = \
            <ZZ_c **>sage_malloc(sizeof(ZZ_c*)*(self.bc_len+1))

        for n in range(self.bc_len+1):
            self.b[n] = ZZ_new()
            self.c[n] = ZZ_new()
            self.b_star[n] = ZZ_new()
            self.c_star[n] = ZZ_new()

        ZZ_set_from_int(self.c[0], 1)
        ZZ_set_from_int(self.c_star[0], 1)

    def __dealloc__(self):
        for i in range(self.num_trace_tables):
            if self.trace_table_sizes[i] != -1:
                sage_free(self.trace_tables[i])
        sage_free(self.trace_tables)
        sage_free(self.trace_table_sizes)

        for i in range(self.bc_len+1):
            ZZ_delete(self.b[i])
            ZZ_delete(self.c[i])
            ZZ_delete(self.b_star[i])
            ZZ_delete(self.c_star[i])

        sage_free(self.b)
        sage_free(self.c)
        sage_free(self.b_star)
        sage_free(self.c_star)

    def _build_euler_table(self, n):
        self._allocate_table(n, self.q**n+1)

        self._surface_init(n)
        __euler_table(self.trace_tables[n], 0, self.q**n, self.d)

        cdef Integer b_n = Integer(0)
        for i in range(self.q**n+1):
            b_n += self.trace_tables[n][i]

        b_n._to_ZZ(self.b[n])

    def _twist_euler_table(self, f, _ellff_EllipticCurve_c E, n):
        cdef zz_pEX_c f_ = _to_zz_pEX_(f, self.phi, self.p)

        assert self.trace_table_sizes[n] == self.q**n+1
        assert E.trace_table_sizes[n]    == self.q**n+1
        __twist_table(f_, self.trace_tables[n], E.trace_tables[n],
            0, self.q**n)

        cdef Integer b_n = Integer(0)
        for i in range(self.q**n+1):
            b_n += E.trace_tables[n][i]

        b_n._to_ZZ(E.b[n])

    def _pullback_euler_table(self, f, _ellff_EllipticCurve_c E, n):
        cdef zz_pEratX_c *f_
        f_ = _to_EratX_(f, self.phi, self.p)

        cdef zz_pEX_c    fin_disc, inf_disc
        fin_disc = _to_zz_pEX_(self.__finite_disc, self.phi, self.p)
        inf_disc = _to_zz_pEX_(self.__infinite_disc, self.phi, self.p)

        __pullback_table(fin_disc, inf_disc, f_[0],
            self.trace_tables[n], E.trace_tables[n], 0, self.q**n)

        cdef Integer b_n = Integer(0)
        for i in range(self.q**n+1):
            b_n += E.trace_tables[n][i]

        b_n._to_ZZ(E.b[n])
        ZZ_sub(E.b_star[n][0], E.b[n][0], self.b[n][0])

    def _allocate_table(self, n, size):
        if self.trace_table_sizes[n] == -1:
            self.trace_tables[n] \
                = <long*>sage_malloc(sizeof(long)*size)
            self.trace_table_sizes[n] = size
        elif self.trace_table_sizes[n] != size:
            raise RuntimeError("table size is wrong")

    def _deallocate_table(self, i):
        assert self.trace_table_sizes[i] > 0
        sage_free(self.trace_tables[i])
        self.trace_table_sizes[i] = -1

    def _num_trace_tables(self):
        return self.num_trace_tables

    def _trace_table_sizes(self, n):
        return self.trace_table_sizes[n]

    def _compute_c_n(self, n):
        __compute_c_n(self.b, self.c, n)

    def _compute_rel_c_n(self, n):
        __compute_c_n(self.b_star, self.c_star, n)

    def __common_functional_equation(self, euler_deg, deg_L, sign, optimal_deg):
        cdef ZZ_c **b, **c
        b = self.__b
        c = self.__c

        cdef ZZ_c fec_c, m_c
        if euler_deg >= optimal_deg:
            if deg_L % 2 == 1:  m = self.q * sign
            else:               m =    1   * sign

            assert not(deg_L % 2 == 0 and sign==1) or euler_deg*2 >= deg_L

            ZZ_set_from_int(&m_c, m)
            for n in range((deg_L+deg_L%2)/2, deg_L+1):
                if n <= euler_deg:
                    if deg_L-n < n:
                        # compare what we calculated against f.e.
                        ZZ_mul(fec_c, m_c, c[deg_L-n][0])
                        if not ZZ_equal(c[n][0], fec_c):
                            e = 'Coefficient', n, \
                                'does not match. By hand:', \
                                ZZ_to_int(&c[n][0]), ', FE:', \
                                ZZ_to_int(&fec_c)
                            raise RuntimeError(e)

                elif deg_L % 2 == 0 and n == deg_L/2:
                    assert sign == -1
                    # f.e. implies middle coefficient vanishes
                    ZZ_set_from_int(c[n], 0)

                else:
                    ZZ_mul(c[n][0], m_c, c[deg_L-n][0])

                ZZ_mul_long(m_c, m_c, self.q * self.q)

    def _apply_functional_equation(self, euler_deg):
        self.__b = self.b
        self.__c = self.c
        self.__common_functional_equation(euler_deg=euler_deg,
            deg_L=self.deg_L,
            sign=self.sign,
            optimal_deg=self.optimal_deg)
        self.b = self.__b
        self.c = self.__c

        # perform a further consistency check, if possible
        if euler_deg > self.deg_L:
            for n in range(self.deg_L+1,euler_deg+1):
                if not ZZ_IsZero(self.c[n][0]):
                    return RuntimeError("fails a consistency check")
            self.c_len = euler_deg
        elif euler_deg >= self.optimal_deg:
            self.c_len = self.deg_L
        else:
            self.c_len = euler_deg

    def _apply_rel_functional_equation(self, euler_deg):
        self.__b = self.b_star
        self.__c = self.c_star
        self.__common_functional_equation(euler_deg=euler_deg,
            deg_L=self.rel_deg_L,
            sign=self.rel_sign,
            optimal_deg=self.rel_optimal_deg)
        self.b_star = self.__b
        self.c_star = self.__c

        # perform a further consistency check, if possible
        if euler_deg > self.rel_deg_L:
            for n in range(self.rel_deg_L+1,euler_deg+1):
                if not ZZ_IsZero(self.c[n][0]):
                    return RuntimeError("fails a consistency check")
            self.c_star_len = euler_deg
        elif euler_deg >= self.rel_optimal_deg:
            self.c_star_len = self.rel_deg_L
        else:
            self.c_star_len = euler_deg

    def _build_L_function(self, int euler_deg):
        if euler_deg < self.optimal_deg: deg_L = euler_deg
        else:                            deg_L = self.deg_L

        cdef Integer tmp
        v = []
        for n in range(deg_L+1):
            tmp = PY_NEW(Integer)
            ZZ_to_mpz(&tmp.value, self.c[n])
            v.append(tmp)

        R = PolynomialRing(ZZ, 'T')
        return R(v)

    def _build_relative_L_function(self, int euler_deg):
        if euler_deg < self.rel_optimal_deg: rel_deg_L = euler_deg
        else:                                rel_deg_L = self.rel_deg_L

        cdef Integer tmp
        v = []
        for n in range(rel_deg_L+1):
            tmp = PY_NEW(Integer)
            ZZ_to_mpz(&tmp.value, self.c_star[n])
            v.append(tmp)

        R = PolynomialRing(ZZ, 'T')
        return R(v)

    def __euler_table(self, n):
        q = self.q
        if self.trace_table_sizes[n] == -1:
            raise RuntimeError("table is empty")
        assert self.trace_table_sizes[n] == q**n+1

        v = []
        for i in range(0, q**n+1):
            v.append(self.trace_tables[n][i])
        return v

    def __set_euler_table(self, n, table, force=False):
        q = self.q
        if self.trace_table_sizes[n] != -1:
            if not force:
                raise RuntimeError("run with force=True to force an overwrite")
            self._deallocate_table(n)

        self._allocate_table(n, q**n+1)
        for i in range(q**n+1):
            self.trace_tables[n][i] = table[i]

        cdef Integer b_n = Integer(0)
        for i in range(self.q**n+1):
            b_n += self.trace_tables[n][i]

        b_n._to_ZZ(self.b[n])

    def _b(self, n):
        cdef Integer tmp
        tmp = PY_NEW(Integer)
        ZZ_to_mpz(&tmp.value, self.b[n])
        return tmp

    def _c(self, n):
        cdef Integer tmp
        tmp = PY_NEW(Integer)
        ZZ_to_mpz(&tmp.value, self.c[n])
        return tmp

    def _b_star(self, n):
        cdef Integer tmp
        tmp = PY_NEW(Integer)
        ZZ_to_mpz(&tmp.value, self.b_star[n])
        return tmp

    def _c_star(self, n):
        cdef Integer tmp
        tmp = PY_NEW(Integer)
        ZZ_to_mpz(&tmp.value, self.c_star[n])
        return tmp


class ellff_EllipticCurve(_ellff_EllipticCurve_c,SageObject):
    r"""
    The \class{ellff_EllipticCurve} class represents an ELLFF
    elliptic curve.
    """
    
    def __init__(self, field, ainvs):
        r"""
        Create the ellff elliptic curve with invariants
        \code{ainvs}, which is a list of 5 \emph{field elements}
        $a_1$, $a_2$, $a_3$, $a_4$, and $a_6$.

        INPUT:
            - field --  ground field
            - ainvs --  a list of 5 integers
                       
        EXAMPLES::

            sage: import psage
            sage: K.<t> = psage.FunctionField(GF(5))
            sage: E = psage.ellff_EllipticCurve(K,[0,0,0,t^2,t+1]); E
            <class 'psage.ellff.ellff.ellff_EllipticCurve'>
        """

        _ellff_EllipticCurve_c.__init__(self)

        if not is_Field(field):
            raise TypeError, "field must be a field."

        if not is_FunctionField(field):
            raise TypeError, "field must be of type FunctionField"

        if not isinstance(field, RationalFunctionField):
            raise NotImplementedError, "K must be the function field of P^1."

        # TODO: make this into a check for K=F_q(t)
        # true/false: the following words *only if* K =F_q(C)
        self.K = field
        self.R = self.K.maximal_order()

        F = self.F = self.K.constant_field()
        p = self.p = self.F.characteristic()
        d = self.d = self.F.degree()
        q = self.q = p**d

        self._E = EllipticCurve(field, ainvs)
        self.L_function_calculated = False;
        self.relative_L_function_calculated = False;

        if not isinstance(ainvs, list) and len(ainvs) == 5:
            raise TypeError, "ainvs must be a list of length 5."

        cdef int sign, deg_L, constant_f
        sign, deg_L, constant_f, phi = self._surface_init(1)

        self.sign    = sign
        self.deg_L   = deg_L
        self.constant_f = constant_f
        self.phi     = phi
        self.__calc_optimal_deg()
        self.euler_deg = 0
        self.rel_euler_deg = 0

    def _retrieve_invariants(self):
        r"""
        Convert the invariants associated to self (local reduction
        information for the finite and infinite models, the
        $a$ invariants (only $a_4$ and $a_6$ for now), the $j$-invariant, and
        the finite and infinite discriminant) from elements of Sage's
        constant field to elements of ELLFF's constant field via a
        field embedding.

        EXAMPLES::

            sage: import psage
            sage: K.<t> = psage.FunctionField(GF(11))
            sage: E = psage.ellff_EllipticCurve(K,[0,0,0,(t+1)*(t+2),t^2+1])
            sage: E._retrieve_invariants()

        """
        
        phi, F2 = _ellff_field_embedding(self.F, self.p, self.d, 1)

        # respective finite fields for sage and ellff
        sage_F  = self.F
        ellff_F = F2

        cdef ZZ_pEX_c **divisors
        divisors = <ZZ_pEX_c**>sage_malloc(sizeof(ZZ_pEX_c*)*9)
        for i in range(9):
            divisors[i] = ZZ_pEX_new()

        # get information about bad reduction away from infinity
        ell_get_reduction(divisors, 1)
        self._finite_reduction = []
        for i in range(9):
            self._finite_reduction.append(
                self.R(__from_ZZ_pEX(sage_F, ellff_F, divisors[i][0], self.R)))

        #  self._finite_reduction[0] = M_sp
        #  self._finite_reduction[1] = M_ns
        #  self._finite_reduction[2] = I^*
        #  self._finite_reduction[3] = II
        #  self._finite_reduction[4] = II^*
        #  self._finite_reduction[5] = III
        #  self._finite_reduction[6] = III^*
        #  self._finite_reduction[7] = IV
        #  self._finite_reduction[8] = IV^*

        # convert into course information about reduction
        self._finite_M = 1
        self._finite_A = 1
        for i in range(0,2):
            self._finite_M *= self._finite_reduction[i]
        for i in range(2,9):
            self._finite_A *= self._finite_reduction[i]

        # get information about reduction for minimal about infinity
        ell_get_reduction(divisors, 0)
        self._infinite_reduction = []
        for i in range(9):
            self._infinite_reduction.append(
                self.R(__from_ZZ_pEX(sage_F, ellff_F, divisors[i][0], self.R)))

        #  self._infinite_reduction[0] = M_sp
        #  self._infinite_reduction[1] = M_ns
        #  self._infinite_reduction[2] = I^*
        #  self._infinite_reduction[3] = II
        #  self._infinite_reduction[4] = II^*
        #  self._infinite_reduction[5] = III
        #  self._infinite_reduction[6] = III^*
        #  self._infinite_reduction[7] = IV
        #  self._infinite_reduction[8] = IV^*

        # convert into course information about reduction
        self._infinite_M = 1
        self._infinite_A = 1
        for i in range(0,2):
            self._infinite_M *= self._infinite_reduction[i]
        for i in range(2,9):
            self._infinite_A *= self._infinite_reduction[i]

        cdef ZZ_pEX_c *a4 = ZZ_pEX_new(), *a6 = ZZ_pEX_new()

        ell_get_an(a4[0], a6[0])
        self.a4 = self.R(__from_ZZ_pEX(self.F, F2, a4[0], self.R))
        self.a6 = self.R(__from_ZZ_pEX(self.F, F2, a6[0], self.R))

        cdef ZZ_pEX_c *j_num = ZZ_pEX_new(), *j_den = ZZ_pEX_new()

        ell_get_j_invariant(j_num[0], j_den[0])
        self.j = self.R(__from_ZZ_pEX(self.F, F2, j_num[0], self.R)) / self.R(__from_ZZ_pEX(self.F, F2, j_den[0], self.R))

        cdef ZZ_pEX_c *finite_disc   = ZZ_pEX_new()
        cdef ZZ_pEX_c *infinite_disc = ZZ_pEX_new()
        ell_get_disc(finite_disc[0],   1)
        ell_get_disc(infinite_disc[0], 0)

        self.__finite_disc   = self.R(__from_ZZ_pEX(self.F, F2, finite_disc[0], self.R))
        self.__infinite_disc = self.R(__from_ZZ_pEX(self.F, F2, infinite_disc[0], self.R))

        # clean up allocated memory
        for i in range(9):
            ZZ_pEX_delete(divisors[i])
        sage_free(divisors)

        ZZ_pEX_delete(a4)
        ZZ_pEX_delete(a6)
        ZZ_pEX_delete(j_num)
        ZZ_pEX_delete(j_den)
        ZZ_pEX_delete(finite_disc)
        ZZ_pEX_delete(infinite_disc)

    def _surface_init(self, n):
        r"""
        Initialize the elliptic surface associated to self over $F_{q^n}$.

        INPUT:

            - n -- an integer

        OUTPUT:

            - info.sign -- the sign of the functional equation
            - info.deg_L -- the degree of the L-function
            - info.constant_f -- True (1) if curve is constant and has no places of bad reduction
            - phi -- an embedding of Sage's F_q to ELLFF's F_q

        EXAMPLES::

            sage: import psage
            sage: K.<t> = psage.FunctionField(GF(11))
            sage: E = psage.ellff_EllipticCurve(K,[0,0,0,(t+1)*(t+2),t^2+1])
            sage: E._surface_init(2)
            [1, 4, 0, Ring morphism:
              From: Finite Field of size 11
              To:   Finite Field in c of size 11^2
              Defn: 1 |--> 1]

        """
        
        p = self.p
        d = self.d

        phi, F2 = _ellff_field_embedding(self.F, p, d, n)

        # convert coefficients to format acceptible for ellff library
        cdef zz_pEratX_c *a1_, *a2_, *a3_, *a4_, *a6_

        a1_ = _to_EratX_(self._E.a1(), phi, p)
        a2_ = _to_EratX_(self._E.a2(), phi, p)
        a3_ = _to_EratX_(self._E.a3(), phi, p)
        a4_ = _to_EratX_(self._E.a4(), phi, p)
        a6_ = _to_EratX_(self._E.a6(), phi, p)

        ell_surface_init(a1_[0], a2_[0], a3_[0], a4_[0], a6_[0])

        zz_pEratX_delete(a1_)
        zz_pEratX_delete(a2_)
        zz_pEratX_delete(a3_)
        zz_pEratX_delete(a4_)
        zz_pEratX_delete(a6_)

        init_NTL_ff(p, d*n, 0, 0, 0, 0)

        cdef ell_surfaceInfoT_c *info
        info = ell_surface_getSurfaceInfo()

        # TODO: add context save (and later restore)
        # cdef ell_surfaceContext_c *context

        if n == 1:
            self._retrieve_invariants()

        return [info.sign, info.deg_L, info.constant_f, phi]

    def __calc_optimal_deg(self):
        r"""
        Calculate and assign the optimal degree to compute the
        $L$-function of self, i.e. what is the maximum degree of the
        extension of $\mathbb{F}_q$ needed to determine the
        $L$-function.

        EXAMPLES::

            sage: import psage
            sage: K.<t> = psage.FunctionField(GF(11))
            sage: E = psage.ellff_EllipticCurve(K,[0,0,0,(t+1)*(t+2),t^2+1])
            sage: E.__calc_optimal_deg()

        """
        
        if self.deg_L % 2 == 1:
            self.optimal_deg = (self.deg_L-1)/2;
        elif self.sign == 1:
            self.optimal_deg = self.deg_L/2
        else:
            self.optimal_deg = self.deg_L/2 - 1

    def quadratic_twist(self, f, tables=False, euler_deg=0, force=False, verbose=False):
        r"""
        Create a new elliptic curve which is a quadratic twist of this
        curve.  If specified, calculates tables of Euler factors for
        twisted curve using tables of this curve.

        INPUT:
            f -- twisting polynomial (twist by extension K(sqrt{f})/K)
            tables -- if True, computes Euler factors of twist.
                (default: False)
            euler_deg -- if positive and tables==True, a bound on the
                degree of the Euler factors to compute.  if zero, becomes
                the optimal Euler degree for computing the L-function.
                (default: 0)
            force -- if True, make tables of Euler factors for this
                curve, if necessary, as well as for new curve.
                otherwise, use brute force.
                (default: False)
            verbose -- if True, outputs progress information.
                (default: False)
                       
        EXAMPLES::

            sage: import psage.ellff.ellff as ellff
            sage: K.<t> = FunctionField(GF(5))
            sage: E = ellff.ellff_EllipticCurve(K,[0,0,0,t^2,t+1]); E
            <class 'ellff.ellff_EllipticCurve'>
            sage: E_twist = E.quadratic_twist(t^2+1); E_twist
            <class 'ellff.ellff_EllipticCurve'>
            sage: E_twist.a4
            t^8 + 2*t^6 + t^4
            sage: E_twist.a6
            t^9 + t^8 + 3*t^7 + 3*t^6 + 3*t^5 + 3*t^4 + t^3 + t^2

        """

        q = self.q

        self._surface_init(1)

        E = ellff_EllipticCurve(self.K, [0, 0, 0, self.a4*f*f, self.a6*f*f*f])

        cdef ntl_ZZ   b_n = PY_NEW(ntl_ZZ)
        if tables:
            if euler_deg == 0:
                euler_deg = E.optimal_deg

            for n in range(1,euler_deg+1):
                if force and self._trace_table_sizes(n) == -1:
                    if verbose:
                        print "rebuilding own Euler table (n = ", n, ")"
                    self._build_euler_table(n)

                if self._trace_table_sizes(n) != -1:
                    if verbose:
                        print "twisting old table into new (n = ", n, ")"
                    E._surface_init(n)
                    E._allocate_table(n, q**n+1)
                    self._twist_euler_table(f, E, n)

                else:
                    if verbose:
                        print "building new Euler table from scratch (n = ", n, ")"
                    E._build_euler_table(n)

        return E

    def pullback(self, f, tables=False, euler_deg=-1, verbose=False):
        r"""
        Create a new elliptic curve which is a pullback of this
        curve.  If specified, calculates tables of Euler factors for
        pullback curve using tables of this curve.  Will also
        compute 'relative L-function', that is, the quotient
        L(T,E_pullback)/L(T,E).

        INPUT:
            f -- Non-constant rational function.  Gives rise to field
                embedding F_q(f)-->F_q(t).   Treat current curve as defined
                over F_q(f) and pullback along field extension F_q(t)/F_q(t)

            tables -- if True, computes Euler factors of pullback.  When
                base curve has tables of Euler factors, also computes
                'relative b_n' (whose 'relative c_n' are coefficients of
                the 'relative L-function').
                (default: False)

            euler_deg -- if non-negative and tables==True, equals
                bound on the degree of Euler factors to compute.  if negative,
                becomes the optimal degree.
                (default: None)

            verbose -- if True, outputs progress information.
                (default: False)
                       
        EXAMPLES::
        
            sage: import psage.ellff.ellff as ellff
            sage: K.<t> = FunctionField(GF(5))
            sage: E = ellff.ellff_EllipticCurve(K,[0,0,0,t^2,t+1]); E
            <class 'ellff.ellff_EllipticCurve'>
            sage: E_pullback = E.pullback(t^3); E_pullback
            <class 'ellff.ellff_EllipticCurve'>
            sage: E_pullback.a4
            t^4
            sage: E_pullback.a6
            t^3 + t^2
            sage: E_pullback = E.pullback(t^3, tables=True, verbose=True); E_pullback
            rebuilding own Euler table  (n =  1 )
            pulling back old table      (n =  1 )
            <class 'ellff.ellff_EllipticCurve'>

        """

        q = self.q

        self._surface_init(1)

        new_a4 = self.a4.numerator().subs(f) / self.a4.denominator().subs(f)
        new_a6 = self.a6.numerator().subs(f) / self.a6.denominator().subs(f)
        E = ellff_EllipticCurve(self.K, [0, 0, 0, new_a4, new_a6])

        E.rel_deg_L = E.deg_L - self.deg_L
        E.rel_sign  = E.sign / self.sign
        E.rel_optimal_deg = E.optimal_deg - self.optimal_deg

        cdef ntl_ZZ   b_n = PY_NEW(ntl_ZZ)
        if tables:
            if euler_deg < 0:
                euler_deg = E.rel_optimal_deg
            E.rel_euler_deg = 0

            for n in range(1,euler_deg+1):
                if self._trace_table_sizes(n) == -1:
                    if verbose:
                        print "rebuilding own Euler table  (n = ", n, ")"
                    self._build_euler_table(n)

                if self._trace_table_sizes(n) != -1:
                    if verbose:
                        print "pulling back old table      (n = ", n, ")"
                    E._surface_init(n)
                    E._allocate_table(n, q**n+1)
                    self._pullback_euler_table(f, E, n)

                    if E.rel_euler_deg == n-1:
                        E.rel_euler_deg = n

                else:
                    if verbose:
                        print "building new Euler table    (n = ", n, ")"
                    E._build_euler_table(n)

        return E

    def relative_L_function(self, verbose=False):
        r"""
        If not done so already, calculate the relative L-function
        (or part of it) of the elliptic curve.

        INPUT:
            euler_deg -- if non-negative, a bound on the degree of the
                Euler factors to compute.  if negative, becomes the optimal
                Euler degree for computing the L-function.
                (default: -1)
                       
        EXAMPLES::
        
            sage: import psage.ellff.ellff as ellff
            sage: K.<t> = FunctionField(GF(5))
            sage: E = ellff.ellff_EllipticCurve(K,[0,0,0,t^2,t+1]); E
            <class 'ellff.ellff_EllipticCurve'>
            sage: E_pullback = E.pullback(t^3, tables=True, verbose=True); E_pullback
            rebuilding own Euler table  (n =  1 )
            pulling back old table      (n =  1 )
            <class 'ellff.ellff_EllipticCurve'>
            sage: E_pullback.relative_L_function()
            125*T^3 - 5*T^2 - T + 1

        """

        if self.relative_L_function_calculated == True:
            return self._relative_L_function

        # retrieve useful constants
        p = self.p
        d = self.d
        q = self.q = p**d

        rel_deg_L = self.rel_deg_L
        rel_sign  = self.rel_sign
        rel_optimal_deg = self.rel_optimal_deg
        rel_euler_deg = self.rel_euler_deg

        # determine relation between rel_euler_deg and rel_optimal_deg
        if rel_optimal_deg > rel_euler_deg:
            print "rel_euler_deg (= ", rel_euler_deg, " ) ", \
                "is too small to compute entire L-function; ", \
                "need rel_euler_deg >= ", rel_optimal_deg

        if rel_euler_deg >= self._num_trace_tables():
            raise RuntimeError("rel_euler_deg is too large")

        # calculate necessary tables
        for n in range(1, rel_euler_deg+1):
            self._compute_rel_c_n(n)

        # take functional equation into account, if possible
        self._apply_rel_functional_equation(self.rel_euler_deg)

        # produce output polynomial
        self._relative_L_function = self._build_relative_L_function(rel_euler_deg)

        if rel_euler_deg >= self.rel_optimal_deg:
            self.relative_L_function_calculated = True

        return self._relative_L_function

    def L_function(self, int euler_deg=-1, force=False, verbose=False):
        r"""
        If not done so already, calculate the L-function (or part of it)
        of the elliptic curve.

        INPUT:
            euler_deg -- if non-negative, a bound on the degree of the
                Euler factors to compute.  if negative, becomes the optimal
                Euler degree for computing the L-function.
                (default: 0)
            force -- if True, forces recalculation of each table of
                Euler factors used.  if False, tables are only
                constructed as needed.
                (default: False)
                       
        EXAMPLES::
        
            sage: import psage.ellff.ellff as ellff
            sage: K.<t> = FunctionField(GF(5))
            sage: E = ellff.ellff_EllipticCurve(K,[0,0,0,t^2,t+1]); E
            <class 'ellff.ellff_EllipticCurve'>
            sage: E.L_function()
            625*T^4 + 125*T^3 + 5*T + 1

            sage: E_pullback = E.pullback(t^3, tables=True, verbose=True); E_pullback
            rebuilding own Euler table  (n =  1 )
            pulling back old table      (n =  1 )
            <class 'ellff.ellff_EllipticCurve'>
            sage: E_pullback.L_function()
            78125*T^7 + 1

            sage: E_twist = E.quadratic_twist(t^2+1,tables=True,force=True); E_twist
            <class 'ellff.ellff_EllipticCurve'>
            sage: E_twist.L_function()
            1953125*T^9 - 312500*T^8 - 31250*T^7 + 14750*T^6 - 10*T^5 - 2*T^4 + 118*T^3 - 10*T^2 - 4*T + 1

        """

        if self.L_function_calculated == True:
            return self._L_function

        if self.constant_f == True:
            raise ValueError, "Elliptic curve must be non-constant or have at least one place of bad reduction"

        # retrieve useful constants
        p = self.p
        d = self.d
        q = self.q = p**d

        deg_L = self.deg_L
        sign  = self.sign
        optimal_deg = self.optimal_deg

        # determine relation between euler_deg and optimal_deg
        if euler_deg < 0:
            euler_deg = optimal_deg
        elif deg_L > euler_deg*2:
            print "euler_deg (= ", euler_deg, " ) ", \
                "is too small to compute entire L-function; ", \
                "need euler_deg >= ", optimal_deg

        if euler_deg >= self._num_trace_tables():
            raise RuntimeError("euler_deg is too large")

        # calculate necessary tables
        for n in range(1, euler_deg+1):
            if force and self._trace_table_sizes(n) != -1:
                self._deallocate_table(n)

            if self._trace_table_sizes(n) == -1:
                if verbose:
                    print "building euler table (n = ", n, ")"
                self._allocate_table(n, q**n+1)
                self._build_euler_table(n)
            elif verbose:
                print "already have euler table (n = ", n, ")"

            self._compute_c_n(n)

        if force or euler_deg >= self.euler_deg:
            self.euler_deg = euler_deg

        # take functional equation into account, if possible
        self._apply_functional_equation(euler_deg)

        # produce output polynomial
        self._L_function = self._build_L_function(euler_deg)

        if euler_deg >= optimal_deg:
            self.L_function_calculated = True

        return self._L_function

    def _euler_table(self, n):
        r"""
        Returns a python list representing the table of Euler factors


        INPUT:
            n -- which degree

        EXAMPLES::
        
            sage: import psage.ellff.ellff as ellff
            sage: K.<t> = FunctionField(GF(5))
            sage: E = ellff.ellff_EllipticCurve(K,[0,0,0,t^2,t+1]); E
            <class 'ellff.ellff_EllipticCurve'>
            sage: E._euler_table(1)
            Traceback (most recent call last):
                ...
            RuntimeError: table is empty
            sage: E_pullback = E.pullback(t^3, tables=True, verbose=True); E_pullback
            rebuilding own Euler table  (n =  1 )
            pulling back old table      (n =  1 )
            <class 'ellff.ellff_EllipticCurve'>
            sage: E._euler_table(1)
            [0, 2, 3, -2, 2, 0]
            sage: E_pullback._euler_table(1)
            [0, 2, 1, -1, 2, 0]
            sage: E_twist = E.quadratic_twist(t^2+1, tables=True, force=True, verbose=True); E_twist
            twisting old table into new (n =  1 )
            twisting old table into new (n =  2 )
            twisting old table into new (n =  3 )
            twisting old table into new (n =  4 )
            <class 'ellff.ellff_EllipticCurve'>
            sage: E._euler_table(1)
            [0, 2, 3, -2, 2, 0]
            sage: E._euler_table(2)
            [-10, -6, -1, -6, -6, -4, 1, 4, -1, 1, -4, -1, 6, 1, -1, -4, -1, 6, 1, -1, -4, 1, 4, -1, 1, 0]
            sage: E._euler_table(3)[0:10]
            [0, -22, -18, 22, -22, 18, -12, -2, 3, 7]

            sage: E_twist._euler_table(1)
            [0, -2, 0, 0, -2, 0]
            sage: E_twist._euler_table(2)
            [0, -6, 0, 0, -6, -4, -1, 4, -1, -1, -4, 1, -6, -1, 1, -4, 1, -6, -1, 1, -4, -1, 4, -1, -1, 0]
            sage: E_twist._euler_table(3)[0:10]
            [0, 22, 0, 0, 22, 18, 12, 2, 3, 7]

        """
        return self.__euler_table(n)

    def _set_euler_table(self, n, table, force=False):
        r"""
        Sets the euler table of self over $F_{q^n}$ to table if table
        not already computed. If euler table has already been
        computed, it can be overwritten with force. Will also compute
        the corresponding exponential coefficients of the
        $L$-function, $b_n$.

        INPUT:

            - n -- the degree of the extenstion of F_q
            - table -- a table of euler factors
            - force -- a boolean that forces overwriting of existing euler table

        EXAMPLES::

            sage: import psage
            sage: K.<t> = psage.FunctionField(GF(11))
            sage: E = psage.ellff_EllipticCurve(K,[0,0,0,(t+1)*(t+2),t^2+1])
            sage: E._euler_table(1)
            Traceback (most recent call last):
                ...
            RuntimeError: table is empty
            sage: E.L_function()
            14641*T^4 + 1
            sage: E._euler_table(1)
            [-4, -2, 1, -4, 3, -6, 3, 5, 4, 0, 0, 0]
            sage: et = E._euler_table(1)
            sage: E._set_euler_table(1,et)
            Traceback (most recent call last):
                ...
            RuntimeError: run with force=True to force an overwrite
            sage: E._set_euler_table(1,table=et,force=True)
            
        """
        return self.__set_euler_table(n, table, force)

    def _save_euler_table(self, n, verbose=False):
        r"""
        Save the euler table for self over the degree n extension of
        $\mathbb{F}_q$ to disk. This is currently implemented with
        sage.database.db, which uses ZODB. If self is the versal
        j-curve, it stores the table in the database
    
        SAGE_ROOT/data/jcurve_euler_tables .

        Otherwise, the tables are stored in the `user` table
        
        SAGE_ROOT/data/local_euler_tables .

        It currently doesn't check if the table already is stored; it
        merely writes over it in that case. This should eventually be
        implemented using MongoDB.

        INPUT:
        
            - n -- the degree of the extension of F_q

        EXAMPLES::

            sage: import psage
            sage: K.<t> = psage.FunctionField(GF(11))
            sage: E = psage.ellff_EllipticCurve(K,[0,0,0,-27*t/(t-1728),54*t/(t-1728)])
            sage: E._build_euler_table(1)
            sage: E._euler_table(1)
            [0, 0, 4, -6, 3, 5, 1, -2, 4, -2, 3, 1]
            sage: E._build_euler_table(2)
            sage: E._build_euler_table(3)
            sage: E._euler_table(1)
            [0, 0, 4, -6, 3, 5, 1, -2, 4, -2, 3, 1]
            sage: E._save_euler_table(1)
            sage: E._save_euler_table(2)
            sage: E._save_euler_table(3)
        
    """
        edb._save_euler_table(self, n, verbose)

    def _load_euler_table(self, n, force=False, verbose=False):
        r"""
        Load the euler table for self over the degree n extension of
        $\mathbb{F}_q$ to disk. If self is the versal j-curve, the
        table is pulled from
    
        SAGE_ROOT/data/jcurve_euler_tables .
    
        Otherwise, the table is pulled from the `user` table

        SAGE_ROOT/data/local_euler_tables .

        This should eventually be implemented using MongoDB.

        It currently doesn't check if the key exist. If the key
        doesn't exist, a RuntimeError is raised by
        sage.database.db. This RuntimeError should be sufficient, so
        key checking may not be necessary.

        INPUT:

            - n -- the degree of the extension of F_q

        EXAMPLES::

            sage: import psage
            sage: K.<t> = psage.FunctionField(GF(11))
            sage: E = psage.ellff_EllipticCurve(K,[0,0,0,-27*t/(t-1728),54*t/(t-1728)])
            sage: E._euler_table(1)
            Traceback (most recent call last):
                ...
            RuntimeError: table is empty
            sage: E._load_euler_table(1)
            sage: E._euler_table(1)
            [0, 0, 4, -6, 3, 5, 1, -2, 4, -2, 3, 1]
            
        """

        edb._load_euler_table(self, n, force, verbose)

    def M(self,fine=False):
        r"""
        Returns the divisor of multiplicative reduction.

        INPUT:
            fine -- if True, returns pair of following outputs,
                one for split-multiplicative reduction and the
                other for non-split multiplicative reduction.

                (default: False)

        OUTPUT:
            ( finite, infinite_degree )

            finite -- polynomial in function field variable.
               the zeros are simple and precisely the places of
               (split or non-split) multiplicative reduction away
               from infinity.

            infinite_degree -- 1 if there is (split or non-split)
              multiplicative reduction at infinity, and 0 otherwise.
                       
        EXAMPLES::

            sage: import psage
            sage: K.<t> = psage.FunctionField(GF(11))
            sage: E = psage.ellff_EllipticCurve(K,[0,0,0,(t+1)*(t+2),t^2+1])
            sage: E.M()
            (t^6 + 9*t^5 + 4*t^4 + 8*t^3 + 8*t^2 + 3*t + 1, 0)

        """
        if fine == False:
            if self._infinite_M(0) == 0:    return (self._finite_M, 1)
            else:                           return (self._finite_M, 0)

        if self._infinite_reduction[0](0) == 0:
            p1 = (self._finite_reduction[0], 1)
        else:
            p1 = (self._finite_reduction[0], 0)

        if self._infinite_reduction[1](0) == 0:
            p2 = (self._finite_reduction[1], 1)
        else:
            p2 = (self._finite_reduction[1], 0)

        return (p1,p2)

    def A(self,fine=False):
        r"""
        Returns the divisor additive reduction as a polynomial,
        for the finite part, and an integer.

        INPUT:
            fine -- if True, returns tuple of following outputs,
                one for each type of bad reduction:
                    I_n^*, II, II^*, III, III^*, IV, IV^*

                (default: False)

        OUTPUT:
            ( finite, infinite_degree )

            finite -- polynomial in function field variable.
               the zeros are simple and precisely the places of
               (special) additive reduction away from infinity.

            infinite_degree -- 1 if there is (special) additive
              reduction at infinity, and 0 otherwise.
                       
        EXAMPLES::

            sage: import psage
            sage: K.<t> = psage.FunctionField(GF(11))
            sage: E = psage.ellff_EllipticCurve(K,[0,0,0,(t+1)*(t+2),t^2+1])
            sage: E.A()
            (1, 1)

        """
        if fine == False:
            if self._infinite_A(0) == 0:    return (self._finite_A, 1)
            else:                           return (self._finite_A, 0)

        if self._infinite_reduction[2](0) == 0:
            p1 = (self._finite_reduction[2], 1)
        else:
            p1 = (self._finite_reduction[2], 0)

        if self._infinite_reduction[3](0) == 0:
            p2 = (self._finite_reduction[3], 1)
        else:
            p2 = (self._finite_reduction[3], 0)

        if self._infinite_reduction[4](0) == 0:
            p3 = (self._finite_reduction[4], 1)
        else:
            p3 = (self._finite_reduction[4], 0)

        if self._infinite_reduction[5](0) == 0:
            p4 = (self._finite_reduction[5], 1)
        else:
            p4 = (self._finite_reduction[5], 0)

        if self._infinite_reduction[6](0) == 0:
            p5 = (self._finite_reduction[6], 1)
        else:
            p5 = (self._finite_reduction[6], 0)

        if self._infinite_reduction[7](0) == 0:
            p6 = (self._finite_reduction[7], 1)
        else:
            p6 = (self._finite_reduction[7], 0)

        if self._infinite_reduction[8](0) == 0:
            p7 = (self._finite_reduction[8], 1)
        else:
            p7 = (self._finite_reduction[8], 0)

        return (p1,p2,p3,p4,p5,p6,p7)

    r"""
    def __repr__(self):
        a = self.ainvs()
        s = "y^2"
        if a[0] == -1:
            s += "- x*y "
        elif a[0] == 1:
            s += "+ x*y "
        elif a[0] != 0:
            s += "+ %s*x*y "%a[0]
        if a[2] == -1:
            s += " - y"
        elif a[2] == 1:
            s += "+ y"
        elif a[2] != 0:
            s += "+ %s*y"%a[2]
        s += " = x^3 "
        if a[1] == -1:
            s += "- x^2 "
        elif a[1] == 1:
            s += "+ x^2 "
        elif a[1] != 0:
            s += "+ %s*x^2 "%a[1]
        if a[3] == -1:
            s += "- x "
        elif a[3] == 1:
            s += "+ x "
        elif a[3] != 0:
            s += "+ %s*x "%a[3]
        if a[4] == -1:
            s += "-1"
        elif a[4] == 1:
            s += "+1"
        elif a[4] != 0:
            s += "+ %s"%a[4]
        s = s.replace("+ -","- ")
        return s
    """

    def conductor(self):
        r"""
        Return the conductor of this curve
        """
        raise NotImplementedError

####################

##  def my_func1(a1):
##     if is_FractionFieldElement(a1):
##         if is_Polynomial(a1.numerator()):
##             return 1
##     elif is_Polynomial(a1):
##         return 1
##     else:
##         raise TypeError, "Input must be rational function"

cdef zz_pE_c _to_zz_pE_(z, phi, p):
    cdef zz_pE_c*  z_
    cdef ntl_zz_pX p_
    i_ = phi(z)
    if i_.category().object().degree() > 1:
        p_ = ntl_zz_pX(i_._vector_(), p)
    else:
        p_ = ntl_zz_pX([i_], p)
    z_ = zz_pE_new()
    zz_pE_conv(z_[0], p_.x)
    return z_[0]

## def _from_zz_pE(ntl_zz_pE z, F):
##     r"""
##     blah
## 
##     INPUT:
##         z -- element to embed.
##         F -- field to embed z into.  must be isomorphic to field of z.
## 
##             (default: 0)
##                    
##     EXAMPLES:
##     """
##     p = F.characteristic()
##     d = F.degree()
## 
##     # void zz_pE_conv "conv"(zz_pE_c x, zz_pX_c a)
##     # zz_pX_c zz_pE_rep "rep"(zz_pE_c z)


cdef zz_pEX_c _to_zz_pEX_(f, phi, p):
    cdef zz_pEX_c* f_

    f_ = zz_pEX_new()
    
    # make sure f is really a polynomial
    if is_FunctionField(f.parent().fraction_field()):
        if f.denominator().degree() > 0:
            raise TypeError("f must be a polynomial")
        # it is, but it looks like a rational function
        f = f.numerator()
 
    for i in range(f.degree()+1):
        zz_pEX_SetCoeff(f_[0], i, _to_zz_pE_(f.coeffs()[i], phi, p))
        
    return f_[0]


cdef zz_pEratX_c* _to_EratX_(r, phi, p):
    cdef zz_pEX_c n, d
    n = _to_zz_pEX_(r.numerator(), phi, p)
    d = _to_zz_pEX_(r.denominator(), phi, p)
    cdef zz_pEratX_c* r_
    # TODO: does this create a memory leak?
    r_ = zz_pEratX_new(n, d)
    return r_


    r"""
    blah

    INPUT:
        blah --

            (default: 0)
                   
    EXAMPLES:
    """

cdef __from_ZZ_pE(sage_F, ellff_F, ZZ_pE_c z):
    phi = ellff_F.Hom(sage_F)([ellff_F.gen().minpoly().change_ring(sage_F).roots()[0][0]])

    cdef ZZ_pX_c  z_
    z_ = ZZ_pE_rep(z);
    d = ZZ_pX_deg(z_)
    v = []
    cdef long l
    for i in range(d+1):
        ZZ_conv_to_long(l, ZZ_p_rep(ZZ_pX_coeff(z_, i)))
        v.append(l)
    if (ZZ_pE_IsZero(z)):
        v = [0]
    if sage_F.degree() > 1:
        R = ZZ[ellff_F.variable_name()]
        v_ = R(v)
        return phi(ellff_F(v_))
    return sage_F(v[0])

cdef __from_ZZ_pEX(sage_F, ellff_F, ZZ_pEX_c f, R):
    r"""
    blah

    INPUT:
        ring -- polynomial ring to embed in

            (default: 0)
                   
    EXAMPLES:
    """
    
    d = ZZ_pEX_deg(f)
    cdef ZZ_pE_c   c
    if not is_FunctionField(R.fraction_field()):
        raise TypeError("R.fraction_field() must be an order in a function field")
        # v = []
        # for i in range(d+1):
            # c = ZZ_pEX_coeff(f, i)
            # C = __from_ZZ_pE(sage_F, ellff_F, c)
            # v.append(C)
            #
        # return v
    s = R.gen()
    poly = R(0)
    for i in reversed(range(d+1)):
        c = ZZ_pEX_coeff(f, i)
        C = __from_ZZ_pE(sage_F, ellff_F, c)
        poly = s*poly + C
    return poly
        
def _ellff_field_embedding(F, p, d, n):
    r"""
    Return an embedding of F, an arbitrary finite field in Sage, to
    F2, a Sage finite field compatible with ELLFF. This embedding
    respects how Sage and ELLFF both treat relative extensions.

    INPUT:
        - F -- a finite field in Sage
        - p -- the characteristic of F
        - d -- F_q = F_{p^d}
        - n -- F2 = F_{q^n}

    OUTPUT:
        - phi -- the homomorphism F --> F2
        - F2 -- the ELLFF finite field as a Sage object

    EXAMPLES::

        sage: import psage
        sage: psage.ellff.ellff._ellff_field_embedding(GF(11), 11, 2, 3)
        [Ring morphism:
          From: Finite Field of size 11
          To:   Finite Field in c of size 11^6
          Defn: 1 |--> 1, Finite Field in c of size 11^6]
        
    """
    cdef ntl_zz_pX pi_1  = ntl_zz_pX(modulus=p)
    cdef ntl_zz_pX pi_2  = ntl_zz_pX(modulus=p)
    cdef ntl_zz_pX alpha = ntl_zz_pX(modulus=p)

    # get modulus from ellff library for representing F_q=F_{p^d},F_{q^n}
    __get_modulus(pi_1.x, pi_2.x, alpha.x, p, d, n)

    # get embedding F_q-->F_{q^r}, r=relative_degree
    if d == 1:
        F2 = GF(p**pi_2.degree(), name='c', modulus=ZZ['x'](pi_2.list()))
        phi = F.Hom(F2)([F.gen().minpoly().change_ring(F2).roots()[0][0]])
    else:
        F1 = GF(p**pi_1.degree(), name='b', modulus=ZZ['x'](pi_1.list()))
        F2 = GF(p**pi_2.degree(), name='c', modulus=ZZ['x'](pi_2.list()))

        phi1 = F.Hom(F1)([F.gen().minpoly().change_ring(F1).roots()[0][0]])
        phi2 = F1.Hom(F2)(ZZ['x'](alpha.list()))
        phi = phi1.post_compose(phi2)

    return [phi, F2]

def jacobi_sum(p, n, d, verbose=False):
    r"""
    Returns a polynomial which represents an abstract Jacobi sum
    over the field F_q=F_{p^n}.

    Suppose d|q-1 and let chi : F_q^* --> Z/d be the composition
    of x|-->x^{(q-1)/d} and a chosen isomorphism
    F_q^*/F_q^{*d} --> Z/d.  let z1,z2 be independent generic
    elements satisfying z1^d=z2^d=1.  The routine calculates the
    'abstract' Jacobi sum

        J(d,q) = sum_{x!=0,1} z1^{chi(x)}*z2^{chi(1-x)}

    and returns result as array A=(a_ij) where
    
        a_ij = coefficient of z1^i*z2^j in J(d,q)

    INPUT:

        - d -- the abstract order of the characters; must divide q-1

    OUTPUT:

        A

    EXAMPLES::

        sage: import psage
        sage: psage.ellff.ellff.jacobi_sum(11,3,4)
        Traceback (most recent call last):
            ...
        RuntimeError: d must divide q-1
        sage: psage.ellff.ellff.jacobi_sum(11,1,2)
        [[2, 2], [2, 3]]
        sage: psage.ellff.ellff.jacobi_sum(11,2,2)
        [[29, 30], [30, 30]]
        sage: psage.ellff.ellff.jacobi_sum(11,3,2)
        [[332, 332], [332, 333]]

    """

    if (not is_Integer(p)) or (not is_prime(p)):
        raise RuntimeError("p must be prime")

    if (not is_Integer(n)) or (n <= 0):
        raise RuntimeError("n must be a positive integer")

    q = p ** n
    if (q-1) % d != 0:
        raise RuntimeError("d must divide q-1")

    # calculate order of q mod d
    x = q % d
    f = 1
    while x != 1:
        x = (x*q) % d
        f = f+1
    if verbose:
        print "q = ", q, ", f = ", f

    # allocate memory for Jacobi sum
    cdef ZZ_c   ***sum
    sum = <ZZ_c***>sage_malloc(sizeof(ZZ_c**)*d)
    for i in range(d):
        sum[i] = <ZZ_c**>sage_malloc(sizeof(ZZ_c*)*d)
        for j in range(d):
            sum[i][j] = ZZ_new()

    if verbose:
        print "setting up finite field"
    init_NTL_ff(p, n, 0, 0, 0, 0)

    if verbose:
        print "calling library's jacobi_sum()"
    __jacobi_sum(sum, d)

    # extract sum
    if verbose:
        print "extracting sum"
    cdef Integer tmp
    A = []
    for i in range(d):
        w = []
        for j in range(d):
            tmp = PY_NEW(Integer)
            ZZ_to_mpz(&tmp.value, sum[i][j])
            w.append(tmp)
        A.append(w)

    # free memory we borrowed
    for i in range(d):
        for j in range(d):
            ZZ_delete(sum[i][j])
        sage_free(sum[i])
    sage_free(sum)

    return A

def j_curve(K):
    r"""
    Returns the following elliptic curve over K:

        y^2 = x^3 - 108*t/(t-1728)*x + 432*j/(j-1728).

    INPUT:

        - K -- function field F_q(t)

    OUTPUT:

        ellff_EllipticCurve

    EXAMPLES::

        sage: import psage
        sage: K.<t> = psage.FunctionField(GF(11))
        sage: E = psage.ellff.ellff.j_curve(K); E
        <class 'psage.ellff.ellff.ellff_EllipticCurve'>
        sage: E.a4, E.a6
        (2*t^4 + 5*t^3 + 6*t^2 + 9*t, 3*t^6 + 7*t^5 + 8*t^4 + 3*t^3 + 4*t^2 + 8*t)

    """

    t = K.gens()[0]
    return ellff_EllipticCurve(K, [0,0,0,-108*t/(t-1728),432*t/(t-1728)])

