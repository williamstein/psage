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
Fast Cython code needed to compute L-series of elliptic curves over F = Q(sqrt(5)).

USES:

   - The fast residue class rings defined in
     psage.modform.hilbert.sqrt5.sqrt5_fast for naive enumeration.

   - Drew Sutherlands smalljac for point counting

   - Lcalc for evaluating the L-series

   - Dokchitser as well.

   - Computes the *root number* in addition to the L-series.
   
"""

####################################################################
# Straightforward Elliptic Curve Pointcount
####################################################################

from psage.modform.hilbert.sqrt5.sqrt5_fast cimport (
    ResidueRingElement, ResidueRing_abstract, residue_element
    )

from psage.modform.hilbert.sqrt5.sqrt5_fast import ResidueRing

cpdef long ap_via_enumeration(ResidueRingElement a4, ResidueRingElement a6) except -1:
    """
    Compute the trace of Frobenius `a_p` on the elliptic curve defined
    by `y^2 = x^3 + a_4 x + a_6` using a straightforward enumeration
    algorithm.  Here `p` must be a prime of good reduction for the
    equation of the curve.

    EXAMPLES::

        sage: import psage.ellcurve.lseries.sqrt5 as s
        sage: K.<a> = NumberField(x^2-x-1)
        sage: E = EllipticCurve([1,-a,a,5*a-4,-3*a+4])
        sage: _,_,_,a4,a6 = E.short_weierstrass_model().a_invariants()
        sage: P = K.ideal(163)
        sage: import psage.modform.hilbert.sqrt5.sqrt5_fast as t
        sage: R = t.ResidueRing(P, 1)
        sage: s.ap_via_enumeration(R(a4), R(a6))
        120
        sage: k = P.residue_field(); E0 = E.change_ring(k); k.cardinality() + 1 - E0.cardinality()
        120
    """
    cdef long i, cnt = 1  # point at infinity
    cdef ResidueRing_abstract R = a4.parent()

    if R.e != 1:
        raise ValueError, "residue ring must be a field"

    cdef long p = R.p
    cdef long n = p*(p-1)/2
    cdef residue_element x, z, w
    for i in range(R.cardinality()):
        R.unsafe_ith_element(x, i)
        R.mul(z, x, x)   # z = x*x
        R.mul(z, z, x)   # z = x^3
        R.mul(w, a4.x, x)  # w = a4*x
        R.add(z, z, w)     # z = z + w = x^3 + a4*x
        R.add(z, z, a6.x)  # z = x^3 + a4*x + a6
        if R.element_is_0(z):
            cnt += 1
        elif R.is_square(z):
            cnt += 2
    return R.cardinality() + 1 - cnt

cdef extern from "pari/pari.h":
    unsigned long Fl_sqrt(unsigned long a, unsigned long p)

from sage.libs.gmp.mpz cimport (
    mpz_t, mpz_set, mpz_init, mpz_clear,
    mpz_fdiv_ui)

cpdef unsigned long sqrt_mod(unsigned long a, unsigned long p):
    return Fl_sqrt(a, p)

from sage.rings.integer cimport Integer

from sage.all import prime_range, pari

from psage.libs.smalljac.wrapper1 import elliptic_curve_ap

include "stdsage.pxi"
include "interrupt.pxi"

# Since we assume computer is 64-bit, this is fine -- this is
# definitely way bigger than any actual a_p we would ever see.
cdef long UNKNOWN = 2**63  

cdef class TracesOfFrobenius:
    cdef long bound, table_size
    cdef long* primes
    cdef long* sqrt5
    cdef long* ap
    cdef mpz_t Ax, Ay, Bx, By
    cdef object a1, a2, a3, a4, a6

    ##########################################################
    # Allocate and de-allocate basic data structures
    ##########################################################
    
    def __cinit__(self):
        self.primes = NULL
        self.sqrt5 = NULL
        self.ap = NULL
        mpz_init(self.Ax); mpz_init(self.Ay); mpz_init(self.Bx); mpz_init(self.By)
        
    def __init__(self, a1, a2, a3, a4, a6, long bound):
        # Make table of primes up to the bound, and uninitialized
        # corresponding a_p (traces of Frobenius)
        self._initialize_coefficients(a1,a2,a3,a4,a6)
        self._initialize_prime_ap_tables(bound)
        #self._compute_split_traces()

    def _initialize_prime_ap_tables(self, long bound):
        self.bound = bound
        cdef long max_size = 2*bound+20
        self.primes = <long*> sage_malloc(sizeof(long)*max_size)
        self.sqrt5 = <long*> sage_malloc(sizeof(long)*max_size)
        self.ap = <long*> sage_malloc(sizeof(long)*max_size)

        cdef long p, t, i=0, sr0, sr1
        for p in prime_range(bound):
            if i >= max_size:
                raise RuntimeError, "memory assumption violated"
            t = p % 5
            if t == 1 or t == 4:
                # split case
                self.primes[i] = p
                self.primes[i+1] = p
                self.ap[i] = UNKNOWN
                self.ap[i+1] = UNKNOWN
                sr0 = sqrt_mod(5, p)
                sr1 = p - sr0
                if sr0 > sr1: # swap
                    sr0, sr1 = sr1, sr0
                self.sqrt5[i] = sr0
                self.sqrt5[i+1] = sr1
                i += 2
            else:
                # inert or ramified cases (ignore primes with too big norm)
                if p == 5 or p*p < bound:
                    self.primes[i] = p
                    self.sqrt5[i] = 0
                    self.ap[i] = UNKNOWN
                    i += 1
                
        self.table_size = i
    
    def _initialize_coefficients(self,a1,a2,a3,a4,a6):
        # Store all coefficients -- may be needed for trace of frob mod 2 and 3 and bad primes.
        self.a1 = a1; self.a2 = a2; self.a3 = a3; self.a4 = a4; self.a6 = a6

        # Compute short Weierstrass form directly (for speed purposes)
        if a1 or a2 or a3:
            b2 = a1*a1 + 4*a2; b4 = a1*a3 + 2*a4; b6 = a3**2 + 4*a6
            if b2:
                A = 8*b4; B = 16*b6
            else:
                c4 = b2**2 - 24*b4; c6 = -b2**3 + 36*b2*b4 - 216*b6
                A = -27*c4; B = -54*c6
        else:
            A = a4; B = a6

        # Now do the change of variables y|-->y/8; x|-->x/4 to get an isomorphic
        # curve so that A and B are in Z[sqrt(5)], and the reduction isn't
        # any worse at 2.
        A *= 16
        B *= 64

        # Store the short Weierstrass form coefficients as a 4-tuple of GMP ints, 
        # where (Ax,Ay) <--> Ax + sqrt(5)*Ay.
        # By far the best way to get these Ax, Ay are to use the parts
        # method, which is very fast, and gives precisely what we want
        # (for elements of a quadratic field).  As a bonus, at least
        # presently, non-quadratic number field elements don't have a
        # parts method, which is a safety check.
        v = A.parts()
        cdef Integer z
        z = Integer(v[0]); mpz_set(self.Ax, z.value)
        z = Integer(v[1]); mpz_set(self.Ay, z.value)
        v = B.parts()
        z = Integer(v[0]); mpz_set(self.Bx, z.value)
        z = Integer(v[1]); mpz_set(self.By, z.value)

    def __dealloc__(self):
        mpz_clear(self.Ax); mpz_clear(self.Ay); mpz_clear(self.Bx); mpz_clear(self.By)
        if self.primes:
            sage_free(self.primes)
        if self.sqrt5:
            sage_free(self.sqrt5)
        if self.ap:
            sage_free(self.ap)
        
    def _tables(self):
        """
        Return Python dictionary with copies of the internal tables.
        This is used mainly for debugging.

        OUTPUT:
            - dictionary {'primes':primes, 'sqrt5':sqrt5, 'ap':ap}
        """
        cdef long i
        primes = [self.primes[i] for i in range(self.table_size)]
        sqrt5 = [self.sqrt5[i] for i in range(self.table_size)]
        ap = [self.ap[i] if self.ap[i] != UNKNOWN else None for i in range(self.table_size)]
        return {'primes':primes, 'sqrt5':sqrt5, 'ap':ap}

    def aplist(self):
        cdef long i
        v = []
        for i in range(self.table_size):
            v.append( (self.primes[i], self.sqrt5[i], self.ap[i] if self.ap[i] != UNKNOWN else None) )
        return v

    ##########################################################
    # Compute traces
    ##########################################################
    def _compute_split_traces(self, algorithm='smalljac'):
        cdef long i, p, a, b
        for i in range(self.table_size):
            if self.sqrt5[i] != 0:
                p = self.primes[i]
                # 1. Reduce A, B modulo the prime to get a curve over GF(p)
                a = mpz_fdiv_ui(self.Ax, p) + mpz_fdiv_ui(self.Ay, p)*self.sqrt5[i]
                b = mpz_fdiv_ui(self.Bx, p) + mpz_fdiv_ui(self.By, p)*self.sqrt5[i]

                # 2. Check that the resulting curve is nonsingular
                if (-64*((a*a) % p) *a - 432*b*b)%p == 0:
                    # curve is singular at p -- leave this to some other algorithm
                    #print "p = %s singular"%p
                    continue

                # 3. Compute the trace using smalljac (or PARI?)
                if algorithm == 'smalljac':
                    self.ap[i] = elliptic_curve_ap(a, b, p)
                elif algorithm == 'pari':
                    self.ap[i] = pari('ellap(ellinit([0,0,0,%s,%s],1),%s)'%(a,b,p))   # the 1 in the elliptic curve constructor is super important!
                else:
                    raise ValueError
                

###############################################################
# Specialized code for computing traces modulo the inert primes
###############################################################

from sage.stats.intlist cimport IntList

cdef class InertTraceCalculator:
    cdef public dict tables
    
    def __init__(self):
        self.tables = {}

    def init_table(self, int p):
        assert p >= 7 and (p%5 == 2 or p%5 == 3)  # inert prime
        if self.tables.has_key(p):
            return
        # create the table for the given prime p.
        from psage.modform.hilbert.sqrt5.sqrt5 import F
        cdef ResidueRing_abstract R = ResidueRing(F.ideal(p), 1)

        cdef IntList ap, c_quo
        ap = IntList(R.cardinality())
        c_quo = IntList(R.cardinality())
        self.tables[p] = [R, ap, c_quo]
        
        cdef long i
        cdef residue_element j, a4, a6
        for i in range(R.cardinality()):
            R.unsafe_ith_element(j, i)
            self.elliptic_curve_from_j(a4, a6, j, R)
            self.ap_via_enumeration(&ap._values[i], &c_quo._values[i], a4, a6, R)

    cdef int elliptic_curve_from_j(self, residue_element a4, residue_element a6,
                                   residue_element j,
                                   ResidueRing_abstract R) except -1:
        cdef residue_element k, m_three, m_two

        # k = 1728
        k[0] = 1728 % R.p; k[1] = 0

        if R.element_is_0(j):   # if j==0
            R.set_element_to_0(a4)
            R.set_element_to_1(a6)
            return 0

        if R.cmp_element(j, k) == 0:   # if j==1728
            R.set_element_to_1(a4)
            R.set_element_to_0(a6)
            return 0

        # -3 and -2               
        m_three[0] = R.p - 3; m_three[1] = 0
        m_two[0] = R.p - 2; m_two[1] = 0

        # k = j-1728        
        R.sub(k, j, k)  

        # a4 = -3*j*k
        R.mul(a4, m_three, j)
        R.mul(a4, a4, k)

        # a6 = -2*j*k^2
        R.mul(a6, m_two, j)   
        R.mul(a6, a6, k)
        R.mul(a6, a6, k)
        

    cdef int ap_via_enumeration(self, int* ap, int* c_quo,
                                residue_element a4, residue_element a6,
                                ResidueRing_abstract R) except -1:
        assert R.p >= 7
        cdef long i, cnt = 1  # start 1 because of point at infinity
        cdef residue_element x, z, w
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
        ap[0] = R.cardinality() + 1 - cnt

        # Now compute c4/c6 = a4/(18*a6)
        if R.element_is_0(a6):
            c_quo[0] = -1   # signifies "infinity"
        else:
            x[0] = 18 % R.p;  x[1] = 0  # x = 18
            R.mul(z, x, a6)             # z = 18*a6
            R.inv(w, z)                 # w = 1/(18*a6)
            R.mul(z, a4, w)             # z = a4/(18*a6)
            c_quo[0] = R.index_of_element(z)
        return 0   # no error occurred
            
            


