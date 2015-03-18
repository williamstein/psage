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
Fast Cython code needed to compute Hilbert modular forms over F = Q(sqrt(5)).

All the code in this file is meant to be highly optimized.
"""

include 'stdsage.pxi'
include 'cdefs.pxi'
include 'python.pxi'

from sage.rings.ideal import is_Ideal
from sage.rings.all import Integers, ZZ, QQ
from sage.matrix.all import MatrixSpace, zero_matrix

from sage.rings.integer cimport Integer

from sqrt5 import F

cdef Integer temp1 = Integer(0)

cdef long SQRT_MAX_LONG = 2**(4*sizeof(long)-1)

#from sqrt5_fast_python import ResidueRing

def is_ideal_in_F(I):
    return is_Ideal(I) and I.number_field() is F

cpdef long ideal_characteristic(I):
    return I.number_field()._pari_().idealtwoelt(I._pari_())[0]

###########################################################################
# Logging
###########################################################################
GLOBAL_VERBOSE=False
def global_verbose(s):
    if GLOBAL_VERBOSE: print s

###########################################################################
# Rings
###########################################################################

def ResidueRing(P, unsigned int e):
    """
    Create the residue class ring O_F/P^e, where P is a prime ideal of the ring
    O_F of integers of Q(sqrt(5)).

    INPUT:
        - P -- prime ideal
        - e -- positive integer

    OUTPUT:
        - a fast residue class ring

    """
    if e <= 0:
        raise ValueError, "exponent e must be positive"
    cdef long p = ideal_characteristic(P)
    if p**e >= SQRT_MAX_LONG:
        raise ValueError, "residue field must have size less than %s"%SQRT_MAX_LONG
    if p == 5:   # ramified
        if e % 2 == 0:
            R = ResidueRing_nonsplit(P, p, e/2)
        else:
            R = ResidueRing_ramified_odd(P, p, e)
    elif p%5 in [2,3]:  # nonsplit
        R = ResidueRing_nonsplit(P, p, e)
    else: # split
        R = ResidueRing_split(P, p, e)
    return R

class ResidueRingIterator:
    def __init__(self, R):
        self.R = R
        self.i = 0

    def __iter__(self):
        return self

    def next(self):
        if self.i < self.R.cardinality():
            self.i += 1
            return self.R[self.i-1]
        else:
            raise StopIteration

cdef class ResidueRing_abstract(CommutativeRing):
    def __init__(self, P, p, e):
        """
        INPUT:

            - P -- prime ideal of O_F
            - e -- positive integer
        """
        self.P = P
        self.e = e
        self.p = p
        self.F = P.number_field()

    # Do not just change the definition of compare.
    # It is used by the ResidueRing_ModN code to ensure that the prime
    # over 2 appears first. If you totally changed the ordering, that
    # would suddenly break.
    def __richcmp__(ResidueRing_abstract left, right, int op):
        if left is right:
            # Make case of easy equality fast
            return _rich_to_bool(op, 0)
        if not isinstance(right, ResidueRing_abstract):
            return _rich_to_bool(op, cmp(ResidueRing_abstract, type(right)))
        cdef ResidueRing_abstract r = right
        # slow
        return _rich_to_bool(op,
                   cmp((left.p,left.e,left.P), (r.p,r.e,r.P)))

    def __hash__(self):
        return hash((self.p,self.e,self.P))

    def residue_field(self):
        if self._residue_field is None:
            self._residue_field = self.P.residue_field()
        return self._residue_field

    cdef object element_to_residue_field(self, residue_element x):
        raise NotImplementedError

    def _moduli(self):
        # useful for testing purposes
        return self.n0, self.n1

    def __call__(self, x):
        cdef NumberFieldElement_quadratic y
        cdef ResidueRingElement z
        try:
            y = x
            z = self.new_element()
            self.coerce_from_nf(z.x, y)
            return z
        except TypeError:
            pass

        if hasattr(x, 'parent'):
            P = x.parent()
        else:
            P = None
        if P is self:
            return x
        elif P is not self.F:
            x = self.F(x)
        return self(x)

    cdef ResidueRingElement new_element(self):
        raise NotImplementedError

    cdef int coefficients(self, long* v0, long* v1, NumberFieldElement_quadratic x) except -1:
        # The attributes of x are x.a, x.b and x.denom, defined by
        #         x = (x.a + x.b*sqrt(5))/x.denom.
        # Let n be the characteristic of self, and assume n is odd. Set
        #       a = x.a (mod n),   b = x.b (mod n),   e = x.denom^(-1) (mod n).
        # all long ints.
        #   We have x = (a+sqrt(5)*b)*e (mod n).
        # The r[0] and r[1] we want are the coefficients of 1 and (1+sqrt(5))/2.
        #    c + d*(1+sqrt(5))/2 = a*e + sqrt(5)*b*e
        #    c + d/2 + (d/2)*sqrt(5) = a*e + sqrt(5)*b*e
        # So
        #    *v0 = c = a*e - d/2 = a*e - b*e
        #    *v1 = d = 2*b*e
        #
        # The above works fine when n is not a power of 2.  When n is a power of 2,
        # it fails, since x.denom is often divisible by 2.  Write
        #         x = (x.a + x.b*sqrt(5))/(2*f), so 2*f=x.denom.
        # If x is 2-integral, then f is not divisible by 2.  We have
        #     e = 1/x.denom = 1/(2*f).
        # Let f' = 1/f (mod n).
        # We have from the algebra above that in case x.denom % 2 == 0:
        #    *v0 = c =(a-b)*e = (a-b)/(2*f) = ((a-b)/2)/f = (a-b)/2 * f'
        #    *v1 = 2*e*b = b/f = b*f'
        # Thus an integrality condition is that a=b(mod 2).
        # If x.denom % 2 != 0, we just proceed as before.

        cdef long a, b, e, t, f, n
        n = self.n0 if (self.n0 > self.n1) else self.n1
        cdef int is_even = n%2==0
        e = mpz_mod_ui(temp1.value, x.denom, n)
        if e == 0 or (is_even and e%2==0):
            if is_even:
                a = mpz_mod_ui(temp1.value, x.a, 2*n)
                b = mpz_mod_ui(temp1.value, x.b, 2*n)
                # Special 2-power case, as described in comment above.
                f = mpz_mod_ui(temp1.value, x.denom, 2*n) / 2
                if f == 0:
                    raise ZeroDivisionError, "x = %s"%x
                else:
                    t = a - b
                    if t%2 != 0:
                        raise ZeroDivisionError, "x = %s"%x
                    else:
                        if t < 0:
                            t += 2*n
                        if f != 1:
                            f = invmod_long(f, n)
                        v0[0] = ((t/2) * f)%n
                        v1[0] = (b * f)%n
                        return 0 # success
            else:
                raise ZeroDivisionError, "x = %s"%x

        a = mpz_mod_ui(temp1.value, x.a, n)  # all are guaranteed >= 0 by GMP docs
        b = mpz_mod_ui(temp1.value, x.b, n)
        if e != 1:
            e = invmod_long(e, n)
        v0[0] = (a*e)%n - (b*e)%n
        if v0[0] < 0:
            v0[0] += n
        v1[0] = (2*b*e)%n


        return 0  # success

    cdef int coerce_from_nf(self, residue_element r, NumberFieldElement_quadratic x) except -1:
        raise NotImplementedError

    def __repr__(self):
        return "Residue class ring of %s^%s of characteristic %s"%(
            self.P._repr_short(), self.e, self.p)

    def lift(self, x): # for consistency with residue fields
        if x.parent() is self:
            return x.lift()
        raise TypeError

    cdef void add(self, residue_element rop, residue_element op0, residue_element op1):
        raise NotImplementedError
    cdef void sub(self, residue_element rop, residue_element op0, residue_element op1):
        raise NotImplementedError
    cdef void mul(self, residue_element rop, residue_element op0, residue_element op1):
        raise NotImplementedError
    cdef int inv(self, residue_element rop, residue_element op) except -1:
        raise NotImplementedError
    cdef bint is_unit(self, residue_element op):
        raise NotImplementedError
    cdef void neg(self, residue_element rop, residue_element op):
        raise NotImplementedError

    cdef bint element_is_1(self, residue_element op):
        return op[0] == 1 and op[1] == 0

    cdef bint element_is_0(self, residue_element op):
        return op[0] == 0 and op[1] == 0

    cdef void set_element_to_1(self, residue_element op):
        op[0] = 1
        op[1] = 0

    cdef void set_element_to_0(self, residue_element op):
        op[0] = 0
        op[1] = 0

    cdef void set_element(self, residue_element rop, residue_element op):
        rop[0] = op[0]
        rop[1] = op[1]

    cdef int set_element_from_tuple(self, residue_element rop, x) except -1:
        rop[0] = x[0]%self.n0; rop[1] = x[1]%self.n1
        return 0

    cdef int cmp_element(self, residue_element left, residue_element right):
        if left[0] < right[0]:
            return -1
        elif left[0] > right[0]:
            return 1
        elif left[1] < right[1]:
            return -1
        elif left[1] > right[1]:
            return 1
        else:
            return 0

    cdef int pow(self, residue_element rop, residue_element op, long e) except -1:
        cdef residue_element op2
        if e < 0:
            self.inv(op2, op)
            return self.pow(rop, op2, -e)
        self.set_element_to_1(rop)
        cdef residue_element z
        self.set_element(z, op)
        while e:
            if e & 1:
                self.mul(rop, rop, z)
            e /= 2
            if e:
                self.mul(z, z, z)

    cdef bint is_square(self, residue_element op):
        raise NotImplementedError
    cdef int sqrt(self, residue_element rop, residue_element op) except -1:
        raise NotImplementedError

    def __getitem__(self, i):
        cdef ResidueRingElement z = PY_NEW(self.element_class)
        z._parent = self
        self.ith_element(z.x, i)
        return z

    def __iter__(self):
        return ResidueRingIterator(self)

    def is_finite(self):
        return True

    cdef int ith_element(self, residue_element rop, long i) except -1:
        if i < 0 or i >= self.cardinality(): raise IndexError
        self.unsafe_ith_element(rop, i)

    cpdef long cardinality(self):
        return self._cardinality

    cdef void unsafe_ith_element(self, residue_element rop, long i):
        # This assumes 0 <=i < self._cardinality.
        pass # no point in raising exception...

    cdef int next_element(self, residue_element rop, residue_element op) except -1:
        # Raises ValueError if there is no next element.
        # op is assumed valid.
        raise NotImplementedError

    cdef bint is_last_element(self, residue_element op):
        raise NotImplementedError

    cdef long index_of_element(self, residue_element op) except -1:
        # Return the index of the given residue element.
        raise NotImplementedError

    cdef long index_of_element_in_P(self, residue_element op) except -1:
        # Return the index of the given residue element, which is
        # assumed to be in the maximal ideal P.  This is the index, as
        # an element of P.
        raise NotImplementedError

    cdef int next_element_in_P(self, residue_element rop, residue_element op) except -1:
        # Sets rop to next element in the enumeration of self that is in P, assuming
        # that op is in P.
        # Raises ValueError if there is no next element.
        # op is assumed valid (i.e., is in P, etc.).
        raise NotImplementedError

    cdef bint is_last_element_in_P(self, residue_element op):
        raise NotImplementedError

    cdef element_to_str(self, residue_element op):
        cdef ResidueRingElement z = PY_NEW(self.element_class)
        z._parent = self
        z.x[0] = op[0]
        z.x[1] = op[1]
        return str(z)


cdef class ResidueRing_split(ResidueRing_abstract):
    cdef ResidueRingElement new_element(self):
        cdef ResidueRingElement_split z = PY_NEW(ResidueRingElement_split)
        z._parent = self
        return z

    def __init__(self, P, p, e):
        ResidueRing_abstract.__init__(self, P, p, e)
        self.element_class = ResidueRingElement_split
        self.n0 = p**e
        self.n1 = 1
        self._cardinality = self.n0
        # Compute the mapping by finding roots mod p^e.
        alpha = P.number_field().gen()
        cdef long z, z0, z1
        z = sqrtmod_long(5, p, e)
        z0 = divmod_long(1+z, 2, p**e)  # z = (1+sqrt(5))/2 mod p^e.
        if alpha - z0 in P:
            self.im_gen0 = z0
        else:
            z1 = divmod_long(1-z, 2, p**e)  # z = (1-sqrt(5))/2 mod p^e.
            if alpha - z1 in P:
                self.im_gen0 = z1
            else:
                raise RuntimeError, "bug -- maybe F is misdefined to have wrong gen"

    cdef int coerce_from_nf(self, residue_element r, NumberFieldElement_quadratic x) except -1:
        r[1] = 0
        cdef long v0, v1
        self.coefficients(&v0, &v1, x)
        r[0] = v0 + (self.im_gen0*v1)%self.n0
        if r[0] >= self.n0:
            r[0] -= self.n0
        return 0 # success

    def __reduce__(self):
        return ResidueRing_split, (self.P, self.p, self.e)

    cdef object element_to_residue_field(self, residue_element x):
        return self.residue_field()(x[0])

    cdef void add(self, residue_element rop, residue_element op0, residue_element op1):
        rop[0] = (op0[0] + op1[0])%self.n0
        rop[1] = 0

    cdef void sub(self, residue_element rop, residue_element op0, residue_element op1):
        rop[0] = (op0[0] - op1[0])%self.n0
        rop[1] = 0

    cdef void mul(self, residue_element rop, residue_element op0, residue_element op1):
        rop[0] = (op0[0] * op1[0])%self.n0
        rop[1] = 0

    cdef int inv(self, residue_element rop, residue_element op) except -1:
        rop[0] = invmod_long(op[0], self.n0)
        rop[1] = 0
        return 0

    cdef bint is_unit(self, residue_element op):
        return gcd_long(op[0], self.n0) == 1

    cdef void neg(self, residue_element rop, residue_element op):
        rop[0] = self.n0 - op[0] if op[0] else 0
        rop[1] = 0

    cdef bint is_square(self, residue_element op):
        cdef residue_element rop
        if op[0] % self.p == 0:
            # TODO: This is inefficient.
            try:
                # do more complicated stuff involving valuation, unit part, etc.
                self.sqrt(rop, op)
            except ValueError:
                return False
            return True
        # p is odd and op[0] is a unit, so we just check if op[0] is a
        # square modulo self.p
        return kronecker_symbol_long(op[0], self.p) != -1

    cdef int sqrt(self, residue_element rop, residue_element op) except -1:
        rop[0] = sqrtmod_long(op[0], self.p, self.e)
        rop[1] = 0
        if (rop[0] * rop[0])%self.n0 != op[0]:
            raise ValueError, "element must be a square"
        return 0

    cdef void unsafe_ith_element(self, residue_element rop, long i):
        rop[0] = i
        rop[1] = 0

    cdef int next_element(self, residue_element rop, residue_element op) except -1:
        rop[0] = op[0] + 1
        rop[1] = 0
        if rop[0] >= self.n0:
            raise ValueError, "no next element"

    cdef bint is_last_element(self, residue_element op):
        return op[0] == self.n0 - 1

    cdef long index_of_element(self, residue_element op) except -1:
        # Return the index of the given residue element.
        return op[0]

    cdef int next_element_in_P(self, residue_element rop, residue_element op) except -1:
        rop[0] = op[0] + self.p
        rop[1] = 0
        if rop[0] >= self.n0:
            raise ValueError, "no next element in P"
        return 0

    cdef bint is_last_element_in_P(self, residue_element op):
        return op[0] == self.n0 - self.p

    cdef long index_of_element_in_P(self, residue_element op) except -1:
        return op[0]/self.p



cdef class ResidueRing_nonsplit(ResidueRing_abstract):
    cdef ResidueRingElement new_element(self):
        cdef ResidueRingElement_nonsplit z = PY_NEW(ResidueRingElement_nonsplit)
        z._parent = self
        return z

    def __init__(self, P, p, e):
        ResidueRing_abstract.__init__(self, P, p, e)
        self.element_class = ResidueRingElement_nonsplit
        self.n0 = p**e
        self.n1 = p**e
        self._cardinality = self.n0 * self.n1

    cdef int coerce_from_nf(self, residue_element r, NumberFieldElement_quadratic x) except -1:
        # As in the __init__ method from ResidueRingElement_nonsplit, we make r
        # by letting v = x._coefficients(), then r[0]=v[0] and r[1]=v[1],
        # reduced mod self.n0 and self.n1.
        self.coefficients(&r[0], &r[1], x)

    def __reduce__(self):
        return ResidueRing_nonsplit, (self.P, self.p, self.e)

    cdef object element_to_residue_field(self, residue_element x):
        k = self.residue_field()
        if self.p != 5:
            return k([x[0], x[1]])

        # OK, now p == 5 and power of prime sqrt(5) is even...
        # We have that self is Z[x]/(5^n, x^2-x-1) --> F_5, where the map sends x to 3, since x^2-x-1=(x+2)^2 (mod 5).
        # Thus a+b*x |--> a+3*b.
        return k(x[0] + 3*x[1])

    cdef void add(self, residue_element rop, residue_element op0, residue_element op1):
        rop[0] = (op0[0] + op1[0])%self.n0
        rop[1] = (op0[1] + op1[1])%self.n1

    cdef void sub(self, residue_element rop, residue_element op0, residue_element op1):
        # cython uses python mod normalization convention, so no need to do that manually
        rop[0] = (op0[0] - op1[0])%self.n0
        rop[1] = (op0[1] - op1[1])%self.n1

    cdef void mul(self, residue_element rop, residue_element op0, residue_element op1):
        # it is significantly faster doing these arrays accesses only once in the line below!
        cdef long a = op0[0], b = op0[1], c = op1[0], d = op1[1]
        rop[0] = (a*c + b*d)%self.n0
        rop[1] = (a*d + b*d + b*c)%self.n1

    cdef int inv(self, residue_element rop, residue_element op) except -1:
        cdef long a = op[0], b = op[1], w
        w = invmod_long(a*a + a*b - b*b, self.n0)
        rop[0] = (w*(a+b)) % self.n0
        rop[1] = (-b*w) % self.n0
        return 0

    cdef bint is_unit(self, residue_element op):
        cdef long a = op[0], b = op[1], w
        return gcd_long(a*a + a*b - b*b, self.n0) == 1

    cdef void neg(self, residue_element rop, residue_element op):
        rop[0] = self.n0 - op[0] if op[0] else 0
        rop[1] = self.n1 - op[1] if op[1] else 0

    cdef bint is_square(self, residue_element op) except -2:
        """
        Return True only if op is a perfect square.

        ALGORITHM:
        We view the residue ring as R[x]/(p^e, x^2-x-1).
        We first compute the power of p that divides our
        element a+bx, which is just v=min(ord_p(a),ord_p(b)).
        If v=infinity, i.e., a+bx=0, then it's a square.
        If this power is odd, then a+bx can't be a square.
        Otherwise, it is even, and we next look at the unit
        part of a+bx = p^v*unit.  Then a+bx is a square if and
        only if the unit part is a square modulo p^(e-v) >= p,
        since a+bx!=0.
        When p>2, the unit part is a square if and only if it is a
        square modulo p, because of Hensel's lemma, and we can check
        whether or not something is a square modulo p quickly by
        raising to the power (p^2-1)/2, since
        "modulo p", means in the ring O_F/p = GF(p^2), since p
        is an inert prime (this is why we can't use the Legendre
        symbol).
        When p=2, it gets complicated to decide if the unit part is a
        square using code, so just call the sqrt function and catch
        the exception.
        """
        cdef residue_element unit, z
        cdef int v = 0
        cdef long p = self.p

        if op[0] == 0 and op[1] == 0:
            return True

        if p == 5:
            try:
                self.sqrt(z, op)
            except ValueError:
                return False
            return True

        # The "unit" stuff below doesn't make sense for p=5, since
        # then the maximal ideal is generated by sqrt(5) rather than 5.

        unit[0] = op[0]; unit[1] = op[1]
        # valuation at p
        while unit[0]%p ==0 and unit[1]%p==0:
            unit[0] = unit[0]/p; unit[1] = unit[1]/p
            v += 1
        if v%2 != 0:
            return False  # can't be a square

        if p == 2:
            if self.e == 1:
                return True  # every element is a square, since squaring is an automorphism
            try:
                self.sqrt(z, op)
            except ValueError:
                return False
            return True

        # Now unit is not a multiple of p, so it is a unit.
        self.pow(z, unit, (self.p*self.p-1)/2)

        # The result z is -1 or +1, so z[1]==0 (mod p), no matter what, so no need to look at it.
        return z[0]%p == 1

    cdef int sqrt(self, residue_element rop, residue_element op) except -1:
        """
        Set rop to a square root of op, if one exists.

        ALGORITHM:
        We view the residue ring as R[x]/(p^e, x^2-x-1), and want to compute
        a square root a+bx of c+dx.  We have (a+bx)^2=c+dx, so
           a^2+b^2=c and 2ab+b^2=d.
        Setting B=b^2 and doing algebra, we get
           5*B^2 - (4c+2d)B + d^2 = 0,
        so B = ((4c+2d) +/- sqrt((4c+2d)^2-20d^2))/10.
        In characteristic 2 or 5, we compute the numerator mod p^(e+1), then
        divide out by p first, then 10/p.
        """
        if op[0] == 0 and op[1] == 0:
            # easy special case that is hard to deal with using general algorithm.
            rop[0] = 0; rop[1] = 0
            return 0

        # TODO: The following is stupid, and doesn't scale up well... except
        # that if we're computing with a given level, we will do many
        # other things that have at least this bad complexity, so this
        # isn't actually going to slow down anything by more than a
        # factor of 2.   Fix later, when test suite fully in place.
        # I decided to just do this stupidly, since I spent _hours_ trying
        # to get a very fast implementation right and failing...  Yuck.
        cdef long a, b
        cdef residue_element t
        for a in range(self.n0):
            rop[0] = a
            for b in range(self.n1):
                rop[1] = b
                self.mul(t, rop, rop)
                if t[0] == op[0] and t[1] == op[1]:
                    return 0
        raise ValueError, "not a square"

##         cdef long m, a, b, B, c=op[0], d=op[1], p=self.p, n0=self.n0, n1=self.n1, s, t, e=self.e
##         if p == 2 or p == 5:
##             # TODO: This is slow, but painful to implement directly.
##             # Find p-adic roots B of 5*B^2 - (4c+2d)B + d^2 = 0 to at least precision e+1.
##             # One of the roots should be the B we seek.
##             f = ZZ['B']([d*d, -(4*c+2*d), 5])
##             t = 4*c+2*d
##             for B in range(n0):
##                 if (5*B*B - t*B + d*d)%n0 == 0:
##                 #for F, exp in f.factor_padic(p, e+1):
##                 #if F.degree() == 1:
##                 #B = (-F[0]).lift()
##                     print "B=", B
##                     try:
##                         rop[1] = sqrtmod_long(B, p, e)
##                         rop[0] = sqrtmod_long((c-B)%n0, p, e)
##                         print "try: %s,%s"%(rop[0],rop[1])
##                         self.mul(z, rop, rop)
##                         if z[0] == op[0] and z[1] == op[1]:
##                             return 0
##                         else: print "fail"
##                         # try other sign for "b":
##                         rop[1] = n0 - rop[1]
##                         print "try: %s,%s"%(rop[0],rop[1])
##                         self.mul(z, rop, rop)
##                         if z[0] == op[0] and z[1] == op[1]:
##                             return 0
##                         else: print "fail"
##                     except ValueError:
##                         pass
##             raise ValueError

##         # general case (p!=2,5):
##         t = 4*c + 2*d
##         s = sqrtmod_long(submod_long(mulmod_long(t,t,n0), 20*mulmod_long(d,d,n0), n0), p, e)
##         cdef residue_element z
##         try:
##             B = divmod_long((t+s)%n0, 10, n0)
##             rop[1] = sqrtmod_long(B, p, e)
##             rop[0] = sqrtmod_long((c-B)%n0, p, e)
##             self.mul(z, rop, rop)
##             if z[0] == op[0] and z[1] == op[1]:
##                 return 0
##             # try other sign for "b":
##             rop[1] = n0 - rop[1]
##             self.mul(z, rop, rop)
##             if z[0] == op[0] and z[1] == op[1]:
##                 return 0
##         except ValueError:
##             pass
##         # Try the other choice of sign (other root s).
##         B = divmod_long((t-s)%n0, 10, n0)
##         rop[1] = sqrtmod_long(B, p, e)
##         rop[0] = sqrtmod_long((c-B)%n0, p, e)
##         self.mul(z, rop, rop)
##         if z[0] == op[0] and z[1] == op[1]:
##             return 0
##         rop[1] = n0 - rop[1]
##         self.mul(z, rop, rop)
##         assert z[0] == op[0] and z[1] == op[1]
##         return 0



##     cdef int sqrt(self, residue_element rop, residue_element op) except -1:
##         k = self.residue_field()
##         if k.degree() == 1:
##             if self.e == 1:
##                 # happens in ramified case with odd exponent
##                 rop[0] = sqrtmod_long(op[0], self.p, 1)
##                 rop[1] = 0
##                 return 0
##             raise NotImplementedError, 'sqrt not implemented in non-prime ramified case...'
##
##         # TODO: the stupid overhead in this step alone is vastly more than
##         # the time to actually compute the sqrt in Givaro or Pari (say)...
##         a = k([op[0],op[1]]).sqrt()._vector_()
##
##         # The rest of this function is fast (hundreds of times
##         # faster than above line)...
##
##         rop[0] = a[0]; rop[1] = a[1]
##         if self.e == 1:
##             return 0
##
##         if rop[0] == 0 and rop[1] == 0:
##             # TODO: Finish sqrt when number is 0 mod P in inert case.
##             raise NotImplementedError
##
##         # Hensel lifting to get square root modulo P^e.
##         # See the code sqrtmod_long below, which is basically
##         # the same, but more readable.
##         cdef int m
##         cdef long pm, ppm, p
##         p = self.p; pm = p
##         # By "x" below, we mean the input op.
##         # By "y" we mean the square root just computed above, currently rop.
##         cdef residue_element d, y_squared, y2, t, inv
##         for m in range(1, self.e):
##             ppm = p*pm
##             # We compute ((x-y^2)/p^m) / (2*y)
##             self.mul(y_squared, rop, rop)  # y2 = y^2
##             self.sub(t, op, y_squared)
##             assert t[0]%pm == 0  # TODO: remove this after some tests
##             assert t[1]%pm == 0
##             t[0] /= pm; t[1] /= pm
##             # Now t = (x-y^2)/p^m.
##             y2[0] = (2*rop[0])%ppm;   y2[1] = (2*rop[1])%ppm
##             self.inv(inv, y2)
##             self.mul(t, t, inv)
##             # Now t = ((x-y^2)/p^m) / (2*y)
##             t[0] *= pm;  t[1] *= pm
##             # Now do y = y + pm * t
##             self.add(rop, rop, t)
##
##         return 0 # success

    cdef void unsafe_ith_element(self, residue_element rop, long i):
        # Definition: If the element is rop = [a,b], then
        # we have i = a + b*self.n0
        rop[0] = i%self.n0
        rop[1] = (i - rop[0])/self.n0

    cdef int next_element(self, residue_element rop, residue_element op) except -1:
        rop[0] = op[0] + 1
        rop[1] = op[1]
        if rop[0] == self.n0:
            rop[0] = 0
            rop[1] += 1
            if rop[1] >= self.n1:
                raise ValueError, "no next element"

    cdef bint is_last_element(self, residue_element op):
        return op[0] == self.n0 - 1 and op[1] == self.n1 - 1

    cdef long index_of_element(self, residue_element op) except -1:
        # Return the index of the given residue element.
        return op[0] + op[1]*self.n0

    cdef int next_element_in_P(self, residue_element rop, residue_element op) except -1:
        """

        For p = 5, we use the following enumeration, for R=O_F/P^(2e).  First, note that we
        represent R as
            R = (Z/5^e)[x]/(x^2-x-1)
        The maximal ideal is the kernel of the natural map to Z/5Z sending x to 3, i.e., it
        has kernel the principal ideal (x+2).
        By considering (a+bx)(x+2) = 2a+b + (a+3b)x, we see that the maximal ideal is
           { a + (3a+5b)x  :  a in Z/5^eZ, b in Z/5^(e-1)Z }.
        We enumerate this by the standard dictionary ordered enumeration
        on (Z/5^eZ) x (Z/5^(e-1)Z).
        """
        if self.p == 5:
            # a = op[0]
            rop[0] = (op[0] + 1)%self.n0
            rop[1] = (op[1] + 3)%self.n1
            if rop[0] == 0: # a wraps around
                # done, i.e., element with a=5^e-1 and b=5^(e-1)-1 last values?
                # We get the formula below from this calculation:
                #    a+b*x = 5^e-1 + (3*a+5*b)*x.  Only the coeff of x matters, which is
                #   3a+5b = 3*5^e - 3 + 5^e - 5 = 4*5^e - 8.
                # And of course self.n0 = 5^e, so we use that below.
                if (rop[1]+5) % self.n1 == 0:
                    raise ValueError, "no next element in P"
                # not done, so wrap around and increase b by 1
                rop[0] = 0
                rop[1] = (rop[1] + 5)%self.n1
            return 0

        # General case (p!=5):
        rop[0] = op[0] + self.p
        rop[1] = op[1]
        if rop[0] >= self.n0:
            rop[0] = 0
            rop[1] += self.p
            if rop[1] >= self.n1:
                raise ValueError, "no next element in P"
        return 0

    cdef bint is_last_element_in_P(self, residue_element op):
        if self.p == 5:
            # see the comment above in "next_element_in_P".
            if self.e == 1:
                # next element got by incrementing a in enumeration
                return op[0] == self.n0 - 1
            else:
                # next element got by incrementing both a and b by 1 in enumeration,
                # and 3a+5b=8.
                return op[0] == self.n0 - 1 and (op[1]+8)%self.n1 == 0
        # General case (p!=5):
        return op[0] == self.n0 - self.p and op[1] == self.n1 - self.p

    cdef long index_of_element_in_P(self, residue_element op) except -1:
        if self.p == 5:
            # see the comment above in "next_element_in_P".
            # The index is a + 5^e*b, so we have to recover a and b in op=a+(3a+5b)x.
            # We have a = op[0], which is easy.  Then b=(op[1]-3a)/5
            return op[0] + ((op[1] - 3*op[0])/5 % (self.n1/self.p)) * self.n0

        # General case (p!=5):
        return op[0]/self.p + op[1]/self.p * (self.n0/self.p)

cdef class ResidueRing_ramified_odd(ResidueRing_abstract):
    cdef ResidueRingElement new_element(self):
        cdef ResidueRingElement_ramified_odd z = PY_NEW(ResidueRingElement_ramified_odd)
        z._parent = self
        return z

    cdef long two_inv
    def __init__(self, P, p, e):
        ResidueRing_abstract.__init__(self, P, p, e)
        # This is assumed in the code for the element_class so it better be true!
        assert self.F.defining_polynomial().list() == [-1,-1,1]
        self.n0 = p**(e//2 + 1)
        self.n1 = p**(e//2)
        self._cardinality = self.n0 * self.n1
        self.two_inv = (Integers(self.n0)(2)**(-1)).lift()
        self.element_class = ResidueRingElement_ramified_odd

    cdef int coerce_from_nf(self, residue_element r, NumberFieldElement_quadratic x) except -1:
        # see docs for __init__ method of
        #    cdef class ResidueRingElement_ramified_odd(ResidueRingElement)
        cdef long v0, v1
        self.coefficients(&v0, &v1, x)
        if v0 == 0 and v1 == 0:
            r[0] = 0; r[1] = 0
        elif v1 == 0:
            r[0] = v0; r[1] = 0
        else:
            r[0] = (v0 + v1*self.two_inv) % self.n0
            r[1] = (v1*self.two_inv) % self.n1
        return 0 # success

    def __reduce__(self):
        return ResidueRing_ramified_odd, (self.P, self.p, self.e)

    cdef object element_to_residue_field(self, residue_element x):
        k = self.residue_field()
        # For odd ramified case, we use a different basis, which is:
        # x[0]+sqrt(5)*x[1], so mapping to residue field mod (sqrt(5))
        # is easy:
        return k(x[0])

    cdef void add(self, residue_element rop, residue_element op0, residue_element op1):
        rop[0] = (op0[0] + op1[0])%self.n0
        rop[1] = (op0[1] + op1[1])%self.n1

    cdef void sub(self, residue_element rop, residue_element op0, residue_element op1):
        rop[0] = (op0[0] - op1[0])%self.n0
        rop[1] = (op0[1] - op1[1])%self.n1

    cdef void mul(self, residue_element rop, residue_element op0, residue_element op1):
        cdef long a = op0[0], b = op0[1], c = op1[0], d = op1[1]
        rop[0] = (a*c + 5*b*d)%self.n0
        rop[1] = (a*d + b*c)%self.n1

    cdef void neg(self, residue_element rop, residue_element op):
        rop[0] = self.n0 - op[0] if op[0] else 0
        rop[1] = self.n1 - op[1] if op[1] else 0

    cdef int inv(self, residue_element rop, residue_element op) except -1:
        """
        We derive by algebra the inverse of an element of the quotient ring.
        We represent the ring as:
          Z[x]/(x^2-5, x^e) = {a+b*x with 0<=a<p^((e/2)+1), 0<=b<p^(e/2)}
        Assume a+b*x is a unit, so a!=0(mod 5).  Let c+d*x be inverse of a+b*x.
        We have (a+b*x)*(c+d*x)=1, so ac+5db + (ad+bc)*x = 1.  Thus
           ac+5bd=1 (mod p^((e/2)+1)) and ad+bc=0 (mod p^(e/2)).
        Doing algebra (multiply both sides by a^(-1) of first, and subtract
        second equation), we arrive at
           d = (a^(-1)*b)*(5*a^(-1)*b^2 - a)^(-1) (mod p^(e/2))
           c = a^(-1)*(1-5*b*d) (mod p^(e/2+1)).
        Notice that the first equation determines d only modulo p^(e/2), and
        the second only requires d modulo p^(e/2)!
        """
        cdef long a = op[0], b = op[1], ainv = invmod_long(a, self.n0), n0=self.n0, n1=self.n1
        # Set rop[1] to "ainv*b*(5*ainv*b*b - a)".  We have to use ugly notation to avoid overflow
        rop[1] = divmod_long(mulmod_long(ainv,b,n1), submod_long(mulmod_long(mulmod_long(mulmod_long(self.p,ainv,n1),b,n1),b,n1), a, n1), n1)
        # Set rop[0] to "ainv*(1-5*b*d)
        rop[0] = mulmod_long(ainv, submod_long(1,mulmod_long(5,mulmod_long(b,rop[1],n0),n0),n0),n0)
        return 0

    cdef bint is_unit(self, residue_element op):
        """
        Return True if the element op=(a,b) is a unit.
        We represent the ring as
            Z[x]/(x^2-5, x^e) = {a+b*x with 0<=a<p^((e/2)+1), 0<=b<p^(e/2)}
        An element is in the maximal ideal (x) if and only if it is of the form:
             (a+b*x)*x = a*x+b*x*x = a*x + 5*b.
        Since a is arbitrary, this means the maximal ideal is the set of
        elements op=(a,b) with a%5==0, so the units are the elements
        with a%5!=0.
        """
        return op[0]%self.p != 0

    cdef bint is_square(self, residue_element op) except -2:
        cdef residue_element rop
        try:
            self.sqrt(rop, op)
        except:  # ick
            return False
        return True

    cdef int sqrt(self, residue_element rop, residue_element op) except -1:
        """
        Algorithm: Given c+dx in Z[x]/(x^2-5,x^e), we seek a+bx with
                  (a+bx)^2=a^2+5b^2 + 2abx = c+dx,
        so a^2+5b^2 = c (mod p^n0)
           2ab = d (mod p^n1)
        Multiply first equation through by A=a^2, to get quadratic in A:
              A^2 - c*A + 5d^2/4 = 0 (mod n0).
        Solve for A = (c +/- sqrt(c^2-5d^2))/2.

        Thus the algorithm is:
              1. Extract sqrt s of c^2-5d^2 modulo n0.  If no sqrt, then not square
              2. Extract sqrt a of (c+s)/2 or (c-s)/2.  Both are squares.
              3. Let b=d/(2*a) for each choice of a.  Square each and see which works.
        """
        # see comment for sqrt above (in inert case).
        cdef long a, b
        cdef residue_element t
        for a in range(self.n0):
            rop[0] = a
            for b in range(self.n1):
                rop[1] = b
                self.mul(t, rop, rop)
                if t[0] == op[0] and t[1] == op[1]:
                    return 0
        raise ValueError, "not a square"

##         cdef long c=op[0], d=op[1], p=self.p, n0=self.n0, n1=self.n1, s, two_inv
##         s = sqrtmod_long((c*c - p*d*d)%n0, p, self.e//2 + 1)
##         two_inv = invmod_long(2, n0)
##         cdef residue_element t
##         try:
##             rop[0] = sqrtmod_long(((c+s)*two_inv)%n0, p, self.e//2 + 1)
##             rop[1] = (divmod_long(d,rop[0],n1) * two_inv)%n1
##             self.mul(t, rop, rop)
##             if t[0] == op[0] and t[1] == op[1]:
##                 return 0
##         except ValueError:
##             pass
##         rop[0] = sqrtmod_long(((c-s)*two_inv)%n0, p, self.e//2 + 1)
##         rop[1] = (divmod_long(d,rop[0],n1) * two_inv)%n1
##         self.mul(t, rop, rop)
##         assert t[0] == op[0] and t[1] == op[1]
##         return 0


    cdef void unsafe_ith_element(self, residue_element rop, long i):
        # Definition: If the element is rop = [a,b], then
        # we have i = a + b*self.n0
        rop[0] = i%self.n0
        rop[1] = (i - rop[0])/self.n0

    cdef int next_element(self, residue_element rop, residue_element op) except -1:
        rop[0] = op[0] + 1
        rop[1] = op[1]
        if rop[0] == self.n0:
            rop[0] = 0
            rop[1] += 1
            if rop[1] >= self.n1:
                raise ValueError, "no next element"

    cdef bint is_last_element(self, residue_element op):
        return op[0] == self.n0 - 1 and op[1] == self.n1 - 1

    cdef long index_of_element(self, residue_element op) except -1:
        # Return the index of the given residue element.
        return op[0] + op[1]*self.n0

    cdef int next_element_in_P(self, residue_element rop, residue_element op) except -1:
        rop[0] = op[0] + self.p
        if rop[0] == self.n0:
            rop[0] = 0
            rop[1] += 1
            if rop[1] >= self.n1:
                raise ValueError, "no next element"

    cdef bint is_last_element_in_P(self, residue_element op):
        return op[0] == self.n0 - self.p and op[1] == self.n1 - 1

    cdef long index_of_element_in_P(self, residue_element op) except -1:
        return op[0]/self.p + op[1] * (self.n0/self.p)


###########################################################################
# Ring elements
###########################################################################
cdef inline bint _rich_to_bool(int op, int r):  # copied from sage.structure.element...
    if op == Py_LT:  #<
        return (r  < 0)
    elif op == Py_EQ: #==
        return (r == 0)
    elif op == Py_GT: #>
        return (r  > 0)
    elif op == Py_LE: #<=
        return (r <= 0)
    elif op == Py_NE: #!=
        return (r != 0)
    elif op == Py_GE: #>=
        return (r >= 0)


cdef class ResidueRingElement:
    cpdef parent(self):
        return self._parent
    cdef new_element(self):
        raise NotImplementedError

    cpdef long index(self):
        """
        Return the index of this element in the enumeration of
        elements of the parent.
        """
        return self._parent.index_of_element(self.x)

    def __add__(ResidueRingElement left, ResidueRingElement right):
        cdef ResidueRingElement z = left.new_element()
        left._parent.add(z.x, left.x, right.x)
        return z
    def __sub__(ResidueRingElement left, ResidueRingElement right):
        cdef ResidueRingElement z = left.new_element()
        left._parent.sub(z.x, left.x, right.x)
        return z
    def __mul__(ResidueRingElement left, ResidueRingElement right):
        cdef ResidueRingElement z = left.new_element()
        left._parent.mul(z.x, left.x, right.x)
        return z
    def __div__(ResidueRingElement left, ResidueRingElement right):
        cdef ResidueRingElement z = left.new_element()
        left._parent.inv(z.x, right.x)
        left._parent.mul(z.x, left.x, z.x)
        return z
    def __neg__(ResidueRingElement self):
        cdef ResidueRingElement z = self.new_element()
        self._parent.neg(z.x, self.x)
        return z
    def __pow__(ResidueRingElement self, e, m):
        cdef ResidueRingElement z = self.new_element()
        self._parent.pow(z.x, self.x, e)
        return z

    def __invert__(ResidueRingElement self):
        cdef ResidueRingElement z = self.new_element()
        self._parent.inv(z.x, self.x)
        return z

    def __richcmp__(ResidueRingElement left, ResidueRingElement right, int op):
        cdef int c
        if left.x[0] < right.x[0]:
            c = -1
        elif left.x[0] > right.x[0]:
            c = 1
        elif left.x[1] < right.x[1]:
            c = -1
        elif left.x[1] > right.x[1]:
            c = 1
        else:
            c = 0
        return _rich_to_bool(op, c)

    def __hash__(self):
        return self.x[0] + self._parent.n0*self.x[1]

    cpdef bint is_unit(self):
        return self._parent.is_unit(self.x)

    cpdef bint is_square(self):
        return self._parent.is_square(self.x)

    cpdef sqrt(self):
        cdef ResidueRingElement z = self.new_element()
        self._parent.sqrt(z.x, self.x)
        return z

cdef class ResidueRingElement_split(ResidueRingElement):
    def __init__(self, ResidueRing_split parent, x):
        self._parent = parent
        assert x.parent() is parent.F
        v = x._coefficients()
        self.x[1] = 0
        if len(v) == 0:
            self.x[0] = 0
            return
        elif len(v) == 1:
            self.x[0] = v[0] % self._parent.n0
            return
        self.x[0] = v[0] + parent.im_gen0*v[1]
        self.x[0] = self.x[0] % self._parent.n0
        self.x[1] = 0

    cdef new_element(self):
        cdef ResidueRingElement_split z = PY_NEW(ResidueRingElement_split)
        z._parent = self._parent
        return z

    def __repr__(self):
        return str(self.x[0])

    def lift(self):
        """
        Return lift of element to number field F.
        """
        return self._parent.F(self.x[0])

cdef class ResidueRingElement_nonsplit(ResidueRingElement):
    """
    EXAMPLES::

        sage: from psage.modform.hilbert.sqrt5.sqrt5_fast import ResidueRing

    Inert example::

        sage: F.<a> = NumberField(x^2 - x - 1)
        sage: P_inert = F.primes_above(3)[0]
        sage: R_inert = ResidueRing(P_inert, 2)
        sage: s = R_inert(F.0); s
        g
        sage: s + s + s
        3*g
        sage: s*s*s*s
        2 + 3*g
        sage: R_inert(F.0^4)
        2 + 3*g

    Ramified example (with even exponent)::

        sage: P = F.primes_above(5)[0]
        sage: R = ResidueRing(P, 2)
        sage: s = R(F.0); s
        g
        sage: s*s*s*s*s
        3
        sage: R(F.0^5)
        3
        sage: a = (s+s - R(F(1))); a*a
        0
        sage: R = ResidueRing(P, 3)
        sage: t = R(2*F.0-1); t   # reduction of sqrt(5)
        s
        sage: t*t*t
        0
    """
    def __init__(self, ResidueRing_nonsplit parent, x):
        self._parent = parent
        assert x.parent() is parent.F
        v = x._coefficients()
        if len(v) == 0:
            self.x[0] = 0; self.x[1] = 0
            return
        elif len(v) == 1:
            self.x[0] = v[0]
            self.x[1] = 0
        else:
            self.x[0] = v[0]
            self.x[1] = v[1]
        self.x[0] = self.x[0] % self._parent.n0
        self.x[1] = self.x[1] % self._parent.n1

    cdef new_element(self):
        cdef ResidueRingElement_nonsplit z = PY_NEW(ResidueRingElement_nonsplit)
        z._parent = self._parent
        return z

    def __repr__(self):
        if self.x[0]:
            if self.x[1]:
                if self.x[1] == 1:
                    return '%s + g'%self.x[0]
                else:
                    return '%s + %s*g'%(self.x[0], self.x[1])
            return str(self.x[0])
        else:
            if self.x[1]:
                if self.x[1] == 1:
                    return 'g'
                else:
                    return '%s*g'%self.x[1]
            return '0'

    def lift(self):
        """
        Return lift of element to number field F.
        """
        # TODO: test...
        return self._parent.F([self.x[0], self.x[1]])

cdef class ResidueRingElement_ramified_odd(ResidueRingElement):
    """
    Element of residue class ring R = O_F / P^(2f-1), where e=2f-1 is
    odd, and P=sqrt(5)O_F is the ramified prime.

    Computing with this ring is trickier than all the rest,
    since it's not a quotient of Z[x] of the form Z[x]/(m,g),
    where m is an integer and g is a polynomial.
    The ring R is the quotient of
        O_F/P^(2f) = O_F/5^f = (Z/5^fZ)[x]/(x^2-x-1),
    by the ideal x^e.  We have
        O_F/P^(2f-2) subset R subset O_F/P^(2f)
    and each successive quotient has order 5 = #(O_F/P).
    Thus R has cardinality 5^(2f-1) and characteristic 5^f.
    The ring R can't be a quotient of Z[x] of the
    form Z[x]/(m,g), since such a quotient has
    cardinality m^deg(g) and characteristic m, and
    5^(2f-1) is not a power of 5^f.

    We thus view R as

        R = (Z/5^fZ)[x] / (x^2 - 5,  5^(f-1)*x).

    The elements of R are pairs (a,b) in (Z/5^fZ) x (Z/5^(f-1)Z),
    which correspond to the class of a + b*x.  The arithmetic laws
    are thus:

       (a,b) + (c,d) = (a+c mod 5^f, b+d mod 5^(f-1))

    and

       (a,b) * (c,d) = (a*c+b*d*5 mod 5^f, a*d+b*c mod 5^(f-1))

    The element omega = F.gen(), which is (1+sqrt(5))/2 maps to
    (1+x)/2 = (1/2, 1/2), which is the generator of this ring.
    """
    def __init__(self, ResidueRing_ramified_odd parent, x):
        self._parent = parent
        assert x.parent() is parent.F
        # We can assume that the defining poly of F is T^2-T-1 (this is asserted
        # in some code above), so the gen is (1+sqrt(5))/2.  We can also
        # assume x has no denom. Then sqrt(5)=2*x-1, so need to find c,d such that:
        #   a + b*x = c + d*(2*x-1)
        # ===> c = a + b/2,  d = b/2
        cdef long a, b
        v = x._coefficients()
        if len(v) == 0:
            self.x[0] = 0; self.x[1] = 0
            return
        if len(v) == 1:
            self.x[0] = v[0]
            self.x[1] = 0
        else:
            a = v[0]; b = v[1]
            self.x[0] = a + b*parent.two_inv
            self.x[1] = b*parent.two_inv
        self.x[0] = self.x[0] % self._parent.n0
        self.x[1] = self.x[1] % self._parent.n1

    cdef new_element(self):
        cdef ResidueRingElement_ramified_odd z = PY_NEW(ResidueRingElement_ramified_odd)
        z._parent = self._parent
        return z

    def __repr__(self):
        if self.x[0]:
            if self.x[1]:
                if self.x[1] == 1:
                    return '%s + s'%self.x[0]
                else:
                    return '%s + %s*s'%(self.x[0], self.x[1])
            return str(self.x[0])
        else:
            if self.x[1]:
                if self.x[1] == 1:
                    return 's'
                else:
                    return '%s*s'%self.x[1]
            return '0'


####################################################################
# misc arithmetic needed elsewhere
####################################################################

cdef mpz_t mpz_temp
mpz_init(mpz_temp)
cdef long kronecker_symbol_long(long a, long b):
    mpz_set_si(mpz_temp, b)
    return mpz_si_kronecker(a, mpz_temp)

cdef inline long mulmod_long(long x, long y, long m):
    return (x*y)%m

cdef inline long submod_long(long x, long y, long m):
    return (x-y)%m

cdef inline long addmod_long(long x, long y, long m):
    return (x+y)%m

#cdef long kronecker_symbol_long2(long a, long p):
#    return powmod_long(a, (p-1)/2, p)

cdef long powmod_long(long x, long e, long m):
    # return x^m mod e, where x and e assumed < SQRT_MAX_LONG
    cdef long y = 1, z = x
    while e:
        if e & 1:
            y = (y * z) % m
        e /= 2
        if e: z = (z * z) % m
    return y

# This is kind of ugly...
from sage.rings.fast_arith cimport arith_llong
cdef arith_llong Arith_llong = arith_llong()
cdef long invmod_long(long x, long m) except -1:
    cdef long long g, s, t
    g = Arith_llong.c_xgcd_longlong(x, m, &s, &t)
    if g != 1:
        raise ZeroDivisionError
    return s%m

cdef long gcd_long(long a, long b) except -1:
    return Arith_llong.c_gcd_longlong(a, b)

cdef long divmod_long(long x, long y, long m) except -1:
    #return quotient x/y mod m, assuming x, y < SQRT_MAX_LONG
    return (x * invmod_long(y, m)) % m

cpdef long sqrtmod_long(long x, long p, int n) except -1:
    cdef long a, r, y, z
    if n <= 0:
        raise ValueError, "n must be positive"

    if x == 0 or x == 1:
        # easy special case that works no matter what p is.
        return x

    if p == 2: # pain in the butt case
        if n == 1:
            return x
        elif n == 2:
            # since x isn't 1 or 2 (see above special case)
            raise ValueError, "not a square"
        elif n == 3:
            if x == 4: return 2
            # since x isn't 1 or 2 or 4 (see above special case)
            raise ValueError, "not a square"
        elif n == 4:
            if x == 4: return 2
            if x == 9: return 3
            raise ValueError, "not a square"
        else:
            # a general slow algorithm -- use pari; this won't get called much
            try:
                # making a string as below is still *much* faster than
                # just doing say:  "Zp(2,20)(16).sqrt()".  Why? Because
                # even evaluating 'Zp(2,20)' in sage takes longer than
                # creating and evaluating the whole string with pari!
                # And then the way Zp in Sage computes the square root is via...
                # making the corresponding Pari object.
                from sage.libs.pari.all import pari, PariError
                cmd = "lift(sqrt(%s+O(2^%s)))%%(2^%s)"%(x, 2*n, n)
                return int(pari(cmd))
            except PariError:
                raise ValueError, "not a square"

    if x%p == 0:
        a = x/p
        y = 1
        r = 1
        while r < n and a%p == 0:
            a /= p
            r += 1
            if r%2 == 0:
                y *= p
        if r % 2 != 0:
            raise ValueError, "not a square"
        return y * sqrtmod_long(a, p, n - r//2)

    # Compute square root of x modulo p^n using Hensel lifting, where
    # p id odd.

    # Step 1. Find a square root modulo p. (See section 4.5 of my Elementary Number Theory book.)
    if p % 4 == 3:
        # Easy case when p is 3 mod 4.
        y = powmod_long(x, (p+1)/4, p)
    elif p == 5: # hardcode useful special case
        z = x % p
        if z == 0 or z == 1:
            y = x
        elif z == 2 or z == 3:
            raise ValueError
        elif z == 4:
            y = 2
        else:
            assert False, 'bug'
    else:
        # Harder case -- no known poly time deterministic algorithm.
        # TODO: There is a fast algorithm we could implement here, but
        # I'll just use pari for now and come back to this later. This
        # is on page 88 of my book, esp. since this is not critical
        # to the speed of the whole algorithm.
        # Another reason to redo -- stupid PariError is hard to trap...
        from sage.all import pari
        y = pari(x).Mod(p).sqrt().lift()

    if n == 1: # done
        return y

    # Step 2. Hensel lift to mod p^n.
    cdef int m
    cdef long t, pm, ppm
    pm = p
    for m in range(1, n):
        ppm = p*pm
        # Compute t = ((x-y^2)/p^m) / (2*y)
        t = divmod_long( ((x - y*y)%ppm) / pm, (2*y)%ppm, ppm)
        # And set y = y + pm*t, which gives next correct approximation for y.
        y += (pm * t) % ppm
        pm = ppm

    return y

###########################################################################
# The quotient ring O_F/N
###########################################################################

# The worse case is the inert prime 2 (of norm 4).  The maximum level possible is 2^32.
cdef enum:
    MAX_PRIME_DIVISORS = 16

ctypedef residue_element modn_element[MAX_PRIME_DIVISORS]

# 2 x 2 matrices over O_F/N
ctypedef modn_element modn_matrix[4]


cdef class ResidueRingModN:
    cdef public list residue_rings
    cdef object N
    cdef int r  # number of prime divisors

    def __init__(self, N):
        self.N = N

        # A guarantee we make for other code is that if 2 divides N,
        # then the first factor below will be 2.  Thus we sort the
        # factors by their residue characteristic in
        # order to ensure this.  It is also ** SUPER IMPORTANT ** that
        # the order *not* depend on the exponent, so we next sort by gens_reduced, and finally by exponent.
        v = [(P.smallest_integer(), P.gens_reduced()[0], e, ResidueRing(P, e)) for P, e in N.factor()]
        v.sort()
        self.residue_rings = [x[-1] for x in v]
        self.r = len(self.residue_rings)

    def __repr__(self):
        return "Residue class ring modulo the ideal %s of norm %s"%(
            self.N._repr_short(), self.N.norm())

    cdef int set(self, modn_element rop, x) except -1:
        self.coerce_from_nf(rop, self.N.number_field()(x), 0)

    cdef int coerce_from_nf(self, modn_element rop, NumberFieldElement_quadratic op, int odd_only) except -1:
        # Given an element op in the field F, try to reduce it modulo
        # N and create the corresponding modn element.
        # If the odd_only flag is set to 1, leave the result locally
        # at the primes over 2 undefined.
        cdef int i
        cdef ResidueRing_abstract R
        for i in range(self.r):
            R = self.residue_rings[i]
            if odd_only and R.p == 2:
                continue
            R.coerce_from_nf(rop[i], op)
        return 0 # success

    cdef void set_element(self, modn_element rop, modn_element op):
        cdef int i
        for i in range(self.r):
            rop[i][0] = op[i][0]
            rop[i][1] = op[i][1]

    cdef void set_element_omitting_ith_factor(self, modn_element rop, modn_element op, int i):
        cdef int j
        for j in range(i):
            rop[j][0] = op[j][0]
            rop[j][1] = op[j][1]
        for j in range(i+1, self.r):
            rop[j-1][0] = op[j][0]
            rop[j-1][1] = op[j][1]

    cdef int set_element_reducing_exponent_of_ith_factor(self, modn_element rop, modn_element op, int i) except -1:
        cdef ResidueRing_abstract R
        cdef int j, e
        cdef long p, temp
        R = self.residue_rings[i]
        e = R.e
        p = R.p
        # No question what to do with everything before i-th factor:
        for j in range(i):
            rop[j][0] = op[j][0]
            rop[j][1] = op[j][1]

        if (e == 1 and p!=5) or (p == 5 and e == 1 and isinstance(R, ResidueRing_ramified_odd)):
            # Easy: just delete the i-th factor completely.
            # Now fill in the rest
            for j in range(i+1, self.r):
                rop[j-1][0] = op[j][0]
                rop[j-1][1] = op[j][1]
        elif p != 5:
            # If p != 5 then the prime is unramified, so the representative
            # a+b*x for op just works so long as we reduce it mod R.n0/p and R.n1/p.
            rop[i][0] = op[i][0] % (R.n0//p)
            if R.n1 > 1:  # otherwise this factor doesn't matter (and we're in the split case)
                rop[i][1] = op[i][1] % (R.n1//p)
            else:
                rop[i][1] = 0

            # And fill in the rest
            for j in range(i+1, self.r):
                rop[j][0] = op[j][0]
                rop[j][1] = op[j][1]
        else:
            assert p == 5
            # Now we're in the ramified case, which is the hardest, since we use
            # completely different representations for O_F / (sqrt(5))^e, for e even
            # and for e odd.  Also, note that the e above (i.e., R.e) in this case
            # isn't even the power of P=(sqrt(5)) in the even case.
            # Recall that we represent O_F / (sqrt(5))^e for e even as
            #    (Z/5^(e/2)Z)[x]/(x^2-x-1)
            # and for e odd as
            #    {a + b*sqrt(5)  : a < p^(e//2 + 1),   b < p^(e//2)}

            if isinstance(R, ResidueRing_ramified_odd):
                # Case 1: odd power >= 3  of sqrt(5) :
                # We have a+b*sqrt(5), and we need to write it as
                # c+d*x, where x=(1+sqrt(5))/2, so 2*x-1 = sqrt(5).
                # Thus  a+b*sqrt(5) = a+b*(2*x-1) = a-b + 2*b*x
                #print "odd ramified case"
                #print op[i][0],op[i][1],
                rop[i][0] = op[i][0]-op[i][1]
                rop[i][1] = 2*op[i][1]
                #print " |--->", rop[i][0], rop[i][1]
            else:
                #print "even ramified case"
                # Case 2: even power >= 2 of sqrt(5)
                # We have a+b*x and have to write in terms of sqrt(5), so
                #    a+b*x = a+b*(1+sqrt(5))/2 = a+b/2 + (b/2)*sqrt(5)
                temp = divmod_long(op[i][1], 2, R.n1)
                rop[i][0] = (op[i][0] + temp)%R.n1
                rop[i][1] = temp
            # Fill in the rest.
            for j in range(i+1, self.r):
                rop[j][0] = op[j][0]
                rop[j][1] = op[j][1]
        return 0

    cdef element_to_str(self, modn_element op):
        s = ','.join([(<ResidueRing_abstract>self.residue_rings[i]).element_to_str(op[i]) for
                          i in range(self.r)])
        if self.r > 1:
            s = '(' + s + ')'
        return s

    cdef int cmp_element(self, modn_element op0, modn_element op1):
        cdef int i, c
        cdef ResidueRing_abstract R
        for i in range(self.r):
            R = self.residue_rings[i]
            c = R.cmp_element(op0[i], op1[i])
            if c: return c
        return 0

    cdef void mul(self, modn_element rop, modn_element op0, modn_element op1):
        cdef int i
        cdef ResidueRing_abstract R
        for i in range(self.r):
            R = self.residue_rings[i]
            R.mul(rop[i], op0[i], op1[i])

    cdef int inv(self, modn_element rop, modn_element op) except -1:
        cdef int i
        cdef ResidueRing_abstract R
        for i in range(self.r):
            R = self.residue_rings[i]
            R.inv(rop[i], op[i])
        return 0

    cdef int neg(self, modn_element rop, modn_element op) except -1:
        cdef int i
        cdef ResidueRing_abstract R
        for i in range(self.r):
            R = self.residue_rings[i]
            R.neg(rop[i], op[i])
        return 0

    cdef void add(self, modn_element rop, modn_element op0, modn_element op1):
        cdef int i
        cdef ResidueRing_abstract R
        for i in range(self.r):
            R = self.residue_rings[i]
            R.add(rop[i], op0[i], op1[i])

    cdef void sub(self, modn_element rop, modn_element op0, modn_element op1):
        cdef int i
        cdef ResidueRing_abstract R
        for i in range(self.r):
            R = self.residue_rings[i]
            R.sub(rop[i], op0[i], op1[i])

    cdef bint is_last_element(self, modn_element x):
        cdef int i
        cdef ResidueRing_abstract R
        for i in range(self.r):
            R = self.residue_rings[i]
            if not R.is_last_element(x[i]):
                return False
        return True

    cdef int next_element(self, modn_element rop, modn_element op) except -1:
        cdef int i, done = 0
        cdef ResidueRing_abstract R
        for i in range(self.r):
            R = self.residue_rings[i]
            if done:
                R.set_element(rop[i], op[i])
            else:
                if R.is_last_element(op[i]):
                    R.set_element_to_0(rop[i])
                else:
                    R.next_element(rop[i], op[i])
                    done = True

    cdef bint is_square(self, modn_element op):
        cdef int i
        cdef ResidueRing_abstract R
        for i in range(self.r):
            R = self.residue_rings[i]
            if not R.is_square(op[i]):
                return False
        return True

    cdef int sqrt(self, modn_element rop, modn_element op) except -1:
        cdef int i
        cdef ResidueRing_abstract R
        for i in range(self.r):
            R = self.residue_rings[i]
            R.sqrt(rop[i], op[i])
        return 0

    #######################################
    # modn_matrices
    #######################################
    cdef matrix_to_str(self, modn_matrix A):
        return '[%s, %s; %s, %s]'%tuple([self.element_to_str(A[i]) for i in range(4)])

    cdef bint matrix_is_invertible(self, modn_matrix x):
        # det = x[0]*x[3] - x[1]*x[2], so we return True
        # when x[0]*x[3] != x[1]*x[2]
        cdef modn_element t1, t2
        self.mul(t1, x[0], x[3])
        self.mul(t2, x[1], x[2])
        return self.cmp_element(t1, t2) != 0

    cdef int matrix_inv(self, modn_matrix rop, modn_matrix x) except -1:
        cdef modn_element d, t1, t2
        self.mul(t1, x[0], x[3])
        self.mul(t2, x[1], x[2])
        self.sub(d, t1, t2)   # x[0]*x[2] - x[1]*x[3]
        self.inv(d, d)
        # Inverse is [x3,-x1,-x2,x0] * d
        self.set_element(rop[0], x[3])
        self.set_element(rop[3], x[0])
        self.neg(t1, x[1])
        self.neg(t2, x[2])
        self.set_element(rop[1], t1)
        self.set_element(rop[2], t2)
        cdef int i
        for i in range(4):
            self.mul(rop[i], rop[i], d)
        return 0  # success

    cdef matrix_mul(self, modn_matrix rop, modn_matrix x, modn_matrix y):
        cdef modn_element t, t2
        self.mul(t, x[0], y[0])
        self.mul(t2, x[1], y[2])
        self.add(rop[0], t, t2)
        self.mul(t, x[0], y[1])
        self.mul(t2, x[1], y[3])
        self.add(rop[1], t, t2)
        self.mul(t, x[2], y[0])
        self.mul(t2, x[3], y[2])
        self.add(rop[2], t, t2)
        self.mul(t, x[2], y[1])
        self.mul(t2, x[3], y[3])
        self.add(rop[3], t, t2)

    cdef matrix_mul_ith_factor(self, modn_matrix rop, modn_matrix x, modn_matrix y, int i):
        cdef ResidueRing_abstract R = self.residue_rings[i]
        cdef residue_element t, t2
        R.mul(t, x[0][i], y[0][i])
        R.mul(t2, x[1][i], y[2][i])
        R.add(rop[0][i], t, t2)
        R.mul(t, x[0][i], y[1][i])
        R.mul(t2, x[1][i], y[3][i])
        R.add(rop[1][i], t, t2)
        R.mul(t, x[2][i], y[0][i])
        R.mul(t2, x[3][i], y[2][i])
        R.add(rop[2][i], t, t2)
        R.mul(t, x[2][i], y[1][i])
        R.mul(t2, x[3][i], y[3][i])
        R.add(rop[3][i], t, t2)

###########################################################################
# The projective line P^1(O_F/N)
###########################################################################
ctypedef modn_element p1_element[2]

cdef class ProjectiveLineModN:
    cdef ResidueRingModN S
    cdef int r
    cdef long _cardinality
    cdef long cards[MAX_PRIME_DIVISORS]  # cardinality of factors
    def __init__(self, N):
        self.S = ResidueRingModN(N)
        self.r = self.S.r  # number of prime divisors
        self._cardinality = 1
        cdef int i
        for i in range(self.r):
            R = self.S.residue_rings[i]
            self.cards[i] = (R.cardinality() + R.cardinality()/R.residue_field().cardinality())
            self._cardinality *= self.cards[i]

    cpdef long cardinality(self):
        return self._cardinality

    def __repr__(self):
        return "Projective line modulo the ideal %s of norm %s"%(
            self.S.N._repr_short(), self.S.N.norm())

##     cdef element_to_str(self, p1_element op):
##         return '(%s : %s)'%(self.S.element_to_str(op[0]), self.S.element_to_str(op[1]))

    cdef element_to_str(self, p1_element op):
        cdef ResidueRing_abstract R
        v = [ ]
        for i in range(self.r):
            R = self.S.residue_rings[i]
            v.append('(%s : %s)'%(R.element_to_str(op[0][i]), R.element_to_str(op[1][i])))
        s = ', '.join(v)
        if self.r > 1:
            s = '(' + s + ')'
        return s

    cdef int matrix_action(self, p1_element rop, modn_matrix op0, p1_element op1) except -1:
        # op0  = [a,b; c,d];    op1=[u;v]
        # product is  [a*u+b*v; c*u+d*v]
        cdef modn_element t0, t1
        self.S.mul(t0, op0[0], op1[0])  # a*u
        self.S.mul(t1, op0[1], op1[1])  # b*v
        self.S.add(rop[0], t0, t1)      # rop[0] = a*u + b*v
        self.S.mul(t0, op0[2], op1[0])  # c*u
        self.S.mul(t1, op0[3], op1[1])  # d*v
        self.S.add(rop[1], t0, t1)      # rop[1] = c*u + d*v

    cdef void set_element(self, p1_element rop, p1_element op):
        self.S.set_element(rop[0], op[0])
        self.S.set_element(rop[1], op[1])

    ######################################################################
    # Reducing elements to canonical form, so can compute index
    ######################################################################

    cdef int reduce_element(self, p1_element rop, p1_element op) except -1:
        # set rop to the result of reducing op to canonical form
        cdef int i
        for i in range(self.r):
            self.ith_reduce_element(rop, op, i)

    cdef int ith_reduce_element(self, p1_element rop, p1_element op, int i) except -1:
        # If the ith factor is (a,b), then, as explained in the next_element
        # docstring, the normalized form of this element is either:
        #      (u,1) or (1,v) with v in P.
        # The code below is careful to allow for the case when op and rop
        # are actually the same object, since this is a standard use case.
        cdef residue_element inv, r0, r1
        cdef ResidueRing_abstract R = self.S.residue_rings[i]
        if R.is_unit(op[1][i]):
            # in first case
            R.inv(inv, op[1][i])
            R.mul(r0, op[0][i], inv)
            R.set_element_to_1(r1)
        else:
            # can't invert b, so must be in second case
            R.set_element_to_1(r0)
            R.inv(inv, op[0][i])
            R.mul(r1, inv, op[1][i])
        R.set_element(rop[0][i], r0)
        R.set_element(rop[1][i], r1)
        return 0

    def test_reduce(self):
        """
        The test passes if this returns the empty list.
        Uses different random numbers each time, so seed the
        generator if you want control.
        """
        cdef p1_element x, y, z
        cdef long a
        self.first_element(x)
        v = []
        from random import randrange
        cdef ResidueRing_abstract R
        while True:
            # get y from x by multiplying through each entry by 3
            for i in range(self.r):
                R = self.S.residue_rings[i]
                a = randrange(1, R.p)
                y[0][i][0] = (a*x[0][i][0])%R.n0
                y[0][i][1] = (a*x[0][i][1])%R.n1
                y[1][i][0] = (a*x[1][i][0])%R.n0
                y[1][i][1] = (a*x[1][i][1])%R.n1
            self.reduce_element(z, y)
            v.append((self.element_to_str(x), self.element_to_str(y), self.element_to_str(z)))
            try:
                self.next_element(x, x)
            except ValueError:
                return [w for w in v if w[0] != w[2]]

    def test_reduce2(self, int m, int n):
        cdef int i
        cdef p1_element x, y, z
        self.first_element(x)
        for i in range(m):
            self.next_element(x, x)
        print self.element_to_str(x)
        cdef ResidueRing_abstract R
        for i in range(self.r):
            R = self.S.residue_rings[i]
            y[0][i][0] = (3*x[0][i][0])%R.n0
            y[0][i][1] = (3*x[0][i][1])%R.n1
            y[1][i][0] = (3*x[1][i][0])%R.n0
            y[1][i][1] = (3*x[1][i][1])%R.n1
        print self.element_to_str(y)
        from sage.all import cputime
        t = cputime()
        for i in range(n):
            self.reduce_element(z, y)
        return cputime(t)


    ######################################################################
    # Standard index of elements
    ######################################################################

    cdef long standard_index(self, p1_element z) except -1:
        # return the standard index of the assumed reduced element z of P^1
        # The global index of z is got from the local indexes at each factor.
        # If we let C_i denote the cardinality of the i-th factor, then the
        # global index I is:
        #         ind(z) = ind(z_0) + ind(z_1)*C_0 + ind(z_2)*C_0*C_1 + ...
        cdef int i
        cdef long C=1, ind = self.ith_standard_index(z, 0)
        for i in range(1, self.r):
            C *= self.cards[i-1]
            ind += C * self.ith_standard_index(z, i)
        return ind

    cdef long ith_standard_index(self, p1_element z, int i) except -1:
        cdef ResidueRing_abstract R = self.S.residue_rings[i]
        # Find standard index of normalized element
        #        (z[0][i], z[1][i])
        # of R x R.  The index is defined by the ordering on
        # normalized elements given by the next_element method (see
        # docs below).

        if R.element_is_1(z[1][i]):
            # then the index is the index of the first entry
            return R.index_of_element(z[0][i])

        # The index is the cardinality of R plus the index of z[1][i]
        # in the list of multiples of P in R.
        return R._cardinality + R.index_of_element_in_P(z[1][i])

    ######################################################################
    # Enumeration of elements
    ######################################################################

    def test_enum(self):
        cdef p1_element x
        self.first_element(x)
        v = []
        while True:
            c = (self.element_to_str(x), self.standard_index(x))
            print c
            v.append(c)
            try:
                self.next_element(x, x)
            except ValueError:
                return v

    def test_enum2(self):
        cdef p1_element x
        self.first_element(x)
        while True:
            try:
                self.next_element(x, x)
            except ValueError:
                return

    cdef int first_element(self, p1_element rop) except -1:
        # set rop to the first standard element, which is 1 in every factor
        cdef int i
        for i in range(self.r):
            # The 2-tuple (0,1) represents the 1 element in all the residue rings.
            rop[0][i][0] = 0
            rop[0][i][1] = 0
            rop[1][i][0] = 1
            rop[1][i][1] = 0
        return 0

    cdef int next_element(self, p1_element rop, p1_element z) except -1:
        # set rop equal to the next standard element after z.
        # The input z is assumed standard!
        # The enumeration is to start with the first ring, enumerate its P^1, and
        # if rolls over, go to the next ring, etc.
        # At the end, raise ValueError

        if rop != z:  # not the same C-level pointer
            self.set_element(rop, z)

        cdef int i=0
        while i < self.r and self.next_element_factor(rop, i):
            # it rolled over
            # Reset i-th one to starting P^1 elt, and increment the (i+1)th
            rop[0][i][0] = 0
            rop[0][i][1] = 0
            rop[1][i][0] = 1
            rop[1][i][1] = 0
            i += 1
        if i == self.r: # we're done
            raise ValueError

        return 0

    cdef bint next_element_factor(self, p1_element op, int i) except -2:
        # modify op in place by replacing the element in the i-th P^1 factor by
        # the next element.  If this involves rolling around, return True; otherwise, False.

        cdef ResidueRing_abstract k = self.S.residue_rings[i]
        # The P^1 local factor is (a,b) where a = op[0][i] and b = op[1][i].
        # Our "abstract" enumeration of the local P^1 with residue ring k is:
        #    [(a,1) for a in k] U [(1,b) for b in P*k]
        if k.element_is_1(op[1][i]):   # b == 1
            # iterate a
            if k.is_last_element(op[0][i]):
                # Then next elt is (1,b) where b is the first element of P*k.
                # So set b to first element in P, which is 0.
                k.set_element_to_0(op[1][i])
                # set a to 1
                k.set_element_to_1(op[0][i])
            else:
                k.next_element(op[0][i], op[0][i])
            return False # definitely didn't loop around whole P^1
        else:
            # case when b != 1
            if k.is_last_element_in_P(op[1][i]):
                # Next element is (1,0) and we return 1 to indicate total loop around
                k.set_element_to_0(op[1][i])
                return True # looped around
            else:
                k.next_element_in_P(op[1][i], op[1][i])
            return False

# Version using exception handling -- it's about 20% percent slower...
##     cdef bint next_element_factor(self, p1_element op, int i):
##         # modify op in place by replacing the element in the i-th P^1 factor by
##         # the next element.  If this involves rolling around, return True; otherwise,
##         # return False.
##         cdef ResidueRing_abstract k = self.S.residue_rings[i]
##         # The P^1 local factor is (a,b) where a = op[0][i] and b = op[1][i].
##         # Our "abstract" enumeration of the local P^1 with residue ring k is:
##         #    [(a,1) for a in k] U [(1,b) for b in P*k]
##         if k.element_is_1(op[1][i]):   # b == 1
##             # iterate a
##             try:
##                 k.next_element(op[0][i], op[0][i])
##             except ValueError:
##                 # Then next elt is (1,b) where b is the first element of P*k.
##                 # So set b to first element in P, which is 0.
##                 k.set_element_to_0(op[1][i])
##                 # set a to 1
##                 k.set_element_to_1(op[0][i])
##             return 0 # definitely didn't loop around whole P^1
##         else:
##             # case when b != 1
##             try:
##                 k.next_element_in_P(op[1][i], op[1][i])
##             except ValueError:
##                 # Next element is (1,0) and we return 1 to indicate total loop around
##                 k.set_element_to_0(op[1][i])
##                 return 1 # looped around
##             return 0


# readable version of the function below (not used):
def quaternion_in_terms_of_icosian_basis2(alpha):
    """
    Write the quaternion alpha in terms of the fixed basis for integer icosians.
    """
    b0 = alpha[0]; b1=alpha[1]; b2=alpha[2]; b3=alpha[3]
    a = F.gen()
    return [F_two*b0,
            (a-F_one)*b0 - b1 + a*b3,
            (a-F_one)*b0 + a*b1 - b2,
            (-a-F_one)*b0 + a*b2 - b3]

cdef NumberFieldElement_quadratic F_one = F(1), F_two = F(2), F_gen = F.gen(), F_gen_m1=F_gen-F_one, F_gen_m2 = -F_gen-F_one

cpdef object quaternion_in_terms_of_icosian_basis(object alpha):
    """
    Write the quaternion alpha in terms of the fixed basis for integer icosians.
    """
    # this is a critical path for speed.
    cdef NumberFieldElement_quadratic b0 = alpha[0], b1=alpha[1], b2=alpha[2], b3=alpha[3],
    return [F_two._mul_(b0), F_gen_m1._mul_(b0)._sub_(b1)._add_(F_gen._mul_(b3)),
            F_gen_m1._mul_(b0)._add_(F_gen._mul_(b1))._sub_(b2),
            F_gen_m2._mul_(b0)._add_(F_gen._mul_(b2))._sub_(b3)]


####################################################################
# Reduction from Quaternion Algebra mod N.
####################################################################
cdef class ModN_Reduction:
    cdef ResidueRingModN S
    cdef modn_matrix G[4]
    cdef bint is_odd
    def __init__(self, N, bint init=True):
        cdef int i

        if not is_ideal_in_F(N):
            raise TypeError, "N must be an ideal of F"

        cdef ResidueRingModN S = ResidueRingModN(N)
        self.S = S
        self.is_odd =  (N.norm() % 2 != 0)

        if init:
            # Set the identity matrix (2 part of this will get partly overwritten when N is even.)
            S.set(self.G[0][0], 1); S.set(self.G[0][1], 0); S.set(self.G[0][2], 0); S.set(self.G[0][3], 1)

            for i in range(S.r):
                self.compute_ith_local_splitting(i)

    def reduce_exponent_of_ith_factor(self, i):
        """
        Let P be the ith prime factor of the modulus N of self. This function
        returns the splitting got by reducing self modulo N/P.  See
        the documentation for the degeneracy_matrix of
        IcosiansModP1ModN for why we wrote this function.
        """
        cdef ResidueRing_abstract R = self.S.residue_rings[i]
        P = R.P
        N2 = self.S.N / P
        cdef ModN_Reduction f = ModN_Reduction(N2, init=False)
        # We have to set the matrices f.G[i] for i=0,1,2,3, using the
        # set_element_reducing_exponent_of_ith_factor method.
        cdef int j, k
        for j in range(4):
            for k in range(4):
                self.S.set_element_reducing_exponent_of_ith_factor(f.G[j][k], self.G[j][k], i)
        return f

    cdef compute_ith_local_splitting(self, int i):
        cdef ResidueRingElement z
        cdef long m, n
        cdef ResidueRing_abstract R = self.S.residue_rings[i]
        F = self.S.N.number_field()

        if R.p != 2:
            # I = [0,-1, 1, 0] works.
            R.set_element_to_0(self.G[1][0][i])
            R.set_element_to_1(self.G[1][1][i]); R.neg(self.G[1][1][i], self.G[1][1][i])
            R.set_element_to_1(self.G[1][2][i])
            R.set_element_to_0(self.G[1][3][i])
        else:
            from sqrt5_fast_python import find_mod2pow_splitting
            w = find_mod2pow_splitting(R.e)
            for n in range(4):
                v = w[n].list()
                for m in range(4):
                    z = v[m]
                    self.G[n][m][i][0] = z.x[0]
                    self.G[n][m][i][1] = z.x[1]
            return

        #######################################################################
        # Now find J:
        # TODO: *IMPORTANT* -- this must get redone in a way that is
        # much faster when the power of the prime is large.  Right now
        # it does a big enumeration, which is very slow.  There is a
        # Hensel lift style approach that is dramatically massively
        # faster.  See the code for find_mod2pow_splitting above,
        # which uses an iterative lifting algorithm.
        # Must implement this... once everything else is done. This is also slow since my sqrt
        # method is ridiculously stupid.  A good test is
        # sage: time ModN_Reduction(F.primes_above(7)[0]^5)
        # CPU times: user 0.65 s, sys: 0.00 s, total: 0.65 s
        # which should be instant.  And don't try exponent 6, which takes minutes (at least)!
        #######################################################################
        cdef residue_element a, b, c, d, t, t2, minus_one
        found_it = False
        R.set_element_to_1(b)
        R.coerce_from_nf(minus_one, F(-1))
        while not R.is_last_element(b):
            # Let c = -1 - b*b
            R.mul(t, b, b)
            R.mul(t, t, minus_one)
            R.add(c, minus_one, t)
            if R.is_square(c):
                # Next set a = -sqrt(c).
                R.sqrt(a, c)
                # Set the matrix self.G[2] to [a,b,(j2-a*a)/b,-a]
                R.set_element(self.G[2][0][i], a)
                R.set_element(self.G[2][1][i], b)
                # Set t to (-1-a*a)/b
                R.mul(t, a, a)
                R.mul(t, t, minus_one)
                R.add(t, t, minus_one)
                R.inv(t2, b)
                R.mul(self.G[2][2][i], t, t2)
                R.mul(self.G[2][3][i], a, minus_one)

                # Check that indeed we found an independent matrices, over the residue field
                good, K = self._matrices_are_independent_mod_ith_prime(i)
                if good:
                    self.S.matrix_mul_ith_factor(self.G[3], self.G[1], self.G[2], i)
                    return
            R.next_element(b, b)

        raise NotImplementedError

    def _matrices_are_independent_mod_ith_prime(self, int i):
        cdef ResidueRing_abstract R = self.S.residue_rings[i]
        k = R.residue_field()
        M = MatrixSpace(k, 2)
        J = M([R.element_to_residue_field(self.G[2][n][i]) for n in range(4)])
        I = M([R.element_to_residue_field(self.G[1][n][i]) for n in range(4)])
        K = I * J
        B = [M.identity_matrix(), I, J, I*J]
        V = k**4
        W = V.span([x.list() for x in B])
        return W.dimension() == 4, K.list()

    def __repr__(self):
        return 'Reduction modulo %s of norm %s defined by sending I to %s and J to %s'%(
            self.S.N._repr_short(), self.S.N.norm(),
            self.S.matrix_to_str(self.G[1]), self.S.matrix_to_str(self.G[2]))

    cdef int quatalg_to_modn_matrix(self, modn_matrix M, alpha) except -1:
        # Given an element alpha in the quaternion algebra, find its image M as a modn_matrix.
        cdef modn_element t
        cdef modn_element X[4]
        cdef int i, j
        cdef ResidueRingModN S = self.S
        cdef ResidueRing_abstract R
        cdef residue_element z[4]
        cdef residue_element t2

        for i in range(4):
            S.coerce_from_nf(X[i], alpha[i], 1)
        for i in range(4):
            S.mul(M[i], self.G[0][i], X[0])
            for j in range(1, 4):
                S.mul(t, self.G[j][i], X[j])
                S.add(M[i], M[i], t)

        if not self.is_odd:
            # The first prime is 2, so we have to fix that the above was
            # totally wrong at 2.
            # TODO: I'm writing this code slowly to work rather
            # than quickly, since I'm running out of time.
            # Will redo this to be fast later.  The algorithm is:
            # 1. Write alpha in terms of Icosian ring basis.
            # 2. Take the coefficients from *that* linear combination
            #    and apply the map, just as above, but of course only
            #    to the first factor of M.
            v = quaternion_in_terms_of_icosian_basis(alpha)
            # Now v is a list of 4 elements of O_F (unless alpha wasn't in R).
            R = self.S.residue_rings[0]
            assert R.p == 2, '%s\nBUG: order of factors wrong? '%(self.S.residue_rings,)
            for i in range(4):
                R.coerce_from_nf(z[i], v[i])
            for i in range(4):
                R.mul(M[i][0], self.G[0][i][0], z[0])
                for j in range(1,4):
                    R.mul(t2, self.G[j][i][0], z[j])
                    R.add(M[i][0], M[i][0], t2)


    def __call__(self, alpha):
        """
        A sort of joke for now.  Reduce alpha using this map, then return
        string representation...
        """
        cdef modn_matrix M
        self.quatalg_to_modn_matrix(M, alpha)
        return self.S.matrix_to_str(M)


####################################################################
# R^* \ P1(O_F/N)
####################################################################

ctypedef modn_matrix icosian_matrices[120]

cdef class IcosiansModP1ModN:
    """
    Create object that represents that set of orbits
          (Icosian Group) \ P^1(O_F / N).

    EXAMPLES::

        sage: from psage.modform.hilbert.sqrt5.sqrt5_fast import F, IcosiansModP1ModN
        sage: I = IcosiansModP1ModN(F.prime_above(31)); I
        The 2 orbits for the action of the Icosian group on Projective line modulo the ideal (5*a - 2) of norm 31
        sage: I = IcosiansModP1ModN(F.prime_above(389)); I
        The 7 orbits for the action of the Icosian group on Projective line modulo the ideal (18*a - 5) of norm 389
        sage: I = IcosiansModP1ModN(F.prime_above(2011)); I
        The 35 orbits for the action of the Icosian group on Projective line modulo the ideal (41*a - 11) of norm 2011

        sage: from psage.modform.hilbert.sqrt5.tables import ideals_of_bounded_norm
        sage: [len(IcosiansModP1ModN(b)) for b in ideals_of_bounded_norm(300)]
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 1, 1, 3, 3, 2, 2, 2, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 3, 4, 4, 3, 3, 3, 3, 3, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 6, 5, 5, 4, 4, 6, 4, 4, 6, 6, 4, 4, 4, 4, 5, 5, 6, 6, 6, 5, 5, 5, 5, 4, 4, 6, 6, 7, 7, 6, 5, 5, 6, 6, 6, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6]
    """
    cdef icosian_matrices G
    cdef ModN_Reduction f
    cdef ProjectiveLineModN P1
    cdef long *std_to_rep_table
    cdef long *orbit_reps
    cdef long _cardinality
    cdef p1_element* orbit_reps_p1elt
    def __init__(self, N, f=None, init=True):
        # compute choice of splitting, unless one is provided
        if f is not None:
            self.f = f
        else:
            self.f = ModN_Reduction(N)
        self.P1 = ProjectiveLineModN(N)
        self.orbit_reps = <long*>0
        self.std_to_rep_table = <long*> sage_malloc(sizeof(long) * self.P1.cardinality())

        # TODO: This is very wasteful of memory, by a factor of about 120 on average?
        self.orbit_reps_p1elt = <p1_element*>sage_malloc(sizeof(p1_element) * self.P1.cardinality())

        # initialize the group G of the 120 mod-N icosian matrices
        from sqrt5 import all_icosians
        X = all_icosians()
        cdef int i
        for i in range(len(X)):
            self.f.quatalg_to_modn_matrix(self.G[i], X[i])
            #print "%s --> %s"%(X[i], self.P1.S.matrix_to_str(self.G[i]))

        if init:
            self.compute_std_to_rep_table()

    def __reduce__(self):
        return unpickle_IcosiansModP1ModN_v1, (self.P1.S.N.gens_reduced()[0], )

    def __len__(self):
        return self._cardinality

    def __dealloc__(self):
        sage_free(self.std_to_rep_table)
        sage_free(self.orbit_reps_p1elt)
        if self.orbit_reps:
            sage_free(self.orbit_reps)

    def __repr__(self):
        return "The %s orbits for the action of the Icosian group on %s"%(self._cardinality, self.P1)

    def level(self):
        """
        Return the level N, which is the ideal we quotiented out by in
        constructing this set.

        EXAMPLES::

        """
        return self.P1.S.N

    cpdef compute_std_to_rep_table(self):
        # Algorithm is pretty obvious.  Just take first element of P1,
        # hit by all icosians, making a list of those that are
        # equivalent to it.  Then find next element of P1 not in the
        # first list, etc.
        cdef long orbit_cnt
        cdef p1_element x, Gx
        self.P1.first_element(x)
        cdef long ind=0, j, i=0, k
        for j in range(self.P1._cardinality):
            self.std_to_rep_table[j] = -1
        self.std_to_rep_table[0] = 0
        orbit_cnt = 1
        reps = []
        while ind < self.P1._cardinality:
            reps.append(ind)
            self.P1.set_element(self.orbit_reps_p1elt[i], x)
            for j in range(120):
                self.P1.matrix_action(Gx, self.G[j], x)
                self.P1.reduce_element(Gx, Gx)
                k = self.P1.standard_index(Gx)
                if self.std_to_rep_table[k] == -1:
                    self.std_to_rep_table[k] = i
                    orbit_cnt += 1
                else:
                    # This is a very good test that we got things right.  If any
                    # arithmetic or reduction is wrong, the orbits for the R^*
                    # "action" are likely to fail to be disjoint.
                    assert self.std_to_rep_table[k] == i, "Bug: orbits not disjoint"
            # This assertion below is also an extremely good test that we got the
            # local splittings, etc., right.  If any of that goes wrong, then
            # the "orbits" have all kinds of random orders.
            assert 120 % orbit_cnt == 0, "orbit size = %s must divide 120"%orbit_cnt
            orbit_cnt = 0
            while ind < self.P1._cardinality and self.std_to_rep_table[ind] != -1:
                ind += 1
                if ind < self.P1._cardinality:
                    self.P1.next_element(x, x)
            i += 1
        self._cardinality = len(reps)
        self.orbit_reps = <long*> sage_malloc(sizeof(long)*self._cardinality)
        for j in range(self._cardinality):
            self.orbit_reps[j] = reps[j]

    def compute_std_to_rep_table_debug(self):
        # Algorithm is pretty obvious.  Just take first element of P1,
        # hit by all icosians, making a list of those that are
        # equivalent to it.  Then find next element of P1 not in the
        # first list, etc.
        cdef long orbit_cnt
        cdef p1_element x, Gx
        self.P1.first_element(x)
        cdef long ind=0, j, i=0, k
        for j in range(self.P1._cardinality):
            self.std_to_rep_table[j] = -1
        self.std_to_rep_table[0] = 0
        orbit_cnt = 1
        reps = []
        while ind < self.P1._cardinality:
            reps.append(ind)
            print "Found representative number %s (which has standard index %s): %s"%(
                i, ind, self.P1.element_to_str(x))
            self.P1.set_element(self.orbit_reps_p1elt[i], x)
            for j in range(120):
                self.P1.matrix_action(Gx, self.G[j], x)
                self.P1.reduce_element(Gx, Gx)
                k = self.P1.standard_index(Gx)
                if self.std_to_rep_table[k] == -1:
                    self.std_to_rep_table[k] = i
                    orbit_cnt += 1
                else:
                    # This is a very good test that we got things right.  If any
                    # arithmetic or reduction is wrong, the orbits for the R^*
                    # "action" are likely to fail to be disjoint.
                    assert self.std_to_rep_table[k] == i, "Bug: orbits not disjoint"
            # This assertion below is also an extremely good test that we got the
            # local splittings, etc., right.  If any of that goes wrong, then
            # the "orbits" have all kinds of random orders.
            assert 120 % orbit_cnt == 0, "orbit size = %s must divide 120"%orbit_cnt
            print "It has an orbit of size %s"%orbit_cnt
            orbit_cnt = 0
            while ind < self.P1._cardinality and self.std_to_rep_table[ind] != -1:
                ind += 1
                if ind < self.P1._cardinality:
                    self.P1.next_element(x, x)
            i += 1
        self._cardinality = len(reps)
        self.orbit_reps = <long*> sage_malloc(sizeof(long)*self._cardinality)
        for j in range(self._cardinality):
            self.orbit_reps[j] = reps[j]

    cpdef long cardinality(self):
        return self._cardinality

    def hecke_matrix(self, P, sparse=True):
        """
        Return matrix of Hecke action, acting from the right, so the
        i-th row is the image of the i-th standard basis vector.

        INPUT:
            - P -- prime ideal
            - ``sparse`` -- bool (default: True) if True, then a sparse
              matrix is returned, otherwise a dense one

        OUTPUT:
            - a dense or sparse matrix over the rational numbers

        NOTE: This function does not cache its results.

        EXAMPLES::

            sage: from psage.modform.hilbert.sqrt5.sqrt5_fast import F, IcosiansModP1ModN
            sage: I = IcosiansModP1ModN(F.primes_above(31)[0]); I
            The 2 orbits for the action of the Icosian group on Projective line modulo the ideal (5*a - 2) of norm 31
            sage: I.hecke_matrix(F.prime_above(2))
            [0 5]
            [3 2]
            sage: I.hecke_matrix(F.prime_above(3))
            [5 5]
            [3 7]
            sage: I.hecke_matrix(F.prime_above(5))
            [1 5]
            [3 3]
            sage: I.hecke_matrix(F.prime_above(3)).fcp()
            (x - 10) * (x - 2)
            sage: I.hecke_matrix(F.prime_above(2)).fcp()
            (x - 5) * (x + 3)
            sage: I.hecke_matrix(F.prime_above(5)).fcp()
            (x - 6) * (x + 2)

        A bigger example::

            sage: I = IcosiansModP1ModN(F.prime_above(389)); I
            The 7 orbits for the action of the Icosian group on Projective line modulo the ideal (18*a - 5) of norm 389
            sage: t2 = I.hecke_matrix(F.prime_above(2)); t2
            [0 3 0 1 1 0 0]
            [3 0 0 0 1 0 1]
            [0 0 2 1 0 1 1]
            [1 0 1 0 1 0 2]
            [1 1 0 1 0 1 1]
            [0 0 2 0 2 1 0]
            [0 1 1 2 1 0 0]
            sage: t2.is_sparse()
            True
            sage: I.hecke_matrix(F.prime_above(2), sparse=False).is_sparse()
            False
            sage: t3 = I.hecke_matrix(F.prime_above(3))
            sage: t5 = I.hecke_matrix(F.prime_above(5))
            sage: t11a = I.hecke_matrix(F.primes_above(11)[0])
            sage: t11b = I.hecke_matrix(F.primes_above(11)[1])
            sage: t2.fcp()
            (x - 5) * (x^2 + 5*x + 5) * (x^4 - 3*x^3 - 3*x^2 + 10*x - 4)
            sage: t3.fcp()
            (x - 10) * (x^2 + 3*x - 9) * (x^4 - 5*x^3 + 3*x^2 + 6*x - 4)
            sage: t5.fcp()
            (x - 6) * (x^2 + 4*x - 1) * (x^2 - x - 4)^2
            sage: t11a.fcp()
            (x - 12) * (x + 3)^2 * (x^4 - 17*x^2 + 68)
            sage: t11b.fcp()
            (x - 12) * (x^2 + 5*x + 5) * (x^4 - x^3 - 23*x^2 + 18*x + 52)
            sage: t2*t3 == t3*t2, t2*t5 == t5*t2, t11a*t11b == t11b*t11a
            (True, True, True)

        Examples involving a prime dividing the level::

            sage: from psage.modform.hilbert.sqrt5.sqrt5_fast import F, IcosiansModP1ModN
            sage: v = F.primes_above(31); I = IcosiansModP1ModN(v[0])
            sage: t31 = I.hecke_matrix(v[1]); t31
            [17 15]
            [ 9 23]
            sage: t31.fcp()
            (x - 32) * (x - 8)
            sage: t5 = I.hecke_matrix(F.prime_above(5))
            sage: t5*t31 == t31*t5
            True
            sage: I.hecke_matrix(v[0])
            Traceback (most recent call last):
            ...
            NotImplementedError: computation of T_P for P dividing the level is not implemented

        Another test involving a prime that divides the norm of the level::

            sage: v = F.primes_above(809); I = IcosiansModP1ModN(v[0])
            sage: t = I.hecke_matrix(v[1])
            sage: I
            The 14 orbits for the action of the Icosian group on Projective line modulo the ideal (-26*a + 7) of norm 809
            sage: t.fcp()
            (x - 810) * (x + 46) * (x^4 - 48*x^3 - 1645*x^2 + 65758*x + 617743) * (x^8 + 52*x^7 - 2667*x^6 - 118268*x^5 + 2625509*x^4 + 55326056*x^3 - 1208759312*x^2 + 4911732368*x - 4279699152)
            sage: s = I.hecke_matrix(F.prime_above(11))
            sage: t*s == s*t
            True
        """
        if ZZ(self.level().norm()).gcd(ZZ(P.norm())) != 1:
            return self._hecke_matrix_badnorm(P, sparse=sparse)
        cdef modn_matrix M
        cdef p1_element Mx
        from sqrt5 import hecke_elements

        T = zero_matrix(QQ, self._cardinality, sparse=sparse)
        cdef long i, j

        for alpha in hecke_elements(P):
            self.f.quatalg_to_modn_matrix(M, alpha)
            for i in range(self._cardinality):
                self.P1.matrix_action(Mx, M, self.orbit_reps_p1elt[i])
                self.P1.reduce_element(Mx, Mx)
                j = self.std_to_rep_table[self.P1.standard_index(Mx)]
                T[i,j] += 1
        return T


    def _hecke_matrix_badnorm(self, P, sparse=True):
        """
        Compute the Hecke matrix for a prime P whose norm divides the
        norm of the level.

        This function doesn't work if P itself divides the level.

        EXAMPLES::

            sage: from psage.modform.hilbert.sqrt5.sqrt5_fast import F, IcosiansModP1ModN; v=F.primes_above(71); I=IcosiansModP1ModN(v[0])
            sage: I._hecke_matrix_badnorm(v[1])
            [61 11]
            [55 17]
            sage: I._hecke_matrix_badnorm(v[1]).fcp()
            (x - 72) * (x - 6)
        """
        cdef modn_matrix M, M0
        cdef p1_element Mx
        from sqrt5 import hecke_elements

        T = zero_matrix(QQ, self._cardinality, sparse=sparse)
        cdef long i, j

        if P.divides(self.level()):
            raise NotImplementedError, "computation of T_P for P dividing the level is not implemented"

        for beta in hecke_elements(P):
            alpha = beta**(-1)
            self.f.quatalg_to_modn_matrix(M0, alpha)
            if self.P1.S.matrix_is_invertible(M0):
                self.P1.S.matrix_inv(M, M0)
                for i in range(self._cardinality):
                    self.P1.matrix_action(Mx, M, self.orbit_reps_p1elt[i])
                    self.P1.reduce_element(Mx, Mx)
                    j = self.std_to_rep_table[self.P1.standard_index(Mx)]
                    T[i,j] += 1
        return T

    def hecke_operator_on_basis_element(self, P, long i):
        """
        Return image of i-th basis vector of self under the Hecke
        operator T_P, for P coprime to the level, computed efficiently
        in the sense of not computing images of all basis vectors.

        INPUT:
            - P -- nonzero prime ideal of ring of integers of Q(sqrt(5))
              that does not divide the level
            - i -- nonnegative integer less than the dimension

        OUTPUT:
            - a dense vector over the integers

        EXAMPLES::

            sage: from psage.modform.hilbert.sqrt5.sqrt5_fast import F, IcosiansModP1ModN
            sage: v = F.primes_above(31); I = IcosiansModP1ModN(v[0])
            sage: I.hecke_matrix(F.prime_above(2))
            [0 5]
            [3 2]
            sage: I.hecke_matrix(v[1])
            [17 15]
            [ 9 23]
            sage: I.hecke_operator_on_basis_element(F.prime_above(2), 0)
            (0, 5)
            sage: I.hecke_operator_on_basis_element(F.prime_above(2), 1)
            (3, 2)
            sage: I.hecke_operator_on_basis_element(v[1], 0)
            (17, 15)
            sage: I.hecke_operator_on_basis_element(v[1], 1)
            (9, 23)
        """
        cdef modn_matrix M, M0
        cdef p1_element Mx

        from sqrt5 import hecke_elements
        cdef list Q = hecke_elements(P)
        v = (ZZ**self._cardinality)(0)  # important -- do not use the vector function to make this: see trac 11657.

        if ZZ(self.level().norm()).gcd(ZZ(P.norm())) == 1:
            for alpha in Q:
                self.f.quatalg_to_modn_matrix(M, alpha)
                self.P1.matrix_action(Mx, M, self.orbit_reps_p1elt[i])
                self.P1.reduce_element(Mx, Mx)
                v[self.std_to_rep_table[self.P1.standard_index(Mx)]] += 1
        else:
            if P.divides(self.level()):
                raise NotImplementedError
            for beta in Q:
                alpha = beta**(-1)
                self.f.quatalg_to_modn_matrix(M0, alpha)
                if self.P1.S.matrix_is_invertible(M0):
                    self.P1.S.matrix_inv(M, M0)
                    self.P1.matrix_action(Mx, M, self.orbit_reps_p1elt[i])
                    self.P1.reduce_element(Mx, Mx)
                    v[self.std_to_rep_table[self.P1.standard_index(Mx)]] += 1
        return v


    def degeneracy_matrix(self, P, sparse=True):
        """
        Return degeneracy matrix from level N of self to level N/P.
        Here P must be a prime ideal divisor of the ideal N.
        """
        #################################################################################
        # Important comment / thought regarding this function!
        # For this to work, the local splitting maps must be chosen
        # in a compatible way.  Let R be the icosian ring.
        # We compute somehow a map
        #            R ---> M_2(O/P^n)
        # and somehow else, we compute a map
        #            R ---> M_2(O/P^(n+1)).
        # If we are trying to compute the degeneracy map from level P^(n+1)
        # to level P^n, and we choose incompatible maps, then our computation
        # will simply be WRONG.
        #   SO WATCH OUT!    (This took me like 3 hours of frustration to realize...)
        ##################################################################################


        N = self.P1.S.N
        # Figure out the index of P as a prime divisor of N
        cdef ResidueRing_abstract R
        cdef int k, ind = -1, m = len(self.P1.S.residue_rings)

        for k, R in enumerate(self.P1.S.residue_rings):
            if R.P == P:
                # Got it
                ind = k
                break
        if ind == -1:
            raise ValueError, "P must be a prime divisor of the level"

        # Compute basis of the P^1 for lower level.
        N2 = N/P

        # CRITICAL: see above comment.  We have to compute the mod-N2
        # local splitting map in terms of the fixed choice of map for
        # self.  If we don't, we will get very subtle bugs.
        g = self.f.reduce_exponent_of_ith_factor(ind)
        cdef IcosiansModP1ModN I2 = IcosiansModP1ModN(N2, g)

        D = zero_matrix(QQ, self._cardinality, I2._cardinality, sparse=sparse)
        cdef Py_ssize_t i, j
        cdef p1_element x, y

        for i in range(self._cardinality):
            # Make the p1_element y from self.orbit_reps_p1elt[i]
            # by reducing at the the "ind" position
            #print "reducing", self.P1.element_to_str(self.orbit_reps_p1elt[i])
            self.P1.S.set_element_reducing_exponent_of_ith_factor(
                                     y[0], self.orbit_reps_p1elt[i][0], ind)
            self.P1.S.set_element_reducing_exponent_of_ith_factor(
                                     y[1], self.orbit_reps_p1elt[i][1], ind)
            I2.P1.reduce_element(x, y)
            #print "... to ", I2.P1.element_to_str(x)
            j = I2.std_to_rep_table[I2.P1.standard_index(x)]
            D[i,j] += 1

        return D

def unpickle_IcosiansModP1ModN_v1(x):
    import sqrt5
    return IcosiansModP1ModN(sqrt5.F.ideal(x))


####################################################################
# fragments/test code below
####################################################################

## def test(v, int n):
##     cdef list w = v
##     cdef int i
##     cdef ResidueRing_abstract R
##     for i in range(n):
##         R = w[i%3]


## def test(ResidueRing_abstract R):
##     cdef modn_element x
##     x[0][0] = 5
##     x[0][1] = 7
##     x[1][0] = 2
##     x[1][1] = 8
##     R.add(x[2], x[0], x[1])
##     print x[2][0], x[2][1]


## def test_kr(long a, long b, long n):
##     cdef long i, j=0
##     for i in range(n):
##         j += kronecker_symbol_si(a,b)
##     return j
## def test1():
##     cdef int x[6]
##     x[2] = 1
##     x[3] = 2
##     x[4] = 1
##     x[5] = 2
##     from sqrt5 import F
##     Pinert = F.primes_above(7)[0]
##     cdef ResidueRing_nonsplit Rinert = ResidueRing(Pinert, 2)
##     print (x+2)[0], (x+2)[1]
##     Rinert.add(x, x+2, x+4)
##     return x[0], x[1], x[2], x[3]

## def bench1(int n):
##     from sqrt5 import F
##     Pinert = F.primes_above(3)[0]
##     Rinert = ResidueRing(Pinert, 2)
##     s_inert = Rinert(F.gen(0)+3)
##     t = s_inert
##     cdef int i
##     for i in range(n):
##         t = t * s_inert
##     return t

## def bench2(int n):
##     from sqrt5 import F
##     Pinert = F.primes_above(3)[0]
##     cdef ResidueRing_nonsplit Rinert = ResidueRing(Pinert, 2)
##     cdef ResidueRingElement_nonsplit2 s_inert = Rinert(F.gen(0)+3)
##     cdef residue_element t, s
##     s[0] = s_inert.x[0]
##     s[1] = s_inert.x[1]
##     t[0] = s[0]
##     t[1] = s[1]
##     cdef int i
##     for i in range(n):
##         Rinert.mul(t, t, s)
##     return t[0], t[1]

