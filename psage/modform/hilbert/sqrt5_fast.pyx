"""
Hilbert modular forms over F = Q(sqrt(5)).

All the code in this file is meant to be highly optimized. 
"""

include 'stdsage.pxi'
include 'cdefs.pxi'


from sage.rings.ring cimport CommutativeRing
from sage.rings.all import Integers

cdef long SQRT_MAX_LONG = 2**(4*sizeof(long))

# Residue element
ctypedef long residue_element[2]

#from sqrt5_fast_python import ResidueRing

cpdef long ideal_characteristic(I):
    return I.number_field()._pari_().idealtwoelt(I._pari_())[0]

###########################################################################
# Rings
###########################################################################

def ResidueRing(P, e):
    cdef long p = ideal_characteristic(P)
    if p**e >= SQRT_MAX_LONG:
        raise ValueError, "Residue field must have size less than %s"%SQRT_MAX_LONG
    if p == 5:   # ramified
        if e % 2 == 0:
            return ResidueRing_nonsplit(P, p, e)
        else:
            return ResidueRing_ramified_odd(P, p, e)
    elif p%5 in [2,3]:  # nonsplit
        return ResidueRing_nonsplit(P, p, e)
    else: # split
        return ResidueRing_split(P, p, e)

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
    cdef object P, e, F
    cdef public object element_class, _residue_field
    cdef long n0, n1, p, _cardinality
    cdef long im_gen0
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

    def residue_field(self):
        if self._residue_field is None:
            self._residue_field = self.P.residue_field()
        return self._residue_field

    def __call__(self, x):
        if hasattr(x, 'parent'):
            P = x.parent()
            if P is self:
                return x
            elif P is self.F:
                # TODO: This denominator is very slow.
                d = x.denominator()
                if d != 1:
                    return self(d*x) * self(self.F(d.inverse_mod(self.p**self.e)))
                return self.element_class(self, x)
        raise TypeError

    def __repr__(self):
        return "Residue class ring of %s^%s of characteristic %s"%(
            self.P._repr_short(), self.e, self.p)

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

    cdef int coerce_from_nf(self, residue_element rop, op) except -1:
        # TODO: this is probably super slow; definitely inefficient
        cdef ResidueRingElement z = self(op)
        rop[0] = z.x[0]
        rop[1] = z.x[1]
        return 0 # success

    cdef bint element_is_1(self, residue_element op):
        return op[0] == 1 and op[1] == 0

    cdef bint element_is_0(self, residue_element op):
        return op[0] == 0 and op[1] == 0
    
    cdef void elt_set_to_1(self, residue_element op):
        op[0] = 1
        op[1] = 0

    cdef void elt_set_to_0(self, residue_element op):
        op[0] = 0
        op[1] = 0

    cdef void elt_set(self, residue_element rop, residue_element op):
        rop[0] = op[0]
        rop[1] = op[1]

    cdef int pow(self, residue_element rop, residue_element op, long e) except -1:
        cdef residue_element op2
        if e < 0:
            self.inv(op2, op)
            return self.pow(rop, op2, -e)
        self.elt_set_to_1(rop)
        cdef residue_element z
        self.elt_set(z, op)
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
        rop[0] = self.n0 - op[0]
        rop[1] = 0        

    cdef bint is_square(self, residue_element op):
        # p is odd, so we just check if op[0] is a square modulo self.p
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
    def __init__(self, P, p, e):
        ResidueRing_abstract.__init__(self, P, p, e)
        self.element_class = ResidueRingElement_nonsplit
        self.n0 = p**e
        self.n1 = p**e
        self._cardinality = self.n0 * self.n1
        
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
        rop[0] = self.n0 - op[0]
        rop[1] = self.n1 - op[1]

    cdef bint is_square(self, residue_element op):
        if self.p == 2 and self.e >= 2:
            raise NotImplementedError, "TODO: implement square testing for p=2 in general"
        cdef residue_element z
        self.pow(z, op, (self.p*self.p-1)/2)
        return z[0] % self.p == 1
    
    cdef int sqrt(self, residue_element rop, residue_element op) except -1:
        k = self.residue_field()
        if k.degree() == 1:
            if self.e == 1:
                # happens in ramified case with odd exponent
                rop[0] = sqrtmod_long(op[0], self.p, 1)
                rop[1] = 0
                return 0
            raise NotImplementedError, 'sqrt not implemented in ramified case...'

        # TODO: the stupid overhead in this step alone is vastly more than
        # the time to actually compute the sqrt in Givaro or Pari (say)...
        a = k([op[0],op[1]]).sqrt()._vector_()
        
        # The rest of this function is fast (hundreds of times
        # faster than above line)...

        rop[0] = a[0]; rop[1] = a[1]
        if self.e == 1:
            return 0

        if rop[0] == 0 and rop[1] == 0:
            # TODO: Finish sqrt when number is 0 mod P in inert case.
            raise NotImplementedError
        
        # Hensel lifting to get square root modulo P^e.
        # See the code sqrtmod_long below, which is basically
        # the same, but more readable.
        cdef int m
        cdef long pm, ppm, p
        p = self.p; pm = p
        # By "x" below, we mean the input op.
        # By "y" we mean the square root just computed above, currently rop.
        cdef residue_element d, y_squared, y2, t, inv
        for m in range(1, self.e):
            ppm = p*pm
            # We compute ((x-y^2)/p^m) / (2*y)
            self.mul(y_squared, rop, rop)  # y2 = y^2
            self.sub(t, op, y_squared)
            assert t[0]%pm == 0  # TODO: remove this after some tests
            assert t[1]%pm == 0
            t[0] /= pm; t[1] /= pm
            # Now t = (x-y^2)/p^m.
            y2[0] = (2*rop[0])%ppm;   y2[1] = (2*rop[1])%ppm
            self.inv(inv, y2)
            self.mul(t, t, inv)
            # Now t = ((x-y^2)/p^m) / (2*y)
            t[0] *= pm;  t[1] *= pm
            # Now do y = y + pm * t
            self.add(rop, rop, t)
            
        return 0 # success

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
        rop[1] = op[1]
        if rop[0] >= self.n0:
            rop[0] = 0
            rop[1] += self.p
            if rop[1] >= self.n1:
                raise ValueError, "no next element in P"
        return 0
        
    cdef bint is_last_element_in_P(self, residue_element op):
        return op[0] == self.n0 - self.p and op[1] == self.n1 - self.p

    cdef long index_of_element_in_P(self, residue_element op) except -1:
        return op[0]/self.p + op[1]/self.p * (self.n0/self.p)
        
cdef class ResidueRing_ramified_odd(ResidueRing_abstract):
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
        rop[0] = self.n0 - op[0]
        rop[1] = self.n1 - op[1]

    cdef int inv(self, residue_element rop, residue_element op) except -1:
        # TODO: implement inverse in ramified case
        raise NotImplementedError

    cdef bint is_unit(self, residue_element op):
        raise NotImplementedError        

    cdef bint is_square(self, residue_element op):
        raise NotImplementedError

    cdef int sqrt(self, residue_element rop, residue_element op) except -1:
        raise NotImplementedError

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
        # TODO!
        raise NotImplementedError, "TODO: implement iteration over P in odd ramified case"

    cdef bint is_last_element_in_P(self, residue_element op):
        raise NotImplementedError, "TODO: implement iteration over P in odd ramified case"

    cdef long index_of_element_in_P(self, residue_element op) except -1:
        raise NotImplementedError, "TODO: implement iteration over P in odd ramified case"        
    

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
    cdef residue_element x
    cdef ResidueRing_abstract _parent
    cpdef parent(self):
        return self._parent
    cdef new(self):
        raise NotImplementedError
    
    def __add__(ResidueRingElement left, ResidueRingElement right):
        cdef ResidueRingElement z = left.new()
        left._parent.add(z.x, left.x, right.x)
        return z
    def __sub__(ResidueRingElement left, ResidueRingElement right):
        cdef ResidueRingElement z = left.new()
        left._parent.sub(z.x, left.x, right.x)
        return z
    def __mul__(ResidueRingElement left, ResidueRingElement right):
        cdef ResidueRingElement z = left.new()
        left._parent.mul(z.x, left.x, right.x)
        return z
    def __div__(ResidueRingElement left, ResidueRingElement right):
        cdef ResidueRingElement z = left.new()
        left._parent.inv(z.x, right.x)
        left._parent.mul(z.x, left.x, z.x)
        return z
    def __neg__(ResidueRingElement self):
        cdef ResidueRingElement z = self.new()        
        self._parent.neg(z.x, self.x)
        return z
    def __pow__(ResidueRingElement self, e, m):
        cdef ResidueRingElement z = self.new()
        self._parent.pow(z.x, self.x, e)
        return z
    
    def __invert__(ResidueRingElement self):
        cdef ResidueRingElement z = self.new()
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
        
    cpdef bint is_square(self):
        return self._parent.is_square(self.x)
    cpdef sqrt(self):
        cdef ResidueRingElement z = self.new()
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
            self.x[0] = v[0]
            return 
        self.x[0] = v[0] + parent.im_gen0*v[1]
        self.x[0] %= self._parent.n0
        self.x[1] = 0

    cdef new(self):
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

        sage: from psage.modform.hilbert.sqrt5_fast import ResidueRing

    Inert example::
    
        sage: F = NumberField(x**2 - x -1, 'a')
        sage: P = F.primes_above(3)[0]
        sage: R = ResidueRing(Pinert, 2)
        sage: s = Rinert(F.0); s
        0 + 1*gamma_bar
        sage: s + s + s
        0 + 3*gamma_bar
        sage: s*s*s*s
        2 + 3*gamma_bar
        sage: Rinert(F.0^4)
        2 + 3*gamma_bar

    Ramified example::

        sage: P = F.primes_above(5)[0]
        sage: R = ResidueRing(P, 2)
        sage: s = R(F.0); s
        0 + 1*gamma_bar
        sage: s*s*s*s*s
        3 + 0*gamma_bar
        sage: R(F.0^5)
        3 + 0*gamma_bar
        sage: a = (s+s - R(F(1))); a*a
        0 + 0*gamma_bar
        sage: R = ResidueRing(P, 3)
        sage: t = R(2*F.0-1)   # reduction of sqrt(5)
        sage: t*t*t
        0 + 0*gamma_bar
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
        self.x[0] %= self._parent.n0
        if self.x[0] < 0:
            self.x[0] += self._parent.n0
        self.x[1] %= self._parent.n0
        if self.x[1] < 0:
            self.x[1] += self._parent.n0

    cdef new(self):
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
        # We can assume that the defining poly of F is x^2-x-1 (this is asserted
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
            self.x[0] = a
            self.x[1] = 0
        else:
            a = v[0]; b = v[1]
            self.x[0] = a + b*parent.two_inv
            self.x[1] = b*parent.two_inv
        self.x[0] %= self._parent.n0
        if self.x[0] < 0:
            self.x[0] += self._parent.n0
        self.x[1] %= self._parent.n1
        if self.x[1] < 0:
            self.x[1] += self._parent.n1

    cdef new(self):
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

cdef long gcd_long(long a, long b):
    return Arith_llong.c_gcd_longlong(a, b)
    
cdef long divmod_long(long x, long y, long m):
    #return quotient x/y mod m, assuming x, y < SQRT_MAX_LONG
    return (x * invmod_long(y, m)) % m
    
cpdef long sqrtmod_long(long x, long p, int n) except -1:
    # Compute square root of x modulo p^n using Hensel lifting, where
    # p id odd.

    # Step 1. Find a square root modulo p. (See section 4.5 of my Elementary Number Theory book.)
    cdef long y
    if p % 4 == 3:
        # Easy case when p is 3 mod 4.
        y = powmod_long(x, (p+1)/4, p)
    else:
        # Harder case -- no known poly time deterministic algorithm.
        # TODO: There is a fast algorithm we could implement here, but
        # I'll just use pari for now and come back to this later. This
        # is on page 88 of my book, esp. since this is not critical
        # to the speed of the whole algorithm. 
        from sage.all import pari
        y = pari(x).Mod(p).sqrt().lift()

    if n == 1: # done
        return y

    if y == 0:
        # TODO: finish this
        print "TODO: Finish sqrtmod_long in case when y==0."
        raise NotImplementedError
    
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

cdef enum:
    MAX_PRIME_DIVISORS = 10

ctypedef residue_element modn_element[MAX_PRIME_DIVISORS]

# 2 x 2 matrices over O_F/N
ctypedef modn_element modn_matrix[4]


cdef class ResidueRingModN:
    cdef list residue_rings
    cdef object N
    cdef int r  # number of prime divisors
    
    def __init__(self, N):
        self.N = N
        self.residue_rings = [ResidueRing(P, e) for P, e in N.factor()]
        self.r = len(self.residue_rings)
            
    def __repr__(self):
        return "Residue class ring modulo the ideal %s of norm %s"%(
            self.N._repr_short(), self.N.norm())

    cdef int set(self, modn_element rop, x) except -1:
        self.coerce_from_nf(rop, self.N.number_field()(x))

    cdef int coerce_from_nf(self, modn_element rop,  op) except -1:
        # Given an element op in the field F, try to reduce it modulo
        # N and create the corresponding modn element.
        cdef int i
        cdef ResidueRing_abstract R
        for i in range(self.r):
            R = self.residue_rings[i]
            R.coerce_from_nf(rop[i], op)
        return 0 # success

    cdef void set_element(self, modn_element rop, modn_element op):
        cdef int i
        for i in range(self.r):
            rop[i][0] = op[i][0]
            rop[i][1] = op[i][1]

    cdef element_to_str(self, modn_element op):
        s = ','.join([(<ResidueRing_abstract>self.residue_rings[i]).element_to_str(op[i]) for
                          i in range(self.r)])
        if self.r > 1:
            s = '(' + s + ')'
        return s

    cdef void mul(self, modn_element rop, modn_element op0, modn_element op1):
        cdef int i
        cdef ResidueRing_abstract R
        for i in range(self.r):
            R = self.residue_rings[i]
            R.mul(rop[i], op0[i], op1[i])
        
    cdef void inv(self, modn_element rop, modn_element op):
        cdef int i
        cdef ResidueRing_abstract R
        for i in range(self.r):
            R = self.residue_rings[i]
            R.inv(rop[i], op[i])
        
    cdef void add(self, modn_element rop, modn_element op0, modn_element op1):
        cdef int i
        cdef ResidueRing_abstract R
        for i in range(self.r):
            R = self.residue_rings[i]
            R.add(rop[i], op0[i], op1[i])

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
                R.elt_set(rop[i], op[i])
            else:
                if R.is_last_element(op[i]):
                    R.elt_set_to_0(rop[i])
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
    # modn_matrices over R
    #######################################    
    cdef matrix_to_str(self, modn_matrix A):
        return '[%s,%s; %s,%s]'%tuple([self.element_to_str(A[i]) for i in range(4)])

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
            v.append('(%s, %s)'%(R.element_to_str(op[0][i]), R.element_to_str(op[1][i])))
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

    cdef void reduce_element(self, p1_element rop, p1_element op):
        # set rop to the result of reducing op to canonical form
        cdef int i
        for i in range(self.r):
            self.ith_reduce_element(rop, op, i)

    cdef void ith_reduce_element(self, p1_element rop, p1_element op, int i):
        # If the ith factor is (a,b), then, as explained in the next_element
        # docstring, the normalized form of this element is either:
        #      (u,1) or (1,v) with v in P.
        cdef residue_element inv
        cdef ResidueRing_abstract R = self.S.residue_rings[i]
        if R.is_unit(op[1][i]):
            # in first case
            R.inv(inv, op[1][i])
            R.mul(rop[0][i], op[0][i], inv)
            R.elt_set_to_1(rop[1][i])
        else:
            # can't invert b, so must be in second case
            R.elt_set_to_1(rop[0][i])
            R.inv(inv, op[0][i])
            R.mul(rop[1][i], inv, op[1][i])

    def test_reduce(self):
        cdef p1_element x, y, z
        self.first_element(x)
        v = []
        cdef ResidueRing_abstract R 
        while True:
            # get y from x by multiplying through each entry by 3
            for i in range(self.r):
                R = self.S.residue_rings[i]
                y[0][i][0] = (3*x[0][i][0])%R.n0
                y[0][i][1] = (3*x[0][i][1])%R.n1
                y[1][i][0] = (3*x[1][i][0])%R.n0
                y[1][i][1] = (3*x[1][i][1])%R.n1
            self.reduce_element(z, y)
            v.append((self.element_to_str(x), self.element_to_str(y), self.element_to_str(z)))
            try:
                self.next_element(x, x)
            except ValueError:
                return v
        
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
            v.append((self.element_to_str(x), self.standard_index(x)))
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

    cdef bint next_element_factor(self, p1_element op, int i):
        # modify op in place by replacing the element in the i-th P^1 factor by
        # the next element.  If this involves rolling around, return True; otherwise,
        # return False.

        cdef ResidueRing_abstract k = self.S.residue_rings[i]
        # The P^1 local factor is (a,b) where a = op[0][i] and b = op[1][i].
        # Our "abstract" enumeration of the local P^1 with residue ring k is:
        #    [(a,1) for a in k] U [(1,b) for b in P*k]
        if k.element_is_1(op[1][i]):   # b == 1
            # iterate a
            if k.is_last_element(op[0][i]):
                # Then next elt is (1,b) where b is the first element of P*k.
                # So set b to first element in P, which is 0.
                k.elt_set_to_0(op[1][i])
                # set a to 1
                k.elt_set_to_1(op[0][i])
            else:
                k.next_element(op[0][i], op[0][i])
            return 0 # definitely didn't loop around whole P^1
        else:
            # case when b != 1
            if k.is_last_element_in_P(op[1][i]):
                # Next element is (1,0) and we return 1 to indicate total loop around
                k.elt_set_to_0(op[1][i])
                return 1 # looped around
            else:
                k.next_element_in_P(op[1][i], op[1][i])
            return 0

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
##                 k.elt_set_to_0(op[1][i])
##                 # set a to 1
##                 k.elt_set_to_1(op[0][i])
##             return 0 # definitely didn't loop around whole P^1
##         else:
##             # case when b != 1
##             try:
##                 k.next_element_in_P(op[1][i], op[1][i])
##             except ValueError:
##                 # Next element is (1,0) and we return 1 to indicate total loop around
##                 k.elt_set_to_0(op[1][i])
##                 return 1 # looped around
##             return 0


####################################################################
# Reduction from Quaternion Algebra mod N.
####################################################################
cdef class ModN_Reduction:
    cdef ResidueRingModN S
    cdef modn_matrix[4] G
    def __init__(self, N):
        cdef ResidueRingModN S = ResidueRingModN(N)
        self.S = S
        if N.norm() % 2 == 0:
            # TODO: finish this
            raise NotImplementedError, "need to implement mod-N reduction for even ideals!"

        # I = [0,i2, 1, 0] works.
        F = N.number_field()
        S.coerce_from_nf(self.G[1][0], F(0))
        S.coerce_from_nf(self.G[1][1], F(-1))
        S.coerce_from_nf(self.G[1][2], F(1))
        S.coerce_from_nf(self.G[1][3], F(0))

        # Now find J.
        cdef modn_element a, b, c, d, t, t2, minus_one
        found_it = False
        S.set(b, 1)
        S.set(minus_one, -1)
        while not S.is_last_element(b):
            # Let c = -1 - b*b
            S.mul(t, b, b)
            S.mul(t, t, minus_one)
            S.add(c, minus_one, t)
            if S.is_square(c):
                # Next set a = -sqrt(c).
                S.sqrt(a, c)  
                found_it = True
                break
            S.next_element(b, b)

        if not found_it:  # sometimes needed
            raise NotImplementedError

        # Set the matrix self.G[2] to [a,b,(j2-a*a)/b,-a]
        S.set_element(self.G[2][0], a)
        S.set_element(self.G[2][1], b)
        # Set t to (-1-a*a)/b
        S.mul(t, a, a)
        S.mul(t, t, minus_one)
        S.add(t, t, minus_one)
        S.inv(t2, b)
        S.mul(self.G[2][2], t, t2)
        S.mul(self.G[2][3], a, minus_one)

        S.matrix_mul(self.G[3], self.G[1], self.G[2])

        # Finally, set the identity matrix
        S.set(self.G[0][0], 1); S.set(self.G[0][1], 0); S.set(self.G[0][2], 0); S.set(self.G[0][3], 1)

    def __repr__(self):
        return 'Reduction modulo %s of norm %s defined by sending I to %s and J to %s'%(
            self.S.N._repr_short(), self.S.N.norm(),
            self.S.matrix_to_str(self.G[1]), self.S.matrix_to_str(self.G[2]))

    cdef int quatalg_to_modn_matrix(self, modn_matrix M, alpha) except -1:
        # Given an element alpha in the quaternion algebra, find its image M as a modn_matrix.
        cdef modn_element t
        cdef modn_element[4] X
        cdef int i, j
        cdef ResidueRingModN S = self.S
        for i in range(4):
            S.coerce_from_nf(X[i], alpha[i])
        for i in range(4):
            S.mul(M[i], self.G[0][i], X[0])
            for j in range(1, 4):
                S.mul(t, self.G[j][i], X[j])
                S.add(M[i], M[i], t)
                

####################################################################
# R^* \ P1(O_F/N)
####################################################################

ctypedef modn_matrix icosian_matrices[120]

cdef class IcosiansModP1ModN:
    cdef icosian_matrices G
    cdef ModN_Reduction f
    cdef ProjectiveLineModN P1
    cdef long *std_to_rep_table
    cdef long *orbit_reps
    cdef long _cardinality
    cdef p1_element* orbit_reps_p1elt
    def __init__(self, N):
        # compute choice of splitting
        self.f = ModN_Reduction(N)
        self.P1 = ProjectiveLineModN(N)
        self.orbit_reps = <long*>0
        self.std_to_rep_table = <long*> sage_malloc(sizeof(long) * self.P1.cardinality())
        
        # This is very wasteful, by a factor of about 120 on average.  REDO.
        self.orbit_reps_p1elt = <p1_element*>sage_malloc(sizeof(p1_element) * self.P1.cardinality())
        
        # initialize the group G of the 120 mod-N icosian matrices
        from sqrt5 import all_icosians
        X = all_icosians()
        cdef int i
        for i in range(len(X)):
            self.f.quatalg_to_modn_matrix(self.G[i], X[i])

        self.compute_std_to_rep_table()

    def __dealloc__(self):
        sage_free(self.std_to_rep_table)
        sage_free(self.orbit_reps_p1elt)
        if self.orbit_reps:
            sage_free(self.orbit_reps)

    def __repr__(self):
        return "The %s orbits for the action of the Icosian group on %s"%(self._cardinality, self.P1)

    cpdef compute_std_to_rep_table(self):
        # Algorithm is pretty obvious.  Just take first element of P1,
        # hit by all icosians, making a list of those that are
        # equivalent to it.  Then find next element of P1 not in the
        # first list, etc.

        cdef p1_element x, Gx
        self.P1.first_element(x)
        cdef long ind=0, j, i=0
        for j in range(self.P1._cardinality):
            self.std_to_rep_table[j] = -1
        self.std_to_rep_table[0] = 0
        reps = []
        while ind < self.P1._cardinality:
            reps.append(ind)
            self.P1.set_element(self.orbit_reps_p1elt[i], x)
            for j in range(120):
                self.P1.matrix_action(Gx, self.G[j], x)
                self.P1.reduce_element(Gx, Gx)
                self.std_to_rep_table[self.P1.standard_index(Gx)] = i
            while self.std_to_rep_table[ind] != -1 and ind < self.P1._cardinality:
                ind += 1
                if ind < self.P1._cardinality:
                    self.P1.next_element(x, x)
            i += 1
        self._cardinality = len(reps)
        self.orbit_reps = <long*> sage_malloc(sizeof(long)*self._cardinality)
        for j in range(self._cardinality):
            self.orbit_reps[j] = reps[j]

    def compute_std_to_rep_table_debug(self):
        cdef p1_element x, Gx
        self.P1.first_element(x)
        cdef long ind=0, j, i=0, k
        for j in range(self.P1._cardinality):
            self.std_to_rep_table[j] = -1
        self.std_to_rep_table[0] = 0
        reps = []
        orbits = []
        while ind < self.P1._cardinality:
            reps.append(ind)
            print "Found representative number %s: %s"%(i, self.P1.element_to_str(x))
            self.P1.set_element(self.orbit_reps_p1elt[i], x)
            orbit = [ind]
            for j in range(120):
                self.P1.matrix_action(Gx, self.G[j], x)
                self.P1.reduce_element(Gx, Gx)
                k = self.P1.standard_index(Gx)
                orbit.append(k)
                self.std_to_rep_table[k] = i
            orbit = list(sorted(set(orbit)))
            orbits.append(orbit)
            print "It has an orbit of size %s"%len(orbit)
            while self.std_to_rep_table[ind] != -1 and ind < self.P1._cardinality:
                ind += 1
                if ind < self.P1._cardinality:
                    self.P1.next_element(x, x)
            i += 1
        self._cardinality = len(reps)
        self.orbit_reps = <long*> sage_malloc(sizeof(long)*self._cardinality)
        for j in range(self._cardinality):
            self.orbit_reps[j] = reps[j]
        return orbits

    cpdef long cardinality(self):
        return self._cardinality

    def hecke_matrix(self, P):
        """
        Return matrix of Hecke action.
        """
        # TODO: deal with when P divides the level...
        
        cdef modn_matrix M
        cdef p1_element Mx
        from sqrt5 import hecke_elements
        from sage.all import zero_matrix, ZZ
        T = zero_matrix(ZZ, self._cardinality)
        cdef long i, j
        
        for alpha in hecke_elements(P):
            # TODO: probably better to invert after???  not at all?
            self.f.quatalg_to_modn_matrix(M, alpha)
            for i in range(self._cardinality):
                self.P1.matrix_action(Mx, M, self.orbit_reps_p1elt[i])
                self.P1.reduce_element(Mx, Mx)
                j = self.std_to_rep_table[self.P1.standard_index(Mx)]
                T[i,j] = T[i,j] + 1

        return T
        


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

