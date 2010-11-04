"""
Hilbert modular forms over F = Q(sqrt(5)).

All the code in this file is meant to be highly optimized. 
"""

include 'stdsage.pxi'

def roots_mod_pn(f, p, n):
    """
    Return the roots of the integer polynomial f modulo p**n, as
    integers.

    NOTE: Though Sage has a roots function that does this, it is dog
    slow if n is at all large.  However, one can do the same quickly
    by factoring over the p-adics (which uses PARI).
    """
    return [(-a[0][0]).lift() for a in f.factor_padic(p, prec=n+1)]

from sage.rings.ring cimport CommutativeRing
from sage.rings.all import Integers

# Residue element
ctypedef int residue_element[2]

#from sqrt5_fast_python import ResidueRing

cpdef int ideal_characteristic(I):
    return I.number_field()._pari_().idealtwoelt(I._pari_())[0]

def ResidueRing(P, e):
    cdef int p = ideal_characteristic(P)
    if p == 5:   # ramified
        if e % 2 == 0:
            return ResidueRing_nonsplit(P, p, e)
        else:
            return ResidueRing_ramified_odd(P, p, e)
    elif p%5 in [2,3]:  # nonsplit
        return ResidueRing_nonsplit(P, p, e)
    else: # split
        return ResidueRing_split(P, p, e)

cdef class ResidueRing_abstract(CommutativeRing):
    cdef object P, e, F
    cdef public object element_class
    cdef int n0, n1, p
    cdef int im_gen0
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
    
    def __call__(self, x):
        if hasattr(x, 'parent'):
            P = x.parent()
            if P is self:
                return x
            elif P is self.F:
                # TODO: This denominator is very slow.
                d = x.denominator()
                if d != 1:
                    return self(d*x) * self(d.inverse_mod(self.p**self.e))
                return self.element_class(self, x)
        raise TypeError

    def __repr__(self):
        return "Residue class ring of %s^%s of characteristic %s"%(
            self.P._repr_short(), self.e, self.p)

    cdef add(self, residue_element rop, residue_element op0, residue_element op1):
        raise NotImplementedError
    cdef sub(self, residue_element rop, residue_element op0, residue_element op1):
        raise NotImplementedError        
    cdef void mul(self, residue_element rop, residue_element op0, residue_element op1):
        raise NotImplementedError
    
cdef class ResidueRing_split(ResidueRing_abstract):
    def __init__(self, P, p, e):
        ResidueRing_abstract.__init__(self, P, p, e)
        self.element_class = ResidueRingElement_split
        self.n0 = p**e
        self.n1 = 1
        # Compute the mapping by finding roots mod p^e.
        alpha = P.number_field().gen()
        for z in roots_mod_pn(self.F.defining_polynomial(), p, e):
            if z - alpha in P:
                self.im_gen0 = z
        assert self.im_gen0 is not None

    cdef add(self, residue_element rop, residue_element op0, residue_element op1):
        rop[0] = (op0[0] + op1[0])%self.n0
        
    cdef sub(self, residue_element rop, residue_element op0, residue_element op1):
        rop[0] = (op0[0] - op1[0])%self.n0
        
    cdef void mul(self, residue_element rop, residue_element op0, residue_element op1):
        rop[0] = (op0[0] * op1[0])%self.n0

        
cdef class ResidueRing_nonsplit(ResidueRing_abstract):
    def __init__(self, P, p, e):
        ResidueRing_abstract.__init__(self, P, p, e)
        self.element_class = ResidueRingElement_nonsplit
        self.n0 = p**e
        self.n1 = p**e

    cdef add(self, residue_element rop, residue_element op0, residue_element op1):
        rop[0] = (op0[0] + op1[0])%self.n0
        rop[1] = (op0[1] + op1[1])%self.n1
        
    cdef sub(self, residue_element rop, residue_element op0, residue_element op1):
        # cython uses python mod normalization convention, so no need to do that manually
        rop[0] = (op0[0] - op1[0])%self.n0
        rop[1] = (op0[1] - op1[1])%self.n1

    cdef void mul(self, residue_element rop, residue_element op0, residue_element op1):
        # it is significantly faster doing these arrays accesses only once in the line below!
        cdef int a = op0[0], b = op0[1], c = op1[0], d = op1[1]
        rop[0] = (a*c + b*d)%self.n0
        rop[1] = (a*d + b*d + b*c)%self.n1

cdef class ResidueRing_ramified_odd(ResidueRing_abstract):
    cdef int two_inv
    def __init__(self, P, p, e):
        ResidueRing_abstract.__init__(self, P, p, e)        
        # This is assumed in the code for the element_class so it better be true!
        assert self.F.defining_polynomial().list() == [-1,-1,1]
        self.n0 = p**(e//2+1)
        self.n1 = p**(e//2)
        self.two_inv = (Integers(self.n0)(2)**(-1)).lift()
        self.element_class = ResidueRingElement_ramified_odd

    cdef add(self, residue_element rop, residue_element op0, residue_element op1):
        rop[0] = (op0[0] + op1[0])%self.n0
        rop[1] = (op0[1] + op1[1])%self.n1
        
    cdef sub(self, residue_element rop, residue_element op0, residue_element op1):
        rop[0] = (op0[0] - op1[0])%self.n0
        rop[1] = (op0[1] - op1[1])%self.n1
        
    cdef void mul(self, residue_element rop, residue_element op0, residue_element op1):
        cdef int a = op0[0], b = op0[1], c = op1[0], d = op1[1]
        rop[0] = (a*c + 5*b*d)%self.n0
        rop[1] = (a*d + b*c)%self.n1
    

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

cdef class ResidueRingElement_split(ResidueRingElement):
    def __init__(self, ResidueRing_split parent, x):
        self._parent = parent
        assert x.parent() is parent.F
        v = x._coefficients()
        if len(v) == 0:
            self.x[0] = 0
            return
        elif len(v) == 1:
            self.x[0] = v[0]
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
            self.x[0] = 0; self.x1 = 0
            return
        elif len(v) == 1:
            self.x[0] = v[0]
            self.x1 = 0
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
        return '%s + %s*gamma_bar'%(self.x[0], self.x[1])

    def lift(self):
        """
        Return lift of element to number field F.
        """
        # TODO: test...
        return self._parent.F([self.x0, self.x1])

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
        cdef int a, b
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
        return '%s + %s*sqrt5'%(self.x[0], self.x[1])




####################################################################
# fragments/test code below
####################################################################

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

