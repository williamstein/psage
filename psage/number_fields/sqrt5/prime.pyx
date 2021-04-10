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
Fast prime ideals of the ring R of integers of Q(sqrt(5)).

This module implements Cython classes for prime ideals of R, and their
enumeration.  The main entry function is primes_of_bounded_norm::

    sage: from psage.number_fields.sqrt5 import primes_of_bounded_norm
    sage: v = primes_of_bounded_norm(50); v
    [2, 5a, 3, 11a, 11b, 19a, 19b, 29a, 29b, 31a, 31b, 41a, 41b, 7]

Note that the output of primes_of_bounded_norm is a list.  Each entry
is a prime ideal, which prints using a simple label consisting of the
characteristic of the prime then "a" or "b", where "b" only appears for
the second split prime.::

    sage: type(v[8])
    <type 'psage.number_fields.sqrt5.prime.Prime'>
    sage: v[8].sage_ideal()
    Fractional ideal (a + 5)

AUTHOR:
    - William Stein 
"""


from cysignals.memory cimport sig_free,sig_malloc
from cysignals.signals cimport sig_on,sig_off
from sage.ext.stdsage cimport PY_NEW

cdef extern from "pari/pari.h":
    unsigned long Fl_sqrt(unsigned long, unsigned long)
    unsigned long Fl_div(unsigned long, unsigned long, unsigned long)

from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal

cdef class Prime:
    """
    Nonzero prime ideal of the ring of integers of Q(sqrt(5)).  This
    is a fast customized Cython class; to get at the corresponding
    Sage prime ideal use the sage_ideal method.
    """
    def __init__(self, p, long r=0, bint first=True):
        """
        Create Prime ideal with given residue characteristic, root,
        and first or not with that characterstic.

        INPUT form 1:
            - `p` -- prime
            - `r` -- root (or 0)
            - ``first`` -- boolean: True if first prime over p

        INPUT form 2:
            - p -- prime ideal of integers of Q(sqrt(5)); validity of the
              input is not checked in any way!

        NOTE: No checking is done to verify that the input is valid.
        
        EXAMPLES::

            sage: from psage.number_fields.sqrt5.prime import Prime
            sage: P = Prime(2,0,True); P
            2
            sage: type(P)
            <type 'psage.number_fields.sqrt5.prime.Prime'>
            sage: P = Prime(5,3,True); P
            5a
            sage: P = Prime(11,8,True); P
            11a
            sage: P = Prime(11,4,False); P
            11b

        We can also set P using a prime ideal of the ring of integers::

            sage: from psage.number_fields.sqrt5.prime import Prime; K.<a> = NumberField(x^2-x-1)
            sage: P1 = Prime(K.primes_above(11)[0]); P2 = Prime(K.primes_above(11)[1]); P1, P2
            (11b, 11a)
            sage: P1 > P2
            True
            sage: Prime(K.prime_above(2))
            2
            sage: P = Prime(K.prime_above(5)); P, P.r
            (5a, 3)
            sage: Prime(K.prime_above(3))
            3

        Test that conversion both ways works for primes up to norm `10^5`::
        
            sage: from psage.number_fields.sqrt5.prime import primes_of_bounded_norm, Prime
            sage: v = primes_of_bounded_norm(10^5)
            sage: w = [Prime(z.sage_ideal()) for z in v]
            sage: v == w
            True
        """
        cdef long t, r1
        if isinstance(p, NumberFieldFractionalIdeal):
            # Set self using a prime ideal of Q(sqrt(5)).
            H = p.pari_hnf()
            self.p = H[0,0]
            self.first = True
            t = self.p % 5
            if t == 1 or t == 4:
                self.r = self.p - H[0,1]
                r1 = self.p + 1 - self.r
                if self.r > r1:
                    self.first = False
            elif t == 0:
                self.r = 3
            else:
                self.r = 0
        else:
            self.p = p
            self.r = r
            self.first = first
        
    def __repr__(self):
        """
        EXAMPLES::
        
            sage: from psage.number_fields.sqrt5.prime import Prime
            sage: Prime(11,4,False).__repr__()
            '11b'
            sage: Prime(11,4,True).__repr__()
            '11a'
            sage: Prime(7,0,True).__repr__()
            '7'
        """
        if self.r:
            return '%s%s'%(self.p, 'a' if self.first else 'b')
        return '%s'%self.p

    def __hash__(self):
        return self.p*(self.r+1) + int(self.first)

    def _latex_(self):
        """
        EXAMPLES::

            sage: from psage.number_fields.sqrt5.prime import Prime
            sage: Prime(11,8,True)._latex_()
            '11a'
            sage: Prime(11,4,False)._latex_()
            '11b'
        """
        return self.__repr__()

    cpdef bint is_split(self):
        """
        Return True if this prime is split (and not ramified).
        
        EXAMPLES::

            sage: from psage.number_fields.sqrt5.prime import Prime
            sage: Prime(11,8,True).is_split()
            True
            sage: Prime(3,0,True).is_split()
            False
            sage: Prime(5,3,True).is_split()
            False
        """
        return self.r != 0 and self.p != 5
    
    cpdef bint is_inert(self):
        """
        Return True if this prime is inert. 

        EXAMPLES::

            sage: from psage.number_fields.sqrt5.prime import Prime
            sage: Prime(11,8,True).is_inert()
            False
            sage: Prime(3,0,True).is_inert()
            True
            sage: Prime(5,3,True).is_inert()
            False
        """
        return self.r == 0
    
    cpdef bint is_ramified(self):
        """
        Return True if this prime is ramified (i.e., the prime over 5).
                
        EXAMPLES::

            sage: from psage.number_fields.sqrt5.prime import Prime
            sage: Prime(11,8,True).is_ramified()
            False
            sage: Prime(3,0,True).is_ramified()
            False
            sage: Prime(5,3,True).is_ramified()
            True
        """
        return self.p == 5

    cpdef long norm(self):
        """
        Return the norm of this ideal.

        EXAMPLES::

            sage: from psage.number_fields.sqrt5.prime import Prime
            sage: Prime(11,4,True).norm()
            11
            sage: Prime(7,0,True).norm()
            49
        """
        return self.p if self.r else self.p*self.p

    def __cmp__(self, Prime right):
        """
        Compare two prime ideals. First sort by the norm, then in the
        (only remaining) split case if the norms are the same, compare
        the residue of (1+sqrt(5))/2 in the interval [0,p).

        WARNING: The ordering is NOT the same as the ordering of
        fractional ideals in Sage.

        EXAMPLES::

            sage: from psage.number_fields.sqrt5 import primes_of_bounded_norm
            sage: v = primes_of_bounded_norm(50); v
            [2, 5a, 3, 11a, 11b, 19a, 19b, 29a, 29b, 31a, 31b, 41a, 41b, 7]
            sage: v[3], v[4]
            (11a, 11b)
            sage: v[3] < v[4]
            True
            sage: v[4] > v[3]
            True

        We test the ordering a bit by sorting::
        
            sage: v.sort(); v
            [2, 5a, 3, 11a, 11b, 19a, 19b, 29a, 29b, 31a, 31b, 41a, 41b, 7]
            sage: v = list(reversed(v)); v
            [7, 41b, 41a, 31b, 31a, 29b, 29a, 19b, 19a, 11b, 11a, 3, 5a, 2]
            sage: v.sort(); v
            [2, 5a, 3, 11a, 11b, 19a, 19b, 29a, 29b, 31a, 31b, 41a, 41b, 7]

        A bigger test::

            sage: v = primes_of_bounded_norm(10^7)
            sage: w = list(reversed(v)); w.sort()
            sage: v == w
            True
        """
        cdef long selfn = self.norm(),  rightn = right.norm()
        if   selfn > rightn: return 1
        elif rightn > selfn: return -1
        elif self.r > right.r: return 1
        elif right.r > self.r: return -1
        else: return 0

    def sage_ideal(self):
        """
        Return the usual prime fractional ideal associated to this
        prime.  This is slow, but provides substantial additional
        functionality.

        EXAMPLES::

            sage: from psage.number_fields.sqrt5 import primes_of_bounded_norm
            sage: v = primes_of_bounded_norm(20)
            sage: v[1].sage_ideal()
            Fractional ideal (2*a - 1)
            sage: [P.sage_ideal() for P in v]
            [Fractional ideal (2), Fractional ideal (2*a - 1), Fractional ideal (3), Fractional ideal (3*a - 1), Fractional ideal (3*a - 2), Fractional ideal (-4*a + 1), Fractional ideal (-4*a + 3)]
        """
        from misc import F
        cdef long p=self.p, r=self.r
        if r: # split and ramified cases
            return F.ideal(p, F.gen()-r)
        else: # inert case
            return F.ideal(p)

from sage.rings.integer import Integer

def primes_above(long p, bint check=True):
    """
    Return ordered list of all primes above p in the ring of integers
    of Q(sqrt(5)).  See the docstring for primes_of_bounded_norm.

    INPUT:
        - p -- prime number in integers ZZ (less than `2^{31}`)
        - check -- bool (default: True); if True, check that p is prime
    OUTPUT:
        - list of 1 or 2 Prime objects

    EXAMPLES::

        sage: from psage.number_fields.sqrt5.prime import primes_above
        sage: primes_above(2)
        [2]
        sage: primes_above(3)
        [3]
        sage: primes_above(5)
        [5a]
        sage: primes_above(11)
        [11a, 11b]
        sage: primes_above(13)
        [13]
        sage: primes_above(17)
        [17]
        sage: primes_above(4)
        Traceback (most recent call last):
        ...
        ValueError: p must be a prime
        sage: primes_above(4, check=False)
        [2]
    """
    if check and not Integer(p).is_pseudoprime():
        raise ValueError("p must be a prime")
    cdef long t = p%5
    if t == 1 or t == 4 or t == 0:
        return prime_range(p, p+1)
    else: # inert
        return [Prime(p, 0, True)]

def primes_of_bounded_norm(bound):
    """
    Return ordered list of all prime ideals of the ring of integers of
    Q(sqrt(5)) of norm less than bound.

    The primes are instances of a special fast Primes class (they are
    *not* usual Sage prime ideals -- use the sage_ideal() method to
    get those).  They are sorted first by norm, then in the remaining
    split case by the integer in the interval [0,p) congruent to
    (1+sqrt(5))/2.  

    INPUT:
        - ``bound`` -- nonnegative integer, less than `2^{31}`

    WARNING: The ordering is NOT the same as the ordering of primes by
    Sage.   Even if you order first by norm, then use Sage's ordering
    for primes of the same norm, then the orderings do not agree.::

    EXAMPLES::

        sage: from psage.number_fields.sqrt5 import primes_of_bounded_norm
        sage: primes_of_bounded_norm(0)
        []
        sage: primes_of_bounded_norm(10)
        [2, 5a, 3]
        sage: primes_of_bounded_norm(50)
        [2, 5a, 3, 11a, 11b, 19a, 19b, 29a, 29b, 31a, 31b, 41a, 41b, 7]
        sage: len(primes_of_bounded_norm(10^6))
        78510

    We grab one of the primes::
    
        sage: v = primes_of_bounded_norm(100)
        sage: P = v[3]; type(P)
        <type 'psage.number_fields.sqrt5.prime.Prime'>

    It prints with a nice label::
    
        sage: P
        11a

    You can get the corresponding fractional ideal as a normal Sage ideal::
    
        sage: P.sage_ideal()
        Fractional ideal (3*a - 1)

    You can also get the underlying residue characteristic::
    
        sage: P.p
        11

    And, the image of (1+sqrt(5))/2 modulo the prime (or 0 in the inert case)::
    
        sage: P.r
        4
        sage: z = P.sage_ideal(); z.residue_field()(z.number_field().gen())
        4

    Prime enumeration is reasonable fast, even when the input is
    relatively large (going up to `10^8` takes a few seconds, and up
    to `10^9` takes a few minutes), and the following should take less
    than a second::

        sage: len(primes_of_bounded_norm(10^7))  # less than a second
        664500

    One limitation is that the bound must be less than `2^{31}`::

        sage: primes_of_bounded_norm(2^31)
        Traceback (most recent call last):
        ...
        ValueError: bound must be less than 2^31
    """
    return prime_range(bound)

def prime_range(long start, stop=None):
    """
    Return ordered list of all prime ideals of the ring of integers of
    Q(sqrt(5)) of norm at least start and less than stop.  If only
    start is given then return primes with norm less than start.

    The primes are instances of a special fast Primes class (they are
    *not* usual Sage prime ideals -- use the sage_ideal() method to
    get those).  They are sorted first by norm, then in the remaining
    split case by the integer in the interval [0,p) congruent to
    (1+sqrt(5))/2.  For optimal speed you can use the Prime objects
    directly from Cython, which provides direct C-level access to the
    underlying data structure.

    INPUT:
        - ``start`` -- nonnegative integer, less than `2^{31}`
        - ``stop``  -- None or nonnegative integer, less than `2^{31}`

    EXAMPLES::

        sage: from psage.number_fields.sqrt5.prime import prime_range
        sage: prime_range(10, 60)
        [11a, 11b, 19a, 19b, 29a, 29b, 31a, 31b, 41a, 41b, 7, 59a, 59b]
        sage: prime_range(2, 11)
        [2, 5a, 3]
        sage: prime_range(2, 12)
        [2, 5a, 3, 11a, 11b]
        sage: prime_range(3, 12)
        [2, 5a, 3, 11a, 11b]
        sage: prime_range(9, 12)
        [3, 11a, 11b]
        sage: prime_range(5, 12)
        [5a, 3, 11a, 11b]
    """
    if start >= 2**31 or (stop and stop >= 2**31):
        raise ValueError("bound must be less than 2^31")

    cdef long p, p2, sr, r0, r1, t, bound
    cdef Prime P
    cdef list v = []

    if stop is None:
        bound = start
        start = 2
    else:
        bound = stop
    
    from sage.all import prime_range as prime_range_ZZ
    
    for p in prime_range_ZZ(bound, py_ints=True):
        t = p % 5
        if t == 1 or t == 4:   # split
            if p >= start:
                # Compute a square root of 5 modulo p.
                sr = Fl_sqrt(5, p)
                # Find the two values of (1+sqrt(5))/2.
                r0 = Fl_div(1+sr, 2, p)
                r1 = p+1-r0
                # Sort
                if r0 > r1: r0, r1 = r1, r0
                # Append each prime to the list
                P = PY_NEW(Prime); P.p = p; P.r = r0; P.first = True; v.append(P)
                P = PY_NEW(Prime); P.p = p; P.r = r1; P.first = False; v.append(P)
        elif p == 5:   # ramified
            if p >= start:  
                v.append(Prime(p, 3, True))
        else:
            p2 = p*p
            if p2 < bound and p2 >= start:
                v.append(Prime(p, 0, True))
                
    v.sort()
    return v

    



