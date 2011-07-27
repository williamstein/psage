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
    [2a, 5a, 3a, 11a, 11b, 19a, 19b, 29a, 29b, 31a, 31b, 41a, 41b, 7a]

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


include "stdsage.pxi"
include "interrupt.pxi"

cdef extern from "pari/pari.h":
    unsigned long Fl_sqrt(unsigned long a, unsigned long p)
    unsigned long Fl_div(unsigned long a, unsigned long b, unsigned long p)

cdef class Prime:
    """
    Nonzero prime ideal of the ring of integers of Q(sqrt(5)).  This
    is a fast customized Cython class; to get at the corresponding
    Sage prime ideal use the sage_ideal method.
    """
    def __init__(self, long p, long r, bint first):
        """
        Create Prime ideal with given residue characteristic, root,
        and first or not with that characterstic.

        INPUT:
            - `p` -- prime
            - `r` -- root (or 0)
            - ``first`` -- boolean: True if first prime over p

        NOTE: No checking is done to verify that the input is valid.
        
        EXAMPLES::

            sage: from psage.number_fields.sqrt5.prime import Prime
            sage: P = Prime(2,0,True); P
            2a
            sage: type(P)
            <type 'psage.number_fields.sqrt5.prime.Prime'>
            sage: P = Prime(5,-2,True); P
            5a
            sage: P = Prime(11,-3,True); P
            11a
            sage: P = Prime(11,4,False); P
            11b
        """
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
        """
        return '%s%s'%(self.p, 'a' if self.first else 'b')

    def _latex_(self):
        """
        EXAMPLES::

            sage: from psage.number_fields.sqrt5.prime import Prime
            sage: Prime(11,-3,True)._latex_()
            '11a'
            sage: Prime(11,4,False)._latex_()
            '11b'
        """
        return self.__repr__()

    cpdef long norm(self):
        """
        Return the norm of this ideal.

        EXAMPLES::

            sage: from psage.number_fields.sqrt5.prime import Prime
            sage: Prime(11,3,True).norm()
            11
            sage: Prime(7,0,True).norm()
            49
        """
        return self.p*self.p if self.r == 0 else self.p

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
            [2a, 5a, 3a, 11a, 11b, 19a, 19b, 29a, 29b, 31a, 31b, 41a, 41b, 7a]
            sage: v[3], v[4]
            (11a, 11b)
            sage: v[3] < v[4]
            True
            sage: v[4] > v[3]
            True
        

        We test the ordering a bit by sorting::
        
            sage: v.sort(); v
            [2a, 5a, 3a, 11a, 11b, 19a, 19b, 29a, 29b, 31a, 31b, 41a, 41b, 7a]
            sage: v = list(reversed(v)); v
            [7a, 41b, 41a, 31b, 31a, 29b, 29a, 19b, 19a, 11b, 11a, 3a, 5a, 2a]
            sage: v.sort(); v
            [2a, 5a, 3a, 11a, 11b, 19a, 19b, 29a, 29b, 31a, 31b, 41a, 41b, 7a]

        A bigger test::

            sage: v = primes_of_bounded_norm(10^7)
            sage: w = list(reversed(v)); w.sort()
            sage: v == w
            True
        """
        cdef int c
        c = cmp(self.norm(), right.norm())
        if c: return c
        return cmp(self.r, right.r)
    
    def sage_ideal(self):
        """
        Return the usual prime fractional ideal associated to this
        prime.  This is slow, but provides substantial additional
        functionality.

        EXAMPLES::

            sage: from psage.number_fields.sqrt5 import primes_of_bounded_norm
            sage: v = primes_of_bounded_norm(20)
            sage: v[1].sage_ideal()
            Fractional ideal (a + 2)
            sage: [P.sage_ideal() for P in v]
            [Fractional ideal (2), Fractional ideal (a + 2), Fractional ideal (3), Fractional ideal (3*a - 1), Fractional ideal (3*a - 2), Fractional ideal (-4*a + 1), Fractional ideal (-4*a + 3)]
        """
        from misc import F
        cdef long p=self.p, r=self.r
        if p == 5:   # ramified
            return F.ideal(F.gen()-r)
        elif r == 0: # inert case
            return F.ideal([p])
        else:        # split case
            return F.ideal([p, F.gen()-r])

def primes_of_bounded_norm(bound):
    """
    Return ordered list of all prime ideals of the ring of integers of
    Q(sqrt(5)) of norm less than bound.

    The primes are instances of a special fast Primes class.  They are
    sorted first by norm, then in the remaining split case by the
    integer in the interval [0,p) congruent to (1+sqrt(5))/2.  For
    optimal speed you can use the Prime objects directly from Cython,
    which provides direct C-level access to the underlying data
    structure.


    INPUT:
        - ``bound`` -- nonnegative integer, less than `2^31`

    WARNING: The ordering is NOT the same as the ordering of primes by
    Sage.   Even if you order first by norm, then use Sage's ordering
    for primes of the same norm, then the orderings do not agree.::

    EXAMPLES::

        sage: from psage.number_fields.sqrt5 import primes_of_bounded_norm
        sage: primes_of_bounded_norm(0)
        []
        sage: primes_of_bounded_norm(10)
        [2a, 5a, 3a]
        sage: primes_of_bounded_norm(50)
        [2a, 5a, 3a, 11a, 11b, 19a, 19b, 29a, 29b, 31a, 31b, 41a, 41b, 7a]
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
    if bound < 0:
        return []
    if bound >= 2**31:
        raise ValueError, "bound must be less than 2^31"

    cdef long p, i = 0, sr, r0, r1
    cdef Prime P
    cdef list v = []

    from sage.all import prime_range
    
    for p in prime_range(bound):
        t = p % 5
        if t == 1 or t == 4:   # split
            # Compute a square root of 5 modulo p.
            sr = Fl_sqrt(5, p)
            # Find the two values of (1+sqrt(5))/2.
            r0 = Fl_div(1+sr, 2, p)
            r1 = Fl_div(1+p-sr, 2, p)
            if r0 > r1: # swap
                r0, r1 = r1, r0
            # Append each prime to the list
            P = PY_NEW(Prime); P.p = p; P.r = r0; P.first = True; v.append(P)
            P = PY_NEW(Prime); P.p = p; P.r = r1; P.first = False; v.append(P)
        elif p == 5:  # ramified
            v.append(Prime(p, -2, True))
        elif p*p < bound:  # inert
            v.append(Prime(p, 0, True))
    v.sort()
    return v

    



