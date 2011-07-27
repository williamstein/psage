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
    sage: v[8].sage_prime()
    Fractional ideal (-a + 6)

This module also includes an object PrimesOfBoundedNorm, which is used
behind the scenes to implement the function primes_of_bounded_norm.
The algorithm used for prime enumeration is straightforward: list each
rational prime, and for the split primes, find the two square roots of
5 modulo p.  Then sort the result.  To make the code fast, we use
simple data structures, and a custom insertion sort that uses memmove
and inserts the inert primes in their proper places.  The code works
for listing the primes up to `10^7` in less then a second, and up to
`10^9` in a few minutes.  Here is a simple example of the
PrimesOfBoundedNorm object::

    sage: import psage.number_fields.sqrt5.prime
    sage: P = psage.number_fields.sqrt5.prime.PrimesOfBoundedNorm(1000)
    sage: P
    Prime ideals of the ring of integers of Q(sqrt(5)) of norm less than 1000
    sage: P[7]
    29a
    sage: len(P)
    163    

TESTS::

    sage: from psage.number_fields.sqrt5 import primes_of_bounded_norm
    sage: v = primes_of_bounded_norm(10^6)
    sage: w = list(v)
    sage: w.sort()
    sage: w == v
    True

AUTHOR:
    - William Stein 
"""


include "stdsage.pxi"
include "interrupt.pxi"

cdef extern from "pari/pari.h":
    unsigned long Fl_sqrt(unsigned long a, unsigned long p)
    unsigned long Fl_div(unsigned long a, unsigned long b, unsigned long p)

cdef extern from "string.h":
    void * memmove(void *s1, void *s2, size_t n)
    void * memcpy(void *s1, void *s2, size_t n)

from misc import F

cdef class Prime:
    """
    Nonzero prime ideal of the ring of integers of Q(sqrt(5)).  This
    is a fast customized Cython class; to get at the corresponding
    Sage prime ideal use the sage_prime method.
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
        the fixed residue of (1+sqrt(5))/2 in the interval (-p/2,p/2).

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
    
    def sage_prime(self):
        """
        Return the usual prime fractional ideal associated to this
        prime.  This is slow, but provides substantial additional
        functionality.

        EXAMPLES::

            sage: from psage.number_fields.sqrt5 import primes_of_bounded_norm
            sage: v = primes_of_bounded_norm(20)
            sage: v[1].sage_prime()
            Fractional ideal (a + 2)
            sage: [P.sage_prime() for P in v]
            [Fractional ideal (2), Fractional ideal (a + 2), Fractional ideal (3), Fractional ideal (3*a - 2), Fractional ideal (3*a - 1), Fractional ideal (-4*a + 3), Fractional ideal (-4*a + 1)]
        """
        cdef long p=self.p, r=self.r
        if p == 5:   # ramified
            return F.ideal(F.gen()-r)
        elif r == 0: # inert case
            return F.ideal([p])
        else:        # split case
            return F.ideal([p, F.gen()-r])

cdef class PrimesOfBoundedNorm:
    """
    Return object representing the prime ideals of Q(sqrt(5)) of norm
    less than a given bound, sorted first by norm, then in the split
    case by the integer in (-p/2,p/2) congruent to (1+sqrt(5))/2.

    EXAMPLES::

        sage: from psage.number_fields.sqrt5.prime import PrimesOfBoundedNorm
        sage: v = PrimesOfBoundedNorm(20); v
        Prime ideals of the ring of integers of Q(sqrt(5)) of norm less than 20
        sage: v.bound
        20
        sage: len(v)
        7
        sage: v[0]
        2a

    You can make a list of the primes as Prime objects::
    
        sage: list(v)
        [2a, 5a, 3a, 11a, 11b, 19a, 19b]

    WARNING: The ordering is NOT the same as the ordering of primes by
    Sage. (Even if you order first by norm, then use Sage's ordering
    for primes of the same norm, then the orderings do not agree.)::

        sage: w = list(v)
        sage: w.sort()
        sage: w == list(v)
        True
    
    For optimal speed you can use the PrimesOfBoundedNorm object
    directly from Cython, which provides direct C-level access to the
    underlying information in the tuples.

    Prime enumeration is reasonable fast, even when the input is
    relatively large (going up to `10^8` takes a few seconds, and up
    to `10^9` takes a few minutes), and the following should take less
    than a second::

        sage: len(PrimesOfBoundedNorm(10^7))  # less than a second
        664500

    One limitation is that the bound must be less than `2^{31}`::

        sage: PrimesOfBoundedNorm(2^31)
        Traceback (most recent call last):
        ...
        NotImplementedError: bound must be less than 2^31
    """
    def __cinit__(self):
        """
        Initialize memory.

        EXAMPLES::

            sage: from psage.number_fields.sqrt5.prime import PrimesOfBoundedNorm
            sage: PrimesOfBoundedNorm(50)   # indirect test
            Prime ideals of the ring of integers of Q(sqrt(5)) of norm less than 50
        """
        # Make sure these are NULL asap, so just in case __dealloc__ is
        # called it won't segfault.
        self.prime = NULL
        self.root = NULL

    def __dealloc__(self):
        """
        Free any memory if it was actually allocated.
        """
        if self.prime:
            sage_free(self.prime)
        if self.root:
            sage_free(self.root)

    cdef long get_prime(self, Py_ssize_t i) except -1:
        """
        Return the i-th prime's residue characteristic.  Bounds are checked.
        """
        if i >= 0 and i < self.table_size:
            return self.prime[i]
        else:
            raise IndexError

    cdef long get_root(self, Py_ssize_t i) except 1099511627776: # bigger than any long that could occur
        """
        Return the root of `x^2-x-1` corresponding to the i-th prime. Bounds are checked.
        """
        if i >= 0 and i < self.table_size:
            return self.root[i]
        else:
            raise IndexError
        
    def __getitem__(self, Py_ssize_t i):
        """
        Get the i-th prime as a Prime object.  Here i must be between
        0 and self.table_size-1; Python negative indexes are *not*
        supported.

        EXAMPLES::

            sage: from psage.number_fields.sqrt5.prime import PrimesOfBoundedNorm
            sage: v = PrimesOfBoundedNorm(50); v
            Prime ideals of the ring of integers of Q(sqrt(5)) of norm less than 50
            sage: v[0]
            2a
            sage: v[3]
            11a
            sage: v[-1]
            Traceback (most recent call last):
            ...
            IndexError
            sage: v[100]
            Traceback (most recent call last):
            ...
            IndexError
        """
        cdef long p, r
        cdef Prime P
        if i >= 0 and i < self.table_size:
            P = PY_NEW(Prime)
            P.p = self.prime[i]; P.r = self.root[i]
            if P.r and self.prime[i-1] == self.prime[i]:
                P.first = False
            else:
                P.first = True
            return P
        raise IndexError

    def __repr__(self):
        """
        EXAMPLES::

            sage: from psage.number_fields.sqrt5.prime import PrimesOfBoundedNorm
            sage: PrimesOfBoundedNorm(50).__repr__()
            'Prime ideals of the ring of integers of Q(sqrt(5)) of norm less than 50'
        """
        return "Prime ideals of the ring of integers of Q(sqrt(5)) of norm less than %s"%self.bound

    def __len__(self):
        """
        Return the number of primes less than the bound.
        
        EXAMPLES::
        
            sage: from psage.number_fields.sqrt5.prime import PrimesOfBoundedNorm
            sage: len(PrimesOfBoundedNorm(50))
            14
            sage: len(PrimesOfBoundedNorm(10^7))
            664500
        """
        return self.table_size
        
    def __init__(self, long bound):
        """
        EXAMPLES::
        
            sage: from psage.number_fields.sqrt5.prime import PrimesOfBoundedNorm
            sage: P = PrimesOfBoundedNorm(50); type(P)
            <type 'psage.number_fields.sqrt5.prime.PrimesOfBoundedNorm'>
        """
        assert sizeof(long) == 8, "this object requires 64-bit"
        if bound >= 2**31:
            raise NotImplementedError, "bound must be less than 2^31"
        if bound < 0:
            raise ValueError, "bound must be nonnegative"
        self.bound = bound
        max_size = self._allocate_memory()
        self._enumerate_primes(max_size)
        self._sort_primes()
        self._reallocate_memory(max_size)

    cdef long _allocate_memory(self) except -1:
        # Allocate memory used to store the primes.  We very slightly
        # overestimate how much memory will be needed by assuming all
        # primes are split.  We give back the extra memory later.
        from sage.all import prime_pi
        max_size = 2*prime_pi(self.bound)
        self.prime = <long*> sage_malloc(sizeof(long)*max_size)
        if self.prime == NULL:
            raise MemoryError
        self.root = <long*> sage_malloc(sizeof(long)*max_size)
        if self.root == NULL:
            sage_free(self.prime)
            raise MemoryError
        return max_size

    cdef int _reallocate_memory(self, long max_size) except -1:
        # Give back the extra over-allocated memory.  This doesn't
        # free up much, since we only overestimated by a small amount.
        # But it is pretty fast (compared to the overall time), so we
        # do it.
        if self.table_size == max_size: return 0
        cdef long i, *prime, *root
        prime = <long*> sage_malloc(sizeof(long)*self.table_size)
        if prime == NULL:  # deal properly with out of memory condition, which may occur
            raise MemoryError
        root  = <long*> sage_malloc(sizeof(long)*self.table_size)
        if root == NULL:
            sage_free(prime)
            raise MemoryError
        memcpy(<void*>prime, <void*>self.prime, sizeof(long)*self.table_size)
        memcpy(<void*>root, <void*>self.root, sizeof(long)*self.table_size)
        sage_free(self.prime)
        sage_free(self.root)
        self.prime = prime
        self.root = root
        return 0

    def _enumerate_primes(self, long max_size):
        """
        This is an internal function that is called once when
        constructing this object.  It does the actual prime
        enumeration, but does not sort the primes.  It is called
        *after* the memory to store primes has been allocated (since
        otherwise there would be segfaults).

        EXAMPLES::

            sage: from psage.number_fields.sqrt5.prime import PrimesOfBoundedNorm
            sage: P = PrimesOfBoundedNorm(80); list(P)
            [2a, 5a, 3a, 11a, 11b, 19a, 19b, 29a, 29b, 31a, 31b, 41a, 41b, 7a, 59a, 59b, 61a, 61b, 71a, 71b, 79a, 79b]
            sage: P._enumerate_primes(len(P))

        Notice that the primes are no longer in the correct order,
        since we re-enumerated them, but didn't sort them::

            sage: list(P)
            [2a, 3a, 5a, 7a, 11a, 11b, 19a, 19b, 29a, 29b, 31a, 31b, 41a, 41b, 59a, 59b, 61a, 61b, 71a, 71b, 79a, 79b]        
        """
        cdef long p, i = 0, sr, r0, r1
        from sage.all import prime_range
        
        _sig_on
        for p in prime_range(self.bound):
            if i >= max_size:
                raise RuntimeError, "memory assumption violated"
            t = p % 5
            if t == 1 or t == 4:
                # split case
                self.prime[i] = p
                self.prime[i+1] = p
                # Compute a square root of 5 modulo p.
                sr = Fl_sqrt(5, p)
                # Find the two values of (1+sqrt(5))/2.
                r0 = Fl_div(1+sr, 2, p)
                r1 = Fl_div(1+p-sr, 2, p)
                # Next we normalize the roots to be between -p/2 and p/2.
                # NOTE: If you wanted to change the sort order to use the root between
                # 0 and p for some reason, delete the following four lines.
                if r0 >= p//2:
                    r0 -= p
                elif r1 >= p//2:
                    r1 -= p
                if r0 > r1: # swap
                    r0, r1 = r1, r0
                self.root[i] = r0
                self.root[i+1] = r1
                i += 2
            else:
                if p == 5:
                    self.prime[i] = p
                    self.root[i] = -2
                    i += 1
                elif p*p < self.bound:  # inert or ramified
                    self.prime[i] = p
                    self.root[i] = 0
                    i += 1
        _sig_off
        self.table_size = i

    def _sort_primes(self):
        """
        Used internally when constructing this object after the primes have been enumerated.
        This sorts them, assuming they are out of order precisely in the way that they
        would be out of order after running _enumerate_primes.

        EXAMPLES::

            sage: from psage.number_fields.sqrt5.prime import PrimesOfBoundedNorm
            sage: P = PrimesOfBoundedNorm(80); list(P)
            [2a, 5a, 3a, 11a, 11b, 19a, 19b, 29a, 29b, 31a, 31b, 41a, 41b, 7a, 59a, 59b, 61a, 61b, 71a, 71b, 79a, 79b]
            sage: P._enumerate_primes(len(P)); list(P)
            [2a, 3a, 5a, 7a, 11a, 11b, 19a, 19b, 29a, 29b, 31a, 31b, 41a, 41b, 59a, 59b, 61a, 61b, 71a, 71b, 79a, 79b]

        Now sorting puts them back in the right order::
        
            sage: P._sort_primes(); list(P)
            [2a, 5a, 3a, 11a, 11b, 19a, 19b, 29a, 29b, 31a, 31b, 41a, 41b, 7a, 59a, 59b, 61a, 61b, 71a, 71b, 79a, 79b]
        """
        # sort: move the inert prime to their proper position
        cdef long i, j, k, p, nrm
        if self.table_size >= 3:  # swap (3) and (sqrt(5))
            self.prime[1] = 5
            self.root[1] = -2
            self.prime[2] = 3
            self.root[2] = 0
        else:
            # already sorted
            return

        import math
        i = self.table_size - 2
        _sig_on
        while i >= 3:
            if self.root[i] == 0:
                j = i+1  # find spot to insert
                nrm = self.prime[i]*self.prime[i]
                while (j < self.table_size and
                       # the following scary thing is a conditional
                       # expression for norm of j-th prime:
                            nrm > (self.prime[j] if self.root[j] else self.prime[j]*self.prime[j])):
                    j += 1
                if j != i+1:
                    # now self.norms[j-1] <= self.norms[i] < self.norms[j]
                    # so we move what is in position i now to position j
                    # whilst moving everything from position i+1 to position j
                    # to the left by 1, using memmove.
                    p = self.prime[i]
                    memmove(<void*>(self.prime+i), <void*>(self.prime+i+1), sizeof(long)*(j-i-1))
                    memmove(<void*>(self.root+i), <void*>(self.root+i+1), sizeof(long)*(j-i-1))
                    self.prime[j-1] = p
                    self.root[j-1] = 0
            i -= 1
        _sig_off


def primes_of_bounded_norm(bound):
    """
    Return ordered list of all primes of the ring of integers of
    Q(sqrt(5)) of norm less than bound.  The primes are instances of a
    special fast Primes class.  They are sorted first by norm, then in
    the split case by the integer in (-p/2,p/2) congruent to
    (1+sqrt(5))/2.

    INPUT:
        - ``bound`` -- nonnegative integer

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
    
        sage: P.sage_prime()
        Fractional ideal (3*a - 2)

    You can also get the underlying residue characteristic::
    
        sage: P.p
        11

    And, the image of (1+sqrt(5))/2 modulo the prime (or 0 in the inert case)::
    
        sage: P.r
        -3
        sage: z = P.sage_prime(); z.residue_field()(z.number_field().gen())
        8

    For optimal speed you can use the Prime objects directly from
    Cython, which provides direct C-level access to the underlying
    data structure.
    """
    return list(PrimesOfBoundedNorm(bound))



