
from sage.stats.intlist cimport IntList

from sage.all import prime_range

def prime_powers_intlist(Py_ssize_t B):
    """
    Return IntList of the prime powers and corresponding primes, up to
    B.  The list is *not* in the usual order; instead the order is like
    this: 2,4,8,...,3,9,27,..., 5,25,..., etc.

    INPUT:

        - B -- positive integer

    EXAMPLES::
    
        sage: from psage.ellcurve.lseries.helper import prime_powers_intlist
        sage: prime_powers_intlist(10)
        ([1, 2, 4, 8, 3, 9, 5, 7], [1, 2, 2, 2, 3, 3, 5, 7])
        sage: prime_powers_intlist(10)
        ([1, 2, 4, 8, 3, 9, 5, 7], [1, 2, 2, 2, 3, 3, 5, 7])
        sage: prime_powers_intlist(20)
        ([1, 2, 4, 8, 16 ... 7, 11, 13, 17, 19], [1, 2, 2, 2, 2 ... 7, 11, 13, 17, 19])
        sage: list(sorted(prime_powers_intlist(10^6)[0])) == prime_powers(10^6)
        True
        sage: set(prime_powers_intlist(10^6)[1][1:]) == set(primes(10^6))
        True
    """
    v = prime_range(B, py_ints=True)
    cdef IntList w = IntList(len(v)*2), w0 = IntList(len(v)*2)
    w[0] = 1
    w0[0] = 1
    # Now fill in prime powers
    cdef Py_ssize_t i=1
    cdef int p
    cdef long long q   # so that the "q < B" test doesn't fail due to overflow.
    for p in v:
        q = p
        while q < B:
            w._values[i] = q
            w0._values[i] = p
            q *= p
            i += 1
    return w[:i], w0[:i]

def extend_multiplicatively(IntList a):
    """
    Given an IntList a such that the a[p^r] is filled in, for all
    prime powers p^r, fill in all the other a[n] multiplicatively.

    INPUT:
        - a -- IntList with prime-power-index entries set; all other
          entries are ignored

    OUTPUT:
        - the input object a is modified to have all entries set
          via multiplicativity.

    EXAMPLES::

        sage: from psage.ellcurve.lseries.helper import extend_multiplicatively
        sage: B = 10^5; E = EllipticCurve('389a'); an = stats.IntList(B)
        sage: for pp in prime_powers(B):
        ...     an[pp] = E.an(pp)
        ...
        sage: extend_multiplicatively(an)
        sage: list(an) == E.anlist(len(an))[:len(an)]
        True
    """
    cdef Py_ssize_t i, B = len(a)
    cdef IntList P, P0
    P, P0 = prime_powers_intlist(B)

    # Known[n] = 1 if a[n] is set.  We use this separate array
    # to avoid using a sentinel value in a.
    cdef IntList known = IntList(B)
    for i in range(len(P)):
        known._values[P[i]] = 1
        
    cdef int k, pp, p, n
    # fill in the multiples of pp = prime power
    for i in range(len(P)):
        pp = P._values[i]; p = P0._values[i]
        k = 2
        n = k*pp
        while n < B:
            # only consider n exactly divisible by pp
            if k % p and known._values[k]:
                a._values[n] = a._values[k] * a._values[pp]
                known._values[n] = 1
            n += pp
            k += 1
        
        
    
    
    
