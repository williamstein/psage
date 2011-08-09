
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
        
        
def extend_multiplicatively_generic(list a):
    """
    Given a list a of numbers such that the a[p^r] is filled in, for
    all prime powers p^r, fill in all the other a[n] multiplicatively.

    INPUT:
        - a -- list with prime-power-index entries set; all other
          entries are ignored

    OUTPUT:
        - the input object a is modified to have all entries set
          via multiplicativity.

    EXAMPLES::

        sage: from psage.ellcurve.lseries.helper import extend_multiplicatively_generic
        sage: B = 10^5; E = EllipticCurve('389a'); an = [0]*B
        sage: for pp in prime_powers(B):
        ...     an[pp] = E.an(pp)
        ...
        sage: extend_multiplicatively_generic(an)
        sage: list(an) == E.anlist(len(an))[:len(an)]
        True


    A test using large integers::
    
        sage: v = [0, 1, 2**100, 3**100, 4, 5, 0]
        sage: extend_multiplicatively_generic(v)
        sage: v
        [0, 1, 1267650600228229401496703205376, 515377520732011331036461129765621272702107522001, 4, 5, 653318623500070906096690267158057820537143710472954871543071966369497141477376]
        sage: v[-1] == v[2]*v[3]
        True

    A test using variables::

        sage: v = list(var(' '.join('x%s'%i for i in [0..30]))); v
        [x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30]
        sage: extend_multiplicatively_generic(v)
        sage: v
        [x0, x1, x2, x3, x4, x5, x2*x3, x7, x8, x9, x2*x5, x11, x3*x4, x13, x2*x7, x3*x5, x16, x17, x2*x9, x19, x4*x5, x3*x7, x11*x2, x23, x3*x8, x25, x13*x2, x27, x4*x7, x29, x2*x3*x5]
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
                a[n] = a[k] * a[pp]
                known._values[n] = 1
            n += pp
            k += 1
    
    
    
