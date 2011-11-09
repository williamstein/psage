cdef extern from "rankbound.h":
    double rankbound(char * curve_string,
                 double logN,
                 long * bad_primes_,
                 int * ap_,
                 int len_bad_primes_,
                 double delta_,
                 int verbose_)

cdef extern from "math.h":
    double log(double)

cdef extern from "stdlib.h":
    void free(void* ptr)
    void* malloc(size_t size)
    void* realloc(void* ptr, size_t size)

from sage.schemes.elliptic_curves.ell_rational_field import EllipticCurve_rational_field

def xxx_rankbound(E, float delta, verbose = 0):
    if not isinstance(E, EllipticCurve_rational_field):
        raise NotImplementedError("Only elliptic curves over Q are implemented for now.")

    L = str(list(E.ainvs()))
    logN = log(E.conductor())

    bad_primes_list = [x.prime().gen() for x in E.local_data()]
    bad_primes_list = [x for x in bad_primes_list if x < 2**63]
    cdef long* bad_primes = <long*>malloc(sizeof(long) * len(bad_primes_list))
    cdef int* bad_ap = <int*>malloc(sizeof(int) * len(bad_primes_list))

    for n in range(len(bad_primes_list)):
        p = bad_primes_list[n]
        bad_primes[n] = p
        bad_ap[n] = E.ap(p)

    bound = rankbound(L, logN, bad_primes, bad_ap, len(bad_primes_list), delta, verbose)
    free(bad_primes)
    free(bad_ap)
    return bound
