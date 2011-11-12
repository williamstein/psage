include 'interrupt.pxi'

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

def xxx_rankbound(E, float delta, bad_primes = None, verbose = 0):
    r"""
    Assuming the Riemann hypothesis for `L(s, E)` and the BSD conjecture,
    this funtion computes and returns an upper bound for the rank of `E`,
    which must be defined over `Q`.

    More precisely, we compute the sum

    `\sum_\gamma f(\gamma, \Delta)`

    where `1/2 + i\gamma` runs over the nontrivial zeros of `L(s, E)`
    (counted with multiplicity) and

    `f(t, \delta) = (sin(\pi \Delta t)/(\pi \Delta t))^2`.

    If all of the `\gamma` are real, then this is an upper bound for the 
    analytic rank of `E`, and so if the BSD conjecture is true, it is an
    upper bound for the rank of `E`.


    INPUT:

     - ``E`` - An Elliptic Curve, whose rank is to be bounded.

     - ``delta`` - a real number between 1 and 5; the `\Delta` described above.

     - ``verbose`` - if > 0, print out some progress information.

    OUTPUT:

     -  The sum described above, which happens to be a bound for the rank if
        RH+BSD holds, expected to be accurate to about `5e-5`.


    EXAMPLES:

    This method works well on some high rank curves. Under BSD, the parity
    of the rank of the following curve is known, so the bound computed
    is tight.

    ::

        sage: from psage.ellcurve.xxx.rankbound import xxx_rankbound
        sage: E = EllipticCurve([1,0,0,-431092980766333677958362095891166,5156283555366643659035652799871176909391533088196]) # rank >= 20
        sage: xxx_rankbound(E, 2.1)
        21.48...

    It is also possible to specify a list of the bad primes, so that computing
    the conductor is fast. (In the following example, a larger value of Delta
    would give a better bound, but a tight bound seems unattainable at the
    moment.)

    ::

        sage: from psage.ellcurve.xxx.rankbound import xxx_rankbound
        sage: E = EllipticCurve([1, -1, 1, -20067762415575526585033208209338542750930230312178956502, 34481611795030556467032985690390720374855944359319180361266008296291939448732243429]) # rank >= 28
        sage: xxx_rankbound(E, 2.0, bad_primes = [2, 3, 5, 7, 11, 13, 17, 19, 48463, 20650099, 315574902691581877528345013999136728634663121, 376018840263193489397987439236873583997122096511452343225772113000611087671413])
        36.13...
        

    It also sometimes works well on some not so high rank curves.

    ::

        sage: from psage.ellcurve.xxx.rankbound import xxx_rankbound
        sage: E = EllipticCurve([0, -1, 1, 109792, 10201568]) # rank 0
        sage: xxx_rankbound(E, 2.0)
        0.53...

    We test the accuracy with some curves where we can actually
    compute the zeros. In both of the following cases, it is
    likely (but not certain) that the value computed by xxx_rankbound()
    is more accurate than the value computed by summing over zeros,
    since we are not finding enough zeros for high accuracy.

    ::

        sage: from psage.ellcurve.xxx.rankbound import xxx_rankbound
        sage: E = EllipticCurve('11a1') # rank 0
        sage: a = xxx_rankbound(E, 2.0); # answer is around .0027
        sage: abs(a - .0027) < 1e-5
        True
        sage: f = lambda t : (sin(RR(t * 2.0 * pi))/RR(t * pi * 2.0))^2 # long time
        sage: zeros = lcalc.zeros(1000, E)                              # long time
        sage: b = 2 * sum( (f(z) for z in zeros) )                      # long time
        sage: abs(a - b) < 1e-4                                         # long time
        True
        sage: E = EllipticCurve('389a1') # rank 2
        sage: a = xxx_rankbound(E, 1.1); # the answer is around 2.03
        sage: abs(a - 2.03) < 1e-4
        True
        sage: f = lambda t : (sin(RR(t * 1.1 * pi))/RR(t * pi * 1.1))^2 # long time
        sage: zeros = lcalc.zeros(4000, E)                              # long time
        sage: b = 2 * sum( (f(z) for z in zeros[2:]) ) + 2              # long time
        sage: abs(a - b) < 5e-4                                         # long time
        True    

    Often it doesn't work so well for curves of small rank and
    large conductor, but I'm not going to write down an example
    of that right now.

    NOTES:

    We use the explicit formula to compute the sum over zeros. (So
    we don't actually find the zeros.) As such, the complexity grows
    exponentially with delta (while, unfortunately, giving only slightly
    better bounds). Guidelines on a fast machine: delta = 2.0 is fast,
    delta = 2.5 is a little slower, delta = 3.0 takes a few minutes,
    delta = 3.5 takes at least a few hours, delta = 4.0 takes over a
    day (a few days?).

    """

    if not isinstance(E, EllipticCurve_rational_field):
        raise NotImplementedError("Only elliptic curves over Q are implemented for now.")

    L = str(list(E.ainvs()))

    if bad_primes is not None:
        N = 1
        for p in bad_primes:
            e = E.local_data(p).conductor_valuation()
            a = p**e
            N = N * a
        logN = float(log(N))
    else:
        logN = float(log(E.conductor()))
        bad_primes = [x.prime().gen() for x in E.local_data()]

    bad_primes = [x for x in bad_primes if x < 2**63]
    cdef long* bad_primes_array = <long*>malloc(sizeof(long) * len(bad_primes))
    cdef int* bad_ap = <int*>malloc(sizeof(int) * len(bad_primes))

    for n in range(len(bad_primes)):
        p = bad_primes[n]
        bad_primes_array[n] = p
        bad_ap[n] = E.ap(p)

    sig_on()
    #
    # We make no attempt to cleanup nicely if there is a keyboard interrupt.
    # Oh well, only a few bytes will be lost...
    #
    bound = rankbound(L, logN, bad_primes_array, bad_ap, len(bad_primes), delta, verbose)
    sig_off()
    free(bad_primes_array)
    free(bad_ap)
    return bound
