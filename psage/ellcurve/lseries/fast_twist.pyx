"""
Numerically computing L'(E/K,1) for E an elliptic curve over Q and K a
quadratic imaginary field satisfying the Heegner hypothesis.
"""

from cysignals.memory cimport sig_free,sig_malloc
from cysignals.signals cimport sig_on,sig_off

from sage.rings.all import ZZ

from psage.libs.smalljac.wrapper import SmallJac

from sage.finance.time_series cimport TimeSeries
from sage.stats.intlist cimport IntList

cdef extern:
    cdef double exp(double)
    cdef double log(double)
    cdef double sqrt(double)

cdef extern from "gsl/gsl_sf_expint.h":
    cdef double gsl_sf_expint_E1(double)

cdef double pi = 3.1415926535897932384626433833

cdef class anlist:
    cdef IntList a
    cdef object E
    cdef Py_ssize_t B
    
    def __init__(self, E, Py_ssize_t B):
        self.a = IntList(B)
        cdef int* a = self.a._values
        self.E = E
        self.B = B

        cdef long UNKNOWN = 10*B   # coefficients not yet known; safely outside Hasse interval
        cdef Py_ssize_t i
        for i in range(B):
            a[i] = UNKNOWN

        # easy cases
        a[0] = 0
        a[1] = 1

        cdef long N = E.conductor()
            
        ##############################################################
        # First compute the prime indexed coefficients using SmallJac
        # (and also Sage for bad primes)
        ##############################################################        
        a1,a2,a3,a4,a6 = E.a_invariants()
        if a1 or a2 or a3:
            _,_,_,a4,a6 = E.short_weierstrass_model().a_invariants()
        R = ZZ['x']
        f = R([a6,a4,0,1])
        C = SmallJac(f)
        sig_on()
        ap_dict = C.ap(0,B)
        sig_off()
        cdef Py_ssize_t p, q, n
        cdef object v
        for p, v in ap_dict.iteritems():
            if v is None:
                a[p] = E.ap(p)
            else:
                a[p] = v
            # prime powers
            # a_{p^r} := a_p * a_{p^{r-1}} - eps(p)*p*a_{p^{r-2}}
            q = p*p
            if N%p == 0:
                while q < B:
                    a[q] = a[p]*a[q//p]
                    q *= p
            else:
                while q < B:
                    a[q] = a[p]*a[q//p] - p*a[q/(p*p)]
                    q *= p
                
        ##############################################################
        # Fill in products of coprime numbers using multiplicativity,
        # since we now have all prime powers.
        ##############################################################
        primes = ap_dict.keys(); primes.sort()  # todo: speedup?
        for p in primes:
            q = p
            while q < B:
                # set a_{qn} = a_q * a_n for each n coprime to p with q*n < B with a_n known.
                n = 2
                while n*q < B:
                    if n%p and a[n] != UNKNOWN:
                        a[q*n] = a[q]*a[n]
                    n += 1
                q *= p
        

    def __repr__(self):
        return "List of coefficients a_n(E)."

    def __getitem__(self, Py_ssize_t i):
        return self.a[i]

cdef class FastHeegnerTwists:
    cdef Py_ssize_t B
    cdef TimeSeries a
    cdef object E
    #cdef object bad_primes
    cdef int root_number
    cdef int N
    
    def __init__(self, E, Py_ssize_t B):
        self.E = E
        self.N = E.conductor()
        #self.bad_primes = E.conductor().prime_divisors()
        self.B = B
        self.a = TimeSeries(B)
        cdef Py_ssize_t i
        cdef anlist v = anlist(E, B)
        for i in range(B):
            self.a[i] = v.a[i]
        self.root_number = E.root_number()

    def __repr__(self):
        return "Fast Heegner twist object associated to %s"%self.E

    def elliptic_curve(self):
        return self.E

    def discs(self, bound):
        return self.E.heegner_discriminants(bound)

    def twist_special_value(self, long D, tol=1e-6):
        # 1. Make floating point table of value of the quadratic
        # symbol (D/-) modulo 4*|D|.
        assert D < 0
        cdef long D4 = abs(D)*4
        chi = TimeSeries(D4)
        cdef int n
        Dz = ZZ(D)
        for n in range(len(chi)):
            chi._values[n] = Dz.kronecker(n)
        
        # 2. Approximate either L(E,chi,1) or L'(E,chi,1) depending on
        # root_number of E.
        # Conductor of E^D = N * |D|^2.
        cdef double  s = 0 # running sum below
        cdef double* a = self.a._values
        cdef double  X = 1/(abs(D) * sqrt(self.N))

        cdef Py_ssize_t nterms

        nterms = <int>(log(tol*(1-exp(-2*pi*X))/2)/(-2*pi*X)) + 3
        if nterms >= self.B:
            raise ValueError, "not enough terms of L-series known (%s needed, but %s known)"%(
                nterms, self.B)

        if self.root_number == -1:
            # compute L(E,chi,1) = 
            #        2 * sum_{n=1}^{k} (chi(n) * a_n / n) * exp(-2*pi*n/sqrt(N*D^2))
            for n in range(1, nterms):
                s += (chi._values[n%D4] * a[n] / n)  * exp(-2*pi*n*X)
        else:
            # compute L'(E,chi,1) = 
            #      2 * sum_{n=1}^{k} (chi(n) * a_n / n) * E_1(2*pi*n/sqrt(N*D^2))
            for n in range(1, nterms):
                s += (chi._values[n%D4] * a[n] / n)  * gsl_sf_expint_E1 (2*pi*n*X)
                
        return 2*s, nterms
    
        
        
