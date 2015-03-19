#
# A new implementation of Dirichlet characters based on the numbering scheme
# devised by Brian Conrey.
#
include "interrupt.pxi"  # ctrl-c interrupt block support
include "stdsage.pxi"  # ctrl-c interrupt block support
include "cdefs.pxi"

from libc.stdlib cimport malloc, free

from sage.all import factor,        \
                     primitive_root,\
                     euler_phi,     \
                     gcd,           \
                     lcm,           \
                     is_prime,      \
                     DirichletGroup,\
                     vector,        \
                     Integer,       \
                     power_mod,     \
                     prod,          \
                     crt,           \
                     mod,           \
                     inverse_mod,   \
                     multiplicative_order, \
                     pi,            \
                     RR,            \
                     CC,            \
                     ZZ,            \
                     diagonal_matrix, \
                     Mod,           \
                     round,         \
                     imag

from sage.modular.dirichlet import DirichletCharacter

import cmath

cdef complex twopii = 3.1415926535897932384626433833 * 2.0 * 1.0j
cdef float PI = float(pi)

cdef class DirichletGroup_conrey:
    #
    # Note: perhaps the discrete log tables should be stored
    # separately for each prime. This will make computation a
    # little slower, but it would be possible to work efficiently
    # when the modulus is huge but divisible only by small primes
    # to small powers.
    #
    # Random question for self: Given a discrete log table mod p,
    # is it easy to solve the discrete log problem mod p^a? If so,
    # it would be possible to remove the last three words from the
    # above paragraph.
    #
    # Pascal >   The answer to the above question is YES. Just
    #        > do a p-adic decomposition (keyword Pohlig-Hellman)
    #        > of the result.
    #        >   In general, there sould be a mix of precomputed
    #        > logs for small primes (to any power) and discrete
    #        > logs computations for large primes.
    #        > This is a TODO.

    cdef long q             # the modulus
                            
    cdef long q_even        # for computation we will strip out the even
    cdef long q_odd         # factors from the modulus. q == q_even * q_odd
                            
    cdef long k             # the number of factors of q_odd
                            
    cdef long * primes      # a list of the odd prime factors of the modulus
    cdef long * exponents   # a list of the exponents of the odd prime factors in the factorization
    cdef long * generators  # a primitive root for each odd prime factor

    cdef long precomp       # shall we precompute logs ?

    cdef long * A           # exponent vectors:
                            # for each m coprime to q_odd we store an array
                            # with the property that
                            #
                            #   m == g[j]**A[m][j] mod p[j]**e[j]
                            #
                            # (where "A[m][k] == A[m * self.k + j]")
                            #
                            # where g[j] is a primitive root mod p[j]**e[j],
                            # and p[j] is the j-th prime factor of q_odd.
                            #
                            # This array is the obstacle that will prevent this
                            # implementation from working reasonably for very
                            # large modulus. We will need something else which
                            # does not use any precomputation for that case.
 
    cdef long * B           # exponents for q_even:
                            # for each odd m, 0 <= m < q_even, we will compute B
                            # so that
                            # 
                            #   m == B[m-1] * 5**B[m] mod q_even,
                            # 
                            # where B[m-1] = +- 1 and 0 <= B[m] < q_even/4

    cdef long * PHI         # PHI[j] = phi(q_odd)/phi(p[j]**e[j]). This will make it easier
                            # to compute the characters.

    cdef long phi_q_odd     # phi(q_odd)
    cdef long phi_q         # phi(q)
    
    cdef complex * zeta_powers_odd  # an array holding powers of a root of unity.
                                    # this should be the only part of the code that
                                    # needs to change in order to work over a cyclotomic field

    cdef complex * zeta_powers_even # for the even part of the character
    
    cdef _standard_dirichlet_group
    cdef _invariants
    cdef _gens

    def __cinit__(self, modulus, basering = None):
        if modulus <= 0:
            raise ArithmeticError("The modulus of a Dirichlet group must be a positive integer.")

        try:
            self.q = modulus
        except OverflowError:
            raise NotImplementedError("Currently this implementation does not allow a modulus that large.")
            
        #self.precomp = (self.q < 100000) # should test to find right value

        self.precomp = 1

        self.q_even = 1
        self.q_odd = self.q
        while self.q_odd % 2 == 0:
            self.q_odd = self.q_odd/2
            self.q_even = self.q_even * 2

        if self.q_odd > 1:
            X = factor(self.q_odd)
            self.k = len(X)
            self.primes = <long *>malloc(self.k * sizeof(long))
            self.exponents = <long *>malloc(self.k * sizeof(long))
            for n in range(self.k):
                self.primes[n] = X[n][0]
                self.exponents[n] = X[n][1]
            self.generators = <long *>malloc(self.k * sizeof(long))
            self.PHI = <long *>malloc(self.k * sizeof(long))
            if self.precomp:
                self.A = <long*>malloc(self.q_odd * self.k * sizeof(long))
                self.zeta_powers_odd = <complex*>malloc(self.q * sizeof(complex))

        if self.precomp and self.q_even > 4:
            # We are only going to use zeta_powers_even if q_even is large enough.
            # When q_even == 2, it would just be {1}, and when q_even == 4, it
            # would just be {1,-1}.
            #
            # This way, it is always the case that, if zeta_powers_even has been
            # initialized, it will be of size q_even/4

            self.B = <long*>malloc(self.q_even * sizeof(long))
            self.zeta_powers_even = <complex*>malloc(self.q_even/4 * sizeof(complex))

    def __init__(self, modulus, basering = None):
        # Once we've hit this stage, all of our arrays are allocated,
        # and both self.prime and self.exponents contain the right things.
        #
        # We now set up the rest of the precomputed arrays.

        self.phi_q_odd = euler_phi(self.q_odd)
     
        if self.q_even > 1:
            self.phi_q = self.phi_q_odd * self.q_even/2
        else:
            self.phi_q = self.phi_q_odd

        cdef long g
        cdef long a
        for j in range(self.k):
            x = self.primes[j]**self.exponents[j]
            g = primitive_root(x)
            self.generators[j] = g
            phi = self.primes[j]**(self.exponents[j] - 1) * (self.primes[j] - 1)
            self.PHI[j] = self.phi_q_odd/phi
            if self.precomp:
                a = 1
                for l in range(phi):
                    for m in range(a, self.q_odd, x):
                        self.A[m * self.k + j] = l
                    a = (a * g) % x
        #
        # Store a flag in A for each m that is not coprime to q_odd.
        # (This will save on expensive gcd computations later.)
        #
        if self.precomp:
            if self.q_odd > 1:
                for m in range(self.q_odd):
                    if gcd(m,self.q_odd) > 1:
                        self.A[m * self.k] = -1

        #
        # Compute a table of powers of the root of unity. This will
        # save on expensive calls to exp() later. It does increase
        # memory usage by an appreciable amount, though, so we might
        # want to add an option to not do this.
        #
        # We will _need_ to not do this later when allowing very
        # large moduli.
        #
        if self.precomp and self.q_odd > 1:
            for n in range(self.phi_q_odd):
                self.zeta_powers_odd[n] = cmath.exp(twopii * n/<double>self.phi_q_odd)

        cdef long pow_five = 1
        if self.precomp and self.q_even > 4:
            for n in range(self.q_even/4):
                self.zeta_powers_even[n] = cmath.exp(twopii * n * 4/<double>self.q_even)

            for e in range(self.q_even/4):
                self.B[pow_five] = e
                self.B[pow_five - 1] = 1
                self.B[self.q_even - pow_five] = e
                self.B[self.q_even - pow_five - 1] = -1
                pow_five = pow_five * 5
                pow_five = pow_five % self.q_even

    cpdef long _chi_odd_exponent(self, long m, long n):
        r"""
        BE CAREFUL CALLING THIS. It implicitly assumes:
            - 1 <= m < self.q_odd
            - 1 <= n < self.q_odd
            - gcd(m, self.q_odd) == 1
            - gcd(n, self.q_odd) == 1
            - (anything else?)
        """
        cdef long x = 0
        if self.precomp:
            for j in range(self.k):
                x += self.A[m * self.k + j]*self.A[n * self.k + j]*self.PHI[j]
                x = x % self.phi_q_odd
            return x;
        else:
            for j in range(self.k):
                pj =  self.primes[j]**self.exponents[j]
                gj = Mod(self.generators[j],pj)
                logm = Mod(m,pj).log(gj)
                logn = Mod(n,pj).log(gj)
                x = x + self.PHI[j]*logm*logn % self.phi_q_odd
            return x
            

    cpdef long _chi_even_exponent(self, long m, long n):
        r"""
        BE CAREFUL CALLING THIS. It implicitly assumes that:
            - 0 < m < self.q_even
            - 0 < n < self.q_even
            - self.q_even > 4
            - m and n are odd
        """
        cdef long exponent = 0
        if self.precomp:
            exponent = self.B[m]*self.B[n]
            if self.B[m-1] == -1 and self.B[n-1] == -1:
                exponent += self.q_even/8
            return exponent % (self.q_even/4)
        else:
            if self.q_even > 2:
                if m % 4 == 3 and n % 4 == 3:
                    exponent = self.q_even//8
                if self.q_even > 4:
                    g2 = Mod(5,self.q_even)
                    if(m % 4 == 3):
                        m = -m
                    if(n % 4 == 3):
                        n = -n
                    logm = Mod(m,self.q_even).log(g2)
                    logn = Mod(n,self.q_even).log(g2)
                    exponent += logn*logn*self.q_even//4
                return exponent % (self.q_even/4)
            else:
                return 0
            
    cpdef complex chi(self, long m, long n):
        if not self.precomp:
            raise NotImplementedError
        cdef complex odd_part = 1
        cdef complex even_part = 1
        if self.q_even > 1:
            if m % 2 == 0 or n % 2 == 0:
                return 0
            elif self.q_even == 2:
                even_part = 1
            elif self.q_even == 4:
                if m % 4 == 3 and n % 4 == 3:
                    even_part = -1
                else:
                    even_part = 1
            else:
                even_part = self.zeta_powers_even[self._chi_even_exponent(m % self.q_even, n % self.q_even)]
        if self.q_odd > 1:
            m = m % self.q_odd
            n = n % self.q_odd
            if self.A[m * self.k] == -1 or self.A[n * self.k] == -1:
                odd_part = 0;
            else:
                odd_part = self.zeta_powers_odd[self._chi_odd_exponent(m,n)]

        return even_part * odd_part

    def __iter__(self):
        cdef long n = 1
        if self.precomp:
            while n <= self.q:
                if self.q_odd == 1 or self.A[(n % self.q_odd) * self.k] != -1:
                    if self.q_even == 1 or n % 2 == 1:
                        yield self._getitem_(n)
                n = n + 1
        else:
            while n <= self.q:
                if gcd(n,self.q) == 1:
                    yield self._getitem_(n)
                n = n + 1

    def primitive_characters(self):
        for chi in self:
            if chi.is_primitive():
                yield chi

    def __dealloc__(self):
        if self.primes != NULL:
            free(self.primes)
        if self.exponents != NULL:
            free(self.exponents)
        if self.generators != NULL:
            free(self.generators)
        if self.PHI != NULL:
            free(self.PHI)
        if self.A != NULL:
            free(self.A)
        if self.zeta_powers_odd != NULL:
            free(self.zeta_powers_odd)
        if self.zeta_powers_even != NULL:
            free(self.zeta_powers_even)
        if self.B != NULL:
            free(self.B)

    def __getitem__(self, n):
        return self._getitem_(n)

    cdef DirichletCharacter_conrey _getitem_(self, long n):
        return DirichletCharacter_conrey(self, n)

    def __repr__(self):
        return "Group of dirichlet characters with modulus %d" % self.q

    def standard_dirichlet_group(self):
        """
        Return the "standard" Sage Dirichlet group with the same modulus,
        when characters taking values in a cyclotomic field. This is only
        computed when it is asked for, but it is cached after being
        computed.

        Maybe this function needs a better name.
        """

        if self._standard_dirichlet_group is None:
            self._standard_dirichlet_group = DirichletGroup(self.q)

        return self._standard_dirichlet_group


    cpdef modulus(self):
        return self.q

    cpdef order(self):
        return self.phi_q

    cpdef _gen_invariants(self):
        """
        compute elementary divisors structure of the group
        first make a list w of cyclic components over primes
        and lift local generators
          gj == g[j] mod pj
          gj == 1    mod qj, qj = q//pj
          uj qj + vj pj = 1
          -> gj = uj qj g[j] + vj pj = 1 + uj qj (g[j]-1)
        """
        w,g = [], []
        for j in range(self.k):
            p,e = self.primes[j], self.exponents[j]
            pj = p**(e - 1)
            w.append(pj*(p-1)) # phi(p**e)
            pj *= p
            qj = self.q//pj
            assert qj*pj == self.q
            gj = self.generators[j]
            uj = inverse_mod(qj,pj)
            gj = 1 + qj*uj*(gj-1) % self.q
            g.append(gj)
        if self.q_even >= 4:
            w.append(2)
            p2, q2 = self.q_even, self.q_odd
            g2 = -1
            u2 = inverse_mod(q2, p2)
            g2 = 1 + q2*u2*(g2-1) % self.q
            g.append(g2)
        if self.q_even > 4:
            w.append(self.q_even//4)
            g2 = 5
            g2 = 1 + q2*u2*(g2-1) % self.q
            g.append(g2)
        """
        then compute the Smith normal form
        and perform base change
        """
        k = len(w)
        W = diagonal_matrix(ZZ,w)
        d,u,v = W.smith_form()
        gens, inv = [], []
        u = u.inverse()
        for j in range(k-1,-1,-1):
            if d[j,j]>1:
                inv.append(d[j,j])
                s = 1
                for l in range(k):
                    s = s * g[l]**u[l,j] % self.q
                gens.append(s)
        self._gens = gens
        self._invariants = inv

    def gens(self):
        if self._gens is None:
            self._gen_invariants()
        return tuple(self._gens)

    def invariants(self):
        if self._invariants is None:
            self._gen_invariants()
        return tuple(self._invariants)

    cpdef zeta_order(self):
        r"""
        Every character in this Dirichlet group takes values in some
        cyclotomic field `\Q(\zeta_n)`. This function returns the smallest
        such `n`.

        Warning: We don't actually use this function internally for our table
        of roots of unity right now, and so, internally, we may use a larger
        zeta order.

        EXAMPLES::

            sage: from dirichlet_conrey import *
            sage: G = DirichletGroup_conrey(37)
            sage: G.zeta_order()
            36

        TESTS::
            
            sage: from dirichlet_conrey import *
            sage: L = [ZZ.random_element(1, 100) for n in range(10)]
            sage: for n in L:
            ...    if DirichletGroup_conrey(n).zeta_order() != DirichletGroup(n).zeta_order():
            ...        print "wrong answer for n =", n
        """

        if self.q > 2:
            return self.invariants()[0]
        else:
            return 1

    cpdef from_sage_character(self, chi):
        r"""
        Given a 'sage character' (that is, a Dirichlet character in the
        standard sage format), return the character that correponds to it.

        EXAMPLES::
            sage: from dirichlet_conrey import *
            sage: q = ZZ.random_element(2, 500)
            sage: q = q - 1 + q%2
            sage: G = DirichletGroup_conrey(q)
            sage: G[1] == G.from_sage_character(G[1].sage_character())
            True
            sage: G[2] == G.from_sage_character(G[2].sage_character())
            True
            sage: x = ZZ.random_element(1,q)
            sage: while gcd(x,q) > 1:
            ...    x = ZZ.random_element(1,q)
            sage: if G[x] == G.from_sage_character(G[x].sage_character()):
            ...    print True
            ... else:
            ...    print q, x
            True
            sage: q = 2 * q
            sage: G = DirichletGroup_conrey(q)
            sage: G[1] == G.from_sage_character(G[1].sage_character())
            True
            sage: x = ZZ.random_element(1,q)
            sage: while gcd(x,q) > 1:
            ...    x = ZZ.random_element(1,q)
            sage: G[x] == G.from_sage_character(G[x].sage_character())
            True
            sage: q = 2 * q
            sage: G = DirichletGroup_conrey(q)
            sage: G[1] == G.from_sage_character(G[1].sage_character())
            True
            sage: x = ZZ.random_element(1,q)
            sage: while gcd(x,q) > 1:
            ...    x = ZZ.random_element(1,q)
            sage: G[x] == G.from_sage_character(G[x].sage_character())
            True
            sage: q = 2 * q
            sage: G = DirichletGroup_conrey(q)
            sage: G[1] == G.from_sage_character(G[1].sage_character())
            True
            sage: x = ZZ.random_element(1,q)
            sage: while gcd(x,q) > 1:
            ...    x = ZZ.random_element(1,q)
            sage: G[x] == G.from_sage_character(G[x].sage_character())
            True
            sage: q = 2 * q
            sage: G = DirichletGroup_conrey(q)
            sage: G[1] == G.from_sage_character(G[1].sage_character())
            True
            sage: x = ZZ.random_element(1,q)
            sage: while gcd(x,q) > 1:
            ...    x = ZZ.random_element(1,q)
            sage: if G[x] == G.from_sage_character(G[x].sage_character()):
            ...    print True
            ... else:
            ...    print q, x
            True
        """
        #
        # At odd prime powers, it is relatively easy to construct a character
        # in our format from the sage character. For even prime powers it is
        # a little bit trickier, but not much worse.
        #
        # So to construct our character, we will decompose the sage character
        # into a product of characters to prime power modulus, do the translation
        # there, and then glue everything back together.
        #
        if chi.modulus() != self.q:
            raise ArithmeticError("The character we are translating from must have the same modulus as this Dirichlet group.")

        decomposition = chi.decomposition()
        n_even = 1

        # We deal with the even part first. First of all, we strip off the
        # even modulus from the decomposition of the character.

        if self.q_even > 1:
            psi, decomposition = decomposition[0], decomposition[1:]

        # Now the even part splits into 3 cases depending on how
        # many times 2 divides q. (i.e. 1 time, 2 times, or more than 2
        # times.)

        # if 2 divides q just once, the relevant character mod 2 is
        # always the trivial one
        if self.q_even == 2:
            n_even = 1

        # if 2 divides q twice, then we just need to check if the
        # part at two is trivial or not.
        elif self.q_even == 4:
            if psi.is_trivial():
                n_even = 1
            else:
                n_even = 3

        # When 8 divides q, we want to somehow find the exponents
        # of the value of the character on 5 and -5, and use these
        # to construct the character that we are interested in.
        #
        # This is going to be ugly...
        elif self.q_even > 4:
            psi5_exponent = round(imag(CC(psi(5)).log() * (self.q_even/4)/(2*pi)))
            n_even = power_mod(5, psi5_exponent, ZZ(self.q_even))
            if psi(5) != psi(-5):
                n_even = self.q_even - n_even

        # I think this might work...
        exponent_vector = [ZZ(psi.element()[0]) for psi in decomposition]
        odd_prime_powers = [ZZ(self.primes[j])**ZZ(self.exponents[j]) for j in range(self.k)]
        generators = [ZZ(self.generators[j]) for j in range(self.k)]
        n_odd = crt(
            [power_mod(g, a, pp) for (g,a,pp) in zip(generators, exponent_vector, odd_prime_powers)],
             odd_prime_powers)

        if n_odd == 0:
            n_odd = 1
        n = crt(n_odd, n_even, ZZ(self.q_odd), ZZ(self.q_even))

        return self[n]

def test_conversion(q):
    G = DirichletGroup_conrey(q)
    for chi in G:
        if not G.from_sage_character(chi.sage_character()) == chi:
            print "Failure for q = {}, chi number {}".format(q, chi.number())
            break

cdef class DirichletCharacter_conrey:
    cdef long _n        # we will store the number used to create this character,
    cdef _number         # e.g., -1, but for all computations we use _n, which is number % q.
    cdef long q         # the modulus of this character

    cdef DirichletGroup_conrey _parent

    def __init__(self, DirichletGroup_conrey parent, long n):
        """
            The nth character for the Dirichlet Group parent.
        """
        self._parent = parent
        self._number = n
        self._n = n % parent.q
        self.q = parent.q

    def __call__(self, long m):
        return self.value(m)

    def __cmp__(self, DirichletCharacter_conrey other):
        r"""
        Compare self to other. Return equality if and only if the moduli and
        the character number are the same. When different, characters are first
        ordered by modulus and then by number.

        EXAMPLES:

            sage: from dirichlet_conrey import *
            sage: G = DirichletGroup_conrey(17)
            sage: G2 = DirichletGroup_conrey(20)
            sage: chi1 = G[3]
            sage: chi2 = G[20]
            sage: chi3 = G[21]
            sage: chi1 == chi2
            True
            sage: chi1 < chi3
            True
            sage: chi2 > chi3
            False
            sage: chi1 < G2[0]
            True
        """
        if(self._parent.q != other._parent.q):
            return cmp(self._parent.q, other._parent.q)
        else:
            return cmp(self._n, other._n)

    def conductor(self):
        r"""
        Return the conductor of this character.

        TESTS::

        sage: from dirichlet_conrey import *
        sage: G = DirichletGroup_conrey(17)
        sage: G[1].conductor()
        1
        sage: G[18].conductor()
        1
        sage: G[2].conductor()
        17

        This is tested more thoroughly in primitive_character().
        """
        odd_parts = [self._primitive_part_at_known_p(k) for k in range(self._parent.k)]
        odd_conductor = prod( [c for (n, c) in odd_parts] )
        _, even_conductor = self._primitive_part_at_two()
        return odd_conductor * even_conductor

    cpdef decomposition(self):
        r"""
        Return the factorization of this character into characters of prime
        power modulus.

        EXAMPLES:

            sage: from dirichlet_conrey import *
            sage: G = DirichletGroup_conrey(8 * 3 * 25)
            sage: (chi1, chi2, chi3) = G[7].decomposition()
            sage: chi1.modulus()
            8
            sage: chi2.modulus()
            3
            sage: chi3.modulus()
            25
            sage: chi1 * chi2 * chi3 == G[7]
            True

        TESTS:
            sage: chi1.sage_character().extend(600).maximize_base_ring() * chi2.sage_character().extend(600).maximize_base_ring() * chi3.sage_character().extend(600).maximize_base_ring() == G[7].sage_character()
            True
        """

        cdef long n = self._n
        cdef long q
        L = []
        if self._parent.q_even > 1:
            L.append( DirichletGroup_conrey(self._parent.q_even)[n % self._parent.q_even])

        for m in range(self._parent.k):
            q = self._parent.primes[m]**self._parent.exponents[m]
            L.append( DirichletGroup_conrey(q)._getitem_(n % q) )

        return L
    
    cpdef long exponent(self, long m):
        r"""
        Return the number a such that chi(m) = e(a/phi(q)).
        """
        cdef long exponent
        cdef long q_even = self._parent.q_even
        cdef long q_odd = self._parent.q_odd

        if q_odd > 1:
            odd_exponent = self._parent._chi_odd_exponent(self._n % q_odd, m % q_odd)
        else:
            odd_exponent = 0

        if q_even > 4:
            even_exponent = self._parent._chi_even_exponent(self._n % q_even, m % q_even)
            even_exponent *= 2  # the function just above computes the exponent of
                                # e(1/ (q_even/4) ), but we want the exponent of
                                # e(1/phi(q_even)) = e(1/(q_even/2))
        elif q_even == 4:
            if (self._n % q_even) == 3 and (m % q_even) == 3:
                even_exponent = 1
            else:
                even_exponent = 0
        else:
            even_exponent = 0

        if q_even == 1: # special case because phi(1) != 1/2.
            exponent = odd_exponent
        else:
            exponent = odd_exponent * q_even/2 + even_exponent * self._parent.phi_q_odd
    
        # we now have the value of chi(m) as e(exponent/phi(q))

        # it could be equal to phi(q), though, and in that case we
        # want it to be zero...
        if exponent == self._parent.phi_q:
            exponent -= self._parent.phi_q

        return exponent

    cpdef extend(self, M):
        r"""
        Return the extension of this character to a character of modulus M,
        where M is a multiple of the modulus of this character.

        EXAMPLES:

            sage: from dirichlet_conrey import *
            sage: G = DirichletGroup_conrey(17)
            sage: chi = G[5].extend(17 * 3); chi
            Dirichlet character with index 22 modulo 51
            sage: chi.sage_character() == G[5].sage_character().extend(17 * 3).maximize_base_ring()
            True
            sage: chi2 = DirichletGroup_conrey(17 * 3)[7]
            sage: chi2.extend(17 * 9).sage_character() == chi2.sage_character().extend(17 * 9).maximize_base_ring()
            True
        """

        if M % self.modulus() != 0:
            raise ArithmeticError("M must be a multiple of the modulus")

        # This is a little tricky with the definition of characters that we
        # are using, and it is easiest to do this a prime at a time. For an
        # odd prime p, to extend a character mod p^a to a character mod p^e
        # we just raise the index to the p^(e-a)th power. We treat independent
        # primes separately and glue things together using the chinese
        # remainder theorem.

        cdef DirichletGroup_conrey G = DirichletGroup_conrey(M)
        cdef long x, y
        y = 0
        indices = []
        moduli = []
        for x in range(self._parent.k):
            while G.primes[y] != self._parent.primes[x]:
                indices.append(1)
                moduli.append(G.primes[y]**G.exponents[y])
                y = y + 1
            if G.exponents[y] == self._parent.exponents[x]:
                q = G.primes[y]**G.exponents[y]
                indices.append(self._n % q)
                moduli.append(q)
            else:
                p = self._parent.primes[y]
                q1 = p**self._parent.exponents[y]
                q2 = p**G.exponents[y]
                indices.append(power_mod(self._n % q1, q2/q1, q2))
                moduli.append(q2)

        # that takes care of the even part of the modulus. still need to deal
        # with the odd part.

        cdef q_even1 = self._parent.q_even
        cdef q_even2 = G.q_even

        if q_even1 == q_even2:
            indices.append(self._n % q_even1)
        else:
            raise NotImplementedError

            # this stuff isn't really written yet...
            if q_even1 <= 2:
                indices.append(1)
            if q_even1 == 4:
                if self._n % 4 == 1:
                    indices.append(1)
                else:
                    indices.append(q_even2/2 - 1)
            if q_even1 == 8:
                pass
                
        moduli.append(q_even2)
        return G[crt(indices, moduli)]

    def galois_orbit(self):
        r"""
        Return the galois conjugates of this character.

        This is going to be rather inefficient for now.

        TESTS::

            sage: from dirichlet_conrey import *
            sage: G = DirichletGroup_conrey(56)
            sage: chi = G[3]
            sage: set([psi.sage_character() for psi in chi.galois_orbit()]) == set(chi.sage_character().galois_orbit())
            True
        """

        L = []
        N = self._parent.zeta_order()
        q = self._parent.modulus()
        for a in range(N):
            if gcd(a,N) == 1:
                L.append(power_mod(self._n, a, q))

        return [self._parent[n] for n in set(L)]




    cpdef complex gauss_sum(self, long a = 1):
        r"""
        Return the Gauss sum

        ``\tau_a(\chi) = \sum_{n=0}^{q-1} \chi(n) exp(2 \pi i an/q)``

        We compute in the complex numbers, so we compute and return this as
        a complex double precision number.

        EXAMPLES:
        sage: from dirichlet_conrey import *
        sage: G = DirichletGroup_conrey(17)
        sage: G[3].gauss_sum()
        (2.325224300372...+3.404898229456...j)
        
        When `\chi` is primitive, and `a` is relatively prime to `q`,
        the Gauss sum always has asbolute value
        `\sqrt(q)`.

        sage: G = DirichletGroup_conrey(17)
        sage: G[4].is_primitive()
        True
        sage: abs(G[4].gauss_sum(5))
        4.12310562561...
        sage: sqrt(17.0)
        4.12310562561766

        TESTS::
        sage: from dirichlet_conrey import *
        sage: G = DirichletGroup_conrey(4 * 11)
        sage: [abs(chi.gauss_sum(3) - chi.sage_character().gauss_sum_numerical(a=3)) < 5e-14 for chi in G] == [True] * euler_phi(4 * 11)
        True
        sage: G = DirichletGroup_conrey(23)
        sage: [abs(chi.gauss_sum() - chi.sage_character().gauss_sum_numerical()) < 5e-14 for chi in G] == [True] * euler_phi(23)
        True
        """
        cdef complex S = 0
        cdef complex x = 0
        cdef long q = self._parent.q
        for n in range(q):
            x = self.value(n)
            if x != 0:
                S = S + x * cmath.exp( (twopii * a * n) / q)

        return S

    cpdef complex gauss_sum_numerical(self, prec = 53, long a = 1):
        r"""
        This is a synonym for gauss_sum(), except that it accepts an extra
        precision argument for uniformity with the other implementation of
        Dirichlet characters. Right now higher precision is not implemented,
        however.
        """

        if prec != 53:
            raise NotImplementedError("Right now we only support a precision of 53 bits.")

        return self.gauss_sum(a)
    
    def number(self):
        return self._number

    def __invert__(self):
        r"""
        Return the inverse of this character. The inverse of a character
        `\chi` is the character `\overline \chi` such that
        `\chi(n) \overline \chi(n)` is the trivial character.

        EXAMPLE:
            
            sage: from dirichlet_conrey import *
            sage: G = DirichletGroup_conrey(23)
            sage: chi = G[3]
            sage: chi2 = chi.__invert__()
            sage: chi * chi2
            Dirichlet character with index 1 modulo 23
            sage: abs(chi2(3) * chi(3) - 1.0) < 1e-15
            True
        """
        return DirichletCharacter_conrey(self._parent, inverse_mod(self._n, self._parent.q))

    cpdef is_even(self):
        r"""
        Return ``true`` if this character is even, ``false`` otherwise.
        
        A character `\chi` is even if `\chi(n) = \chi(-n)` for all n, and odd if
        `\chi(n) = -\chi(-n)` for all ``n``. Every character is either even or odd.

        TESTS::

        sage: from dirichlet_conrey import *
        sage: G = DirichletGroup_conrey(11)
        sage: [chi.is_even() for chi in G] == [chi.sage_character().is_even() for chi in G]
        True
        """
        return self.exponent(-1) == 0

    cpdef is_odd(self):
        r"""
        Return ``true`` if this character is odd, ``false`` otherwise.
        
        A character chi is even if `\chi(n) = \chi(-n)` for all n, and odd if
        `\chi(n) = -\chi(-n)` for all ``n``. Every character is either even or odd.

        TESTS::

        sage: from dirichlet_conrey import *
        sage: G = DirichletGroup_conrey(16)
        sage: [chi.is_odd() for chi in G] == [chi.sage_character().is_odd() for chi in G]
        True
        """
        return self.exponent(-1) != 0

    def is_primitive(self):
        """
        Return whether or not this character is primitive.

        EXAMPLES::

        sage: from dirichlet_conrey import *
        sage: G = DirichletGroup_conrey(5*5*5)
        sage: G[1].is_primitive()
        False
        sage: G[2].is_primitive()
        True

        TESTS::
        sage: from dirichlet_conrey import *
        sage: [chi.is_primitive() for chi in G] == [chi.sage_character().is_primitive() for chi in G]
        True
        """

        for j in range(self._parent.k):
            if not self._is_primitive_at_known_p(j):
                return False

        return self._is_primitive_at_two()

    cdef int _is_primitive_at_known_p(self, long j):
        r"""
        Return whether or not this character is primitive at the jth prime
        factor of the odd part of the modulus.
        """

        cdef long p = self._parent.primes[j]
        n = self._n % self._parent.q_odd
        cdef long dlog = self._parent.A[n * self._parent.k + j]
        return dlog % p != 0

    cdef _is_primitive_at_two(self):
        cdef long q_even = self._parent.q_even
        cdef long * B = self._parent.B
        cdef long n = self._n % q_even
        if q_even == 1:
            return True
        elif q_even == 2:
            return False
        elif q_even == 4:
            return n == 3
        else:
            return B[n] % 2 == 1

    def is_trivial(self):
        r"""
        Return ``true`` if this character is trivial, ``false`` otherwise.

        A character `\chi` is trivial if `\chi(n) = 1` whenever it is nonzero.

        TESTS::

        sage: from dirichlet_conrey import *
        sage: G = DirichletGroup_conrey(17)
        sage: G[1].is_trivial()
        True
        sage: G[18].is_trivial()
        True
        sage: G[2].is_trivial()
        False
        """
        return self._n == 1

    def kernel(self, outtype=None):
        r"""
        Return the kernel of this character as a list. By default, this is
        a list of Python integers, but ``outtype`` can be specified to
        change this to something else.
        
        TESTS::

        sage: from dirichlet_conrey import *
        sage: G = DirichletGroup_conrey(4 * 25)
        sage: [chi.kernel() for chi in G] == [chi.sage_character().kernel() for chi in G]
        True
        """
        if outtype is not None:
            return [outtype(n) for n in range(self._parent.q) if gcd(n,self._parent.q) == 1 and self.exponent(n) == 0]
        else:
            return [n for n in range(self._parent.q) if gcd(n,self._parent.q) == 1 and self.exponent(n) == 0]

    def level(self):
        r"""
        synonym for conductor().
        """
        return self.conductor()

    def logvalue(self, long m):
        r"""
        Return log(chi(m))/(2 pi i) as a rational number; i.e., return a/b
        so that chi(m) = e(a/b). 
        """
        if gcd(m,self._parent.q) != 1:
            return -1
        cdef long exponent = self.exponent(m)
        return Integer(exponent)/Integer(self._parent.phi_q) # TODO: there is probably
                                                             # a better way to construct
                                                             # a rational number.

    def Lone(self):
        r"""
        Return the L-series of this character evaluated at 1. If the character
        is trivial, raise an error.
        """
        #
        # This implementation is not optimal, but is simple. When chi is even,
        # we use the formula
        #
        # L(1, chi) = -tau(chi)/q sum_{a=1}^{q-1} chibar(a) log(sin(pi a/q))
        #
        # (which is eq 9.8, page 288 in Montgomery-Vaughan)
        #
        # while for odd chi we use the formula
        #
        # L(1, chi) = pi * i/(tau(chibar) * (2 - chibar(2))) * sum_{a=1}^{q/2} chibar(a)

        cdef complex S = 0
        cdef long p

        if self.is_trivial():
            raise ArithmeticError

        if not self.is_primitive():
            psi = self.primitive_character()
            S = psi.Lone()
            for n in range(self._parent.k):
                p = self._parent.primes[n]
                S = S * (1 - psi(p)/p)
            if self._parent.q % 2 == 0:
                S = S * (1 - psi(2)/2)
            return S

        chibar = self.__invert__()
        cdef long q = self._parent.q
        if self.is_even():
            for a in range(1, q):
                S = S + chibar(a) * cmath.log(cmath.sin(PI * a/float(q)))
            return -S * self.gauss_sum()/q
        else:
            for a in range(1, q):
                S = S + a * chibar(a)
            return S * PI * 1.0j * self.gauss_sum()/(q*q)
                
    cpdef max_sum(self, return_location = False):
        """
        return the maximum of the character sum and its location.
        """
        
        cdef complex S = 0
        cdef float absmax = 0
        cdef complex max = 0
        cdef long m = 0
        for n in range(self.q/2 + 1):
            S = S + self.value(n)
            if abs(S) > absmax:
                max = S
                absmax = abs(S)
                m = n

        return max, m
            
    def modulus(self):
        r"""
        Return the modulus of this charater as a Python integer.

        EXAMPLES:

            sage: from dirichlet_conrey import *
            sage: G = DirichletGroup_conrey(68)
            sage: G[7].modulus()
            68
        """
        return self._parent.q

    def __mul__(self, other):
        r"""
        Return the product of this character and another.

        If the two characters have the same modulus, then we just multiply
        the indices; otherwise we map the characters to a common modulus and
        do the multiplication there. (Second part not implemented yet.)

        EXAMPLES:

            sage: from dirichlet_conrey import *
            sage: G = DirichletGroup_conrey(21)
            sage: chi = G[4]
            sage: chi2 = G[5]
            sage: chi * chi2
            Dirichlet character with index 20 modulo 21
            sage: G2 = DirichletGroup_conrey(17)
            sage: chi * chi2 * G2[7]
            Dirichlet character with index 41 modulo 357
        
        TESTS:

            sage: chi.sage_character().extend(357).maximize_base_ring() * chi2.sage_character().extend(357).maximize_base_ring() * G2[7].sage_character().extend(357).maximize_base_ring() == (chi * chi2 * G2[7]).sage_character()
            True
        """
        return self._mul_(other)
     
    cpdef DirichletCharacter_conrey _mul_(self, DirichletCharacter_conrey other):
        cdef long q, q1, q2, n1, n2, n
        q1 = self._parent.q
        n1 = self._n
        q2 = other._parent.q
        n2 = other._n
        if q1 == q2:
            return DirichletCharacter_conrey(self._parent, (n1 * n2) % q1)
        elif gcd(q1, q2) == 1:
            q = q1 * q2
            n = crt(n1, n2, q1, q2)
            return DirichletGroup_conrey(q)[n]
        else:
            raise NotImplementedError

    def multiplicative_order(self):
        r"""
        Return the order of this character as an element of the Dirichlet
        group that it is a member of. The set of characters modulo q forms an
        abelian group which is isomorphic to `(Z/qZ)^*`, so every element has
        finite order.

        EXAMPLES::

        sage: from dirichlet_conrey import *
        sage: G = DirichletGroup_conrey(11)
        sage: G[1].multiplicative_order()
        1
        sage: G[-1].multiplicative_order()
        2

        TESTS::

        sage: from dirichlet_conrey import *
        sage: G = DirichletGroup_conrey(4 * 49)
        sage: [chi.multiplicative_order() for chi in G] == [chi.sage_character().multiplicative_order() for chi in G]
        True
        """
        
        return multiplicative_order(mod(self._n, self._parent.q))
        
    cpdef parent(self):
        return self._parent

    def primitive_character(self):
        r"""
        Return the primitive character that induces this character.

        TESTS::

        sage: from dirichlet_conrey import *
        sage: G = DirichletGroup_conrey(17)
        sage: G[1].primitive_character()
        Dirichlet character with index 0 modulo 1
        sage: G[18].primitive_character()
        Dirichlet character with index 0 modulo 1
        sage: G[2].primitive_character()
        Dirichlet character with index 2 modulo 17
        sage: G = DirichletGroup_conrey(1)
        sage: [chi.primitive_character().sage_character() for chi in G] == [chi.sage_character().primitive_character() for chi in G]
        True
        sage: G = DirichletGroup_conrey(2)
        sage: [chi.primitive_character().sage_character() for chi in G] == [chi.sage_character().primitive_character() for chi in G]
        True
        sage: G = DirichletGroup_conrey(4)
        sage: [chi.primitive_character().sage_character() for chi in G if not chi.is_trivial()] == [chi.sage_character().primitive_character() for chi in G if not chi.is_trivial()] # omitting trivial character because of a bug in sage
        True
        sage: G = DirichletGroup_conrey(8)
        sage: [chi.primitive_character().sage_character() for chi in G if not chi.is_trivial()] == [chi.sage_character().primitive_character() for chi in G if not chi.is_trivial()] # omitting trivial character because of a bug in sage
        True
        sage: G = DirichletGroup_conrey(16)
        sage: [chi.primitive_character().sage_character() for chi in G] == [chi.sage_character().primitive_character() for chi in G]
        True
        sage: G = DirichletGroup_conrey(8 * 5)
        sage: [chi.primitive_character().sage_character() for chi in G] == [chi.sage_character().primitive_character() for chi in G]
        True
        sage: G = DirichletGroup_conrey(8 * 5 * 5)
        sage: [chi.primitive_character().sage_character() for chi in G] == [chi.sage_character().primitive_character() for chi in G]
        True
        sage: G = DirichletGroup_conrey(16 * 3 * 3)
        sage: [chi.primitive_character().sage_character() for chi in G] == [chi.sage_character().primitive_character() for chi in G]
        True
        """
        odd_parts = [self._primitive_part_at_known_p(k) for k in range(self._parent.k)]
        even_index, even_conductor = self._primitive_part_at_two()
        indices = [Integer(even_index)] + [Integer(n) for (n,m) in odd_parts]
        moduli = [Integer(even_conductor)] + [Integer(m) for (n,m) in odd_parts]
        q = prod(moduli)
        index = crt(indices, moduli)
        return DirichletGroup_conrey(q)[index]

    cdef tuple _primitive_part_at_known_p(self, long j):
        r"""
        Return the conductor and the index of the primitive character
        associated to the p-part of the character for the j-th odd prime
        factor p of the modulus.
        """

        cdef long e, conductor, index
        cdef long p = self._parent.primes[j]
        n = self._n % self._parent.q_odd
        cdef long dlog = self._parent.A[n * self._parent.k + j]
        if dlog == 0:
            return (1, 1)
        else:
            e = 0
            while dlog % p == 0:
                dlog /= p
                e = e + 1
            conductor = p**(self._parent.exponents[j] - e)
            index = power_mod(self._parent.generators[j], dlog, conductor)
            return (index, conductor)

    cdef _primitive_part_at_two(self):
        r"""
        Return the conductor and the index of the primitive
        character associated to the even part of the modulus.

        We return a pair (n, q), where q is the even part of the
        conductor of chi and n in the index of the even part of the
        primitive character inducing chi.
        """
        cdef long q_even = self._parent.q_even
        cdef long * B = self._parent.B
        cdef long n = self._n % q_even
        cdef long e, dlog

        # We basically do everything by hand.

        if q_even == 1:                 # If q is odd, there is no even
            return (1,1)                # part of the character,
        
        elif q_even == 2:               # while if 2 is only divisible by 2
            return (1,1)                # the inducing character must be
                                        # trivial.

        elif q_even == 4:               # When q_even is 4, the conductor
            if n == 3:                  # just depends on n mod 4
                return (3,4)
            else:
                return (1,1)
        elif q_even == 8:
            if n == 1:
                return (1,1)
            elif n == 3:
                return (3,8)
            elif n == 5:
                return (5,8)
            else:
                return (3,4)
        else:
            #if n == 1:
            #    return (1,1)
            #elif n == q_even - 1:
            #   return (7, 8) # special case for primitive character mod 8
            #elif n == q_even/2 - 1:
            #    return (3,4)  # special case for primitive character mod 4
            alpha = B[n]
            epsilon = B[n-1]
            if alpha == 0:
                if epsilon == 1:
                    return (1,1)
                else:
                    return (3,4)
            f = 0
            m = alpha
            while m % 2 == 0:
                m /= 2
                f = f + 1
            conductor = q_even/(2**f)
            index = power_mod(5, m, conductor)
            if epsilon == -1:
                index = conductor - index
            return (index, conductor)

    def __repr__(self):
        return "Dirichlet character with index %d modulo %d" % (self._n, self._parent.q)

    def restrict(self):
        r"""
        Return the "restriction" of this character to a smaller modulus.

        If `\chi` is a character modulo `q` and `q'` divides `q`, then we
        define the restriction as...
        """

        raise NotImplementedError

    def sage_character(self):
        """
        return the sage.modular.dirichlet.DirichletCharacter that corresponds
        to this character.

        This function has a stupid name, because eventually this code should
        be part of Sage, so there will be two available implementations. I
        don't know what to call it right now.
        """

        G = self._parent.standard_dirichlet_group()

        gens = G.unit_gens()  # grabbing the generators this way
                              # ensures that they will be the same
                              # generators used by the
                              # DirichletGroup_class

        # We can construct a DirichletCharacter by giving it
        # a list of the exponents on the generators, so we
        # compute these.
        
        # Because we are computing the odd and even parts of
        # the character separately, we have to properly combine
        # the exponents.

        cdef long zeta_order = G.zeta_order()
        exponents = []
        for a in gens:
            exponent = (self.exponent(a) * zeta_order)/self._parent.phi_q

            exponents.append(exponent)

        # To make sure that the exponent vector has the right type, I'm
        # mimicking a bit what is done in modular/dirichlet.pyx, without
        # necessarily understanding if what I'm doing it right.
        # 
        # Thus I put the XXX here to mark this as a possible trouble
        # spot if there are problems in the future...

        M = self._parent.standard_dirichlet_group()._module

        exponents = M(exponents)
        return DirichletCharacter(self._parent.standard_dirichlet_group(), exponents)

    cpdef sum(self, long m):
        """
        return the sum of chi up to n
        """
        cdef complex S = 0
        for n in range(m + 1):
            S = S + self.value(n)

        return S

    cpdef complex value(self, long m):
        return self._parent.chi(self._n, m)


    def values(self):
        return [self.value(n) for n in range(self._parent.q)]

    #cdef complex __call__unsafe(self, long m):
    #    return self.parent._chi_unsafe(self._n, m)
