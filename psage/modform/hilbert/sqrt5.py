"""
EXAMPLES:

LEVEL 31::

    sage: from psage.modform.hilbert.sqrt5 import *    
    sage: F.<a> = QuadraticField(5)
    sage: B.<i,j,k> = QuaternionAlgebra(F,-1,-1)
    sage: c = F.factor(31)[1][0]
    sage: P = F.primes_above(5)[0]
    sage: TH = THETA(B, 20)   # currently about a minute
    pi = [2, 3, -1/2*a + 5/2, -a + 4, -1/2*a + 7/2, -1/2*a + 9/2, -3/2*a + 11/2]
    tr = [4, 6, 5, 8, 7, 9, 11]
    N = 11
    Sorting through 22440 elements
    sage: T = hecke_ops(B, c, TH)
    sage: T
    [(5, -1/2*a + 5/2, [3 3]
    [5 1]), (9, 3, [7 3]
    [5 5]), (11, -1/2*a + 7/2, [9 3]
    [5 7]), (11, -a + 4, [ 6  6]
    [10  2]), (19, -1/2*a + 9/2, [11  9]
    [15  5]), (19, -3/2*a + 11/2, [14  6]
    [10 10])]
    sage: for nm,p,t in T:
    ...       print nm, p, t.charpoly().factor()
    5 -1/2*a + 5/2 (x - 6) * (x + 2)
    9 3 (x - 10) * (x - 2)
    11 -1/2*a + 7/2 (x - 12) * (x - 4)
    11 -a + 4 (x - 12) * (x + 4)
    19 -1/2*a + 9/2 (x - 20) * (x + 4)
    19 -3/2*a + 11/2 (x - 20) * (x - 4)


LEVEL 41::

    sage: from psage.modform.hilbert.sqrt5 import * 
    sage: F.<a> = QuadraticField(5)
    sage: B.<i,j,k> = QuaternionAlgebra(F,-1,-1)
    sage: F.primes_above(41)
    [Fractional ideal (1/2*a - 13/2), Fractional ideal (1/2*a + 13/2)]
    sage: c = F.primes_above(41)[0]
    sage: TH = THETA(B, 11)   # about 30 seconds
    pi = [2, 3, -1/2*a + 5/2, -a + 4, -1/2*a + 7/2]
    tr = [4, 6, 5, 8, 7]
    N = 8
    Sorting through 6660 elements
    sage: T = hecke_ops(B, c, TH)
    sage: T
    [(5, -1/2*a + 5/2, [4 2]
    [5 1]), (9, 3, [ 6  4]
    [10  0]), (11, -a + 4, [10  2]
    [ 5  7]), (11, -1/2*a + 7/2, [ 8  4]
    [10  2])]
    sage: for nm,p,t in T:
    ...         print nm, p, t.charpoly().factor()
    5 -1/2*a + 5/2 (x - 6) * (x + 1)
    9 3 (x - 10) * (x + 4)
    11 -a + 4 (x - 12) * (x - 5)
    11 -1/2*a + 7/2 (x - 12) * (x + 2)

LEVEL 389!:

This relies on having TH from above (say from the level 31 block above)::
    
    sage: c = F.primes_above(389)[0]
    sage: T = hecke_ops(B, c, TH)
    sage: for nm,p,t in T:
    ...       print nm, p, t.charpoly().factor()
    5 -1/2*a + 5/2 (x - 6) * (x^2 + 4*x - 1) * (x^2 - x - 4)^2
    9 3 (x - 10) * (x^2 + 3*x - 9) * (x^4 - 5*x^3 + 3*x^2 + 6*x - 4)
    11 -1/2*a + 7/2 (x - 12) * (x^2 + 5*x + 5) * (x^4 - x^3 - 23*x^2 + 18*x + 52)
    11 -a + 4 (x - 12) * (x + 3)^2 * (x^4 - 17*x^2 + 68)
    19 -1/2*a + 9/2 (x - 20) * (x^2 + 3*x - 9) * (x^4 + x^3 - 23*x^2 + 16*x + 52)
    19 -3/2*a + 11/2 (x - 20) * (x^2 + 3*x - 9) * (x^4 + 5*x^3 - 65*x^2 - 278*x + 404)
    sage: F.primes_above(389)
    [Fractional ideal (9*a + 4), Fractional ideal (-9*a + 4)]
"""

def modp_splitting(B, p):
    """
    INPUT:
        
        - B -- quaternion algebra of the form K<i,j> where i^2=a, j^2=b.
        - p -- ideal of the number field K = B.base() with ring O of integers.
        
    OUTPUT:
        
        - matrices I, J in M_2(O/p) such that i |--> I and j |--> J defines 
          an algebra morphism, i.e., I^2=a, J^2=b, I*J=-J*I.

    EXAMPLES::    
                        
        sage: F.<a> = QuadraticField(5); F
        Number Field in a with defining polynomial x^2 - 5
        sage: B.<i,j,k> = QuaternionAlgebra(F,-1,-1); B
        Quaternion Algebra (-1, -1) with base ring Number Field in a with defining polynomial x^2 - 5
        sage: c = F.factor(31)[0][0]
        sage: from psage.modform.hilbert.sqrt5 import modp_splitting
        sage: modp_splitting(B, c)
        (
        [ 0 30]  [18  4]
        [ 1  0], [ 4 13]
        )
        sage: c = F.factor(37)[0][0]; c
        Fractional ideal (37)
        sage: I, J = modp_splitting(B, c); I, J
        (
        [ 0 36]  [ 7*abar + 15 21*abar + 32]
        [ 1  0], [21*abar + 32 30*abar + 22]
        )
        sage: I^2
        [36  0]
        [ 0 36]
        sage: J^2
        [36  0]
        [ 0 36]
        sage: I*J == -J*I
        True

    AUTHOR: William Stein
    """
    # Inspired by the code in the function
    # modp_splitting_data in algebras/quatalg/quaternion_algebra.py
    if p.number_field() != B.base():
        raise ValueError, "p must be a prime ideal in the base field of the quaternion algebra"
    if not p.is_prime():
        raise ValueError, "p must be prime"
    F = p.residue_field()
    from sage.all import MatrixSpace
    M = MatrixSpace(F, 2)
    i2, j2 = B.invariants()
    i2 = F(i2); j2 = F(j2)    
    if F.characteristic() == 2:
        if i2 == 0 or j2 == 0:
            raise NotImplementedError
        return M([0,1,1,0]), M([1,1,0,1])
    # Find I -- just write it down
    I = M([0,i2,1,0])
    # Find J -- I figured this out by just writing out the general case
    # and seeing what the parameters have to satisfy
    i2inv = 1/i2
    a = None
    for b in list(F):
        if not b: continue
        c = j2 + i2inv * b*b
        if c.is_square():
            a = -c.sqrt()
            break        
    if a is None:
        # do a fallback search; needed in char 3 sometimes.
        for J in M:
            K = I*J
            if J*J == j2 and K == -J*I:
                return I, J, K
            
    J = M([a,b,(j2-a*a)/b, -a])
    K = I*J
    assert K == -J*I, "bug in that I,J don't skew commute"    
    return I, J

def modp_splitting_map(B, p):
    """
    Return a map from subset of B to 2x2 matrix space isomorphic
    to R tensor OF/p.

    INPUT:
        - `B` -- quaternion algebra over F=Q(sqrt(5))
          with invariants -1, -1.
        - `p` -- prime ideal of F=Q(sqrt(5))

    EXAMPLES::

        sage: from psage.modform.hilbert.sqrt5 import *    
        sage: attach /tmp/a.sage
        sage: F.<a> = QuadraticField(5)
        sage: B.<i,j,k> = QuaternionAlgebra(F,-1,-1)
        sage: theta = modp_splitting_map(B, F.primes_above(5)[0])
        sage: theta(i+j-k)
        [2 1]
        [3 3]
        sage: s = 2 + 3*i - 2*j - 2*k
        sage: theta(s)
        [1 3]
        [4 3]
        sage: s.reduced_characteristic_polynomial()
        x^2 - 4*x + 21
        sage: theta(s).charpoly()
        x^2 + x + 1
        sage: s.reduced_characteristic_polynomial().change_ring(GF(5))
        x^2 + x + 1
        sage: theta = modp_splitting_map(B, F.primes_above(3)[0])
        sage: smod = theta(s); smod
        [  abar + 2     2*abar]
        [    2*abar 2*abar + 2]
        sage: smod^2 - 4*smod + 21
        [0 0]
        [0 0]    
    """
    I, J = modp_splitting(B,p)
    F = I.parent().base_ring()
    def f(x):
        return F(x[0]) + I*F(x[1]) + J*F(x[2]) + I*J*F(x[3])
    return f

def icosian_gens(B):
    """
    Return generators of the icosian group, as elements of the 
    Hamilton quaternion algebra B over Q(sqrt(5)).
    
    AUTHOR: William Stein

    EXAMPLES::
    
        sage: F.<a> = QuadraticField(5)
        sage: B.<i,j,k> = QuaternionAlgebra(F,-1,-1)
        sage: from psage.modform.hilbert.sqrt5 import icosian_gens
        sage: icosian_gens(B)
        [i, j, k, -1/2 + 1/2*i + 1/2*j + 1/2*k, 1/2*i + (-1/4*a + 1/4)*j + (1/4*a + 1/4)*k]
        sage: [a.reduced_norm() for a in icosian_gens(B)]
        [1, 1, 1, 1, 1]
    """
    F = B.base()
    sqrt5 = F.gen(); assert sqrt5*sqrt5==5
    sigma = (1-sqrt5)/2
    tau = (1+sqrt5)/2
    return [B(v)/2 for v in [(0,2,0,0), (0,0,2,0), 
                 (0,0,0,2), (-1,1,1,1), (0,1,tau,sigma)]]
    # NOTE: In a previous version of this code (0,1,tau,sigma) was
    # accidentally replaced by (0,1,sigma,tau), which is wrong (odd
    # permutation!), and led to much work to fix.

def icosian_ring_gens(B):
    """
    Return generators for the icosian ring (a maximal order) in the
    quaternion algebra ramified only at infinity over F=Q(sqrt(5)).
    These are generators over the ring of integers of F.
    """
    # See page 6 of Dembele.
    sqrt5 = B.base().gen(); assert sqrt5*sqrt5==5
    omega = (1+sqrt5)/2
    omega_bar = (1-sqrt5)/2
    return [B(v)/2 for v in [(1,-omega_bar,omega,0),
                             (0,-omega_bar,1,omega),
                             (0,omega,-omega_bar,1),
                             (0,1,omega,-omega_bar)]]

def icosian_ring_gens_over_ZZ(B):
    """
    Return basis over ZZ for the icosian ring, which has ZZ-rank 8.
    """
    I = icosian_ring_gens(B)
    sqrt5 = B.base().gen(); assert sqrt5*sqrt5==5    
    omega = (1+sqrt5)/2
    return I + [omega*x for x in I]

def tensor_over_QQ_with_RR(B, prec=53):
    """
    Return map from the quaternion algebra B to the tensor product of
    B over QQ with RR, viewed as an 8-dimensional real vector space.
    """
    from sage.all import RealField
    RR = RealField(prec=prec)
    V = RR**8
    F = B.base()
    S = F.embeddings(RR)
    def f(x):
        return V(sum([[sigma(a) for a in x] for sigma in S],[]))
    return f

def modp_icosians(B, p):
    """
    Return matrices of images of all 120 icosians mod p.
    """
    I, J = modp_splitting(B, p); K = I*J
    k = p.residue_field()
    G = [k(g[0]) + k(g[1])*I + k(g[2])*J + k(g[3])*K for g in icosian_gens(B)]
    from sage.all import MatrixGroup
    return [g.matrix() for g in MatrixGroup(G)]

class P1ModList(object):
    def __init__(self, c):
        """
        INPUT:
           - c -- a prime of O_F, where F is the totally real field Q(sqrt(5)).
        """
        self._c = c

        F = c.residue_field()
        V = F**2
        self._V = V
        self._F = F
        self._list = [ V([0,1]) ] + [V([1,a]) for a in F]
        for a in self._list:
            a.set_immutable()

    def random_element(self):
        import random
        return random.choice(self._list)
        
    def normalize(self, uv):
        w = self._V(uv)
        if w[0]:
            w = (~w[0]) * w
            w.set_immutable()
            #assert w in self._list
            return w
        else:
            return self._list[0]

    def __len__(self):
        return len(self._list)
    
    def __getitem__(self, i):
        return self._list[i]

    def __call__(self, x):
        return self.normalize(x)

    def __repr__(self):
        return 'Projective line over %s'%self._F
    

def P1_orbits(B, p):
    """
    INPUT:

       - B -- quaternion algebra
       - p -- a prime of O_F, where F is the totally real field Q(sqrt(5)).

    AUTHOR: William Stein

    EXAMPLES::

    """
    P1 = P1ModList(p)
    ICO = modp_icosians(B, p)
    
    def act(u, t):
        return P1(u*t)

    cur = P1.random_element()
    reps = [cur]
    orbits = {cur:cur}
    while len(orbits) < len(P1):
        for u in ICO:
            s = act(u, cur)
            if not orbits.has_key(s):
                orbits[s] = cur
        if len(orbits) < len(P1):
            # choose element of P1 not a key
            while True:
                c = P1.random_element()
                if c not in orbits:
                    cur = c
                    reps.append(cur)
                    orbits[cur] = cur
                    break
    # done
    return orbits, reps, P1

def P1_orbits2(B, p):
    """
    INPUT:

       - B -- quaternion algebra
       - p -- a prime of O_F, where F is the totally real field Q(sqrt(5)).

    AUTHOR: William Stein

    EXAMPLES::

    """
    P1 = P1ModList(p)
    ICO = modp_icosians(B, p)
    
    orbits = []
    while sum(len(x) for x in orbits) < len(P1):
        v = P1.random_element()
        skip = False
        for O in orbits:
            if v in O:
                skip = True
                break
        if skip: continue
        O = set([P1(g*v) for g in ICO])
        orbits.append(O)

    # Now make a dictionary
    return orbits
    
    return orbits, reps, P1

def totally_pos_gen(p):
    """
    Given a prime ideal p of a narrow class number 1 real quadratic
    field, return a totally positive generator.

    INPUT:
        - p -- prime ideal of narrow class number 1 real
          quadratic field

    OUTPUT:
        - generator of p that is totally positive

    AUTHOR: William Stein

    EXAMPLES::
    
        sage: F.<a> = QuadraticField(5)
        sage: from psage.modform.hilbert.sqrt5 import totally_pos_gen
        sage: g = totally_pos_gen(F.factor(19)[0][0]); g
        -1/2*a + 9/2
        sage: g.complex_embeddings()
        [5.61803398874989, 3.38196601125011]        

        sage: for p in primes(14):
        ...       for P, e in F.factor(p):
        ...           g = totally_pos_gen(P)
        ...           print P, g, g.complex_embeddings()
        Fractional ideal (2) 2 [2.00000000000000, 2.00000000000000]
        Fractional ideal (3) 3 [3.00000000000000, 3.00000000000000]
        Fractional ideal (a) -1/2*a + 5/2 [3.61803398874989, 1.38196601125011]
        Fractional ideal (7) 7 [7.00000000000000, 7.00000000000000]
        Fractional ideal (-3/2*a + 1/2) -a + 4 [6.23606797749979, 1.76393202250021]
        Fractional ideal (-3/2*a - 1/2) -1/2*a + 7/2 [4.61803398874989, 2.38196601125011]
        Fractional ideal (13) 13 [13.0000000000000, 13.0000000000000]
    """
    F = p.number_field()
    assert F.degree() == 2 and F.discriminant() > 0
    G = p.gens_reduced()
    if len(G) != 1:
        raise ValueError, "ideal not principal"
    g = G[0]

    from sage.all import RR
    sigma = F.embeddings(RR)
    
    e = [s(g) for s in sigma]
    u = F.unit_group().gen(1)
    ue = [s(u) for s in sigma]
    if ue[0] > 0 and ue[1] < 0:
        u *= -1
    if e[0] < 0 and e[1] < 0:
        return -g
    elif e[0] < 0 and e[1] > 0:
        if ue[0] < 0 and ue[1] > 0:
            return u*g
        else:
            raise ValueError, "no totally positive generator"
    elif e[0] > 0 and e[1] > 0:
        return g
    elif e[0] > 0 and e[1] < 0:
        if ue[0] < 0 and ue[1] > 0:
            return -u*g
        else:
            raise ValueError, "no totally positive generator"
    assert False, "bug"

def gram_matrix_of_maximal_order(R):
    """
    Return 8x8 Gram matrix of maximal order R.
    """
    G = [[(R[i]*R[j].conjugate()).reduced_trace().trace()
          for i in range(8)] for j in range(8)]
    from sage.all import matrix, ZZ
    return matrix(ZZ, G)

def bounded_elements(B, N):
    """
    Return elements in maximal order of B that have reduced norm
    whose trace is at most N.

    EXAMPLES::

        sage: from psage.modform.hilbert.sqrt5 import *    
        sage: F.<a> = QuadraticField(5)
        sage: B.<i,j,k> = QuaternionAlgebra(F,-1,-1)
        sage: X = bounded_elements(B,3)
        sage: len(X)
        180
        sage: rnX = [a.reduced_norm() for a in X]
        sage: set([a.trace() for a in rnX])
        set([2, 3])
        sage: set([a.norm() for a in rnX])
        set([1])
        sage: X = bounded_elements(B,5)
        sage: len(X)
        1200
        sage: rnX = [a.reduced_norm() for a in X]
        sage: set([a.trace() for a in rnX])
        set([2, 3, 4, 5])
        sage: set([a.norm() for a in rnX])
        set([1, 4, 5])    
    """
    # Get our maximal order
    R = icosian_ring_gens_over_ZZ(B)
    
    # Compute Gram matrix of R
    G = gram_matrix_of_maximal_order(R)

    # Make PARI quadratic form
    from sage.all import pari
    qf = pari(G)

    # Get the vectors of norm up to N.
    # The 2 is because we had to scale by 2 to get
    # rid of denominator in Gram matrix. 
    Z = qf.qfminim(2*N, 2**32)  # TODO: not sure about 2^32...?
    Z2 = Z[2].sage().transpose()

    # For each vector, make the corresponding element of B.
    # TODO: This step massively dominates the runtime, and can be
    # easily made trivial with careful thought.
    V = []
    for i in range(Z2.nrows()):
        w = Z2[i]
        V.append(sum(w[j]*R[j] for j in range(8)))
    return V

def primes_of_bounded_norm(F, N):
    """
    Return the primes of the quadratic field F = Q(sqrt(5))
    of norm up to N.

    EXAMPLES::

        sage: from psage.modform.hilbert.sqrt5 import *
        sage: F.<a> = QuadraticField(5)
        sage: primes_of_bounded_norm(F, 5)
        [Fractional ideal (2), Fractional ideal (a)]         
        sage: primes_of_bounded_norm(F,25)
        [Fractional ideal (2), Fractional ideal (3), Fractional ideal (-a),
         Fractional ideal (-3/2*a + 1/2), Fractional ideal (-3/2*a - 1/2),
         Fractional ideal (-2*a - 1), Fractional ideal (-2*a + 1)]
    """
    # The answer is the set of primes over primes p of ZZ
    # with p<=N with p split (or ramified) or p^2<=N with p inert.
    X = []
    from sage.all import primes
    for p in primes(N+1):
        fac = F.factor(p)
        if len(fac) == 2: # split
            X.append(fac[0][0])
            X.append(fac[1][0])
        elif len(fac) == 1 and fac[0][1]==2: # ramified
            X.append(fac[0][0])
        elif p*p <= N: # inert
            X.append(fac[0][0])
    return X

def THETA(B, N):
    r"""
    Return representative elements of the maximal order of `R` of norm
    `\pi_p` up to `N` modulo the left action of the units of `R` (the
    icosians).  Here `\pi_p` runs through totally positive generators
    of the prime ideals with norm up to `N`.

    INPUT:
       - `B` -- quaternion algebra
       - `N` -- a positive integer

    AUTHOR: William Stein

    EXAMPLES::
    """
    # ** NOTE: This algorithm will not scale up well, because there
    #    are so many vectors of bounded norm. 
    ####################################################################
    # Algorithm:
    #   * Enumerate set S of primes of norm <= N.
    #   * For each prime p in S:
    #        - Find a totally positive generator pi_p for p.
    #        - Compute mod-p local splitting theta_p.
    #   * Compute set X of elements in maximal order of B of norm <= N
    #     using the Gram matrix of the icosian ring (maximal order).
    #   * For each element z of X, compute the reduced norm of z.
    #     If it equals pi_p for one of the pi_p, compute the
    #     top row v of the reduced row echelon form of theta_p(z).
    #     Store v:z if there isn't already something for v.
    #     Also, of course, store this stuff separately for each p.
    ####################################################################

    # Get primes of norm up to N.
    F = B.base_ring()
    S = primes_of_bounded_norm(F, N)

    # Find totally positive generators pi_p
    pi = [totally_pos_gen(p) for p in S]
    print "pi =",pi

    # Compute traces of the generators, since that's what
    # the bounded_elements command computes up to.
    tr = [abs(x.trace()) for x in pi]
    print "tr =", tr
    N = max(tr)
    print "N =", N

    # A list that at least includes all elements (up to -1) whose
    # reduced norm has trace at most N.
    X = bounded_elements(B, N)
    
    # Compute mod-p local splitting maps
    theta_map = {}
    for i, p in enumerate(S):
        theta_map[pi[i]] = modp_splitting_map(B, p)

    # Sort through the elements of X.
    pi_set = set(pi)
    
    # TODO: We skip the prime 2, since the mod-p splitting map is
    # broken there.
    pi_set.remove(2)
    
    Theta = {}
    for pi_p in pi_set:
        Theta[pi_p] = {}

    # The dictionary Theta will have keys the pi_p and
    # the dictionary for pi_p has keys reduced vectors
    # in (OF/p)^2.  Here "reduced" just means "reduced
    # row echelon form", so scaled so first entry is 1.
    print "Sorting through %s elements"%len(X)
    for a in X:
        nrm = a.reduced_norm()
        if nrm in pi_set:
            # this is: mod right action of R^* acting on the right,
            # so column echelon form
            v = theta_map[nrm](a).transpose().echelon_form()[0]

            ## for reference, this is: mod left action of R^*,
            ## which is wrong, I guess:
            # v = theta_map[nrm](a).echelon_form()[0]
            
            z = Theta[nrm]
            if z.has_key(v):
                pass
            else:
                z[v] = a
                
    return Theta

def hecke_ops(B, c, X):
    orbits, reps, P1 = P1_orbits(B, c)
    theta_c = modp_splitting_map(B, c)    
    def Tp(pi):
        z = X[pi]
        mat = []
        for x in reps:
            row = [0]*len(reps)
            for _, w in z.iteritems():
                w_c = theta_c(w)
                y = w_c**(-1) * x
                y_red = orbits[P1.normalize(y)]
                row[reps.index(y_red)] += 1
            mat.append(row)
        from sage.all import ZZ, matrix
        return matrix(ZZ, mat)
    ans = [(pi.norm(), pi, Tp(pi)) for pi in X.keys()]
    ans.sort()
    return ans


def hecke_ops2(B, c, X):
    reduce, reps, P1 = P1_orbits(B, c)
    theta_c = modp_splitting_map(B, c)    
    def Tp(pi):
        z = X[pi]
        mat = []
        for x in reps:
            print "x = %s,  card = %s"%(x, len([M for M in reduce.keys() if reduce[M]==x]))
            row = [0]*len(reps)
            for _, w in z.iteritems():
                w_c = theta_c(w)
                y = w_c**(-1) * x
                print "y =", y
                y_red = reduce[P1(y)]
                row[reps.index(y_red)] += 1
                print "y_red =", y_red
            mat.append(row)
        from sage.all import ZZ, matrix
        return matrix(ZZ, mat)
    ans = [(pi.norm(), pi, Tp(pi)) for pi in X.keys()]
    ans.sort()
    return ans

