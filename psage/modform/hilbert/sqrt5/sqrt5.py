#################################################################################
#
# (c) Copyright 2010 William Stein
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
Toy Implementation of Hilbert Modular Forms

This file contains an incredibly slow naive toy implementation of
Dembele's quaternion algebra algorithm for computing Hilbert modular
forms of weight (2,2) and ramified or split prime level.  This is for
testing and educational purposes only.  The file sqrt5_fast.pyx
contains a dramatically faster version.  That said, figuring out the
content of this file based on the contents of Dembele's paper
"Explicit Computations of Hilbert Modular Forms on Q(sqrt(5))"
was a timing consuming and very painful task.

EXAMPLES:

LEVEL 31::

    sage: from psage.modform.hilbert.sqrt5.sqrt5 import THETA, hecke_ops, F
    sage: B.<i,j,k> = QuaternionAlgebra(F,-1,-1)
    sage: c = F.factor(31)[1][0]
    sage: P = F.primes_above(5)[0]
    sage: TH = THETA(20)        # about 1 minute
    pi = [2, a + 2, 3, 2*a + 3, a + 3, 3*a + 4, a + 4]
    Sorting through 22440 elements
    sage: T = hecke_ops(c, TH); T   # random output do to choice of basis
    [(5, a + 2, [1 5]
    [3 3]), (9, 3, [5 5]
    [3 7]), (11, a + 3, [ 2 10]
    [ 6  6]), (11, 2*a + 3, [7 5]
    [3 9]), (19, a + 4, [10 10]
    [ 6 14]), (19, 3*a + 4, [ 5 15]
    [ 9 11])]
    sage: for nm,p,t in T:
    ...       print nm, p, t.charpoly().factor()
    5 a + 2 (x - 6) * (x + 2)
    9 3 (x - 10) * (x - 2)
    11 a + 3 (x - 12) * (x + 4)
    11 2*a + 3 (x - 12) * (x - 4)
    19 a + 4 (x - 20) * (x - 4)
    19 3*a + 4 (x - 20) * (x + 4)

LEVEL 41::

    sage: from psage.modform.hilbert.sqrt5.sqrt5 import THETA, hecke_ops, F
    sage: B.<i,j,k> = QuaternionAlgebra(F,-1,-1)
    sage: F.primes_above(41)
    [Fractional ideal (a - 7), Fractional ideal (a + 6)]
    sage: c = F.primes_above(41)[0]
    sage: TH = THETA(11)    # about 30 seconds
    pi = [2, a + 2, 3, 2*a + 3, a + 3]
    Sorting through 6660 elements
    sage: T = hecke_ops(c, TH); T   # random output do to choice of basis
    [(5, a + 2, [4 2]
    [5 1]), (9, 3, [ 6  4]
    [10  0]), (11, a + 3, [10  2]
    [ 5  7]), (11, 2*a + 3, [ 8  4]
    [10  2])]
    sage: for nm,p,t in T:
    ...         print nm, p, t.charpoly().factor()
    5 a + 2 (x - 6) * (x + 1)
    9 3 (x - 10) * (x + 4)
    11 a + 3 (x - 12) * (x - 5)
    11 2*a + 3 (x - 12) * (x + 2)
    

LEVEL 389!:

This relies on having TH from above (say from the level 31 block above)::
    
    sage: F.primes_above(389)
    [Fractional ideal (18*a - 5), Fractional ideal (-18*a + 13)]
    sage: c = F.primes_above(389)[0]
    sage: T = hecke_ops(c, TH)
    sage: for nm,p,t in T:
    ...       print nm, p, t.charpoly().factor()
    5 a + 2 (x - 6) * (x^2 + 4*x - 1) * (x^2 - x - 4)^2
    9 3 (x - 10) * (x^2 + 3*x - 9) * (x^4 - 5*x^3 + 3*x^2 + 6*x - 4)
    11 a + 3 (x - 12) * (x + 3)^2 * (x^4 - 17*x^2 + 68)
    11 2*a + 3 (x - 12) * (x^2 + 5*x + 5) * (x^4 - x^3 - 23*x^2 + 18*x + 52)
"""


from sage.all import NumberField, polygen, QQ, ZZ, QuaternionAlgebra, cached_function, disk_cached_function

x = polygen(QQ,'x')
F = NumberField(x**2 - x -1, 'a')
O_F = F.ring_of_integers()
B = QuaternionAlgebra(F,-1,-1,'i,j,k')

def modp_splitting(p):
    """
    INPUT:
        
        - p -- ideal of the number field K = B.base() with ring O of integers.
        
    OUTPUT:
        
        - matrices I, J in M_2(O/p) such that i |--> I and j |--> J defines 
          an algebra morphism, i.e., I^2=a, J^2=b, I*J=-J*I.

    EXAMPLES::    

        sage: from psage.modform.hilbert.sqrt5.sqrt5 import F, B, modp_splitting
        sage: c = F.factor(31)[0][0]
        sage: modp_splitting(c)
        (
        [ 0 30]  [18  4]
        [ 1  0], [ 4 13]
        )
        sage: c = F.factor(37)[0][0]; c
        Fractional ideal (37)
        sage: I, J = modp_splitting(c); I, J
        (
        [ 0 36]  [23*abar + 21  36*abar + 8]
        [ 1  0], [ 36*abar + 8 14*abar + 16]
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
    global B, F
    
    # Inspired by the code in the function
    # modp_splitting_data in algebras/quatalg/quaternion_algebra.py
    if p.number_field() != B.base():
        raise ValueError, "p must be a prime ideal in the base field of the quaternion algebra"
    if not p.is_prime():
        raise ValueError, "p must be prime"
    
    if F is not p.number_field():
        raise ValueError, "p must be a prime of %s"%F
    
    k = p.residue_field()

    from sage.all import MatrixSpace
    M = MatrixSpace(k, 2)
    i2, j2 = B.invariants()
    i2 = k(i2); j2 = k(j2)    
    if k.characteristic() == 2:
        if i2 == 0 or j2 == 0:
            raise NotImplementedError
        return M([0,1,1,0]), M([1,1,0,1])
    # Find I -- just write it down
    I = M([0,i2,1,0])
    # Find J -- I figured this out by just writing out the general case
    # and seeing what the parameters have to satisfy
    i2inv = 1/i2
    a = None
    for b in list(k):
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

def modp_splitting_map(p):
    """
    Return a map from subset of B to 2x2 matrix space isomorphic
    to R tensor OF/p.

    INPUT:
        - `B` -- quaternion algebra over F=Q(sqrt(5))
          with invariants -1, -1.
        - `p` -- prime ideal of F=Q(sqrt(5))

    EXAMPLES::

        sage: from psage.modform.hilbert.sqrt5.sqrt5 import F, B, modp_splitting_map
        sage: i,j,k = B.gens()
        sage: theta = modp_splitting_map(F.primes_above(5)[0])
        sage: theta(i + j - k)
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
        sage: theta = modp_splitting_map(F.primes_above(3)[0])
        sage: smod = theta(s); smod
        [2*abar + 1   abar + 1]
        [  abar + 1       abar]
        sage: smod^2 - 4*smod + 21
        [0 0]
        [0 0]    
    """
    I, J = modp_splitting(p)
    F = p.residue_field()
    def f(x):
        return F(x[0]) + I*F(x[1]) + J*F(x[2]) + I*J*F(x[3])
    return f

def icosian_gens():
    """
    Return generators of the icosian group, as elements of the 
    Hamilton quaternion algebra B over Q(sqrt(5)).
    
    AUTHOR: William Stein

    EXAMPLES::
    
        sage: from psage.modform.hilbert.sqrt5.sqrt5 import icosian_gens
        sage: icosian_gens()
        [i, j, k, -1/2 + 1/2*i + 1/2*j + 1/2*k, 1/2*i + 1/2*a*j + (-1/2*a + 1/2)*k]
        sage: [a.reduced_norm() for a in icosian_gens()]
        [1, 1, 1, 1, 1]
    """
    global B, F
    
    omega = F.gen()  # (1+sqrt(5))/2
    omega_bar = 1 - F.gen() # (1-sqrt(5))/2 = 1 - (1+sqrt(5))/2
    return [B(v)/2 for v in [(0,2,0,0), (0,0,2,0), 
                 (0,0,0,2), (-1,1,1,1), (0,1,omega,omega_bar)]]



def compute_all_icosians():
    """
    Return a list of the elements of the Icosian group of order 120,
    which we compute by generating enough products of icosian
    generators.
    
    EXAMPLES::

        sage: from psage.modform.hilbert.sqrt5.sqrt5 import compute_all_icosians, all_icosians
        sage: v = compute_all_icosians()
        sage: len(v)
        120
        sage: v
        [1/2 + 1/2*a*i + (-1/2*a + 1/2)*k, 1/2 + (-1/2*a + 1/2)*i + 1/2*a*j,..., -k, i, j, -i]
        sage: assert set(v) == set(all_icosians())  # double check
    """
    from sage.all import permutations, cartesian_product_iterator
    Icos = []
    ig = icosian_gens()
    per = permutations(range(5))
    exp = cartesian_product_iterator([range(1,i) for i in [5,5,5,4,5]])
    for f in exp:
        for p in per:
            e0 = ig[p[0]]**f[0]
            e1 = ig[p[1]]**f[1]
            e2 = ig[p[2]]**f[2]
            e3 = ig[p[3]]**f[3]
            e4 = ig[p[4]]**f[4]
            elt = e0*e1*e2*e3*e4
            if elt not in Icos:
                Icos.append(elt)
        if len(Icos) == 120:
            return Icos

@cached_function
def all_icosians():
    """
    Return a list of all 120 icosians, from a precomputed table.

    EXAMPLES::

    sage: from psage.modform.hilbert.sqrt5.sqrt5 import all_icosians
    sage: v = all_icosians()
    sage: len(v)
    120
    """
    s = '[1+a*i+(-a+1)*k,1+(-a+1)*i+a*j,-a+i+(a-1)*j,1+(-a)*i+(-a+1)*k,-a-j+(-a+1)*k,1+(a-1)*i+(-a)*j,-1+(-a)*i+(a-1)*k,-1+(a-1)*i+(-a)*j,-a+1+(-a)*j-k,-1+a*i+(-a+1)*k,-a+1+i+(-a)*k,-1+(-a+1)*i+(-a)*j,(-a+1)*i+j+(-a)*k,a-1+(-a)*j+k,(a-1)*i-j+a*k,a+i+(-a+1)*j,1+a*i+(a-1)*k,-1+(-a)*i+(-a+1)*k,a*i+(-a+1)*j-k,a-1+i+a*k,(-a)*i+(a-1)*j+k,a+j+(-a+1)*k,1+(-a+1)*i+(-a)*j,-1+(a-1)*i+a*j,a-i+(-a+1)*j,-1+a*i+(a-1)*k,a+j+(a-1)*k,-1+(-a+1)*i+a*j,(-a+1)*i+j+a*k,a-1+a*j+k,-a+i+(-a+1)*j,1+(-a)*i+(a-1)*k,a*i+(-a+1)*j+k,a-1-i+a*k,-a+j+(a-1)*k,1+(a-1)*i+a*j,(a-1)*i+j+a*k,-a+1+a*j+k,(-a)*i+(-a+1)*j+k,-a+1+i+a*k,-a-i+(-a+1)*j,-a+1+(-a)*j+k,(-a+1)*i-j+a*k,(a-1)*i+j+(-a)*k,a+i+(a-1)*j,a-1+a*j-k,-a+j+(-a+1)*k,-a+1-i+a*k,a*i+(a-1)*j+k,(-a)*i+(-a+1)*j-k,a-j+(a-1)*k,a-1+i+(-a)*k,-1+i+j+k,-1-i-j+k,1-i-j-k,1+i-j+k,a-j+(-a+1)*k,-1+i-j-k,1-i+j+k,(a-1)*i-j+(-a)*k,-a-j+(a-1)*k,1+i+j-k,a*i+(a-1)*j-k,-1-i+j-k,a-1+(-a)*j-k,-a+1+a*j-k,(-a)*i+(a-1)*j-k,-a+1-i+(-a)*k,-a-i+(a-1)*j,a-1-i+(-a)*k,(-a+1)*i-j+(-a)*k,a-i+(a-1)*j,-i+(-a)*j+(a-1)*k,i+a*j+(a-1)*k,i+a*j+(-a+1)*k,-i+a*j+(a-1)*k,-i+a*j+(-a+1)*k,i+(-a)*j+(a-1)*k,-i+(-a)*j+(-a+1)*k,i+(-a)*j+(-a+1)*k,-1-i-j-k,1+i+j+k,2,-a+1+a*i-j,-2,-a+(-a+1)*i-k,a-1+a*i+j,a+(-a+1)*i+k,a-1+(-a)*i+j,1+(-a+1)*j+(-a)*k,-a+1+a*i+j,-1+(-a+1)*j+a*k,a+(a-1)*i+k,-1+(a-1)*j+a*k,-a+(-a+1)*i+k,1+(-a+1)*j+a*k,-a+1+(-a)*i+j,-a+(a-1)*i+k,a-1+a*i-j,1+(a-1)*j+a*k,a+(-a+1)*i-k,-1+(-a+1)*j+(-a)*k,-a+1+(-a)*i-j,-a+(a-1)*i-k,a-1+(-a)*i-j,1+(a-1)*j+(-a)*k,a+(a-1)*i-k,-1+(a-1)*j+(-a)*k,-1+i+j-k,1-i+j-k,1-i-j+k,-1-i+j+k,-1+i-j+k,1+i-j-k,2*k,(-2)*j,(-2)*k,2*i,2*j,(-2)*i]'
    v = eval(s, {'a':F.gen(), 'i':B.gen(0), 'j':B.gen(1), 'k':B.gen(0)*B.gen(1)})
    return [B(x)/B(2) for x in v]


def icosian_ring_gens():
    """
    Return ring generators for the icosian ring (a maximal order) in the
    quaternion algebra ramified only at infinity over F=Q(sqrt(5)).
    These are generators over the ring of integers of F.

    OUTPUT:

       - list

    EXAMPLES::

        sage: from psage.modform.hilbert.sqrt5.sqrt5 import icosian_ring_gens
        sage: icosian_ring_gens()
        [1/2 + (1/2*a - 1/2)*i + 1/2*a*j, (1/2*a - 1/2)*i + 1/2*j + 1/2*a*k, 1/2*a*i + (1/2*a - 1/2)*j + 1/2*k, 1/2*i + 1/2*a*j + (1/2*a - 1/2)*k]
    """
    global B, F
    # See page 6 of Dembele.
    # DO NOT CHANGE THESE!  You'll break, e.g., the optimized code for
    # writing elements in terms of these...
    omega = F.gen()
    omega_bar = 1 - F.gen()
    return [B(v)/2 for v in [(1,-omega_bar,omega,0),
                             (0,-omega_bar,1,omega),
                             (0,omega,-omega_bar,1),
                             (0,1,omega,-omega_bar)]]

def icosian_ring_gens_over_ZZ():
    """
    Return basis over ZZ for the icosian ring, which has ZZ-rank 8.

    EXAMPLES::
    
        sage: from psage.modform.hilbert.sqrt5.sqrt5 import icosian_ring_gens_over_ZZ
        sage: v = icosian_ring_gens_over_ZZ(); v
        [1/2 + (1/2*a - 1/2)*i + 1/2*a*j, (1/2*a - 1/2)*i + 1/2*j + 1/2*a*k, 1/2*a*i + (1/2*a - 1/2)*j + 1/2*k, 1/2*i + 1/2*a*j + (1/2*a - 1/2)*k, 1/2*a + 1/2*i + (1/2*a + 1/2)*j, 1/2*i + 1/2*a*j + (1/2*a + 1/2)*k, (1/2*a + 1/2)*i + 1/2*j + 1/2*a*k, 1/2*a*i + (1/2*a + 1/2)*j + 1/2*k]
        sage: len(v)
        8
    """
    global B, F
    I = icosian_ring_gens()
    omega = F.gen()
    return I + [omega*x for x in I]

def tensor_over_QQ_with_RR(prec=53):
    """
    Return map from the quaternion algebra B to the tensor product of
    B over QQ with RealField(prec), viewed as an 8-dimensional real
    vector space.

    INPUT:
        - prec -- (integer: default 53); bits of real precision 

    OUTPUT:
        - a Python function

    EXAMPLES::

        sage: from psage.modform.hilbert.sqrt5.sqrt5 import tensor_over_QQ_with_RR, B
        sage: f = tensor_over_QQ_with_RR()
        sage: B.gens()
        [i, j, k]
        sage: f(B.0)
        (0.000000000000000, 1.00000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 1.00000000000000, 0.000000000000000, 0.000000000000000)
        sage: f(B.1)
        (0.000000000000000, 0.000000000000000, 1.00000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 1.00000000000000, 0.000000000000000)
        sage: f = tensor_over_QQ_with_RR(20)
        sage: f(B.0 - (1/9)*B.1)
        (0.00000, 1.0000, -0.11111, 0.00000, 0.00000, 1.0000, -0.11111, 0.00000)
    """
    global B, F
    from sage.all import RealField
    RR = RealField(prec=prec)
    V = RR**8
    S = F.embeddings(RR)
    def f(x):
        return V(sum([[sigma(a) for a in x] for sigma in S],[]))
    return f

def modp_icosians(p):
    """
    Return matrices of images of all 120 icosians modulo p.

    INPUT:
        - p -- *split* or ramified prime ideal of real quadratic field F=Q(sqrt(5))

    OUTPUT:
        - list of matrices modulo p.

    EXAMPLES::

        sage: from psage.modform.hilbert.sqrt5.sqrt5 import F, modp_icosians
        sage: len(modp_icosians(F.primes_above(5)[0]))
        120
        sage: v = modp_icosians(F.primes_above(11)[0])
        sage: len(v)
        120
        sage: v[0]
        [10  8]
        [ 8  1]
        sage: v[-1]
        [0 3]
        [7 0]
    """
    I, J = modp_splitting(p); K = I*J
    k = p.residue_field()
    G = [k(g[0]) + k(g[1])*I + k(g[2])*J + k(g[3])*K for g in icosian_gens()]
    from sage.all import MatrixGroup
    return [g.matrix() for g in MatrixGroup(G)]

class P1ModList(object):
    """
    Object the represents the elements of the projective line modulo
    a nonzero *prime* ideal of the ring of integers of Q(sqrt(5)).

    Elements of the projective line are represented by elements of a 2-dimension
    vector space over the residue field.

    EXAMPLES::

    We construct the projective line modulo the ideal (2), and illustrate
    all the standard operations with it::
    
        sage: from psage.modform.hilbert.sqrt5.sqrt5 import F, P1ModList
        sage: P1 = P1ModList(F.primes_above(2)[0]); P1
        Projective line over Residue field in abar of Fractional ideal (2)
        sage: len(P1)
        5
        sage: list(P1)
        [(0, 1), (1, 0), (1, abar), (1, abar + 1), (1, 1)]
        sage: P1.random_element()   # random output
        (1, abar + 1)
        sage: z = P1.random_element(); z   # random output
        (1, 0)
        sage: z[0].parent()
        Residue field in abar of Fractional ideal (2)
        sage: g = z[0].parent().gen()
        sage: P1.normalize((g,g))
        (1, 1)
        sage: P1((g,g))
        (1, 1)
    """
    def __init__(self, c):
        """
        INPUT:
           - c -- a nonzero prime of the ring of integers of Q(sqrt(5))

        EXAMPLES::

            sage: from psage.modform.hilbert.sqrt5.sqrt5 import F, P1ModList
            sage: P1ModList(F.primes_above(3)[0])
            Projective line over Residue field in abar of Fractional ideal (3)
            sage: P1ModList(F.primes_above(11)[1])
            Projective line over Residue field of Fractional ideal (3*a - 1)
            sage: list(P1ModList(F.primes_above(5)[0]))
            [(0, 1), (1, 0), (1, 1), (1, 2), (1, 3), (1, 4)]
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
        """
        Return a random element of this projective line.

            sage: from psage.modform.hilbert.sqrt5.sqrt5 import F, P1ModList
            sage: P1 = P1ModList(F.primes_above(13)[0]); P1
            Projective line over Residue field in abar of Fractional ideal (13)
            sage: P1.random_element()   # random output
            (1, 10*abar + 5)
        """
        import random
        return random.choice(self._list)
        
    def normalize(self, uv):
        """
        Normalize a representative element so it is either of the form
        (1,*) if the first entry is nonzero, or of the form (0,1)
        otherwise.

        EXAMPLES::

            sage: from psage.modform.hilbert.sqrt5.sqrt5 import F, P1ModList
            sage: p = F.primes_above(13)[0]
            sage: P1 = P1ModList(p)
            sage: k = p.residue_field()
            sage: g = k.gen()
            sage: P1.normalize([3,4])
            (1, 10)
            sage: P1.normalize([g,g])
            (1, 1)
            sage: P1.normalize([0,g])
            (0, 1)
        """
        w = self._V(uv)
        if w[0]:
            w = (~w[0]) * w
            w.set_immutable()
            #assert w in self._list
            return w
        else:
            return self._list[0]

    def __len__(self):
        """
        Return number of elements of this P1.

        EXAMPLES::

            sage: from psage.modform.hilbert.sqrt5.sqrt5 import F, P1ModList
            sage: len(P1ModList(F.primes_above(3)[0]))
            10
            sage: len(P1ModList(F.primes_above(5)[0]))
            6
            sage: len(P1ModList(F.primes_above(19)[0]))
            20
            sage: len(P1ModList(F.primes_above(19)[1]))
            20
        """
        return len(self._list)
    
    def __getitem__(self, i):
        """
        Return i-th element.

        EXAMPLES::

            sage: from psage.modform.hilbert.sqrt5.sqrt5 import F, P1ModList
            sage: P = P1ModList(F.primes_above(3)[0]); list(P)
            [(0, 1), (1, 0), (1, 2*abar), (1, abar + 1), (1, abar + 2), (1, 2), (1, abar), (1, 2*abar + 2), (1, 2*abar + 1), (1, 1)]
            sage: P[2]
            (1, 2*abar)
        """
        return self._list[i]

    def __call__(self, x):
        """
        Coerce x into this P1 list. Here x is anything that coerces to
        the 2-dimensional vector space over the residue field.  The
        result is normalized (in fact this function is just an alias
        for the normalize function).

        EXAMPLES::

            sage: from psage.modform.hilbert.sqrt5.sqrt5 import F, P1ModList
            sage: p = F.primes_above(13)[0]
            sage: k = p.residue_field(); g = k.gen()
            sage: P1 = P1ModList(p)
            sage: P1([3,4])
            (1, 10)
            sage: P1([g,g])
            (1, 1)
            sage: P1(P1([g,g]))
            (1, 1)
        """
        return self.normalize(x)

    def __repr__(self):
        """
        EXAMPLES::
        
            sage: from psage.modform.hilbert.sqrt5.sqrt5 import F, P1ModList
            sage: P1ModList(F.primes_above(19)[1]).__repr__()
            'Projective line over Residue field of Fractional ideal (-4*a + 3)'
        """
        return 'Projective line over %s'%self._F
    

def P1_orbits(p):
    """
    INPUT:
       - p -- a split or ramified prime of the integers of Q(sqrt(5)).

    OUTPUT:
        - ``orbits`` -- dictionary mapping elements of P1 to a choice of orbit rep
        - ``reps`` -- list of representatives for the orbits
        - ``P1`` -- the P1ModList object

    AUTHOR: William Stein

    EXAMPLES::

        sage: from psage.modform.hilbert.sqrt5.sqrt5 import P1_orbits, F
        sage: orbits, reps, P1 = P1_orbits(F.primes_above(5)[0])
        sage: orbits   # random output
        {(1, 2): (1, 0), (0, 1): (1, 0), (1, 3): (1, 0), (1, 4): (1, 0), (1, 0): (1, 0), (1, 1): (1, 0)}
        sage: reps   # random output
        [(1, 1)]
        sage: len(reps)
        1
        sage: P1
        Projective line over Residue field of Fractional ideal (2*a - 1)
        sage: orbits, reps, P1 = P1_orbits(F.primes_above(41)[0])
        sage: reps   # random output
        [(1, 40), (1, 5)]
        sage: len(reps)
        2
    """
    global B
    
    P1 = P1ModList(p)
    ICO = modp_icosians(p)
    
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

def P1_orbits2(p):
    """
    INPUT:
        - p -- a split or ramified prime of the integers of Q(sqrt(5)).

    OUTPUT:
        - list of disjoint sets of elements of P1 that are orbits

    AUTHOR: William Stein

    EXAMPLES::

        sage: from psage.modform.hilbert.sqrt5.sqrt5 import P1_orbits2, F
        sage: P1_orbits2(F.primes_above(5)[0])  # random output
        [set([(1, 2), (0, 1), (1, 3), (1, 4), (1, 0), (1, 1)])]
        sage: P1_orbits2(F.primes_above(11)[0])  # random output
        [set([(0, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 8), (1, 6), (1, 9), (1, 7), (1, 10), (1, 0), (1, 1)])]
        sage: len(P1_orbits2(F.primes_above(41)[0]))
        2
    """
    global B
    
    P1 = P1ModList(p)
    ICO = modp_icosians(p)
    
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

def totally_pos_gen(p):
    """
    Given a prime ideal p of a narrow class number 1 real quadratic
    field, return a totally positive generator for p.

    INPUT:
        - p -- prime ideal of narrow class number 1 real
          quadratic field

    OUTPUT:
        - generator of p that is totally positive

    AUTHOR: William Stein

    EXAMPLES::
    
        sage: from psage.modform.hilbert.sqrt5.sqrt5 import F, totally_pos_gen
        sage: g = totally_pos_gen(F.factor(19)[0][0]); g
        3*a + 4
        sage: g.norm()
        19
        sage: g.complex_embeddings()
        [2.14589803375032, 8.85410196624968]

        sage: for p in primes(14):
        ...       for P, e in F.factor(p):
        ...           g = totally_pos_gen(P)
        ...           print P, g, g.complex_embeddings()
        Fractional ideal (2) 2 [2.00000000000000, 2.00000000000000]
        Fractional ideal (3) 3 [3.00000000000000, 3.00000000000000]
        Fractional ideal (2*a - 1) a + 2 [1.38196601125011, 3.61803398874989]
        Fractional ideal (7) 7 [7.00000000000000, 7.00000000000000]
        Fractional ideal (3*a - 2) a + 3 [2.38196601125011, 4.61803398874989]
        Fractional ideal (3*a - 1) 2*a + 3 [1.76393202250021, 6.23606797749979]
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
    Return 8x8 Gram matrix of maximal order defined by R, which is assumed to be
    a basis for maximal order over ZZ.

    EXAMPLES::

        sage: from psage.modform.hilbert.sqrt5.sqrt5 import icosian_ring_gens_over_ZZ, gram_matrix_of_maximal_order
        sage: R = icosian_ring_gens_over_ZZ()
        sage: gram_matrix_of_maximal_order(R)
        [4 2 2 1 2 1 1 3]
        [2 4 1 1 1 2 3 3]
        [2 1 4 1 1 3 2 3]
        [1 1 1 4 3 3 3 2]
        [2 1 1 3 6 3 3 4]
        [1 2 3 3 3 6 4 4]
        [1 3 2 3 3 4 6 4]
        [3 3 3 2 4 4 4 6]
    """
    G = [[(R[i]*R[j].conjugate()).reduced_trace().trace()
          for i in range(8)] for j in range(8)]
    from sage.all import matrix
    return matrix(ZZ, G)

def qfminim(qf, N):
    """
    Call the PARI qfminim method on qf and 2*N, with smaller and
    and smaller search range, until a MemoryError is *not* raised.
    On a large-memory machine this will succeed the first time.

    EXAMPLES::

        sage: from psage.modform.hilbert.sqrt5.sqrt5 import icosian_ring_gens_over_ZZ, gram_matrix_of_maximal_order
        sage: R = icosian_ring_gens_over_ZZ()
        sage: G = gram_matrix_of_maximal_order(R)
        sage: qf = pari(G)
        sage: from psage.modform.hilbert.sqrt5.sqrt5 import icosian_ring_gens_over_ZZ, gram_matrix_of_maximal_order, qfminim
        sage: n, m, v = qfminim(qf, 2)
        sage: n
        120
        sage: m
        4
        sage: v[0]
        [0, 0, 0, -1, 1, 0, 0, 0]~
    """
    i = 32
    while i>10:
        try:
            return qf.qfminim(2*N, 2**i)
        except MemoryError:
            i -= 1
       

def bounded_elements(N):
    """
    Return elements in maximal order of B that have reduced norm
    whose trace is at most N.

    EXAMPLES::

        sage: from psage.modform.hilbert.sqrt5.sqrt5 import bounded_elements
        sage: X = bounded_elements(3)
        sage: len(X)
        180
        sage: rnX = [a.reduced_norm() for a in X]
        sage: set([a.trace() for a in rnX])
        set([2, 3])
        sage: set([a.norm() for a in rnX])
        set([1])
        sage: X = bounded_elements(5)
        sage: len(X)
        1200
        sage: rnX = [a.reduced_norm() for a in X]
        sage: set([a.trace() for a in rnX])
        set([2, 3, 4, 5])
        sage: set([a.norm() for a in rnX])
        set([1, 4, 5])    
    """
    # Get our maximal order
    R = icosian_ring_gens_over_ZZ()
    
    # Compute Gram matrix of R
    G = gram_matrix_of_maximal_order(R)

    # Make PARI quadratic form
    from sage.all import pari
    qf = pari(G)

    # Get the vectors of norm up to N.
    # The 2 is because we had to scale by 2 to get
    # rid of denominator in Gram matrix. 
    Z = qfminim(qf, N)
    Z2 = Z[2].sage().transpose()

    # For each vector, make the corresponding element of B.
    # TODO: This step massively dominates the runtime, and can be
    # easily made trivial with careful thought.
    V = []
    for i in range(Z2.nrows()):
        w = Z2[i]
        V.append(sum(w[j]*R[j] for j in range(8)))
    return V

from tables import primes_of_bounded_norm

def THETA(N):
    r"""
    Return representative elements of the maximal order of `R` of
    reduced norm `\pi_p` up to `N` modulo the left action of the units
    of `R` (the icosians).  Here `\pi_p` runs through totally positive
    generators of the odd prime ideals with norm up to `N`.

    INPUT:
        - `N` -- a positive integer >= 4.

    OUTPUT:
        - dictionary with keys the totally positive generators of the
          odd prime ideals with norm up to and including `N`, and values
          a dictionary with values the actual elements of reduced norm
          `\pi_p`.

    AUTHOR: William Stein

    EXAMPLES::

        sage: from psage.modform.hilbert.sqrt5.sqrt5 import THETA
        sage: d = THETA(9)
        pi = [2, a + 2, 3]
        Sorting through 2400 elements
        sage: d.keys()
        [a + 2, 3]
        sage: d[3]
        {(0, 1): a + 1/2*i + (-1/2*a + 1)*j + (-1/2*a + 1/2)*k, (1, 2): a - 1 + a*j, (1, abar): a - 1/2 + (1/2*a + 1/2)*i + (-1/2*a + 1)*j, (1, abar + 1): a - 1/2 + 1/2*i + 1/2*j + (-a + 1/2)*k, (1, abar + 2): a - 1 + (1/2*a + 1/2)*i + 1/2*j + (-1/2*a)*k, (1, 2*abar + 2): a - 1 + 1/2*a*i + (1/2*a + 1/2)*j + 1/2*k, (1, 2*abar): a - 1/2 + i + (-1/2*a + 1/2)*j + (-1/2*a)*k, (1, 2*abar + 1): a - 1/2 + (-1/2*a)*i + j + (-1/2*a + 1/2)*k, (1, 0): a + (-a + 1)*k, (1, 1): a - 1/2 + 1/2*a*i + j + (-1/2*a + 1/2)*k}
        sage: k = d[3].keys(); k
        [(0, 1), (1, 2), (1, abar), (1, abar + 1), (1, abar + 2), (1, 2*abar + 2), (1, 2*abar), (1, 2*abar + 1), (1, 0), (1, 1)]
        sage: d[3][k[0]]
        a + 1/2*i + (-1/2*a + 1)*j + (-1/2*a + 1/2)*k
        sage: d[3][k[0]].reduced_norm()
        3
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
    global B, F

    # Get primes of norm up to N.
    S = primes_of_bounded_norm(N)

    # Find totally positive generators pi_p
    pi = [totally_pos_gen(p) for p in S]
    print "pi =",pi

    # Compute traces of the generators, since that's what
    # the bounded_elements command computes up to.
    tr = [abs(x.trace()) for x in pi]
    N = max(tr)

    # A list that at least includes all elements (up to -1) whose
    # reduced norm has trace at most N.
    X = bounded_elements(N)
    
    # Compute mod-p local splitting maps
    theta_map = {}
    for i, p in enumerate(S):
        theta_map[pi[i]] = modp_splitting_map(p)

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

def hecke_ops(c, X):
    orbits, reps, P1 = P1_orbits(c)
    theta_c = modp_splitting_map(c)    
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
        from sage.all import matrix
        return matrix(ZZ, mat)
    ans = [(pi.norm(), pi, Tp(pi)) for pi in X.keys()]
    ans.sort()
    return ans


def hecke_ops2(c, X):
    reduce, reps, P1 = P1_orbits(c)
    theta_c = modp_splitting_map(c)    
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
        from sage.all import matrix
        return matrix(ZZ, mat)
    ans = [(pi.norm(), pi, Tp(pi)) for pi in X.keys()]
    ans.sort()
    return ans

class AlphaZ:
    def __init__(self, P):
        """
        Computing elements with norm pi, where P=(pi).

        INPUT:
            - P - a prime of O_F

        OUTPUT:
            - element alpha in B with norm pi
              whose image via the mod-p splitting map
              has column echelon form with first column z
        """
        self.R_Z = icosian_ring_gens_over_ZZ()
        self.P = P
        self.p = P.smallest_integer()
        self.pi = totally_pos_gen(P)
        self.deg = P.residue_class_degree()
        f = modp_splitting_map(self.P)
        self.n = [f(g) for g in self.R_Z]

    def local_map(self):
        n = self.n
        if self.deg == 1:
            k = n[0].parent().base_ring()
            W = k**4
            V = k**8
            return V.hom([W(x.list()) for x in n])
        else:
            # P is an inert prime
            from sage.all import GF
            k = GF(self.p)
            W = k**8
            V = k**8
            return V.hom([W(sum([y._vector_().list() for y in x.list()],[])) for x in n])

    def ideal_mod_p(self, z):
        """
        INPUT:
            - z - an element of P^1(O_F/p).
        """
        A = self.local_map()
        phi = self.local_map()
        V = phi.codomain()
        if self.deg == 1:
            g0 = V([z[0],0,z[1],0])
            g1 = V([0,z[0],0,z[1]])
            W = V.span([g0,g1])
        else:
            n = self.n
            M2 = n[0].parent()
            a = M2.base_ring().gen()
            g0 = M2([z[0],0,z[1],0])
            g1 = a*g0
            g2 = M2([0,z[0],0,z[1]])
            g3 = a*g2
            z = [g0,g1,g2,g3]
            W = V.span([V(sum([y._vector_().list() for y in x.list()],[])) for x in z])
        return phi.inverse_image(W)

    def ideal(self, z):
        Imod = self.ideal_mod_p(z)
        A = Imod.basis_matrix().lift()
        from sage.all import identity_matrix
        p = self.p
        B = A.stack(p*identity_matrix(ZZ,8))
        V = B.row_module()
        return V

    def ideal_basis(self, z):
        J = self.ideal(z)
        R = self.R_Z
        return [sum(g[i]*R[i] for i in range(8)) for g in J.gens()]

    def ideal_gram(self, W):
        G = [[(W[i]*W[j].conjugate()).reduced_trace().trace()
              for i in range(8)] for j in range(8)]
        from sage.all import matrix
        return matrix(ZZ, G)

    def alpha(self, z):
        """
        INPUT:
            - z - an element of P^1(O_F/P).
        """
        W = self.ideal_basis(z)
        G = self.ideal_gram(W)
        from sage.all import pari
        qf = pari(G)
        t = self.pi.trace()
        c = qfminim(qf, t)
        #print "number of vectors", c[0]
        for r in c[2].sage().transpose():
            a = sum([W[i]*r[i] for i in range(8)])
            if a.reduced_norm() == self.pi:
                return a
        raise ValueError, "bug"

    def all_alpha(self):
        return [self.alpha(z) for z in P1ModList(self.P)]

#@cached_function

import os
path = '/tmp/hmf-%s'%os.environ['USER']
if not os.path.exists(path):
    os.makedirs(path)
@disk_cached_function(path, memory_cache=True)
def hecke_elements(P):
    if P.norm() == 4:
        # hardcode this special case.
        return [~a for a in hecke_elements_2()]
    else:
        return [~a for a in AlphaZ(P).all_alpha()]


# Dumb code to get this special case.  The answer turns out to be:
# The following elements, divided by 2, where coordinates are in terms of 1,i,j,k, and a=(1+sqrt(5))/2:
#  [[-a-1,0,a-2,1], [-a-1,a-1,-a+1,a-1], [-a-1,-a+1,-a+1,a-1], [-a-1,-a+2,-1,0], [-a-1,a-2,-1,0]]
def hecke_elements_2():
    P = F.primes_above(2)[0]
    from sqrt5_fast import ModN_Reduction
    from sage.matrix.all import matrix
        
    f = ModN_Reduction(P)
    G = icosian_ring_gens()
    k = P.residue_field()
    g = k.gen()
    def pi(z):
        # Note that f(...) returns a string right now, since it's all done at C level.
        # This will prboably change, breaking this code. 
        M = matrix(k,2,eval(f(z).replace(';',','), {'g':g})).transpose()
        v = M.echelon_form()[0]
        v.set_immutable()
        return v
        
    # now just run through elements of icosian ring until we find enough...
    ans = {}
    a = F.gen()
    B = 1
    X = [i + a*j for i in range(-B,B+1) for j in range(-B,B+1)]
    from sage.misc.all import cartesian_product_iterator
    for v in cartesian_product_iterator([X]*4):
        z = sum(v[i]*G[i] for i in range(4))
        if z.reduced_norm() == 2:
            t = pi(z)
            if not ans.has_key(t):
                ans[t] = z
            if len(ans) == 5:
                return [x for _, x in ans.iteritems()]
    raise RuntimeError

class HMF:
    def __init__(self, N):
        from sage.all import QuadraticField, QuaternionAlgebra
        self._N = N
        red, reps, P1 = P1_orbits(N)
        self._reduce = red
        self._reps = reps
        self._P1 = P1
        self._theta_N = modp_splitting_map(N)

    def __repr__(self):
        return "Space of Hilbert modular forms over Q(sqrt(5)) of level %s (norm=%s) and dimension %s"%(
            self._N, self._N.norm(), self.dimension())

    def dimension(self):
        return len(self._reps)

    def hecke_matrix(self, P):
        mat = []
        alpha = AlphaZ(P)
        theta = self._theta_N
        Z = [theta(x)**(-1) for x in alpha.all_alpha()]
        P1 = self._P1
        red = self._reduce
        for x in self._reps:
            row = [0]*len(self._reps)
            for w in Z:
                y_red = red[P1(w*x)]
                row[self._reps.index(y_red)] += 1
            mat.append(row)
        from sage.all import matrix
        return matrix(ZZ, mat)

    Tp = hecke_matrix

    


########################################
# Generalizing to arbitrary level
########################################

def residue_ring(N):
    fac = N.factor()
    if len(fac) != 1:
        raise NotImplementedError, "P must be a prime power"
    from sage.rings.all import kronecker_symbol

    P, e = fac[0]
    p = P.smallest_integer()
    s = kronecker_symbol(p, 5)
    if e == 1:
        return P.residue_field()
    if s == 1:
        return ResidueRing_split(N, P, p, e)
    elif s == -1:
        return ResidueRing_inert(N, P, p, e)
    else:
        if e % 2 == 0:
            # easy case
            return ResidueRing_ramified_even(N, P, p, e)
        else:
            # hardest case
            return ResidueRing_ramified_odd(N, P, p, e)

class ResidueRing_base:
    def __call__(self, x):
        if x.parent() is not self._F:
            x = self._F(x)
        return self._to_ring(x)

    def list(self):
        return self._ring.list()
    
    def _to_ring(self, x):
        return self._ring(x[0]) + self._im_gen*self._ring(x[1])

    def ring(self):
        return self._ring

    def __repr__(self):
        return "Residue class ring of ZZ[(1+sqrt(5))/2] modulo %s of characteristic %s"%(
            self._N, self._p)
                                                                                         
    
class ResidueRing_split(ResidueRing_base):
    """
    Quotient of ring of integers of F by a prime power N.
    """
    def __init__(self, N, P, p, e):
        self._N = N
        self._F = P.number_field()
        self._p = p
        # Figure out the map to Z/p^e.
        fac = self._F.defining_polynomial().factor_padic(p, prec=e+1)
        assert len(fac) == 2
        roots = [(-a[0][0]).lift() for a in fac]
        gen = self._F.gen()
        if gen - roots[0] in N:
            im_gen = roots[0]
        elif gen - roots[1] in N:
            im_gen = roots[1]
        else:
            raise RuntimError, 'bug'
        self._ring = ZZ.quotient(p**e)
        self._im_gen = self._ring(im_gen)

    def lift(self, x):
        return self._F(x.lift())
    
class ResidueRing_inert(ResidueRing_base):
    def __init__(self, N, P, p, e):
        self._N = N
        self._F = P.number_field()
        self._p = p
        from sage.rings.all import Integers
        R = Integers(p**e)
        modulus = self._F.defining_polynomial().change_ring(R)
        S = R['x']
        self._base = R
        self._ring = S.quotient_by_principal_ideal(modulus)
        self._im_gen = self._ring.gen()

    def lift(self, x):
        f = x.lift().change_ring(ZZ)
        # f is a linear poly in generator of field
        return self._F(f)

    def list(self):
        R = self._ring
        x = R.gen()
        return [R(a + b*x) for a in self._base for b in self._base]
    

class ResidueRing_ramified_even(ResidueRing_base):
    def __init__(self, N, P, p, e):
        self._N = N
        self._F = P.number_field()
        self._p = p
        from sage.rings.all import Integers
        assert e%2 == 0
        R = Integers(p**(e//2))
        modulus = self._F.defining_polynomial().change_ring(R)
        S = R['x']
        self._ring = S.quotient_by_principal_ideal(modulus)
        self._im_gen = self._ring.gen()

    def lift(self, x):
        f = x.lift().change_ring(ZZ)
        return self._F(f)
        
class ResidueRing_ramified_odd(ResidueRing_base):
    """
    Residue class ring R = O_F / P^(2f-1), where e=2f-1
    is odd, and P=sqrt(5)O_F is the ramified prime. 

    Computing with this ring is trickier than all the rest,
    since it's not a quotient of Z[x] of the form Z[x]/(m,g),
    where m is an integer and g is a polynomial.
    The ring R is the quotient of
        O_F/P^(2f) = O_F/5^f = (Z/5^fZ)[x]/(x^2-x-1),
    by the ideal x^e.  We have
        O_F/P^(2f-2) subset R subset O_F/P^(2f)
    and each successive quotient has order 5 = #(O_F/P).
    Thus R has cardinality 5^(2f-1) and characteristic 5^f.
    The ring R can't be a quotient of Z[x] of the
    form Z[x]/(m,g), since such a quotient has
    cardinality m^deg(g) and characteristic m, and
    5^(2f-1) is not a power of 5^f.

    We thus view R as

        R = (Z/5^fZ)[x] / (x^2 - 5,  5^(f-1)*x).

    The elements of R are pairs (a,b) in (Z/5^fZ) x (Z/5^(f-1)Z),
    which correspond to the class of a + b*x.  The arithmetic laws
    are thus:

       (a,b) + (c,d) = (a+c mod 5^f, b+d mod 5^(f-1))

    and

       (a,b) * (c,d) = (a*c+b*d*5 mod 5^f, a*d+b*c mod 5^(f-1))

    The element omega = F.gen(), which is (1+sqrt(5))/2 maps to
    (1+x)/2 = (1/2, 1/2), which is the generator of this ring.

    EXAMPLES::

        sage: from psage.modform.hilbert.sqrt5.sqrt5 import F, residue_ring
        sage: P = F.primes_above(5)[0]; P
        Fractional ideal (2*a - 1)
        sage: R = residue_ring(P^5)
        sage: a = R(F.gen()); a
        [0 + 1*sqrt(5)]
        sage: a*a
        [5 + 0*sqrt(5)]
        sage: a*a*a
        [0 + 5*sqrt(5)]
        sage: a*a*a*a
        [0 + 0*sqrt(5)]
    """
    def __init__(self, N, P, p, e):
        self._N = N
        self._F = P.number_field()
        self._p = p
        from sage.rings.all import Integers
        f = e//2 + 1
        assert f*2 - 1 == e
        R0 = Integers(p**f)
        R1 = Integers(p**(f-1))
        self._ring = RamifiedProductRing(R0, R1)
        self._im_gen = self._ring.gen()
        self._sqrt5 = (self._F.gen()*2-1)**2

    def lift(self, x):
        return x[0].lift() + self._sqrt5 * x[1].lift()

class RamifiedProductRingElement:
    def __init__(self, parent, x, check=True):
        self._parent = parent
        if isinstance(x, RamifiedProductRingElement) and x._parent is parent:
            self._x = x
        else:
            if check:
                if isinstance(x, (list, tuple, RamifiedProductRingElement)):
                    self._x = (parent._R0(x[0]), parent._R1(x[1]))
                else:
                    self._x = (parent._R0(x), parent._R1(0))
            else:
                self._x = (x[0], x[1])

    def __getitem__(self, i):
        return self._x[i]
                
    def __repr__(self):
        return '[%s + %s*sqrt(5)]'%self._x

    def __add__(left, right):
        a, b = left._x
        c, d = right._x
        return RamifiedProductRingElement(left._parent, [a+b, c+d], check=False)

    def __sub__(left, right):
        a, b = left._x
        c, d = right._x
        return RamifiedProductRingElement(left._parent, [a-b, c-d], check=False)

    def __mul__(left, right):
        a, b = left._x
        c, d = right._x
        return RamifiedProductRingElement(left._parent, [a*c+b*d*5, a*d+b*c], check=False)

class RamifiedProductRing:
    def __init__(self, R0, R1):
        self._R0 = R0
        self._R1 = R1
        self._gen = self([ZZ(1)/2, ZZ(1)/2])

    def __call__(self, x):
        return RamifiedProductRingElement(self, x)

    def gen(self):
        return self._gen

        
class Mod_P_reduction_map:
    def __init__(self, P):
        FAC = P.factor()
        assert len(FAC) == 1
        self._p, self._e = FAC[0]
        self._I, self._J, self._residue_ring = self._compute_IJ(self._p, self._e)
        self._IJ = self._I*self._J

    def __repr__(self):
        return "(Partial) homomorphism from %s onto 2x2 matrices modulo %s^%s"%(
            B, self._p._repr_short(), self._e)

    def domain(self):
        return B

    def codomain(self):
        return self._I.parent()

    def __call__(self, x):
        R = self._residue_ring
        if x.parent() is not B:
            x = B(x)
        return R(x[0]) + self._I*R(x[1]) + self._J*R(x[2]) + self._IJ*R(x[3])

    def _compute_IJ(self, p, e):
        global B, F
        if p.number_field() != B.base():
            raise ValueError, "p must be a prime ideal in the base field of the quaternion algebra"
        if not p.is_prime():
            raise ValueError, "p must be prime"

        if p.number_field() != B.base():
            raise ValueError, "p must be a prime ideal in the base field of the quaternion algebra"
        if not p.is_prime():
            raise ValueError, "p must be prime"
        if F is not p.number_field():
            raise ValueError, "p must be a prime of %s"%F

        R = residue_ring(p**e)
        if isinstance(R, ResidueRing_base):
            k = R.ring()
        else:
            k = R

        from sage.all import MatrixSpace
        M = MatrixSpace(k, 2)

        i2, j2 = B.invariants()
        i2 = R(i2); j2 = R(j2)

        if k.characteristic() == 2:
            raise NotImplementedError

        # Find I -- just write it down
        I = M([0,i2,1,0])
        # Find J -- I figured this out by just writing out the general case
        # and seeing what the parameters have to satisfy
        i2inv = k(1)/i2
        a = None
        for b in R.list():
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
        return I, J, R
    
