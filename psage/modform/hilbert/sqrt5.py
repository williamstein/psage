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
    if F.characteristic() == 2:
        raise NotImplementedError
    from sage.all import MatrixSpace
    M = MatrixSpace(F, 2)
    i2, j2 = B.invariants()
    i2 = F(i2); j2 = F(j2)    
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
    sqrt5 = F(5).sqrt()
    sigma = (1-sqrt5)/2
    tau = (1+sqrt5)/2
    return [B(v)/2 for v in [(0,2,0,0), (0,0,2,0), 
                 (0,0,0,2), (-1,1,1,1), (0,1,sigma,tau)]]


def modp_icosians(B, p):
    """
    Return matrices of images of all 168 icosians mod p.
    """
    I, J = modp_splitting(B, p); K = I*J
    k = p.residue_field()
    G = [k(g[0]) + k(g[1])*I + k(g[2])*J + k(g[3])*K for g in icosian_gens(B)]
    from sage.all import MatrixGroup
    return [g.matrix() for g in MatrixGroup(G)]

from sage.all import cached_function
@cached_function
def P1_orbits(B, p):
    """
    INPUT:

       - B -- quaternion algebra
       - p -- a prime of O_F, where F is the totally real field Q(sqrt(5)).

    AUTHOR: William Stein

    EXAMPLES::

    """
    from sage.all import P1NFList
    
    P1 = P1NFList(p)
    ICO = modp_icosians(B, p)
    
    def act(g, t):
        return P1.normalize(tuple(t[0]*g[0] + t[1]*g[1])).tuple()

    # This implementation is *horrendously* idiotically slow and
    # stupid.  But it works.
        
    orbits = {}
    cur = P1[0].tuple()
    while len(orbits) < len(P1):
        for g in ICO:
            s = act(g, cur)
            if not orbits.has_key(s):
                orbits[s] = cur
        if len(orbits) < len(P1):
            # choose element of P1 not a key
            for c in P1:
                if c.tuple() not in orbits:
                    cur = c.tuple()
                    break
    # done
    return orbits, act

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

def Theta(B, p):
    r"""
    Return representative elements of the maximal order of `R` of norm
    `\pi_p` up to the units of `R` (the icosians).  Here `\pi_p` is a
    totally positive generator of `p`.

    INPUT:
       - `B` -- quaternion algebra
       - `p` -- a prime of `\mathcal{O}_F`, where `F` is the totally
         real field `\QQ(\sqrt{5})`.

    AUTHOR: William Stein

    EXAMPLES::
    """
    raise NotImplementedError

def hecke_operator(B, P1, Theta_p):
    raise NotImplementedError
