import sqrt5, sqrt5_fast
from sage.misc.all import cputime
from sage.rings.all import Integer

F = sqrt5.F

def ideals_of_bounded_norm(B):
    return sum([v for n, v in F.ideals_of_bdd_norm(B).iteritems() if n != 1], [])
    
def ideals_of_norm(v):
    try:
        v = list(v)
    except TypeError:
        v = [Integer(v)]
    z = F.ideals_of_bdd_norm(max(v))
    return sum([z[n] for n in v if n>1],[])

def canonical_gen(I):
    return I.gens_reduced()[0]

def no_space(s):
    return str(s).replace(' ', '')

def dimensions(v, filename=None):
    """
    Compute dimensions of spaces of Hilbert modular forms for all the levels in v.
    The format is:

        Norm   dimension  generator  time
    """
    F = open(filename,'a') if filename else None
    for N in ideals_of_norm(v):
        t = cputime()
        H = sqrt5_fast.IcosiansModP1ModN(N)
        tm = cputime(t)
        s = '%s %s %s %s'%(N.norm(), H.cardinality(), no_space(canonical_gen(N)), tm)
        print s
        if F:
            F.write(s+'\n')
            F.flush()

def charpolys(v, B, filename=None):
    """
    Compute characteristic polynomials of T_P for primes P with norm <= B
    for spaces of Hilbert modular forms for all the levels in v.
    """
    F = open(filename,'a') if filename else None
    P = [p for p in ideals_of_bounded_norm(B) if p.is_prime()]
    for N in ideals_of_norm(v):
        t = cputime()
        H = sqrt5_fast.IcosiansModP1ModN(N)
        T = [(p.smallest_integer(),H.hecke_matrix(p).fcp()) for p in P if
             gcd(Integer(p.norm()), Integer(N.norm())) == 1]
        tm = cputime(t)
        s = '%s %s %s %s'%(N.norm(), no_space(canonical_gen(N)), tm, no_space(T))
        print s
        if F:
            F.write(s+'\n')
            F.flush()

def one_charpoly(v, filename=None):
    """
    Compute and factor one characteristic polynomials for all the
    levels in v.  Always compute the charpoly of T_P where P is the
    smallest prime not dividing the level.
    """
    F = open(filename,'a') if filename else None
    P = [p for p in ideals_of_bounded_norm(100) if p.is_prime()]
    for N in ideals_of_norm(v):
        NN = Integer(N.norm())
        t = cputime()
        H = sqrt5_fast.IcosiansModP1ModN(N)
        t0 = cputime(t)
        for p in P:
            if Integer(p.norm()).gcd(NN) == 1:
                break
        t = cputime()            
        T = H.hecke_matrix(p)
        t1 = cputime(t)
        t = cputime()
        f = T.fcp()
        t2 = cputime(t)
        s = '%s\t%s\t%s\t%s\t%s\t(%.1f,%.1f,%.1f)'%(N.norm(), no_space(canonical_gen(N)), 
                                                    p.smallest_integer(), no_space(canonical_gen(p)), no_space(f),
                                                    t0, t1, t2,)
        print s
        if F:
            F.write(s+'\n')
            F.flush()
