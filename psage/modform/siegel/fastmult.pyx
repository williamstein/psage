"""
Low level functions for coefficients of Siegel modular forms
"""


# coeffs_dict* are dictionaries of the form (a,b,c) -> Integer, where
# (a,b,c) is GL(2,Z)-reduced, i.e. |b| <= a <= c, and b^2-4ac <= 0,
# and where the keys consist either of all reduced triples in a box of the form  
# (0,A)x(0,B)x(0,C) or else of all reduced triples (a,b,c) with 4ac-b^2 below
# a given bound (and c <=  if 4ac-b^2=0).

include 'gmp.pxi'

from sage.structure.element cimport Element
import operator
include 'sage/structure/coerce.pxi'
from sage.rings.integer_ring import ZZ
from sage.rings.integer cimport Integer

cdef struct int_triple:
    int a
    int b
    int c

cdef struct long_triple:
    long a
    long b
    long c

cdef struct int_septuple:
    int a
    int b
    int c
    int O
    int o
    int U
    int u


cdef struct int_quadruple:
    int a
    int b
    int c
    int s


cpdef  mult_coeff_int(int a, int b, int c, coeffs_dict1, coeffs_dict2):
    """
    Returns the value at (a,b,c) of the coefficient dict of the product
    of the two forms with dict coeffs_dict1 and coeffs_dict2.
    It is assumed that (a,b,c) is a key in coeffs_dict1 and coeffs_dict2.

    EXAMPLES::

        sage: from sage.modular.siegel.fastmult import mult_coeff_int
        sage: A, B, C, D = SiegelModularFormsAlgebra().gens()
        sage: A[(0,0,0)]
        1
        sage: mult_coeff_int(0, 0, 0, A.coeffs(), A.coeffs())
        1
        sage: mult_coeff_int(2, 2, 2, C.coeffs(), D.coeffs())
        1
        sage: mult_coeff_int(2, 1, 2, C.coeffs(), D.coeffs())
        8
    """
    cdef int a1, a2
    cdef int b1, b2
    cdef int c1, c2
    cdef int B1, B2

    cdef Integer result
    
    cdef mpz_t tmp
    cdef mpz_t left, right
    cdef mpz_t mult_coeff

    cdef mpz_t mpz_zero

    cdef long_triple tmp_triple

    result = PY_NEW(Integer)

    mpz_init(tmp)
    mpz_init(left)
    mpz_init(right)
    mpz_init(mult_coeff)
    mpz_init(mpz_zero)

    ## We assume that (a,b,c) is semipositive definite! (Havoc otherwise.)
    ## Reduce the form first so that our search space is smaller
    tmp_triple = _reduce_GL(a, b, c)

    sig_on()

    mpz_set_si( mpz_zero, 0)
    a = tmp_triple.a
    b = tmp_triple.b
    c = tmp_triple.c
    for a1 from 0 <= a1 < a+1:
        a2 = a - a1
        for c1 from 0 <= c1 < c+1:
            c2 = c - c1
            mpz_set_si(tmp, 4*a1*c1)
            mpz_sqrt(tmp, tmp)
            B1 = mpz_get_si(tmp)

            mpz_set_si(tmp, 4*a2*c2)
            mpz_sqrt(tmp, tmp)
            B2 = mpz_get_si(tmp)

            for b1 from max(-B1, b - B2) <= b1 < min(B1 + 1, b + B2 + 1) :
                ## Guarantees that both (a1,b1,c1) and (a2,b2,c2) are
                ## positive semidefinite                
                b2 = b - b1
                
                set_coeff(left, a1, b1, c1, coeffs_dict1)
                if 0 == mpz_cmp( left, mpz_zero): continue
                
                set_coeff(right, a2, b2, c2, coeffs_dict2)
                if  0 == mpz_cmp( right, mpz_zero): continue
                
                mpz_mul(tmp, left, right)
                mpz_add(mult_coeff, mult_coeff, tmp)
                
    sig_off()
    
    mpz_set(result.value, mult_coeff)
    mpz_clear(tmp)
    mpz_clear(left)
    mpz_clear(right)
    mpz_clear(mult_coeff)
    mpz_clear(mpz_zero)
    
    return result


cdef inline void set_coeff(mpz_t dest, int a, int b, int c, coeffs_dict):
    """
    Return the value of coeffs_dict at the triple obtained from
    reducing (a,b,c).
    It is assumed that the latter is a valid key
    and that (a,b,c) is positive semidefinite.

    """
    cdef long_triple tmp_triple

    mpz_set_si(dest, 0)

    tmp_triple = _reduce_GL(a, b, c)
    triple = (tmp_triple.a, tmp_triple.b, tmp_triple.c)
    try:
        mpz_set(dest, (<Integer>(coeffs_dict[triple])).value)
    except KeyError:
        pass
    return


cpdef mult_coeff_generic(int a, int b, int c, coeffs_dict1, coeffs_dict2, Ring R):
    """
    Returns the value at (a,b,c) of the coefficient dict of the product
    of the two forms with dict coeffs_dict1 and coeffs_dict2.
    It is assumed that (a,b,c) is a key in coeffs_dict1 and coeffs_dict2.

    EXAMPLES::

        sage: from sage.modular.siegel.fastmult import mult_coeff_generic
        sage: A, B, C, D = SiegelModularFormsAlgebra().gens()
        sage: AB = A.satoh_bracket(B)
        sage: CD = C.satoh_bracket(D)
        sage: R = CD.base_ring()
        sage: mult_coeff_generic(3, 2, 7, AB.coeffs(), CD.coeffs(), R)
        -7993474656/5*x^4 - 18273219840*x^3*y - 134223332976/5*x^2*y^2 - 120763645536/5*x*y^3 - 17770278912/5*y^4
    """
    cdef int a1, a2
    cdef int b1, b2
    cdef int c1, c2
    cdef int B1, B2
    cdef Integer nt = PY_NEW(Integer)    
    nmc = R(0)
    cdef long_triple tmp_triple

    cdef mpz_t tmp

    mpz_init(tmp)

    ## We assume that (a,b,c) is semipositive definite! (Havoc otherwise.)
    ## Reduce the form first so that our search space is smaller
    tmp_triple = _reduce_GL(a, b, c)

    sig_on()
    
    a = tmp_triple.a
    b = tmp_triple.b
    c = tmp_triple.c
    for a1 from 0 <= a1 < a+1:
        a2 = a - a1
        for c1 from 0 <= c1 < c+1:
            c2 = c - c1
            mpz_set_si(tmp, 4*a1*c1)
            mpz_sqrt(tmp, tmp)
            B1 = mpz_get_si(tmp)

            mpz_set_si(tmp, 4*a2*c2)
            mpz_sqrt(tmp, tmp)
            B2 = mpz_get_si(tmp)

            for b1 from max(-B1, b - B2) <= b1 < min(B1 + 1, b + B2 + 1) :
                ## Guarantees that both (a1,b1,c1) and (a2,b2,c2) are
                ## positive semidefinite                
                b2 = b - b1

                nl = get_coeff(a1, b1, c1, coeffs_dict1)
                if nl is None or nl.is_zero() : continue
                	
                nr = get_coeff(a2, b2, c2, coeffs_dict2)
                if nr is None or nr.is_zero() : continue
                
                nmc += nl*nr
                
    sig_off()

    mpz_clear(tmp)

    return nmc


cdef get_coeff(int a, int b, int c, coeffs_dict):
    """
    Return the value of coeffs_dict at the triple obtained from
    reducing (a,b,c).
    It is assumed that the latter is a valid key
    and that (a,b,c) is positive semidefinite. 


    """
    cdef long_triple tmp_triple

    tmp_triple = _reduce_GL(a, b, c)
    triple = (tmp_triple.a, tmp_triple.b, tmp_triple.c)
    try:
        return coeffs_dict[triple]
    except KeyError:
        return None


cpdef mult_coeff_generic_with_action(int a, int b, int c, coeffs_dict1, coeffs_dict2, Ring R):
    """
    Returns the value at (a,b,c) of the coefficient dict of the product
    of the two forms with dict coeffs_dict1 and coeffs_dict2.
    It is assumed that (a,b,c) is a key in coeffs_dict1 and coeffs_dict2.

    TO DO:

        I'm not sure what this does and it's not working

    """
    cdef int a1, a2
    cdef int b1, b2
    cdef int c1, c2
    cdef int B1, B2
    cdef Integer nt = PY_NEW(Integer)    
    nmc = R(0)
    cdef long_triple tmp_triple

    cdef mpz_t tmp

    mpz_init(tmp)

    ## We assume that (a,b,c) is semipositive definite! (Havoc otherwise.)
    ## Reduce the form first so that our search space is smaller
    tmp_triple = _reduce_GL(a, b, c)

    sig_on()
    
    a = tmp_triple.a
    b = tmp_triple.b
    c = tmp_triple.c
    for a1 from 0 <= a1 < a+1:
        a2 = a - a1
        for c1 from 0 <= c1 < c+1:
            c2 = c - c1
            mpz_set_si(tmp, 4*a1*c1)
            mpz_sqrt(tmp, tmp)
            B1 = mpz_get_si(tmp)

            mpz_set_si(tmp, 4*a2*c2)
            mpz_sqrt(tmp, tmp)
            B2 = mpz_get_si(tmp)

            for b1 from max(-B1, b - B2) <= b1 < min(B1 + 1, b + B2 + 1) :
                ## Guarantees that both (a1,b1,c1) and (a2,b2,c2) are
                ## positive semidefinite                
                b2 = b - b1

                nl = get_coeff_with_action(a1, b1, c1, coeffs_dict1, R)
                if nl is None or nl.is_zero() : continue
                	
                nr = get_coeff_with_action(a2, b2, c2, coeffs_dict2, R)
                if nr is None or nr.is_zero() : continue
                
                nmc += nl*nr
                
    sig_off()

    mpz_clear(tmp)

    return nmc



cdef get_coeff_with_action(int a0, int b0, int c0, coeffs_dict, R):
    """
    Return the value of coeffs_dict at the triple obtained from
    reducing (a,b,c).
    It is assumed that the latter is a valid key
    and that (a,b,c) is positive semidefinite. 
    """
    cdef int_septuple tmp_septuple
    cdef int a,b,c,d
    
    tmp_septuple = _xreduce_GL(a0, b0, c0)
    triple = (tmp_septuple.a, tmp_septuple.b, tmp_septuple.c)
    a = tmp_septuple.O
    b = tmp_septuple.o
    c = tmp_septuple.U
    d = tmp_septuple.u
    
    try:
        v = coeffs_dict[triple]
    except KeyError:
        return None

    from sage.algebras.all import GroupAlgebra
    if isinstance( R, GroupAlgebra):
        B = R.base_ring()
    else:
        B = R

    from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
    if is_PolynomialRing( B):
        x = B.gen(0)
        y = B.gen(1)
        xp = x*a + y*c
        yp = x*b + y*d

    if is_PolynomialRing( R):
        v = v( xp, yp) 
        return v

    res = R(0)
    for p,chi in v._fs:
        if is_PolynomialRing( B): p = p(xp,yp)
        G = B.group()
        det = G.gen(0)
        sig = G.gen(1)
        if chi == G.identity(): e = 1
        elif chi == det: e = a*d-b*c
        elif chi == sig: e = _sig(a,b,c,d)
        else: e  = (a*d-b*c) * _sig(a,b,c,d)
        res += p*e*R(chi)
    return res


cdef _sig( int a, int b, int c, int d):
    a = a%2
    b = b%2
    c = c%2
    d = d%2
    if 0 == b and 0 == c: return 1
    if 1 == a+d: return 1
    return -1
    

def reduce_GL(a, b, c):
    """
    Return the GL2(ZZ)-reduced form equivalent to (positive semidefinite)
    quadratic form ax^2+bxy+cy^2.

    EXAMPLES::

        sage: from sage.modular.siegel.fastmult import reduce_GL
        sage: reduce_GL(2,2,4)
        (2, 2, 4)
        sage: reduce_GL(3,5,3)
        (1, 1, 3)
        sage: reduce_GL(1,2,1)
        (0, 0, 1)
    """
    cdef long_triple res

    # We want to check that (a,b,c) is semipositive definite
    # since otherwise we might end up in an infinite loop.
    # TODO: the discriminant can become to big
    if b*b-4*a*c > 0 or a < 0 or c < 0:
        raise NotImplementedError, "only implemented for nonpositive discriminant"
            
    res = _reduce_GL(a, b, c)
    return (res.a, res.b, res.c)


cdef long_triple _reduce_GL(long a, long b, long c):
    """
    Return the GL2(ZZ)-reduced form equivalent to (positive semidefinite)
    quadratic form ax^2+bxy+cy^2.
    (Following algorithm 5.4.2 found in Cohen's Computational Algebraic Number Theory.)
    """
    cdef long_triple res
    cdef long twoa, q, r
    cdef long tmp

    if b < 0:
        b = -b
        
    while not (b<=a and a<=c):  ## while not GL-reduced
        twoa = 2 * a
        #if b not in range(-a+1,a+1):
        if b > a:
            ## apply Euclidean step (translation)
            q = b / twoa #(2*a) 
            r = b % twoa  #(2*a)
            if r > a:
                ## want (quotient and) remainder such that -a<r<=a
                r = r - twoa  # 2*a;
                q = q + 1
            c = c-b*q+a*q*q
            b = r
            
        ## apply symmetric step (reflection) if necessary
        if a > c:
            #(a,c) = (c,a)
            tmp = c
            c = a
            a = tmp

        #b=abs(b)
        if b < 0:
            b = -b
        ## iterate
    ## We're finished! Return the GL2(ZZ)-reduced form (a,b,c):

    res.a = a
    res.b = b
    res.c = c
            
    return res


def xreduce_GL(a, b, c):
    """
    Return the GL2(ZZ)-reduced form equivalent to (positive semidefinite)
    quadratic form ax^2+bxy+cy^2.

    EXAMPLES::

        sage: from sage.modular.siegel.fastmult import xreduce_GL
        sage: xreduce_GL(3,1,1)
        ((1, 1, 3), (0, 1, 1, 0))
        sage: xreduce_GL(1,2,1)
        ((0, 0, 1), (0, 1, 1, 1))
        sage: xreduce_GL(3,8,2)
        Traceback (most recent call last):
        NotImplementedError: only implemented for nonpositive discriminant
        sage: xreduce_GL(1,1,1)
        ((1, 1, 1), (1, 0, 0, 1))
        sage: xreduce_GL(3,2,9)
        ((3, 2, 9), (1, 0, 0, 1))
        sage: xreduce_GL(3,7,9)
        ((3, 1, 5), (1, 0, 1, 1))


    """
    cdef int_septuple res

    # We want to check that (a,b,c) is semipositive definite
    # since otherwise we might end up in an infinite loop.
    if b*b-4*a*c > 0 or a < 0 or c < 0:
        raise NotImplementedError, "only implemented for nonpositive discriminant"
            
    res = _xreduce_GL(a, b, c)
    return ((res.a, res.b, res.c), (res.O, res.o, res.U, res.u))


cdef int_septuple _xreduce_GL(int a, int b, int c):
    """
    Return the GL2(ZZ)-reduced form f and a matrix M such that
    ax^2+bxy+cy^2 = f((x,y)*M).

    TODO
         _xreduce_GL is a factor 20 slower than xreduce_GL.
         How can we improve this factor ?

    """
    cdef int_septuple res
    cdef int twoa, q, r
    cdef int tmp
    cdef int O,o,U,u
    
    # If [A,B,C] is the form to be reduced, then after each reduction
    # step we have [A,B,C] =[a,b,c]( (x,y)*matrix(2,2,[O,o, U,u]) )
    O = u = 1
    o = U = 0

    if b < 0:
        b = -b
        O = -1
        
    while not (b<=a and a<=c):  ## while not GL-reduced
        twoa = 2 * a
        #if b not in range(-a+1,a+1):
        if b > a:
            ## apply Euclidean step (translation)
            q = b / twoa #(2*a) 
            r = b % twoa  #(2*a)
            if r > a:
                ## want (quotient and) remainder such that -a<r<=a
                r = r - twoa  # 2*a;
                q = q + 1
            c = c-b*q+a*q*q
            b = r
            O += q*o
            U += q*u
            
        ## apply symmetric step (reflection) if necessary
        if a > c:
            #(a,c) = (c,a)
            tmp = c; c = a; a = tmp
            tmp = O; O = o; o = tmp
            tmp = U; U = u; u = tmp
            
        #b=abs(b)
        if b < 0:
            b = -b
            O = -O
            U = -U

        ## iterate
    ## We're finished! Return the GL2(ZZ)-reduced form (a,b,c) and the O,o,U,u-matrix:

    res.a = a
    res.b = b
    res.c = c
    res.O = O
    res.o = o
    res.U = U
    res.u = u
            
    return res





cdef inline int poly(int e, int f, int g, int h,
                     int i, int j, int k, int l,
                     int m, int n, int o, int p) :
## reverse the linebreaks; this won't work
    return (  4*(f*(k*p - l*o) - g*(j*p - l*n) + h*(j*o - k*n))
              - 6*(e*(k*p - l*o) - g*(i*p - l*m) + h*(i*o - k*m))
              + 10*(e*(j*p - l*n) - f*(i*p - l*m) + h*(i*n - j*m))
              - 12*(e*(j*o - k*n) - f*(i*o - k*m) + g*(i*n - j*m)) )


def chi35(disc, A,B,C,D) :
    r"""
    Returns the odd weight generator the ring of Siegel modular forms of degree 2
    as in Igusa.  Uses a formula from Ibukiyama [I].

    EXAMPLES::

        sage: A, B, C, D = SiegelModularFormsAlgebra().gens()
        sage: from sage.modular.siegel.fastmult import chi35
        sage: chi = chi35(51,A,B,C,D)
        sage: chi[(0,0,0)]
        Traceback (most recent call last):
        KeyError: (0, 0, 0)
        sage: chi[(2,1,3)]
        1
        sage: chi[(3,1,4)]
        -129421

    """
    cdef dict Acoeffs, Bcoeffs, Ccoeffs, Dcoeffs
    cdef int a1, a2, a3, a4
    cdef int c1, c2, c3, c4
    cdef int b1, b2, b3, b4
    cdef int a, b, c
    cdef int B1, B1r, B2, B2r, B3, B4

    cdef mpz_t coeff, acc, tmp

    cdef mpz_t mpz_zero

    mpz_init(coeff)
    mpz_init(acc)
    mpz_init(tmp)
    mpz_init(mpz_zero)


    sig_on()
    
    mpz_set_si( mpz_zero, 0)
    coeffs = PY_NEW( dict )

    prec = min([ disc, A.prec(), B.prec(), C.prec(), D.prec()])

    if prec != disc:
        raise ValueError, "Insufficient precision of Igusa forms"
    from siegel_modular_form_prec import SiegelModularFormPrecision
    prec_class = SiegelModularFormPrecision( prec)

    Acoeffs = A.coeffs()
    Bcoeffs = B.coeffs()
    Ccoeffs = C.coeffs()
    Dcoeffs = D.coeffs()

    for trip in prec_class.positive_forms():
        (a, b, c) = trip
        mpz_set_si(coeff, 0)
        for c1 from 0 <= c1 < c-1 :
            for c2 from 0 <= c2 < c-1-c1 :
                for c3 from 1 <= c3 < c-c1-c2 :
                    c1r = c - c1
                    c2r = c1r - c2
                    c4 = c2r - c3
                    for a1 from 0 <= a1 < a-1 :
                        for a2 from 0<= a2 < a-1-a1 :
                            for a3 from 1 <= a3 < a-a1-a2 :
                                a1r = a - a1
                                a2r = a1r - a2
                                a4 = a2r - a3

                                mpz_set_si(tmp, 4*a1*c1)
                                mpz_sqrt(tmp, tmp)
                                B1 = mpz_get_si(tmp)

                                mpz_set_si(tmp, 4*a1r*c1r)
                                mpz_sqrt(tmp, tmp)
                                B1r = mpz_get_si(tmp)

                                mpz_set_si(tmp, 4*a2*c2)
                                mpz_sqrt(tmp, tmp)
                                B2 = mpz_get_si(tmp)

                                mpz_set_si(tmp, 4*a2r*c2r)
                                mpz_sqrt(tmp, tmp)
                                B2r = mpz_get_si(tmp)

                                mpz_set_si(tmp, 4*a3*c3)
                                mpz_sqrt(tmp, tmp)
                                B3 = mpz_get_si(tmp)

                                mpz_set_si(tmp, 4*a4*c4)
                                mpz_sqrt(tmp, tmp)
                                B4 = mpz_get_si(tmp)

                                for b1 from max(-B1, b - B1r) <= b1 < min(B1 + 1, b + B1r + 1) :
                                    for b2 from max(-B2, b - b1 - B2r) <= b2 < min(B2 + 1, b - b1 + B2r + 1) :
                                        for b3 from max(-B3, b - b1 - b2 - B4) <= b3 < min(B3 + 1, b - b1 - b2 + B4 + 1) :
                                            b4 = b-(b1+b2+b3)

                                            set_coeff(acc, a1, b1, c1, Acoeffs)
                                            set_coeff(tmp, a2, b2, c2, Bcoeffs)
                                            mpz_mul(acc, acc, tmp)
                                            set_coeff(tmp, a3, b3, c3, Ccoeffs)
                                            mpz_mul(acc, acc, tmp)
                                            set_coeff(tmp, a4, b4, c4, Dcoeffs)
                                            mpz_mul(acc, acc, tmp)

                                            mpz_set_si(tmp, poly(a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4))
                                            mpz_mul(acc, acc, tmp)

                                            mpz_add(coeff, coeff, acc)
        if mpz_cmp( coeff, mpz_zero) != 0:
            coeffs[trip] = PY_NEW( Integer)
            mpz_set((<Integer>coeffs[trip]).value, coeff)

    for t in coeffs:
        coeffs[t] = Integer(coeffs[t]/41472)

    sig_off()
    
    mpz_clear(coeff)
    mpz_clear(acc)
    mpz_clear(tmp)
    mpz_clear(mpz_zero)
    
    return coeffs

