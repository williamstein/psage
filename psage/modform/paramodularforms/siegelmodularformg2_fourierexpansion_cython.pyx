r"""
Functions for reduction of fourier indices and for multiplication of
Siegel modular forms.

AUTHORS:

- Craig Citro, Nathan Ryan Initian version
- Martin Raum (2009 - 04) Refined multiplication by Martin Raum
- Martin Raum (2009 - 07 - 28) Port to new framework
"""

#===============================================================================
# 
# Copyright (C) 2009 Martin Raum
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

## TODO : Clean this file up.

include "interrupt.pxi"
include "stdsage.pxi"
include "cdefs.pxi"

#include 'sage/ext/gmp.pxi'

from sage.rings.integer cimport Integer
from sage.rings.ring cimport Ring

cdef struct int_triple:
    int a
    int b
    int c

cdef struct int_quinruple:
    int a
    int b
    int c
    int d # determinant of reducing g \in GL2
    int s # sign of reducing g \in GL2
    
cdef struct int_septuple:
    int a
    int b
    int c
    int O # upper left entry of reducing g \in GL2
    int o # upper right entry of reducing g \in GL2
    int U # lower left entry of reducing g \in GL2
    int u # lower right entry of reducing g \in GL2


#####################################################################
#####################################################################
#####################################################################

#===============================================================================
# reduce_GL
#===============================================================================

def reduce_GL(tripel) :
    r"""
    Return the `GL_2(\ZZ)`-reduced form equivalent to (positive semidefinite)
    quadratic form `a x^2 + b x y + c y^2`.
    """
    cdef int a, b, c
    cdef int_triple res

    # We want to check that (a,b,c) is semipositive definite
    # since otherwise we might end up in an infinite loop.
    # TODO: the discriminant can become to big
    if b*b-4*a*c > 0 or a < 0 or c < 0:
        raise NotImplementedError, "only implemented for nonpositive discriminant"

    (a, b, c) = tripel            
    res = _reduce_GL(a, b, c)

    return ((res.a, res.b, res.c), 0)

#===============================================================================
# _reduce_GL
#===============================================================================

cdef int_triple _reduce_GL(int a, int b, int c) :
    r"""
    Return the `GL_2(\ZZ)`-reduced form equivalent to (positive semidefinite)
    quadratic form `a x^2 + b x y + c y^2`.
    (Following Algorithm 5.4.2 found in Cohen's Computational Algebraic Number Theory.)
    """
    cdef int_triple res
    cdef int twoa, q, r
    cdef int tmp

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

##TODO: This is a mess. In former implementations the sign was the determinant
##      of the reducing matrix. This has to be replaced by d and the sign is the
##      GL2(F_2) \cong S_3 character

#===============================================================================
# #==============================================================================
# # sreduce_GL
# #==============================================================================
# 
# def sreduce_GL(a, b, c):
#    """
#    Return the GL2(ZZ)-reduced form f and a determinant d and a sign s such
#    that d is the determinant of the matrices M in GL2(ZZ) such that
#    ax^2+bxy+cy^2 = f((x,y)*M).
#    """
#    cdef int_quadruple res
# 
#    # We want to check that (a,b,c) is semipositive definite
#    # since otherwise we might end up in an infinite loop.
#    if b*b-4*a*c > 0 or a < 0 or c < 0:
#        raise NotImplementedError, "only implemented for nonpositive discriminant"
#            
#    res = _sreduce_GL(a, b, c)
#    return ((res.a, res.b, res.c), (res.d, res.s))
# 
# #===============================================================================
# # int_quintuple _sreduce_GL
# #===============================================================================
# 
# cdef int_quintuple _sreduce_GL(int a, int b, int c):
#    """
#    Return the GL2(ZZ)-reduced form equivalent to (positive semidefinite)
#    quadratic form ax^2+bxy+cy^2.
#    (Following algorithm 5.4.2 found in Cohen's Computational Algebraic Number Theory.)
# 
#    TODO
#         _xreduce_GL is a factor 20 slower than xreduce_GL.
#         How can we improve this factor ?
# 
#    """
#    cdef int_quadruple res
#    cdef int twoa, q, r
#    cdef int tmp
#    cdef int s
#    
#    # If [A,B,C] is the form to be reduced, then after each reduction
#    # step we have s = det(M), where M is the GL2(ZZ)-matrix
#    # such that [A,B,C] =[a,b,c]( (x,y)*M)
#    s = 1
#    if b < 0:
#        b = -b
#        s = -1
#        
#    while not (b<=a and a<=c):  ## while not GL-reduced
#        twoa = 2 * a
#        #if b not in range(-a+1,a+1):
#        if b > a:
#            ## apply Euclidean step (translation)
#            q = b / twoa #(2*a) 
#            r = b % twoa  #(2*a)
#            if r > a:
#                ## want (quotient and) remainder such that -a<r<=a
#                r = r - twoa  # 2*a;
#                q = q + 1
#            c = c-b*q+a*q*q
#            b = r
#            
#        ## apply symmetric step (reflection) if necessary
#        if a > c:
#            #(a,c) = (c,a)
#            tmp = c; c = a; a = tmp
#            s *= -1
#            
#        #b=abs(b)
#        if b < 0:
#            b = -b
#            s *= -1
# 
#        ## iterate
#    ## We're finished! Return the GL2(ZZ)-reduced form (a,b,c) and the O,o,U,u-matrix:
# 
#    res.a = a
#    res.b = b
#    res.c = c
#    res.s = s
#            
#    return res
#===============================================================================

#===============================================================================
# xreduce_GL
#===============================================================================

def xreduce_GL(tripel) :
    r"""
    Return the `GL_2(\ZZ)`-reduced form equivalent to (positive semidefinite)
    quadratic form `a x^2 + b x y + c y^2`.
    """
    cdef int a, b, c
    cdef int_septuple res

    # We want to check that (a,b,c) is semipositive definite
    # since otherwise we might end up in an infinite loop.
#    if b*b-4*a*c > 0 or a < 0 or c < 0:
#        raise NotImplementedError, "only implemented for nonpositive discriminant"

    (a, b, c) = tripel            
    res = _xreduce_GL(a, b, c)
    return ((res.a, res.b, res.c), (res.O, res.o, res.U, res.u))

#===============================================================================
# _xreduce_GL
#===============================================================================

cdef int_septuple _xreduce_GL(int a, int b, int c) :
    r"""
    Return the `GL_2(\ZZ)`-reduced form `f` and a matrix `M` such that
    `a x^2 + b x y + c y^2 = f( (x,y)*M )`.

    TODO:
    
         ``_xreduce_GL`` is a factor 20 slower than ``xreduce_GL``.
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

#####################################################################
#####################################################################
#####################################################################

#===============================================================================
# mult_coeff_int_without_character
#===============================================================================

cpdef mult_coeff_int_without_character(result_key,
       coeffs_dict1, coeffs_dict2, ch1, ch2, result) :
    r"""
    Returns the value at `(a, b, c)` of the coefficient dictionary of the product
    of the two forms with dictionaries ``coeffs_dict1`` and ``coeffs_dict2``.
    It is assumed that `(a, b, c)` is a key in ``coeffs_dict1`` and ``coeffs_dict2``.
    """
    cdef int a, b, c
    cdef int a1, a2
    cdef int b1, b2
    cdef int c1, c2
    cdef int B1, B2
    
    cdef mpz_t tmp, mpz_zero
    cdef mpz_t left, right
    cdef mpz_t acc

    mpz_init(tmp)
    mpz_init(mpz_zero)
    mpz_init(left)
    mpz_init(right)
    mpz_init(acc)

    (a, b, c) = result_key

    sig_on()
    for a1 from 0 <= a1 < a+1 :
        a2 = a - a1
        for c1 from 0 <= c1 < c+1 :
            c2 = c - c1
            mpz_set_si(tmp, 4*a1*c1)
            mpz_sqrt(tmp, tmp)
            B1 = mpz_get_si(tmp)

            mpz_set_si(tmp, 4*a2*c2)
            mpz_sqrt(tmp, tmp)
            B2 = mpz_get_si(tmp)

            for b1 from max(-B1, b - B2) <= b1 < min(B1 + 1, b + B2 + 1) :
                ## Guarantes that both (a1,b1,c1) and (a2,b2,c2) are
                ## positive semidefinite                
                b2 = b - b1
                
                get_coeff_int(left, a1, b1, c1, coeffs_dict1)
                if mpz_cmp(left, mpz_zero) == 0 : continue
                
                get_coeff_int(right, a2, b2, c2, coeffs_dict2)
                if mpz_cmp(right, mpz_zero) == 0 : continue
                
                mpz_mul(tmp, left, right)
                mpz_add(acc, acc, tmp)
    sig_off()
    
    mpz_set((<Integer>result).value, acc)
    
    mpz_clear(tmp)
    mpz_clear(mpz_zero)
    mpz_clear(left)
    mpz_clear(right)
    mpz_clear(acc)
    
    return result

#===============================================================================
# mult_coeff_generic_without_character
#===============================================================================

cpdef mult_coeff_generic_without_character(result_key,
       coeffs_dict1, coeffs_dict2, ch1, ch2, result) :
    r"""
    Returns the value at `(a, b, c)`of the coefficient dictionary of the product
    of the two forms with dictionaries ``coeffs_dict1`` and ``coeffs_dict2``.
    It is assumed that `(a, b, c)` is a key in ``coeffs_dict1`` and ``coeffs_dict2``.
    """
    cdef int a, b, c
    cdef int a1, a2
    cdef int b1, b2
    cdef int c1, c2
    cdef int B1, B2

    cdef mpz_t tmp

    mpz_init(tmp)

    (a, b, c) = result_key

    sig_on()
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
                ## Guarantes that both (a1,b1,c1) and (a2,b2,c2) are
                ## positive semidefinite                
                b2 = b - b1

                left = get_coeff_generic(a1, b1, c1, coeffs_dict1)
                if left is None : continue
                
                right = get_coeff_generic(a2, b2, c2, coeffs_dict2)
                if right is None : continue
                
                result += left*right
    sig_off()

    mpz_clear(tmp)

    return result

#####################################################################

#===============================================================================
# get_coeff_int
#===============================================================================

cdef inline void get_coeff_int(mpz_t dest, int a, int b, int c, coeffs_dict):
    r"""
    Return the value of ``coeffs_dict`` at the triple obtained from
    reducing `(a, b, c)`.
    It is assumed that the latter is a valid key
    and that `(a, b, c)` is positive semi-definite. 
    """
    cdef int_triple tmp_triple

    mpz_set_si(dest, 0)

    tmp_triple = _reduce_GL(a, b, c)
    triple = (tmp_triple.a, tmp_triple.b, tmp_triple.c)
    try :
        mpz_set(dest, (<Integer>(coeffs_dict[triple])).value)
    except KeyError :
        pass

#===============================================================================
# get_coeff_generic
#===============================================================================

cdef get_coeff_generic(int a, int b, int c, coeffs_dict):
    r"""
    Return the value of ``coeffs_dict`` at the triple obtained from
    reducing `(a, b, c)`.
    It is assumed that the latter is a valid key
    and that `(a, b, c)` is positive semi-definite. 
    """
    cdef int_triple tmp_triple

    tmp_triple = _reduce_GL(a, b, c)
    triple = (tmp_triple.a, tmp_triple.b, tmp_triple.c)
    
    try :
        return coeffs_dict[triple]
    except KeyError :
        return None
