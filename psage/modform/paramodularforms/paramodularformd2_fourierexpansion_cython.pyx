r"""
Functions for reduction of Fourier indices and for multiplication of
paramodular modular forms.

AUTHORS:

- Martin Raum (2010 - 04 - 09) Intitial version
"""

#===============================================================================
# 
# Copyright (C) 2010 Martin Raum
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

include "interrupt.pxi"
include "stdsage.pxi"
include "cdefs.pxi"

include 'sage/ext/gmp.pxi'

from sage.rings.arith import gcd
from sage.rings.integer cimport Integer

## [a,b,c] a quadratic form; (u,x) an element of P1(ZZ/pZZ)
cdef struct para_index :
    int a
    int b
    int c
    int u
    int x 

#####################################################################
#####################################################################
#####################################################################

cdef dict _inverse_mod_N_cache = dict()

#####################################################################
#####################################################################
#####################################################################

#===============================================================================
# apply_GL_to_form
#===============================================================================

cpdef apply_GL_to_form(g, s) :
    """
    Return g s g^\tr if s is written as a matrix [a, b/2, b/2, c].
    """
    cdef int a, b, c,  u, x

    (u, x) = g
    (a,b,c) = s
    
    if u == 0 :
        return (c, b, a)
    else :
        return (a, b + 2*x*a, c + b*x + x**2 * a)
    
#===============================================================================
# reduce
#===============================================================================

def reduce_GL(index, p1list) :
    r"""
    A reduced form `((a,b,c), l)` satisfies `0 \le b \le a \le c`.
    """
    cdef int a, b, c, l, u, x, N
    cdef para_index res
    
    ((a,b,c), l) = index
    (u, x) = p1list[l]
    N = p1list.N()

    res = _reduce_ux(a, b, c, u, x, N)
    
    ## We use an exceptional rule for N = 1, since otherwise l = -1
    if N == 1 :
        return (((res.a, res.b, res.c), 0), 0)
    else :
        return (((res.a, res.b, res.c), p1list.index(res.u,res.x)), 0)

#===============================================================================
# _reduce_ux
#===============================================================================

cdef para_index _reduce_ux(int a, int b, int c, int u, int x, int N) :
    r"""
    We reduce the pair `(t, v)` by right action of `GL_2(ZZ)`.
    v is a left coset representative with respect to `(GL_2(ZZ))_0 (N)`.
    `(x, u)` is the second row of `v`. 
    
    NOTE:
        The quadratic form is reduced following Algorithm 5.4.2 found in Cohen's
        Computational Algebraic Number Theory.
    """
    cdef para_index res
    cdef int twoa, q, r
    cdef int tmp

    
    global _inverse_mod_N_cache
    if not N in _inverse_mod_N_cache :
        invN = PY_NEW(dict)
        
        for i in xrange(N) :
            if gcd(i,N) == 1 :
                invN[i] = Integer(i).inverse_mod(N)
        _inverse_mod_N_cache[N] = invN
        
    inverse_mod_N = _inverse_mod_N_cache[N] 

    if b < 0 :
        b = -b
        
        if u != 0 :
            x = (-x) % N
        
    while ( b > a or a > c ) :
        twoa = 2 * a
        
        #if b not in range(-a+1,a+1):
        if b > a :
            ## apply Euclidean step (translation)
            q = b // twoa #(2*a)
            r = b % twoa  #(2*a)
            if r > a :
                ## want (quotient and) remainder such that -a < r <= a
                r = r - twoa  # 2*a;
                q = q + 1
            c = c - b*q + a*q*q
            b = r
            
            if u != 0 :
                x = (x + q) % N
            
        ## apply symmetric step (reflection) if necessary
        if a > c:
            #(a,c) = (c,a)
            tmp = c
            c = a
            a = tmp
            
            if u == 0 :
                u = 1
                x = 0
            elif x == 0 :
                u = 0
                x = 1
            else :
                x = (-inverse_mod_N[x]) % N

        #b = abs(b)
        if b < 0:
            b = -b
            
            if u != 0 :
                x = (-x) % N
        ## iterate
    ## We're finished! Return the GL2(ZZ)-reduced form (a,b,c):

    res.a = a
    res.b = b
    res.c = c
    res.u = u
    res.x = x
            
    return res
