r"""
Various functions needed for some calculations in the package.

AUTHORS:

- Martin Raum (2009 - 08 - 27) Initial version
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

include "cysignals/signals.pxi"
include 'sage/ext/stdsage.pxi'
include "sage/ext/cdefs.pxi"
#include 'sage/ext/gmp.pxi'

from sage.rings.integer cimport Integer

cpdef divisor_dict(int precision) :
    r"""
    Return a dictionary of assigning to each `k <` ``precision`` a list of its divisors.
    
    INPUT :
    
    - precision    -- a positive integer
    """
    cdef int k
    cdef int l
    cdef int maxl
    cdef int divisor
    
    cdef mpz_t tmp
    
    mpz_init(tmp)
    
    mpz_set_si(tmp, precision)
    mpz_sqrt(tmp, tmp)
    bound = mpz_get_si(tmp) + 1
    
    mpz_clear(tmp)
    
    div_dict = PY_NEW(dict)
    
    # these are the trivial divisors
    div_dict[1] = [1]
    for k from 2 <= k < precision :
        div_dict[k] = [1,k]
        
    for k from 2 <= k < bound // 2 :
        divisor = 2*k
        
        for l from 2 <= l < bound // k :
            (<list>(div_dict[divisor])).append(k)
            divisor += k

    return div_dict
    
cpdef negative_fundamental_discriminants(int precision) :
    r"""
    Return a list of all negative fundamental discriminants `> -precision`
    """
    cdef int *markers = <int *>sage_malloc(precision * sizeof(int))
    cdef int k
    cdef int maxh
    cdef int hs

    cdef int tmp
    cdef mpz_t mpz_tmp
            
    mpz_init(mpz_tmp)
    
    sig_on()
    ## consider congruences mod 4
    for k from 1 <= k < precision :
        tmp = k % 16
        if tmp == 3 or tmp == 7 or tmp == 11 or tmp == 15 or tmp == 4 or tmp == 8 :
            markers[k] = 1
        else :
            markers[k] = 0

    mpz_set_si(mpz_tmp, precision)
    mpz_sqrt(mpz_tmp, mpz_tmp)
    maxh = mpz_get_si(mpz_tmp)
        
    ## consider square divisors
    for h from 3 <= h < maxh :
        hs = h*h 
        tmp = 0
        
        for k from 1 <= k < precision // hs :
            tmp += hs
            markers[tmp] = 0
        
    fund_discs = PY_NEW(list)
    for k from 2 <= k < precision :
        if markers[k] :
            fund_discs.append(-k)
    sig_off()
    
    mpz_clear(mpz_tmp)
                
    return fund_discs
 
