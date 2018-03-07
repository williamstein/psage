"""
Functions to reduce indices of Jacobi forms and to multiply expansions
of Jacobi forms.

AUTHORS:

    - Martin Raum (2010 - 04 - 05) Initial version
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


#include "cysignals/signals.pxi"
#include "stdsage.pxi"
include "cdefs.pxi"
from stdsage cimport PY_NEW
from cysignals.signals cimport sig_on,sig_off


from sage.rings.integer cimport Integer
from sage.rings.ring cimport Ring

cdef struct jac_index_sgn :
    int n
    int r
    int s

#####################################################################
#####################################################################
#####################################################################

#===============================================================================
# creduce
#===============================================================================

def creduce(k, m) :
    cdef jac_index_sgn res
    
    (n, r) = k
    res = _creduce(n, r, m)
    
    return ((res.n, res.r), res.s)
    
#===========================================================================
# _creduce
#===========================================================================

cdef jac_index_sgn _creduce(int n, int r, int m) :
    cdef int r_red, n_red, s
    cdef jac_index_sgn res
    
    r_red = r % (2 * m)
    if r_red > m :
        res.s = -1
        r_red = 2 * m - r_red
    else :
        res.s = 1
    
    n_red = n - (r**2 - r_red**2)/(4 * m)
    
    while n_red < 0 :
        r_red = r_red + 2 * m
        n_red = n + (4*m*r_red - 4*m**2)

    res.n = n_red
    res.r = r_red

    return res

#####################################################################
#####################################################################
#####################################################################

#===============================================================================
# mult_coeff_int
#===============================================================================

cpdef mult_coeff_int(result_key, coeffs_dict1, coeffs_dict2, ch1, ch2, result, int m) :
    r"""
    Returns the value at ``result_key`` of the coefficient dictionary of the product
    of the two Jacobi forms with dictionary ``coeffs_dict1`` and ``coeffs_dict2``.
    """
    cdef int n,r
    cdef int n1, r1
    cdef int n2, r2
    cdef int fm
    
    cdef mpz_t sqrt1, sqrt2, tmp, mpz_zero
    cdef mpz_t left, right
    cdef mpz_t acc
    
    mpz_init(sqrt1)
    mpz_init(sqrt2)
    mpz_init(tmp)
    mpz_init(mpz_zero)
    mpz_init(left)
    mpz_init(right)
    mpz_init(acc)
    
    (n,r) = result_key
    
    sig_on()
    fm = 4 * m
    for n1 in range(n + 1) :
        n2 = n - n1
        
        mpz_set_si(sqrt1, fm * n1)
        mpz_sqrt(sqrt1, sqrt1)
        
        mpz_set_si(sqrt2, fm * n2)
        mpz_sqrt(sqrt2, sqrt2)
        
        for r1 in range( r1 - mpz_get_si(sqrt2),
                         min( r1 + mpz_get_si(sqrt2) + 1,
                              mpz_get_si(sqrt1) + 1 ) ) :
            r2 = r - r1
            get_coeff_int(left, n1, r1, ch1, coeffs_dict1, m)
            if mpz_cmp(left, mpz_zero) == 0 : continue
            
            get_coeff_int(right, n2, r2, ch2, coeffs_dict2, m)
            if mpz_cmp(right, mpz_zero) == 0 : continue
            
            mpz_mul(tmp, left, right)
            mpz_add(acc, acc, tmp)
    sig_off()
    
    mpz_set((<Integer>result).value, acc)
    
    mpz_clear(sqrt1)
    mpz_clear(sqrt2)
    mpz_clear(tmp)
    mpz_clear(mpz_zero)
    mpz_clear(left)
    mpz_clear(right)
    mpz_clear(acc)
    
    return result
    
#===============================================================================
# mult_coeff_int_weak
#===============================================================================

cpdef mult_coeff_int_weak(result_key, coeffs_dict1, coeffs_dict2, ch1, ch2, result, m) :
    r"""
    Returns the value at ``result_key`` of the coefficient dictionary of the product
    of the two weak Jacobi forms with dictionary ``coeffs_dict1`` and ``coeffs_dict2``.
    """
    cdef int n,r
    cdef int n1, r1
    cdef int n2, r2
    cdef int fm, msq
    
    cdef mpz_t sqrt1, sqrt2, tmp, mpz_zero
    cdef mpz_t left, right
    cdef mpz_t acc
    
    mpz_init(sqrt1)
    mpz_init(sqrt2)
    mpz_init(tmp)
    mpz_init(mpz_zero)
    mpz_init(left)
    mpz_init(right)
    mpz_init(acc)
    
    (n,r) = result_key
    
    sig_on()
    fm = 4 * m
    msq = m**2
    for n1 in range(n + 1) :
        n2 = n - n1
        
        mpz_set_si(sqrt1, fm * n1 + msq)
        mpz_sqrt(sqrt1, sqrt1)
        
        mpz_set_si(sqrt2, fm * n2 + msq)
        mpz_sqrt(sqrt2, sqrt2)
        
        for r1 in range( r1 - mpz_get_si(sqrt2),
                         min( r1 + mpz_get_si(sqrt2) + 1,
                              mpz_get_si(sqrt1) + 1 ) ) :
            r2 = r - r1
            get_coeff_int(left, n1, r1, ch1, coeffs_dict1, m)
            if mpz_cmp(left, mpz_zero) == 0 : continue
            
            get_coeff_int(right, n2, r2, ch2, coeffs_dict2, m)
            if mpz_cmp(right, mpz_zero) == 0 : continue
            
            mpz_mul(tmp, left, right)
            mpz_add(acc, acc, tmp)
    sig_off()
    
    mpz_set((<Integer>result).value, acc)
    
    mpz_clear(sqrt1)
    mpz_clear(sqrt2)
    mpz_clear(tmp)
    mpz_clear(mpz_zero)
    mpz_clear(left)
    mpz_clear(right)
    mpz_clear(acc)
    
    return result

#==============================================================================
# mult_coeff_generic
#==============================================================================

cpdef mult_coeff_generic(result_key, coeffs_dict1, coeffs_dict2, ch1, ch2, result, int m) :
    r"""
    Returns the value at ``result_key`` of the coefficient dictionary of the product
    of the two Jacobi forms with dictionary ``coeffs_dict1`` and ``coeffs_dict2``.
    """
    cdef int n,r
    cdef int n1, r1
    cdef int n2, r2
    cdef int fm
    
    cdef mpz_t sqrt1, sqrt2, tmp, mpz_zero
    
    mpz_init(sqrt1)
    mpz_init(sqrt2)
    
    (n,r) = result_key
    
    sig_on()
    fm = 4 * m
    for n1 in range(n + 1) :
        n2 = n - n1
        
        mpz_set_si(sqrt1, fm * n1)
        mpz_sqrt(sqrt1, sqrt1)
        
        mpz_set_si(sqrt2, fm * n2)
        mpz_sqrt(sqrt2, sqrt2)
        
        for r1 in range( r1 - mpz_get_si(sqrt2),
                         min( r1 + mpz_get_si(sqrt2) + 1,
                              mpz_get_si(sqrt1) + 1 ) ) :
            r2 = r - r1
            left = get_coeff_generic(n1, r1, ch1, coeffs_dict1, m)
            if left is None : continue
            
            right = get_coeff_generic(n2, r2, ch2, coeffs_dict2, m)
            if right is None : continue
            
            result += left*right
    sig_off()

    mpz_clear(sqrt1)
    mpz_clear(sqrt2)
    
    return result

#===============================================================================
# mult_coeff_generic_weak
#===============================================================================

cpdef mult_coeff_generic_weak(result_key, coeffs_dict1, coeffs_dict2, ch1, ch2, result, int m) :
    r"""
    Returns the value at ``result_key`` of the coefficient dictionary of the product
    of the two Jacobi forms with dictionary ``coeffs_dict1`` and ``coeffs_dict2``.
    """
    cdef int n,r
    cdef int n1, r1
    cdef int n2, r2
    cdef int fm, msq
    
    cdef mpz_t sqrt1, sqrt2, tmp, mpz_zero
    
    mpz_init(sqrt1)
    mpz_init(sqrt2)
    
    (n,r) = result_key
    
    sig_on()
    fm = 4 * m
    msq = m**2
    for n1 in range(n + 1) :
        n2 = n - n1
        
        mpz_set_si(sqrt1, fm * n1 + msq)
        mpz_sqrt(sqrt1, sqrt1)
        
        mpz_set_si(sqrt2, fm * n2 + msq)
        mpz_sqrt(sqrt2, sqrt2)
        
        for r1 in range( r1 - mpz_get_si(sqrt2),
                         min( r1 + mpz_get_si(sqrt2) + 1,
                              mpz_get_si(sqrt1) + 1 ) ) :
            r2 = r - r1
            left = get_coeff_generic(n1, r1, ch1, coeffs_dict1, m)
            if left is None : continue
            
            right = get_coeff_generic(n2, r2, ch2, coeffs_dict2, m)
            if right is None : continue
            
            result += left*right
    sig_off()
    
    mpz_clear(sqrt1)
    mpz_clear(sqrt2)
    
    return result

#####################################################################
#####################################################################
#####################################################################

#===============================================================================
# get_coeff_int
#===============================================================================

cdef inline void get_coeff_int(mpz_t dest, int n, int r, int ch, coeffs_dict, int m) :
    cdef jac_index_sgn tmp_ind

    mpz_set_si(dest, 0)

    tmp_ind = _creduce(n, r, m)
    try :
        mpz_set(dest, (<Integer>(coeffs_dict[(tmp_ind.n, tmp_ind.r)])).value)
        if ch == -1 and tmp_ind.s == -1 :
            mpz_neg(dest, dest)
    except KeyError :
        pass

#===============================================================================
# get_coeff_generic
#===============================================================================

cdef inline get_coeff_generic(int n, int r, int ch, coeffs_dict, int m):
    cdef jac_index_sgn tmp_ind

    tmp_ind = _creduce(n, r, m)
    try :
        if ch == -1 and tmp_ind.s == -1 :
            return -coeffs_dict[(tmp_ind.n, tmp_ind.r)]
        else :
            return coeffs_dict[(tmp_ind.n, tmp_ind.r)]
    
    except KeyError :
        return None

