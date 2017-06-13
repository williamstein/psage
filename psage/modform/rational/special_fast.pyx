#################################################################################
#
# (c) Copyright 2011 William Stein
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
Cython code for very fast high precision computation of certain specific modular forms of interest.
"""


########################################

# The following functions are here mainly because computing f(q^d),
# for f(q) a power series, is by default quite slow in Sage, since
# FLINT's polynomial composition code is slow in this special case.

from sage.rings.rational cimport Rational
from sage.rings.integer cimport Integer
from sage.rings.polynomial.polynomial_rational_flint cimport Polynomial_rational_flint
from sage.libs.flint.fmpq_poly cimport (fmpq_poly_get_coeff_mpq, fmpq_poly_set_coeff_mpq,
                                        fmpq_poly_length)

from sage.libs.flint.fmpz_poly cimport fmpz_poly_set_coeff_mpz
from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint
                                                                   

from sage.libs.gmp.mpq cimport mpq_numref

from sage.rings.all import ZZ

def _change_ring_ZZ(Polynomial_rational_flint f):
    """
    Return the polynomial of numerators of coefficients of f.  
    
    INPUT:
        - f -- a polynomial over the rational numbers with no denominators
    OUTPUT:
        - a polynomial over the integers

    EXAMPLES::

        sage: import psage.modform.rational.special_fast as s
        sage: R.<q> = QQ[]
        sage: f = 3 + q + 17*q^2 - 4*q^3 - 2*q^5
        sage: g = s._change_ring_ZZ(f); g
        -2*q^5 - 4*q^3 + 17*q^2 + q + 3
        sage: g.parent()
        Univariate Polynomial Ring in q over Integer Ring

    Notice that the denominators are just uniformly ignored::

        sage: f = 3/2 + -5/8*q + 17/3*q^2
        sage: s._change_ring_ZZ(f)
        17*q^2 - 5*q + 3
    """
    cdef Polynomial_integer_dense_flint res = ZZ[f.parent().variable_name()](0)
    cdef Rational x = Rational(0)
    cdef unsigned long i
    for i in range(fmpq_poly_length(f.__poly)):
        fmpq_poly_get_coeff_mpq(x.value, f.__poly, i)
        fmpz_poly_set_coeff_mpz(res.__poly, i, mpq_numref(x.value))
    return res

def _evaluate_poly_at_power_of_gen(Polynomial_rational_flint f, unsigned long n, bint truncate):
    """
    INPUT:
        - f -- a polynomial over the rational numbers
        - n -- a positive integer
        - truncate -- bool; if True, truncate resulting polynomial so
          it has coefficients up to deg(f) and possibly slightly more,
          e.g., this is natural if we are interested in the power
          series f+O(x^(d+1)).
          
    OUTPUT:
        - a polynomial over the rational numbers

    EXAMPLES::

        sage: import psage.modform.rational.special_fast as s
        sage: R.<q> = QQ[]
        sage: f = 2/3 + q + 17*q^2 - 4*q^3 - 2*q^5
        sage: s._evaluate_poly_at_power_of_gen(f, 2, False)
        -2*q^10 - 4*q^6 + 17*q^4 + q^2 + 2/3
        sage: s._evaluate_poly_at_power_of_gen(f, 2, True)
        -4*q^6 + 17*q^4 + q^2 + 2/3
    """
    if n == 0:
        raise ValueError, "n must be positive"
    cdef Polynomial_rational_flint res = f._new()
    cdef unsigned long k, m
    cdef Rational z = Rational(0)
    m = fmpq_poly_length(f.__poly)
    if truncate:
        m = m // n
    for k in range(m+1):
        fmpq_poly_get_coeff_mpq(z.value, f.__poly, k)
        fmpq_poly_set_coeff_mpq(res.__poly, n*k, z.value)
    return res

def _evaluate_series_at_power_of_gen(h, unsigned long n, bint truncate):
    """
    INPUT:
        - h -- a power series over the rational numbers
        - n -- a positive integer
        - if  
        
    OUTPUT:
        - a power series over the rational numbers

    EXAMPLES::
    
        sage: import psage.modform.rational.special_fast as s
        sage: R.<q> = QQ[[]]
        sage: f = 2/3 + 3*q + 14*q^2 + O(q^3)
        sage: s._evaluate_series_at_power_of_gen(f, 2, True)
        2/3 + 3*q^2 + O(q^3)
        sage: s._evaluate_series_at_power_of_gen(f, 2, False)
        2/3 + 3*q^2 + 14*q^4 + O(q^5)    
    """
    # Return h(q^n) truncated to the precision of h massively quicker
    # than writing "h(q^n)" in Sage.
    f = _evaluate_poly_at_power_of_gen(h.polynomial(), n, truncate=truncate)
    prec = h.prec() if truncate else (h.prec()-1) * n + 1
    return h.parent()(f, prec)



def Integer_list_to_polynomial(list v, var='q'):
    cdef Polynomial_integer_dense_flint res = ZZ[var](0)
    cdef Py_ssize_t i
    cdef Integer a
    for i in range(len(v)):
        a = v[i]
        fmpz_poly_set_coeff_mpz(res.__poly, i, a.value)
    return res
    
