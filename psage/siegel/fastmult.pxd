include 'sage/ext/interrupt.pxi'
include 'sage/ext/stdsage.pxi'
include 'sage/ext/cdefs.pxi'

from sage.rings.ring cimport Ring
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
from sage.algebras.group_algebra import GroupAlgebra

cpdef mult_coeff_int(int a, int b, int c, coeffs_dict1, coeffs_dict2)
cpdef mult_coeff_int(int a, int b, int c, coeffs_dict1, coeffs_dict2)
cpdef mult_coeff_generic(int a, int b, int c, coeffs_dict1, coeffs_dict2, Ring R)
