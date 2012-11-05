include "sage/ext/interrupt.pxi" 
include "sage/ext/stdsage.pxi"  
include "sage/ext/cdefs.pxi"
include "sage/rings/mpc.pxi"
from psage.matrix.matrix_complex_dense cimport Matrix_complex_dense
from sage.rings.complex_mpc cimport MPComplexField_class
from sage.rings.real_mpfr cimport RealField_class,RealNumber
from sage.rings.rational cimport Rational
cdef class LinearSystem(object):
    cdef int _verbose
    cdef int _prec
    cdef int _M0
    cdef Rational _weight
    cdef Matrix_complex_dense _matrix
    cdef object _multiplier
    cdef MPComplexField_class _CF
    cdef RealField_class _RF
    cpdef  setup_system_vv(self,int M0,RealNumber Y)
