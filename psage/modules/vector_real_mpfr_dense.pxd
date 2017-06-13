include '../ext/stdsage.pxi'
include '../ext/cdefs.pxi'
#include "../rings/mpc.pxi"
#include "../ext/gmp.pxi"

from sage.modules.free_module_element cimport FreeModuleElement

#from sage.rings.complex_mpc cimport *
from sage.rings.real_mpfr cimport *
from sage.structure.element cimport Vector,Element

cdef class Vector_real_mpfr_dense(FreeModuleElement):
        cdef mpfr_t* _entries
        cdef int _prec
        cdef RealField_class _base_ring
        cpdef _add_(self, other)
        cpdef _sub_(self, other)
        cdef _new_c(self,v=*)
        cdef _init(self, Py_ssize_t degree, parent)
        cpdef RealNumber norm(self,int ntype=?)
        cpdef base_ring(self)
        cdef RealNumber _scalar_product_(self,Vector right)
        #def scalar_product(self,right)
        cdef int _cmp_c_impl(left, Element right) except -2