cimport psage.modform.maass.common_defs
from psage.modform.maass.common_defs cimport *

from sage.modules.free_module_element cimport FreeModuleElement
from psage.matrix.matrix_complex_dense cimport Matrix_complex_dense
from sage.rings.complex_mpc cimport MPComplexField_class,MPComplexNumber
from sage.rings.real_mpfr cimport RealNumber
from psage.rings.mpfr_nogil cimport *
from sage.structure.element cimport Vector

cdef class Vector_complex_dense(FreeModuleElement):
	cdef mpc_t* _entries
	cdef int _prec
	cdef MPComplexField_class _base_ring

	cdef _new_c(self,v=*)
	cdef _init(self, Py_ssize_t degree, parent)
	cpdef RealNumber norm(self,int ntype=?)
	cpdef base_ring(self)
	cdef MPComplexNumber _scalar_product_(self,Vector right)
	#def scalar_product(self,right)
