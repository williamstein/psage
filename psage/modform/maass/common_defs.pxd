
include 'sage/ext/stdsage.pxi' 
include "sage/ext/cdefs.pxi"
include 'sage/ext/interrupt.pxi'
from sage.libs.mpfr cimport *

from sage.rings.complex_mpc cimport * 
from sage.rings.complex_mpc cimport MPComplexNumber
from sage.rings.real_mpfr cimport RealNumber,RealField_class
from sage.rings.complex_number cimport ComplexNumber
cimport sage.structure.element
from sage.structure.element cimport Element, ModuleElement, RingElement, Vector
from sage.rings.integer cimport Integer
from sage.rings.complex_mpc cimport MPComplexNumber
from sage.rings.complex_mpc cimport MPComplexField_class
from sage.modules.free_module_element cimport FreeModuleElement
from sage.rings.real_mpfr cimport RealNumber

cdef mpc_rnd_t rnd
cdef mpfr_rnd_t rnd_re
rnd = MPC_RNDNN
rnd_re = MPFR_RNDN
