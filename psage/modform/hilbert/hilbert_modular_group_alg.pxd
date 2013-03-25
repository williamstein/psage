#include "sage/ext/cdefs.pxi"
#include "sage/libs/ntl/decl.pxi"

#include "../maass/math_inc.pxi"
#from psage.rings.double_prec_math cimport *
#from psage.modform.hilbert.hn_class cimport Hn
#cpdef get_closest_cusp(Hn z,G,int denom_max=?,int verbose=?)

#clang c++

#include "sage/ext/cdefs.pxi"
include "sage/libs/ntl/decl.pxi"
include 'sage/ext/stdsage.pxi'
include 'sage/ext/interrupt.pxi'
include "sage/rings/mpc.pxi"
cimport sage.structure.element
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.structure.element cimport Element, FieldElement, RingElement, ModuleElement
from sage.structure.parent_base cimport ParentWithBase
from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ

from sage.libs.ntl.ntl_ZZX cimport ZZX_c
from sage.libs.ntl.ntl_ZZ cimport ZZ_c

from sage.rings.number_field.number_field_element cimport NumberFieldElement, NumberFieldElement_absolute

