include "sage/ext/interrupt.pxi" 
include "sage/ext/stdsage.pxi"  
include "sage/ext/cdefs.pxi"
include "sage/rings/mpc.pxi"
#include "sage/ext/gmp.pxi"

from sage.libs.mpfr cimport mpfr_t

cdef incgamma_pint_c(mpfr_t res,int n,mpfr_t x,int verbose=*)
cdef incgamma_hint_c(mpfr_t res,int n,mpfr_t x,int verbose=*)


cdef incgamma_nint_c(mpfr_t res,int n,mpfr_t x,int verbose=*)


