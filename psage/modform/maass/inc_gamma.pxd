include "sage/ext/interrupt.pxi" 
include "sage/ext/stdsage.pxi"  
include "sage/ext/cdefs.pxi"
#include "sage/rings/mpc.pxi"
#include "sage/ext/gmp.pxi"

from sage.libs.mpfr cimport *
include "../../rings/mpfr_loc.pxi" #ugly -- redefining things. but without the cimport above I get problems with mpfr_get_prec() and such

cdef int incgamma_pint_c(mpfr_t res,int n,mpfr_t x,int verbose=*) nogil
cdef int incgamma_hint_c(mpfr_t res,int n,mpfr_t x,int verbose=*) nogil
cdef int incgamma_phint_c(mpfr_t res, int n,mpfr_t x,int verbose=*) nogil
cdef int incgamma_nhint_c(mpfr_t res, int n,mpfr_t x,int verbose=*) nogil
cdef int incgamma_nint_c(mpfr_t res,int n,mpfr_t x,int verbose=*) nogil


