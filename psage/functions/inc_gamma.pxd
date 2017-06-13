#include "sage/ext/interrupt.pxi" 
include "sage/ext/stdsage.pxi"  
include "sage/ext/cdefs.pxi"
from psage.rings.mpfr_nogil cimport *

cdef int incgamma_pint_c(mpfr_t res,int n,mpfr_t x,int verbose=*) nogil
cdef int incgamma_hint_c(mpfr_t res,int n,mpfr_t x,int verbose=*) nogil
cdef int incgamma_phint_c(mpfr_t res, int n,mpfr_t x,int verbose=*) nogil
cdef int incgamma_nhint_c(mpfr_t res, int n,mpfr_t x,int verbose=*) nogil
cdef int incgamma_nint_c(mpfr_t res,int n,mpfr_t x,int verbose=*) nogil


