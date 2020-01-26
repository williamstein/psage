r"""
Collection of includes of double precision mathematical functions

"""
cdef extern from "pari/paricfg.h":
    pass
cdef extern from "pari/pari.h":
    pass
cdef extern from "pari/paripriv.h":
    pass


cdef extern from "math.h":
    cdef double fabs(double)
    cdef double fmax(double,double)
    cdef int ceil(double)
    cdef int floor(double) 
    cdef double sqrt(double)
    cdef double sin(double)
    cdef double cos(double)
    cdef double log(double)
    cdef double power(double,double)

cdef extern from "complex.h":
    ### "Hack" suggested by Robert Bradshaw in 2009.
    ### TODO: Is there now a better way for pure c double complex?
    ctypedef double cdouble "double complex"
    cdef double creal(double complex)
    cdef double cimag(double complex)
    cdef double complex _Complex_I
    cdef double carg(double complex)
    cdef double cabs(double complex)
    cdef double complex cexp(double complex)
    cdef double complex csqrt(double complex)
    cdef double complex cpow(double complex,double complex)

cdef double complex _I = _Complex_I
