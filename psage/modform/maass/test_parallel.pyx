include "stdsage.pxi"  
include "cdefs.pxi"

cdef extern from "mpfr.h" nogil:
    ctypedef short mpfr_prec_t
    ctypedef unsigned short mpfr_uprec_t
    ctypedef int mpfr_sign_t
    ctypedef short mpfr_exp_t
    ctypedef unsigned short mpfr_uexp
    
    struct __mpfr_struct:
        mpfr_prec_t  _mpfr_prec
        mpfr_sign_t  _mpfr_sign
        mpfr_exp_t   _mpfr_exp
        mp_limb_t   *_mpfr_d

    ctypedef enum mpfr_rnd_t:
        MPFR_RNDN=0,
        MPFR_RNDZ,
        MPFR_RNDU,
        MPFR_RNDD,
        MPFR_RNDA,
        MPFR_RNDF,
        MPFR_RNDNA=-1
        
    ctypedef __mpfr_struct mpfr_t[1]
    ctypedef __mpfr_struct *mpfr_ptr
    
    ctypedef __mpfr_struct *mpfr_srcptr
    
    void mpfr_init2 (mpfr_ptr, mpfr_prec_t)
    void mpfr_init (mpfr_ptr)
    void mpfr_clear (mpfr_ptr)
    int mpfr_const_pi (mpfr_ptr, mpfr_rnd_t)
    int mpfr_mul_ui (mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t)
    int mpfr_sqrt (mpfr_ptr, mpfr_srcptr, mpfr_rnd_t)
    int mpfr_add (mpfr_ptr, mpfr_srcptr,mpfr_srcptr, mpfr_rnd_t)
    mpfr_prec_t mpfr_get_prec(mpfr_srcptr)
    int mpfr_add_ui (mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t)
    void mpfr_set4 (mpfr_ptr, mpfr_srcptr, mpfr_rnd_t, int)
    ctypedef int mpfr_sign_t
    void mpfr_set (mpfr_ptr, mpfr_srcptr, mpfr_rnd_t)
    unsigned long mpfr_get_ui (mpfr_srcptr, mpfr_rnd_t)
    int mpfr_set_ui (mpfr_ptr, unsigned long, mpfr_rnd_t)
    int mpfr_cos (mpfr_ptr, mpfr_srcptr,mpfr_rnd_t)
    int mpfr_sin (mpfr_ptr, mpfr_srcptr,mpfr_rnd_t)
    int mpfr_exp (mpfr_ptr, mpfr_srcptr,mpfr_rnd_t)
    double mpfr_get_d (mpfr_srcptr, mpfr_rnd_t)
    
import cython
from libc.stdlib cimport abort, malloc, free
    
from cython.parallel cimport parallel, prange
from sage.rings.real_mpfr cimport RealNumber, RealField_class
from sage.rings.real_mpfr import RealField

cdef mpfr_rnd_t rnd_re
rnd_re = MPFR_RNDN

cdef mpfr_t * ef1cosv = NULL
#cpdef int * ef1cosv = NULL
ef1cosv = <mpfr_t *> malloc( sizeof(mpfr_t) * (100000) )

cimport openmp
openmp.omp_set_dynamic(1)
#openmp.omp_set_num_threads(5)

cpdef test_par_loop(int m):
    global ef1cosv
    par_loop(m,ef1cosv)
    #for n in range(m):
    #    print mpfr_get_d(ef1cosv[n],rnd_re)

cpdef free_e1():
    free(ef1cosv)

cpdef test_loop_seq(m):
    global e1cosv
    return seq_loop(m, ef1cosv)

cdef int par_loop(int m, mpfr_t * l):
    cdef int n = 0
    cdef int k = 0
    #for n in range(m):
        #mpfr_init(ef1cosv[n])
        #mpfr_set_ui(ef1cosv[n],n,rnd_re)
        #test.append(0)
    #    l[n]=int(0)
    #for n in range(m):
    #    print n, ef1cosv[n]
    with nogil, parallel():
        for n in prange(m):
            mpfr_init(l[n])
            mpfr_set_ui(l[n],n+openmp.omp_get_thread_num(),rnd_re)
            mpfr_cos(l[n],l[n],rnd_re)
            for k in range(1000):
                mpfr_cos(l[n],l[n],rnd_re)
                mpfr_sin(l[n],l[n],rnd_re)
                mpfr_exp(l[n],l[n],rnd_re)  
    return 0

cdef int seq_loop(int m, mpfr_t * l):
    cdef int n = 0
    cdef int k = 0
    for n in range(m):
        mpfr_init(l[n])
        mpfr_set_ui(l[n],n+openmp.omp_get_thread_num(),rnd_re)
        mpfr_cos(l[n],l[n],rnd_re)
        for k in range(1000):
            mpfr_cos(l[n],l[n],rnd_re)
            mpfr_sin(l[n],l[n],rnd_re)
            mpfr_exp(l[n],l[n],rnd_re)
    return 0
