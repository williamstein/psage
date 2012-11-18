include "stdsage.pxi"  
include "cdefs.pxi"
#include "mpfr_loc.pxi"

from sage.libs.mpfr cimport *

include "../../rings/mpfr_loc.pxi"

import cython
from libc.stdlib cimport abort, malloc, free
    
from cython.parallel cimport parallel, prange

cdef mpfr_rnd_t rnd_re
rnd_re = MPFR_RNDN

cdef mpfr_t * ef1cosv = NULL
#cpdef int * ef1cosv = NULL
ef1cosv = <mpfr_t *> malloc( sizeof(mpfr_t) * (100000) )

cdef int nmax=10000

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

## Gamma(n,x) with n<0 and real x>0
cdef int incgamma_nint_c(mpfr_t res, int n, mpfr_t x):
    r"""
    Incomplete Gamma, Gamma(n,x) with n negative integer.
    We use Gamma(0,x)=-Ei(-x) for x>0
    """
    if n!=0:
        return -1
    cdef mpfr_prec_t prec
    cdef mpfr_prec_t wp
    prec=mpfr_get_prec(x)
    cdef int ok = 1
    cdef mpfr_t xabs_t, wp_t
    wp=prec+20
    mpfr_init2(xabs_t, prec)
    mpfr_init2(wp_t, prec)
    mpfr_set_ui(wp_t,wp,rnd_re)
    mpfr_abs(xabs_t, x, rnd_re)
    mpfr_get_ui(x,rnd_re)
    #cdef int xabs = mpfr_get_ui(xabs_t,rnd_re)
    cdef mpfr_t xnew
    cdef int rh
    mpfr_mul_d(wp_t, wp_t,0.693, rnd_re)
    rh = mpfr_get_ui(wp_t,rnd_re)
    if mpfr_cmp_ui(xabs_t,rh) > 0:
        mpfr_init2(xnew,wp)
        mpfr_neg(xnew,x,rnd_re)
        ok = ei_asymp_c(res,xnew)
    else:
        mpfr_mul_ui(xabs_t,xabs_t,2,rnd_re)
        mpfr_set_ui(wp_t,wp,rnd_re)
        mpfr_add(wp_t,wp_t,xabs_t,rnd_re)
        mpfr_init2(xnew,mpfr_get_ui(wp_t,rnd_re))
        mpfr_set_prec(res,mpfr_get_ui(wp_t,rnd_re))
        mpfr_neg(xnew,x,rnd_re)
        ok = ei_taylor_c(res,xnew)
    mpfr_neg(res,res,rnd_re)
    mpfr_clear(xabs_t)
    mpfr_clear(wp_t)
    mpfr_clear(xnew)
    return ok


cdef int ei_asymp_c(mpfr_t res, mpfr_t x):
    r"""
    Compute the exponential integral of x via asymptotic formula.
    """
    cdef mpfr_t tmp,tmp2,summa,r,eps
    cdef int k
    cdef mpfr_prec_t prec
    #cdef RealField_class RF
    cdef int ok = 1
    prec = mpfr_get_prec(x)
    #RF = RealField(prec)
    #tmp=RF(1); summa=RF(1); r=RF(1); tmp2=RF(0)
    mpfr_init2(eps,prec)
    mpfr_init2(tmp,prec)
    mpfr_set_ui(tmp,1,rnd_re)
    mpfr_init2(summa,prec)
    mpfr_init2(r,prec)
    mpfr_init2(tmp2,prec)
    #eps = RF((2.**-(prec+1)))
    mpfr_set_si_2exp(eps,1,-prec-1,rnd_re)
    mpfr_set(summa,tmp,rnd_re)
    mpfr_div(r,tmp,x,rnd_re)
    mpfr_exp(tmp2,x,rnd_re)
    mpfr_mul(tmp2,tmp2,r,rnd_re)
    #if abs(tmp2)<eps:
    #    mpfr_set(res,tmp2.value,rnd_re)
    #    return
    mpfr_div(eps,eps,tmp2,rnd_re)
    mpfr_abs(eps,eps,rnd_re)
    #if verbose>0:
    #    print "r = 1/x = ", r
    #    print "eps = ", eps
    #    print "exp(x)/x = ", tmp2
    #    print "tmp = ", tmp
    for k in range(1,nmax): #from 1 <= k <= nmax:
        mpfr_mul_ui(tmp,tmp,k,rnd_re)
        mpfr_mul(tmp,tmp,r,rnd_re)
        mpfr_add(summa,summa,tmp,rnd_re)
        #if verbose>0:
        #    print "k = ", k
        #    print "r = ", r
        #    print "tmp = ", tmp
        #    print "summa = ", summa
        if  mpfr_cmpabs(tmp,eps)<0:
            #if verbose>0:
            #    print 'break at k=', k
            break
    if k>= nmax:
        mpfr_set(summa,x,rnd_re)
    else:
        ok = 0
        #raise ArithmeticError,"k>nmax! in Ei(x)!, error of order {0} for x={1}".format(tmp,summa)
    mpfr_mul(res,summa,tmp2,rnd_re)
    mpfr_clear(tmp);mpfr_clear(tmp2)
    mpfr_clear(summa);mpfr_clear(r)
    mpfr_clear(eps)
    return ok


cdef int ei_taylor_c(mpfr_t res, mpfr_t x):
    cdef int ok=1
    cdef mpfr_t lnx
    cdef mpfr_prec_t prec = mpfr_get_prec(x)
    mpfr_init2(lnx, prec)
    ok = Ei_ml_c(res, x)
    #if verbose>0:
    #    print  "Ei(x)-log(x)={0}, prec={1}".format(mpfr_get_ld(res,rnd_re), prec)
    mpfr_abs(x,x,rnd_re)
    mpfr_log(lnx,x,rnd_re)
    #if verbose>0:
    #    print  "ln(x)={0}, prec={1}".format(mpfr_get_ld(lnx,rnd_re), prec)
    mpfr_add(res,res,lnx,rnd_re)
    mpfr_clear(lnx)
    return ok


cdef int Ei_ml_c(mpfr_t res,mpfr_t x):
    r"""
    Compute the exponential integral of x  - ln|x|
    """
    cdef mpfr_t tmp,summa,eps
    cdef int k
    cdef mpfr_prec_t prec
    cdef mpfr_t ec
    cdef int ok = 1
    prec=mpfr_get_prec(x)
    mpfr_init2(ec,prec)
    #RF = RealField(prec)
    #tmp=RF(1); summa=RF(0)
    mpfr_init2(tmp,prec)
    mpfr_init2(eps,prec)
    mpfr_init2(summa,prec)
    mpfr_set_si(summa,0,rnd_re)
    mpfr_set_si_2exp(eps,1,-prec-20,rnd_re)
    #eps = RF(2.0**-((prec+20)))
    #print "eps={0}, prec={1}".format(eps, prec)
    #call set_eulergamma()
    mpfr_set(tmp,x,rnd_re)
    #summa=tmp
    mpfr_set(summa,tmp,rnd_re)
    for k in range(2,nmax+1): #from 2 <= k <= nmax:
        #kk=RF(k)
        ##tmp=tmp*(kk-RF(1))/kk**2
        mpfr_mul_ui(tmp,tmp,k-1,rnd_re)
        mpfr_div_ui(tmp,tmp,k*k,rnd_re)
        mpfr_mul(tmp,tmp,x,rnd_re)
        mpfr_add(summa,summa,tmp,rnd_re)
        #print "tmp=",tmp
        if mpfr_cmpabs(tmp,eps)<0:
            break
    if k>= nmax:
        mpfr_set(summa,x,rnd_re)
    else:
        ok = 0 #return 0
        #raise ArithmeticError,"k>nmax! in Ei(x)!, error of order {0} for x={1}".format(tmp,summa)
    mpfr_const_euler(ec,rnd_re)
    mpfr_add(summa,summa,ec,rnd_re)
    mpfr_set(res,summa,rnd_re)
    mpfr_clear(ec)
    mpfr_clear(summa)
    mpfr_clear(tmp)
    mpfr_clear(eps)
    #summa=summa+RF.euler_constant()
    #print 'Ei(',x,')=',summa
    return ok

