# encoding: utf-8
# cython: profile=False
# -*- coding=utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2012
#
#  Authors:
#  Fredrik Strömberg <stroemberg@mathematik.tu-darmstadt.de>
#  Stephan Ehlen <ehlen@mathematik.tu-darmstadt.de>
#  
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
r"""
Algorithms for computing the incomplete gamma function (and related functions)



Conventions for return codes:
 '0' - success
 '1' - failure of reaching precision
 '-1'- error in input parameters
 
"""

## For multiprecision support
cdef mpfr_rnd_t rnd_re
rnd_re = MPFR_RNDN
cdef int nmax = 10000
from sage.rings.real_mpfr cimport RealNumber
from sage.rings.real_mpfr cimport RealField_class
from sage.rings.real_mpfr import RealField
from sage.functions.special import error_fcn as erfc


cdef extern from "math.h":
    double fabs(double)
    double fmax(double,double)
    int ceil(double)
    double exp(double)
    double M_LN10
    double log(double)

import cython

####
####
#### The incomplete Gamma function of integer parameter
####

cpdef RealNumber incgamma_int(int n,RealNumber x,int verbose=0):
    cdef RealNumber res
    cdef int ok = 1
    res = x.parent()(0)
    if n>0:
        if verbose>2:
            print "In incgamma_int! Calling pint!"
        ok = incgamma_pint_c(res.value,n,x.value,verbose)
    else:
        if verbose>2:
            print "In incgamma_int! Calling nint!"
        ok = incgamma_nint_c(res.value,n,x.value,verbose)
    if verbose>2:
        print "Got ok={0}!".format(ok)
    if ok == 0:
        return res
    else:
        raise ArithmeticError,"Incomplete Gamma failed with code:{0}".format(ok)
    
####  The incomplete gamma function of half-integral parameter
cpdef RealNumber incgamma_hint(int n,RealNumber x,int verbose=0):
    r"""
    Incomplete Gamma function of half-integer parameter, Gamma(n+1/2,x).

    INPUT:
    -`n` -- integer
    -`x` -- real number
    -`verbose` -- integer

    """
    cdef RealNumber res
    cdef int ok = 1
    res = x.parent()(0)
    ok = incgamma_hint_c(res.value,n,x.value)
    if ok == 0:        
        return res
    else:
        raise ArithmeticError,"Incomplete Gamma failed with code:{0}".format(ok)
    
## The exponential integral
cpdef RealNumber Ei_ml(RealNumber x):
    r"""
    Compute Ei(x)-ln|x| where Ei(x) is the exponential integral of x.    
    """
    cdef RealNumber res
    cdef int ok = 1
    res = x.parent()(0)
    ok = Ei_ml_c(res.value,x.value)
    if ok == 0:        
        return res
    else:
        raise ArithmeticError,"Exponential integral failed with code:{0}".format(ok)

cpdef RealNumber ei(RealNumber x, int verbose=0):
    cdef int wp = x.prec() + 20
    cdef int xabs= abs(x.integer_part())
    cdef RealField_class RF_orig = x.parent()
    cdef RealField_class RF
    cdef int rh = RF_orig(wp*0.693).integer_part() + 10
    cdef int ok = 1
    if xabs > rh:
        if verbose>0:
            print "can use asymptotic"
        RF=RealField(wp)
        return RF_orig(ei_asymp(RF(x),verbose))
    else:
        if verbose>0:
            print "can NOT use asymptotic"
        RF=RealField(wp+2*xabs)
        return RF_orig(ei_taylor(RF(x), verbose))

cpdef RealNumber ei_taylor(RealNumber x, int verbose=0):
    r"""
    Compute the exponential integral of x  - ln|x|
    """
    cdef RealNumber res
    cdef int ok = 1
    res = x.parent()(0)
    ok = ei_taylor_c(res.value,x.value,verbose)
    if ok == 0:
        return res
    else:
        raise ArithmeticError,"Taylor expansion method for exponential integral failed with code:{0}".format(ok)
    
cpdef RealNumber ei_asymp(RealNumber x, int verbose=0):
    r"""
    Compute the exponential integral of x via asymptotic formula
    """
    cdef RealNumber res
    cdef ok = 1
    res = x.parent()(0)
    ok = ei_asymp_c(res.value,x.value,verbose)
    if ok == 0:
        return res
    else:
        raise ArithmeticError,"Asymptotic method for the exponential integral failed with code:{0}".format(ok)
    
###
###  cdef'd functions
### 
#cdef RealNumber incgamma_pint(int n,RealNumber x):
## Changed that these return 1 on success and 0 on fail (to avoid try -- catch clauses)
cdef int incgamma_pint_c(mpfr_t res, int n,mpfr_t x,int verbose=0) nogil:
    r"""
    Incomplete Gamma, Gamma(n,x) with n positive integer.
    Uses the formula:
    Gamma(n,x)=(n-1)!e^(-x) \sum_{k=0}^{n-1} x^k/k!
    This is just a finite sum so the only error should come from precision of the mpfr-functions and operations.
    """
    cdef mpfr_t tmp,tmp2,tmp_fak
    cdef int j
    cdef mpfr_prec_t prec
    prec = mpfr_get_prec(x)
    #res = RF(-x)
    mpfr_neg(res,x,rnd_re)
    #print "prec=",prec
    mpfr_init2(tmp,prec)
    mpfr_init2(tmp2,prec)
    mpfr_init2(tmp_fak,prec)
    mpfr_set_ui(tmp,1,rnd_re)
    mpfr_set_ui(tmp2,1,rnd_re)
    mpfr_set_ui(tmp_fak,1,rnd_re)
    for j from 1 <=j <= n-1:
        #jj=<double>j
        mpfr_mul(tmp2,tmp2,x,rnd_re)
        mpfr_div_ui(tmp2,tmp2,j,rnd_re)
        #tmp2=tmp2*x/jj
        mpfr_add(tmp,tmp,tmp2,rnd_re)
        #tmp=tmp+tmp2
        mpfr_mul_ui(tmp_fak,tmp_fak,j,rnd_re)
        #tmp_fak=tmp_fak*jj
    # res = exp(-x)*tmp_fak*tmp
    # res = -x
    #mpfr_neg(res.value,x,rnd_re)
    #mpfr_exp(res.value,res.value,rnd_re)
    mpfr_exp(res,res,rnd_re)
    mpfr_mul(tmp,tmp,tmp_fak,rnd_re)
    mpfr_mul(res,res,tmp,rnd_re)
    mpfr_clear(tmp)
    mpfr_clear(tmp2)
    mpfr_clear(tmp_fak)
    return 0

## Gamma(n,x) with n<0 and real x>0
cdef int incgamma_nint_c(mpfr_t res, int n,mpfr_t x,int verbose=0) nogil:
    r"""
    Incomplete Gamma, Gamma(n,x) with n negative integer.
    We use Gamma(0,x)=-Ei(-x) for x>0
    """
    if n!=0:
        return -1
    cdef mpfr_prec_t prec, wp
    cdef int ok = 1
    cdef mpfr_t xabs_t, wp_t
    prec = mpfr_get_prec(x)
    wp=prec+20
    mpfr_init2(xabs_t, prec)
    mpfr_init2(wp_t, prec)
    mpfr_set_ui(wp_t,wp,rnd_re)
    mpfr_abs(xabs_t, x, rnd_re)
    #cdef int xabs = mpfr_get_ui(xabs_t,rnd_re)
    cdef mpfr_t xnew
    cdef int rh
    mpfr_mul_d(wp_t, wp_t,0.693, rnd_re)
    rh = mpfr_get_ui(wp_t,rnd_re)
    if mpfr_cmp_ui(xabs_t,rh) > 0:
        mpfr_init2(xnew,wp)
        mpfr_neg(xnew,x,rnd_re)
        ok = ei_asymp_c(res,xnew,verbose)
    else:
        mpfr_mul_ui(xabs_t,xabs_t,2,rnd_re)
        mpfr_set_ui(wp_t,wp,rnd_re)
        mpfr_add(wp_t,wp_t,xabs_t,rnd_re)
        mpfr_init2(xnew,mpfr_get_ui(wp_t,rnd_re))
        mpfr_set_prec(res,mpfr_get_ui(wp_t,rnd_re))
        mpfr_neg(xnew,x,rnd_re)
        ok = ei_taylor_c(res,xnew,verbose)
    mpfr_neg(res,res,rnd_re)
    mpfr_clear(xabs_t)
    mpfr_clear(wp_t)
    mpfr_clear(xnew)
    return ok


cdef int ei_asymp_c(mpfr_t res, mpfr_t x, int verbose=0) nogil:
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


cdef int ei_taylor_c(mpfr_t res, mpfr_t x, int verbose=0) nogil:
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


cdef int Ei_ml_c(mpfr_t res,mpfr_t x) nogil:
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

cdef int incgamma_hint_c(mpfr_t res,int n,mpfr_t x,int verbose=0) nogil:
    cdef mpfr_t sqpi,sqx
    cdef int ok = 1
    cdef mpfr_prec_t prec
    if n > 0:
        #return
        ok = incgamma_phint_c(res,n,x,verbose)
    elif n<0:
        ok = incgamma_nhint_c(res,-n,x,verbose)
    else:
        prec = mpfr_get_prec(x)
        mpfr_init2(sqpi,prec)
        mpfr_init2(sqx,prec)
        mpfr_const_pi(sqpi,rnd_re)
        mpfr_sqrt(sqpi,sqpi,rnd_re)
        mpfr_sqrt(sqx,x,rnd_re)
        mpfr_erfc(res,sqx,rnd_re)
        mpfr_mul(res,res,sqpi,rnd_re)
        ok = 0
    return ok
### sqrt(pi)* erfc(sqrt(x))  is used at several places

cdef int mpfr_sqpi_erfc_sqx(mpfr_t res,mpfr_t x,int verbose=0) nogil:
    cdef mpfr_prec_t prec = mpfr_get_prec(x)
    cdef mpfr_t sqpi
    mpfr_init2(sqpi,prec)
    mpfr_const_pi(sqpi,rnd_re)
    mpfr_sqrt(sqpi,sqpi,rnd_re)
    mpfr_sqrt(res,x,rnd_re)
    mpfr_erfc(res,res,rnd_re)
    mpfr_mul(res,res,sqpi,rnd_re)
    mpfr_clear(sqpi)
    return 0


cdef int incgamma_phint_c(mpfr_t res, int n,mpfr_t x,int verbose=0) nogil:
    cdef int j,m
    cdef mpfr_prec_t prec
    cdef mpfr_t half_m_n,summa,tmp,term,tmp2,sqx
    cdef int ok = 0
    #assert n>=0
    if n < 0:
        return -1
    #RF = x._parent
    prec = mpfr_get_prec(x) #RF.__prec
    #summa=RF(0); tmp=RF(1); term=RF(0)
    #tmp2=RF(0)
    #nn=RF(n)
    mpfr_init2(half_m_n,prec)
    mpfr_init2(tmp,prec)
    mpfr_init2(sqx,prec)
    mpfr_init2(tmp2,prec)
    mpfr_init2(term,prec)
    mpfr_init2(summa,prec)
    #!! do the sum first
    mpfr_set_ui(tmp,1,rnd_re)
    mpfr_set_ui(term,0,rnd_re)
    mpfr_set_ui(summa,0,rnd_re)
    mpfr_set_d(half_m_n,0.5,rnd_re) # = RF(0.5)
    mpfr_sub_si(half_m_n,half_m_n,n,rnd_re)
    for j from 0 <= j <= n-1:
        _mppochammer_mpfr(term,half_m_n,n-j-1)
        mpfr_mul(term,term,tmp,rnd_re)
        mpfr_add(summa,summa,term,rnd_re)
        mpfr_mul(tmp,tmp,x,rnd_re)
        mpfr_neg(tmp,tmp,rnd_re)
    #m=2*n+1
    mpfr_set_ui(tmp,1,rnd_re)
    #for j from 1 <= j <= 2*n-1:
    for j from 1 <= j <= n:
        #1,3,...,2n-1
        #for j in xrange(1,m,2): #  do j=1,m-2,2
        m=2*j-1
        mpfr_mul_si(tmp,tmp,m,rnd_re) #=tmp*RF(2*j-1)
    #      tmp=tmp*sqrt(mppic)/mp_two**(mp_half*(mm-mp_one))
    #tmp2 = RF.pi().sqrt()/RF(2)**nn
    mpfr_sqpi_erfc_sqx(tmp2,x)
    #mpfr_const_pi(tmp2,rnd_re)
    #mpfr_sqrt(tmp2,tmp2,rnd_re)
    mpfr_div_si(tmp2,tmp2,2**n,rnd_re)
    #mpfr_mul(tmp,tmp,tmp2,rnd_re)
    mpfr_sqrt(sqx,x,rnd_re)
    #mpfr_erfc(tmp2,sqx,rnd_re)
    #tmp2 = erfc(x.sqrt())
    #tmp=tmp*erfc(x.sqrt())
    mpfr_mul(tmp,tmp,tmp2,rnd_re)
    mpfr_neg(tmp2,x,rnd_re)
    mpfr_exp(tmp2,tmp2,rnd_re)
    mpfr_mul(tmp2,tmp2,sqx,rnd_re)
    #tmp2=(-x).exp()*x.sqrt()
    if n % 2 ==0: #( mod(n,2).eq.0) then
        mpfr_neg(tmp2,tmp2,rnd_re) #-(-x).exp()*x.sqrt()
    mpfr_mul(tmp2,tmp2,summa,rnd_re)
    #tmp2=tmp2*summa
    mpfr_add(tmp2,tmp2,tmp,rnd_re)
    #tmp=tmp+tmp2
    mpfr_set(res,tmp2,rnd_re)
    mpfr_clear(half_m_n)
    mpfr_clear(summa)
    mpfr_clear(tmp)
    mpfr_clear(sqx)
    mpfr_clear(tmp2)
    return ok
    #print 'incG(',n,'+0.5,',x,'=',tmp
    #return tmp2

#  !!
# cpdef RealNumber incgamma_nhint(int n,RealNumber x,int verbose=0):
#     cdef RealNumber res
#     cdef int ok = 1
#     res = x.parent()(1)
#     ok = incgamma_nhint_c(res.value,n,x.value,verbose)
#     if ok == 0:
#         return res
#     else:
#         raise ArithmeticError,"Did not achieve sought precision!"        

#!!  incgamma(-n+1/2,x)
#      !! for integer n>0 and real x>0
cdef int incgamma_nhint_c(mpfr_t res,int n,mpfr_t x,int verbose=0) nogil:
    cdef int j
    cdef mpfr_prec_t prec
    cdef mpfr_t tmp,summa,tmp2,tmp3,half,half_m_n,lnx,tmpr
    cdef int ok = 0
    #assert n>=0
    if n < 0:
        return -1
    #write(*,*) 'ERROR: incgamma_nhint for n>0! n=',n
    #RF=x._parent
    prec = mpfr_get_prec(x)
    mpfr_init2(tmpr,prec)
    mpfr_init2(lnx,prec)
    mpfr_init2(tmp,prec)
    mpfr_init2(tmp2,prec)
    mpfr_init2(tmp3,prec)
    mpfr_init2(half,prec)
    mpfr_init2(half_m_n,prec)
    mpfr_init2(summa,prec)

    #mpfr_set_ui(tmp,1,rnd_re)
    mpfr_set_si(summa,0,rnd_re)
    mpfr_set_ui(half,1,rnd_re)
    mpfr_div_ui(half,half,2,rnd_re)
    mpfr_sub_si(half_m_n,half,n,rnd_re)
    mpfr_pow(tmp,x,half_m_n,rnd_re)

    # half_m_n = 1/2 - n
    # !! do the sum first
    for j from 0 <= j <=n-1:
        _mppochammer_mpfr(tmp2,half_m_n,j+1)
        mpfr_div(tmp2,tmp,tmp2,rnd_re)
        mpfr_add(summa,summa,tmp2,rnd_re)
        mpfr_mul(tmp,tmp,x,rnd_re)
    #tmp2=mppochammer(RF(0.5),n)
    _mppochammer_mpfr(tmp2,half,n)
    mpfr_sqpi_erfc_sqx(tmp3,x)
    mpfr_div(tmp2,tmp3,tmp2,rnd_re)
    #tmp2=RF.pi().sqrt()/tmp2*erfc(x.sqrt())
    if n%2==1:
        mpfr_neg(tmp2,tmp2,rnd_re)
    #mpfr_log(lnx,x,rnd_re)
    #mpfr_mul(tmp,lnx,half_m_n,rnd_re)
    mpfr_neg(tmp,x,rnd_re)
    mpfr_exp(tmp3,tmp,rnd_re)
    # tmp3 = x**(1/2-n)*exp(-x)
    mpfr_neg(summa,summa,rnd_re)
    mpfr_mul(tmp3,tmp3,summa,rnd_re)
    mpfr_add(tmp2,tmp2,tmp3,rnd_re)
    #print 'incG(',n,'+0.5,',x,'=',tmp
    mpfr_set(res,tmp2,rnd_re)
    mpfr_clear(tmp)
    mpfr_clear(summa)
    mpfr_clear(tmp2)
    mpfr_clear(tmp3)
    mpfr_clear(half)
    mpfr_clear(lnx)
    mpfr_clear(half_m_n)
    return ok
    #return tmp2


def incgamma_nhint_test(n,x,verbose=0):
    #RF=RealField(53)
    RF = RealField(x.parent().prec())
    #RF =
    #prec = RF.prec()
    tmp = RF(1)
    summa=RF(0)
    half_m_n=RF(0.5)-RF(n)
    b=RF(1)/RF(0.5-n)
    summa=b
    for j from 1 <= j <= n-1:
        #tmp2 = mppochammer(half_m_n,j+1)
        #term = tmp/tmp2
        b = b*(x/RF(0.5-n+j))
        summa = summa + b
        #if verbose>0:
        #    print "term=",term
        #tmp=tmp*x
    if verbose>0:
        print "sum=",summa

    RF = RealField(103)
    tmp2 = pochammer(RF(0.5),n)
    tmp = RF.pi().sqrt()*(RF(x).sqrt()).erfc()/tmp2
    RF = RealField(x.parent().prec())
    if verbose>0:
        print "tmp=",tmp
    if (n%2) == 1:
        tmp=-tmp
    #tmp3 = (x**half_m_n)*(-x).exp()
    tmp3 = ( half_m_n*(x.log())-x).exp()
    tmp2 = tmp3*summa
    if verbose>0:
        print "tmp=",tmp
        print "tmp2=",tmp2

    return tmp - tmp2

cpdef RealNumber pochammer(RealNumber a,int k):
    cdef RealNumber res
    res = a._parent(0)
    _mppochammer_mpfr(res.value,a.value, k)
    return res

cdef void _mppochammer_mpfr(mpfr_t res, mpfr_t a,int k) nogil:
    r"""
    res should be initialized outside
    """
    cdef int j
    cdef mpfr_prec_t prec
    cdef mpfr_t tmp
    prec = mpfr_get_prec(a)
    #mpfr_init2(res,prec)
    mpfr_init2(tmp,prec)
    if k == 0:
        mpfr_set_ui(res,1,rnd_re)
    else:
        mpfr_set(res,a,rnd_re)
        for j from 1<=j<=k-1:
            mpfr_add_ui(tmp,a,j,rnd_re)
            mpfr_mul(res,res,tmp,rnd_re)
            # res=res*(a+<double>j)
    mpfr_clear(tmp)
