# -*- coding=utf-8 -*-
#*****************************************************************************
#ö  Copyright (C) 2010  Fredrik Strömberg <stroemberg@mathematik.tu-darmstadt.de>
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


#include "sage/ext/cdefs.pxi"
#include "sage/ext/interrupt.pxi"  # ctrl-c interrupt block support
#include "sage/ext/stdsage.pxi"  # ctrl-c interrupt block support
#include "sage/ext/gmp.pxi"
#include "sage/rings/mpc.pxi"


## For multiprecision support
from sage.libs.mpfr cimport *
cdef mpc_rnd_t rnd
cdef mpfr_rnd_t rnd_re
#rnd = MPC_RNDNN
rnd_re = GMP_RNDN
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
#### The incomplete Gamma function of integer parameter


cpdef RealNumber incgamma_int(int n,RealNumber x,int verbose=0):
    cdef RealNumber res
    res = x.parent()(0)
    if n>0:
        incgamma_pint_c(res.value,n,x.value,verbose)
    else:
        incgamma_nint_c(res.value,n,x.value,verbose)
    return res
## Integer a=n>=0
## Very crude implementations...

#cdef RealNumber incgamma_pint(int n,RealNumber x):
cdef incgamma_pint_c(mpfr_t res, int n,mpfr_t x,int verbose=0):
    r"""
    Incomplete Gamma, Gamma(n,x) with n positive integer.
    """
    cdef mpfr_t tmp,tmp2,tmp_fak
    cdef int j,prec
    cdef RealField_class RF
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
    #return res

## Gamma(n,x) with n<0 and real x>0
#cpdef RealNumber incgamma_nint(int n,RealNumber x):
cdef incgamma_nint_c(mpfr_t res, int n,mpfr_t x,int verbose=0):
    r"""
    Incomplete Gamma, Gamma(n,x) with n negative integer.

    We use Gamma(0,x)=-Ei(-x) for x>0
    """
    cdef RealNumber summa,tmp,tmp2,fak
    cdef RealField_class RF
    cdef int nn,k,prec
    if n > 0:
        raise ValueError, "ERROR: n<=0 in this function!"
    if verbose>0:
        print "In incgamma_nint_c with n={0}".format(n)
    prec = mpfr_get_prec(x)
    RF = RealField(prec)
    summa=RF(0); tmp=RF(1); tmp2=RF(1); fak=RF(1)
    if mpfr_cmp_si(x,0)<0:
        mpfr_set(tmp.value,x,rnd_re)
        raise ValueError, "ERROR: we assume that x>0. Got:{0}".format(tmp)
    nn=-n
    if verbose>0:
        print "-n=",nn
        print "-x=",tmp
    if n==0:
        ## Gamma(0,x)=-Ei(-x)
        mpfr_neg(tmp.value,x,rnd_re)
        if verbose>0:
            print "n=0"
        mpfr_neg(res,tmp.value,rnd_re)
        #mpfr_sub(res,res,tmp.value,rnd_re)
        return
    Ei_ml_c(tmp.value,tmp.value) # = Ei(x) - ln|x|
    #!! Make the Ei -sum here
    #tmp=Ei_ml_c(-x) = Ei(-x) - ln|x|
    #
    if verbose>0:
        print "Ei(-x)-ln(|x|)=",tmp
    raise NotImplementedError,"Doesn't work right now..."
    #tmp2=(x**n*RF(nn))**-1
    mpfr_pow_si(tmp2.value,x,n,rnd_re)
    mpfr_mul_si(tmp2.value,tmp2.value,nn,rnd_re)
    mpfr_pow_si(tmp2.value,tmp2.value,-1,rnd_re)
    #  write(*,*) 'here2!'
    mpfr_set(summa.value,tmp2.value,rnd_re)
    if verbose>0:
        print "Summa1=",summa
    for k in range(1,nn+1): #from 1 <= k <= -n:
        #write(*,*) 'here! k=',k
        mpfr_set_si(fak.value,k-1-n,rnd_re)
        # fak=RF(-n+k-1)
        #tmp2=tmp2*x/fak
        mpfr_mul(tmp2.value,tmp2.value,x,rnd_re)
        mpfr_div(tmp2.value,tmp2.value,fak.value,rnd_re)
        #summa=summa+tmp2
        mpfr_add(summa.value,summa.value,tmp2.value,rnd_re)
        # write(*,*) 'summa:',dble(summa)
    #fak=RF( (-1)**(n-1))
    if verbose>0:
        print "Summa2=",summa
    if n-1 % 2 == 0:
        mpfr_set_si(fak.value,1,rnd_re)
    else:
        mpfr_set_si(fak.value,-1,rnd_re)
    mpfr_fac_ui(tmp2.value,nn,rnd_re)  #=RF.factorial(nn) #gamma(n+1)
    mpfr_div(res,fak.value,tmp2.value,rnd_re)
    mpfr_mul(res,res,tmp.value,rnd_re)
    if verbose>0:
        print "fak=",fak
        print "res=",tmp
    #tmp2=tmp*fak/tmp2
    #print 'Gamma(',n,',',x,')=',tmp2
    #return tmp2

cpdef RealNumber ei(RealNumber x, int verbose=0):
    cdef int wp = x.prec() + 20
    cdef int xabs= abs(x.integer_part())
    cdef RealField_class RF_orig = x.parent()
    cdef RealField_class RF
    cdef int rh = RF_orig(wp*0.693).integer_part() + 10
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

cpdef RealNumber ei_asymp(RealNumber x, int verbose=0):
    r"""
    Compute the exponential integral of x via asymptotic formula
    """
    cdef RealNumber res
    res = x.parent()(0)
    ei_asymp_c(res.value,x.value,verbose)
    return res

cdef ei_asymp_c(mpfr_t res, mpfr_t x, int verbose=0):
    r"""
    Compute the exponential integral of x via asymptotic formula.
    Only valid for x <= -40.
    """
    cdef RealNumber tmp,tmp2,summa,r
    cdef double eps
    cdef int k,prec
    cdef RealField_class RF
    prec = mpfr_get_prec(x)
    RF = RealField(prec)
    tmp=RF(1); summa=RF(1); r=RF(1); tmp2=RF(0)
    eps = 2.**-(prec+1)
    if verbose>0:
        print "prec = ", prec
    mpfr_set(summa.value,tmp.value,rnd_re)
    mpfr_div(r.value,tmp.value,x,rnd_re)
    mpfr_exp(tmp2.value,x,rnd_re)
    mpfr_mul(tmp2.value,tmp2.value,r.value,rnd_re)
    eps=abs(eps/tmp2)
    if verbose>0:
        print "r = 1/x = ", r
        print "eps = ", eps
        print "exp(x)/x = ", tmp2
        print "tmp = ", tmp
    for k in range(1,nmax): #from 1 <= k <= nmax:
        mpfr_mul_ui(tmp.value,tmp.value,k,rnd_re)
        mpfr_mul(tmp.value,tmp.value,r.value,rnd_re)
        mpfr_add(summa.value,summa.value,tmp.value,rnd_re)
        if verbose>0:
            print "k = ", k
            print "r = ", r
            print "tmp = ", tmp
            print "summa = ", summa
        if mpfr_cmp_d(tmp.value,eps)<0 and mpfr_cmp_d(tmp.value,-eps)>0:
            # if  abs(tmp) < eps:
            if verbose>0:
                print 'break at k=', k
            break
    if k>= nmax:
        mpfr_set(summa.value,x,rnd_re)
        raise ArithmeticError,"k>nmax! in Ei(x)!, error of order {0} for x={1}".format(tmp,summa)
    mpfr_mul(res,summa.value,tmp2.value,rnd_re)
    return

    #!! Ei(x)-ln(|x|)  for real x
cpdef RealNumber Ei_ml(RealNumber x):
    r"""
    Compute the exponential integral of x  - ln|x|
    """
    cdef RealNumber res
    res = x.parent()(0)
    Ei_ml_c(res.value,x.value)
    return res

cpdef RealNumber ei_taylor(RealNumber x, int verbose=0):
    r"""
    Compute the exponential integral of x  - ln|x|
    """
    cdef RealNumber res
    cdef mpfr_t lnx
    cdef int prec=x.prec()
    mpfr_init2(lnx,prec)
    res = x.parent()(0)
    Ei_ml_c(res.value,x.value)
    if verbose>0:
        print "Ei(-x)-ln|x| = ", res
    if x<0:
        x=-x
    mpfr_log(lnx,x.value,rnd_re)
    if verbose>0:
        print "ln|x| = ", x
    mpfr_add(res.value,res.value,lnx,rnd_re)
    return res

cdef Ei_ml_c(mpfr_t res,mpfr_t x):
    r"""
    Compute the exponential integral of x  - ln|x|
    Only good for x>-40.
    """
    cdef RealNumber tmp,summa
    cdef double eps
    cdef int k,prec
    cdef RealField_class RF
    cdef mpfr_t ec
    prec=mpfr_get_prec(x)
    mpfr_init2(ec,prec)
    RF = RealField(prec)
    tmp=RF(1); summa=RF(0)
    eps = 2.0**-(prec+1)
    #call set_eulergamma()
    mpfr_set(tmp.value,x,rnd_re)
    #summa=tmp
    mpfr_set(summa.value,tmp.value,rnd_re)
    for k in range(2,nmax+1): #from 2 <= k <= nmax:
        #kk=RF(k)
        ##tmp=tmp*(kk-RF(1))/kk**2
        mpfr_mul_ui(tmp.value,tmp.value,k-1,rnd_re)
        mpfr_div_ui(tmp.value,tmp.value,k*k,rnd_re)

        mpfr_mul(tmp.value,tmp.value,x,rnd_re)
        mpfr_add(summa.value,summa.value,tmp.value,rnd_re)
        if mpfr_cmp_d(tmp.value,eps)<0 and mpfr_cmp_d(tmp.value,-eps)>0:
            # if  abs(tmp) < eps:
            break
    if k>= nmax:
        mpfr_set(summa.value,x,rnd_re)
        raise ArithmeticError,"k>nmax! in Ei(x)!, error of order {0} for x={1}".format(tmp,summa)
    mpfr_const_euler(ec,rnd_re)
    mpfr_add(summa.value,summa.value,ec,rnd_re)
    mpfr_clear(ec)
    #summa=summa+RF.euler_constant()
    #print 'Ei(',x,')=',summa
    mpfr_set(res,summa.value,rnd_re)
    #return summa


#      !!  incgamma(n+1/2,x) integer n and real x>0
cpdef RealNumber incgamma_hint(int n,RealNumber x,int verbose=0):
    r"""
    Incomplete Gamma function of half-integer parameter, Gamma(n+1/2,x).

    INPUT:
    -`n` -- integer
    -`x` -- real number
    -`verbose` -- integer

    """
    cdef RealNumber res
    res = x.parent()(0)
    incgamma_hint_c(res.value,n,x.value)
    return res
    #cdef RealField_class RF
    #if n > 0:
    #    return incgamma_phint(n,x)
    #elif n<0:
    #    return incgamma_nhint(-n,x)
    #else:
    #   RF=x._parent
    #    return RF.pi().sqrt()*erfc(x.sqrt())

cdef incgamma_hint_c(mpfr_t res,int n,mpfr_t x,int verbose=0):
    cdef RealField_class RF
    cdef mpfr_t sqpi,sqx
    cdef int prec
    if n > 0:
        #return
        incgamma_phint_c(res,n,x,verbose)
    elif n<0:
        incgamma_nhint_c(res,-n,x,verbose)
    else:
        prec = mpfr_get_prec(x)
        mpfr_init2(sqpi,prec)
        mpfr_init2(sqx,prec)
        mpfr_const_pi(sqpi,rnd_re)
        mpfr_sqrt(sqpi,sqpi,rnd_re)
        mpfr_sqrt(sqx,x,rnd_re)
        mpfr_erfc(res,sqx,rnd_re)
        mpfr_mul(res,res,sqpi,rnd_re)


### sqrt(pi)* erfc(sqrt(x))  is used at several places

cdef mpfr_sqpi_erfc_sqx(mpfr_t res,mpfr_t x,int verbose=0):
    cdef int prec = mpfr_get_prec(x)
    cdef mpfr_t sqpi
    mpfr_init2(sqpi,prec)
    mpfr_const_pi(sqpi,rnd_re)
    mpfr_sqrt(sqpi,sqpi,rnd_re)
    mpfr_sqrt(res,x,rnd_re)
    mpfr_erfc(res,res,rnd_re)
    mpfr_mul(res,res,sqpi,rnd_re)
    mpfr_clear(sqpi)

#    !!  incgamma(n+1/2,x)
# !! for integer n>=0 and real x>0
cpdef RealNumber incgamma_phint(int n,RealNumber x,int verbose=0):
    cdef RealNumber res
    res = x.parent()(1)
    incgamma_phint_c(res.value,n,x.value,verbose)
    return res

cdef void incgamma_phint_c(mpfr_t res, int n,mpfr_t x,int verbose=0):
    cdef RealNumber jj,nn,mm
    cdef RealField_class RF
    cdef int j,m,prec
    cdef mpfr_t half_m_n,summa,tmp,term,tmp2,sqx
    assert n>=0
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
    #print 'incG(',n,'+0.5,',x,'=',tmp
    #return tmp2

#  !!
cpdef RealNumber incgamma_nhint(int n,RealNumber x,int verbose=0):
    cdef RealNumber res
    res = x.parent()(1)
    incgamma_nhint_c(res.value,n,x.value,verbose)
    return res

#!!  incgamma(-n+1/2,x)
#      !! for integer n>0 and real x>0
cdef void incgamma_nhint_c(mpfr_t res,int n,mpfr_t x,int verbose=0):
    cdef RealNumber tmpr
    cdef int j,prec
    cdef RealField_class RF
    cdef mpfr_t tmp,summa,tmp2,tmp3,half,half_m_n,lnx
    assert n>=0
    #write(*,*) 'ERROR: incgamma_nhint for n>0! n=',n
    #RF=x._parent
    prec = mpfr_get_prec(x)
    if verbose>0:
        tmpr=RealField(prec)(1)
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
    if verbose>0:
        mpfr_set(tmpr.value,summa,rnd_re)
        print "sum=",tmpr
    _mppochammer_mpfr(tmp2,half,n)
    mpfr_sqpi_erfc_sqx(tmp3,x)
    mpfr_div(tmp2,tmp3,tmp2,rnd_re)
    if verbose>0:
        mpfr_set(tmpr.value,tmp2,rnd_re)
        print "tmp=",tmpr
    #tmp2=RF.pi().sqrt()/tmp2*erfc(x.sqrt())
    if n%2==1:
        mpfr_neg(tmp2,tmp2,rnd_re)
    #mpfr_log(lnx,x,rnd_re)
    #mpfr_mul(tmp,lnx,half_m_n,rnd_re)
    mpfr_neg(tmp,x,rnd_re)
    mpfr_exp(tmp3,tmp,rnd_re)
    # tmp3 = x**(1/2-n)*exp(-x)
    if verbose>0:
        mpfr_set(tmpr.value,tmp2,rnd_re)
        print "tmp2=",tmpr
    mpfr_neg(summa,summa,rnd_re)
    mpfr_mul(tmp3,tmp3,summa,rnd_re)
    if verbose>0:
        mpfr_set(tmpr.value,tmp3,rnd_re)
        print "tmp2=",tmpr

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

cdef void _mppochammer_mpfr(mpfr_t res, mpfr_t a,int k):
    r"""
    res should be initialized outside
    """
    cdef int j,prec
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




