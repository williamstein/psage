# -*- coding=utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2010 Fredrik Strömberg <stroemberg@mathematik.tu-darmstadt.de>
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
Cython algorithms for scalar-valued Poincare series.

TODO: add proven error estimates

"""
include 'interrupt.pxi'
from psage.rings.mp_cimports cimport *

from sage.rings.complex_mpc import MPComplexField
from sage.rings.real_mpfr import RealField
#from sage.all import Bessel_J
import mpmath
import cython

cpdef Bplus_triv(int m,int k,int n,int N,int prec=53, int NMAX=50,int verbose=0):
    r"""
    Get the coefficient b^+(m,n,N) of Q^+(m,k,N) the m-th weight k Poincare series of Gammm0(N) with m<=-1. Sum up to a bound of NMAX

    """
    cdef MPComplexNumber res,tmp
    cdef mpc_t term,kval
    cdef mpfr_t factor,factor2
    CF =MPComplexField(prec)
    res=CF(0) #; kval=CF(0)
    tmp=CF(0)
    cdef int c,cn,cnk
    mpc_init2(term,prec)
    mpc_init2(kval,prec)
    mpfr_init2(factor,prec)
    mpfr_init2(factor2,prec)
    if n==0:
        for c from 1<=c<=NMAX:
            cn=N*c
            Ktriv0(m,n,cn,kval)
            mpc_set(tmp.value,kval,rnd)            
            #print "K,m,n(",cn,")=",tmp
            mpc_div_ui(term,kval,cn**k,rnd)
            #term=kval/cn**2
            mpc_add(res.value,res.value,term,rnd)
            mpc_set(tmp.value,term,rnd)
            #print "term(",cn,")=",tmp

    mpfr_const_pi(factor,rnd_re)
    mpfr_mul_ui(factor,factor,2,rnd_re)
    mpfr_pow_si(factor,factor,k,rnd_re)
    mpfr_mul_ui(factor,factor,(-m)**(k-1),rnd_re)
    mpfr_fac_ui(factor2,k-1,rnd_re)
    mpfr_mul(factor,factor,factor2,rnd_re)
    if k % 2==0:
        if k % 4 == 2:
            mpfr_neg(factor,factor,rnd_re)
        mpc_set_fr(term,factor,rnd)
    else:
        mpc_set_fr(term,factor,rnd)
        if k % 4 == 1:
            mpc_mul_i(term,term,1,rnd)
        else:
            mpc_mul_i(term,term,-1,rnd)
    mpc_mul(res.value,res.value,term,rnd)
    mpfr_clear(factor)
    mpfr_clear(factor2)
    mpc_clear(term)
    mpc_clear(kval)
    return res


cpdef Ktriv(int N,int m, int n,int c,int prec):
    cdef MPComplexNumber res
    CF =MPComplexField(prec)
    res=CF(0)
    if c % N==0:
        Ktriv0(m,n,c,res.value)
    return res

cdef Ktriv0(int m, int n,int c,mpc_t res,int verbose=0):
    r"""
    K(m,n,c) = sum_{d (c)} e((md+n\bar{d})/c)
    """
    cdef int prec
    prec = mpc_get_prec(res)
    CF = MPComplexField(prec)
    mpc_set_si(res,0,rnd)
    #z=CyclotomicField(c).gen()
    cdef mpz_t dd,dinv,cc,arg,mm,tmp
    cdef mpc_t arg1
    mpc_init2(arg1,prec)
    mpz_init(dd);    mpz_init(cc)
    mpz_init(mm);    mpz_init(dinv)
    mpz_init(arg); mpz_init(tmp)
    mpz_set_ui(cc,c)
    mpz_set_si(mm,m)
    cdef int t
    cdef MPComplexNumber term,arg0,tmpc
    #cdef MPComplexField CF
    if c==1:
        mpc_set_ui(res,1,rnd)
        return 
    tmpc=CF(0)
    term=CF(0)
    arg0=CF(0,2)*CF.base_ring().pi() #/CF(c)
    mpc_div_ui(arg0.value,arg0.value,c,rnd)
    #arg1=CF(0)
    #print "arg0=",arg0
    for d from 1<=d<c:
        mpz_set_ui(dd,d)
        t = mpz_invert(dinv,dd,cc)
        if t<>0:
            mpz_mul(arg,dinv,mm)
            mpz_set_si(tmp,n*d)
            mpz_add(arg,arg,tmp)
            mpc_mul_si(arg1,arg0.value,mpz_get_si(arg),rnd)
            #print "arg=",mpz_get_si(arg)
            mpc_exp(arg1,arg1,rnd)
            #mpc_set(tmpc.value,arg1,rnd)
            #print "exp(arg/c)=",tmpc
            mpc_add(res,res,arg1,rnd)
    mpc_clear(arg1)
    #mpc_set(term.value,res,rnd)
    #print "res=",term



cpdef Aplus_triv(int m,int k,int n,int N,int prec=53, int NMAX=50,int verbose=0):
    r"""
    Get the coefficient a^+(m,n,N) of P^+(m,k,N) the m-th weight k Poincare series of Gammm0(N) with m>0. Sum up to a bound of NMAX

    """
    assert n>0
    cdef MPComplexNumber res,tmp
    cdef double eps,epscn
    cdef mpc_t term,kval
    cdef mpfr_t factor,factor2
    cdef RealNumber arg,fourpi,kone,besj,tmpr,konehalf
    if verbose>0:
        print "m,n,k,N=",m,n,k,N
    eps = 2.0**-prec
    CF =MPComplexField(prec)
    RF = RealField(prec)
    fourpi = RF.pi()*RF(4)
    besj=RF(0)
    tmpr=RF(0)
    kone=RF(k-1)
    konehalf=kone/RF(2)
    arg = RF(0)
    res=CF(0) #; kval=CF(0)
    tmp=CF(0)
    mpmath.mp.prec = prec
    cdef int c,cn,cnk
    mpc_init2(term,prec)
    mpc_init2(kval,prec)
    mpfr_init2(factor,prec)
    mpfr_init2(factor2,prec)
    mpfr_set_si(arg.value,m,rnd_re)
    mpfr_mul_si(arg.value,arg.value,n,rnd_re)
    #if verbose>0:
    #    print "sqrt(mn)*4pi={0}".format(arg)
    if m<0:
        mpfr_abs(arg.value,arg.value,rnd_re)
    mpfr_sqrt(arg.value,arg.value,rnd_re)    
    mpfr_mul(arg.value,arg.value,fourpi.value,rnd_re)
    if verbose>0:
        print "sqrt(mn)*4pi={0}".format(arg)
        print "eps=",eps
    for c from 1<=c<=NMAX:
        cn=N*c
        Ktriv0(m,n,cn,kval,verbose-1)
        epscn=eps*cn
        if mpfr_cmp_d(mpc_realref(kval),epscn)<0 and mpfr_cmp_d(mpc_realref(kval),-epscn)>0:
            mpc_set(tmp.value,kval,rnd)
            if mpfr_cmp_d(mpc_imagref(kval),epscn)<0 and mpfr_cmp_d(mpc_imagref(kval),-epscn)>0:            
                if verbose>1:
                    print "Skipping kval({0})={1}".format(cn,tmp)
                continue

        mpc_set(tmp.value,kval,rnd)
        if verbose>1:
            print "K({0},{1},{2})={3}".format(m,n,cn,tmp)
        mpc_div_ui(tmp.value,tmp.value,cn,rnd)
        mpfr_div_si(tmpr.value,arg.value,cn,rnd_re)
        if m<0:
            besj = RF(mpmath.mp.besseli(kone,tmpr))
        else:
            besj = RF(mpmath.mp.besselj(kone,tmpr))
        if verbose>1:
            print "J_{0}({1})={2}".format(kone,tmpr,besj)
        mpc_mul_fr(tmp.value,tmp.value,besj.value,rnd)
        mpc_add(res.value,res.value,tmp.value,rnd)
        if verbose>1:
            #print "J_{0}({1})={2}".format(kone,tmpr,besj)
            print "term(",cn,")=",tmp
    if verbose>0:
        print "res=",res
    
    mpfr_const_pi(factor,rnd_re)
    mpfr_mul_ui(factor,factor,2,rnd_re)
    #mpfr_pow_si(factor,factor,k,rnd_re)
    mpfr_set_si(tmpr.value,n,rnd_re)
    if m<0:
        mpfr_div_si(tmpr.value,tmpr.value,-m,rnd_re)
    else:
        mpfr_div_si(tmpr.value,tmpr.value,m,rnd_re)
    mpfr_pow(tmpr.value,tmpr.value,konehalf.value,rnd_re)
    mpfr_mul(factor,factor,tmpr.value,rnd_re)
    if k % 2==0:
        if k % 4 == 2:
            mpfr_neg(factor,factor,rnd_re)
        mpc_set_fr(term,factor,rnd)
    else:
        mpc_set_fr(term,factor,rnd)
        if k % 4 == 1:
            mpc_mul_i(term,term,1,rnd)
        else:
            mpc_mul_i(term,term,-1,rnd)
    mpc_mul(res.value,res.value,term,rnd)
    mpfr_clear(factor)
    mpfr_clear(factor2)
    mpc_clear(term)
    mpc_clear(kval)
    return res
