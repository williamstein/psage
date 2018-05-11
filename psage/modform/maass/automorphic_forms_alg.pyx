# cython: profile=False
# -*- coding: utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2010 Fredrik Strömberg <stroemberg@mathematik.tu-darmstadt.de>,
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
Cython algorithms for Harmoic weak Maass forms.
Used by routines in atomorphic_forms.py

"""
from libc.stdint cimport uint64_t

from psage.rings.mpfr_nogil cimport *
from cysignals.memory cimport sig_free,sig_malloc,check_allocarray
from cysignals.signals cimport sig_on,sig_off

import logging
log = logging.getLogger(__name__)
cdef extern from "stdio.h":
    cdef extern void printf(char *fmt,...) nogil
    
#from sage.libs.mpfr cimport *
cdef mpc_rnd_t rnd
cdef mp_rnd_t rnd_re
import cython
from cython.parallel cimport parallel, prange
rnd = MPC_RNDNN
rnd_re = MPFR_RNDN
from sage.rings.complex_mpc cimport MPComplexNumber,MPComplexField_class
from sage.rings.complex_mpc import MPComplexField
from sage.rings.real_mpfr cimport RealNumber,RealField_class
from sage.rings.real_mpfr import RealField
from sage.rings.complex_number cimport ComplexNumber
import sage.structure.element
#from sage.libs.pari.gen cimport GEN
#cimport sage.structure.element
from psage.modules.vector_real_mpfr_dense cimport Vector_real_mpfr_dense
from psage.functions.inc_gamma cimport incgamma_hint_c,incgamma_nint_c,incgamma_pint_c
from psage.rings.mpc_extras cimport _mpc_div,_mpc_mul,_mpc_set,_mpc_sub,_mpc_mul_fr,_mpc_add

from sage.structure.element cimport Element, ModuleElement, RingElement
from sage.all import MatrixSpace,vector,ComplexField
cdef extern from "math.h" nogil:
    cdef double fabs(double)
    cdef double fmax(double,double)
    cdef int ceil(double) 

from sage.functions.all import ceil as pceil
from sage.functions.all import floor
#load "
from maass_forms_alg import *
#from maass_forms_alg cimport *

from psage.modform.arithgroup.mysubgroups_alg import normalize_point_to_cusp_mpfr as normalize_point_to_cusp
from psage.modform.arithgroup.mysubgroups_alg cimport _apply_sl2z_map_mpfr
from psage.modform.arithgroup.mysubgroups_alg import apply_sl2z_map,pullback_to_psl2z_mat #,pullback_to_G
#from maass_forms_alg import smallest_inf_norm
from sage.modular.arithgroup.congroup_sl2z import SL2Z
import mpmath
from sage.all import SageObject,copy,deepcopy
from pullback_algorithms import pullback_pts_mpc,j_fak_int_mpfr,pullback_pts_mp
from psage.matrix.matrix_complex_dense cimport Matrix_complex_dense
from psage.modules.vector_complex_dense cimport Vector_complex_dense
from pullback_algorithms cimport pullback_pts_mpc_new_c
from pullback_algorithms cimport pullback_pts_mpc_new_c_sym

cimport openmp

openmp.omp_set_dynamic(1)#openmp.omp_set_num_threads(2)

def get_Y_and_M_for_hwmf(G,PP,weight,ndig):
    r"""
    Find a good Y and M for computing coefficents with precison 10^-ndig

    """
     # generalized_level
    Y0=G.minimal_height()*mpmath.mpf(0.98)
    # then get corresponding M
    #M0=get_M_for_hwmf(Y0,Kmax,Cmax,weight,ndig)
    M0=get_M_for_hwmf(Y0,weight,ndig,PP)
    return [Y0,M0]


def get_M_for_holom(Y,weight,prec=10):
    r""" Get the correct M needed to approximate a  holomorphic
    modular form. 
    """
    # to get the range of M we check we make some crude
    # estimates to begin with
    if weight <= 0:
        raise NotImplementedError,"We have not implemented holomorphic forms of weight <= 0!"
    elif(abs(weight-mpmath.mp.mpf(1))>mpmath.mp.eps()):
        t1=mpmath.mp.pi()*mpmath.mp.mpf(2)*mpmath.mpf(Y)/mpmath.mp.mpf(weight-1)
    else:
        t1=mpmath.mp.mpf(0.2)
        #print "t1=",t1
    x=mpmath.mp.pi()*mpmath.mpf(Y)
    q=mpmath.exp(-x)
    t0=(prec-mpmath.mp.ln(mpmath.mp.mpf(1)-q))/x
    mstart=max(ceil(t0)+1,10) # A crude estimate of good M
    #print "mstart=",mstart
    eps=mpmath.mp.mpf(0.5)*mpmath.power(10,-prec)
    for M in range(mstart,3*mstart):
        t0=mpmath.mp.ln(M)/mpmath.mp.mpf(M)
        #print "t0(",M,")=",t0
        if(t0 > t1):
            continue
        t2=err_est_holo(Y,M,weight)
        #print "t2(",M,")=",t2
        if(t2<eps):
            #print "t2:%s < eps:%s for M:%s" %(t2,eps,M)
            break
    if(M>=3*mstart-1):
        t2=err_est_holo(Y,M,weight)
        if(t2>eps):
            raise ArithmeticError,"Could not find good M for Y:%s!" % Y
    #print "Got M=",M
    return M
    
def err_est_holo(Y,M,weight):
    r"""
    The theoretical error estimate for the truncation error.
    """
    arg=mpmath.mp.pi()*mpmath.mp.mpf(2)*mpmath.mp.mpf(Y)
    mm=mpmath.mp.mpf(M+1)
    wt=mpmath.mpf(weight+1)/mpmath.mp.mpf(2)
    res=mpmath.mp.power(arg,-wt)
    res=2*res*mpmath.mp.gammainc(wt,arg*mm)
    return res
    
cpdef get_M_for_hwmf(Y,weight,ndig,PP={},verbose=0):
    r""" Computes truncation point for Harmonic Maass waveform.
    """
    ## Use low precision
    dold=mpmath.mp.dps
    #mpmath.mp.dps=int(mpmath.ceil(abs(mpmath.log10(eps))))+5
    mpmath.mp.dps=ndig+5
    twopi=2*mpmath.pi()
    twopiy=twopi*mpmath.mpf(Y)
    # an extra two for the accumulation of errors
    eps=mpmath.mpf(10)**mpmath.mpf(-ndig)
    minm=mpmath.ceil(10)
    do_hint = 0;do_int = 0

    if hasattr(Y,"parent"):        
        RF = RealField(Y.parent().prec())
    else:
        RF = RealField(mpmath.mp.prec)
    YY = RF(Y); 
    if floor(weight)==ceil(weight):
        do_int = 1
        wt = int(weight)-1 
    elif floor(weight-0.5)==ceil(weight-0.5):
        do_hint = 1
        wt = int(weight-0.5) 
    if(len(PP.values())==0):
        Cmax=1
        Kmax=0        
    else:
        Kmax=0        
        Cmax=RF(max(PP.values()))
        for (c,l) in PP.keys():
            if(abs(l)>Kmax):
                Kmax=abs(l)
    rweight = RF(weight)
    reps = RF(eps)
    if verbose>0:
        print "Kmax,cmax=",Kmax,Cmax
        print "do_int=",do_int
        print "do_hint=",do_hint
        print "eps=",eps
    try:
        if do_int==1:
            for m in range(minm,10000,5):
                errest1=err_est_hwmf_pos(YY,m,rweight,Kmax,Cmax)
                #if verbose>0:
                #    print "e1=",errest1
                errest2=err_est_hwmf_neg_int(YY,m,wt,Kmax,Cmax)
                #if verbose>0:
                #     print "e2=",errest2
                if max(abs(errest1),abs(errest2))<reps:
                    raise StopIteration()
        elif do_hint==1:
            for m in range(minm,10000):
                errest1=err_est_hwmf_pos(YY,m,rweight,Kmax,Cmax)
                errest2=err_est_hwmf_neg_hint(YY,m,wt,Kmax,Cmax)
                #if verbose>0:
                #     print "e1=",errest1
                #     print "e2=",errest2
                if(max(errest1,errest2)<eps):
                    raise StopIteration()
        else:
            for m in range(minm,10000):
                errest1=err_est_hwmf_pos_mpmath(YY,m,weight,Kmax,Cmax)
                errest2=err_est_hwmf_mpmath(YY,m,weight,Kmax,Cmax)
                if(max(errest1,errest2)<eps):
                    raise StopIteration()
    except StopIteration:
        if verbose>0:
            print "er +=",errest1
            print "er -=",errest2
            print "m=",m
        mpmath.mp.dps=dold
    return m
    

cpdef err_est_hwmf_pos(RealNumber Y,int M,RealNumber k,int K0,RealNumber K1):
    r"""
    Estimate the positive part of the truncated series

    CAVEATS:
    - we assume the Ramanujan bound for the coefficients $a(n)$, i.e. $|a(n)|\le2$.
    - error estimate holds for 2piMY >> R

    INPUT:

    - ``Y`` -- real > 0
    - ``M``  -- integer <>0
    - ``k``  -- real
    - ``K0``  -- real (-smallest degree of princ. part)
    - ``K1``  -- real (largest coeff. of princ. part)

    OUTPUT:

    - ``r``  -- error estimate:

    EXAMPLES::


        sage: err_est_hwmf_pos(0.85,20,0.5,1)
        mpf('6.0466581020500474e-23')
        sage: err_est_hwmf_pos(0.85,20,0.5,2) 
        mpf('2.5677901894222784e-13')
        sage: err_est_hwmf_pos(0.5,20,0.5,1) 
        mpf('4.3209202868264023e-6')
        sage: err_est_hwmf_pos(0.85,100,0.5,1)
        mpf('2.4812998619046392e-178')
        sage: err_est_hwmf_pos(0.85,20,2.5,2)
        mpf('3.8398588394365332e-14')


    """
    cdef RealNumber YY,res
    cdef int prec = Y.parent().prec()
    #cdef mpfr_t pi,f1,sqM,sqK,B,Bsq,sqrt2pi,sqpiY,f2,f3
    cdef mpfr_t f1,sqM,sqK,B,Bsq,sqpiY,f2,f3
    cdef mpfr_t pi[1]
    cdef mpfr_t sqrt2pi[1]
    
    #mpfr_init2(t,prec);
    #mpfr_init2(c,prec)
    mpfr_init2(pi[0],prec); mpfr_init2(sqrt2pi[0],prec)
    mpfr_init2(f1,prec); mpfr_init2(sqM,prec)
    mpfr_init2(sqK,prec); mpfr_init2(B,prec)
    mpfr_init2(Bsq,prec);mpfr_init2(f2,prec)
    mpfr_init2(sqpiY,prec);mpfr_init2(f3,prec)
    RF = RealField(prec)
    YY = RF(Y); res = RF(0)
    #mpfr_set_ui(c,1,rnd_re) ## Should be some "true" bound
    mpfr_const_pi(pi[0],rnd_re)
    mpfr_mul_ui(sqrt2pi[0],pi[0],2,rnd_re)
    mpfr_sqrt(sqrt2pi[0],sqrt2pi[0],rnd_re)
    mpfr_mul(sqpiY,pi[0],YY.value,rnd_re)
    mpfr_sqrt(sqpiY,sqpiY,rnd_re)
    #mpfr_set_pi=mpmath.mp.pi()
    mpfr_mul_si(f1,pi[0],2,rnd_re)
    mpfr_mul_si(f1,f1,K0,rnd_re)
    mpfr_div(f1,f1,YY.value,rnd_re)
    mpfr_exp(f1,f2,rnd_re)
    #f1=mpmath.mp.exp(2.0*pi*K0/Y)
    mpfr_sqrt_ui(sqM,M,rnd_re)
    #sqM=mpmath.mp.sqrt(mpmath.mp.mpf(M))
    mpfr_set_si(sqK,K0,rnd_re)
    mpfr_sqrt(sqK,sqK,rnd_re)
    #sqK=mpmath.mp.sqrt(mpmath.mp.mpf(K0))
    mpfr_div(B,sqK,YY.value,rnd_re)
    mpfr_neg(B,B,rnd_re)
    mpfr_add(B,B,sqM,rnd_re)
    mpfr_mul(B,B,sqrt2pi[0],rnd_re)
    mpfr_mul(Bsq,B,B,rnd_re)
    mpfr_neg(Bsq,Bsq,rnd_re)
    #B=mpmath.mp.sqrt(2*pi)*(sqM-sqK/Y)
    #Bsq=B*B
    mpfr_exp(f2,Bsq,rnd_re)
    mpfr_mul(f2,f2,f2,rnd_re)
    #f2=f1*mpmath.mp.exp(-Bsq)
    if k<=2:
        mpfr_set_d(B,<double>1.5-<double>k,rnd_re)
        mpfr_pow(f3,YY.value,B,rnd_re)
        mpfr_mul(f3,f3,f2,rnd_re)
        #f3=f2*mpmath.mp.power(Y,mpmath.mpf(1.5)-k)
        if K0>0:
            mpfr_set_d(B,<double>(0.5*k-1.0),rnd_re)
            mpfr_set_si(f2,K0,rnd_re)
            mpfr_pow(f2,f2,B,rnd_re)
            mpfr_mul(f3,f3,f2,rnd_re)
            #fak=f3*mpmath.mp.power(K0,mpmath.mpf(0.5)*k-mpmath.mpf(1))
        #else
        #mpfr_set(fak,f3,rnd_re)
        #print "f2=",f2,"f3=",f3,

    else:
        mpfr_set_d(B,<double>(k-2.5),rnd_re)
        mpfr_set_d(Bsq,<double>(2),rnd_re)
        mpfr_pow(f3,Bsq,B,rnd_re)
        mpfr_mul(f3,f3,f2,rnd_re)
        #f3=f2*mpmath.mp.power(mpmath.mpf(2),k-mpmath.mpf(2.5))
        mpfr_set_si(B,k-1,rnd_re)
        mpfr_mul(f3,f3,B,rnd_re)
        mpfr_div_si(f3,f3,2,rnd_re)
        #f4=f3*mpmath.mp.mpf(k-1)/mpmath.mp.mpf(2)
        mpfr_div(B,sqK,YY.value,rnd_re)
        mpfr_neg(B,B,rnd_re)
        mpfr_add(B,B,sqM,rnd_re)
        #B=(sqM-sqK/Y)
        mpfr_set_d(Bsq,<double>(k-3),rnd_re)
        mpfr_pow(B,B,Bsq,rnd_re)
        mpfr_mul(f3,f3,B,rnd_re)
        #f5=f4*mpmath.mp.power(B,k-3)
        mpfr_div(f3,f3,sqpiY,rnd_re)
        #fak=f5/mpmath.sqrt(pi*Y)
    mpfr_set(res.value,f3,rnd_re)
    res = res*RF(K1)
    #mpfr_clear(c)
    mpfr_clear(pi[0]); mpfr_clear(sqrt2pi[0])
    mpfr_clear(f1); mpfr_clear(sqM)
    mpfr_clear(sqK); mpfr_clear(B)
    mpfr_clear(Bsq);mpfr_clear(f2)
    mpfr_clear(sqpiY); mpfr_clear(f3)
    return res
    #return mpmath.mpf(K1)*fak


def err_est_hwmf_pos_mpmath(Y,M,k,K0,K1):
    r"""
    Estimate the positive part of the truncated series

    CAVEATS:
    - we assume the Ramanujan bound for the coefficients $a(n)$, i.e. $|a(n)|\le2$.
    - error estimate holds for 2piMY >> R

    INPUT:

    - ``Y`` -- real > 0
    - ``M``  -- integer <>0
    - ``k``  -- real
    - ``K0``  -- real (-smallest degree of princ. part)
    - ``K1``  -- real (largest coeff. of princ. part)

    OUTPUT:

    - ``r``  -- error estimate:

    EXAMPLES::


        sage: err_est_hwmf_pos(0.85,20,0.5,1)
        mpf('6.0466581020500474e-23')
        sage: err_est_hwmf_pos(0.85,20,0.5,2) 
        mpf('2.5677901894222784e-13')
        sage: err_est_hwmf_pos(0.5,20,0.5,1) 
        mpf('4.3209202868264023e-6')
        sage: err_est_hwmf_pos(0.85,100,0.5,1)
        mpf('2.4812998619046392e-178')
        sage: err_est_hwmf_pos(0.85,20,2.5,2)
        mpf('3.8398588394365332e-14')


    """
    c=mpmath.mp.mpf(1) ## Should be some "true" bound
    pi=mpmath.mp.pi()
    f1=mpmath.mp.exp(2.0*pi*K0/Y)
    sqM=mpmath.mp.sqrt(mpmath.mp.mpf(M))
    sqK=mpmath.mp.sqrt(mpmath.mp.mpf(K0))
    B=mpmath.mp.sqrt(2*pi)*(sqM-sqK/Y)
    Bsq=B*B
    f2=f1*mpmath.mp.exp(-Bsq)
    if(k <=2):
        f3=f2*mpmath.mp.power(Y,mpmath.mpf(1.5)-k)
        if(K0>0):
            fak=f3*mpmath.mp.power(K0,mpmath.mpf(0.5)*k-mpmath.mpf(1))
        else:
            fak=f3
        #print "f2=",f2,"f3=",f3,

    else:
        f3=f2*mpmath.mp.power(mpmath.mpf(2),k-mpmath.mpf(2.5))
        f4=f3*mpmath.mp.mpf(k-1)/mpmath.mp.mpf(2)
        B=(sqM-sqK/Y)
        f5=f4*mpmath.mp.power(B,k-3)
        fak=f5/mpmath.sqrt(pi*Y)
    return mpmath.mpf(K1)*fak

def err_est_hwmf_neg_mpmath(Y,M,k,K0,K1):
    r"""
    Estimate the negative part of the truncated series.

    CAVEATS:
    - we assume the Ramanujan bound for the coefficients $a(n)$, i.e. $|a(n)|\le2$.
    - error estimate holds for 2piMY >> R

    INPUT:

    - ``Y`` -- real > 0
    - ``M``  -- integer <>0
    - ``N``  -- integer <>0
    - ``k``  -- real
    - ``K0``  -- real (-smallest degree of princ. part)
    - ``K1``  -- real (largest coeff. of princ. part)

    OUTPUT:

    - ``r``  -- error estimate

    EXAMPLES::


 


    """
    c=mpmath.mp.mpf(1) ## Should be some "true" bound
    pi=mpmath.mp.pi()
    two=mpmath.mp.mpf(2)
    f1=mpmath.power(two,-k-mpmath.mpf(1))
    f2=f1/pi/Y
    X=mpmath.mpf(2*M)*Y*pi
    fak=f2*mpmath.gammainc(k-1,X)
    #if(k>=0):
    ##    fak=f1*mpmath.power(X,-k)*mpmath.exp(-X)
    #äelse:
    #    fak=f1*mpmath.power(X,-k)*mpmath.exp(-X)
    #fak=f2mpmath.mp.exp(-two*pi*M*y)
    return mpmath.mpf(K1)*fak

cpdef err_est_hwmf_neg_int(RealNumber Y,int M,int km1,int K0,RealNumber K1):
    r"""
    Estimate the negative part of the truncated series.

    CAVEATS:
    - we assume the Ramanujan bound for the coefficients $a(n)$, i.e. $|a(n)|\le2$.
    - error estimate holds for 2piMY >> R

    INPUT:

    - ``Y`` -- real > 0
    - ``M``  -- integer <>0
    - ``N``  -- integer <>0
    - ``k``  -- real
    - ``K0``  -- real (-smallest degree of princ. part)
    - ``K1``  -- real (largest coeff. of princ. part)

    OUTPUT:

    - ``r``  -- error estimate

    EXAMPLES::


 


    """
    cdef mpfr_t pi,f1,X,ig
    cdef RealNumber res
    cdef int prec = Y.parent().prec()
    RF = RealField(prec)
    res = RF(1)
    #c=mpmath.mp.mpf(1) ## Should be some "true" boun
    mpfr_init2(pi,prec);    mpfr_init2(X,prec)
    mpfr_init2(f1,prec);    mpfr_init2(ig,prec)
    mpfr_const_pi(pi,rnd_re)
    #pi=mpmath.mp.pi()
    #two=mpmath.mp.mpf(2)
    #mpfr_add_ui(f1,k,1,rnd_re)
    #mpfr_neg(f1,f1,rnd_re)
    #mpfr_exp2(f1,f1,rnd_re)
    mpfr_set_si_2exp(f1,1,-km1-2,rnd_re)
  
    #f1=mpmath.power(two,-k-mpmath.mpf(1))
    mpfr_div(f1,f1,Y.value,rnd_re)
    mpfr_div(f1,f1,pi,rnd_re)
    #f2=f1/pi/Y
    mpfr_mul(X,Y.value,pi,rnd_re)
    mpfr_mul_si(X,X,2*M,rnd_re)
    #X=mpmath.mpf(2*M)*Y*pi
    cdef int ok = 1
    if km1>0:
        ok = incgamma_pint_c(ig,km1,X)
    else:
        ok = incgamma_nint_c(ig,km1,X)
    mpfr_set(res.value,ig,rnd_re)
    mpfr_mul(res.value,res.value,f1,rnd_re)
    mpfr_mul(res.value,res.value,K1.value,rnd_re)
    mpfr_clear(pi);    mpfr_clear(X)
    mpfr_clear(f1);    mpfr_clear(ig)
    return res
    #fak=f2*mpmath.gammainc(k-1,X)
    #if(k>=0):
    ##    fak=f1*mpmath.power(X,-k)*mpmath.exp(-X)
    #äelse:
    #    fak=f1*mpmath.power(X,-k)*mpmath.exp(-X)
    #fak=f2mpmath.mp.exp(-two*pi*M*y)
    #return mpmath.mpf(K1)*fak

cpdef err_est_hwmf_neg_hint(RealNumber Y,int M,int kmh,int K0,RealNumber K1):
    r"""
    Estimate the negative part of the truncated series.
    km1 = k - 1/2 is an integer  
    CAVEATS:
    - we assume the Ramanujan bound for the coefficients $a(n)$, i.e. $|a(n)|\le2$.
    - error estimate holds for 2piMY >> R

    INPUT:

    - ``Y`` -- real > 0
    - ``M``  -- integer <>0
    - ``N``  -- integer <>0
    - ``k``  -- real
    - ``K0``  -- real (-smallest degree of princ. part)
    - ``K1``  -- real (largest coeff. of princ. part)

    OUTPUT:

    - ``r``  -- error estimate

    EXAMPLES::


 


    """
    cdef mpfr_t pi,f1,X,ig
    cdef int prec = Y.parent().prec()
    cdef RealNumber res
    #c=mpmath.mp.mpf(1) ## Should be some "true" boun
    RF = RealField(prec)
    res = RF(1)
    mpfr_init2(pi,prec);    mpfr_init2(X,prec)
    mpfr_init2(f1,prec);    mpfr_init2(ig,prec)
    mpfr_const_pi(pi,rnd_re)

    mpfr_set_d(f1,kmh+0.5,rnd_re)
    mpfr_neg(f1,f1,rnd_re)
    mpfr_exp2(f1,f1,rnd_re)
    #mpfr_set_si_2exp(f1,1,-k-0.5,rnd_re)
    #f1=mpmath.power(two,-k-mpmath.mpf(1))
    mpfr_div(f1,f1,Y.value,rnd_re)
    mpfr_div(f1,f1,pi,rnd_re)
    #f2=f1/pi/Y
    mpfr_mul(X,Y.value,pi,rnd_re)
    mpfr_mul_si(X,X,2*M,rnd_re)
    #X=mpmath.mpf(2*M)*Y*pi
    cdef int ok = 1
    ok = incgamma_hint_c(ig,kmh,X)
    mpfr_set(res.value,ig,rnd_re)
    mpfr_mul(res.value,res.value,f1,rnd_re)
    mpfr_mul(res.value,res.value,K1.value,rnd_re)
    mpfr_clear(pi);    mpfr_clear(X)
    mpfr_clear(f1);    mpfr_clear(ig)
    return res
    #fak=f2*mpmath.gammainc(k-1,X)
    #if(k>=0):
    ##    fak=f1*mpmath.power(X,-k)*mpmath.exp(-X)
    #äelse:
    #    fak=f1*mpmath.power(X,-k)*mpmath.exp(-X)
    #fak=f2mpmath.mp.exp(-two*pi*M*y)
    #return mpmath.mpf(K1)*fak

cpdef setup_matrix_for_harmonic_Maass_waveforms(H,Y_in,int M,int Q,principal_parts,use_sym=1,version=1,threads=1):
    if H.group().ncusps()<=2 and use_sym==1:
        return setup_matrix_for_harmonic_Maass_waveforms_sym(H,Y_in,M,Q,principal_parts,version,threads)
    else:
        return setup_matrix_for_harmonic_Maass_waveforms_no_sym(H,Y_in,M,Q,principal_parts,version,threads)
    
@cython.cdivision(True)
@cython.boundscheck(False)
cpdef setup_matrix_for_harmonic_Maass_waveforms_sym(H,RealNumber Y_in,int M,int Q,principal_parts,version=1,threads=1):
    r"""

    Set up the matrix for the system of equations giving the Fourier coefficients of a Harmonic Maass waveforms.

    INPUT:
    
        - ``H`` -- Space of harmonic weak Maass forms
        - ``Y`` -- height of horocycle used for sampling (mpmath real)
        - ``k`` -- weight (mpmath real) 
        - ``M`` -- integer
        - ``Q`` -- integer > M
        - ``PP``-- dict : Principal parts at cusp nr. j = \Sum_ PP(j,n) q^[-n]

    OUTPUT:
    
        - ``W`` -- dictionary
        - ``W['Mf']`` -- M start
        - ``W['nc']`` -- number of cusps
        - ``W['Ms']`` -- M stop
        - ``W['V']``  -- ((Ms-Mf+1)*num_cusps)**2 matrix 
        - ``W['RHS']``  -- ((Ms-Mf+1)*num_cusps) matrix
    
    
    EXAMPLES::

        sage: setup_matrix_for_harmonic_Maass_waveforms_sv(MySubgroup(Gamma0(1)),mpmath.mpf(0.5),mpmath.mpf(0),10,20,{(0,-1):1})

    

    """
    cdef int l,i,jj,j,k,icusp,jcusp,jjcusp,n,ni,li,Ml,Ms,Mf,Qs,Qf,Ql,s,nc,lj
    cdef int prec
    prec = <int>H._prec
    nc=int(H._group.ncusps())
    #cdef mpc_t tmpc
    cdef MPComplexNumber iargpb,iargm,iargpb2
    #cdef RealField_class RF
    #cdef MPComplexField_class CF
    RF = RealField(prec)
    CF = MPComplexField(prec)
    iargpb=CF(0);iargm=CF(0)
    iargpb2=CF(0)
    weight=H._weight
    cdef RealNumber pi,one,two,zero,p,Qfak,ypb,xpb,trn
    cdef mpfr_t twopi,fourpi,nrfourpi,twopiY,fourpiY, kint_t, Y
    cdef mpfr_t tmpr_t,tmpar,tmpar1,tmpab,tmpcos,tmpsin, nr, tmpr, tmpr2
    cdef mpc_t iargpb_t,tmpc_t, tmp1, tmp2
    cdef mpfr_t lr,arg,nrY2pi,kbes
    cdef mpc_t ppc,summa,ppc_minus,summa_minus

    cdef int kinti
    cdef Matrix_complex_dense V
    cdef MPComplexNumber tmpc,ch
    cdef int verbose=int(H._verbose)
    if nc > 2:
        raise ValueError,"Can not use symmetrized algo. with more than 2 cusps!"
    pi=RF.pi() #mpmath_ctx.pi()
    one=RF(1) #mpmath_ctx.mpf(1)
    two=RF(2) #mpmath_ctx.mpf(2)
    zero=RF(0) #mpmath_ctx.mpf(0)
    ypb=RF(0); xpb=RF(0); ch=CF(0); tmpc=CF(0)
    mpfr_init2(twopi,prec)
    mpfr_init2(fourpi,prec)
    mpfr_init2(Y,prec)
    mpfr_init2(twopiY,prec)
    mpfr_init2(fourpiY,prec)
    mpfr_init2(nrfourpi, prec)
    mpfr_init2(kint_t,prec)
    mpc_init2(tmp1,prec)
    mpc_init2(tmp2,prec)
    mpfr_init2(nr, prec)
    mpfr_init2(tmpr, prec)
    mpfr_init2(tmpr2, prec)
    mpfr_init2(tmpr_t,prec)
    mpfr_init2(tmpcos,prec)
    mpfr_init2(tmpsin,prec)
    mpfr_init2(tmpar,prec)
    mpfr_init2(tmpar1,prec)
    mpfr_init2(tmpab,prec)
    mpc_init2(iargpb_t,prec)
    mpc_init2(tmpc_t,prec)
    mpfr_init2(lr,prec)
    mpfr_init2(nrY2pi,prec)
    mpfr_init2(kbes,prec)
    mpfr_init2(arg,prec)
    mpc_init2(ppc,prec)
    mpc_init2(ppc_minus,prec)
    mpc_init2(summa,prec)
    mpc_init2(summa_minus,prec)
    mpfr_mul_ui(twopi,pi.value,2,rnd_re)
    mpfr_mul_ui(fourpi,twopi,2,rnd_re)
    mpfr_set(Y,Y_in.value,rnd_re)
    mpfr_mul(twopiY,twopi,Y,rnd_re)
    mpfr_mul_ui(fourpiY,twopiY,2,rnd_re)
    p=(weight-one)/two
    cdef RealNumber kint=(one-weight)
    cdef RealNumber tmpreal1,tmpreal2
    cdef MPComplexNumber tmpcplx
    tmpcplx = CF(0)
    tmpreal1 = RF(0)
    tmpreal2 = RF(0)
    cdef int is_int=0
    cdef int is_half_int=0
    ## Test if the weight is integral
    if floor(kint)==pceil(kint):
        kinti = int(kint); is_int = 1
    if verbose>0:
        print "is_int=",is_int
        print "kint=",kint
        print "kinti=",kinti
    mpfr_set(kint_t,kint.value,rnd_re)
    if is_int==0:
        ## Check if kint is half-integral.
        if floor(2*kint)==pceil(2*kint):
            is_half_int = 1
            kinti = int(kint-RF(0.5))
    cdef MPComplexNumber ckint
    ckint = CF(kint)
    #print "kint=",kint
    if Q<M:
        Q=M+20
    if H._holomorphic or H._almost_holomorphic:
        Ms=0; Mf=M; Ml=Mf-Ms+1
    else:
        Ms=-M; Mf=M; Ml=Mf-Ms+1
    Qs=1; Qf=Q; Ql=Qf-Qs+1
    Qfak=RF(2*Q)
    cdef int do_mpmath = 0
    if hasattr(H,"_do_mpmath"):
        #if H._do_mpmath<-1:     
        do_mpmath=int(H._do_mpmath) ## By default we want to use mpmath for now...
        if verbose>0:
            print 'NOT using mpmath for incgamma'
    cdef int Qfaki
    Qfaki=2*Q
    # Pullpack points
    if verbose>0:
        print "In setup_matrix_for_harmonic_Maass_waveforms_sym"
        print "Qs,Qf=",Qs,Qf
    cdef mpfr_t* Xm
    cdef mpfr_t*** Xpb=NULL
    cdef mpfr_t*** Ypb=NULL
    cdef mpfr_t**** RCvec=NULL
    cdef int*** CSvec=NULL
    Xm = <mpfr_t*> check_allocarray( sizeof(mpfr_t), Ql )
    if Xm==NULL: raise MemoryError
    for n in range(Ql):
        mpfr_init2(Xm[n],prec)
    Xpb = <mpfr_t***> check_allocarray( sizeof(mpfr_t** ), nc )
    if Xpb==NULL: raise MemoryError
    Ypb = <mpfr_t***> check_allocarray( sizeof(mpfr_t** ), nc )
    if Ypb==NULL: raise MemoryError
    for i in range(nc):
        Xpb[i] = <mpfr_t**>check_allocarray(sizeof(mpfr_t*), nc )
        Ypb[i] = <mpfr_t**>check_allocarray(sizeof(mpfr_t*), nc )
        if Ypb[i]==NULL or Xpb[i]==NULL:
            raise MemoryError
        for j in range(nc):
            Xpb[i][j] = <mpfr_t*>check_allocarray(sizeof(mpfr_t), Ql )
            Ypb[i][j] = <mpfr_t*>check_allocarray(sizeof(mpfr_t), Ql )
            if Ypb[i][j]==NULL or Xpb[i][j]==NULL:
                raise MemoryError
            for n in range(Ql):
                mpfr_init2(Xpb[i][j][n],prec) 
                mpfr_init2(Ypb[i][j][n],prec) 
                mpfr_set_si(Xpb[i][j][n],0,rnd_re)
                mpfr_set_si(Ypb[i][j][n],0,rnd_re)
                #Ypb[i][j][n]=<double>0
    CSvec = <int***>check_allocarray(sizeof(int**), nc )
    if CSvec==NULL: raise MemoryError
    for i from 0<=i<nc:
        CSvec[i] = <int**>check_allocarray(sizeof(int*), nc )
        if CSvec[i]==NULL:
            raise MemoryError
        for j from 0<=j<nc:
            CSvec[i][j] = <int*>check_allocarray(sizeof(int), Ql )
            if CSvec[i][j]==NULL: raise MemoryError
    RCvec = <mpfr_t****>check_allocarray(sizeof(Xpb),nc )
    if RCvec==NULL: raise MemoryError
    for i from 0<=i<nc:
        RCvec[i] = <mpfr_t***>check_allocarray(sizeof(mpfr_t**),nc )
        if RCvec[i]==NULL:
            raise MemoryError
        for j from 0<=j<nc:
            RCvec[i][j] = <mpfr_t**>check_allocarray(sizeof(mpfr_t*),Ql )
            if RCvec[i][j]==NULL:
                raise MemoryError
            for n from 0<=n<Ql:
                RCvec[i][j][n] = <mpfr_t*>check_allocarray(sizeof(mpfr_t), 3 )
                mpfr_init2(RCvec[i][j][n][0],prec)
                mpfr_init2(RCvec[i][j][n][1],prec)
                mpfr_init2(RCvec[i][j][n][2],prec)
                mpfr_set_si(RCvec[i][j][n][0],0,rnd_re)
                mpfr_set_si(RCvec[i][j][n][1],0,rnd_re)
                mpfr_set_si(RCvec[i][j][n][2],0,rnd_re)
                #Cvec[i][j][n]=<double complex>0
    pullback_pts_mpc_new_c_sym(H,Qs,Qf,Y,Xm,Xpb,Ypb,RCvec,CSvec)
    #pb=pullback_pts_mp(H,Qs,Qf,Y,weight)
    #Xm=pb['xm']; Xpb=pb['xpb']; Ypb=pb['ypb']; Cv=pb['cvec']
    s=nc*Ml
    MS = MatrixSpace(CF,s,s)
    V = Matrix_complex_dense(MS,0,True,True)
    cdef list PPplus,PPminus
    cdef dict pp_info
    pp_info = check_principal_parts(H,principal_parts)    
    if verbose>0:
        print "pp_info=",pp_info
        log.debug("pp_info={0}".format(pp_info))
    #if verbose>1:
        print "Yin=",mpfr_get_d(Y,rnd_re)
        #mpc_set(ch.value,Cvec[0][1][0],rnd)
        #print "Cvec[0,1,0]=",ch
    #    mpfr_set(tmpr,Ypb[0][1][0],rnd_re)
    #    print "Ypb[0,1,0]=",mpfr_get_d(tmpr,rnd_re)
    PPplus = pp_info['PPplus']; PPminus = pp_info['PPminus']
    cdef int d = len(PPplus)
    cdef int **variable_a0_plus
    cdef int **variable_a0_minus
    variable_a0_plus = <int**> check_allocarray(sizeof(int*),d)
    variable_a0_minus = <int**> check_allocarray(sizeof(int*),d)
    for j in range(d):
        variable_a0_plus[j] = <int*> check_allocarray(sizeof(int),nc)
        variable_a0_minus[j] = <int*> check_allocarray(sizeof(int),nc)
        for l in range(nc):
            variable_a0_minus[j][l]=int(pp_info['variable_a0_minus'].get(j,{}).get(l,0))
            variable_a0_plus[j][l]=int(pp_info['variable_a0_plus'].get(j,{}).get(l,0))
#        for l in range(nc):
#            variable_a0_minus[j][l]=int(pp_info['variable_a0_minus'].get(j,{}).get(l,0))
#            variable_a0_plus[j][l]=int(pp_info['variable_a0_plus'].get(j,{}).get(l,0))


    log.debug("HERE0")
    cdef int **PPplus_cusp=NULL
    cdef int **PPplus_n=NULL
    cdef int **PPminus_cusp=NULL
    cdef int **PPminus_n=NULL
    cdef mpc_t **PPplus_values=NULL
    cdef mpc_t **PPminus_values=NULL
    cdef mpfr_t **PPplus_lal=NULL
    cdef int num_ppplus = len(pp_info['PPplus'][0])
    cdef int num_ppminus = len(pp_info['PPminus'][0])
    PPplus_cusp = <int **>check_allocarray(sizeof(int*),d)
    if PPplus_cusp==NULL: raise MemoryError
    PPplus_n = <int **>check_allocarray(sizeof(int*),d)
    if PPplus_n==NULL: raise MemoryError
    PPminus_cusp = <int **>check_allocarray(sizeof(int*),d)
    if PPminus_cusp==NULL: raise MemoryError
    PPminus_n = <int **>check_allocarray(sizeof(int*),d)
    if PPminus_n==NULL: raise MemoryError
    PPminus_values = <mpc_t**>check_allocarray(sizeof(mpc_t*),d)
    if PPminus_values==NULL: raise MemoryError
    PPplus_values = <mpc_t**>check_allocarray(sizeof(mpc_t*),d)
    if PPplus_values==NULL: raise MemoryError
    PPplus_lal = <mpfr_t**>check_allocarray(sizeof(mpfr_t*),d)
    if PPplus_lal==NULL: raise MemoryError
    PPminus_lal = <mpfr_t**>check_allocarray(sizeof(mpfr_t*),d)
    log.debug("HERE1")    
    for j in range(d):
        log.debug("HERE2 {0}".format(d))
        PPplus_cusp[j]=NULL;PPplus_n[j]=NULL;PPminus_cusp[j]=NULL;PPminus_n[j]=NULL
        PPplus_values[j]=NULL;PPminus_values[j]=NULL;PPplus_lal[j]=NULL
        PPplus_cusp[j] = <int *>check_allocarray(sizeof(int),num_ppplus)
        if PPplus_cusp[j]==NULL: raise MemoryError
        PPplus_n[j] = <int *>check_allocarray(sizeof(int),num_ppplus)
        if PPplus_n[j]==NULL: raise MemoryError
        PPminus_cusp[j] = <int *>check_allocarray(sizeof(int),num_ppminus)
        if PPminus_cusp[j]==NULL: raise MemoryError
        PPminus_n[j] = <int *>check_allocarray(sizeof(int),num_ppminus)
        if PPplus_n[j]==NULL: raise MemoryError
        PPminus_values[j] = <mpc_t*>check_allocarray(sizeof(mpc_t),num_ppminus)
        if PPminus_values[j]==NULL: raise MemoryError
        PPplus_values[j] = <mpc_t*>check_allocarray(sizeof(mpc_t),num_ppplus)
        if PPplus_values[j]==NULL: raise MemoryError
        PPplus_lal[j] =  <mpfr_t*>check_allocarray(sizeof(mpfr_t),num_ppplus)
        if PPplus_lal[j]==NULL: raise MemoryError
        PPminus_lal[j] =  <mpfr_t*>check_allocarray(sizeof(mpfr_t),num_ppminus)
        if PPminus_lal[j]==NULL: raise MemoryError        
        l = 0
        log.debug("HERE2a")
        for i,jj in pp_info['PPplus'][j].keys():
            if l>=num_ppplus:
                log.critical("Too large l={0} > num_ppplus={1}!".format(l,num_ppplus))
            
            tmpc = CF(pp_info['PPplus'][j][(i,jj)])        
            PPplus_cusp[j][l]=int(i)
            PPplus_n[j][l]=int(jj)
            mpc_init2(PPplus_values[j][l],prec)
            mpc_set(PPplus_values[j][l],tmpc.value,rnd)
            trn = RF(jj)+RF(H.alpha(i)[0])
            mpfr_set(tmpr,trn.value,rnd_re)
            #print "tmpr=",tmpr
            mpfr_init2(PPplus_lal[j][l],prec)
            mpfr_set(PPplus_lal[j][l],tmpr,rnd_re)
            #print "tmpr=",tmpr
            l+=1
        l = 0
        log.debug("HERE2b")        
        for i,jj in pp_info['PPminus'][j].keys():
            if l>=num_ppminus:
                log.critical("Too large l={0} > num_ppminus={1}!".format(l,num_ppminus))
            PPminus_cusp[j][l]=int(i)
            PPminus_n[j][l]=int(jj)
            tmpc = CF(pp_info['PPminus'][j][(i,jj)])
            mpc_init2(PPminus_values[j][l],prec)
            mpc_set(PPminus_values[j][l],tmpc.value,rnd)
            trn = RF(jj)+RF(H.alpha(i)[0])
            mpfr_set(tmpr,trn.value,rnd_re)
            mpfr_init2(PPminus_lal[j][l],prec)
            mpfr_set(PPminus_lal[j][l],tmpr,rnd_re)
            l+=1

    log.debug("HERE3")

    cdef int has_key = 0
    MSRHS = MatrixSpace(CF,s,d)
    RHS = Matrix_complex_dense(MSRHS,0,True,True)

    cdef mpfr_t **nvec=NULL
    nvec = <mpfr_t**> check_allocarray( sizeof(mpfr_t* ), nc )
    cdef RealNumber alpha_tmp
    alpha_tmp = RF(0)
    for icusp in range(nc):
        nvec[icusp]=<mpfr_t*> check_allocarray( sizeof(mpfr_t), Ml )
        alpha_tmp = RF(H.alpha(icusp)[0])
        for l in range(Ml):
            mpfr_init2(nvec[icusp][l],prec)
            mpfr_set(tmpr,alpha_tmp.value,rnd_re)
            mpfr_set_si(nvec[icusp][l],l+Ms,rnd_re)
            mpfr_add(nvec[icusp][l],nvec[icusp][l],tmpr,rnd_re)
    cdef mpfr_t ***ef2cosv=NULL
    cdef mpfr_t ***ef2sinv=NULL
    ef2cosv = <mpfr_t***> check_allocarray( sizeof(mpfr_t** ), nc )
    for icusp in range(nc):
        ef2cosv[icusp]=<mpfr_t**> check_allocarray( sizeof(mpfr_t* ), Ml )
        for n in range(Ml):
            ef2cosv[icusp][n]=<mpfr_t*> check_allocarray( sizeof(mpfr_t ), Ql )
    ef2sinv = <mpfr_t***> check_allocarray( sizeof(mpfr_t** ), nc )
    for icusp in range(nc):
        ef2sinv[icusp]=<mpfr_t**> check_allocarray( sizeof(mpfr_t* ), Ml )
        for n in range(Ml):
            ef2sinv[icusp][n]=<mpfr_t*> check_allocarray( sizeof(mpfr_t ), Ql )
            
    # cdef mpfr_t ****ef1=NULL
    # ef1 = <mpfr_t****> check_allocarray( sizeof(ef2), nc )
    # for icusp in range(nc):
    #     ef1[icusp]=<mpfr_t***> check_allocarray( sizeof(mpfr_t** ), nc )
    #     for jcusp in range(nc):
    #         ef1[icusp][jcusp]=<mpfr_t**> check_allocarray( sizeof(mpfr_t* ), Ml )
    #         for n in range(Ml):
    #             ef1[icusp][jcusp][n]=<mpfr_t*> check_allocarray( sizeof(mpfr_t ), Ql )
    #             for j in range(Ql):
    #                 mpfr_init2(ef1[icusp][jcusp][n][j],prec)
    cdef mpfr_t **** besv=NULL
    cdef mpc_t **** besv_minus=NULL    
    # Include the special functions etc. for non-principal parts and holom. (+) principal parts
    besv = <mpfr_t****> check_allocarray( sizeof(ef2cosv), nc )
    for icusp in range(nc):
        besv[icusp]=<mpfr_t***> check_allocarray( sizeof(mpfr_t** ), nc )
        for jcusp in range(nc):
            besv[icusp][jcusp]=<mpfr_t**> check_allocarray( sizeof(mpfr_t* ), Ml )
            for n in range(Ml):
                besv[icusp][jcusp][n]=<mpfr_t*> check_allocarray( sizeof(mpfr_t ), Ql )
                for j in range(Ql):
                    mpfr_init2(besv[icusp][jcusp][n][j],prec)
    # Include the special functions etc. for non-holom. (-) principal parts                    
    besv_minus = <mpc_t****> check_allocarray( sizeof(ef2cosv), nc )
    for icusp in range(nc):
        besv_minus[icusp]=<mpc_t***> check_allocarray( sizeof(mpc_t** ), nc )
        for jcusp in range(nc):
            besv_minus[icusp][jcusp]=<mpc_t**> check_allocarray( sizeof(mpc_t* ), Ml )
            for n in range(Ml):
                besv_minus[icusp][jcusp][n]=<mpc_t*> check_allocarray( sizeof(mpc_t ), Ql )
                for j in range(Ql):
                    mpc_init2(besv_minus[icusp][jcusp][n][j],prec)
    cdef mpfr_t **** ef1cosv=NULL
    cdef mpfr_t **** ef1sinv=NULL
    ef1cosv = <mpfr_t****> check_allocarray( sizeof(ef2cosv), nc )
    for icusp in range(nc):
        ef1cosv[icusp]=<mpfr_t***> check_allocarray( sizeof(mpfr_t** ), nc )
        for jcusp in range(nc):
            ef1cosv[icusp][jcusp]=<mpfr_t**> check_allocarray( sizeof(mpfr_t* ), Ml )
            for n in range(Ml):
                ef1cosv[icusp][jcusp][n]=<mpfr_t*> check_allocarray( sizeof(mpfr_t ), Ql )
                for j in range(Ql):
                    mpfr_init2(ef1cosv[icusp][jcusp][n][j],prec)
    ef1sinv = <mpfr_t****> check_allocarray( sizeof(ef2cosv), nc )
    for icusp in range(nc):
        ef1sinv[icusp]=<mpfr_t***> check_allocarray( sizeof(mpfr_t** ), nc )
        for jcusp in range(nc):
            ef1sinv[icusp][jcusp]=<mpfr_t**> check_allocarray( sizeof(mpfr_t* ), Ml )
            for n in range(Ml):
                ef1sinv[icusp][jcusp][n]=<mpfr_t*> check_allocarray( sizeof(mpfr_t ), Ql )
                for j in range(Ql):
                    mpfr_init2(ef1sinv[icusp][jcusp][n][j],prec)
    cdef double eps
    eps = 2.0**float(1-H._dprec)
    cdef int not_holom = int(not H._holomorphic)
    cdef int is_weak = int(H._weak)
    cdef int ok = 1
    if verbose>0:
        print "Ml, Ql =",Ml, Ql
    for n in range(Ml):
        for jcusp in range(nc):
                mpfr_set(nr,nvec[jcusp][n],rnd_re)
                setcossin2(ef2cosv[jcusp][n],ef2sinv[jcusp][n],Xm,nr,Ql,prec)
                if verbose>2:
                    printf("done with ef2cosv[%d][%d], ef2sinv[%d][%d]\n",jcusp,n,jcusp,n)
                    printf("ef2cosv[%d][%d][0]=%f\n",jcusp,n,mpfr_get_d(ef2cosv[jcusp][n][0],rnd_re))
                #mpfr_set(nr,nvec[jcusp][n],rnd_re)
                #if verbose>0 and n==0:
                #    print "n=",n
                #    printf("nn=%f",mpfr_get_d(nr,rnd_re))
                #nr=nvec[n,jcusp]
                mpfr_mul(nrfourpi,nr,fourpi,rnd_re)
                #nrfourpi=nr*fourpi
                if True:    
                    for icusp in xrange(nc):
                        if verbose>2:
                            printf("jcusp=%i\n",jcusp)
                        #setcossin_besv(ef1cosv[icusp][jcusp][n], ef1sinv[icusp][jcusp][n], besv[icusp][jcusp][n], Ypb[icusp][jcusp][n], variable_a0_minus[0][] nr, eps, kinti, kint_t, is_int, is_half_int, do_mpmath,)
                        for j in xrange(Ql):
                            ## ef1 contains -Xpb*l
                            if mpfr_zero_p(Ypb[icusp][jcusp][j])<>0:
                                #mpfr_set_si(ef1[icusp][jcusp][xn][j],0,rnd_re) 
                                continue
                            if verbose>2:
                                printf("n,icusp,jcusp,j=%d,%d,%d,%d\n",n,icusp,jcusp,j)
                            mpfr_mul(tmpar,Xpb[icusp][jcusp][j],nr,rnd_re)
                            mpfr_add(tmpar,tmpar,RCvec[icusp][jcusp][j][1],rnd_re)
                            mpfr_cos(ef1cosv[icusp][jcusp][n][j],tmpar,rnd_re)
                            mpfr_sin(ef1sinv[icusp][jcusp][n][j],tmpar,rnd_re)
                            mpfr_mul(tmpr_t,twopi,Ypb[icusp][jcusp][j],rnd_re)
                            mpfr_mul(tmpr_t,tmpr_t,nr,rnd_re)
                            mpfr_neg(tmpr_t,tmpr_t,rnd_re)
                            if verbose>2:
                                printf("arg=%f nr=%f",mpfr_get_d(tmpr_t,rnd_re),mpfr_get_d(nr,rnd_re))
                            #if verbose>1 and n+Ms==0:
                            #    printf("arg=%f \n n=%f \n",mpfr_get_d(tmpr_t,rnd_re),mpfr_get_d(nr,rnd_re))
                            mpfr_exp(besv[icusp][jcusp][n][j],tmpr_t,rnd_re)
                            #if verbose>1 and n+Ms==0:
                            #    mpfr_set(tmpreal1.value,besv[icusp][jcusp][n][j],rnd_re)
                            #    print "exp(arg)=",tmpreal1
                            if mpfr_cmp_d(nr,eps)>0:
                                pass
                            elif mpfr_cmp_d(nr,-eps)<0 and not_holom==1 and is_weak==1:
                                mpfr_abs(tmpr,nrfourpi,rnd_re)
                                mpfr_mul(tmpr,tmpr,Ypb[icusp][jcusp][j],rnd_re)
                                ok = 1
                                if (is_int==1 or is_half_int==1) and do_mpmath==0:
                                    #try:
                                    if is_int==1:
                                        if kinti>0:
                                            ok = incgamma_pint_c(tmpr2,kinti,tmpr,verbose)
                                        else:
                                            if verbose>2:
                                                print "doing incgamma_nint_c"
                                            if kinti <>0:
                                                mpfr_set(tmpreal1.value,tmpr,rnd_re)
                                                tmpreal2 = RF(mpmath.mp.gammainc(kint,tmpreal1).real)
                                                mpfr_set(tmpr2,tmpreal2.value,rnd_re)
                                                ok = 0
                                            else:
                                                ok = incgamma_nint_c(tmpr2,kinti,tmpr,verbose)
                                    elif is_half_int==1:
                                        ok = incgamma_hint_c(tmpr2,kinti,tmpr)                                
                                        #except ArithmeticError: ## In case we can not achieve the required error
                                        #tmpr2 = RF(mpmath.mp.gammainc(kint,tmpr).real)                    
                                if (verbose>0 and ok<>0) or verbose>2:
                                    print "ok={0}, tmpr={1} Gamma={2} do_mpmath={3}, kint={4}".format(ok,mpfr_get_d(tmpr,rnd_re),mpfr_get_d(tmpr2,rnd_re),do_mpmath,kinti)
                                if ok <> 0:
                                    #tmpr2 = RF(mpmath.mp.gammainc(kint,mpfr_get_d(tmpr,rnd_re)).real).value
                                    if verbose>0:
                                        printf('return value from incomplete gamma not good')
                                        return -1
                                #if verbose>1 and icusp==1 and n==0:
                                    #print "f1[{0},{1},{2},{3}]={4}".format(icusp,jcusp,n,j,iargpb)
                                #    print "be[{0},{1},{2},{3}]={4}".format(icusp,jcusp,n,j,mpfr_get_d(tmpr2,rnd_re))
                                mpfr_mul(besv[icusp][jcusp][n][j],besv[icusp][jcusp][n][j],tmpr2,rnd_re)
                            elif mpfr_cmp_d(nr,-eps)<0:
                                mpfr_set_si(besv[icusp][jcusp][n][j],0,rnd_re)
                            else:
                                ## Note that principal parts have to be determined
                                ## and will appear in the right hand side
                                if verbose>0:
                                    print "we are here with n+Ms=",0
                                if variable_a0_plus[0][jcusp]==1:
                                    mpfr_set_si(besv[icusp][jcusp][n][j],1,rnd_re)
                                elif variable_a0_minus[0][jcusp]==1:
                                    mpfr_set(besv[icusp][jcusp][n][j],Ypb[icusp][jcusp][j],rnd_re)
                                    if kinti==0:
                                        mpfr_log(besv[icusp][jcusp][n][j],besv[icusp][jcusp][n][j],rnd_re)
                                    else:
                                        mpfr_pow(besv[icusp][jcusp][n][j],besv[icusp][jcusp][n][j],kint_t,rnd_re)
                                        if verbose>2:
                                            mpfr_set(tmpreal1.value,besv[icusp][jcusp][n][j],rnd_re)
                                            print "set besv[",icusp,jcusp,n,j,"=",tmpreal1

                                else:
                                    ## If none of them are variables this means
                                    ## that they are both set in the principal parts
                                    mpfr_set_si(besv[icusp][jcusp][n][j],0,rnd_re)
    for i in range(num_ppminus):
        jcusp = PPminus_cusp[k][i]
        n = PPminus_n[k][i]
        if n < 0:
            raise ValueError,"invalid - principal part"
        if verbose>0:
            print "PP- at cusp={0} and n={1}".format(jcusp,n)
        for icusp in range(nc):
            mpfr_set(nr,nvec[jcusp][n-Ms],rnd_re)
            mpfr_mul(nrfourpi,nr,fourpi,rnd_re)
            #nrfourpi=nr*fourpi
            for j in xrange(Ql):
                ## ef1 contains -Xpb*l
                if mpfr_zero_p(Ypb[icusp][jcusp][j])<>0:
                    #mpfr_set_si(ef1[icusp][jcusp][n][j],0,rnd_re) 
                    continue
                if n==0:
                    if verbose>0:
                        printf("Ypb=%f \n",mpfr_get_d(Ypb[icusp][jcusp][j],rnd_re))
                        printf("kint_t=%f \n",mpfr_get_d(kint_t,rnd_re))
                    mpfr_pow(tmpr_t,Ypb[icusp][jcusp][j],kint_t,rnd_re)
                    mpc_set_fr(besv_minus[icusp][jcusp][n-Ms][j],tmpr_t,rnd)
                    if verbose>0:
                        mpc_set(tmpcplx.value,besv_minus[icusp][jcusp][n-Ms][j],rnd)
                        print "besv_minus0[",icusp,jcusp,n-Ms,j,"]=",tmpcplx
                       
                else:
                    # Gamma(1-k,-4piny)+(-1)^{1-k}\pi i/Gamma(k)
                    #  x e^(-2*pi*Ypb*n)
                    ##
                    mpfr_mul(tmpr_t,twopi,Ypb[icusp][jcusp][j],rnd_re) # 2pi y
                    mpfr_mul(tmpr_t,tmpr_t,nr,rnd_re) # 2pi y*n
                    mpfr_neg(tmpr_t,tmpr_t,rnd_re)  # -2pi y*n
                    mpfr_exp(tmpr_t,tmpr_t,rnd_re) # e^(-2pi y*n)

                    mpfr_mul(tmpr,nrfourpi,Ypb[icusp][jcusp][j],rnd_re)
                    mpfr_neg(tmpr,tmpr,rnd_re)
                    mpfr_set(tmpreal1.value,tmpr,rnd_re)
                    tmpcplx = CF(kint).gamma_inc(tmpreal1)   ## gamma(1-k,4pi |n| y)
                    tmpcplx2 = CF(weight).gamma() #(mptmp.real,mptmp.imag)
                    tmpcplx = tmpcplx+CF(0,1)**(3-2*weight)*pi/tmpcplx2
                    mpc_set(besv_minus[icusp][jcusp][n-Ms][j],tmpcplx.value,rnd)
                    mpc_mul_fr(besv_minus[icusp][jcusp][n-Ms][j],besv_minus[icusp][jcusp][n-Ms][j],tmpr_t,rnd)
                    #mpc_mul_fr(besv_minus[icusp][jcusp][n-Ms][j],besv_minus[icusp][jcusp][n-Ms][j],tmpr_t,rnd)
                    if verbose>0:
                        print "tmpcplx=",tmpcplx
                        print "tmpcplx2=",tmpcplx2
                        mpc_set(tmpcplx.value,besv_minus[icusp][jcusp][n-Ms][j],rnd)
                        print "besv_minus1[",icusp,jcusp,n-Ms,j,"]=",tmpcplx
            
                                    
    cdef int nrows,ncols
    nrows = int(V.nrows()); ncols = int(V.ncols())
    if verbose>1:
        printf("nrows,ncols=%d,%d\n",nrows,ncols)
        printf("Ml,nc,Ql=%d,%d,%d\n",Ml,nc,Ql)
    for n in xrange(nrows):
        for l in xrange(ncols):
            mpc_init2(V._matrix[n][l],prec)
            mpc_set_ui(V._matrix[n][l],0,rnd)
    for l in prange(Ml, nogil=True):
        setV(V._matrix, RCvec,CSvec,besv,Ypb,ef1cosv,ef1sinv,ef2cosv,ef2sinv,nc, Ql, Ml, l, prec)                                        
    if verbose>1:
        print "Mi,Ms=",Ms,Mf
        print "V1(",0,0,")=",V[0,0]
        print "V1(",1,0,")=",V[1,0]
        print "V1(",0,1,")=",V[0,1]
        print "V1(",1,1,")=",V[1,1]
        #if V.ncols()>=11:
        #    print "V1(",0,11,")=",V[0,11]
    for n in range(nrows):
        for l in range(ncols):        
            mpc_div_ui(V._matrix[n][l],V._matrix[n][l],Qfaki,rnd)
            #V[n,l]=V[n,l]/Qfak
    if verbose>1: 
        print "V1(",0,0,")=",V[0,0]
        print "V1(",1,0,")=",V[1,0]            
        print "V1(",0,1,")=",V[0,1]
        print "V1(",1,1,")=",V[1,1]
        #if V.ncols()>=11:
        #    print "V1(",0,11,")=",V[0,11]
    cdef MPComplexNumber f1,f2
    f1 = CF(0); f2=CF(0)
    #verbose=3
    for n in range(Ml):
        for icusp in range(nc):
            sig_check()
            mpfr_set(nr,nvec[icusp][n],rnd_re)
            mpfr_mul(nrY2pi,nr,twopiY,rnd_re)
            mpfr_mul(nrfourpi,nr,fourpi,rnd_re)
            ni=Ml*icusp+n
            if mpfr_cmp_d(nr,eps)>0:
                mpfr_neg(kbes,nrY2pi,rnd_re)
                mpfr_exp(kbes,kbes,rnd_re)      #kbes=e^-2piy
            elif mpfr_cmp_d(nr,-eps)<0 and not_holom==1 and is_weak==1:
                #kbes=RF(mpmath.mp.gammainc(kint,abs(nr)*fourpiY))
                #kbes=ckint.gamma_inc(abs(nr)*fourpiY).real()                
                mpfr_mul(tmpr,nr,fourpiY,rnd_re)
                mpfr_abs(tmpr,tmpr,rnd_re)
                ok = 1
                if (is_int==1 or is_half_int==1) and  do_mpmath==0:
                    #try:
                    if is_int==1:
                        if kinti>0:
                            ok = incgamma_pint_c(kbes,kinti,tmpr)
                        else:
                            ok = incgamma_nint_c(kbes,kinti,tmpr)
                    elif is_half_int==1:
                        ok = incgamma_hint_c(kbes,kinti,tmpr)                                
                    #except ArithmeticError: ## In case we can not achieve the required error
                    #    kbes = RF(mpmath.mp.gammainc(kint,tmpr).real)                    
                #if (verbose>0 and ok<>0) or verbose>2:
                #    print "ok={0}, tmpr={1} Gamma={2} do_mpmath={3}".format(ok,mpfr_get_d(tmpr,rnd_re),mpfr_get_d(kbes,rnd_re),do_mpmath)
                #    print "ok={0}, tmpr={1}".format(ok,mpfr_get_d(tmpr, rnd_re))
                if ok <> 0:
                    #kbes = RF(mpmath.mp.gammainc(kint,tmpr).real) #TODO: Fix call (tmpr is mpfr_t now)
                    ret=-1
                    break
                mpfr_abs(tmpr2,nrY2pi,rnd_re)
                mpfr_exp(tmpr2,tmpr2,rnd_re)
                mpfr_mul(kbes,kbes,tmpr2,rnd_re)
                # Gamma(1-k,4pi|n|y)e^(-2piny)
                #if verbose>1 and ni==0:
                #    print "expfax=",mpfr_get_d(tmpr2,rnd_re)
                #    print "Igamma(",kinti,",Arg)*exp()=",mpfr_get_d(tmpr,rnd_re)
            elif mpfr_cmp_d(nr,-eps)<0:
                mpfr_set_ui(kbes,0,rnd_re) #=zero
            else:
                if variable_a0_plus[0][icusp]==1:
                    mpfr_set_ui(kbes,1,rnd_re) #=one
                elif variable_a0_minus[0][icusp]==1:
                    if kinti<>0:
                        mpfr_pow(kbes,Y,kint_t,rnd_re)
                    else:
                        mpfr_log(kbes,Y,rnd_re)
                    #kbes=Y**kint
                else:
                    #kbes=one
                    #raise ValueError,"Need to specify at least one of the n=0 terms!"
                    mpfr_set_ui(kbes,0,rnd_re)
                    #kbes=zero
            #with gil:
            mpc_sub_fr(V._matrix[ni][ni],V._matrix[ni][ni],kbes,rnd)
            ##
            ## Setting the right hand side using the principal parts.
            
            for k in range(d):
                mpc_set_ui(RHS._matrix[ni][k],0,rnd) #=CF(0)
                ## first the holomorphic ones
                for i in range(num_ppplus):
                    mpc_set(ppc,PPplus_values[k][i],rnd)
                    jcusp = PPplus_cusp[k][i]
                    l = PPplus_n[k][i]
                    ### Note: The constant term is treated below
                    if mpc_is_zero(ppc)==1: #==zero or PPplus[k][(jcusp,l)]==0:
                        continue
                    mpc_set_ui(summa,0,rnd) #=CF(0)
                    for j in range(Ql):   
                        if mpfr_zero_p(Ypb[icusp][jcusp][j])<>0:
                            continue
                        mpfr_mul(tmpr_t,Ypb[icusp][jcusp][j],twopi,rnd_re)
                        mpfr_neg(tmpr_t,tmpr_t,rnd_re)
                        mpfr_mul(tmpr_t,tmpr_t,PPplus_lal[k][i],rnd_re)
                        mpfr_exp(tmpr_t,tmpr_t,rnd_re)
                        mpfr_mul(tmpab,tmpr_t,RCvec[icusp][jcusp][j][0],rnd_re)
                        # e^(-2pi*Ypb*(n+alpha_i))*v(icusp,jcusp,j)
                        mpfr_mul(tmpar,PPplus_lal[k][i],Xpb[icusp][jcusp][j],rnd_re)
                        mpfr_mul(tmpr_t,nr,Xm[j],rnd_re)
                        mpfr_sub(tmpar,tmpar,tmpr_t,rnd_re)
                        mpfr_add(tmpar,tmpar,RCvec[icusp][jcusp][j][1],rnd_re)
                        # e^(2pi i *(Xpb(n+alpha_i)-Xm[j]*n + argv)
                        #if verbose>3:
                        #    print "f1(",j,")=",f1
                        #    print "arg=",lr,"*",xpb,"-",nr,"*",Xm[j],"=",arg
                        #    print "f2(",j,")=",f2
                        if mpfr_get_si(RCvec[icusp][jcusp][j][2],rnd_re) % 2 == 0:
                            mpfr_cos(tmpar1,tmpar,rnd_re)
                            mpfr_mul_ui(tmpar1,tmpar1,2,rnd_re)
                            mpc_set_fr(tmp2,tmpar1,rnd)                                     
                        else:
                            mpfr_sin(tmpar1,tmpar,rnd_re)
                            mpc_set_si_si(tmp2,0,2,rnd)
                            mpc_mul_fr(tmp2,tmp2,tmpar1,rnd)
                        mpc_mul_fr(tmp2,tmp2,tmpab,rnd)
                        #mpc_mul_fr(tmp2.value,tmp2.value,RCvec[icusp][jcusp][j][0]):
                        mpc_add(summa,summa,tmp2,rnd)
                        #if verbose>3:
                        #    print "tmpc=",tmpc
                        #    print "summa=",summa                 
#                        summa=summa+ch*tmpc
                    #if verbose>2:
                    #    print "summa=",RR(summa.real()),RR(summa.imag())
                    #RHS[ni,k]=RHS[ni,k]+summa/Qfak
                    mpc_div_ui(summa,summa,Qfaki,rnd)
                    mpc_mul(summa,summa,ppc,rnd)
                    mpc_add(RHS._matrix[ni][k],RHS._matrix[ni][k],summa,rnd)
                #if verbose>2:
                #    print "RHS0[",ni,k,"]=",RHS[ni][k]
                #    print "icusp,n+Ms=",icusp,n+Ms
                #    print "PPplus.keys=",PPplus[k].keys()
                ## Since 0 can also be a variable or set here we subtract the contribution.
                if n+Ms==0:
                    has_key = 0
                    for j in range(num_ppplus):
                        if PPplus_cusp[k][j]==icusp and PPplus_n[k][j]==n+Ms:
                            has_key = 1
                            mpc_set(ppc,PPplus_values[k][j],rnd)
                            break
                    if has_key==1:
                        print "has key! + n+Ms=",n+Ms
                        mpfr_abs(tmpr_t,nrY2pi,rnd_re)
                        mpfr_exp(tmpr_t,tmpr_t,rnd_re)
                        mpc_mul_fr(tmpc_t,ppc,tmpr_t,rnd)
                        mpc_sub(RHS._matrix[ni][k],RHS._matrix[ni][k],tmpc_t,rnd)
                        #RHS[ni,k]=RHS[ni,k]-ppc*(-nrY2pi).exp()
                        if verbose>0:
                            print "num_ppm=",num_ppminus
                #
                # Add contributions of non-holomorphic principal parts.
                #
                for i in range(num_ppminus):
                    jcusp = PPminus_cusp[k][i]
                    l = PPminus_n[k][i] - Ms
                    mpc_set(ppc_minus,PPminus_values[k][i],rnd)
                    if verbose>0:
                        mpc_set(tmpcplx.value,ppc_minus,rnd)                        
                        print "ppart_minus=",tmpcplx
                    if mpc_is_zero(ppc_minus)==1:
                        continue
                    mpc_set_ui(summa_minus,0,rnd) 
                    for j in range(Ql):
                        sig_check()
                        if mpfr_zero_p(Ypb[icusp][jcusp][j])<>0:
                            continue
                        # RCvec = [|v|, arg(v)]
                        if verbose>0:
                            mpc_set(tmpcplx.value,besv_minus[icusp][jcusp][l][j],rnd)
                            print "besv_minus2[",icusp,jcusp,l,j,"]=",tmpcplx
                            mpfr_set(tmpreal1.value,RCvec[icusp][jcusp][j][0],rnd_re)
                            print "Rcvec[",icusp,jcusp,j,0,"]=",tmpreal1
                            mpfr_set(tmpreal1.value,RCvec[icusp][jcusp][j][1],rnd_re)
                            print "Rcvec[",icusp,jcusp,j,1,"]=",tmpreal1
                            mpfr_set(tmpreal1.value,RCvec[icusp][jcusp][j][2],rnd_re)
                            print "Rcvec[",icusp,jcusp,j,2,"]=",tmpreal1                                                        
                        mpfr_mul(tmpar,PPminus_lal[k][i],Xpb[icusp][jcusp][j],rnd_re)
                        if verbose>0:
                            mpfr_set(tmpreal1.value,tmpar,rnd_re)
                            print "tmpar=PPminus_lal[",k,i,"]*Xpb[",icusp,jcusp,j,"]=",tmpreal1
                            mpfr_set(tmpreal1.value,Xm[j],rnd_re)
                            print "Xm[",j,"]=",tmpreal1
                        mpfr_mul(tmpr_t,nr,Xm[j],rnd_re)
                        mpfr_sub(tmpar,tmpar,tmpr_t,rnd_re)
                        #mpfr_add(tmpar,tmpar,ef2[icusp][n][j],rnd_re)
                        mpfr_add(tmpar,tmpar,RCvec[icusp][jcusp][j][1],rnd_re)
                        ## tmpar = (l+alpha(icusp))*Xpb-(n+alpha(jcusp))*Xm+arg(v)
                        ## v includes (cz+d)^(-k) and multiplier
                        if verbose>0:
                            mpfr_set(tmpreal1.value,tmpar,rnd_re)
                            print "tmpar+Rvec[",icusp,jcusp,j,1,"]=",tmpreal1                

                        if mpfr_get_si(RCvec[icusp][jcusp][j][2],rnd_re) % 2 == 0:
                            mpfr_cos(tmpr_t,tmpar,rnd_re)
                            mpfr_mul_ui(tmpr_t,tmpr_t,2,rnd_re)
                            mpc_set_fr(tmp2,tmpr_t,rnd)                                     
                        else:
                            mpfr_sin(tmpr_t,tmpar,rnd_re)
                            mpc_set_si_si(tmp2,0,2,rnd)
                            mpc_mul_fr(tmp2,tmp2,tmpr_t,rnd)
                        if verbose>3:
                            mpc_set(tmpcplx.value,tmp2,rnd)
                            print "tmp2(",j,")=",tmpcplx
                        mpc_mul_fr(tmp2,tmp2,RCvec[icusp][jcusp][j][0],rnd)
                        mpc_mul(tmp2,tmp2,besv_minus[icusp][jcusp][l][j],rnd)
                        # tmp2 = |v|*sin/cos(tmpar)* [exponential / gamma / bessel]
                        mpc_add(summa_minus,summa_minus,tmp2,rnd)
                        #summa_minus=summa_minus+ch*tmpc
                        #summa=summa+ef2[n(2):,j,icusp]*tmpc*ch*ppc
                        #if verbose>2:
                        #    print "j=",j
                        #    print "ypb=",ypb
                        #    print "ch=",ch
                        #    #print "tmpc_minus=",tmpc_minus
                        #    print "tmpc=",tmpc
                        #    print "summa-(",j,")=",summa_minus
                    if verbose>0:
                        mpc_set(tmpcplx.value,summa_minus,rnd)
                        print "Summa-(",ni,")=",tmpcplx
                    mpc_div_ui(summa_minus,summa_minus,Qfaki,rnd)
                    mpc_mul(summa_minus,summa_minus,ppc_minus,rnd)
                    mpc_add(RHS._matrix[ni][k],RHS._matrix[ni][k],summa_minus,rnd)
                    #RHS[ni,k]=RHS[ni,k]+ppc_minus*summa_minus/Qfak
                #print "summa(",ni,")=",RHS[ni,k]
                #if verbose>2:
                #    print "RHS2[{0},{1}]={2}".format(ni,k,RHS[ni,k])

                ## The only free variable we allow in the non-holomorphic part
                ##is a0 and this
                #if n+Ms==0: #variable_a0_minus[0][icusp]==1 and n+Ms==0:
                has_key = 0
                for j in range(num_ppminus):
                    if PPminus_cusp[k][j]==icusp and PPminus_n[k][j]==n+Ms:
                        has_key = 1
                        mpc_set(ppc_minus,PPminus_values[k][j],rnd)
                        break
                if has_key==1:
                    print "has key! - n+Ms=",n+Ms
                    if n+Ms == 0:
                        if kinti==0:
                            mpfr_log(tmpr_t,Y,rnd_re)
                        else:
                            mpfr_pow(tmpr_t,Y,kint_t,rnd_re)
                    else:

                        mpfr_neg(tmpr_t,nrY2pi,rnd_re)  # -2pi Y*n
                        mpfr_exp(tmpr_t,tmpr_t,rnd_re) # e^(-2pi Y*n)

                        mpfr_mul(tmpr,nrfourpi,Y,rnd_re)
                        mpfr_neg(tmpr,tmpr,rnd_re)
                        mpfr_set(tmpreal1.value,tmpr,rnd_re)
                        print "tmprel=",tmpreal1
                        tmpcplx = CF(kint).gamma_inc(tmpreal1)   ## gamma(1-k,- 4pi n y)
                        tmpcplx2 = CF(weight).gamma() #(mptmp.real,mptmp.imag)
                        tmpcplx = tmpcplx+CF(0,1)**(3-2*weight)*pi/tmpcplx2
                        mpc_mul(ppc_minus,ppc_minus,tmpcplx.value,rnd)
                    mpc_mul_fr(ppc_minus,ppc_minus,tmpr_t,rnd)
                    mpc_sub(RHS._matrix[ni][k],RHS._matrix[ni][k],ppc_minus,rnd)
                    #RHS[ni,k]=RHS[ni,k]-ppc_minus*Y**kint
                    #print "subtracting:",ppc_minus*mpmath.power(Y,kint)
                #if verbose>2:
                #    print "RHS3[{0},{1}]={2}".format(ni,k,RHS[ni,k])
                    #print "RHS(",ni,")-=",ppc*mpmath_ctx.exp(-nrY2pi)
    if verbose>1: 
        for n in range(2):
            print "V2(",n,n,")=",V[n,n]
        #if V.ncols()>=11:
        #    print "V2(",0,11,")=",V[0,11]
            # Clearing up allocated variables

    for j in range(d):
        #printf("a0jm=%p",variable_a0_minus[j])
        sig_free(variable_a0_minus[j])
        sig_free(variable_a0_plus[j])        
        sig_free(PPplus_cusp[j])
        sig_free(PPplus_n[j])
        sig_free(PPminus_cusp[j])
        sig_free(PPminus_n[j])
        for l in range(num_ppplus):
            mpc_clear(PPplus_values[j][l])
            mpfr_clear(PPplus_lal[j][l])
        sig_free(PPplus_lal[j])
        sig_free(PPplus_values[j])
        for l in range(num_ppminus):
            mpc_clear(PPminus_values[j][l])
            mpfr_clear(PPminus_lal[j][l])
        sig_free(PPminus_lal[j])
        sig_free(PPminus_values[j])
        
    sig_free(PPminus_lal)
    sig_free(PPplus_lal)
    sig_free(PPplus_cusp)
    sig_free(PPplus_n)
    sig_free(PPplus_values)
    sig_free(PPminus_cusp)
    sig_free(PPminus_n)
    sig_free(PPminus_values)

    sig_free(variable_a0_minus)
    sig_free(variable_a0_plus)


    if Ypb<>NULL:
        for i in range(nc):
            if Ypb[i]<>NULL:
                for j in range(nc):
                    if Ypb[i][j]<>NULL:
                        sig_free(Ypb[i][j])
                sig_free(Ypb[i])
        sig_free(Ypb)
    if Xpb<>NULL:
        for i in range(nc):
            if Xpb[i]<>NULL:
                for j in range(nc):
                    if Xpb[i][j]<>NULL:
                        sig_free(Xpb[i][j])
                sig_free(Xpb[i])
        sig_free(Xpb)
    if CSvec<>NULL:
        for i in range(nc):
            if CSvec[i]<>NULL:
                for j in range(nc):
                    if CSvec[i][j]<>NULL:
                        sig_free(CSvec[i][j])
                sig_free(CSvec[i])
        sig_free(CSvec)
    if RCvec<>NULL:
        for i in range(nc):
            if RCvec[i]<>NULL:
                for j in range(nc):
                    if RCvec[i][j]<>NULL:
                        for n in range(Ql):
                            if RCvec[i][j][n]<>NULL:
                                mpfr_clear(RCvec[i][j][n][0])
                                mpfr_clear(RCvec[i][j][n][1])
                                mpfr_clear(RCvec[i][j][n][2])
                                sig_free(RCvec[i][j][n])
                sig_free(RCvec[i])
        sig_free(RCvec)
    if nvec<>NULL:
        for i in range(nc):
            if nvec[i]<>NULL:
                for l in range(Ml):
                    mpfr_clear(nvec[i][l])
                sig_free(nvec[i])
        sig_free(nvec)
    if Xm<>NULL:
        for l in range(Ql):
            mpfr_clear(Xm[l])
        sig_free(Xm)
    if ef2cosv<>NULL:
        for icusp in range(nc):
            if ef2cosv[icusp]<>NULL:
                for n in range(Ml):
                    if ef2cosv[icusp][n]<>NULL:
                        for j in range(Ql):
                            mpfr_clear(ef2cosv[icusp][n][j])
                        sig_free(ef2cosv[icusp][n])
                sig_free(ef2cosv[icusp])
        sig_free(ef2cosv)
    if ef2sinv<>NULL:
        for icusp in range(nc):
            if ef2sinv[icusp]<>NULL:
                for n in range(Ml):
                    if ef2sinv[icusp][n]<>NULL:
                        for j in range(Ql):
                            mpfr_clear(ef2sinv[icusp][n][j])
                        sig_free(ef2sinv[icusp][n])
                sig_free(ef2sinv[icusp])
        sig_free(ef2sinv)
        
    if ef1cosv<>NULL:
        for icusp in range(nc):
            if ef1cosv[icusp]<>NULL:
                for jcusp in range(nc):
                    if ef1cosv[icusp][jcusp]<>NULL:
                        for n in range(Ml):
                            if ef1cosv[icusp][jcusp][n]<>NULL:
                                for j in range(Ql):
                                    mpfr_clear(ef1cosv[icusp][jcusp][n][j])
                                sig_free(ef1cosv[icusp][jcusp][n])
                        sig_free(ef1cosv[icusp][jcusp])                        
                sig_free(ef1cosv[icusp])
        sig_free(ef1cosv)
    if ef1sinv<>NULL:
        for icusp in range(nc):
            if ef1sinv[icusp]<>NULL:
                for jcusp in range(nc):
                    if ef1sinv[icusp][jcusp]<>NULL:
                        for n in range(Ml):
                            if ef1sinv[icusp][jcusp][n]<>NULL:
                                for j in range(Ql):
                                    mpfr_clear(ef1sinv[icusp][jcusp][n][j])
                                sig_free(ef1sinv[icusp][jcusp][n])
                        sig_free(ef1sinv[icusp][jcusp])                        
                sig_free(ef1sinv[icusp])
        sig_free(ef1sinv)
    if besv<>NULL:
        for icusp in range(nc):
            if besv[icusp]<>NULL:
                for jcusp in range(nc):
                    if besv[icusp][jcusp]<>NULL:
                        for n in range(Ml):
                            if besv[icusp][jcusp][n]<>NULL:
                                for j in range(Ql):
                                    mpfr_clear(besv[icusp][jcusp][n][j])
                                sig_free(besv[icusp][jcusp][n])
                        sig_free(besv[icusp][jcusp])                        
                sig_free(besv[icusp])
        sig_free(besv)
    if besv_minus<>NULL:
        for icusp in range(nc):
            if besv_minus[icusp]<>NULL:
                for jcusp in range(nc):
                    if besv_minus[icusp][jcusp]<>NULL:
                        for n in range(Ml):
                            if besv_minus[icusp][jcusp][n]<>NULL:
                                for j in range(Ql):
                                    mpc_clear(besv_minus[icusp][jcusp][n][j])
                                sig_free(besv_minus[icusp][jcusp][n])
                        sig_free(besv_minus[icusp][jcusp])                        
                sig_free(besv_minus[icusp])
        sig_free(besv_minus)
    mpc_clear(iargpb_t)
    mpc_clear(tmpc_t)
    mpc_clear(tmp1)
    mpc_clear(tmp2)
    mpc_clear(ppc)
    mpc_clear(summa)
    mpc_clear(summa_minus)
    mpfr_clear(tmpr_t)
    mpfr_clear(tmpar)
    mpfr_clear(tmpar1)
    mpfr_clear(tmpab)
    mpfr_clear(tmpcos)
    mpfr_clear(tmpsin)
    mpfr_clear(tmpr)
    mpfr_clear(tmpr2)
    mpfr_clear(Y)
    mpfr_clear(twopi)
    mpfr_clear(fourpi)
    mpfr_clear(twopiY)
    mpfr_clear(fourpiY)
    mpfr_clear(nrfourpi)
    W=dict()
    W['var_a+']=pp_info['variable_a0_plus']
    W['var_a-']=pp_info['variable_a0_minus']
    #if H._verbose>0:
    #    print "alphas=",H.alphas()
    W['alphas']=H.alphas()
    W['V']=V
    #W['cv']=Cv
    W['RHS']=RHS
    W['Ms']=Ms
    W['Mf']=Mf
    W['Ml']=Ml
    W['nc']=nc
    W['PP']=principal_parts
    W['space']=H
    W['rdim']=H._rdim
    return W

cdef void setcossin2(mpfr_t * lcos, mpfr_t * lsin, mpfr_t * Xm, mpfr_t nr, int Ql, mpfr_prec_t prec) nogil:
    cdef int j = 0
    cdef mpfr_t tmpar
    mpfr_init2(tmpar,prec)
    for j in range(Ql):
        mpfr_init2(lcos[j],prec)
        mpfr_init2(lsin[j],prec)
        mpfr_mul(tmpar,Xm[j],nr,rnd_re)
        mpfr_neg(tmpar,tmpar,rnd_re)
        mpfr_cos(lcos[j],tmpar,rnd_re)
        mpfr_sin(lsin[j],tmpar,rnd_re)
    mpfr_clear(tmpar)

cdef void setV(mpc_t **Vmat, mpfr_t ****RCvec,int ***CSvec, mpfr_t **** besv, mpfr_t *** Ypb, mpfr_t ****ef1cosv, mpfr_t ****ef1sinv, mpfr_t ***ef2cosv, mpfr_t ***ef2sinv, int nc, int Ql, int Ml, int l, mpfr_prec_t prec) nogil:
    cdef mpfr_t tmpar,tmpar1,tmpab,tmpcos,tmpsin
    cdef mpc_t tmpc_t
    cdef int jcusp,lj,icusp,j,n, ni
    cdef mpc_t t[2]
    mpc_init2(t[0],prec)
    mpc_init2(t[1],prec)
    mpc_init2(tmpc_t,prec)
    mpfr_init2(tmpcos,prec)
    mpfr_init2(tmpsin,prec)
    mpfr_init2(tmpar,prec)
    mpfr_init2(tmpar1,prec)
    mpfr_init2(tmpab,prec)
    for jcusp in xrange(nc):
        lj=Ml*jcusp+l
        for icusp in xrange(nc):
            for j in xrange(Ql):
                if mpfr_zero_p(Ypb[icusp][jcusp][j])<>0:
                    continue
                #if mpfr_get_si(RCvec[icusp][jcusp][j][2],rnd_re) % 2 == 0:
                #    mpfr_set(
                #mpfr_add(tmpar,ef1[icusp][jcusp][l][j],RCvec[icusp][jcusp][j][1],rnd_re)
                mpfr_mul(tmpab,besv[icusp][jcusp][l][j],RCvec[icusp][jcusp][j][0],rnd_re)
                for n in xrange(Ml): 
                    ni=icusp*Ml+n
                    #mpfr_add(tmpar1,ef2[icusp][n][j],tmpar,rnd_re)
                    #mpc_mul(tmp2.value,tmp1.value,ef2[icusp][n][j],rnd)
                    #printf("ef2cosv[%d][%d][%d]=%f\n",icusp,n,j,mpfr_get_flt(ef2cosv[icusp][n][j],rnd_re))
                    #if mpfr_get_si(RCvec[icusp][jcusp][j][2],rnd_re) % 2 == 0:
                    if CSvec[icusp][jcusp][j] == 0:
                        #mpfr_cos(tmpar1,tmpar1,rnd_re)
                        mpfr_mul(tmpcos,ef2cosv[icusp][n][j],ef1cosv[icusp][jcusp][l][j],rnd_re)
                        mpfr_mul(tmpsin,ef2sinv[icusp][n][j],ef1sinv[icusp][jcusp][l][j],rnd_re)
                        mpfr_sub(tmpar1,tmpcos,tmpsin,rnd_re)
                        mpfr_mul_ui(tmpar1,tmpar1,2,rnd_re)
                        mpc_set_fr(tmpc_t,tmpar1,rnd)                                     
                    else:
                        #mpfr_sin(tmpar1,tmpar1,rnd_re)
                        mpfr_mul(tmpcos,ef2sinv[icusp][n][j],ef1cosv[icusp][jcusp][l][j],rnd_re)
                        mpfr_mul(tmpsin,ef2cosv[icusp][n][j],ef1sinv[icusp][jcusp][l][j],rnd_re)
                        mpfr_add(tmpar1,tmpcos,tmpsin,rnd_re)
                        mpc_set_si_si(tmpc_t,0,2,rnd)
                        _mpc_mul_fr(&tmpc_t,tmpc_t,tmpar1,rnd,rnd_re)             
                    _mpc_mul_fr(&tmpc_t,tmpc_t,tmpab,rnd,rnd_re)                         
                    #tmp2=tmp1*ef2[n,j,icusp]
                    _mpc_add(&Vmat[ni][lj],Vmat[ni][lj],tmpc_t,rnd_re)
                        #V[ni,lj]=V[ni,lj]+tmp1*
                        # if verbose > 1 and ni==0 and lj==10:
                        #     print "-------------------"
                        #     print "V[{0},{1}]({2})={3}".format(ni,lj,j,V[ni,lj]) #CC(V[ni,lj].real(),V[ni,lj].imag())
                        #     mpfr_set(tmpr.value,ef1[icusp][jcusp][n][j],rnd_re)
                        #     print "tmp2=",tmp2 #CC(tmpc.real(),tmpc.imag())
                        #     print "sym=",mpfr_get_si(RCvec[icusp][jcusp][j][2],rnd_re)
                        #     print "ef1(",j,")=",tmpr #CC(tmpc.real(),tmpc.imag())
                        #     #mpc_set(ch.value,Cvec[icusp][jcusp][j],rnd)
                        #     #print "cv(",j,")=",ch #CC(ch.real(),ch.imag())
                        #     mpfr_set(tmpr.value,ef2[icusp][n][j],rnd_re)
                        #     print "ef2(",j,")=",tmpr
    mpc_clear(t[0])
    mpc_clear(t[1])
    mpc_clear(tmpc_t)
    mpfr_clear(tmpcos)
    mpfr_clear(tmpsin)
    mpfr_clear(tmpar)
    mpfr_clear(tmpar1)
    mpfr_clear(tmpab)   
    
### Version to use when we can not use symmetry

cpdef setup_matrix_for_harmonic_Maass_waveforms_no_sym(H,Y_in,int M,int Q,principal_parts,version=1,threads=1):
    r"""

    Set up the matrix for the system of equations giving the Fourier coefficients of a Harmonic Maass waveforms.

    INPUT:
    
        - ``H`` -- Space of harmonic weak Maass forms
        - ``Y`` -- height of horocycle used for sampling (mpmath real)
        - ``k`` -- weight (mpmath real) 
        - ``M`` -- integer
        - ``Q`` -- integer > M
        - ``PP``-- dict : Principal parts at cusp nr. j = \Sum_ PP(j,n) q^[-n]

    OUTPUT:
    
        - ``W`` -- dictionary
        - ``W['Mf']`` -- M start
        - ``W['nc']`` -- number of cusps
        - ``W['Ms']`` -- M stop
        - ``W['V']``  -- ((Ms-Mf+1)*num_cusps)**2 matrix 
        - ``W['RHS']``  -- ((Ms-Mf+1)*num_cusps) matrix
    
    
    EXAMPLES::

        sage: setup_matrix_for_harmonic_Maass_waveforms_sv(MySubgroup(Gamma0(1)),mpmath.mpf(0.5),mpmath.mpf(0),10,20,{(0,-1):1})

    

    """
    cdef int l,i,j,k,icusp,jcusp,n,ni,li,Ml,Ms,Mf,Qs,Qf,Ql,s,nc
    cdef int prec 
    prec = <int>H._prec
    #cdef mpc_t tmpc
    cdef MPComplexNumber tmp1,tmp2,iargpb,iargm,iargpb2
    cdef RealNumber nr,tmpr,nrfourpi,Y,tmpr2
    #cdef RealField_class RF
    #cdef MPComplexField_class CF
    RF = RealField(prec)
    CF = MPComplexField(prec)
    nr = RF(0); tmpr=RF(0); nrfourpi=RF(0); tmpr2=RF(0)
    iargpb=CF(0);iargm=CF(0); tmp1=CF(0); tmp2=CF(0)
    iargpb2=CF(0)
    Y = RF(Y_in)
    weight=RF(H._weight)
    cdef RealNumber pi,one,two,zero,twopi,fourpi,twopiY,fourpiY,p,kint,Qfak,ypb,xpb
    cdef int kinti
    cdef Matrix_complex_dense V
    cdef MPComplexNumber tmpc,ch
    cdef int verbose=int(H._verbose)
    pi=RF.pi() #mpmath_ctx.pi()
    one=RF(1) #mpmath_ctx.mpf(1)
    two=RF(2) #mpmath_ctx.mpf(2)
    zero=RF(0) #mpmath_ctx.mpf(0)
    ypb=RF(0); xpb=RF(0); ch=CF(0); tmpc=CF(0)
    twopi=two*pi
    fourpi=two*twopi
    twopiY=twopi*Y
    fourpiY=two*twopiY
    p=(weight-one)/two
    kint=one-weight
    cdef int is_int=0
    cdef int is_half_int=0
    ## Test if the weight is integral
    if floor(kint)==pceil(kint):
        kinti = int(kint); is_int = 1
    if is_int==0:
        ## Check if kint is half-integral.
        if floor(2*kint)==pceil(2*kint):
            is_half_int = 1
            kinti = int(kint-RF(0.5))
    if verbose>0:
        print "is_int=",is_int
        print "kint=",kint
        print "kinti=",kinti
            
    cdef MPComplexNumber ckint
    ckint = CF(kint)
    #print "kint=",kint
    nc=int(H._group.ncusps())
    if Q<M:
        Q=M+20
    if H._holomorphic or H._almost_holomorphic:
        Ms=0; Mf=M; Ml=Mf-Ms+1
    else:
        Ms=-M; Mf=M; Ml=Mf-Ms+1
    Qs=1-Q; Qf=Q; Ql=Qf-Qs+1
    Qfak=RF(2*Q)
    cdef int do_mpmath = 1
    if hasattr(H,"_do_mpmath"):
         if H._do_mpmath<-1:     
             do_mpmath=0  ## By default we want to use mpmath for now...
             if verbose>0:
                 print 'NOT using mpmath for incgamma'
    cdef int Qfaki
    Qfaki=2*Q
    # Pullpack points
    if verbose>0:
        print "In setup_matrix_for_harmonic_Maass_waveforms_nosym"
        print "Qs,Qf=",Qs,Qf
    cdef mpfr_t tmpr_t
    cdef mpc_t iargpb_t,tmpc_t
    mpfr_init2(tmpr_t,prec)
    mpc_init2(iargpb_t,prec)
    mpc_init2(tmpc_t,prec)
    cdef mpfr_t* Xm
    cdef mpfr_t*** Xpb=NULL
    cdef mpfr_t*** Ypb=NULL
    cdef mpc_t*** Cvec=NULL
    Xm = <mpfr_t*> check_allocarray( sizeof(mpfr_t), Ql )
    for n in range(Ql):
        mpfr_init2(Xm[n],prec)
    Xpb = <mpfr_t***> check_allocarray( sizeof(mpfr_t** ), nc )
    if Xpb==NULL: raise MemoryError
    Ypb = <mpfr_t***> check_allocarray( sizeof(mpfr_t** ), nc )
    if Ypb==NULL: raise MemoryError
    for i in range(nc):
        Xpb[i] = <mpfr_t**>check_allocarray(sizeof(mpfr_t*), nc )
        Ypb[i] = <mpfr_t**>check_allocarray(sizeof(mpfr_t*), nc )
        if Ypb[i]==NULL or Xpb[i]==NULL:
            raise MemoryError
        for j in range(nc):
            Xpb[i][j] = <mpfr_t*>check_allocarray(sizeof(mpfr_t), Ql )
            Ypb[i][j] = <mpfr_t*>check_allocarray(sizeof(mpfr_t), Ql )
            if Ypb[i][j]==NULL or Xpb[i][j]==NULL:
                raise MemoryError
            for n in range(Ql):
                mpfr_init2(Xpb[i][j][n],prec) 
                mpfr_init2(Ypb[i][j][n],prec) 
                mpfr_set_si(Xpb[i][j][n],0,rnd_re)
                mpfr_set_si(Ypb[i][j][n],0,rnd_re)
                #Ypb[i][j][n]=<double>0
    Cvec = <mpc_t***>check_allocarray(sizeof(mpc_t**), nc )
    if Cvec==NULL: raise MemoryError
    for i from 0<=i<nc:                        #if verbose>1:
                        #    mpfr_set(ypb.value,Ypb[icusp][jcusp][j],rnd_re)
                        #    print "ypb[{0}][{1}{2}={3}".format(icusp,jcusp,j,ypb)
                        #    mpc_set(ch.value,Cvec[icusp][jcusp][j],rnd)
                        #    print "ch[{0}][{1}{2}={3}".format(icusp,jcusp,j,ch)

        Cvec[i] = <mpc_t**>check_allocarray(sizeof(mpc_t*), nc )
        if Cvec[i]==NULL:
            raise MemoryError
        for j from 0<=j<nc:
            Cvec[i][j] = <mpc_t*>check_allocarray(sizeof(mpc_t), Ql )
            if Cvec[i][j]==NULL:
                raise MemoryError
            for n from 0<=n<Ql:
                mpc_init2(Cvec[i][j][n],prec)
                mpc_set_si(Cvec[i][j][n],0,rnd)
                #Cvec[i][j][n]=<double complex>0
    pullback_pts_mpc_new_c(H,Qs,Qf,Y.value,Xm,Xpb,Ypb,Cvec)
    #pb=pullback_pts_mp(H,Qs,Qf,Y,weight)
    #Xm=pb['xm']; Xpb=pb['xpb']; Ypb=pb['ypb']; Cv=pb['cvec']
    s=nc*Ml
    MS = MatrixSpace(CF,s,s)
    V = Matrix_complex_dense(MS,0,True,True)
    cdef list PPplus,PPminus
    cdef dict pp_info
    pp_info = check_principal_parts(H,principal_parts)    
    if verbose>0:
        print "pp_info=",pp_info
        log.debug("pp_info={0}".format(pp_info))
    if verbose>1:
        print "Yin=",Y
        mpc_set(ch.value,Cvec[0][0][0],rnd)
        print "Cvec[0,0,0]=",ch
        mpfr_set(tmpr.value,Ypb[0][0][0],rnd_re)
        print "Ypb[0,0,0]=",tmpr
    PPplus = pp_info['PPplus']; PPminus = pp_info['PPminus']
    cdef int d = len(PPplus)
    cdef int **variable_a0_plus
    cdef int **variable_a0_minus
    variable_a0_plus = <int**> check_allocarray(sizeof(int*),d)
    variable_a0_minus = <int**> check_allocarray(sizeof(int*),d)
    for j in range(d):
        variable_a0_plus[j] = <int*> check_allocarray(sizeof(int), nc)
        variable_a0_minus[j] = <int*> check_allocarray(sizeof(int), nc)
        for l in range(nc):
            variable_a0_minus[j][l]=int(pp_info['variable_a0_minus'].get(j,{}).get(l,0))
            variable_a0_plus[j][l]=int(pp_info['variable_a0_plus'].get(j,{}).get(l,0))
    cdef int **PPplus_cusp=NULL
    cdef int **PPplus_n=NULL
    cdef int **PPminus_cusp=NULL
    cdef int **PPminus_n=NULL
    cdef mpc_t **PPplus_values=NULL
    cdef mpc_t **PPminus_values=NULL
    cdef mpfr_t **PPplus_lal=NULL
    cdef mpfr_t **PPminus_lal=NULL
    cdef int num_ppplus = len(pp_info['PPplus'][0])
    cdef int num_ppminus = len(pp_info['PPminus'][0])
    PPplus_cusp = <int **>check_allocarray(sizeof(int*), d)
    if PPplus_cusp==NULL: raise MemoryError
    PPplus_n = <int **>check_allocarray(sizeof(int*), d)
    if PPplus_n==NULL: raise MemoryError
    PPminus_cusp = <int **>check_allocarray(sizeof(int*), d)
    if PPminus_cusp==NULL: raise MemoryError
    PPminus_n = <int **>check_allocarray(sizeof(int*), d)
    if PPminus_n==NULL: raise MemoryError
    PPminus_values = <mpc_t**>check_allocarray(sizeof(mpc_t*), d)
    if PPminus_values==NULL: raise MemoryError
    PPplus_values = <mpc_t**>check_allocarray(sizeof(mpc_t*), d)
    if PPplus_values==NULL: raise MemoryError
    PPplus_lal = <mpfr_t**>check_allocarray(sizeof(mpfr_t*), d)
    if PPplus_lal==NULL: raise MemoryError
    PPminus_lal = <mpfr_t**>check_allocarray(sizeof(mpfr_t*), d)
    if PPminus_lal==NULL: raise MemoryError
    cdef int jj
    for j in range(d):
        PPplus_cusp[j]=NULL;PPplus_n[j]=NULL;PPminus_cusp[j]=NULL;PPminus_n[j]=NULL
        PPplus_values[j]=NULL;PPminus_values[j]=NULL;PPplus_lal[j]=NULL
        PPplus_cusp[j] = <int *>check_allocarray(sizeof(int), num_ppplus)
        if PPplus_cusp[j]==NULL: raise MemoryError
        PPplus_n[j] = <int *>check_allocarray(sizeof(int), num_ppplus)
        if PPplus_n[j]==NULL: raise MemoryError
        PPminus_cusp[j] = <int *>check_allocarray(sizeof(int), num_ppminus)
        if PPminus_cusp[j]==NULL: raise MemoryError
        PPminus_n[j] = <int *>check_allocarray(sizeof(int), num_ppminus)
        if PPplus_n[j]==NULL: raise MemoryError
        PPminus_values[j] = <mpc_t*>check_allocarray(sizeof(mpc_t), num_ppminus)
        if PPminus_values[j]==NULL: raise MemoryError
        PPplus_values[j] = <mpc_t*>check_allocarray(sizeof(mpc_t), num_ppplus)
        if PPplus_values[j]==NULL: raise MemoryError
        PPplus_lal[j] =  <mpfr_t*>check_allocarray(sizeof(mpfr_t), num_ppplus)
        if PPplus_lal[j]==NULL: raise MemoryError
        PPminus_lal[j] =  <mpfr_t*>check_allocarray(sizeof(mpfr_t), num_ppminus)
        if PPminus_lal[j]==NULL: raise MemoryError
        l = 0
        for i,jj in pp_info['PPplus'][j].keys():
            tmpc = CF(pp_info['PPplus'][j][(i,jj)])        
            PPplus_cusp[j][l]=int(i)
            PPplus_n[j][l]=int(jj)
            mpc_init2(PPplus_values[j][l],prec)
            mpc_set(PPplus_values[j][l],tmpc.value,rnd)
            tmpr = RF(jj)+RF(H.alpha(i)[0])
            #print "tmpr=",tmpr
            mpfr_init2(PPplus_lal[j][l],prec)
            mpfr_set(PPplus_lal[j][l],tmpr.value,rnd_re)
            #print "tmpr=",tmpr
            l+=1
        l = 0
        for i,jj in pp_info['PPminus'][j].keys():
            if l>=num_ppminus:
                log.critical("Too large l={0} > num_ppminus={1}!".format(l,num_ppminus))           
            PPminus_cusp[j][l]=int(i)
            PPminus_n[j][l]=int(jj)
            tmpc = CF(pp_info['PPminus'][j][(i,jj)])
            mpc_init2(PPminus_values[j][l],prec)
            mpc_set(PPminus_values[j][l],tmpc.value,rnd)
            tmpr = RF(jj)+RF(H.alpha(i)[0])
            #mpfr_set(tmpr,trn.value,rnd_re)
            mpfr_init2(PPminus_lal[j][l],prec)
            mpfr_set(PPminus_lal[j][l],tmpr.value,rnd_re)            
            l+=1


    cdef int has_key = 0
    MSRHS = MatrixSpace(CF,s,d)
    RHS = Matrix_complex_dense(MSRHS,0,True,True)

    cdef mpfr_t **nvec=NULL
    nvec = <mpfr_t**> check_allocarray( sizeof(mpfr_t* ), nc )
    cdef RealNumber alpha_tmp
    alpha_tmp = RF(0)
    for icusp in range(nc):
        nvec[icusp]=<mpfr_t*> check_allocarray( sizeof(mpfr_t), Ml )
        for l in range(Ml):
            mpfr_init2(nvec[icusp][l],prec)
            alpha_tmp = RF(H.alpha(icusp)[0])
            mpfr_set(tmpr.value,alpha_tmp.value,rnd_re)
            mpfr_set_si(nvec[icusp][l],l+Ms,rnd_re)
            mpfr_add(nvec[icusp][l],nvec[icusp][l],tmpr.value,rnd_re)
    cdef mpc_t ***ef2=NULL
    ef2 = <mpc_t***> check_allocarray( sizeof(mpc_t** ), nc )
    for icusp in range(nc):
        ef2[icusp]=<mpc_t**> check_allocarray( sizeof(mpc_t* ), Ml )
        for n in range(Ml):
            ef2[icusp][n]=<mpc_t*> check_allocarray( sizeof(mpc_t ), Ql )
            
    cdef mpc_t ****ef1=NULL
    #ef1 = <mpc_t****> check_allocarray( sizeof(mpc_t*** ), nc )
    ef1 = <mpc_t****> check_allocarray( sizeof(ef2), nc )
    for icusp in range(nc):
        ef1[icusp]=<mpc_t***> check_allocarray( sizeof(mpc_t** ), nc )
        for jcusp in range(nc):
            ef1[icusp][jcusp]=<mpc_t**> check_allocarray( sizeof(mpc_t* ), Ml )
            for n in range(Ml):
                ef1[icusp][jcusp][n]=<mpc_t*> check_allocarray( sizeof(mpc_t ), Ql )
                
    cdef double eps
    eps = 2.0**float(1-H._dprec)
    cdef int not_holom = int(not H._holomorphic)
    cdef int is_weak = int(H._weak)

    for n in range(Ml): 
        for icusp in range(nc):
            mpfr_set(nr.value,nvec[icusp][n],rnd_re)
            for j in range(Ql):
                mpfr_mul(tmpr.value,Xm[j],nr.value,rnd_re)
                mpfr_neg(tmpr.value,tmpr.value,rnd_re)
                mpc_set_fr(iargm.value,tmpr.value,rnd)
                mpfr_swap(mpc_realref(iargm.value),mpc_imagref(iargm.value))
                mpc_exp(iargm.value,iargm.value,rnd)
                mpc_init2(ef2[icusp][n][j],prec)
                mpc_set(ef2[icusp][n][j],iargm.value,rnd)
        for jcusp in range(nc):
            mpfr_set(nr.value,nvec[jcusp][n],rnd_re)
            mpfr_mul(nrfourpi.value,nr.value,fourpi.value,rnd_re)
            #nrfourpi=nr*fourpi
            for icusp in range(nc):
                for j in range(Ql): 
                    mpc_init2(ef1[icusp][jcusp][n][j],prec)
                    mpc_set_si(ef1[icusp][jcusp][n][j],0,rnd)
                    mpfr_set(ypb.value,Ypb[icusp][jcusp][j],rnd_re)
                    if mpfr_zero_p(ypb.value)<>0:
                        continue
                    mpfr_set(xpb.value,Xpb[icusp][jcusp][j],rnd_re)
                    mpfr_mul(tmpr_t,twopi.value,ypb.value,rnd_re)
                    mpfr_neg(tmpr_t,tmpr_t,rnd_re)
                    mpc_set_fr_fr(iargpb_t,tmpr_t,Xpb[icusp][jcusp][j],rnd)                    
                    mpc_mul_fr(iargpb_t,iargpb_t,nr.value,rnd)
                    mpc_exp(iargpb_t,iargpb_t,rnd)
                    mpc_set(iargpb.value,iargpb_t,rnd)
                    #iargpb=(nr*CF(-twopi*ypb,xpb)).exp()
                    #iargpb=iargpb.exp()

                    if mpfr_cmp_d(nr.value,eps)>0:
                        ## ef1 = e(2pi*i*n(-xpb+i*ypb))
                        mpc_set(ef1[icusp][jcusp][n][j],iargpb.value,rnd) 
                        #ef1[j,icusp,jcusp,n]=one*iargpb
                    #elif mpfr_cmp_d(nr,-eps)<0 and (not H._holomorphic) and H._weak:
                    elif mpfr_cmp_d(nr.value,-eps)<0 and not_holom==1 and is_weak==1:
                        ## ef1 = gamma(1-k,4*pi*|n|)*e(2pi*n*i*(-x+i*ypb))
                        mpfr_abs(tmpr.value,nrfourpi.value,rnd_re)
                        mpfr_mul(tmpr.value,tmpr.value,ypb.value,rnd_re)
                        if (is_int==1 or is_half_int==1) and do_mpmath==0:
                            #print "tmpr=",tmpr                            
                            try:
                                if is_int==1:
                                    if kinti>0:
                                        incgamma_pint_c(tmpr2.value,kinti,tmpr.value)
                                    else:
                                        incgamma_nint_c(tmpr2.value,kinti,tmpr.value)
                                        # print "incgamma_{0}:{1}={2}".format(is_int,kinti,tmpr2)
                                elif is_half_int==1:
                                    incgamma_hint_c(tmpr2.value,kinti,tmpr.value)
                            except ArithmeticError: ## In case we can not achieve the required error
                                tmpr2 = RF(mpmath.mp.gammainc(kint,tmpr).real)
                        else:
                            if verbose>1 and icusp==1 and n==0:
                                 print "arg[{0},{1},{2},{3}]={4}".format(icusp,jcusp,n,j,tmpr)
                    
                            tmpr2 = RF(mpmath.mp.gammainc(kint,tmpr).real)
                        #if verbose>1 and icusp==1 and n==Ms:
                        #    print "Gamma({0},{1})={2}".format(kint,tmpr,tmpr2)
                        if verbose>1 and icusp==1 and n==0:
                            print "f1[{0},{1},{2},{3}]={4}".format(icusp,jcusp,n,j,iargpb)
                            print "be[{0},{1},{2},{3}]={4}".format(icusp,jcusp,n,j,tmpr2)
                            
                        mpc_mul_fr(iargpb.value,iargpb.value,tmpr2.value,rnd)
                        mpc_set(ef1[icusp][jcusp][n][j],iargpb.value,rnd)

                    elif mpfr_cmp_d(nr.value,-eps)<0:
                        mpc_set_si(ef1[icusp][jcusp][n][j],0,rnd)
                    else:
                        ## This is the constant terms.
                        ## Note that principal parts have to be determined
                        ## and will appear in the right hand side
                        #print "n+alpha=0 (jcusp,n)=",jcusp,n
                        if variable_a0_plus[0][jcusp]==1:
                            mpc_set_si(ef1[icusp][jcusp][n][j],1,rnd)
                            #ef1[j,icusp,jcusp,n]=one
                        elif variable_a0_minus[0][jcusp]==1:
                            #mpfr_pow(tmpr.value,ypb.value,kint.value,rnd_re)
                            mpc_set_fr(ef1[icusp][jcusp][n][j],Ypb[icusp][jcusp][j],rnd)
                            if kint==0:
                                mpc_log(ef1[icusp][jcusp][n][j],ef1[icusp][jcusp][n][j],rnd)
                            else:
                                mpc_pow_fr(ef1[icusp][jcusp][n][j],ef1[icusp][jcusp][n][j],kint.value,rnd)
                            #ef1[j,icusp,jcusp,n]=ypb**kint
                        else:
                            ## If none of them are variables this means
                            ## that they are both set in the principal parts
                            mpc_set_si(ef1[icusp][jcusp][n][j],0,rnd)
                            #ef1[j,icusp,jcusp,n]=one
                        #if verbose>1:
                        #print "ef(nr=0)[",j,icusp,jcusp,n,"]=",ef1[j,icusp,jcusp,n]
                    if verbose > 3 and n-Ms==1:
                        mpc_set(tmpc.value,ef1[icusp][jcusp][n][j],rnd)
                        print "ef1[",n,j,icusp,"]=",CC(tmpc.real(),tmpc.imag())
    for n in range(s):
        for l in range(s):        
            mpc_set_ui(V._matrix[n][l],0,rnd)
    for l in range(Ml): 
        for jcusp in range(nc):
            lj=Ml*jcusp+l
            for icusp in range(nc):
                for j in range(Ql): 

                    if mpfr_zero_p(mpc_realref(Cvec[icusp][jcusp][j]))<>0 and mpfr_zero_p(mpc_imagref(Cvec[icusp][jcusp][j]))<>0:
                        continue
                    #if mpfr_zero_p(Ypb[icusp][jcusp][j])<>0:
                    #    continue                        
                    #mpc_set(tmp1.value,rnd)
                    mpc_mul(tmp1.value,ef1[icusp][jcusp][l][j],Cvec[icusp][jcusp][j],rnd)
                    #tmp1=ef1[j,icusp,jcusp,l]*ch
                    for n in range(Ml): 
                        ni=Ml*icusp+n
                        mpc_mul(tmp2.value,tmp1.value,ef2[icusp][n][j],rnd)
                        #tmp2=tmp1*ef2[n,j,icusp]
                        mpc_add(V._matrix[ni][lj],V._matrix[ni][lj],tmp2.value,rnd)
                        #V[ni,lj]=V[ni,lj]+tmp1*
                        if verbose > 1 and ni==0 and lj==10:
                            print "-------------------"
                            #print "V[1,1](",j,")=",V[ni,lj] #CC(V[ni,lj].real(),V[ni,lj].imag())
                            print "V[{0},{1}]({2})={3}".format(ni,lj,j,V[ni,lj])
                            print "tmp1=",tmp1
                            print "tmp2=",tmp2
                            mpc_set(tmpc.value,ef1[icusp][jcusp][n][j],rnd)
                            print "ef1(",j,")=",tmpc #CC(tmpc.real(),tmpc.imag())
                            mpc_set(ch.value,Cvec[icusp][jcusp][j],rnd)
                            print "cv(",j,")=",ch #CC(ch.real(),ch.imag())
                            mpc_set(tmpc.value,ef2[icusp][n][j],rnd)
                            print "ef2(",j,")=",tmpc
                        
    if verbose>1:
        print "Mi,Ms=",Ms,Mf
        print "V1(",0,0,")=",V[0,0]
        print "V1(",1,0,")=",V[1,0]
        print "V1(",0,1,")=",V[0,1]
        print "V1(",1,1,")=",V[1,1]
        if V.ncols()>=11:
            print "V1(",0,11,")=",V[0,11]
    cdef int nrows,ncols
    nrows = int(V.nrows()); ncols = int(V.ncols())
    for n in range(nrows):
        for l in range(ncols):        
            mpc_div_ui(V._matrix[n][l],V._matrix[n][l],Qfaki,rnd)
            #V[n,l]=V[n,l]/Qfak
    if verbose>1: 
        print "V1(",0,0,")=",V[0,0]
        print "V1(",1,0,")=",V[1,0]            
        print "V1(",0,1,")=",V[0,1]
        print "V1(",1,1,")=",V[1,1]
        if V.ncols()>=11:
            print "V1(",0,11,")=",V[0,11]
    cdef MPComplexNumber f1,f2,ppc,summa,ppc_minus,summa_minus
    cdef RealNumber lr,arg,nrY2pi,kbes
    lr = RF(0); arg=RF(0);nrY2pi=RF(0); kbes=RF(0)
    f1 = CF(0); f2=CF(0); ppc=CF(0); summa=CF(0)
    ppc_minus = CF(0); summa_minus=CF(0)
    for n in range(Ml):
        for icusp in range(nc):
            mpfr_set(nr.value,nvec[icusp][n],rnd_re)
            mpfr_mul(nrY2pi.value,nr.value,twopiY.value,rnd_re)
            #nrY2pi=nr*twopiY
            ni=Ml*icusp+n
            if mpfr_cmp_d(nr.value,eps)>0:
                mpfr_neg(kbes.value,nrY2pi.value,rnd_re)
                mpfr_exp(kbes.value,kbes.value,rnd_re)
                #kbes=(-nrY2pi).exp()
            elif mpfr_cmp_d(nr.value,-eps)<0 and not_holom==1 and is_weak==1:
                #kbes=RF(mpmath.mp.gammainc(kint,abs(nr)*fourpiY))
                #kbes=ckint.gamma_inc(abs(nr)*fourpiY).real()                
                mpfr_mul(tmpr.value,nr.value,fourpiY.value,rnd_re)
                mpfr_abs(tmpr.value,tmpr.value,rnd_re)
                if is_int==1 or is_half_int==1 and  do_mpmath==0:
                    try:
                        if is_int==1:
                            if kinti>0:
                                incgamma_pint_c(kbes.value,kinti,tmpr.value)
                            else:
                                incgamma_nint_c(kbes.value,kinti,tmpr.value)
                        elif is_half_int==1:
                            incgamma_hint_c(kbes.value,kinti,tmpr.value)                                
                    except ArithmeticError: ## In case we can not achieve the required error
                        kbes = RF(mpmath.mp.gammainc(kint,tmpr).real)                    
                else:
                    kbes = RF(mpmath.mp.gammainc(kint,tmpr).real)

                mpfr_abs(tmpr2.value,nrY2pi.value,rnd_re)
                mpfr_exp(tmpr2.value,tmpr2.value,rnd_re)
                if verbose>1 and ni==0:
                    print "Arg=",tmpr
                    print "Igamma(",kinti,",Arg)=",kbes
                mpfr_mul(kbes.value,kbes.value,tmpr2.value,rnd_re)
                if verbose>1 and ni==0:
                    print "expfax=",tmpr2
                    print "Igamma(",kinti,",Arg)*exp()=",kbes
                #kbes=kbes*(-nrY2pi).exp()
            elif mpfr_cmp_d(nr.value,-eps)<0:
                mpfr_set_ui(kbes.value,0,rnd_re) #=zero
            else:
                if variable_a0_plus[0][icusp]==1:
                    mpfr_set_ui(kbes.value,1,rnd_re) #=one
                elif variable_a0_minus[0][icusp]==1:
                    if kint<>0:
                        mpfr_pow(kbes.value,Y.value,kint.value,rnd_re)
                    else:
                        mpfr_log(kbes.value,Y.value,rnd_re)
                    #kbes=Y**kint
                else:
                    #kbes=one
                    #raise ValueError,"Need to specify at least one of the n=0 terms!"
                    mpfr_set_ui(kbes.value,0,rnd_re)
                    #kbes=zero
            mpc_sub_fr(V._matrix[ni][ni],V._matrix[ni][ni],kbes.value,rnd)
            if verbose>1:

                if ni==0:
                    print "Arg=",tmpr
                    print "kbes(",0,")=",kbes
                    print "V[1,1]=",V[0,0]
                if ni==1:
                    print "nr=",nr
                    print "kint,kinti=",kint,kinti
                    print "Arg=",tmpr
                    print "kbes(",1,")=",kbes
                    print "V[1,1]=",V[1,1]
                if ni == 11 and V.ncols()>=11:
                    print "V1(",0,11,")=",V[0,11]
            # setting the right hand side of the system
            for k in range(d):
                mpc_set_ui(RHS._matrix[ni][k],0,rnd) #=CF(0)
                for i in range(num_ppplus):
                #for (jcusp,l) in PPplus[k].keys():
                    mpc_set(ppc.value,PPplus_values[k][i],rnd)
                    jcusp = PPplus_cusp[k][i]
                    l = PPplus_n[k][i]
                    #ppc=CF(PPplus[k][(jcusp,l)])
                    #mpfr_set(lr.value,nvec[jcusp][l],rnd_re)
                    #lr=RF(l+H.alpha(jcusp)[0])
                    #if(abs(lr) <= mpmath.eps()):
                    #    continue
                    ### Note: The constant term is treated below
                    if mpc_is_zero(ppc.value)==1: #==zero or PPplus[k][(jcusp,l)]==0:
                        continue
                    mpc_set_ui(summa.value,0,rnd) #=CF(0)
                    if verbose>2:
                        mpfr_set(lr.value,PPplus_lal[k][i],rnd_re)
                        print "l=",lr
                        print "icusp,jcusp=",icusp,jcusp
                        print "alpha(",jcusp,")=",H.alpha(jcusp)[0]
                        print "n=",n,nr," n=",n+Ms
                        print "ppc=",ppc
                    for j in range(Ql):   
                        if mpc_is_zero(Cvec[icusp][jcusp][j])==1:
                            continue
                        if mpfr_zero_p(Ypb[icusp][jcusp][j])<>0:
                            continue
                        mpfr_mul(tmpr_t,Ypb[icusp][jcusp][j],twopi.value,rnd_re)
                        if verbose>3:
                            mpfr_set(tmpr.value,Ypb[icusp][jcusp][j],rnd_re)
                            print "ypb=",tmpr
                            print "lr=",lr
                            print "twopi=",twopi
                        mpfr_neg(tmpr_t,tmpr_t,rnd_re)
                        mpfr_mul(tmpr_t,tmpr_t,PPplus_lal[k][i],rnd_re)
                        mpfr_exp(tmpr_t,tmpr_t,rnd_re)
                        # tmpr_t = exp(-2pi*l*ypb)
                        if verbose>2:
                            mpfr_set(tmpr.value,tmpr_t,rnd_re)
                            print "-2pi*l*ypb=",tmpr
                        mpc_set_fr(f1.value,tmpr_t,rnd)
                        mpfr_set(xpb.value,Xpb[icusp][jcusp][j],rnd_re)
                        #arg = lr*xpb-nr*Xm[j]
                        mpfr_mul(tmpr_t,PPplus_lal[k][i],Xpb[icusp][jcusp][j],rnd_re)
                        mpfr_mul(arg.value,nr.value,Xm[j],rnd_re)
                        mpfr_sub(arg.value,tmpr_t,arg.value,rnd_re)
                        mpc_set_fr_fr(iargpb_t,zero.value,arg.value,rnd)
                        mpc_exp(f2.value,iargpb_t,rnd)
                        #f2 = CF(0,arg).exp()
                        if verbose>3:
                            print "f1(",j,")=",f1
                            print "arg=",lr,"*",xpb,"-",nr,"*",mpfr_get_d(Xm[j],rnd_re),"=",arg
                            print "f2(",j,")=",f2
                        mpc_mul(tmpc.value,f1.value,f2.value,rnd)
                        mpc_mul(tmpc.value,tmpc.value,ppc.value,rnd)
                        mpc_mul(tmpc.value,tmpc.value,Cvec[icusp][jcusp][j],rnd)
                        mpc_add(summa.value,summa.value,tmpc.value,rnd)
                        if verbose>3:
                            print "tmpc=",tmpc
                            print "summa=",summa                 
#                        summa=summa+ch*tmpc
                    if verbose>2:
                        print "summa=",RR(summa.real()),RR(summa.imag())
                    #RHS[ni,k]=RHS[ni,k]+summa/Qfak
                    mpc_div_ui(summa.value,summa.value,Qfaki,rnd)
                    mpc_add(RHS._matrix[ni][k],RHS._matrix[ni][k],summa.value,rnd)
                if verbose>2:
                    print "RHS0[",ni,k,"]=",RHS[ni][k]
                    print "icusp,n+Ms=",icusp,n+Ms
                    print "PPplus.keys=",PPplus[k].keys()
                has_key = 0
                for j in range(num_ppplus):
                    if PPplus_cusp[k][j]==icusp and PPplus_n[k][j]==n+Ms:
                        has_key = 1
                        mpc_set(ppc.value,PPplus_values[k][j],rnd)
                        break
                if has_key==1:
                    #if PPplus[k].has_key((icusp,n+Ms)):
                    #if( abs(nr) > mpmath.eps()):
                    #ppc=CF(PPplus[k][icusp,n+Ms])
                    mpfr_abs(tmpr_t,nrY2pi.value,rnd_re)
                    mpfr_exp(tmpr_t,tmpr_t,rnd_re)
                    mpc_mul_fr(tmpc_t,ppc.value,tmpr_t,rnd)
                    mpc_sub(RHS._matrix[ni][k],RHS._matrix[ni][k],tmpc_t,rnd)
                    if verbose>2:
                        print "n=",n
                        print "icusp=",icusp
                        print "ppc=",ppc
                        print "nrY2pi=",nrY2pi
                    #RHS[ni,k]=RHS[ni,k]-ppc*(-nrY2pi).exp()
                if verbose>2:
                    print "RHS1[",ni,k,"]=",RHS[ni][k]
                for i in range(num_ppminus):
                    #for (jcusp,l) in PPminus[k].keys():
                    jcusp = PPminus_cusp[k][i]
                    l = PPminus_n[k][i]
                    mpc_set(ppc_minus.value,PPminus_values[k][i],rnd)
                    #ppc_minus = CF(PPminus[k][(jcusp,l)])
                    #if ppc_minus==zero or PPminus[k][(jcusp,l)]==0:
                    if mpc_is_zero(ppc_minus.value)==1:
                        continue
                    if verbose>2:
                        print "ppart_minus=",ppc_minus
                    mpc_set_ui(summa_minus.value,0,rnd) #=zero #summa_minus
                    for j in range(Ql):
                        #if mpfr_zero_p(mpc_realref(Cvec[icusp][jcusp][j]))<>0 and mpfr_zero_p(mpc_imagref(Cvec[icusp][jcusp][j]))<>0:
                        if mpc_is_zero(Cvec[icusp][jcusp][j])==1:
                            continue
                        if mpfr_zero_p(Ypb[icusp][jcusp][j])<>0:
                            continue
                        mpc_set(ch.value,Cvec[icusp][jcusp][j],rnd)
                        mpfr_set(ypb.value,Ypb[icusp][jcusp][j],rnd_re)
                        if kint==0:
                            mpfr_log(tmpr.value,ypb.value,rnd_re)
                        else:
                            mpfr_pow(tmpr.value,ypb.value,kint.value,rnd_re)
                        mpc_set_fr(tmpc_t,tmpr.value,rnd) #ypb**kint

                        mpfr_mul(tmpr.value,nr.value,Xm[j],rnd_re)
                        mpfr_neg(tmpr.value,tmpr.value,rnd_re)
                        mpc_set_fr(iargm.value,tmpr.value,rnd)
                        mpfr_swap(mpc_realref(iargm.value),mpc_imagref(iargm.value))
                        mpc_exp(iargm.value,iargm.value,rnd)
                        mpc_mul(tmpc.value,tmpc_t,iargm.value,rnd)
                        #tmpc= tmpc_minus*CF(0,-nr*Xm[j]).exp()
                        # tmpc = ypb^{k}e(-2*pi*i*Xm[j])
                        mpc_mul(tmpc.value,tmpc.value,ch.value,rnd)
                        mpc_add(summa_minus.value,summa_minus.value,tmpc.value,rnd)
                        #summa_minus=summa_minus+ch*tmpc
                        #summa=summa+ef2[n(2):,j,icusp]*tmpc*ch*ppc
                        if verbose>2:
                            print "j=",j
                            print "ypb=",ypb
                            print "ch=",ch
                            #print "tmpc_minus=",tmpc_minus
                            print "tmpc=",tmpc
                            print "summa-(",j,")=",summa_minus
                    if verbose>2:
                        print "Summa-(",ni,")=",summa_minus
                    mpc_div_ui(summa_minus.value,summa_minus.value,Qfaki,rnd)
                    mpc_mul(summa_minus.value,summa_minus.value,ppc_minus.value,rnd)
                    mpc_add(RHS._matrix[ni][k],RHS._matrix[ni][k],summa_minus.value,rnd)
                    #RHS[ni,k]=RHS[ni,k]+ppc_minus*summa_minus/Qfak
                #print "summa(",ni,")=",RHS[ni,k]
                if verbose>2:
                    print "RHS2[{0},{1}]={2}".format(ni,k,RHS[ni,k])
                has_key = 0
                for j in range(num_ppminus):
                    if PPminus_cusp[k][j]==icusp and PPminus_n[k][j]==n+Ms:
                        has_key = 1
                        mpc_set(ppc_minus.value,PPminus_values[k][j],rnd)
                        break
                if has_key==1:
                    #if PPminus[k].has_key((icusp,n+Ms)):
                    #ppc_minus = CF(PPminus[k][(icusp,n+Ms)])
                    if kinti==0:
                        mpfr_log(tmpr.value,Y.value,rnd_re)
                    else:
                        mpfr_pow(tmpr.value,Y.value,kint.value,rnd_re)
                    mpc_mul_fr(ppc_minus.value,ppc_minus.value,tmpr.value,rnd)
                    mpc_sub(RHS._matrix[ni][k],RHS._matrix[ni][k],ppc_minus.value,rnd)
                    #RHS[ni,k]=RHS[ni,k]-ppc_minus*Y**kint
                    #print "subtracting:",ppc_minus*mpmath.power(Y,kint)
                if verbose>2:
                    print "RHS3[{0},{1}]={2}".format(ni,k,RHS[ni,k])
                    #print "RHS(",ni,")-=",ppc*mpmath_ctx.exp(-nrY2pi)
    if verbose>1: 
        for n in range(2):
            print "V2(",n,n,")=",V[n,n]
        if V.ncols()>=11:
            print "V2(",0,11,")=",V[0,11]
            # Clearing up allocated variables

    for j in range(d):
        sig_free(variable_a0_minus[j])
        sig_free(variable_a0_plus[j])        
        sig_free(PPplus_cusp[j])
        sig_free(PPplus_n[j])
        sig_free(PPminus_cusp[j])
        sig_free(PPminus_n[j])
        for l in range(num_ppplus):
            mpc_clear(PPplus_values[j][l])
            mpfr_clear(PPplus_lal[j][l])
        sig_free(PPplus_lal[j])
        sig_free(PPplus_values[j])
        for l in range(num_ppminus):
            mpc_clear(PPminus_values[j][l])
        sig_free(PPminus_values[j])
        
    sig_free(PPplus_lal)
    sig_free(PPplus_cusp)
    sig_free(PPplus_n)
    sig_free(PPplus_values)
    sig_free(PPminus_cusp)
    sig_free(PPminus_n)
    sig_free(PPminus_values)

    sig_free(variable_a0_minus)
    sig_free(variable_a0_plus)




    if Ypb<>NULL:
        for i in range(nc):
            if Ypb[i]<>NULL:
                for j in range(nc):
                    if Ypb[i][j]<>NULL:
                        sig_free(Ypb[i][j])
                sig_free(Ypb[i])
        sig_free(Ypb)
    if Xpb<>NULL:
        for i in range(nc):
            if Xpb[i]<>NULL:
                for j in range(nc):
                    if Xpb[i][j]<>NULL:
                        sig_free(Xpb[i][j])
                sig_free(Xpb[i])
        sig_free(Xpb)
    if Cvec<>NULL:
        for i in range(nc):
            if Cvec[i]<>NULL:
                for j in range(nc):
                    if Cvec[i][j]<>NULL:
                        sig_free(Cvec[i][j])
                sig_free(Cvec[i])
        sig_free(Cvec)
    if nvec<>NULL:
        for i in range(nc):
            if nvec[i]<>NULL:
                for l in range(Ml):
                    mpfr_clear(nvec[i][l])
                sig_free(nvec[i])
        sig_free(nvec)
    if Xm<>NULL:
        for n in range(Ql):
            mpfr_clear(Xm[n])
        sig_free(Xm)
    if ef2<>NULL:
        for icusp in range(nc):
            if ef2[icusp]<>NULL:
                for n in range(Ml):
                    if ef2[icusp][n]<>NULL:
                        for j in range(Ql):
                            mpc_clear(ef2[icusp][n][j])
                        sig_free(ef2[icusp][n])
                sig_free(ef2[icusp])
        sig_free(ef2)
    if ef1<>NULL:
        for icusp in range(nc):
            if ef1[icusp]<>NULL:
                for jcusp in range(nc):
                    if ef1[icusp][jcusp]<>NULL:
                        for n in range(Ml):
                            if ef1[icusp][jcusp][n]<>NULL:
                                for j in range(Ql):
                                    mpc_clear(ef1[icusp][jcusp][n][j])
                                sig_free(ef1[icusp][jcusp][n])
                        sig_free(ef1[icusp][jcusp])                        
                sig_free(ef1[icusp])
        sig_free(ef1)
        
        #    print "V2(",44,n,")=",V[V.rows-1,n]
    mpc_clear(iargpb_t)
    mpc_clear(tmpc_t)
    mpfr_clear(tmpr_t)
    W=dict()
    W['var_a+']=pp_info['variable_a0_plus']
    W['var_a-']=pp_info['variable_a0_minus']
    #if H._verbose>0:
    #    print "alphas=",H.alphas()
    W['alphas']=H.alphas()
    W['V']=V
    #W['cv']=Cv
    W['RHS']=RHS
    W['Ms']=Ms
    W['Mf']=Mf
    W['Ml']=Ml
    W['nc']=nc
    W['PP']=principal_parts
    W['space']=H
    W['rdim']=H._rdim
    return W

## cpdef setup_matrix_for_harmonic_Maass_waveforms_sv_orig(H,Y_in,int M,int Q,principal_parts,version=1):
##     r"""

##     Set up the matrix for the system of equations giving the Fourier coefficients of a Harmonic Maass waveforms.

##     INPUT:
    
##         - ``H`` -- Space of harmonic weak Maass forms
##         - ``Y`` -- height of horocycle used for sampling (mpmath real)
##         - ``k`` -- weight (mpmath real) 
##         - ``M`` -- integer
##         - ``Q`` -- integer > M
##         - ``PP``-- dict : Principal parts at cusp nr. j = \Sum_ PP(j,n) q^[-n]

##     OUTPUT:
    
##         - ``W`` -- dictionary
##         - ``W['Mf']`` -- M start
##         - ``W['nc']`` -- number of cusps
##         - ``W['Ms']`` -- M stop
##         - ``W['V']``  -- ((Ms-Mf+1)*num_cusps)**2 matrix 
##         - ``W['RHS']``  -- ((Ms-Mf+1)*num_cusps) matrix
    
    
##     EXAMPLES::

##         sage: setup_matrix_for_harmonic_Maass_waveforms_sv(MySubgroup(Gamma0(1)),mpmath.mpf(0.5),mpmath.mpf(0),10,20,{(0,-1):1})

    

##     """
##     cdef int l,i,j,k,icusp,jcusp,n,ni,li,Ml,Ms,Mf,Qs,Qf,Ql,s,nc
##     cdef int prec 
##     prec = <int>H._prec
##     #cdef mpc_t tmpc
##     cdef MPComplexNumber tmp1,tmp2,iargpb,iargm,iargpb2
##     cdef RealNumber nr,tmpr,nrfourpi,Y,tmpr2
##     cdef RealField_class RF
##     cdef MPComplexField_class CF
##     RF = RealField_class(prec)
##     CF = MPComplexField_class(prec)
##     nr = RF(0); tmpr=RF(0); nrfourpi=RF(0); tmpr2=RF(0)
##     iargpb=CF(0);iargm=CF(0); tmp1=CF(0); tmp2=CF(0)
##     iargpb2=CF(0)
##     Y = RF(Y_in)
##     weight=RF(H._weight)
##     cdef RealNumber pi,one,two,zero,twopi,fourpi,twopiY,fourpiY,p,kint,Qfak,ypb,xpb
##     cdef int kinti
##     cdef Matrix_complex_dense V
##     cdef MPComplexNumber tmpc,ch
##     cdef int verbose=int(H._verbose)
##     pi=RF.pi() #mpmath_ctx.pi()
##     one=RF(1) #mpmath_ctx.mpf(1)
##     two=RF(2) #mpmath_ctx.mpf(2)
##     zero=RF(0) #mpmath_ctx.mpf(0)
##     ypb=RF(0); xpb=RF(0); ch=CF(0); tmpc=CF(0)
##     twopi=two*pi
##     fourpi=two*twopi
##     twopiY=twopi*Y
##     fourpiY=two*twopiY
##     p=(weight-one)/two
##     kint=one-weight
##     cdef int is_int=0
##     cdef int is_half_int=0
##     ## Test if the weight is integral
##     if floor(kint)==pceil(kint):
##         kinti = int(kint); is_int = 1
##     if verbose>0:
##         print "is_int=",is_int
##         print "kint=",kint
##         print "kinti=",kinti
##     if is_int==0:
##         ## Check if kint is half-integral.
##         if floor(2*kint)==pceil(2*kint):
##             is_half_int = 1
##             kinti = int(kint-RF(0.5))
##     cdef MPComplexNumber ckint
##     ckint = CF(kint)
##     #print "kint=",kint
##     nc=int(H._group.ncusps())
##     if Q<M:
##         Q=M+20
##     if H._holomorphic:
##         Ms=0; Mf=M; Ml=Mf-Ms+1
##     else:
##         Ms=-M; Mf=M; Ml=Mf-Ms+1
##     Qs=1-Q; Qf=Q; Ql=Qf-Qs+1
##     Qfak=RF(2*Q)
##     cdef int do_mpmath = 0
##     if hasattr(H,"_do_mpmath"):
##         do_mpmath=H._do_mpmath
##     cdef int Qfaki
##     Qfaki=2*Q
##     # Pullpack points
##     if verbose>0:
##         print "In setup_matrix_for_harmonic_Maass_waveforms_sv_orig"
##         print "Qs,Qf=",Qs,Qf
##     cdef mpfr_t tmpr_t
##     cdef mpc_t iargpb_t,tmpc_t
##     mpfr_init2(tmpr_t,prec)
##     mpc_init2(iargpb_t,prec)
##     mpc_init2(tmpc_t,prec)
##     cdef Vector_real_mpfr_dense Xm
##     cdef mpfr_t*** Xpb=NULL
##     cdef mpfr_t*** Ypb=NULL
##     cdef mpc_t*** Cvec=NULL
##     Xm=Vector_real_mpfr_dense(vector(RF,Ql).parent(),0)
##     Xpb = <mpfr_t***> check_allocarray( sizeof(mpfr_t** ) * nc )
##     if Xpb==NULL: raise MemoryError
##     Ypb = <mpfr_t***> check_allocarray( sizeof(mpfr_t** ), nc )
##     if Ypb==NULL: raise MemoryError
##     for i in range(nc):
##         Xpb[i] = <mpfr_t**>check_allocarray(sizeof(mpfr_t*), nc )
##         Ypb[i] = <mpfr_t**>check_allocarray(sizeof(mpfr_t*), nc )
##         if Ypb[i]==NULL or Xpb[i]==NULL:
##             raise MemoryError
##         for j in range(nc):
##             Xpb[i][j] = <mpfr_t*>check_allocarray(sizeof(mpfr_t), Ql )
##             Ypb[i][j] = <mpfr_t*>check_allocarray(sizeof(mpfr_t), Ql )
##             if Ypb[i][j]==NULL or Xpb[i][j]==NULL:
##                 raise MemoryError
##             for n in range(Ql):
##                 mpfr_init2(Xpb[i][j][n],prec) 
##                 mpfr_init2(Ypb[i][j][n],prec) 
##                 mpfr_set_si(Xpb[i][j][n],0,rnd_re)
##                 mpfr_set_si(Ypb[i][j][n],0,rnd_re)
##                 #Ypb[i][j][n]=<double>0
##     Cvec = <mpc_t***>check_allocarray(sizeof(mpc_t**), nc )
##     if Cvec==NULL: raise MemoryError
##     for i from 0<=i<nc:                        #if verbose>1:
##                         #    mpfr_set(ypb.value,Ypb[icusp][jcusp][j],rnd_re)
##                         #    print "ypb[{0}][{1}{2}={3}".format(icusp,jcusp,j,ypb)
##                         #    mpc_set(ch.value,Cvec[icusp][jcusp][j],rnd)
##                         #    print "ch[{0}][{1}{2}={3}".format(icusp,jcusp,j,ch)

##         Cvec[i] = <mpc_t**>check_allocarray(sizeof(mpc_t*), nc )
##         if Cvec[i]==NULL:
##             raise MemoryError
##         for j from 0<=j<nc:
##             Cvec[i][j] = <mpc_t*>check_allocarray(sizeof(mpc_t), Ql )
##             if Cvec[i][j]==NULL:
##                 raise MemoryError
##             for n from 0<=n<Ql:
##                 mpc_init2(Cvec[i][j][n],prec)
##                 mpc_set_si(Cvec[i][j][n],0,rnd)
##                 #Cvec[i][j][n]=<double complex>0
##     pullback_pts_mpc_new_c(H,Qs,Qf,Y,Xm,Xpb,Ypb,Cvec)
##     #pb=pullback_pts_mp(H,Qs,Qf,Y,weight)
##     #Xm=pb['xm']; Xpb=pb['xpb']; Ypb=pb['ypb']; Cv=pb['cvec']
##     s=nc*Ml
##     MS = MatrixSpace(CF,s,s)
##     V = Matrix_complex_dense(MS,0,True,True)
##     cdef list PPplus,PPminus
##     cdef dict pp_info
##     pp_info = check_principal_parts(H,principal_parts)    
##     if verbose>0:
##         print "pp_info=",pp_info
##     PPplus = pp_info['PPplus']; PPminus = pp_info['PPminus']
##     cdef int d = len(PPplus)
##     cdef int **variable_a0_plus,**variable_a0_minus
##     variable_a0_plus = <int**> check_allocarray(d*sizeof(int*))
##     variable_a0_minus = <int**> check_allocarray(d*sizeof(int*))
##     for j in range(d):
##         variable_a0_plus[j] = <int*> check_allocarray(nc*sizeof(int))
##         variable_a0_minus[j] = <int*> check_allocarray(nc*sizeof(int))
##         for l in range(nc):
##             variable_a0_minus[j][l]=int(pp_info['variable_a0_minus'][j][l])
##             variable_a0_plus[j][l]=int(pp_info['variable_a0_plus'][j][l])
##     cdef int **PPplus_cusp=NULL , **PPplus_n=NULL,**PPminus_cusp=NULL , **PPminus_n=NULL
##     cdef mpc_t **PPplus_values=NULL,**PPminus_values=NULL
##     cdef mpfr_t **PPplus_lal=NULL
##     cdef int num_ppplus = len(pp_info['PPplus'][0])
##     cdef int num_ppminus = len(pp_info['PPminus'][0])
##     PPplus_cusp = <int **>check_allocarray(d*sizeof(int*))
##     if PPplus_cusp==NULL: raise MemoryError
##     PPplus_n = <int **>check_allocarray(d*sizeof(int*))
##     if PPplus_n==NULL: raise MemoryError
##     PPminus_cusp = <int **>check_allocarray(d*sizeof(int*))
##     if PPminus_cusp==NULL: raise MemoryError
##     PPminus_n = <int **>check_allocarray(d*sizeof(int*))
##     if PPminus_n==NULL: raise MemoryError
##     PPminus_values = <mpc_t**>check_allocarray(d*sizeof(mpc_t*))
##     if PPminus_values==NULL: raise MemoryError
##     PPplus_values = <mpc_t**>check_allocarray(d*sizeof(mpc_t*))
##     if PPplus_values==NULL: raise MemoryError
##     PPplus_lal = <mpfr_t**>check_allocarray(d*sizeof(mpfr_t*))
##     if PPplus_lal==NULL: raise MemoryError
##     cdef int jj
##     for j in range(d):
##         PPplus_cusp[j]=NULL;PPplus_n[j]=NULL;PPminus_cusp[j]=NULL;PPminus_n[j]=NULL
##         PPplus_values[j]=NULL;PPminus_values[j]=NULL;PPplus_lal[j]=NULL
##         PPplus_cusp[j] = <int *>check_allocarray(num_ppplus*sizeof(int))
##         if PPplus_cusp[j]==NULL: raise MemoryError
##         PPplus_n[j] = <int *>check_allocarray(num_ppplus*sizeof(int))
##         if PPplus_n[j]==NULL: raise MemoryError
##         PPminus_cusp[j] = <int *>check_allocarray(num_ppminus*sizeof(int))
##         if PPminus_cusp[j]==NULL: raise MemoryError
##         PPminus_n[j] = <int *>check_allocarray(num_ppminus*sizeof(int))
##         if PPplus_n[j]==NULL: raise MemoryError
##         PPminus_values[j] = <mpc_t*>check_allocarray(num_ppminus*sizeof(mpc_t))
##         if PPminus_values[j]==NULL: raise MemoryError
##         PPplus_values[j] = <mpc_t*>check_allocarray(num_ppplus*sizeof(mpc_t))
##         if PPplus_values[j]==NULL: raise MemoryError
##         PPplus_lal[j] =  <mpfr_t*>check_allocarray(num_ppplus*sizeof(mpfr_t))
##         if PPplus_lal[j]==NULL: raise MemoryError
##         l = 0
##         for i,jj in pp_info['PPplus'][j].keys():
##             tmpc = CF(pp_info['PPplus'][j][(i,jj)])        
##             PPplus_cusp[j][l]=int(i)
##             PPplus_n[j][l]=int(jj)
##             mpc_init2(PPplus_values[j][l],prec)
##             mpc_set(PPplus_values[j][l],tmpc.value,rnd)
##             tmpr = RF(jj)+RF(H.alpha(i)[0])
##             #print "tmpr=",tmpr
##             mpfr_init2(PPplus_lal[j][l],prec)
##             mpfr_set(PPplus_lal[j][l],tmpr.value,rnd_re)
##             #print "tmpr=",tmpr
##             l+=1
##         l = 0
##         for i,jj in pp_info['PPminus'][j].keys():
##             PPminus_cusp[j][l]=int(i)
##             PPminus_n[j][l]=int(jj)
##             tmpc = CF(pp_info['PPminus'][j][(i,jj)])
##             mpc_init2(PPminus_values[j][l],prec)
##             mpc_set(PPminus_values[j][l],tmpc.value,rnd)
##             l+=1


##     cdef int has_key = 0
##     MSRHS = MatrixSpace(CF,s,d)
##     RHS = Matrix_complex_dense(MSRHS,0,True,True)

##     cdef mpfr_t **nvec=NULL
##     nvec = <mpfr_t**> check_allocarray( sizeof(mpfr_t* ), nc )
##     cdef RealNumber alpha_tmp
##     alpha_tmp = RF(0)
##     for icusp in range(nc):
##         nvec[icusp]=<mpfr_t*> check_allocarray( sizeof(mpfr_t), Ml )
##         for l in range(Ml):
##             mpfr_init2(nvec[icusp][l],prec)
##             alpha_tmp = RF(H.alpha(icusp)[0])
##             mpfr_set(tmpr.value,alpha_tmp.value,rnd_re)
##             mpfr_set_si(nvec[icusp][l],l+Ms,rnd_re)
##             mpfr_add(nvec[icusp][l],nvec[icusp][l],tmpr.value,rnd_re)
##     cdef mpc_t ***ef2=NULL
##     ef2 = <mpc_t***> check_allocarray( sizeof(mpc_t** ), nc )
##     for icusp in range(nc):
##         ef2[icusp]=<mpc_t**> check_allocarray( sizeof(mpc_t* ), Ml )
##         for n in range(Ml):
##             ef2[icusp][n]=<mpc_t*> check_allocarray( sizeof(mpc_t ), Ql )
            
##     cdef mpc_t ****ef1=NULL
##     #ef1 = <mpc_t****> check_allocarray( sizeof(mpc_t*** ), nc )
##     ef1 = <mpc_t****> check_allocarray( sizeof(ef2), nc )
##     for icusp in range(nc):
##         ef1[icusp]=<mpc_t***> check_allocarray( sizeof(mpc_t** ), nc )
##         for jcusp in range(nc):
##             ef1[icusp][jcusp]=<mpc_t**> check_allocarray( sizeof(mpc_t* ), Ml )
##             for n in range(Ml):
##                 ef1[icusp][jcusp][n]=<mpc_t*> check_allocarray( sizeof(mpc_t ), Ql )
                
##     cdef double eps
##     eps = 2.0**float(1-H._dprec)
##     cdef int not_holom = int(not H._holomorphic)
##     cdef int is_weak = int(H._weak)

##     for n in range(Ml): 
##         for icusp in range(nc):
##             mpfr_set(nr.value,nvec[icusp][n],rnd_re)
##             for j in range(Ql):
##                 mpfr_mul(tmpr.value,Xm._entries[j],nr.value,rnd_re)
##                 mpfr_neg(tmpr.value,tmpr.value,rnd_re)
##                 mpc_set_fr(iargm.value,tmpr.value,rnd)
##                 mpfr_swap(mpc_realref(iargm.value),mpc_imagref(iargm.value))
##                 mpc_exp(iargm.value,iargm.value,rnd)
##                 #iargm=-inr*Xm[j]
##                 mpc_init2(ef2[icusp][n][j],prec)
##                 mpc_set(ef2[icusp][n][j],iargm.value,rnd)
##                 #ef2[n,j,icusp]=one*iargm #iargm.exp() #mpmath_ctx.exp(-iargm)
##         for jcusp in range(nc):
##             mpfr_set(nr.value,nvec[jcusp][n],rnd_re)
##             #nr=nvec[n,jcusp]
##             mpfr_mul(nrfourpi.value,nr.value,fourpi.value,rnd_re)
##             #nrfourpi=nr*fourpi
##             for icusp in range(nc):
##                 for j in range(Ql): 
##                     mpc_init2(ef1[icusp][jcusp][n][j],prec)
##                     mpc_set_si(ef1[icusp][jcusp][n][j],0,rnd)
##                     mpfr_set(ypb.value,Ypb[icusp][jcusp][j],rnd_re)
##                     if mpfr_zero_p(ypb.value)<>0:
##                         continue
##                     mpfr_set(xpb.value,Xpb[icusp][jcusp][j],rnd_re)
##                     mpfr_mul(tmpr_t,twopi.value,ypb.value,rnd_re)
##                     mpfr_neg(tmpr_t,tmpr_t,rnd_re)
##                     mpc_set_fr_fr(iargpb_t,tmpr_t,Xpb[icusp][jcusp][j],rnd)                    
##                     mpc_mul_fr(iargpb_t,iargpb_t,nr.value,rnd)
##                     mpc_exp(iargpb_t,iargpb_t,rnd)
##                     mpc_set(iargpb.value,iargpb_t,rnd)
##                     #iargpb=(nr*CF(-twopi*ypb,xpb)).exp()
##                     #iargpb=iargpb.exp()

##                     if mpfr_cmp_d(nr.value,eps)>0:
##                         mpc_set(ef1[icusp][jcusp][n][j],iargpb.value,rnd)
##                         #ef1[j,icusp,jcusp,n]=one*iargpb
##                     #elif mpfr_cmp_d(nr,-eps)<0 and (not H._holomorphic) and H._weak:
##                     elif mpfr_cmp_d(nr.value,-eps)<0 and not_holom==1 and is_weak==1:
##                         mpfr_abs(tmpr.value,nrfourpi.value,rnd_re)
##                         mpfr_mul(tmpr.value,tmpr.value,ypb.value,rnd_re)
##                         #tmpr=abs(nrfourpi)*ypb
##                         #tmp2=iargpb
##                         # pari: iargpb=iargpb*ckint.gamma_inc(tmpr)
##                         #iargpb=iargpb*ckint.gamma_inc(tmpr)
##                         #tmpr2 = RF(mpmath.mp.gammainc(kint,tmpr).real)
##                         #tmpr2 = RF(mpmath.mp.gammainc(kint,tmpr).real)
##                         if is_int==1 or is_half_int==1 and do_mpmath==0:
##                             #print "tmpr=",tmpr                            
##                             try:
##                                 if is_int==1:
##                                     if kinti>0:
##                                         incgamma_pint_c(tmpr2.value,kinti,tmpr.value)
##                                     else:
##                                         incgamma_nint_c(tmpr2.value,kinti,tmpr.value)
##                                         # print "incgamma_{0}:{1}={2}".format(is_int,kinti,tmpr2)
##                                 elif is_half_int==1:
##                                     incgamma_hint_c(tmpr2.value,kinti,tmpr.value)
##                             except ArithmeticError: ## In case we can not achieve the required error
##                                 tmpr2 = RF(mpmath.mp.gammainc(kint,tmpr).real)
##                         else:
##                             tmpr2 = RF(mpmath.mp.gammainc(kint,tmpr).real)
##                             #tmpr = RF(mpmath.mp.gammainc(kint,tmpr).real)
##                         #if tmpr2<>tmpr:
##                         #    print "is_int=",is_int
##                         #    print "is_half_int=",is_half_int
##                         #    print "mpmath.gammainc {0}={1}".format(kint,tmpr2)
##                         #    raise ArithmeticError
                        
                        
##                         mpc_mul_fr(iargpb.value,iargpb.value,tmpr2.value,rnd)
##                         #mpmath: iargpb=iargpb*RF(mpmath.mp.gammainc(kint,tmpr).real)
##                         #ef1[j,icusp,jcusp,n]=iargpb*ckint.gamma_inc(tmpr)
##                         mpc_set(ef1[icusp][jcusp][n][j],iargpb.value,rnd)
##                     elif mpfr_cmp_d(nr.value,-eps)<0:
##                         mpc_set_si(ef1[icusp][jcusp][n][j],0,rnd)
##                     else:
##                         ## Note that principal parts have to be determined
##                         ## and will appear in the right hand side
##                         #print "n+alpha=0 (jcusp,n)=",jcusp,n
##                         if variable_a0_plus[0][jcusp]==1:
##                             mpc_set_si(ef1[icusp][jcusp][n][j],1,rnd)
##                             #ef1[j,icusp,jcusp,n]=one
##                         elif variable_a0_minus[0][jcusp]==1:
##                             #mpfr_pow(tmpr.value,ypb.value,kint.value,rnd_re)
##                             mpc_set_fr(ef1[icusp][jcusp][n][j],Ypb[icusp][jcusp][j],rnd)
##                             if kint==0:
##                                 mpc_log(ef1[icusp][jcusp][n][j],ef1[icusp][jcusp][n][j],rnd)
##                             else:
##                                 mpc_pow_fr(ef1[icusp][jcusp][n][j],ef1[icusp][jcusp][n][j],kint.value,rnd)
##                             #ef1[j,icusp,jcusp,n]=ypb**kint
##                         else:
##                             ## If none of them are variables this means
##                             ## that they are both set in the principal parts
##                             mpc_set_si(ef1[icusp][jcusp][n][j],0,rnd)
##                             #ef1[j,icusp,jcusp,n]=one
##                         #if verbose>1:
##                         #print "ef(nr=0)[",j,icusp,jcusp,n,"]=",ef1[j,icusp,jcusp,n]
##                     if verbose > 3 and n-Ms==1:
##                         mpc_set(tmpc.value,ef1[icusp][jcusp][n][j],rnd)
##                         print "ef1[",n,j,icusp,"]=",CC(tmpc.real(),tmpc.imag())
##     for n in range(s):
##         for l in range(s):        
##             mpc_set_ui(V._matrix[n][l],0,rnd)
##     for l in range(Ml): 
##         for jcusp in range(nc):
##             lj=Ml*jcusp+l
##             for icusp in range(nc):
##                 for j in range(Ql): 

##                     if mpfr_zero_p(mpc_realref(Cvec[icusp][jcusp][j]))<>0 and mpfr_zero_p(mpc_imagref(Cvec[icusp][jcusp][j]))<>0:
##                         continue
##                     #if mpfr_zero_p(Ypb[icusp][jcusp][j])<>0:
##                     #    continue                        
##                     #mpc_set(tmp1.value,rnd)
##                     mpc_mul(tmp1.value,ef1[icusp][jcusp][l][j],Cvec[icusp][jcusp][j],rnd)
##                     #tmp1=ef1[j,icusp,jcusp,l]*ch
##                     for n in range(Ml): 
##                         ni=Ml*icusp+n
##                         mpc_mul(tmp2.value,tmp1.value,ef2[icusp][n][j],rnd)
##                         #tmp2=tmp1*ef2[n,j,icusp]
##                         mpc_add(V._matrix[ni][lj],V._matrix[ni][lj],tmp2.value,rnd)
##                         #V[ni,lj]=V[ni,lj]+tmp1*
##                         if verbose > 2 and ni==1 and lj==1:
##                             print "-------------------"
##                             print "V[1,1](",j,")=",CC(V[ni,lj].real(),V[ni,lj].imag())
##                             mpc_set(tmpc.value,ef1[icusp][jcusp][n][j],rnd)
##                             print "ef1(",j,")=",CC(tmpc.real(),tmpc.imag())
##                             mpc_set(ch.value,Cvec[icusp][jcusp][j],rnd)
##                             print "cv(",j,")=",CC(ch.real(),ch.imag())
##                             mpc_set(tmpc.value,ef2[icusp][n][j],rnd)
##                             print "ef2(",j,")=",tmpc
                        
##     if verbose>1:
##         print "Mi,Ms=",Ms,Mf
##         print "V1(",0,0,")=",V[0,0]
##         print "V1(",1,0,")=",V[1,0]
##         print "V1(",0,1,")=",V[0,1]
##         print "V1(",1,1,")=",V[1,1]
##     cdef int nrows,ncols
##     nrows = int(V.nrows()); ncols = int(V.ncols())
##     for n in range(nrows):
##         for l in range(ncols):        
##             mpc_div_ui(V._matrix[n][l],V._matrix[n][l],Qfaki,rnd)
##             #V[n,l]=V[n,l]/Qfak
##     if verbose>1: 
##         print "V1(",0,0,")=",V[0,0]
##         print "V1(",1,0,")=",V[1,0]            
##         print "V1(",0,1,")=",V[0,1]
##         print "V1(",1,1,")=",V[1,1]
##     cdef MPComplexNumber f1,f2,ppc,summa,ppc_minus,summa_minus
##     cdef RealNumber lr,arg,nrY2pi,kbes
##     lr = RF(0); arg=RF(0);nrY2pi=RF(0); kbes=RF(0)
##     f1 = CF(0); f2=CF(0); ppc=CF(0); summa=CF(0)
##     ppc_minus = CF(0); summa_minus=CF(0)
##     for n in range(Ml):
##         for icusp in range(nc):
##             mpfr_set(nr.value,nvec[icusp][n],rnd_re)
##             mpfr_mul(nrY2pi.value,nr.value,twopiY.value,rnd_re)
##             #nrY2pi=nr*twopiY
##             ni=Ml*icusp+n
##             if mpfr_cmp_d(nr.value,eps)>0:
##                 mpfr_neg(kbes.value,nrY2pi.value,rnd_re)
##                 mpfr_exp(kbes.value,kbes.value,rnd_re)
##                 #kbes=(-nrY2pi).exp()
##             elif mpfr_cmp_d(nr.value,-eps)<0 and not_holom==1 and is_weak==1:
##                 #kbes=RF(mpmath.mp.gammainc(kint,abs(nr)*fourpiY))
##                 #kbes=ckint.gamma_inc(abs(nr)*fourpiY).real()                
##                 mpfr_mul(tmpr.value,nr.value,fourpiY.value,rnd_re)
##                 mpfr_abs(tmpr.value,tmpr.value,rnd_re)
##                 if is_int==1 or is_half_int==1 and  do_mpmath==0:
##                     try:
##                         if is_int==1:
##                             if kinti>0:
##                                 incgamma_pint_c(kbes.value,kinti,tmpr.value)
##                             else:
##                                 incgamma_nint_c(kbes.value,kinti,tmpr.value)
##                         elif is_half_int==1:
##                             incgamma_hint_c(kbes.value,kinti,tmpr.value)                                
##                     except ArithmeticError: ## In case we can not achieve the required error
##                         kbes = RF(mpmath.mp.gammainc(kint,tmpr).real)                    
##                 else:
##                     kbes = RF(mpmath.mp.gammainc(kint,tmpr).real)

##                 mpfr_abs(tmpr2.value,nrY2pi.value,rnd_re)
##                 mpfr_exp(tmpr2.value,tmpr2.value,rnd_re)
##                 if verbose>1 and ni==0:
##                     print "Arg=",tmpr
##                     print "Igamma(",kinti,",Arg)=",kbes
##                 mpfr_mul(kbes.value,kbes.value,tmpr2.value,rnd_re)
##                 if verbose>1 and ni==0:
##                     print "expfax=",tmpr2
##                     print "Igamma(",kinti,",Arg)*exp()=",kbes
##                 #kbes=kbes*(-nrY2pi).exp()
##             elif mpfr_cmp_d(nr.value,-eps)<0:
##                 mpfr_set_ui(kbes.value,0,rnd_re) #=zero
##             else:
##                 if variable_a0_plus[0][icusp]==1:
##                     mpfr_set_ui(kbes.value,1,rnd_re) #=one
##                 elif variable_a0_minus[0][icusp]==1:
##                     if kint<>0:
##                         mpfr_pow(kbes.value,Y.value,kint.value,rnd_re)
##                     else:
##                         mpfr_log(kbes.value,Y.value,rnd_re)
##                     #kbes=Y**kint
##                 else:
##                     #kbes=one
##                     #raise ValueError,"Need to specify at least one of the n=0 terms!"
##                     mpfr_set_ui(kbes.value,0,rnd_re)
##                     #kbes=zero
##             mpc_sub_fr(V._matrix[ni][ni],V._matrix[ni][ni],kbes.value,rnd)
##             if verbose>1:
##                 if ni==0:
##                     print "Arg=",tmpr
##                     print "kbes(",0,")=",kbes
##                     print "V[1,1]=",V[0,0]
##                 if ni==1:
##                     print "nr=",nr
##                     print "kint,kinti=",kint,kinti
##                     print "Arg=",tmpr
##                     print "kbes(",1,")=",kbes
##                     print "V[1,1]=",V[1,1]
##             # setting the right hand side of the system
##             for k in range(d):
##                 mpc_set_ui(RHS._matrix[ni][k],0,rnd) #=CF(0)
##                 for i in range(num_ppplus):
##                 #for (jcusp,l) in PPplus[k].keys():
##                     mpc_set(ppc.value,PPplus_values[k][i],rnd)
##                     jcusp = PPplus_cusp[k][i]
##                     l = PPplus_n[k][i]
##                     #ppc=CF(PPplus[k][(jcusp,l)])
##                     #mpfr_set(lr.value,nvec[jcusp][l],rnd_re)
##                     #lr=RF(l+H.alpha(jcusp)[0])
##                     #if(abs(lr) <= mpmath.eps()):
##                     #    continue
##                     ### Note: The constant term is treated below
##                     if mpc_is_zero(ppc.value)==1: #==zero or PPplus[k][(jcusp,l)]==0:
##                         continue
##                     mpc_set_ui(summa.value,0,rnd) #=CF(0)
##                     if verbose>2:
##                         mpfr_set(lr.value,PPplus_lal[k][i],rnd_re)
##                         print "l=",lr
##                         print "icusp,jcusp=",icusp,jcusp
##                         print "alpha(",jcusp,")=",H.alpha(jcusp)[0]
##                         print "n=",n,nr," n=",n+Ms
##                         print "ppc=",ppc
##                     for j in range(Ql):   
##                         if mpc_is_zero(Cvec[icusp][jcusp][j])==1:
##                             continue
##                         if mpfr_zero_p(Ypb[icusp][jcusp][j])<>0:
##                             continue
##                         mpfr_mul(tmpr_t,Ypb[icusp][jcusp][j],twopi.value,rnd_re)
##                         if verbose>3:
##                             mpfr_set(tmpr.value,Ypb[icusp][jcusp][j],rnd_re)
##                             print "ypb=",tmpr
##                             print "lr=",lr
##                             print "twopi=",twopi
##                         mpfr_neg(tmpr_t,tmpr_t,rnd_re)
##                         mpfr_mul(tmpr_t,tmpr_t,PPplus_lal[k][i],rnd_re)
##                         mpfr_exp(tmpr_t,tmpr_t,rnd_re)
##                         # tmpr_t = exp(-2pi*l*ypb)
##                         if verbose>2:
##                             mpfr_set(tmpr.value,tmpr_t,rnd_re)
##                             print "-2pi*l*ypb=",tmpr
##                         mpc_set_fr(f1.value,tmpr_t,rnd)
##                         mpfr_set(xpb.value,Xpb[icusp][jcusp][j],rnd_re)
##                         #arg = lr*xpb-nr*Xm[j]
##                         mpfr_mul(tmpr_t,PPplus_lal[k][i],Xpb[icusp][jcusp][j],rnd_re)
##                         mpfr_mul(arg.value,nr.value,Xm._entries[j],rnd_re)
##                         mpfr_sub(arg.value,tmpr_t,arg.value,rnd_re)
##                         mpc_set_fr_fr(iargpb_t,zero.value,arg.value,rnd)
##                         mpc_exp(f2.value,iargpb_t,rnd)
##                         #f2 = CF(0,arg).exp()
##                         if verbose>3:
##                             print "f1(",j,")=",f1
##                             print "arg=",lr,"*",xpb,"-",nr,"*",Xm[j],"=",arg
##                             print "f2(",j,")=",f2
##                         mpc_mul(tmpc.value,f1.value,f2.value,rnd)
##                         mpc_mul(tmpc.value,tmpc.value,ppc.value,rnd)
##                         mpc_mul(tmpc.value,tmpc.value,Cvec[icusp][jcusp][j],rnd)
##                         mpc_add(summa.value,summa.value,tmpc.value,rnd)
##                         if verbose>3:
##                             print "tmpc=",tmpc
##                             print "summa=",summa                 
## #                        summa=summa+ch*tmpc
##                     if verbose>2:
##                         print "summa=",RR(summa.real()),RR(summa.imag())
##                     #RHS[ni,k]=RHS[ni,k]+summa/Qfak
##                     mpc_div_ui(summa.value,summa.value,Qfaki,rnd)
##                     mpc_add(RHS._matrix[ni][k],RHS._matrix[ni][k],summa.value,rnd)
##                 if verbose>2:
##                     print "RHS0[",ni,k,"]=",RHS[ni][k]
##                     print "icusp,n+Ms=",icusp,n+Ms
##                     print "PPplus.keys=",PPplus[k].keys()
##                 has_key = 0
##                 for j in range(num_ppplus):
##                     if PPplus_cusp[k][j]==icusp and PPplus_n[k][j]==n+Ms:
##                         has_key = 1
##                         mpc_set(ppc.value,PPplus_values[k][j],rnd)
##                         break
##                 if has_key==1:
##                     #if PPplus[k].has_key((icusp,n+Ms)):
##                     #if( abs(nr) > mpmath.eps()):
##                     #ppc=CF(PPplus[k][icusp,n+Ms])
##                     mpfr_abs(tmpr_t,nrY2pi.value,rnd_re)
##                     mpfr_exp(tmpr_t,tmpr_t,rnd_re)
##                     mpc_mul_fr(tmpc_t,ppc.value,tmpr_t,rnd)
##                     mpc_sub(RHS._matrix[ni][k],RHS._matrix[ni][k],tmpc_t,rnd)
##                     if verbose>2:
##                         print "n=",n
##                         print "icusp=",icusp
##                         print "ppc=",ppc
##                         print "nrY2pi=",nrY2pi
##                     #RHS[ni,k]=RHS[ni,k]-ppc*(-nrY2pi).exp()
##                 if verbose>2:
##                     print "RHS1[",ni,k,"]=",RHS[ni][k]
##                 for i in range(num_ppminus):
##                     #for (jcusp,l) in PPminus[k].keys():
##                     jcusp = PPminus_cusp[k][i]
##                     l = PPminus_n[k][i]
##                     mpc_set(ppc_minus.value,PPminus_values[k][i],rnd)
##                     #ppc_minus = CF(PPminus[k][(jcusp,l)])
##                     #if ppc_minus==zero or PPminus[k][(jcusp,l)]==0:
##                     if mpc_is_zero(ppc_minus.value)==1:
##                         continue
##                     if verbose>2:
##                         print "ppart_minus=",ppc_minus
##                     mpc_set_ui(summa_minus.value,0,rnd) #=zero #summa_minus
##                     for j in range(Ql):
##                         #if mpfr_zero_p(mpc_realref(Cvec[icusp][jcusp][j]))<>0 and mpfr_zero_p(mpc_imagref(Cvec[icusp][jcusp][j]))<>0:
##                         if mpc_is_zero(Cvec[icusp][jcusp][j])==1:
##                             continue
##                         if mpfr_zero_p(Ypb[icusp][jcusp][j])<>0:
##                             continue
##                         mpc_set(ch.value,Cvec[icusp][jcusp][j],rnd)
##                         mpfr_set(ypb.value,Ypb[icusp][jcusp][j],rnd_re)
##                         if kint==0:
##                             mpfr_log(tmpr.value,ypb.value,rnd_re)
##                         else:
##                             mpfr_pow(tmpr.value,ypb.value,kint.value,rnd_re)
##                         mpc_set_fr(tmpc_t,tmpr.value,rnd) #ypb**kint

##                         mpfr_mul(tmpr.value,nr.value,Xm._entries[j],rnd_re)
##                         mpfr_neg(tmpr.value,tmpr.value,rnd_re)
##                         mpc_set_fr(iargm.value,tmpr.value,rnd)
##                         mpfr_swap(mpc_realref(iargm.value),mpc_imagref(iargm.value))
##                         mpc_exp(iargm.value,iargm.value,rnd)
##                         mpc_mul(tmpc.value,tmpc_t,iargm.value,rnd)
##                         #tmpc=tmpc_minus*CF(0,-nr*Xm[j]).exp()
##                         mpc_mul(tmpc.value,tmpc.value,ch.value,rnd)
##                         mpc_add(summa_minus.value,summa_minus.value,tmpc.value,rnd)
##                         #summa_minus=summa_minus+ch*tmpc
##                         #summa=summa+ef2[n(2):,j,icusp]*tmpc*ch*ppc
##                         if verbose>2:
##                             print "j=",j
##                             print "ypb=",ypb
##                             print "ch=",ch
##                             #print "tmpc_minus=",tmpc_minus
##                             print "tmpc=",tmpc
##                             print "summa-(",j,")=",summa_minus
##                     if verbose>2:
##                         print "Summa-(",ni,")=",summa_minus
##                     mpc_div_ui(summa_minus.value,summa_minus.value,Qfaki,rnd)
##                     mpc_mul(summa_minus.value,summa_minus.value,ppc_minus.value,rnd)
##                     mpc_add(RHS._matrix[ni][k],RHS._matrix[ni][k],summa_minus.value,rnd)
##                     #RHS[ni,k]=RHS[ni,k]+ppc_minus*summa_minus/Qfak
##                 #print "summa(",ni,")=",RHS[ni,k]
##                 if verbose>2:
##                     print "RHS2[{0},{1}]={2}".format(ni,k,RHS[ni,k])
##                 has_key = 0
##                 for j in range(num_ppminus):
##                     if PPminus_cusp[k][j]==icusp and PPminus_n[k][j]==n+Ms:
##                         has_key = 1
##                         mpc_set(ppc_minus.value,PPminus_values[k][j],rnd)
##                         break
##                 if has_key==1:
##                     #if PPminus[k].has_key((icusp,n+Ms)):
##                     #ppc_minus = CF(PPminus[k][(icusp,n+Ms)])
##                     if kinti==0:
##                         mpfr_log(tmpr.value,Y.value,rnd_re)
##                     else:
##                         mpfr_pow(tmpr.value,Y.value,kint.value,rnd_re)
##                     mpc_mul_fr(ppc_minus.value,ppc_minus.value,tmpr.value,rnd)
##                     mpc_sub(RHS._matrix[ni][k],RHS._matrix[ni][k],ppc_minus.value,rnd)
##                     #RHS[ni,k]=RHS[ni,k]-ppc_minus*Y**kint
##                     #print "subtracting:",ppc_minus*mpmath.power(Y,kint)
##                 if verbose>2:
##                     print "RHS3[{0},{1}]={2}".format(ni,k,RHS[ni,k])
##                     #print "RHS(",ni,")-=",ppc*mpmath_ctx.exp(-nrY2pi)
##     if verbose>1: 
##         for n in range(2):
##             print "V2(",n,n,")=",V[n,n]
##     # Clearing up allocated variables

##     for j in range(d):
##         sig_free(variable_a0_minus[j])
##         sig_free(variable_a0_plus[j])        
##         sig_free(PPplus_cusp[j])
##         sig_free(PPplus_n[j])
##         sig_free(PPminus_cusp[j])
##         sig_free(PPminus_n[j])
##         for l in range(num_ppplus):
##             mpc_clear(PPplus_values[j][l])
##             mpfr_clear(PPplus_lal[j][l])
##         sig_free(PPplus_lal[j])
##         sig_free(PPplus_values[j])
##         for l in range(num_ppminus):
##             mpc_clear(PPminus_values[j][l])
##         sig_free(PPminus_values[j])
        
##     sig_free(PPplus_lal)
##     sig_free(PPplus_cusp)
##     sig_free(PPplus_n)
##     sig_free(PPplus_values)
##     sig_free(PPminus_cusp)
##     sig_free(PPminus_n)
##     sig_free(PPminus_values)

##     sig_free(variable_a0_minus)
##     sig_free(variable_a0_plus)




##     if Ypb<>NULL:
##         for i in range(nc):
##             if Ypb[i]<>NULL:
##                 for j in range(nc):
##                     if Ypb[i][j]<>NULL:
##                         sig_free(Ypb[i][j])
##                 sig_free(Ypb[i])
##         sig_free(Ypb)
##     if Xpb<>NULL:
##         for i in range(nc):
##             if Xpb[i]<>NULL:
##                 for j in range(nc):
##                     if Xpb[i][j]<>NULL:
##                         sig_free(Xpb[i][j])
##                 sig_free(Xpb[i])
##         sig_free(Xpb)
##     if Cvec<>NULL:
##         for i in range(nc):
##             if Cvec[i]<>NULL:
##                 for j in range(nc):
##                     if Cvec[i][j]<>NULL:
##                         sig_free(Cvec[i][j])
##                 sig_free(Cvec[i])
##         sig_free(Cvec)
##     if nvec<>NULL:
##         for i in range(nc):
##             if nvec[i]<>NULL:
##                 for l in range(Ml):
##                     mpfr_clear(nvec[i][l])
##                 sig_free(nvec[i])
##         sig_free(nvec)
##     if ef2<>NULL:
##         for icusp in range(nc):
##             if ef2[icusp]<>NULL:
##                 for n in range(Ml):
##                     if ef2[icusp][n]<>NULL:
##                         for j in range(Ql):
##                             mpc_clear(ef2[icusp][n][j])
##                         sig_free(ef2[icusp][n])
##                 sig_free(ef2[icusp])
##         sig_free(ef2)
##     if ef1<>NULL:
##         for icusp in range(nc):
##             if ef1[icusp]<>NULL:
##                 for jcusp in range(nc):
##                     if ef1[icusp][jcusp]<>NULL:
##                         for n in range(Ml):
##                             if ef1[icusp][jcusp][n]<>NULL:
##                                 for j in range(Ql):
##                                     mpc_clear(ef1[icusp][jcusp][n][j])
##                                 sig_free(ef1[icusp][jcusp][n])
##                         sig_free(ef1[icusp][jcusp])                        
##                 sig_free(ef1[icusp])
##         sig_free(ef1)
        
##         #    print "V2(",44,n,")=",V[V.rows-1,n]
##     mpc_clear(iargpb_t)
##     mpc_clear(tmpc_t)
##     mpfr_clear(tmpr_t)
##     W=dict()
##     W['var_a+']=pp_info['variable_a0_plus']
##     W['var_a-']=pp_info['variable_a0_minus']
##     #if H._verbose>0:
##     #    print "alphas=",H.alphas()
##     W['alphas']=H.alphas()
##     W['V']=V
##     #W['cv']=Cv
##     W['RHS']=RHS
##     W['Ms']=Ms
##     W['Mf']=Mf
##     W['Ml']=Ml
##     W['nc']=nc
##     W['PP']=principal_parts
##     W['space']=H
##     W['rdim']=H._rdim
##     return W

from sage.all import RR
cpdef setup_matrix_for_sl2z(H,RealNumber Y,int M, int Q,int verbose=2):
    r"""
    Optimized algorithm for SL(2,Z).
    Also: for now, for holomorphic forms only.
    """
    if H.group().index()<>1:
        raise ValueError,"Use only for SL(2,Z)"
    RF = Y.parent()
    cdef int prec =RF.prec()
    CF = MPComplexField(prec)
    cdef RealNumber pi,one,two,twopi,p,twopiY,fourpiY,kint,Qfak,minus_twopiy,half
    cdef RealNumber weight
    cdef MPComplexNumber twopii
    Qfak=RF(2*Q)
    weight = RF(H._weight)
    pi=RF.pi() #mpmath_ctx.pi()
    one=RF(1) #mpmath_ctx.mpf(1)
    two=RF(2) #mpmath_ctx.mpf(2)
    zero=RF(0) #mpmath_ctx.mpf(0)
    half = RF(0.5)
    twopi=two*pi
    fourpi=two*twopi
    twopiY=twopi*Y
    minus_twopiy=-twopiY
    fourpiY=two*twopiY
    p=(weight-one)/two
    kint=one-weight

    cdef Matrix_complex_dense V
    cdef Vector_complex_dense Cv
    #cdef vector Xm,Xpb,Ypb
    MS = MatrixSpace(CF,M+1,M+1)
    V = Matrix_complex_dense(MS,0,True,True)
    MS = MatrixSpace(CF,2*Q,M+1)
    fak1 = Matrix_complex_dense(MS,0,True,True)
    fak2 = Matrix_complex_dense(MS,0,True,True)
    cdef int Ql,n,l,j,Qs,Qf
    cdef int a,b,c,d
    cdef double xdb,ydb
    ydb = <double>Y
    cdef MPComplexNumber tmpcc,jfak
    cdef RealNumber x1,y1,kap
    cdef RealNumber alpha,ll,nn,ar1
    alpha=RF(H.alpha(0)[0])
    jfak = CF(1); tmpcc=CF(1)
    x1=RF(0); y1=RF(0)
    cdef RealNumber fourQ,y11,sqy
    fourQ=RF(4*Q)
    cdef mpfr_t nnn,kappa,mphalf,xm,tmpr,weight_half,sqrtY,twopiypb,sqypb
    cdef mpc_t arg,tmpc
    mpfr_init2(nnn,prec)
    mpfr_init2(tmpr,prec)
    mpfr_init2(sqrtY,prec)
    mpfr_init2(kappa,prec)
    mpfr_init2(mphalf,prec)
    mpfr_init2(weight_half,prec)
    mpfr_init2(twopiypb,prec)
    mpfr_init2(sqypb,prec)
    mpfr_init2(xm,prec)
    mpc_init2(arg,prec)
    mpc_init2(tmpc,prec)
    mpfr_div_ui(weight_half,weight.value,2,rnd_re)
    mpfr_sqrt(sqrtY,Y.value,rnd_re)
    #Xm=vector(RF,2*Q)
    #Xpb=vector(RF,2*Q)
    #Ypb=vector(RF,2*Q)
    Cv=Vector_complex_dense(vector(CF,2*Q).parent(),0)
    mpfr_set_ui(mphalf,1,rnd_re)
    mpfr_div_ui(mphalf,mphalf,2,rnd_re)
    cdef int Q4 = 4*Q
    cdef int jj

    for j from 0 <= j < 2*Q:
        jj = 2*j+1
        mpfr_set_si(xm,jj,rnd_re)
        mpfr_div_ui(xm,xm,Q4,rnd_re)
        mpfr_sub(xm,xm,mphalf,rnd_re)
        xdb = mpfr_get_d(xm,rnd_re) 
        a,b,c,d=pullback_to_psl2z_mat(xdb,ydb)
        mpfr_set(y1.value,Y.value,rnd_re)
        mpfr_set(x1.value,xm,rnd_re)
        _apply_sl2z_map_mpfr(x1.value,y1.value,a,b,c,d)
        jfak =  j_fak_int_mpfr(-c,a,x1,y1,weight,H._unitary_action)
        if H.multiplier().is_trivial():
            mpc_set(Cv._entries[j],jfak.value,rnd)
        else:
            v = H.multiplier()(SL2Z([d,-b,-c,a]))
            if hasattr(v,"complex_embedding"):
                Cv[j]=v.complex_embedding(prec)*CF(jfak)
            else:
                Cv[j]=CF(v)*CF(jfak)
        #print "Cv=",Cv[j]
        mpfr_mul(xm,xm,twopi.value,rnd_re)
        if verbose>2:
            print "Xm[",j,"]=",mpfr_get_d(xm,rnd_re)
            
            print "xpb[",j,"]=",x1
            print "ypb[",j,"]=",y1
            print "Cv[",j,"]=",Cv[j]
        x1=x1*twopi
        mpfr_mul(twopiypb,twopi.value,y1.value,rnd_re)
        mpfr_neg(twopiypb,twopiypb,rnd_re)
        mpfr_sqrt(sqypb,y1.value,rnd_re)
        for n from 0 <= n<= M:
            mpfr_add_si(nnn,alpha.value,n,rnd_re) 
            mpfr_mul(tmpr,nnn,xm,rnd_re)
            mpfr_neg(tmpr,tmpr,rnd_re)
            mpc_set_fr(arg,tmpr,rnd)
            mpfr_swap(mpc_realref(arg),mpc_imagref(arg))
            mpc_exp(arg,arg,rnd)
            mpc_set(fak2._matrix[j][n],arg,rnd)
            if n==-1 and verbose>3:
                print "Xm[",j,"]=",mpfr_get_d(xm,rnd_re)
                mpc_set(tmpcc.value,arg,rnd)
                print "iargm[",j,"]=",tmpcc
                print "fak2[",n,j,0,"]=",fak2[j][n]
            if n == 0 and alpha==0:
                if H._unitary_action == 0:
                    mpc_set_ui(fak1._matrix[j][n],1,rnd)
                else:
                    mpfr_pow(tmpr,y1.value,weight_half,rnd_re)
                    mpc_set_fr(fak1._matrix[j][n],tmpr,rnd)
                    #fak1[j,n]= (CF(-y11*nn,Xpb[j]*nn)).exp()
            else:
                mpc_set_fr_fr(arg,twopiypb,x1.value,rnd)
                mpc_mul_fr(arg,arg,nnn,rnd)
                mpc_exp(arg,arg,rnd)  # arg = ex(-nn(ypb,xpb))
                if H._unitary_action == 1:
                    mpfr_mul(tmpr,y1.value,nnn,rnd_re)
                    mpfr_pow(tmpr,tmpr,p.value,rnd_re)
                    mpfr_mul(tmpr,tmpr,sqypb,rnd_re)
                    mpc_mul_fr(arg,arg,tmpr,rnd)
                mpc_set(fak1._matrix[j][n],arg,rnd)
                if n==1 and verbose>3:
                    print "Ypb[",j,"]=",y1 
                    mpc_set(tmpcc.value,arg,rnd)
                    print "iargpb[",j,"]=",tmpcc
                    print "ef1[",n,j,0,"]=",fak1[j][n]
                
    for n from 0 <= n <= M:
        for l from 0 <= l <= M:
            mpc_set_ui(V._matrix[n][l],0,rnd)
    for l from 0 <= l <= M:
        for j from 0 <= j < 2*Q:
            for n from 0 <= n <= M:
                mpc_mul(tmpc,fak1._matrix[j][l],fak2._matrix[j][n],rnd)
                mpc_mul(tmpc,tmpc,Cv._entries[j],rnd)
                #mpc_mul_fr(tmpc,tmpc,kappa,rnd)
                mpc_add(V._matrix[n][l],V._matrix[n][l],tmpc,rnd)
                if H._verbose > 2 and n==1 and l==1:
                    print "-------------------"
                    #print "ni,lj=",ni,lj
                    #print "icusp,jcusp=",icusp,jcusp
                    print "V[1,1](",j,")=",CC(V[n,l].real(),V[n,l].imag())
                    #print "xpb=",RR(Xpb[icusp,jcusp,j])
                    print "ef1(",j,")=",CC(fak1[j][l].real(),fak1[j][l].imag())
                    print "cv(",j,")=",CC(Cv[j].real(),Cv[j].imag())
                    print "ef2(",j,")=",CC(fak2[j][n].real(),fak2[j][n].imag())
         

                
                #mpfr_set(x1.value,kappa,rnd_re)
                #V[n,l]=V[n,l]+Cv[j]*fak1[j,l]*fak2[j,n]
    if H._verbose>1: 
        for n from 0<=n<2:
            print "V0(",n,n,")=",V[n,n]
    for n from 0 <= n <= M:
        for l from 0 <= l <= M:
            #mpc_div_fr(V._matrix[n][l],V._matrix[n][l],Qfak.value,rnd)
            V[n,l]=V[n,l]/Qfak
    if H._verbose>1: 
        for n from 0<=n<2:
            print "V1(",n,n,")=",V[n,n]
    for n from 0 <= n <= M:
        #nn = RF(n)+alpha
        mpfr_add_si(nnn,alpha.value,n,rnd_re) 
        if n==0 and alpha==0:
            if H._unitary_action==1:
                #kap = Y**(weight*half)
                mpfr_pow(kappa,Y.value,weight_half,rnd_re)
            else:
                mpfr_set_ui(kappa,1,rnd_re)
        else:            
            mpfr_mul(tmpr,minus_twopiy.value,nnn,rnd_re)
            mpfr_exp(kappa,tmpr,rnd_re)
            if H._unitary_action==1:
                mpfr_mul(tmpr,Y.value,nnn,rnd_re)
                mpfr_pow(tmpr,tmpr,p.value,rnd_re)
                mpfr_mul(tmpr,tmpr,sqrtY,rnd_re)
                #kap = Y.sqrt()*(Y*nn)**p*(-twopi*nn*Y).exp() 
                mpfr_mul(kappa,kappa,tmpr,rnd_re)
            #V[n,n]=V[n,n]-kap
        mpc_sub_fr(V._matrix[n][n],V._matrix[n][n],kappa,rnd)
    if H._verbose>1: 
        for n from 0<=n<2:
            print "V2(",n,n,")=",V[n,n]
        
    mpfr_clear(nnn)
    mpfr_clear(tmpr)
    mpfr_clear(kappa)
    mpfr_clear(mphalf)
    mpfr_clear(sqrtY)
    mpfr_clear(weight_half)
    mpfr_clear(xm)
    mpfr_clear(twopiypb)
    mpfr_clear(sqypb)
    mpc_clear(arg)
    mpc_clear(tmpc)
    W = dict()
    W['V']=V
    #W['xm']=Xm
    #W['xpb']=Xpb
    #W['ypb']=Ypb
    #W['cv']=Cv
    #W['f1']=fak1
    #W['f2']=fak2
    W['Ms']=0
    W['Mf']=M
    W['Ml']=M+1
    W['nc']=1
    W['PP']=[]
    W['space']=H
    W['alphas']=[H.alpha(0)]
    W['var_a+']=[]
    W['var_a-']=[]
 
    return W


def solve_for_sl2z(W):
    V=W['V']
    if W['space'].alpha(0)[0]==0:
        V.delete_row(0)
        V.delete_row(0)
        b = V.column(1)
        V.delete_column(0)
        V.delete_column(0)    
    else:
        V.delete_row(1)
        b = V.column(1)
        V.delete_column(1)
    c = V.solve(-b)
    return c

### The oldest backup...

## def setup_matrix_for_harmonic_Maass_waveforms_sv_bak_22(H,Y_in,int M,int Q,principal_parts):
##     r"""

##     Set up the matrix for the system of equations giving the Fourier coefficients of a Harmonic Maass waveforms.

##     INPUT:
    
##         - ``H`` -- Space of harmonic weak Maass forms
##         - ``Y`` -- height of horocycle used for sampling (mpmath real)
##         - ``k`` -- weight (mpmath real) 
##         - ``M`` -- integer
##         - ``Q`` -- integer > M
##         - ``PP``-- dict : Principal parts at cusp nr. j = \Sum_ PP(j,n) q^[-n]

##     OUTPUT:
    
##         - ``W`` -- dictionary
##         - ``W['Mf']`` -- M start
##         - ``W['nc']`` -- number of cusps
##         - ``W['Ms']`` -- M stop
##         - ``W['V']``  -- ((Ms-Mf+1)*num_cusps)**2 matrix 
##         - ``W['RHS']``  -- ((Ms-Mf+1)*num_cusps) matrix
    
    
##     EXAMPLES::

##         sage: setup_matrix_for_harmonic_Maass_waveforms_sv(MySubgroup(Gamma0(1)),mpmath.mpf(0.5),mpmath.mpf(0),10,20,{(0,-1):1})

    

##     """
##     cdef int l,j,icusp,jcusp,n,ni,li,iMl,iMs,iMf,iQs,iQf,iQl,s,nc
##     cdef double RR,YY
##     mpmath_ctx=H._mp_ctx
##     mpmath_ctx=mpmath.mp
##     #if( not test_ctx(Y,mpmath_ctx)):
##     #    raise TypeError," Need Mpmath reals as input!"
##     #GG=H._G
##     if(not isinstance(Y_in,type(mpmath_ctx.mpf(1.0)))):        
##         Y = mpmath_ctx.mpf(Y_in)
##     else:
##         Y = Y_in
##     weight=H._weight
##     pi=mpmath_ctx.pi()
##     one=mpmath_ctx.mpf(1)
##     two=mpmath_ctx.mpf(2)
##     zero=mpmath_ctx.mpf(0)
##     twopi=two*pi
##     fourpi=two*twopi
##     twopiY=twopi*Y
##     fourpiY=two*twopiY
##     p=(weight-mpmath_ctx.mpf(1))/mpmath_ctx.mpf(2)
##     kint=mpmath_ctx.mpf(1)-weight
##     #print "kint=",kint
##     nc=int(H._group.ncusps())
##     if(Q<M):
##         Q=M+20
##     if(H._holomorphic):
##         Ms=0; Mf=M; Ml=Mf-Ms+1
##     else:
##         Ms=-M; Mf=M; Ml=Mf-Ms+1
##     Qs=1-Q; Qf=Q
##     Qfak=mpmath_ctx.mpf(2*Q)
##     iMl=int(Ml); iMs=int(Ms); iMf=int(Mf)
##     iQf=int(Qf); iQs=int(Qs); iQl=int(Qf-Qs+1)
##     # Pullpack points
##     pb=pullback_pts_mpmath(H,Qs,Qf,Y)
##     Xm=pb['xm']; Xpb=pb['xpb']; Ypb=pb['ypb']; Cv=pb['cvec']
##     #[Xm,Xpb,Ypb,Cv]=pullback_pts_hw(H,Qs,Qf,Y)
##     s=int(nc*iMl)
##     #print "Q=",Q
##     #print "kint=",kint
##     V=mpmath_ctx.matrix(int(s),int(s))

##     PPplus = list(); PPminus = list()
##     for pp in principal_parts:
##         PPplus.append(pp['+'])
##         PPminus.append(pp['-'])
##     # cehcking correctness of porincipal parts
##     for pp in PPplus:
##         for icusp,n in pp:
##             if((icusp not in range(nc)) or n>0):
##                 raise ValueError,"Principal parts are in wrong format!"
##     for pp in PPminus:
##         for icusp,n in pp:
##             if((icusp not in range(nc)) or n>0):
##                 raise ValueError,"Principal parts are in wrong format!"            

##     # if we have a holomorphic function we  don't want any non-holomorphic principal parts
##     if(H._holomorphic):
##         for i in range(len(PPminus)):
##             for j in range(0,nc):
##                 PPminus[i][(j,0)]=0


##     #if(principal_parts.has_key('+')):
##     #    PP = principal_parts['+']
##     #else:
##     #    PP = [{}]
##     #if(principal_parts.has_key('-')):
##     #    PPminus = principal_parts['-']
##     #else:
##     #    PPminus = [{}] 
##     if(not isinstance(PPplus,list)):
##         raise ValueError,"Need a list of (+) principal parts!"
##     if(not isinstance(PPminus,list)):
##         raise ValueError,"Need a list of (-) principal parts!"
##     if(len(PPminus)<>len(PPplus)):
##         raise ValueError," Need equal number of holomorphic and non-holomorphic principal parts! Got:%s" % principal_parts
##     d = len(PPplus)
##     if(H._verbose>0): 
##         print "PPplus=",PPplus
##         print "PPminus=",PPminus
##     variable_a0_plus=dict()
##     variable_a0_minus=dict()
##     for j in range(nc):
##         if(H.alpha(j)[1]<>1):
##             #variable_a0_plus[j]=False
##             #variable_a0_minus[j]=False
##             continue
##         variable_a0_plus[j]=True
##         variable_a0_minus[j]=True
##         if(PPplus[0].has_key((j,0))):
##             variable_a0_plus[j]=False
##         if(PPminus[0].has_key((j,0))):
##             variable_a0_minus[j]=False        
##         ## They can not both be variable...
##         if(variable_a0_minus[j] and variable_a0_plus[j]):
##             print "Variables a(0)^+=",variable_a0_plus
##             print "Variables a(0)^-=",variable_a0_minus
##             raise ValueError,"Need to specify constant terms of the principal parts! Got: %s" % principal_parts

##     if(H._verbose>0): 
##         print "Variables a(0)^+=",variable_a0_plus
##         print "Variables a(0)^-=",variable_a0_minus
##     RHS=mpmath_ctx.matrix(int(s),int(d))
##     nvec=dict()
##     for icusp from int(0)<=icusp<nc:
##         for l in range(iMs,iMf+1):
##             nvec[l,icusp]=mpmath_ctx.mpf(l)+H.alpha(icusp)[0]
##         #print "nvec(",icusp,0,")=",nvec[0,icusp]
##     ef1=dict(); ef2=dict()
    
##     for n from iMs <=n <= iMf:        
##         for icusp in range(nc):
##             inr=mpmath_ctx.mpc(0,nvec[n,icusp])
##             for j from iQs <= j <= iQf:
##                 iargm=inr*Xm[j]
##                 ef2[n,j,icusp]=mpmath_ctx.exp(-iargm)
##         for jcusp in range(nc):
##             nr=nvec[n,jcusp]
##             nrfourpi=nr*fourpi
##             for icusp in range(nc):
##                 for j from iQs <= j <= iQf:
##                     if(not Xpb.has_key((icusp,jcusp,j))):
##                         continue
##                     iargpb=nr*mpmath_ctx.mpc(-twopi*Ypb[icusp,jcusp,j],Xpb[icusp,jcusp,j])
##                     if(nr>mpmath_ctx.eps()):
##                         ef1[j,icusp,jcusp,n]=mpmath_ctx.exp(iargpb)
##                     elif(nr<-mpmath_ctx.eps() and not H._holomorphic and H._weak):
##                         tmp=abs(nrfourpi)*Ypb[icusp,jcusp,j]
##                         tmp2=mpmath_ctx.exp(iargpb)
##                         ef1[j,icusp,jcusp,n]=mpmath_ctx.exp(iargpb)*mpmath_ctx.gammainc(kint,tmp)
##                     elif(nr<-mpmath_ctx.eps()):
##                         #print "here!"
##                         ef1[j,icusp,jcusp,n]=mpmath_ctx.mpf(0)
##                     else:
##                         ## here nr = 0
##                         if(variable_a0_plus[jcusp]):
##                             ef1[j,icusp,jcusp,n]=mpmath_ctx.mpf(1)
##                         elif(variable_a0_minus[jcusp]):
##                             ef1[j,icusp,jcusp,n]=mpmath.mp.power(Ypb[icusp,jcusp,j],kint)
##                         else:
##                             ef1[j,icusp,jcusp,n]=mpmath_ctx.mpf(0)

##     for l from  iMs <= l <= iMf:
##         for j from iQs <= j <= iQf:
##             for jcusp from int(0) <= jcusp < int(nc):
##                 lj=iMl*jcusp+l-iMs
##                 for icusp from int(0) <= icusp < int(nc):
##                     if(not Cv.has_key((icusp,jcusp,j))):
##                         continue
##                     ch= Cv[icusp,jcusp,j]  #    continue
##                     tmp1=ef1[j,icusp,jcusp,l]*ch
##                     for n from iMs <= n <= iMf:
##                         ni=iMl*icusp+n-iMs
##                         V[ni,lj]=V[ni,lj]+tmp1*ef2[n,j,icusp]
##                         if(H._verbose > 2 and n==0 and l==0):
##                             print "-------------------"
##                             print "ni,lj=",ni,lj
##                             print "icusp,jcusp=",icusp,jcusp
##                             print "V[0,0](",j,")=",V[ni,lj]
##                             print "xpb=",Xpb[icusp,jcusp,j]
##                             print "ef1(",j,")=",ef1[j,icusp,jcusp,l]
##                             print "cv(",j,")=",Cv[icusp,jcusp,j] 
##                             print "ef2(",j,")=",ef2[n,j,icusp]
##                         #if(ni==44):
##                         #    print l,j,tmp1,ef2[n,j,icusp]
##                         #    print V[ni,lj]

##     if(H._verbose>1):
##         print "Mi,Ms=",iMs,iMf
##         print "V1(",0,0,")=",V[0,0]
##     for n from 0<=n< V.rows:
##         for l from 0 <= l < V.cols:        
##             V[n,l]=V[n,l]/Qfak
##     if(H._verbose>2): 
##         for n from 0<=n< V.rows:   
##             print "V1(",n,n,")=",V[n,n]
            
##     for n from iMs<=n <=iMf:
##         for icusp from int(0)<=icusp<nc:
##             nr=nvec[n,icusp]
##             #nr=mpmath_ctx.mpf(n)
##             nrY2pi=nr*twopiY
##             ni=iMl*icusp+n-iMs
##             if(nr>mpmath_ctx.eps()):
##                 kbes=mpmath_ctx.exp(-nrY2pi)
##             elif(nr<-mpmath_ctx.eps() and (not H._holomorphic) and H._weak):
##                 kbes=mpmath_ctx.gammainc(kint,abs(nr)*fourpiY)*mpmath_ctx.exp(-nrY2pi)
##             elif(nr<-mpmath_ctx.eps()):
##                 kbes=mpmath_ctx.mpf(0)
##             else:
##                 #print "n,icusp=",n,icusp
##                 #kbes=mpmath_ctx.mpf(1)
##                 if(variable_a0_plus[icusp]):
##                     #print "1"
##                     kbes=mpmath_ctx.mpf(1)
##                 elif(variable_a0_minus[icusp]):
##                     #print "2"
##                     kbes=mpmath.mp.power(Y,kint)
##                 else:
##                     #print "3"
##                     kbes=mpmath_ctx.mpf(0)
##             V[ni,ni]=V[ni,ni]-kbes
##             # setting the right hand side of the system
##             for k in range(d):
##                 RHS[ni,k]=zero
##                 for (jcusp,l) in PPplus[k].keys():
##                     ppc=mpmath_ctx.mpc(PPplus[k][(jcusp,l)])
##                     #print "ppart_plus=",PP[k][jcusp,l]
##                     lr=mpmath_ctx.mpf(l+H.alpha(jcusp)[0])
##                     #if(abs(lr) <= mpmath.eps()):
##                     #    continue
##                     ### Note: The constant term is treated below
##                     if (ppc.ae(0)):
##                         continue
##                     summa=mpmath.mpf(0)
##                     for j from Qs <= j <= Qf:
##                         if(not Ypb.has_key((icusp,jcusp,j))):
##                             continue
##                         ch=Cv[icusp,jcusp,j]
##                         tmp=mpmath_ctx.mpc(-twopi*lr*Ypb[icusp,jcusp,j],lr*Xpb[icusp,jcusp,j]-nr*Xm[j])
##                         tmpc=ppc*mpmath_ctx.exp(tmp)
##                         summa=summa+ch*tmpc
##                     RHS[ni,k]=RHS[ni,k]+summa/Qfak
##                 if(PPplus[k].has_key((icusp,n))):
##                     #if( abs(nr) > mpmath.eps()):
##                     ppc=mpmath_ctx.mpc(PPplus[k][icusp,n])
##                     RHS[ni,k]=RHS[ni,k]-ppc*mpmath_ctx.exp(-nrY2pi)
##                 for (jcusp,l) in PPminus[k].keys():
##                     ppc_minus = mpmath_ctx.mpc(PPminus[k][(jcusp,l)])
##                     if(ppc_minus.ae(0)):
##                         continue
##                     #print "ppart_minus=",PPminus[k][jcusp,l]
##                     summa_minus=zero #summa_minus
##                     for j from Qs <= j <= Qf:
##                         if(not Ypb.has_key((icusp,jcusp,j))):
##                             continue
##                         ch=Cv[icusp,jcusp,j]
##                         tmpc_minus=mpmath.mp.power(Ypb[icusp,jcusp,j],kint)
##                         tmpc=tmpc_minus*mpmath_ctx.exp(mpmath_ctx.mpc(0,-nr*Xm[j]))
##                         summa_minus=summa_minus+ch*tmpc
##                         #summa=summa+ef2[n,j,icusp]*tmpc*ch*ppc
##                     #print "summa-(",ni,")=",summa_minus
##                     RHS[ni,k]=RHS[ni,k]+ppc_minus*summa_minus/Qfak
##                 #print "summa(",ni,")=",RHS[ni,k]
##                 if(PPminus[k].has_key((icusp,n))):
##                     ppc_minus = mpmath_ctx.mpc(PPminus[k][(icusp,n)])
##                     RHS[ni,k]=RHS[ni,k]-ppc_minus*mpmath.mp.power(Y,kint)
##                     #print "subtracting:",ppc_minus*mpmath.power(Y,kint)
##                     #print "RHS(",ni,")-=",ppc*mpmath_ctx.exp(-nrY2pi)
##     if(H._verbose>2): 
##         for n from 0<=n< V.rows:   
##             print "V2(",n,n,")=",V[n,n]
##     #for n from 0<=n< V.rows:   
##     #    print "V2(",44,n,")=",V[V.rows-1,n]
##     W=dict()
##     W['var_a+']=variable_a0_plus
##     W['var_a-']=variable_a0_minus
##     W['alphas']=H.alphas()
##     W['V']=V
##     W['space']=H
##     W['RHS']=RHS
##     W['Ms']=Ms
##     W['Mf']=Mf
##     W['Ml']=Ml
##     W['nc']=nc
##     W['rdim']=H._rdim
##     return W





## def setup_matrix_for_harmonic_Maass_waveforms_sv_bak(H,Y,int M,int Q,principal_parts):
##     r"""

##     Set up the matrix for the system of equations giving the Fourier coefficients of a Harmonic Maass waveforms.

##     INPUT:
    
##         - ``H`` -- Space of harmonic weak Maass forms
##         - ``Y`` -- height of horocycle used for sampling (mpmath real)
##         - ``k`` -- weight (mpmath real) 
##         - ``M`` -- integer
##         - ``Q`` -- integer > M
##         - ``PP``-- dict : Principal parts at cusp nr. j = \Sum_ PP(j,n) q^[-n]

##     OUTPUT:
    
##         - ``W`` -- dictionary
##         - ``W['Mf']`` -- M start
##         - ``W['nc']`` -- number of cusps
##         - ``W['Ms']`` -- M stop
##         - ``W['V']``  -- ((Ms-Mf+1)*num_cusps)**2 matrix 
##         - ``W['RHS']``  -- ((Ms-Mf+1)*num_cusps) matrix
    
    
##     EXAMPLES::

##         sage: setup_matrix_for_harmonic_Maass_waveforms_sv(MySubgroup(Gamma0(1)),mpmath.mpf(0.5),mpmath.mpf(0),10,20,{(0,-1):1})

    

##     """
##     cdef int l,j,icusp,jcusp,n,ni,li,iMl,iMs,iMf,iQs,iQf,iQl,s,nc
##     cdef double RR,YY
##     mpmath_ctx=H._mp_ctx
##     mpmath_ctx=mpmath.mp
##     #if( not test_ctx(Y,mpmath_ctx)):
##     #    raise TypeError," Need Mpmath reals as input!"
##     #GG=H._G
##     weight=H._weight
##     pi=mpmath_ctx.pi()
##     one=mpmath_ctx.mpf(1)
##     two=mpmath_ctx.mpf(2)
##     zero=mpmath_ctx.mpf(0)
##     twopi=two*pi
##     fourpi=two*twopi
##     twopiY=twopi*Y
##     fourpiY=two*twopiY
##     p=(weight-mpmath_ctx.mpf(1))/mpmath_ctx.mpf(2)
##     kint=mpmath_ctx.mpf(1)-weight
##     nc=int(H._group.ncusps())
##     two_const_terms=False
##     #if(H._weak and not H._holomorphic and not H.two_constant_terms):
##     #   
##     two_const_terms=True
##     if(Q<M):
##         Q=M+20
##     if(H._holomorphic):
##         Ms=0; Mf=M; Ml=Mf-Ms+1
##     else:
##         ## In general we have two zero-terms in the weak harmonic forms
##         if(two_const_terms):
##             Ms=-M; Mf=M; Ml=Mf-Ms+1
##         else:
##             Ms=-M; Mf=M; Ml=Mf-Ms+1
##     Qs=1-Q; Qf=Q
##     Qfak=mpmath_ctx.mpf(2*Q)
##     iMl=int(Ml); iMs=int(Ms); iMf=int(Mf)
##     iQf=int(Qf); iQs=int(Qs); iQl=int(Qf-Qs+1)
##     # Pullpack points
##     RF=RealField(mpmath.mp.prec)
##     pb=pullback_pts_mpc_new(H,Qs,Qf,RF(Y))
##     Xm=pb['xm']; Xpb=pb['xpb']; Ypb=pb['ypb']; Cv=pb['cvec']
##     #[Xm,Xpb,Ypb,Cv]=pullback_pts_hw(H,Qs,Qf,Y)
##     s=int(nc*iMl)
##     #print "Q=",Q
##     #print "kint=",kint
##     print "s=",s
##     V=mpmath_ctx.matrix(int(s),int(s))
##     if(isinstance(principal_parts,dict)):
##         if(principal_parts.has_key('+')):
##             PP = principal_parts['+']
##         else:
##             PP = [{}]
##         if(principal_parts.has_key('-')):
##             PPminus = principal_parts['-']
##         else:
##             PPminus = list()
##             for i in range(len(PP)):
##                 PPminus.append({})
##     else:
##         PP = list()
##         PPminus = list()
##         for i in range(len(principal_parts)):
##             PP.append(principal_parts[i]['+'])
##             PPminus.append(principal_parts[i]['-'])
            
##             # if we have a holomorphic function we  don't want any non-holomorphic principal parts
##     if(H._holomorphic):
##         print "H is holomorphic!"
##         for i in range(len(PP)):
##             for j in range(0,nc):
##                 PPminus[i][(j,0)]=0

##     if(not isinstance(PP,list)):
##         raise ValueError,"Need a list of (+) principal parts!"
##     if(not isinstance(PPminus,list)):
##         raise ValueError,"Need a list of (-) principal parts!"
##     if(len(PPminus)<>len(PP)):
##         raise ValueError," Need equal number of holomorphic and non-holomorphic principal parts! Got:%s" % principal_parts
##     d = len(PP)
##     if(H._verbose):
##         print "s=",s
##         print "PP=",PP
##         print "PPminus=",PPminus
##     ## From the principal parts we can figure out which of a^+(0) and a^-(0) are unknowns
##     ## Version one: we need to set at least one of them for each cusp.
##     variable_a0_plus=dict()
##     variable_a0_minus=dict()
##     for j in range(nc):
##         if(H.alpha(j)[1]<>1):
##             #variable_a0_plus[j]=False
##             #variable_a0_minus[j]=False
##             continue
##         variable_a0_plus[j]=True
##         variable_a0_minus[j]=True
##         if(PP[0].has_key((j,0))):
##             variable_a0_plus[j]=False
##         if(PPminus[0].has_key((j,0))):
##             variable_a0_minus[j]=False
            
##         ## They can not both be variable...
##         if(variable_a0_minus[j] and variable_a0_plus[j]):
##             print "Variables a(0)^+=",variable_a0_plus
##             print "Variables a(0)^-=",variable_a0_minus
##             raise ValueError,"Need to specify constant terms of the principal parts! Got: %s" % principal_parts
##     if H._verbose>0:
##         print "Variables a(0)^+=",variable_a0_plus
##         print "Variables a(0)^-=",variable_a0_minus
##     RHS=mpmath_ctx.matrix(int(s),int(d))
##     nvec=dict()
##     for icusp from int(0)<=icusp<nc:
##         for l in range(iMs,iMf+1):
##             nvec[l,icusp]=mpmath_ctx.mpf(l)+H.alpha(icusp)[0]
##             #print "nvec(",icusp,0,")=",nvec[0,icusp]
##     ef1=dict() #; ef1m=dict();
##     ef2=dict() #; ef2m=dict()
    
##     for n from iMs <=n <= iMf:        
##         for icusp in range(nc):
##             inr=mpmath_ctx.mpc(0,nvec[n,icusp])
##             for j from iQs <= j <= iQf:
##                 iargm=inr*Xm[j]
##                 ef2[n,j,icusp]=mpmath_ctx.exp(-iargm)
##         for jcusp in range(nc):
##             nr=nvec[n,jcusp]
##             nrfourpi=nr*fourpi
##             for icusp in range(nc):
##                 for j from iQs <= j <= iQf:
##                     if(not Xpb.has_key((icusp,jcusp,j))):
##                         continue
##                     iargpb=nr*mpmath_ctx.mpc(-twopi*Ypb[icusp,jcusp,j],Xpb[icusp,jcusp,j])
##                     if(nr>mpmath_ctx.eps()):
##                         ef1[j,icusp,jcusp,n]=mpmath_ctx.exp(iargpb)
##                     elif(nr<-mpmath_ctx.eps() and not H._holomorphic and H._weak):
##                         tmp=abs(nrfourpi)*Ypb[icusp,jcusp,j]                        
##                         ef1[j,icusp,jcusp,n]=mpmath_ctx.exp(iargpb)*mpmath_ctx.gammainc(kint,tmp)
##                     elif(nr<-mpmath_ctx.eps()):
##                         ef1[j,icusp,jcusp,n]=zero
##                     else: # nr = 0 
##                         ## WAs commented out
##                         if(variable_a0_plus[jcusp]):
##                             ef1[j,icusp,jcusp,n]=one
##                         elif(variable_a0_minus[jcusp]):
##                             ef1[j,icusp,jcusp,n]=mpmath_ctx.power(Ypb[icusp,jcusp,j],kint)
##                         else:
##                             ef1[j,icusp,jcusp,n]=zero #one #zero

##     for l from  iMs <= l <= iMf:
##         for j from iQs <= j <= iQf:
##             for jcusp from int(0) <= jcusp < int(nc):
##                 lj = iMl*jcusp+iMf+l
##                 #print "lj_m,lj_p=",lj_m,lj_p
##                 for icusp from int(0) <= icusp < int(nc):
##                     if(not Cv.has_key((icusp,jcusp,j))):
##                         continue
##                     ch= mpmath.mp.mpc(Cv[icusp,jcusp,j].real(),Cv[icusp,jcusp,j].imag())
##                     if ch==0:
##                         continue
##                     tmp1=ef1[j,icusp,jcusp,l]*ch
##                     if(tmp1.ae(zero)):
##                         continue
##                     for n from iMs <= n <= iMf:                        
##                         ni=iMl*icusp+iMf+n #-iMs
##                         V[ni,lj] = V[ni,lj] + tmp1*ef2[n,j,icusp]
##                         if(H._verbose > 1 and ((n==0 and icusp==2) or (l==0 and  jcusp==2)) and (j>0 and j<5)):
##                             print "-------------------"
##                             print "ni,lj=",ni,lj
##                             print "nrp=",nvec[n,jcusp]
##                             #print "nrm=",nvec[-n,jcusp]
##                             print "icusp,jcusp=",icusp,jcusp
##                             print "V[0,0](",j,")=",V[ni,lj]
##                             #print "xpb=",Xpb[icusp,jcusp,j]
##                             print "tmp1(",j,")=",tmp1
##                             #print "tmp1p(",j,")=",tmp1p
##                             print "cv(",j,")=",Cv[icusp,jcusp,j] 
##                             print "ef2(",j,")=",ef2[n,j,icusp]
##                             #print "ef2p",j,")=",ef2[n,j,icusp]
##                         #if(ni==44):
##                         #    print l,j,tmp1,ef2[n,j,icusp]
##                         #    print V[ni,lj]
##     if H._verbose>0:
##         print "Mi,Ms=",iMs,iMf
##     if(H._verbose>1): 
##         print "V1(",0,0,")=",V[0,0]
##     for n from 0<=n< V.rows:
##         for l from 0 <= l < V.cols:        
##             V[n,l]=V[n,l]/Qfak
##     if(H._verbose>2): 
##         for n from 0<=n< V.rows:   
##             print "V1(",n,n,")=",V[n,n]
            
##     for n from iMs<= n <=iMf:
##         for icusp from int(0)<=icusp<nc:
##             nr=nvec[n,icusp]
##             nrY2pi=nr*twopiY
##             ni=iMl*icusp+iMf+n #-iMs
##             if(nr>mpmath_ctx.eps()):
##                 kbes=mpmath_ctx.exp(-nrY2pi)
##             elif(nr<-mpmath_ctx.eps() and not H._holomorphic and H._weak):                
##                 kbes=mpmath_ctx.gammainc(kint,abs(nr)*fourpiY)*mpmath_ctx.exp(-nrY2pi)
##             elif(nr<-mpmath_ctx.eps()):
##                 kbes=zero
##             else: # nr = 0
##                 if(variable_a0_plus[icusp]):
##                     kbes=one
##                 elif(variable_a0_minus[icusp]):
##                     kbes=mpmath.power(Y,kint)
##                 else:
##                     kbes=zero
##             V[ni,ni]=V[ni,ni]-kbes
##             #V[ni_p,ni_p]=V[ni_p,ni_p]-kbes_p
##             ## Now we have to add the principal part, including any non-variable a^+(0) and a^-(0)
##             for k in range(d):
##                 RHS[ni,k]=zero
##                 summa=zero
##                 for (jcusp,l) in PP[k].keys():
##                     ppc=mpmath_ctx.mpc(PP[k][(jcusp,l)])
##                     lr=mpmath_ctx.mpf(l+H.alpha(jcusp)[0])
##                     if(ppc.ae(0)):
##                         continue
##                     for j from Qs <= j <= Qf:
##                         if(not Ypb.has_key((icusp,jcusp,j))):
##                             continue
##                         ch=Cv[icusp,jcusp,j]
##                         tmp=mpmath_ctx.mpc(-twopi*lr*Ypb[icusp,jcusp,j],lr*Xpb[icusp,jcusp,j]-nr*Xm[j])
##                         summa=summa+ch*ppc*mpmath_ctx.exp(tmp)
##                 RHS[ni,k]=RHS[ni,k]+summa/Qfak
##                 if(PP[k].has_key((icusp,n))):
##                     if( abs(nr) > mpmath.eps()):
##                         ppc=mpmath_ctx.mpc(PP[k][icusp,n])
##                         RHS[ni,k]=RHS[ni,k]-ppc*mpmath_ctx.exp(-nrY2pi)
##                 ## Was commented out
##                 summa=zero
##                 for jcusp in range(nc):
##                     if(H.alpha(jcusp)[1]<>1):
##                         continue
##                     ppc_minus = zero                    
##                     ppc_plus  = zero                    
##                     if( (not variable_a0_minus[jcusp]) and PPminus[k].has_key((jcusp,0))):
##                         ppc_minus = mpmath_ctx.mpc(PPminus[k][(jcusp,0)])
##                         if(not ppc_minus.ae(zero)):
##                             summa_minus=zero
##                             for j from Qs <= j <= Qf:
##                                 if(not Ypb.has_key((icusp,jcusp,j))):
##                                     continue
##                                 ch=mpmath.mp.mpc(Cv[icusp,jcusp,j].real(),Cv[icusp,jcusp,j].imag())
##                                 tmpc_minus=mpmath.power(Ypb[icusp,jcusp,j],kint)
##                                 tmpc=tmpc_minus*mpmath_ctx.exp(mpmath_ctx.mpc(0,-nr*Xm[j]))
##                                 summa_minus=summa_minus+ch*tmpc
##                             summa=summa+summa_minus*ppc_minus
##                     if( (not variable_a0_plus[jcusp]) and PP[k].has_key((jcusp,0))):
##                         ppc_plus = mpmath_ctx.mpc(PP[k][(jcusp,0)])                            
##                         summa_plus=zero
##                         if(not ppc_plus.ae(zero)):
##                             for j from Qs <= j <= Qf:
##                                 if(not Ypb.has_key((icusp,jcusp,j))):
##                                     continue
##                                 ch=Cv[icusp,jcusp,j]
##                                 tmpc=mpmath_ctx.exp(mpmath_ctx.mpc(0,-nr*Xm[j]))
##                                 summa_plus=summa_plus+ch*tmpc
##                             summa=summa+summa_plus*ppc_plus
##                          #summa=summa+ef2[n,j,icusp]*tmpc*ch*ppc
##                 RHS[ni,k]=RHS[ni,k]+summa/Qfak
##                 if(n==0 and H.alpha(icusp)[1]==1):
##                     if(not variable_a0_plus[icusp]):
##                         if(PP[k].has_key((icusp,0))):
##                             RHS[ni,k]=RHS[ni,k]-mpmath_ctx.mpc(PP[k][(icusp,0)])
##                     if(not variable_a0_minus[icusp]):
##                         if(PPminus[k].has_key((icusp,0))):
##                             ppc_minus = mpmath_ctx.mpc(PPminus[k][(icusp,0)])
##                             RHS[ni,k]=RHS[ni,k]-ppc_minus*mpmath.power(Y,kint)


                        
##                 # print "summa(",ni,")=",RHS[ni,k]
##     if(PPminus[k].has_key((icusp,n))):
##                     ppc_minus = mpmath_ctx.mpc(PPminus[k][(icusp,n)])
##                     RHS[ni,k]=RHS[ni,k]-ppc_minus*mpmath.power(Y,kint)
##                     print "subtracting:",ppc_minus*mpmath.power(Y,kint)
##                     #print "RHS(",ni,")-=",ppc*mpmath_ctx.exp(-nrY2pi)
##                 ## until here
##     if(H._verbose>2): 
##         for n from 0<=n< V.rows:   
##             print "V2(",n,n,")=",V[n,n]
##     #for n from 0<=n< V.rows:   
##     #    print "V2(",44,n,")=",V[V.rows-1,n]
##     W=dict()
##     W['space']=H
##     W['var_a+']=variable_a0_plus
##     W['var_a-']=variable_a0_minus
##     W['alphas']=H.alphas()
##     W['V']=V
##     W['RHS']=RHS
##     W['Ms']=Ms 
##     W['Mf']=Mf
##     W['Ml']=Ml
##     W['nc']=nc
##     W['rdim']=H._rdim
##     return W

cpdef check_principal_parts(H,principal_parts):
    r"""
    Make sure that the principal parts are of correct format and extract necessary information.
"""
    res = {}
    PPplus=[]; PPminus=[]
    cdef int nc = H.group().ncusps()
    for pp in principal_parts:
        PPplus.append(pp['+'])
        PPminus.append(pp['-'])

    if H._verbose>1:
        print "PPplus=",PPplus
        print "PPminus=",PPminus
    for pp in PPplus:
        for icusp,n in pp:
            if((icusp not in range(nc)) or n>0):
                raise ValueError,"Principal parts are in wrong format! PP+={0}".format(PPplus)
    for pp in PPminus:
        for icusp,n in pp:
            if((icusp not in range(nc)) or n<0):
                raise ValueError,"Principal parts are in wrong format! PP-={0}".format(PPminus)            

    # We include functions which are almost holomorphic and
    # only non-holomorphic parts are the principal parts
    #if H._holomorphic:
    #    for i in range(len(PPminus)):
    #        for j in range(0,nc):
    #            PPminus[i][(j,0)]=0
    if not isinstance(PPplus,list):
        raise ValueError,"Need a list of (+) principal parts!"
    if not isinstance(PPminus,list):
        raise ValueError,"Need a list of (-) principal parts!"
    if len(PPminus)<>len(PPplus):
        raise ValueError," Need equal number of holomorphic and non-holomorphic principal parts! Got:%s" % principal_parts
    if H._verbose>0: 
        print "PPplus=",PPplus
        print "PPminus=",PPminus
    variable_a0_plus=dict()
    variable_a0_minus=dict()
    for j in range(nc):
        if H.alpha(j)[1]<>1:
            #variable_a0_plus[j]=False
            #variable_a0_minus[j]=False
            continue
        variable_a0_plus[j]=True
        variable_a0_minus[j]=True
        if PPplus[0].has_key((j,0)):
            variable_a0_plus[j]=False
        if PPminus[0].has_key((j,0)):
            variable_a0_minus[j]=False        
        ## They can not both be variable...
        if variable_a0_minus[j] and variable_a0_plus[j]:
            print "Variables a(0)^+=",variable_a0_plus
            print "Variables a(0)^-=",variable_a0_minus
            raise ValueError,"Need to specify constant terms of the principal parts! Got: %s" % principal_parts

    if H._verbose>0: 
        print "Variables a(0)^+=",variable_a0_plus
        print "Variables a(0)^-=",variable_a0_minus
        print "alphas=",H.alphas()
    ### In the future I want to be able to use many right-hand sides at once.
    res['variable_a0_plus']={0:variable_a0_plus}
    res['variable_a0_minus']={0:variable_a0_minus}
    res['PPplus']=PPplus
    res['PPminus']=PPminus
    return res


cdef int mpc_is_zero(mpc_t z):
    if mpfr_zero_p(mpc_realref(z))==0:
        return 0
    if mpfr_zero_p(mpc_imagref(z))==0:
        return 0
    return 1


cpdef setup_and_solve_for_harmonic_Maass_waveforms(H,Y_in,int M,int Q,principal_parts,version=1,cset=[]):
    cdef Matrix_complex_dense V1,RHS1
    W = setup_matrix_for_harmonic_Maass_waveforms(H,Y_in,M,Q,principal_parts,version)
    V1 = W['V']
    RHS1 = W['RHS']
    N = H.set_norm(cset)
    N['Ms']=W['Ms']
    N['Mf']=W['Mf']
    N['nc']=W['nc']
    N['PP']=W['PP']
    N['space']=H
    N['alphas']=W['alphas']
    N['Ml']=W['Ml']
    N['var_a+']=W['var_a+']
    N['var_a-']=W['var_a-']
    C = solve_system_for_harmonic_weak_Maass_waveforms_mp(N,V1,RHS1) 
    return C
    
cpdef solve_system_for_harmonic_weak_Maass_waveforms_mp(dict N, Matrix_complex_dense V1,Matrix_complex_dense RHS1,gr=0):
    r"""
    Solve the linear system to obtain the Fourier coefficients of Maass forms

    INPUT:

    - ``W`` --   (system) dictionary
        - ``W['Ms']``  -- M start
        - ``W['Mf']``  -- M stop
        - ``W['nc']``  -- number of cusps
        - ``W['space']``  -- space of automorphic forms
        - ``W['V']``   -- matrix of size ((Ms-Mf+1)*nc)**2
        - ``W['RHS']`` -- right hand side (for inhomogeneous system) matrix of size ((Ms-Mf+1)*nc)*(dim)
    - ``N`` -- normalisation (dictionary, output from the set_norm_for_maass function)
        - ``N['SetCs']``   -- Which coefficients are set and their values
        - ``N['comp_dim']``-- How large is the assumed dimension of the solution space
        - ``N['num_set']`` -- Number of coefficients which are set
        

    OUTPUT:
    
    - ``C`` -- Fourier coefficients

    EXAMPLES::

        sage: G=MySubgroup(Gamma0(1))
        sage: mpmath.mp.dps=20
        sage: R=mpmath.mpf(9.533695261353557554344235235928770323821256395107251982375790464135348991298347781769255509975435366)
        sage: Y=mpmath.mpf(0.5)
        sage: W=setup_matrix_for_Maass_waveforms(G,R,Y,12,22)
        sage: N=set_norm_maass(1)
        sage: C=solve_system_for_Maass_waveforms(W,N)
        sage: C[0][2]*C[0][3]-C[0][6]
        mpc(real='-1.8055426724989656270259e-14', imag='1.6658248366482944572967e-19')

    If M is too large and the precision is not high enough the matrix might be numerically singular

        W=setup_matrix_for_Maass_waveforms(G,R,Y,20,40)  
        sage: C=solve_system_for_Maass_waveforms(W,N)
        Traceback (most recent call last)
        ...
        ZeroDivisionError: Need higher precision! Use > 23 digits!

    Increasing the precision helps
    
        sage: mpmath.mp.dps=25
        sage: R=mpmath.mpf(9.533695261353557554344235235928770323821256395107251982375790464135348991298347781769255509975435366)
        sage: C=solve_system_for_Maass_waveforms(W,N)
        sage: C[0][2]*C[0][3]-C[0][6]
        mpc(real='3.780824715556976438911480324e-25', imag='2.114746048869188750991752872e-99')


        """
    cdef Matrix_complex_dense RHS,LHS,Q,R
    cdef Vector_complex_dense v,TMP
    cdef int j,r,k,n,nr,Ml,cr,Ms,Mf,prec,use_sym,num_set,roffs,coffs,fn_j,maxit,comp_dim,c
    cdef MPComplexNumber zero,tmp,vv
    #V=N['V']
    Ms=N['Ms']
    Mf=N['Mf']
    nc=N['nc']
    PP=N['PP']
    H = N['space']
    cdef int verbose = H._verbose
    alphas=N['alphas']
    Ml=N['Ml'] #Mf-Ms+1
    variable_a_plus=N['var_a+']
    variable_a_minus=N['var_a-']
    cdef int nrows = int(V1.nrows())
    cdef int ncols = int(V1.ncols())
    if ncols<>Ml*nc or nrows<>Ml*nc:
        raise Exception," Wrong dimension of input matrix!"
    # we have to assume that all normalizations use the same coefficients
    prec = H.prec()
    maxit=1000
    SetCs=N['SetCs']
    SetCs_neg=N.get('SetCs_neg',{})
    CF = MPComplexField(prec)
    zero = CF(0)
    tmp=CF(0)
    comp_dim=N['comp_dim']
    use_sym=0



    SetClist=dict()
    for j in range(0,comp_dim):
        SetClist[j]=dict()
    if len(PP)>0 and ((comp_dim<>len(SetCs.keys()) and comp_dim<>len(PP))):
        print "comp_dim=",comp_dim
        print "SetC=",SetCs
        print "PP=",PP
        raise ValueError," Inconsistent normalization SetCs:%s" % SetCs
    num_set=0
    for j in range(0,comp_dim):
        # # First we treat set values of coefficients not corresponsing to the principal part
        for (r,n) in SetCs[j].keys():
            for j in range(comp_dim):
                nr = r*Ml+n
                if nr>=0 or not H.is_holomorphic():
                    SetClist[j][nr]=SetCs[j][(r,n)]
        if verbose>0:
            print "SetClist_pos=",SetClist
        ## Then we check the zeroth coefficients
        ## Note that these should already be in the right hand side.    
        for r in range(nc):
            if(alphas[r][1]==1):
                if( (not variable_a_plus[j][r]) and (not variable_a_minus[j][r])):
                    nr = r*Ml
                    #if(SetCs_neg.get(j,{}).has_key((r,0))):
                    SetClist[j][nr]=CF(0) #SetCs_neg[j][(r,0)]) 
        num_set=len(SetClist[0].keys())
    if verbose>0:
        print "SetClist_tot=",SetClist
    #use_symmetry=False

    MS = MatrixSpace(CF,int(Ml*nc-num_set),int(comp_dim))
    RHS = Matrix_complex_dense(MS,0,True,True)

    # We allow for either a variation of principal parts or of set coefficients
    #if(W.has_key('RHS')):
    l=RHS1.ncols()
    if(l>1 and l<>comp_dim):
        raise ValueError,"Incorrect number of right hand sides!"
        
    MS2 = MatrixSpace(CF,int(Ml*nc-num_set),int(Ml*nc-num_set))
    LHS = Matrix_complex_dense(MS2,0,True,True)
    #LHS=mpmath_ctx.matrix(int(Ml*nc-num_set),int(Ml*nc-num_set))
    roffs=0
    cdef MPComplexNumber tmpc
    tmpc = CF(0)
    cdef mpc_t ** setc_values=NULL
    cdef int** setc_n=NULL
    setc_values = <mpc_t**>check_allocarray(sizeof(mpc_t*), comp_dim)
    setc_n =<int**>check_allocarray(sizeof(int*), comp_dim)
    if setc_values==NULL: raise MemoryError
    if setc_n==NULL: raise MemoryError        
    for r in range(comp_dim):
        setc_values[r]=NULL
        setc_values[r] = <mpc_t*>check_allocarray(sizeof(mpc_t), num_set)
        if setc_values[r]==NULL: raise MemoryError
        setc_n[r]=NULL
        setc_n[r] = <int*>check_allocarray(sizeof(int), num_set)
        if setc_n[r]==NULL: raise MemoryError        
        i = 0
        for j in SetClist[r].keys():
            setc_n[r][i]=int(j)
            i+=1            
        assert i==num_set
        for j in range(num_set):            
            mpc_init2(setc_values[r][j],prec)
            tmpc = CF(SetClist[r][setc_n[r][j]])
            mpc_set(setc_values[r][j],tmpc.value,rnd)
    if verbose>0:
        print "Ml=",Ml
        print "num_set=",num_set
        print "SetCs=",SetCs
        print "SetClist=",SetClist
        #print "Valslist=",Valslist
        print "V.rows=",nrows
        print "V.cols=",ncols
        print "LHS.rows=",LHS.nrows()
        print "LHS.cols=",LHS.ncols()
        print "RHS.rows=",RHS1.nrows()
        print "RHS.cols=",RHS1.ncols()
        print "use_sym=",use_sym
    cdef int isin = 0
    for r in range(nrows):
        cr=r+Ms
        isin = 0
        for j in range(num_set):
            if setc_n[0][j]==r+Ms:
                isin=1
                break       
        #if(SetClist[0].keys().count(r+Ms)>0):
        if isin==1:
            roffs=roffs+1
            continue
        for fn_j in range(comp_dim):
            ## Assume well-formed...
            mpc_set(RHS._matrix[r-roffs][fn_j],RHS1._matrix[r][fn_j],rnd)
            mpc_neg(RHS._matrix[r-roffs][fn_j],RHS._matrix[r-roffs][fn_j],rnd)
            
            #for c in SetClist[fn_j].keys():
            for j in range(num_set):
                #vv=CF(SetClist[fn_j][c])
                #tmp=vv*V[r,c-Ms]
                #tmp = CF(V[r,c-Ms])
                c = setc_n[fn_j][j]
                mpc_mul(tmp.value,V1._matrix[r][c-Ms],setc_values[fn_j][j],rnd)
                mpc_sub(RHS._matrix[r-roffs][fn_j],RHS._matrix[r-roffs][fn_j],tmp.value,rnd)
                #RHS[r-roffs,fn_j]=RHS[r-roffs,fn_j]-tmp
        coffs=0
        for k in range(ncols):            
            isin=0
            for j in range(num_set):
                if setc_n[0][j]==k+Ms:
                    isin=1
                    break           
            #if(SetClist[0].keys().count(k+Ms)>0):
            if isin==1:
                coffs=coffs+1
                continue
            #try:
            #tmp = CF(V[r][k])
            mpc_set(LHS._matrix[r-roffs][k-coffs],V1._matrix[r][k],rnd)
    if gr==1:
        return LHS,RHS
    smin=smallest_inf_norm(LHS)
    if verbose>0:
        print "sminfn=",smin
    dps0=CF.prec()
    done=False
    i=1
    while (not done and i<=maxit):
        try:
            Q,R = LHS.qr_decomposition()
            #A, p = mpmath_ctx.LU_decomp(LHS)
            done=True
        except ZeroDivisionError:
            #t=int(mpmath_ctx.ceil(-mpmath_ctx.log10(smallest_inf_norm(LHS))))
            t=int(ceil(-log_b(smallest_inf_norm(LHS),10)))
            dps=t+5*i; i=i+1
            if verbose>-1:
                print "raising number of digits to:",dps
            LHS.set_prec(dps)
            # raise ZeroDivisionError,"Need higher precision! Use > %s digits!" % t
    if(i>=maxit):
        raise ZeroDivisionError,"Can not raise precision enough to solve system! Should need > %s digits! and %s digits was not enough!" % (t,dps)
    X=dict()
    for fn_j in range(comp_dim):
        X[fn_j] = dict() 
        X[fn_j][0] = dict() 
        v = RHS.column(fn_j)
        if verbose>0:
            print "len(B)=",len(v)
            #print "RHS=",v
        #b = mpmath_ctx.L_solve(A, RHS.column(fn_j), p)
        TMP = LHS.solve(v) #mpmath_ctx.U_solve(A, b)
        roffs=0
        res = (LHS*TMP-v).norm()
        if verbose>0:
            print "res(",fn_j,")=",res
        #res = mpmath_ctx.norm(mpmath_ctx.residual(LHS, TMP, RHS.column(fn_j)))
        #print "res(",fn_j,")=",res
        for i in range(0,nc):
            X[fn_j][0][i]=dict()
        for i in range(nc):
            roffs2=0
            for n in range(Ml):
                nn=i*Ml+n+Ms
                key=n+Ms
                #if(i==1):
                #    print n,key
                if(SetClist[fn_j].keys().count(nn)>0):
                    if verbose>1:
                        print "We have set ",nn
                    roffs=roffs+1
                    X[fn_j][0][i][key]=SetClist[fn_j][nn] 
                    if verbose>0:
                        print "X[",fn_j,",",i,",",key,"]=",SetClist[fn_j][nn]
                        print "nn=",nn
                    continue
                try:
                    #X[fn_j][0][i][n-roffs2+Ms]=TMP[nn-Ms-roffs,0]
                    X[fn_j][0][i][key]=TMP[nn-Ms-roffs]
                except IndexError:
                    print "n*Mli-roffs=",n,'+',Ml,'*',i,'-',roffs,"=",n+Ml*i-roffs
                ## We also insert the principal part if it is applicable
    #mpmath.mp.dps=dpold
    # return x
    for r in range(comp_dim):
        if setc_values[r]<>NULL:
            for j in range(num_set):
                mpc_clear(setc_values[r][j])
            sig_free(setc_values[r])
        if setc_n[r]<>NULL:
            sig_free(setc_n[r])
    if setc_n<>NULL:
        sig_free(setc_n)
    if setc_values<>NULL:
        sig_free(setc_values)
    return X
