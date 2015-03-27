#cython: profile=False
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
#include "interrupt.pxi"  # ctrl-c interrupt block support
#include "stdsage.pxi"  # ctrl-c interrupt block support
include 'sage/ext/stdsage.pxi'
include "sage/ext/cdefs.pxi"
include 'sage/ext/interrupt.pxi'
#include "sage/ext/gmp.pxi"
#include "sage/rings/mpc.pxi"
from psage.rings.mpfr_nogil cimport *

cdef mpc_rnd_t rnd
cdef mpfr_rnd_t rnd_re
rnd = MPC_RNDNN
rnd_re = MPFR_RNDN
from sage.rings.complex_mpc cimport MPComplexNumber
from sage.rings.complex_mpc import MPComplexField
from sage.rings.real_mpfr cimport RealNumber,RealField_class
from sage.rings.real_mpfr import RealField
from sage.rings.complex_number cimport ComplexNumber
from sage.rings.rational_field import QQ #,Rational
from sage.rings.rational cimport Rational
from sage.rings.integer cimport Integer as cInteger
from sage.rings.integer import Integer,is_Integer

from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.arith import gcd,kronecker,valuation
from sage.modules.free_module import FreeModule
from sage.misc.functional import numerator,denominator,is_odd
from sage.functions.other import floor,ceil as pceil
from sage.modular.arithgroup.congroup_sl2z import SL2Z
from sage.all import MatrixSpace,pi,ComplexField, vector,matrix,copy

##
##
from sage.libs.mpmath.ext_impl cimport MPF_to_tuple
#from sage.libs.mpmath.utils cimport mpfr_from_mpfval

#from maass_forms_alg import *
import sys
#from sys.stdout import.flush
from psage.rings.mpc_extras cimport *
from psage.matrix.matrix_complex_dense cimport Matrix_complex_dense
from psage.modform.arithgroup.mysubgroups_alg cimport pullback_to_psl2z_mat_c,pullback_to_psl2z_mat,_apply_sl2z_map_mpfr
from psage.modules.vector_complex_dense cimport Vector_complex_dense
from psage.modules.vector_real_mpfr_dense cimport Vector_real_mpfr_dense
#from inc_gamma cimport incgamma_hint_c
from psage.functions.inc_gamma import incgamma_hint  ## gamma(n+1/2,x)
from psage.functions.inc_gamma cimport incgamma_hint_c,incgamma_nint_c,incgamma_pint_c
import mpmath
#from vv_harmonic_weak_maass_forms import rn_from_D
#from sage.all import cached_method,cached_function
import cython
from Cython import Utils
from Cython.Utils import cached_function, cached_method
import mpmath.matrices.matrices
mpmath.matrices.matrices.matrix = mpmath.matrix
  
cdef extern from "math.h":
    double fabs(double)
    double fmax(double,double)
    int ceil(double) 

def vv_harmonic_wmwf_setupV(WR,PP,Y_in,M,Q,kappa,sym_type=1,verbose=0):
    r"""
    Set up the matrix for a vector-valued harmonic weak maass form for the Weil representation.
    
    EXAMPLES::

    sage: WR=WeilRepDiscriminantForm(1,False)
    sage: Y=mpmath.mp.mpf(0.5); k=mpmath.mp.mpf(0.5)
    sage: M=10; Q=20; m=-1; beta=1/2
    sage: W= vv_harmonic_wmwf_setupV(WR,beta,m,Y,M,Q,k)


    """
    raise NotImplementedError,"This method is obsolete!"


def vv_harmonic_wmwf_setupV_ef(H,PP,Y_in,M,Q,kappa,sym_type=1,verbose=0):
    r"""
    Set up the matrix for a vector-valued harmonic weak maass form for the Weil representation.
    
    EXAMPLES::

    sage: WR=WeilRepDiscriminantForm(1,False)
    sage: Y=mpmath.mp.mpf(0.5); k=mpmath.mp.mpf(0.5)
    sage: M=10; Q=20; m=-1; beta=1/2
    sage: W= vv_harmonic_wmwf_setupV(WR,beta,m,Y,M,Q,k)


    """
    X = H.multiplier()
    WM = X.weil_module()
    cdef list D=X.D()   #; N=WR.rank()/2
    cdef list oddD = WM.odd_submodule(indices=1)
    cdef int nc = H.group().ncusps()
    cdef int orderD = WM.finite_quadratic_module().order()
    cdef int tt=-1
    betai=dict();mbeta=dict(); mbetai=dict(); mm=dict(); mmm=dict()
    dual_rep = H.multiplier().is_dual()
    if verbose>0:
        print "in vv_harmonic_wmwf_setupV_ef!"
        print "args=",WM,PP,Y_in,M,Q,kappa,sym_type,verbose
        print "D=",D
    PP0=dict()
    
    for t in PP:
        if isinstance(t,tuple):
            (beta,m) = t
        else:
            (beta,m) = rn_from_D(H.multiplier(),t)
        if isinstance(beta,int) or isinstance(beta,Integer):
            i = beta
            beta = D[i]
            betai[beta] = i
        else:
            betai[beta]=D.index(beta)
        PP0[(beta,m)]=PP[t]
        p=numerator(1-beta); q=denominator(1-beta)
        mb=QQ(p % q)/QQ(q)
        mbeta[beta]=mb
        if verbose>0:
            print "b=",beta
            print "1-b=",mb
        mbetai[beta]=D.index(mb)
        if not beta in D:
            raise Exception,"Need beta=%s in D=" %(beta,D)
        mm [(beta,m)]=(m+H.multiplier().Qv[betai[beta]])
        mmm[(beta,m)]=mpmath.mp.mpf(mm[(beta,m)])
        if mm[(beta,m)]>0:
            raise ValueError,"Need principal part with negative exponents! Got:%s for PP:%s" %(mm[(beta,m)],PP)
        if verbose>0:
            print "-beta=",mbeta[beta]
            print "mm=",mm[(beta,m)]
        if mm[(beta,m)]>tt:
            tt=mm[(beta,m)]
    if tt>=0 or H._holomorphic: # Holomorphic
        Ms=int(0); Mf=int(M); Ml=int(Mf-Ms+1)
    else:
        Ms=int(-M); Mf=int(M); Ml=int(Mf-Ms+1)
    abD=len(D)
    if verbose>0:
        print "Ms,Mf=",Ms,Mf
        print "betai=",betai
        print "mbeta=",mbeta
        print "mbetai=",mbetai
        print "mm=",mm
        print "mmm=",mmm
    

    Y=mpmath.mp.mpf(Y_in)
    mpone=mpmath.mp.mpf(1)
    mptwo=mpmath.mp.mpf(2); mpfour=mpmath.mp.mpf(4)
    twopi=mptwo*mpmath.pi(); twopii=mpmath.mp.mpc(0,twopi)
    fourpi=mptwo*twopi
    weight=mpmath.mp.mpf(kappa); weight_over_two=weight/mpmath.mp.mpf(2)
    mp0=mpmath.mpf(0)
    if Q==0:
        Q=M+20
    Qs=int(1-Q); Qf=int(Q)
    Qfak=mpmath.mp.mpf(2*Q)
    fak1=dict(); fak2=dict()
    #besarg_old=-1
    kint=mpone-weight
    ## Using symmetry and we also figure out the symmetry type
    sym_type=get_sym_type(X,kappa)
    if verbose>0:
        print "Symmetry type=",sym_type
    ## Recall that we have different index sets depending on the symmetry
    if sym_type==1:
        Dstart=D[0]; Dfinish=D[-1] # 0,1,...,N
    elif sym_type==-1:
        Dstart=D[0]; Dfinish=D[-1] #len(D)+1 # 1,2,...,N-1  (since -0=0 and -N=N)
    # Pullpack points
    [Xm,Xpb,Ypb,Cv,Cfak]=pullback_pts_weil_rep_ef(X,Q,Y,weight,Dstart,Dfinish,dual_rep)
    Dl=len(D) #Dfinish-Dstart+1
    #setD=range(Dstart,Dfinish+1)
    s=Dl*Ml
    if verbose>0:
        print "Dstart, Dfinish=",Dstart,Dfinish,Dl
        print "Ml=",Ml
        print "D=",D
        print "abD=",abD
        print "matrix size=",s
    V=mpmath.matrix(int(s),int(s))#,force_type=mpc)
    RHS=mpmath.mp.matrix(int(s),int(1))
    #p=(weight-mpone)/mptwo
    if verbose>0:
        print "k/2=",weight_over_two
        print "p=",p

    pbiZ=dict()
    cdef int j,l,n,ai,bi
    for j from 1 <= j <= Q:
        pbiZ[j]=mpmath.mp.mpc(-Ypb[j],Xpb[j])
        #pbiZ[1-j]=pbiZ[j].conjugate()



    for n from Ms <= n <= Mf:
        #for ai from Dstart <= ai <= Dfinish:
        for a in D:
            #a=D[ai]
            nr=mpmath.mp.mpf(n+H.multiplier().Qv[a])
            nri=mpmath.mp.mpc(0,nr)
            #print "nri[",n,ai,"]=",nri
            nrtwopi=nr*twopi; nrtwopii=mpmath.mp.mpc(0,nrtwopi)
            for j from 1 <= j <= Q:
                fak2[n,j,a]=mpmath.mp.exp(-nri*Xm[j])
            if nr>=0:
                for j from 1 <= j <= Q:
                    tmp1=nr*pbiZ[j]
                    #if ai==Dstart and j==1:
                    #    print "i*n*Zpb[",n,ai,j,"]=",tmp1
                    fak1[n,j,a]=mpmath.exp(tmp1)
                    #if ai==Dstart and j==1:
                    #    print "i*n*Zpb[",n,ai,j,"]=",mpmath.exp(tmp1)
            else: 
                nrtwo=abs(nr*mptwo)
                for j from 1 <= j <= Q:
                    tmp1=mpmath.mp.exp(pbiZ[j]*nr)
                    fak1[n,j,a]=tmp1*mpmath.gammainc(kint,nrtwo*Ypb[j])
    #maxv=0
    #return fak1,fak2
    #for l from Ms <= l <= Mf:
    #    for ai from Dstart <= ai <= Dfinish:
    #        a=D[ai]
    #        for j from 1 <= j <= Qf:
    #            print "fak1[",l,ai,j,"]=",fak1[l,j,a]
    if verbose>0:
        for j from 1 <= j <= Q:
            #for ai from Dstart <= ai <= Dfinish:
            #    for bi from Dstart <= bi <= Dfinish:
            #        #if Cv[j].has_key((ai,bi)):
            #        #print "Cv[",j,"][",ai,bi,"]=",Cv[j][ai,bi]
            print "Cfak[",j,"][0]=",Cfak[j][0]
            print "Cfak[",1-j,"][1]=",Cfak[1-j][0]
        print "here0"
    for l from Ms <= l <= Mf:
        #for ai from Dstart <= ai <= Dfinish:
        for a in D:
            #a=D[ai]
            #for bi from Dstart <= bi <= Dfinish:
            for b in D:
                #b=D[bi]
                for j from 1 <= j <= Qf:
                    lj=Ml*(b-Dstart)+l-Ms
                    if b not in oddD: #(2*b) % orderD ==0: #is_int(2*b): #if (b == -b):
                        ch=Cv[j][a,b]
                        bii=(b*Cfak[1-j][1]) % abD
                        if Cv[j][a,bii]==0:
                            ch2=0
                        else:
                            ch2=1/Cv[j][a,bii]  #*Cfak[1-j][0]  #Cv[1-j][ai,bi]
                    else:
                        mbi= -b % abD 
                        ch=(Cv[j][a,b]+sym_type*Cv[j][a,mbi])
                        ch21=0; ch22=0
                        if Cv[j][a,b]<>0:
                            ch21=1/Cv[j][a,b]
                        if Cv[j][a,mbi]<>0:
                            ch22=1/Cv[j][a,mbi]
                        if sym_type==1:
                            ch2=(ch22+ch21)
                        else:
                            ch2=(ch22-ch21)
                    if ch<>0:
                        tmp1=fak1[l,j,b]*ch*Cfak[j][0]
                    else:
                        tmp1=0
                    if ch2<>0:
                        tmp2=fak1[l,j,b].conjugate()*ch2*Cfak[1-j][0] 
                    else:
                        tmp2=0
                    if ch==0 and ch2==0:
                        continue
                    nis=int(Ml*(a-Dstart)-Ms)
                    ## if lj==5 and j==1:
                    ##     print "ai,bi,mbi=",ai,bi,mbi
                    ##     print "sym=",sym_type
                    ##     print "Cv=",Cv[j][ai,bi]
                    ##     print "Cv[-bi]=",Cv[j][ai,mbi]
                    ##     print "ch=",ch
                    ##     print "ch2=",ch2
                    ##     print "tmp1=",tmp1
                    ##     print "tmp2=",tmp2 
                        
                    for n from Ms <= n <= Mf:
                        ni=nis+n
                        tmp=fak2[n,j,a]
                        tmpV=(tmp1*tmp+tmp2*tmp.conjugate())
                        #if lj==5 and n==1:
                        #    print "fak2[",n,ai,j,"]=",fak2[n,j,a]
                        #    print "tmpV1=",j,tmpV
                        V[ni,lj]=V[ni,lj]+tmpV
    #print "maxv=",maxv
    #return Cv
    #if verbose>0:
    #    print "V[9,9]=",V[9,9]
    #    print "V[30,30]=",V[30,30]

    twopiY=twopi*Y
    fourpiY=fourpi*Y
    for n from 0 <= n <= V.rows-1:
        for l from 0 <= l <= V.cols-1:
            V[n,l]=V[n,l]/Qfak
    maxb=0
    maxw=0
    for n from Ms <= n <= Mf:
        #for ai from Dstart <= ai <= Dfinish:
        for a in D:
            #a=D[ai]
            nr=mpmath.mp.mpf(n+H.multiplier().Qv[a])
            nri=mpmath.mp.mpc(0,nr)
            ni=Ml*(a-Dstart)+n-Ms
            nrY2pi=nr*twopiY
            nrY4pi=abs(nr*fourpiY)
            if nr>=0:
                kbes=mpmath.exp(-nrY2pi)
            else:
                kbes=mpmath.gammainc(kint,nrY4pi)*mpmath.exp(-nrY2pi)
            if abs(kbes)>maxb:
                maxb=abs(kbes)
            #if verbose>0:
            #    print "kbes(",a,n,")=",kbes
            V[ni,ni]=V[ni,ni]-kbes
            summa=0
            #ni=Ml*ai+n-Ms
            for (beta,m) in PP0:
                # We only need the first half
                if betai[beta] not in range(Dstart,Dfinish+1):
                    print "Skipping beta=",beta," by symmetry!"
                    print "index=",betai[beta]
                    print "range=",range(Dstart,Dfinish+1)
                    continue
                #print "beta=",beta
                #print "m=",m
                an=mpmath.mp.mpc(PP0[(beta,m)])
                if an.ae(mp0):
                    continue
                if beta == mbeta[beta]:
                    if sym_type==-1:
                        summa=0
                    else:
                        for j from 1 <= j <= Qf:
                            tmp=mmm[(beta,m)]*pbiZ[j]-nri*Xm[j]
                            tmpc=mpmath.exp(tmp)
                            ch=Cv[j][a,betai[beta]]*Cfak[j][0]
                            #ch2=Cv[1-j][ai,betai[beta]]
                            summa=summa+tmpc*ch
                            bii=(betai[beta]*Cfak[1-j][1]) % abD
                            if Cv[j][a,bii]<>0:
                                ch2=Cfak[1-j][0]/Cv[j][a,bii]
                                summa=summa+ch2*tmpc.conjugate()

                else:
                    for j from 1 <= j <= Qf:
                        ch=(Cv[j][a,betai[beta]]+sym_type*Cv[j][a,mbetai[beta]])
                        #ch2=Cv[1-j][ai,betai[beta]]+sym_type*Cv[1-j][ai,mbetai[beta]]
                        ch21=0; ch22=0
                        if Cv[j][a,betai[beta]]<>0:
                            ch21=1/Cv[j][a,betai[beta]]
                        if Cv[j][a,mbetai[beta]]<>0:
                            ch22=1/Cv[j][a,mbetai[beta]]
                        if sym_type==1:
                            ch2=(ch21+ch22)
                        else:
                            ch2=(ch22-ch21)
                        if ch==0 and ch2==0:
                            continue
                        tmp=mmm[(beta,m)]*pbiZ[j]-nri*Xm[j]
                        ch=ch*Cfak[j][0]
                        ch2=ch2*Cfak[1-j][0] 
                        tmpc=mpmath.exp(tmp)
                        summa=summa+tmpc*ch+ch2*tmpc.conjugate()
                #print "summa=",n,a,summa
                RHS[ni,0]=RHS[ni,0]+an*summa/Qfak
                if verbose>0:
                    print "summa=",n,a,summa
                    print "RHS[",ni,"]=",RHS[ni,0]
                if n==m:
                    if mbeta[beta]==beta:
                        if a==beta and sym_type==1:
                            RHS[ni,0]=RHS[ni,0]-an*mpmath.exp(-mmm[(beta,m)]*twopiY)
                    elif a==beta:
                        RHS[ni,0]=RHS[ni,0]-an*mpmath.exp(-mmm[(beta,m)]*twopiY)
                    elif a == mbeta[beta]:
                        RHS[ni,0]=RHS[ni,0]-an*sym_type*mpmath.exp(-mmm[(beta,m)]*twopiY)
                    if verbose>0:
                        print "subtracting for n==m!"
                        print "RHS[",ni,"]=",RHS[ni,0]
    #print "maxb=",maxb
    #print "maxw=",maxw
    W=dict()
    W['space']=H
    W['V']=V
    W['Ms']=Ms
    W['Mf']=Mf
    W['r']=len(D)
    W['nc']=nc
    W['RHS']=RHS
    W['sym_type']=sym_type
    return W


cpdef vv_harmonic_wmwf_setupV_mpc(H,dict PP,RealNumber Y,int M,int Q):
    r"""
    Set up the matrix for a vector-valued harmonic weak maass form for the Weil representation.
    Using MPC

    INPUT:
    - H -- instance of VVHarmonicWeakMaassforms
    
    EXAMPLES::

    sage: WR=WeilRepDiscriminantForm(1,False)
    sage: Y=mpmath.mp.mpf(0.5); k=mpmath.mp.mpf(0.5)
    sage: M=10; Q=20; m=-1; beta=1/2
    sage: W= vv_harmonic_wmwf_setupV(WR,beta,m,Y,M,Q,k)


    """
    cdef int prec
    prec = max(Y.base_ring().prec(),H._setupV_prec)
    if H._verbose>0:
        print "prec in setupV=",prec
    RF=RealField(prec)
    CF=MPComplexField(prec)
    CCF=ComplexField(prec)
    if prec==53:
        mpctx = mpmath.fp
    else:
        mpctx = mpmath.mp
    mpctx = mpmath.mp
    cdef RealNumber mpone,mptwo,twopi,fourpi,weight,mp0,Qfak 
    cdef RealNumber nrtwopi,nrtwo,tmp_r
    cdef MPComplexNumber tmp,ztmp,nrtwopii,twopii,summa
    weight = RF(H._weight)
    p0=RF(0); mpone=RF(1); mptwo=RF(2); mpfour=RF(4)
    twopi=mptwo*RF.pi();  fourpi=mptwo*twopi
    weight_over_two=weight/mptwo
    twopii=CF(0,twopi); summa=CF(0)
    ztmp=CF(1);  nrtwopii=CF(1); tmp=CF(0); tmp_r=RF(0)
    cdef int j,l,n,ai,bi,bii,mbi
    cdef int N,sym_type,verbose,len_of_pp,Qs,Qf,Ql,Ms,Mf,Ml,Qfaki
    cdef list D,PP0,beta_rat,pp_c,mbeta_rat #,mm_real,mmi_cplx
    cdef mpfr_t tmpr_t[1]
    mpfr_init2(tmpr_t[0],prec)
    D=H.index_set() #; N=H.rank()/2
    sym_type = H.sym_type(); verbose = H._verbose
    beta_rat=[]; pp_c=[]; mbeta_rat=[] #; mm_real=[]; mmi_cplx=[]
    cdef mpfr_t *mm_real
    cdef mpc_t *mmi_cplx
    if verbose>0:
        print "in vv_harmonic_wmwf_setupV_mpc!"
        print "args=",H,PP,Y,M,Q
        print "D=",D
    # PP is a dictionary containing the principal part in the form:
    # PP[(j,m)]= c or
    # PP[ D ]  = c or
    # corresponding to a term of the form c*e( Qv( D[j]+m )
    # or equivalently (in the second case) c*e( D/4N )
    # The list beta_int contains the list of all j's
    # The list beta_rat contains the list of all D[j]'s
    # The list pp_c contains the list of all c's
    ## Count how large the principal part is
    len_of_pp = len(PP.keys())
    cdef int  *mbeta_int, *beta_int
    beta_int = <int *>sage_malloc(sizeof(int)*len_of_pp)
    mbeta_int = <int *>sage_malloc(sizeof(int)*len_of_pp)    
    mm_real = <mpfr_t *>sage_malloc(sizeof(mpfr_t)*len_of_pp)
    mmi_cplx = <mpc_t *>sage_malloc(sizeof(mpc_t)*len_of_pp)        
    for i from 0 <= i <len_of_pp:
        mpfr_init2(mm_real[i],prec)
        mpc_init2(mmi_cplx[i],prec)
    for i from 0 <= i <len_of_pp:
        ttt = PP.keys()[i]
        if isinstance(ttt,tuple):
            (beta,m) = ttt
        else:
            (beta,m) = rn_from_D(H.multiplier(),ttt)
        pp_c.append(PP[ttt])
        if isinstance(beta,int) or isinstance(beta,Integer):
            beta_int[i]=int(beta)
            beta_rat.append(D[beta])
        elif isinstance(beta,Rational):
            beta_rat.append(beta)
            beta_int[i]=int(D.index(beta))            
    if verbose>0:
        print "beta_rat=",beta_rat
        print "pp_c    =",pp_c
    cdef Rational tt,mm
    tt = QQ(-1); mm=QQ(0)
    for i from 0 <= i < len_of_pp:
        p=numerator(1-beta_rat[i]); q=denominator(1-beta_rat[i])
        mb=QQ(p % q)/QQ(q) # reduce mod 1
        mbeta_rat.append(mb)
        mbeta_int[i]=int(D.index(mb))
        if verbose>0:
            print "b=",beta
            print "1-b=",mb
            print "mb=",mb
        if not mb in D:
            raise Exception,"Need beta=%s in D=" %(mb,D)
        mm = m+H.multiplier().Qv[beta_int[i]]
        if mm>0:
            raise ValueError,"Need principal part with negative exponents! Got:%s for PP:%s" % (mm,PP)
        tmp_r = RF(mm)
        mpfr_set(mm_real[i],tmp_r.value,rnd_re)
        #mm_real.append(RF(mm))
        #mmi_cplx.append(CF(0,mm))
        ztmp = CF(0,mm)
        mpc_set(mmi_cplx[i],ztmp.value,rnd)
        if verbose>0:
            print "mm=",mm
        if mm>tt:
            tt=mm
    if tt>=0 or H._holomorphic: # Holomorphic
        Ms=int(0); Mf=int(M); Ml=int(Mf-Ms+1)
    else:
        Ms=int(-M); Mf=int(M); Ml=int(Mf-Ms+1)
    cdef unsigned int abD
    abD=len(D)
    assert abD>0
    if verbose>0:
        print "Ms,Mf=",Ms,Mf
        print "beta=",beta_rat
        print "betai="
        for i from 0 <= i < len_of_pp:
            print ",bi[",i,"]=",beta_int[i]
            print "-beta[",i,"]=",mbeta_rat[i]
        print "mbetai="
        print "mm=",mm
        for i from 0 <= i < len_of_pp:
            print ",mbi[",i,"]=",mbeta_int[i]
            mpfr_set(tmp_r.value,mm_real[i],rnd_re)
            print "mmm_real=",tmp_r
            mpc_set(ztmp.value,mmi_cplx[i],rnd)
            print "mmi_cplx=",ztmp
    if Q==0:
        Q=M+20
    # We use symmetry here
    if verbose>0:
        print "Symmetry type=",sym_type
    if sym_type<>0:
        Qs=1; Qf=Q; Ql=Qf-Qs; Qfaki=2*Q
    else:
        Qs=1-Q; Qf=Q; Ql=Qf-Qs; Qfaki=4*Q
    kint=mpmath.mp.mpf(mpone-weight)
    ## Recall that we have different index sets depending on the symmetry
    cdef int Ds,Df,Dl,s
    Ds = int(H.index_set()[0])
    Df = int(H.index_set()[-1])
    Dl = len(H.index_set())
    # Pullpack points
    #cdef Vector_real_mpfr_dense Xm
    #cdef Vector_complex_dense Zpb
    cdef Vector_real_mpfr_dense Qvv
    cdef mpc_t ***Cv
    Cv = NULL
    Cv = <mpc_t***> sage_malloc ( sizeof(mpc_t**)*Ql)
    if Cv==NULL: raise MemoryError
    #cdef int nr = Df -Ds + 1
    cdef int nc = H.multiplier().ambient_rank()
    for j from 0 <= j < Ql: 
        Cv[j] = NULL
        Cv[j]=<mpc_t**>sage_malloc(sizeof(mpc_t*)*Dl)
        if Cv[j]==NULL: raise MemoryError
        for n from 0 <= n <  Dl:
            Cv[j][n] = NULL
            Cv[j][n]=<mpc_t*>sage_malloc(sizeof(mpc_t)*nc) 
            if Cv[j][n]==NULL: raise MemoryError
            for l from 0<= l < nc:
                mpc_init2(Cv[j][n][l],prec) #=<mpc_t>sage_malloc(sizeof(mpc_t)*nc) 
                #print "set Cv[",j,n,l,"]=",printf("%p ",Cv[j][n][l])
                #print "Set Cv",j,n,l
    cdef Vector_real_mpfr_dense Xm
    cdef Vector_complex_dense Zpb,Cfak_plus,Cfak_minus
    Xm = Vector_real_mpfr_dense(FreeModule(RF,Ql),0)
    Zpb= Vector_complex_dense(FreeModule(CF,Ql),0)
    Cfak_plus = Vector_complex_dense(FreeModule(CF,Ql),0)
    Cfak_minus = Vector_complex_dense(FreeModule(CF,Ql),0)
    pullback_pts_weil_rep_mpc2(H,Q,Y,weight,Xm,Zpb,Cfak_plus,Cfak_minus,Cv,Ds,Df)
    Dl=Df-Ds+1
    s=Dl*Ml
    if verbose>0:
        print "Ds, Df=",Ds,Df,Dl
        print "Ml=",Ml
        print "abD=",abD
        print "matrix size=",s
    MS=MatrixSpace(CF,s,s)
    cdef Matrix_complex_dense V,RHS 
    V=Matrix_complex_dense(MS,0)
    MS=MatrixSpace(CF,s,1)
    RHS=Matrix_complex_dense(MS,0)
    MS = MatrixSpace(CF,int(Df-Ds+1),int(Q))
    if verbose>0:
        print "k/2=",weight_over_two
        print "p=",p

    VS = vector(RF,Df-Ds+1).parent()
    Qvv=Vector_real_mpfr_dense(VS,0)
    nrtwopi=RF(0); nrtwo=RF(0)
    cdef mpc_t tmp1[1],tmp2[1]
    mpc_init2(tmp1[0],prec)
    mpc_init2(tmp2[0],prec)

    for j from Ds<= j <= Df:
        Qvv[j-Ds]=RF(H.multiplier().Qv[j])
    cdef mpfr_t nr[1]
    cdef mpc_t nri[1],tmpc_t[1]
    cdef Rational a
    a = QQ(0)
    mpfr_init2(nr[0],prec)
    mpc_init2(nri[0],prec)
    mpc_init2(tmpc_t[0],prec)
    cdef RealNumber tmparg
    tmparg=RF(1)
    cdef RealNumber nrY2pi,nrY4pi,kbes
    kbes=RF(0)
    nrY2pi=RF(0)
    nrY4pi=RF(0)
    cdef MPComplexNumber tmpc
    tmpc=CF(0)
    cdef mpc_t ***fak1,***fak2
    fak1=NULL; fak2=NULL
    fak1 =<mpc_t ***> sage_malloc( Ml*sizeof(mpc_t**))
    for n from 0 <= n < Ml:
        fak1[n]=<mpc_t **> sage_malloc( Dl*sizeof(mpc_t*))
        for ai from 0 <= ai < Dl:
            fak1[n][ai]=<mpc_t *>sage_malloc(Q*sizeof(mpc_t))
            for j from 0 <= j < Ql:
                mpc_init2(fak1[n][ai][j],prec)
    fak2 =<mpc_t ***> sage_malloc( Ml*sizeof(mpc_t**))
    for n from 0 <= n < Ml:
        fak2[n]=<mpc_t **> sage_malloc( Dl*sizeof(mpc_t*))
        for ai from 0 <= ai < Dl:
            fak2[n][ai]=<mpc_t *>sage_malloc(Q*sizeof(mpc_t))
            for j from 0 <= j < Q:
                mpc_init2(fak2[n][ai][j],prec)

    cdef tuple tu
    mparg = mpmath.mp.mpf(0)
    for n from 0 <= n < Ml:
        for ai from 0 <= ai < Dl:
            mpfr_set_si(nr[0],n+Ms,rnd_re)
            mpfr_add(nr[0],nr[0],Qvv._entries[ai],rnd_re)
            mpc_set_fr(nri[0],nr[0],rnd)
            mpfr_swap(mpc_realref(nri[0]),mpc_imagref(nri[0]))
            mpfr_set(tmp_r.value,nr[0],rnd_re)
            mpc_set(ztmp.value,nri[0],rnd)
            for j from 0 <= j < Ql:
                mpc_mul_fr(tmpc_t[0],nri[0],Xm._entries[j],rnd)
                mpc_neg(tmpc_t[0],tmpc_t[0],rnd)
                mpc_exp(fak2[n][ai][j],tmpc_t[0],rnd)
                #fak2[n,j,a]=CF(-nri[0]*Xm[j]).exp()
            if mpfr_sgn(nr[0])>0:
                for j from 0 <= j < Ql:
                    mpc_mul(tmpc_t[0],nri[0],Zpb._entries[j],rnd)
                    mpc_exp(tmpc_t[0],tmpc_t[0],rnd)
                    mpc_set(fak1[n][ai][j],tmpc_t[0],rnd)
            elif mpfr_sgn(nr[0])==0:
                for j from 0 <= j < Ql:
                    mpc_set_ui(fak1[n][ai][j],1,rnd)
            elif not H._holomorphic: 
                mpfr_mul_ui(nrtwo.value,nr[0],2,rnd_re)
                mpfr_abs(nrtwo.value,nrtwo.value,rnd_re)
                mpc_set_si(tmp1[0],1,rnd) 
                for j from 0 <= j < Ql:
                    mpc_mul(tmpc_t[0],nri[0],Zpb._entries[j],rnd)
                    mpc_exp(tmp1[0],tmpc_t[0],rnd)
                    mpfr_mul(tmparg.value,mpc_imagref(Zpb._entries[j]),nrtwo.value,rnd_re)
                    tu = mpfr_to_mpfval(tmparg.value)
                    mparg._set_mpf(tu)
                    testmp= mpctx.gammainc(kint,mparg)
                    mpfr_from_mpfval(kbes.value,testmp._mpf_)
                    #mpgamma = mpctx.gammainc(kint,tmparg)
                    mpc_mul_fr(tmp1[0],tmp1[0],kbes.value,rnd)
                    mpc_set(fak1[n][ai][j],tmp1[0],rnd)
                    #if n==Ms and j==1 and ai==1:
                    #    print "ERROR!"
            #else:
                
    if verbose>5:
        for j from 0 <= j < Ql:
            for ai from 0 <= ai < Dl:
                for bi from 0 <= bi < nc:
                    mpc_set(ztmp.value,Cv[j][ai][bi],rnd)
                    print "Cv[",j,"][",ai,bi,"]=",ztmp
                    #print "Cfak[",j,"][0]=",Cfak_plus[j-Qs]
                    #print "Cfak[",1-j,"][0]=",Cfak_minus[j-Qs]
                    #print "Cfak[",j,"][1]=",Cfak[j][1]
    cdef int t_ch_r,t_ch_i,t_ch2_r,t_ch2_i,ni,lj,nis
    cdef mpc_t ch[1],ch2[1],ch21[1],ch22[1],tmpV[1],tmpV2[1]
    mpc_init2(tmpV2[0],prec)
    mpc_init2(tmpV[0],prec)
    mpc_init2(ch[0],prec)
    mpc_init2(ch2[0],prec)
    mpc_init2(ch21[0],prec)
    mpc_init2(ch22[0],prec)
    cdef MPComplexNumber tmpV1
    cdef int* minus_bi
    minus_bi =<int*>sage_malloc(Df*sizeof(int))
    for i in range(Df):
        minus_bi[i]=int(H.multiplier().weil_module().neg_index(i))
    cdef int t
    tmpV1=CF(0)
    for l from 0 <= l < Ml:
        for ai from 0 <= ai < Dl:
            for bi from Ds <= bi <= Df:
                lj=Ml*(bi-Ds)+l
                if l==0 and Qvv[bi-Ds]<0 and H._holomorphic:
                    if verbose>2:
                        print "Skipping l=",l," b=",bi-Ds, " lj=",lj
                    continue
                for j from 0 <= j < Ql:
                    mbi = minus_bi[bi]
                    if mbi==bi:
                        mpc_set(ch[0],Cv[j][ai][bi],rnd)
                        # ch = Cv[j][ai,bi]
                        #mpc_set(ztmp.value,Cv[j][ai][bii],rnd)
                        if is_zero(Cv[j][ai][mbi]): #t_ch_r<>0 and t_ch_i<>0:
                            mpc_set_si(ch2[0],0,rnd) # ch2=CF(0)
                        else:
                            #print "non-zero!"
                            mpc_set_ui(ch2[0],1,rnd)
                            mpc_div(ch2[0],ch2[0],Cv[j][ai][mbi],rnd)
                            #ch2=1/Cv[j][ai,bii]  #*Cfak[1-j][0]  #Cv[1-j][ai,bi]
                    else:
                        #mbi= abD -bi # else -bi = 2N-bi
                        ##ch = (Cv[j][ai,bi]+sym_type*Cv[j][ai,mbi])
                        mpc_set(ch[0],Cv[j][ai][mbi],rnd)
                        mpc_mul_si(ch[0],ch[0],sym_type,rnd)
                        mpc_add(ch[0],ch[0],Cv[j][ai][bi],rnd)
                        mpc_set_si(ch21[0],0,rnd)
                        mpc_set_si(ch22[0],0,rnd)
                        #ch21=0; ch22=0
                        if is_zero(Cv[j][ai][bi])==0: # i.e. is non-zero
                            mpc_set_ui(ch21[0],1,rnd)
                            mpc_div(ch21[0],ch21[0],Cv[j][ai][bi],rnd)
                            #ch21=Cv[j][ai,bi]**-1
                        if is_zero(Cv[j][ai][mbi])==0:
                            mpc_set_ui(ch22[0],1,rnd)
                            mpc_div(ch22[0],ch22[0],Cv[j][ai][mbi],rnd)
                            #ztmp=Cv[j][ai,mbi]**-1
                        if sym_type==1:
                            mpc_add(ch2[0],ch22[0],ch21[0],rnd)
                        else:
                            mpc_sub(ch2[0],ch22[0],ch21[0],rnd)
                    if is_zero(ch[0])==0:
                        mpc_set(tmp1[0],fak1[l][bi-Ds][j],rnd)
                        mpc_mul(tmp1[0],tmp1[0],ch[0],rnd)
                        mpc_mul(tmp1[0],tmp1[0],Cfak_plus._entries[j],rnd)
                        #tmp1[0]=fak1[l][j-1,bi-Ds]*ch*Cfak[j][0]
                    else:
                        mpc_set_si(tmp1[0],0,rnd)
                    if is_zero(ch2[0])==0:
                        mpc_set(tmp2[0],fak1[l][bi-Ds][j],rnd)
                        mpc_conj(tmp2[0],tmp2[0],rnd)
                        mpc_mul(tmp2[0],tmp2[0],ch2[0],rnd)
                        mpc_mul(tmp2[0],tmp2[0],Cfak_minus._entries[j],rnd)
                    else:
                        mpc_set_si(tmp2[0],0,rnd)
                    if is_zero(ch[0]) and is_zero(ch2[0]):
                        continue
                    nis= Ml*ai
                    ## Do the case n=0 separately
                    if not H._holomorphic or Qvv[ai]>=0:
                        mpc_set(tmpV[0],fak2[0][ai][j],rnd)
                        mpc_mul(tmpV2[0],tmp1[0],tmpV[0],rnd)
                        mpc_conj(tmpV[0],tmpV[0],rnd)
                        mpc_mul(tmpV[0],tmpV[0],tmp2[0],rnd)                        
                        mpc_add(tmpV1.value,tmpV[0],tmpV2[0],rnd)
                        mpc_add(V._matrix[nis][lj],V._matrix[nis][lj],tmpV1.value,rnd)                        
                    for n from 1 <= n < Ml:
                        ni=nis+n
                        mpc_set(tmpV[0],fak2[n][ai][j],rnd)
                        mpc_mul(tmpV2[0],tmp1[0],tmpV[0],rnd)
                        mpc_conj(tmpV[0],tmpV[0],rnd)
                        mpc_mul(tmpV[0],tmpV[0],tmp2[0],rnd)                        
                        mpc_add(tmpV1.value,tmpV[0],tmpV2[0],rnd)
                        #tmpV1=tmp1[0]*ztmp+tmp2*ztmp.conjugate()
                        #tmpV1=(tmp1[0]*ztmp+tmp2*ztmp.conjugate())
                        mpc_add(V._matrix[ni][lj],V._matrix[ni][lj],tmpV1.value,rnd)
    if verbose>0:
        print "here1"
    #print "maxv=",maxv
    #return Cv
    if verbose>0:
        print "V[9,9]=",V[9,9]
        print "V[30,30]=",V[30,30]
    cdef mpfr_t twopiY[1]
    mpfr_init2(twopiY[0],prec)
    mpfr_mul(twopiY[0],twopi.value,Y.value,rnd_re)
    ni = V.nrows()-1
    lj = V.ncols()-1
    for n from 0 <= n <= ni:
        for l from 0 <= l <= lj:
            mpc_div_ui(V._matrix[n][l],V._matrix[n][l],Qfaki,rnd)
    maxb=0
    maxw=0
    cdef MPComplexNumber an
    an = CF(0)

    for ai from 0 <= ai < Dl:
        for n from 0 <= n < Ml:
            #a=D[ai+Ds]
            mpfr_set(nr[0],Qvv._entries[ai],rnd_re)
            mpfr_add_si(nr[0],nr[0],n+Ms,rnd_re)
            mpc_set_fr(nri[0],nr[0],rnd)
            mpfr_swap(mpc_realref(nri[0]),mpc_imagref(nri[0]))
            ni=Ml*ai+n
            mpfr_mul(nrY2pi.value,nr[0],twopiY[0],rnd_re)
            if verbose>1:
                print "ni=",ni
                mpfr_set(tmp_r.value,nr[0],rnd_re)
                print "nr=",tmp_r
            if mpfr_sgn(nr[0])>=0:
                kbes=(-nrY2pi).exp()
            elif not H._holomorphic:
                mpfr_mul_si(nrY4pi.value,nrY2pi.value,-2,rnd_re)
                #if verbose>0:
                #    print "nrY4pi=",nrY4pi
                tu = mpfr_to_mpfval(nrY4pi.value)
                mparg._set_mpf(tu)
                testmp=mpctx.gammainc(kint,nrY4pi)
                mpfr_from_mpfval(kbes.value,testmp._mpf_)
                kbes=kbes*(-nrY2pi).exp()
            else:
                #continue
                if verbose>1:
                    print "kbes=0 (maybe skipping?)"
                kbes = RF(0) #continue
            mpc_sub_fr(V._matrix[ni][ni],V._matrix[ni][ni],kbes.value,rnd)
            mpc_set_si(summa.value,0,rnd) #=CF(0)
            #ni=Ml*ai+n-Ms
            for bi in range(len_of_pp):
                #if verbose>0 :
                #    print "bi=",bi
                # We only need the first half
                if beta_int[bi]<Ds or beta_int[bi]>Df:
                    if verbose>1:
                        print "Skipping beta=",beta," by symmetry!"
                        print "index=",beta_int[bi]
                        print "range=",range(Ds,Df+1)
                    continue
                an=CF(pp_c[bi])
                if an==0:
                    if verbose>1:
                        print "an=0, skipping!"
                    continue
                else:
                    if verbose>1:
                        print "an=",an
                mpc_set_si(summa.value,0,rnd) #=CF(0)
                if beta_int[bi] == mbeta_int[bi]:
                    if sym_type==-1:
                        if verbose>1:
                            print "cancel due to symmetry!"
                        
                        #summa=CF(0)
                    else:
                        for j from 0 <= j < Ql:
                            mpc_mul_fr(tmp1[0],nri[0],Xm._entries[j],rnd)
                            #mpc_mul
                            mpc_set(tmp.value,tmp1[0],rnd)
                            #tmp=mmi_cplx[bi]*Zpb[j]-tmp  
                            mpc_mul(tmp1[0],mmi_cplx[bi],Zpb._entries[j],rnd)
                            mpc_sub(tmp.value,tmp1[0],tmp.value,rnd)
                            #tmp=mmmi[(beta,m)]*Zpb[j]-nri[0]*Xm[j]
                            mpc_exp(tmpc.value,tmp.value,rnd) #=tmp.exp()
                            mpc_set(tmp1[0],Cv[j][ai][beta_int[bi]],rnd)
                            #tmp1[0]=tmp1[0]*Cfak[j][0] TEST
                            mpc_mul(tmp1[0],tmp1[0],Cfak_plus._entries[j],rnd)
                            #ch2=Cv[1-j][ai,betai[beta]]
                            mpc_mul(tmp1[0],tmp1[0],tmpc.value,rnd)
                            mpc_add(summa.value,summa.value,tmp1[0],rnd)
                            #rtmp1[0].summa=summa+tmpc*tmp1[0]
                            #if beta_int[bi] == 0 or beta_int[bi]==N:
                            mbi = minus_bi[bi]
                            #bii = beta_int[bi]
                            #else:
                            #    bii= abD - beta_int[bi] #
                            if is_zero(Cv[j][ai][mbi])==0:
                                #t_ch_r==0 or t_ch_i==0:
                                # if Cv[j][ai,bii]<>0:
                                mpc_div(tmp1[0],Cfak_minus._entries[j],Cv[j][ai][mbi],rnd)
                                mpc_set(tmp2[0],tmpc.value,rnd)
                                mpc_conj(tmp2[0],tmp2[0],rnd)
                                mpc_mul(tmp1[0],tmp1[0],tmp2[0],rnd)

                                mpc_add(summa.value,summa.value,tmp1[0],rnd)
                                #summa=summa+tmp1[0]*tmpc.conjugate()

                else:
                    for j from 0 <= j < Ql:
                        #ch=(Cv[j][ai,beta_int[bi]]+sym_type*Cv[j][ai,mbeta_int[bi]])
                        mpc_set(ch[0],Cv[j][ai][mbeta_int[bi]],rnd)
                        mpc_mul_si(ch[0],ch[0],sym_type,rnd)
                        mpc_add(ch[0],ch[0],Cv[j][ai][beta_int[bi]],rnd)
                        if is_zero(Cv[j][ai][beta_int[bi]])==0:
                            mpc_set_si(ch21[0],1,rnd)
                            #ch21=1/Cv[j][ai,beta_int[bi]]
                            mpc_div(ch21[0],ch21[0],Cv[j][ai][beta_int[bi]],rnd)
                        if is_zero(Cv[j][ai][mbeta_int[bi]])==0:
                            mpc_set_si(ch22[0],1,rnd)
                            mpc_div(ch22[0],ch22[0],Cv[j][ai][mbeta_int[bi]],rnd)
                            #ch22=1/Cv[j][ai,mbeta_int[bi]]
                        if sym_type==1:
                            mpc_add(ch2[0],ch21[0],ch22[0],rnd)
                            #ch2=(ch21+ch22)
                        else:
                            mpc_sub(ch2[0],ch22[0],ch21[0],rnd)
                            #ch2=(ch22-ch21)
                        if is_zero(ch[0]) and is_zero(ch2[0]):
                            continue
                        mpc_mul_fr(tmp1[0],nri[0],Xm._entries[j],rnd)
                        #tmp=mmmi[(beta,m)]*Zpb[j]-nri[0]*Xm[j]
                        mpc_set((<MPComplexNumber>tmp).value,tmp1[0],rnd)
                        
                        #tmp=mmi_cplx[bi]*Zpb[j]-tmp # -tmp1[0]
                        #tmpc=tmp.exp()
                        mpc_mul(tmp1[0],mmi_cplx[bi],Zpb._entries[j],rnd)
                        mpc_sub(tmp.value,tmp1[0],tmp.value,rnd)
                        mpc_exp(tmpc.value,tmp.value,rnd)
                        
                        mpc_mul(ch[0],ch[0],Cfak_plus._entries[j],rnd)
                        mpc_mul(ch2[0],ch2[0],Cfak_minus._entries[j],rnd)
                        mpc_mul(ch[0],tmpc.value,ch[0],rnd)
                        mpc_conj(tmpc.value,tmpc.value,rnd)
                        mpc_mul(ch2[0],tmpc.value,ch2[0],rnd)
                        mpc_add(ch[0],ch[0],ch2[0],rnd)
                        mpc_add(summa.value,summa.value,ch[0],rnd)
                        #summa=summa+tmpc*ch+ch2*tmpc.conjugate()

                mpc_mul(summa.value,summa.value,an.value,rnd)
                mpc_div_ui(summa.value,summa.value,Qfaki,rnd)
                mpc_add(RHS._matrix[ni][0],RHS._matrix[ni][0],summa.value,rnd)
                #RHS[ni,0]=RHS[ni,0]+an*summa/Qfak
                if verbose>2:
                    print "summa=",n,ai,bi,summa
                    print "RHS[",ni,"]=",RHS[ni,0]
                    #print "ai=",ai
                    print "beta_int[bi]=",beta_int[bi]
                    #print "mbeta_int[bi]=",mbeta_int[bi]
                if n+Ms == m:
                    if mbeta_int[bi]==beta_int[bi]:
                        if ai+Ds==beta_int[bi] and sym_type==1:
                            mpfr_mul(tmpr_t[0],twopiY[0],mm_real[bi],rnd_re)
                            mpfr_neg(tmpr_t[0],tmpr_t[0],rnd_re)
                            mpfr_exp(tmpr_t[0],tmpr_t[0],rnd_re)
                            mpc_mul_fr(an.value,an.value,tmpr_t[0],rnd)
                            #an = an*(-mm_real[bi]*twopiY).exp()
                            RHS[ni,0]=RHS[ni,0]- an
                    elif ai+Ds==beta_int[bi]:
                        mpfr_mul(tmpr_t[0],twopiY[0],mm_real[bi],rnd_re)
                        mpfr_neg(tmpr_t[0],tmpr_t[0],rnd_re)
                        mpfr_exp(tmpr_t[0],tmpr_t[0],rnd_re)
                        mpc_mul_fr(an.value,an.value,tmpr_t[0],rnd)
                        #an = an*(-mm_real[bi]*twopiY).exp()
                        RHS[ni,0]=RHS[ni,0]-an #*(-mm_real[bi]*twopiY).exp()
                    elif ai+Ds == mbeta_int[bi]:
                        mpfr_mul(tmpr_t[0],twopiY[0],mm_real[bi],rnd_re)
                        mpfr_neg(tmpr_t[0],tmpr_t[0],rnd_re)
                        mpfr_exp(tmpr_t[0],tmpr_t[0],rnd_re)
                        mpc_mul_fr(an.value,an.value,tmpr_t[0],rnd)
                        mpc_mul_si(an.value,an.value,sym_type,rnd)
                        #an = sym_type*an*(-mm_real[bi]*twopiY).exp()
                        RHS[ni,0]=RHS[ni,0]-an  #*(-mm_real[bi]*twopiY).exp()


    # Clear variables
    mpc_clear(ch[0])
    mpc_clear(tmpV[0])
    mpc_clear(tmp1[0])
    mpc_clear(tmp2[0])
    mpc_clear(tmpV2[0])
    mpc_clear(ch2[0])
    mpc_clear(ch21[0])
    mpc_clear(ch22[0])
    mpfr_clear(nr[0])
    mpfr_clear(tmpr_t[0])
    mpfr_clear(twopiY[0])
    mpc_clear(nri[0])
    mpc_clear(tmpc_t[0])
    sage_free(beta_int)
    sage_free(mbeta_int)
    sage_free(mm_real)
    sage_free(mmi_cplx)
    if Cv<>NULL:
        for j from 0 <= j < Ql: 
            if Cv[j]<>NULL:
                for n from 0 <= n <  Dl:
                    if Cv[j][n]<>NULL:
                        for l from 0<= l < nc:
                            #print "Cv[",j,n,l,"]=",printf("%p ",Cv[j][n][l])
                            #mpc_set(ztmp.value,Cv[j][n][l],rnd)
                            #print ztmp
                            mpc_clear(Cv[j][n][l])
                        sage_free(Cv[j][n])
                sage_free(Cv[j])
        sage_free(Cv)

    if fak1<>NULL:
        for n from 0<=n < Ml:
            if fak1[n]<>NULL:
                for ai from 0<=ai < Dl:
                    if fak1[n][ai]<>NULL:
                        for j from 0<=j<=Q-1:
                            mpc_clear(fak1[n][ai][j])
                        sage_free(fak1[n][ai])
                sage_free(fak1[n])            
        sage_free(fak1)
    if fak2<>NULL:
        for n from 0<=n < Ml:
            if fak2[n]<>NULL:
                for ai from 0<=ai < Dl:
                    if fak2[n][ai]<>NULL:
                        for j from 0<=j<=Q-1:
                            mpc_clear(fak2[n][ai][j])
                        sage_free(fak2[n][ai])
                sage_free(fak2[n])            
        sage_free(fak2)
    W=dict()
    W['space']=H
    W['V']=V
    W['Ms']=Ms
    W['Mf']=Mf
    W['r']=Dl
    W['nc']=nc
    W['RHS']=RHS
    W['sym_type']=sym_type
    return W

cpdef vv_harmonic_wmwf_setupV_mpc2(H,dict PP,RealNumber Y,int M,int Q):
    r"""
    Set up the matrix for a vector-valued harmonic weak maass form for the Weil representation.
    Using MPC

    INPUT:
    - H -- instance of VVHarmonicWeakMaassforms
        if has
    EXAMPLES::

    sage: WR=WeilRepDiscriminantForm(1,False)
    sage: Y=mpmath.mp.mpf(0.5); k=mpmath.mp.mpf(0.5)
    sage: M=10; Q=20; m=-1; beta=1/2
    sage: W= vv_harmonic_wmwf_setupV(WR,beta,m,Y,M,Q,k)


    """
    ### Only implemented for Weil representation
    X = H.multiplier()
    WM = X.weil_module()
    cdef list D = X.D()
    
    if WM == None:
        if H._verbose>0:
            print "WM=",WM
        raise ValueError,"Need a space with a Weil representation! Got:{0}".format(H)
    
    cdef int prec
    if hasattr(H,"_setupV_prec"):
        prec = int(max(Y.base_ring().prec(),H._setupV_prec))
    else:
        prec = int(Y.base_ring().prec())
    if H._verbose>0:
        print "prec in setupV=",prec
    RF=RealField(prec)
    CF=MPComplexField(prec)
    CCF=ComplexField(prec)
    if H._verbose>0:
        print "here!"
    if prec==53:
        mpctx = mpmath.fp
    else:
        mpctx = mpmath.mp
    mpctx = mpmath.mp
    cdef RealNumber mpone,mptwo,twopi,fourpi,weight,mp0,Qfak 
    cdef RealNumber nrtwopi,nrtwo,tmp_r
    cdef MPComplexNumber tmp,ztmp,nrtwopii,twopii,summa
    weight = RF(H._weight)
    p0=RF(0); mpone=RF(1); mptwo=RF(2); mpfour=RF(4)
    twopi=mptwo*RF.pi();  fourpi=mptwo*twopi
    weight_over_two=weight/mptwo
    twopii=CF(0,twopi); summa=CF(0)
    ztmp=CF(1);  nrtwopii=CF(1); tmp=CF(0); tmp_r=RF(0)
    cdef int j,l,n,ai,bi,bii,mbi
    cdef int N,sym_type,verbose,len_of_pp,Qs,Qf,Ql,Ms,Mf,Ml,Qfaki
    cdef list PP0,beta_rat,pp_c,mbeta_rat #,mm_real,mmi_cplx
    cdef mpfr_t tmpr_t[1]
    mpfr_init2(tmpr_t[0],prec)
    #D=WM.basis(); N=WM.rank()/2
    #D = 
    sym_type = H.sym_type(); verbose = H._verbose
    beta_rat=[]; pp_c=[]; mbeta_rat=[] #; mm_real=[]; mmi_cplx=[]
    cdef mpfr_t *mm_real
    cdef mpc_t *mmi_cplx
    if verbose>0:
        print "in vv_harmonic_wmwf_setupV_mpc2!"
        print "args=",H,PP,Y,M,Q
        print "D=",D
    # PP is a dictionary containing the principal part in the form:
    # PP[(j,m)]= c or
    # PP[ D ]  = c or
    # corresponding to a term of the form c*e( Qv( D[j]+m )
    # or equivalently (in the second case) c*e( D/4N )
    # The list beta_int contains the list of all j's
    # The list beta_rat contains the list of all D[j]'s
    # The list pp_c contains the list of all c's
    ## Count how large the principal part is
    if PP.has_key('+'):
        len_of_pp = len(PP['+'].keys())
    else:
        len_of_pp = 0
    cdef int  *mbeta_int, *beta_int
    beta_int = <int *>sage_malloc(sizeof(int)*len_of_pp)
    mbeta_int = <int *>sage_malloc(sizeof(int)*len_of_pp)    
    mm_real = <mpfr_t *>sage_malloc(sizeof(mpfr_t)*len_of_pp)
    mmi_cplx = <mpc_t *>sage_malloc(sizeof(mpc_t)*len_of_pp)        
    for i from 0 <= i <len_of_pp:
        mpfr_init2(mm_real[i],prec)
        mpc_init2(mmi_cplx[i],prec)
    for i from 0 <= i <len_of_pp:
        ttt = PP['+'].keys()[i]
        if isinstance(ttt,tuple):
            (beta,m) = ttt
        elif isinstance(ttt,(int,Integer)):
            (beta,m) = rn_from_D(H.multiplier(),ttt)
        else:
            raise ValueError,"principal part is not in corect form! Got:{0}".format(PP)
        pp_c.append(PP['+'][ttt])
        if isinstance(beta,int) or isinstance(beta,Integer):
            beta_int[i]=int(beta)
        #elif isinstance(beta,Rational):
        #    beta_rat.append(beta)
        #    beta_int[i]=int(D.index(beta))            
    if verbose>0:
        print "pp_c    =",pp_c
    cdef Rational tt,mm
    tt = QQ(-1); mm=QQ(0)
    #cdef int mbi
    for i from 0 <= i < len_of_pp:
        mbeta_int[i] = WM._neg_index(beta_int[i])
        mm = m+H.multiplier().Qv[beta_int[i]]
        if mm>0:
            raise ValueError,"Need principal part with non-positive exponents! Got:%s for PP:%s" % (mm,PP)
        tmp_r = RF(mm)
        mpfr_set(mm_real[i],tmp_r.value,rnd_re)
        #mm_real.append(RF(mm))
        #mmi_cplx.append(CF(0,mm))
        ztmp = CF(0,mm)
        mpc_set(mmi_cplx[i],ztmp.value,rnd)
        if verbose>0:
            print "mm=",mm
        if mm>tt:
            tt=mm
    if tt>=0 or H._holomorphic: # Holomorphic
        Ms=int(0); Mf=int(M); Ml=int(Mf-Ms+1)
    else:
        Ms=int(-M); Mf=int(M); Ml=int(Mf-Ms+1)
    cdef unsigned int abD
    #abD=len(D)
    abD = WM.rank()
    assert abD>0
    if verbose>0:
        print "Ms,Mf=",Ms,Mf
        print "betai="
        for i from 0 <= i < len_of_pp:
            print ",bi[",i,"]=",beta_int[i]
        print "mbetai="
        print "mm=",mm
        for i from 0 <= i < len_of_pp:
            print ",mbi[",i,"]=",mbeta_int[i]
            mpfr_set(tmp_r.value,mm_real[i],rnd_re)
            print "mmm_real=",tmp_r
            mpc_set(ztmp.value,mmi_cplx[i],rnd)
            print "mmi_cplx=",ztmp
    if Q==0:
        Q=M+20
    # We use symmetry here
    if verbose>0:
        print "Symmetry type=",sym_type
        print "D=",D
    if sym_type<>0:
        Qs=1; Qf=Q; Ql=Qf-Qs+1; Qfaki=2*Q
    else:
        Qs=1-Q; Qf=Q; Ql=Qf-Qs+1; Qfaki=4*Q
    cdef int kinti
    cdef int is_int=0
    cdef int is_half_int=0
    cdef RealNumber kint
    kint = mpone - weight
    if floor(kint)==pceil(kint):
        kinti = int(kint); is_int = 1
    if is_int==0:
        ## Check if kint is half-integral.
        if floor(2*kint)==pceil(2*kint):
            is_half_int = 1
            kinti = int(kint-RF(0.5))   
 
    ## Recall that we have different index sets depending on the symmetry
    cdef int Ds,Df,Dl,s
    if sym_type==1:
        Ds =D[0]; Df=D[-1] #len(D)-1  #int(0); Df=int(N) # 0,1,...,N
    elif sym_type==-1:
        Ds=D[0]; Df=D[-1] #len(D)-1 #int(N-1) # 1,2,...,N-1  (since -0=0 and -N=N)
    if verbose>0:
        print "Ds=",Ds,"Df=",D
    Dl = Df - Ds + 1
    assert Dl==len(D)
    # Pullpack points
    #cdef Vector_real_mpfr_dense Xm
    #cdef Vector_complex_dense Zpb
    cdef Vector_real_mpfr_dense Qvv
    cdef mpc_t ***Cv
    Cv = NULL
    Cv = <mpc_t***> sage_malloc ( sizeof(mpc_t**)*Ql)
    if Cv==NULL: raise MemoryError
    #cdef int nr = Df -Ds + 1
    #cdef int nc = len(WM.basis())
    ## Not used. assumed to be 1
    cdef int nc = H.group().ncusps() #len(WM.basis())
    cdef int rank = X.rank()
    cdef int full_rank = len(X.Qv)
    if verbose>0:
        print "precBBBB:",prec
    for j from 0 <= j < Ql: 
        Cv[j] = NULL
        Cv[j]=<mpc_t**>sage_malloc(sizeof(mpc_t*)*Dl)
        if Cv[j]==NULL: raise MemoryError
        for n from 0 <= n <  Dl:
            Cv[j][n] = NULL
            Cv[j][n]=<mpc_t*>sage_malloc(sizeof(mpc_t)*full_rank) 
            if Cv[j][n]==NULL: raise MemoryError
            for l from 0<= l < full_rank:
                #if verbose>0:
                #    print "j,n,l=",j,n,l
                mpc_init2(Cv[j][n][l],prec) #=<mpc_t>sage_malloc(sizeof(mpc_t)*nc) 
                #print "set Cv[",j,n,l,"]=",printf("%p ",Cv[j][n][l])
                #print "Set Cv",j,n,l
    cdef Vector_real_mpfr_dense Xm
    cdef Vector_complex_dense Zpb,Cfak_plus,Cfak_minus
    Xm = Vector_real_mpfr_dense(FreeModule(RF,Ql),0)
    Zpb= Vector_complex_dense(FreeModule(CF,Ql),0)
    Cfak_plus = Vector_complex_dense(FreeModule(CF,Ql),0)
    Cfak_minus = Vector_complex_dense(FreeModule(CF,Ql),0)
    if verbose>0:
        print "Dl=",Dl
        print "rank=",rank
    pullback_pts_weil_rep_mpc2(H,Q,Y,weight,Xm,Zpb,Cfak_plus,Cfak_minus,Cv,Ds,Df)
    Dl=Df-Ds+1
    s=Dl*Ml
    if verbose>0:
        print "Ds, Df=",Ds,Df,Dl
        print "Ml=",Ml
        print "abD=",abD
        print "matrix size=",s
    MS=MatrixSpace(CF,s,s)
    cdef Matrix_complex_dense V,RHS 
    V=Matrix_complex_dense(MS,0)
    MS=MatrixSpace(CF,s,1)
    RHS=Matrix_complex_dense(MS,0)
    MS = MatrixSpace(CF,int(Df-Ds+1),int(Q))
    if verbose>0:
        print "k/2=",weight_over_two
        #print "p=",p

    VS = vector(RF,Df-Ds+1).parent()
    Qvv=Vector_real_mpfr_dense(VS,0)
    nrtwopi=RF(0); nrtwo=RF(0)
    cdef mpc_t tmp1[1],tmp2[1]
    mpc_init2(tmp1[0],prec)
    mpc_init2(tmp2[0],prec)

    for j from Ds<= j <= Df:
        Qvv[j-Ds]=RF(H.multiplier().Qv[j])
    cdef mpfr_t nr[1]
    cdef mpc_t nri[1],tmpc_t[1]
    cdef Rational a
    a = QQ(0)
    mpfr_init2(nr[0],prec)
    mpc_init2(nri[0],prec)
    mpc_init2(tmpc_t[0],prec)
    cdef RealNumber tmparg
    tmparg=RF(1)
    cdef RealNumber nrY2pi,nrY4pi,kbes
    kbes=RF(0)
    nrY2pi=RF(0)
    nrY4pi=RF(0)
    cdef MPComplexNumber tmpc
    tmpc=CF(0)
    cdef mpc_t ***fak1,***fak2
    fak1=NULL; fak2=NULL
    fak1 =<mpc_t ***> sage_malloc( Ml*sizeof(mpc_t**))
    for n from 0 <= n < Ml:
        fak1[n]=<mpc_t **> sage_malloc( Dl*sizeof(mpc_t*))
        for ai from 0 <= ai < Dl:
            fak1[n][ai]=<mpc_t *>sage_malloc(Q*sizeof(mpc_t))
            for j from 0 <= j < Ql:
                mpc_init2(fak1[n][ai][j],prec)
    fak2 =<mpc_t ***> sage_malloc( Ml*sizeof(mpc_t**))
    for n from 0 <= n < Ml:
        fak2[n]=<mpc_t **> sage_malloc( Dl*sizeof(mpc_t*))
        for ai from 0 <= ai < Dl:
            fak2[n][ai]=<mpc_t *>sage_malloc(Q*sizeof(mpc_t))
            for j from 0 <= j < Q:
                mpc_init2(fak2[n][ai][j],prec)

    cdef tuple tu
    cdef int ok=0
    mparg = mpmath.mp.mpf(0)
    for n from 0 <= n < Ml:
        for ai from 0 <= ai < Dl:
            mpfr_set_si(nr[0],n+Ms,rnd_re)
            mpfr_add(nr[0],nr[0],Qvv._entries[ai],rnd_re)
            mpc_set_fr(nri[0],nr[0],rnd)
            mpfr_swap(mpc_realref(nri[0]),mpc_imagref(nri[0]))
            mpfr_set(tmp_r.value,nr[0],rnd_re)
            mpc_set(ztmp.value,nri[0],rnd)
            for j from 0 <= j < Ql:
                mpc_mul_fr(tmpc_t[0],nri[0],Xm._entries[j],rnd)
                mpc_neg(tmpc_t[0],tmpc_t[0],rnd)
                mpc_exp(fak2[n][ai][j],tmpc_t[0],rnd)
                #fak2[n,j,a]=CF(-nri[0]*Xm[j]).exp()
            if mpfr_sgn(nr[0])>0:
                for j from 0 <= j < Ql:
                    mpc_mul(tmpc_t[0],nri[0],Zpb._entries[j],rnd)
                    mpc_exp(tmpc_t[0],tmpc_t[0],rnd)
                    mpc_set(fak1[n][ai][j],tmpc_t[0],rnd)
            elif mpfr_sgn(nr[0])==0:
                for j from 0 <= j < Ql:
                    mpc_set_ui(fak1[n][ai][j],1,rnd)
            elif not H._holomorphic: 
                mpfr_mul_ui(nrtwo.value,nr[0],2,rnd_re)
                mpfr_abs(nrtwo.value,nrtwo.value,rnd_re)
                mpc_set_si(tmp1[0],1,rnd) 
                for j from 0 <= j < Ql:
                    mpc_mul(tmpc_t[0],nri[0],Zpb._entries[j],rnd)
                    mpc_exp(tmp1[0],tmpc_t[0],rnd)
                    mpfr_mul(tmparg.value,mpc_imagref(Zpb._entries[j]),nrtwo.value,rnd_re)
                    ok = 0
                    if is_int==1:
                        if kinti>0:
                            ok = incgamma_pint_c(kbes.value,kinti,tmparg.value,verbose)
                        else:
                            ok = incgamma_nint_c(kbes.value,kinti,tmparg.value,verbose)
                    elif is_half_int==1:
                        ok = incgamma_hint_c(kbes.value,kinti,tmparg.value,verbose)
                    if not ok:
                        tu = mpfr_to_mpfval(tmparg.value)
                        mparg._set_mpf(tu)
                        testmp=mpctx.gammainc(kint,tmparg)
                        mpfr_from_mpfval(kbes.value,testmp._mpf_)
                    #tu = mpfr_to_mpfval(tmparg.value)
                    #mparg._set_mpf(tu)                    
                    #testmp= mpctx.gammainc(kint,mparg)
                    #mpfr_from_mpfval(kbes.value,testmp._mpf_)
                    #mpgamma = mpctx.gammainc(kint,tmparg)
                    mpc_mul_fr(tmp1[0],tmp1[0],kbes.value,rnd)
                    mpc_set(fak1[n][ai][j],tmp1[0],rnd)
                    #if n==Ms and j==1 and ai==1:
                    #    print "ERROR!"
            #else:
                
    if verbose>5:
        for j from 0 <= j < Ql:
            for ai from 0 <= ai < Dl:
                for bi from 0 <= bi < rank:
                    mpc_set(ztmp.value,Cv[j][ai][bi],rnd)
                    print "Cv[",j,"][",ai,bi,"]=",ztmp
                    #print "Cfak[",j,"][0]=",Cfak_plus[j-Qs]
                    #print "Cfak[",1-j,"][0]=",Cfak_minus[j-Qs]
                    #print "Cfak[",j,"][1]=",Cfak[j][1]
    cdef int t_ch_r,t_ch_i,t_ch2_r,t_ch2_i,ni,lj,nis
    cdef mpc_t ch[1],ch2,ch21,ch22,tmpV,tmpV2
    if verbose>0:
        print "precAAA:",prec
    mpc_init2(tmpV2,prec)
    mpc_init2(tmpV,prec)
    mpc_init2(ch[0],prec)
    mpc_init2(ch2,prec)
    mpc_init2(ch21,prec)
    mpc_init2(ch22,prec)
    cdef MPComplexNumber tmpV1
    cdef int order = WM.finite_quadratic_module().order()
    cdef int* minus_ix
    minus_ix = <int*>sage_malloc(sizeof(int)*order)
    #QM = WM.finite_quadratic_module()
    for i in range(order):
        minus_ix[i] = WM._neg_index(i)
        
        
    cdef int t
    tmpV1=CF(0)
    for l from 0 <= l < Ml:
        #for ai from Ds <= ai <= Df:
        for ai from 0 <= ai < Dl:
            for bi from Ds <= bi <= Df:
            #for bi from 0 <= bi < Dl:
                lj=Ml*(bi-Ds)+l
                if l==0 and Qvv[bi-Ds]<0 and H._holomorphic:
                    if verbose>2:
                        print "Skipping l=",l," b=",bi-Ds, " lj=",lj
                    continue
                for j from 0 <= j < Ql:
                    if minus_ix[bi]==bi: # bi==0 or bi==N:
                        mpc_set(ch[0],Cv[j][ai][bi],rnd)
                        # ch = Cv[j][ai,bi]
                        bii= bi # this is true for bi=0 and bi=N
                        #mpc_set(ztmp.value,Cv[j][ai][bii],rnd)
                        if is_zero(Cv[j][ai][bii]): #t_ch_r<>0 and t_ch_i<>0:
                            mpc_set_si(ch2,0,rnd) # ch2=CF(0)
                        else:
                            #print "non-zero!"
                            mpc_set_ui(ch2,1,rnd)
                            #if verbose>0:
                            #    print "j,ai,bii=",j,ai,bii
                            #    mpc_set(ztmp.value,Cv[j][ai][bii],rnd)
                            #    print "Cv[",j,"][",ai,bii,"]=",ztmp
                            mpc_div(ch2,ch2,Cv[j][ai][bii],rnd)
                            #ch2=1/Cv[j][ai,bii]  #*Cfak[1-j][0]  #Cv[1-j][ai,bi]
                    else:
                        #mbi= abD -bi # else -bi = 2N-bi
                        mbi = minus_ix[bi]
                        ##ch = (Cv[j][ai,bi]+sym_type*Cv[j][ai,mbi])
                        #if verbose>0:
                        #    print "j,ai,mbi=",j,ai,mbi
                        #    mpc_set(ztmp.value,Cv[j][ai][mbi],rnd)
                        #    print "Cv[",j,"][",ai,mbi,"]=",ztmp
                        mpc_set(ch[0],Cv[j][ai][mbi],rnd)
                        mpc_mul_si(ch[0],ch[0],sym_type,rnd)
                        mpc_add(ch[0],ch[0],Cv[j][ai][bi],rnd)
                        mpc_set_si(ch21,0,rnd)
                        mpc_set_si(ch22,0,rnd)
                        #ch21=0; ch22=0
                        if is_zero(Cv[j][ai][bi])==0: # i.e. is non-zero
                            mpc_set_ui(ch21,1,rnd)
                            mpc_div(ch21,ch21,Cv[j][ai][bi],rnd)
                            #ch21=Cv[j][ai,bi]**-1
                        if is_zero(Cv[j][ai][mbi])==0:
                            mpc_set_ui(ch22,1,rnd)
                            mpc_div(ch22,ch22,Cv[j][ai][mbi],rnd)
                            #ztmp=Cv[j][ai,mbi]**-1
                        if sym_type==1:
                            mpc_add(ch2,ch22,ch21,rnd)
                        else:
                            mpc_sub(ch2,ch22,ch21,rnd)
                    if is_zero(ch[0])==0:
                        mpc_set(tmp1[0],fak1[l][bi-Ds][j],rnd)
                        mpc_mul(tmp1[0],tmp1[0],ch[0],rnd)
                        mpc_mul(tmp1[0],tmp1[0],Cfak_plus._entries[j],rnd)
                        #tmp1[0]=fak1[l][j-1,bi-Ds]*ch*Cfak[j][0]
                    else:
                        mpc_set_si(tmp1[0],0,rnd)
                    if is_zero(ch2)==0:
                        mpc_set(tmp2[0],fak1[l][bi-Ds][j],rnd)
                        mpc_conj(tmp2[0],tmp2[0],rnd)
                        mpc_mul(tmp2[0],tmp2[0],ch2,rnd)
                        mpc_mul(tmp2[0],tmp2[0],Cfak_minus._entries[j],rnd)
                    else:
                        mpc_set_si(tmp2[0],0,rnd)
                    if is_zero(ch[0]) and is_zero(ch2):
                        continue
                    nis= Ml*ai
                    ## Do the case n=0 separately
                    if not H._holomorphic or Qvv[ai]>=0:
                        mpc_set(tmpV,fak2[0][ai][j],rnd)
                        mpc_mul(tmpV2,tmp1[0],tmpV,rnd)
                        mpc_conj(tmpV,tmpV,rnd)
                        mpc_mul(tmpV,tmpV,tmp2[0],rnd)                        
                        mpc_add(tmpV1.value,tmpV,tmpV2,rnd)
                        mpc_add(V._matrix[nis][lj],V._matrix[nis][lj],tmpV1.value,rnd)                        
                    for n from 1 <= n < Ml:
                        ni=nis+n
                        mpc_set(tmpV,fak2[n][ai][j],rnd)
                        mpc_mul(tmpV2,tmp1[0],tmpV,rnd)
                        mpc_conj(tmpV,tmpV,rnd)
                        mpc_mul(tmpV,tmpV,tmp2[0],rnd)                        
                        mpc_add(tmpV1.value,tmpV,tmpV2,rnd)
                        #tmpV1=tmp1[0]*ztmp+tmp2[0]*ztmp.conjugate()
                        #tmpV1=(tmp1[0]*ztmp+tmp2[0]*ztmp.conjugate())
                        mpc_add(V._matrix[ni][lj],V._matrix[ni][lj],tmpV1.value,rnd)
    if verbose>0:
        print "here1"
    #print "maxv=",maxv
    #return Cv
    if verbose>0:
        if V.nrows()>=9:
            print "V[9,9]=",V[9,9]
        if V.nrows()>=30:
            print "V[30,30]=",V[30,30]
    cdef mpfr_t twopiY[1]
    mpfr_init2(twopiY[0],prec)
    mpfr_mul(twopiY[0],twopi.value,Y.value,rnd_re)
    ni = V.nrows()-1
    lj = V.ncols()-1
    for n from 0 <= n <= ni:
        for l from 0 <= l <= lj:
            mpc_div_ui(V._matrix[n][l],V._matrix[n][l],Qfaki,rnd)
    maxb=0
    maxw=0
    cdef MPComplexNumber an
    an = CF(0)

    for ai from 0 <= ai < Dl:
        for n from 0 <= n < Ml:
            #a=D[ai+Ds]
            mpfr_set(nr[0],Qvv._entries[ai],rnd_re)
            mpfr_add_si(nr[0],nr[0],n+Ms,rnd_re)
            mpc_set_fr(nri[0],nr[0],rnd)
            mpfr_swap(mpc_realref(nri[0]),mpc_imagref(nri[0]))
            ni=Ml*ai+n
            mpfr_mul(nrY2pi.value,nr[0],twopiY[0],rnd_re)
            if verbose>1:
                print "ni=",ni
                mpfr_set(tmp_r.value,nr[0],rnd_re)
                print "nr=",tmp_r
            if mpfr_sgn(nr[0])>=0:
                kbes=(-nrY2pi).exp()
            elif not H._holomorphic:
                mpfr_mul_si(nrY4pi.value,nrY2pi.value,-2,rnd_re)
                #if verbose>0:
                #    print "nrY4pi=",nrY4pi
                ok = 0
                if is_int==1:
                    if kinti>0:
                        ok = incgamma_pint_c(kbes.value,kinti,nrY4pi.value,verbose)
                    else:
                        ok = incgamma_nint_c(kbes.value,kinti,nrY4pi.value,verbose)
                elif is_half_int==1:
                    ok = incgamma_hint_c(kbes.value,kinti,nrY4pi.value,verbose)
                if not ok:
                    tu = mpfr_to_mpfval(nrY4pi.value)
                    mparg._set_mpf(tu)
                    testmp=mpctx.gammainc(kint,nrY4pi)
                    mpfr_from_mpfval(kbes.value,testmp._mpf_)
                kbes=kbes*(-nrY2pi).exp()
            else:
                #continue
                if verbose>1:
                    print "kbes=0 (maybe skipping?)"
                kbes = RF(0) #continue
            mpc_sub_fr(V._matrix[ni][ni],V._matrix[ni][ni],kbes.value,rnd)
            mpc_set_si(summa.value,0,rnd) #=CF(0)
            #ni=Ml*ai+n-Ms
            for bi from 0<= bi < len_of_pp:
                #if verbose>0 :
                #    print "bi=",bi
                # We only need the first half
                if beta_int[bi]<Ds or beta_int[bi]>Df:
                    if verbose>1:
                        print "Skipping beta=",beta," by symmetry!"
                        print "index=",beta_int[bi]
                        print "range=",range(Ds,Df+1)
                    continue
                an=CF(pp_c[bi])
                if an==0:
                    if verbose>1:
                        print "an=0, skipping!"
                    continue
                else:
                    if verbose>1:
                        print "an=",an
                mpc_set_si(summa.value,0,rnd) #=CF(0)
                if beta_int[bi] == mbeta_int[bi]:
                    if sym_type==-1:
                        if verbose>1:
                            print "cancel due to symmetry!"
                        
                        #summa=CF(0)
                    else:
                        for j from 0 <= j < Ql:
                            mpc_mul_fr(tmp1[0],nri[0],Xm._entries[j],rnd)
                            mpc_set(tmp.value,tmp1[0],rnd)
                            #tmp=mmi_cplx[bi]*Zpb[j]-tmp  
                            mpc_mul(tmp1[0],mmi_cplx[bi],Zpb._entries[j],rnd)
                            mpc_sub(tmp.value,tmp1[0],tmp.value,rnd)
                            #tmp=mmmi[(beta,m)]*Zpb[j]-nri[0]*Xm[j]
                            mpc_exp(tmpc.value,tmp.value,rnd) #=tmp.exp()
                            mpc_set(tmp1[0],Cv[j][ai][beta_int[bi]],rnd)
                            #tmp1[0]=tmp1[0]*Cfak[j][0] TEST
                            mpc_mul(tmp1[0],tmp1[0],Cfak_plus._entries[j],rnd)
                            #ch2=Cv[1-j][ai,betai[beta]]
                            mpc_mul(tmp1[0],tmp1[0],tmpc.value,rnd)
                            mpc_add(summa.value,summa.value,tmp1[0],rnd)
                            #rtmp1[0].summa=summa+tmpc*tmp1[0]
                            bii = minus_ix[beta_int[bi]]
                            #if beta_int[bi] == 0 or beta_int[bi]==N:
                            #    bii = beta_int[bi]
                            #else:
                            #    bii= abD - beta_int[bi] #
                            if is_zero(Cv[j][ai][bii])==0:
                                #t_ch_r==0 or t_ch_i==0:
                                # if Cv[j][ai,bii]<>0:
                                mpc_div(tmp1[0],Cfak_minus._entries[j],Cv[j][ai][bii],rnd)
                                mpc_set(tmp2[0],tmpc.value,rnd)
                                mpc_conj(tmp2[0],tmp2[0],rnd)
                                mpc_mul(tmp1[0],tmp1[0],tmp2[0],rnd)

                                mpc_add(summa.value,summa.value,tmp1[0],rnd)
                                #summa=summa+tmp1[0]*tmpc.conjugate()

                else:
                    for j from 0 <= j < Ql:
                        #ch=(Cv[j][ai,beta_int[bi]]+sym_type*Cv[j][ai,mbeta_int[bi]])
                        mpc_set(ch[0],Cv[j][ai][mbeta_int[bi]],rnd)
                        mpc_mul_si(ch[0],ch[0],sym_type,rnd)
                        mpc_add(ch[0],ch[0],Cv[j][ai][beta_int[bi]],rnd)
                        if is_zero(Cv[j][ai][beta_int[bi]])==0:
                            mpc_set_si(ch21,1,rnd)
                            #ch21=1/Cv[j][ai,beta_int[bi]]
                            mpc_div(ch21,ch21,Cv[j][ai][beta_int[bi]],rnd)
                        if is_zero(Cv[j][ai][mbeta_int[bi]])==0:
                            mpc_set_si(ch22,1,rnd)
                            mpc_div(ch22,ch22,Cv[j][ai][mbeta_int[bi]],rnd)
                            #ch22=1/Cv[j][ai,mbeta_int[bi]]
                        if sym_type==1:
                            mpc_add(ch2,ch21,ch22,rnd)
                            #ch2=(ch21+ch22)
                        else:
                            mpc_sub(ch2,ch22,ch21,rnd)
                            #ch2=(ch22-ch21)
                        if is_zero(ch[0]) and is_zero(ch2):
                            continue
                        mpc_mul_fr(tmp1[0],nri[0],Xm._entries[j],rnd)
                        #tmp=mmmi[(beta,m)]*Zpb[j]-nri[0]*Xm[j]
                        mpc_set((<MPComplexNumber>tmp).value,tmp1[0],rnd)
                        
                        #tmp=mmi_cplx[bi]*Zpb[j]-tmp # -tmp1[0]
                        #tmpc=tmp.exp()
                        mpc_mul(tmp1[0],mmi_cplx[bi],Zpb._entries[j],rnd)
                        mpc_sub(tmp.value,tmp1[0],tmp.value,rnd)
                        mpc_exp(tmpc.value,tmp.value,rnd)
                        
                        mpc_mul(ch[0],ch[0],Cfak_plus._entries[j],rnd)
                        mpc_mul(ch2,ch2,Cfak_minus._entries[j],rnd)
                        mpc_mul(ch[0],tmpc.value,ch[0],rnd)
                        mpc_conj(tmpc.value,tmpc.value,rnd)
                        mpc_mul(ch2,tmpc.value,ch2,rnd)
                        mpc_add(ch[0],ch[0],ch2,rnd)
                        mpc_add(summa.value,summa.value,ch[0],rnd)
                        #summa=summa+tmpc*ch+ch2*tmpc.conjugate()

                mpc_mul(summa.value,summa.value,an.value,rnd)
                mpc_div_ui(summa.value,summa.value,Qfaki,rnd)
                mpc_add(RHS._matrix[ni][0],RHS._matrix[ni][0],summa.value,rnd)
                #RHS[ni,0]=RHS[ni,0]+an*summa/Qfak
                if verbose>2:
                    print "summa=",n,ai,bi,summa
                    print "RHS[",ni,"]=",RHS[ni,0]
                    #print "ai=",ai
                    print "beta_int[bi]=",beta_int[bi]
                    #print "mbeta_int[bi]=",mbeta_int[bi]
                if n+Ms == m:
                    if mbeta_int[bi]==beta_int[bi]:
                        if ai+Ds==beta_int[bi] and sym_type==1:
                            mpfr_mul(tmpr_t[0],twopiY[0],mm_real[bi],rnd_re)
                            mpfr_neg(tmpr_t[0],tmpr_t[0],rnd_re)
                            mpfr_exp(tmpr_t[0],tmpr_t[0],rnd_re)
                            mpc_mul_fr(an.value,an.value,tmpr_t[0],rnd)
                            #an = an*(-mm_real[bi]*twopiY).exp()
                            RHS[ni,0]=RHS[ni,0]- an
                    elif ai+Ds==beta_int[bi]:
                        mpfr_mul(tmpr_t[0],twopiY[0],mm_real[bi],rnd_re)
                        mpfr_neg(tmpr_t[0],tmpr_t[0],rnd_re)
                        mpfr_exp(tmpr_t[0],tmpr_t[0],rnd_re)
                        mpc_mul_fr(an.value,an.value,tmpr_t[0],rnd)
                        #an = an*(-mm_real[bi]*twopiY).exp()
                        RHS[ni,0]=RHS[ni,0]-an #*(-mm_real[bi]*twopiY).exp()
                    elif ai+Ds == mbeta_int[bi]:
                        mpfr_mul(tmpr_t[0],twopiY[0],mm_real[bi],rnd_re)
                        mpfr_neg(tmpr_t[0],tmpr_t[0],rnd_re)
                        mpfr_exp(tmpr_t[0],tmpr_t[0],rnd_re)
                        mpc_mul_fr(an.value,an.value,tmpr_t[0],rnd)
                        mpc_mul_si(an.value,an.value,sym_type,rnd)
                        #an = sym_type*an*(-mm_real[bi]*twopiY).exp()
                        RHS[ni,0]=RHS[ni,0]-an  #*(-mm_real[bi]*twopiY).exp()


    # Clear variables
    mpc_clear(ch[0])
    mpc_clear(tmpV)
    mpc_clear(tmp1[0])
    mpc_clear(tmp2[0])
    mpc_clear(tmpV2)
    mpc_clear(ch2)
    mpc_clear(ch21)
    mpc_clear(ch22)
    mpfr_clear(nr[0])
    mpfr_clear(tmpr_t[0])
    mpfr_clear(twopiY[0])
    mpc_clear(nri[0])
    mpc_clear(tmpc_t[0])
    sage_free(beta_int)
    sage_free(mbeta_int)
    sage_free(mm_real)
    sage_free(mmi_cplx)
    sage_free(minus_ix)
    if Cv<>NULL:
        for j from 0 <= j < Ql: 
            if Cv[j]<>NULL:
                for n from 0 <= n <  Dl:
                    if Cv[j][n]<>NULL:
                        for l from 0<= l < rank:
                            #print "Cv[",j,n,l,"]=",printf("%p ",Cv[j][n][l])
                            #mpc_set(ztmp.value,Cv[j][n][l],rnd)
                            #print ztmp
                            mpc_clear(Cv[j][n][l])
                        sage_free(Cv[j][n])
                sage_free(Cv[j])
        sage_free(Cv)

    if fak1<>NULL:
        for n from 0<=n < Ml:
            if fak1[n]<>NULL:
                for ai from 0<=ai < Dl:
                    if fak1[n][ai]<>NULL:
                        for j from 0<=j<=Q-1:
                            mpc_clear(fak1[n][ai][j])
                        sage_free(fak1[n][ai])
                sage_free(fak1[n])            
        sage_free(fak1)
    if fak2<>NULL:
        for n from 0<=n < Ml:
            if fak2[n]<>NULL:
                for ai from 0<=ai < Dl:
                    if fak2[n][ai]<>NULL:
                        for j from 0<=j<=Q-1:
                            mpc_clear(fak2[n][ai][j])
                        sage_free(fak2[n][ai])
                sage_free(fak2[n])            
        sage_free(fak2)
    W=dict()
    #W['WR']=WR
    W['V']=V
    W['Ms']=Ms
    W['Mf']=Mf
    W['r']=Dl
    W['nc']=nc
    W['RHS']=RHS
    W['sym_type']=sym_type
    return W

cpdef vv_holomorphic_setupV_mpc(H,RealNumber Y,int M,int Q):
    r"""
    Set up the matrix for a vector-valued holomorphic modular form for the Weil representation.
    Using MPC

    INPUT:
    - H -- instance of VVHarmonicWeakMaassforms
        if has
    EXAMPLES::

    sage: WR=WeilRepDiscriminantForm(1,False)
    sage: Y=mpmath.mp.mpf(0.5); k=mpmath.mp.mpf(0.5)
    sage: M=10; Q=20; m=-1; beta=1/2
    sage: W= vv_harmonic_wmwf_setupV(WR,beta,m,Y,M,Q,k)


    """
    ### Only implemented for Weil representation
    X = H.multiplier()
    WM = X.weil_module()
    cdef list D = X.D()    
    if WM == None:
        if H._verbose>0:
            print "WM=",WM
        raise ValueError,"Need a space with a Weil representation! Got:{0}".format(H)    
    cdef int prec
    if hasattr(H,"_setupV_prec"):
        prec = int(max(Y.base_ring().prec(),H._setupV_prec))
    else:
        prec = int(Y.base_ring().prec())
    if H._verbose>0:
        print "prec in setupV=",prec
    RF=RealField(prec)
    CF=MPComplexField(prec)
    CCF=ComplexField(prec)
    if H._verbose>0:
        print "here!"
    if prec==53:
        mpctx = mpmath.fp
    else:
        mpctx = mpmath.mp
    mpctx = mpmath.mp
    cdef RealNumber mpone,mptwo,twopi,fourpi,weight,mp0,Qfak 
    cdef RealNumber nrtwopi,nrtwo,tmp_r
    cdef MPComplexNumber tmp,ztmp,nrtwopii,twopii,summa
    weight = RF(H._weight)
    p0=RF(0); mpone=RF(1); mptwo=RF(2); mpfour=RF(4)
    twopi=mptwo*RF.pi();  fourpi=mptwo*twopi
    weight_over_two=weight/mptwo
    twopii=CF(0,twopi); summa=CF(0)
    ztmp=CF(1);  nrtwopii=CF(1); tmp=CF(0); tmp_r=RF(0)
    cdef int j,l,n,ai,bi,bii,mbi
    cdef int N,sym_type,verbose,Qs,Qf,Ql,Ms,Mf,Ml,Qfaki
    #cdef list PP0,beta_rat,pp_c,mbeta_rat #,mm_real,mmi_cplx
    cdef mpfr_t tmpr_t[1]
    mpfr_init2(tmpr_t[0],prec)
    sym_type = H.sym_type(); verbose = H._verbose
    beta_rat=[]; pp_c=[]; mbeta_rat=[] #; mm_real=[]; mmi_cplx=[]
    cdef mpfr_t *mm_real
    cdef mpc_t *mmi_cplx
    if verbose>0:
        print "in vv_holomorphic_setupV_mpc!"
        print "args=",H,Y,M,Q
        print "D=",D
    Ms=int(0); Mf=int(M); Ml=int(Mf-Ms+1)
    cdef unsigned int abD
    abD = WM.rank()
    assert abD>0
    if verbose>0:
        print "Ms,Mf=",Ms,Mf
        print "beta=",beta_rat
        print "betai="
    if Q==0:
        Q=M+20
    # We use symmetry here
    if verbose>0:
        print "Symmetry type=",sym_type
        print "D=",D
    if sym_type<>0:
        Qs=1; Qf=Q; Ql=Qf-Qs; Qfaki=2*Q
    else:
        Qs=1-Q; Qf=Q; Ql=Qf-Qs; Qfaki=4*Q
    kint=mpmath.mp.mpf(mpone-weight)
    ## Recall that we have different index sets depending on the symmetry
    cdef int Ds,Df,Dl,s
    if sym_type==1:
        Ds =D[0]; Df=D[-1] #len(D)-1  #int(0); Df=int(N) # 0,1,...,N
    elif sym_type==-1:
        Ds=D[0]; Df=D[-1] #len(D)-1 #int(N-1) # 1,2,...,N-1  (since -0=0 and -N=N)
    if verbose>0:
        print "Ds=",Ds,"Df=",D
    Dl = Df - Ds + 1
    assert Dl==len(D)
    cdef Vector_real_mpfr_dense Qvv
    cdef mpc_t ***Cv
    Cv = NULL
    Cv = <mpc_t***> sage_malloc ( sizeof(mpc_t**)*Ql)
    if Cv==NULL: raise MemoryError
    cdef int nc = H.group().ncusps() #len(WM.basis())
    cdef int rank = X.rank()
    cdef int full_rank = len(X.Qv)
    if verbose>0:
        print "precBBBB:",prec
    for j from 0 <= j < Ql: 
        Cv[j] = NULL
        Cv[j]=<mpc_t**>sage_malloc(sizeof(mpc_t*)*Dl)
        if Cv[j]==NULL: raise MemoryError
        for n from 0 <= n <  Dl:
            Cv[j][n] = NULL
            Cv[j][n]=<mpc_t*>sage_malloc(sizeof(mpc_t)*full_rank) 
            if Cv[j][n]==NULL: raise MemoryError
            for l from 0<= l < full_rank:
                mpc_init2(Cv[j][n][l],prec) 
    cdef Vector_real_mpfr_dense Xm
    cdef Vector_complex_dense Zpb,Cfak_plus,Cfak_minus
    Xm = Vector_real_mpfr_dense(FreeModule(RF,Ql),0)
    Zpb= Vector_complex_dense(FreeModule(CF,Ql),0)
    Cfak_plus = Vector_complex_dense(FreeModule(CF,Ql),0)
    Cfak_minus = Vector_complex_dense(FreeModule(CF,Ql),0)
    if verbose>0:
        print "Ql=",Ql
        print "Dl=",Dl
        print "rank=",rank
    pullback_pts_weil_rep_mpc2(H,Q,Y,weight,Xm,Zpb,Cfak_plus,Cfak_minus,Cv,Ds,Df)
    Dl=Df-Ds+1
    s=Dl*Ml
    if verbose>0:
        print "Ds, Df=",Ds,Df,Dl
        print "Ml=",Ml
        print "abD=",abD
        print "matrix size=",s
    MS=MatrixSpace(CF,s,s)
    cdef Matrix_complex_dense V #,RHS 
    V=Matrix_complex_dense(MS,0)
    MS=MatrixSpace(CF,s,1)
    #RHS=Matrix_complex_dense(MS,0)
    MS = MatrixSpace(CF,int(Df-Ds+1),int(Q))
    if verbose>0:
        print "k/2=",weight_over_two
    VS = vector(RF,Df-Ds+1).parent()
    Qvv=Vector_real_mpfr_dense(VS,0)
    nrtwopi=RF(0); nrtwo=RF(0)
    cdef mpc_t tmp1[1],tmp2[1]
    mpc_init2(tmp1[0],prec)
    mpc_init2(tmp2[0],prec)
    for j from Ds<= j <= Df:
        Qvv[j-Ds]=RF(H.multiplier().Qv[j])
    cdef mpfr_t nr[1]
    cdef mpc_t nri[1],tmpc_t[1]
    cdef Rational a
    a = QQ(0)
    mpfr_init2(nr[0],prec)
    mpc_init2(nri[0],prec)
    mpc_init2(tmpc_t[0],prec)
    cdef RealNumber tmparg
    tmparg=RF(1)
    cdef RealNumber nrY2pi,nrY4pi,kbes
    kbes=RF(0)
    nrY2pi=RF(0)
    nrY4pi=RF(0)
    cdef MPComplexNumber tmpc
    tmpc=CF(0)
    cdef mpc_t ***fak1,***fak2
    fak1=NULL; fak2=NULL
    fak1 =<mpc_t ***> sage_malloc( Ml*sizeof(mpc_t**))
    for n from 0 <= n < Ml:
        fak1[n]=<mpc_t **> sage_malloc( Dl*sizeof(mpc_t*))
        for ai from 0 <= ai < Dl:
            fak1[n][ai]=<mpc_t *>sage_malloc(Q*sizeof(mpc_t))
            for j from 0 <= j < Ql:
                mpc_init2(fak1[n][ai][j],prec)
    fak2 =<mpc_t ***> sage_malloc( Ml*sizeof(mpc_t**))
    for n from 0 <= n < Ml:
        fak2[n]=<mpc_t **> sage_malloc( Dl*sizeof(mpc_t*))
        for ai from 0 <= ai < Dl:
            fak2[n][ai]=<mpc_t *>sage_malloc(Q*sizeof(mpc_t))
            for j from 0 <= j < Q:
                mpc_init2(fak2[n][ai][j],prec)

    cdef tuple tu
    mparg = mpmath.mp.mpf(0)
    for n from 0 <= n < Ml:
        for ai from 0 <= ai < Dl:
            mpfr_set_si(nr[0],n+Ms,rnd_re)
            mpfr_add(nr[0],nr[0],Qvv._entries[ai],rnd_re)
            mpc_set_fr(nri[0],nr[0],rnd)
            mpfr_swap(mpc_realref(nri[0]),mpc_imagref(nri[0]))
            mpfr_set(tmp_r.value,nr[0],rnd_re)
            mpc_set(ztmp.value,nri[0],rnd)
            for j from 0 <= j < Ql:
                mpc_mul_fr(tmpc_t[0],nri[0],Xm._entries[j],rnd)
                mpc_neg(tmpc_t[0],tmpc_t[0],rnd)
                mpc_exp(fak2[n][ai][j],tmpc_t[0],rnd)
            if mpfr_sgn(nr[0])>0:
                for j in range(Ql): 
                    mpc_mul(tmpc_t[0],nri[0],Zpb._entries[j],rnd)
                    mpc_exp(tmpc_t[0],tmpc_t[0],rnd)
                    mpc_set(fak1[n][ai][j],tmpc_t[0],rnd)
            elif mpfr_sgn(nr[0])==0:
                for j in range(Ql):
                    mpc_set_ui(fak1[n][ai][j],1,rnd)                
    if verbose>5:
        for j from 0 <= j < Ql:
            for ai from 0 <= ai < Dl:
                for bi from 0 <= bi < rank:
                    mpc_set(ztmp.value,Cv[j][ai][bi],rnd)
                    print "Cv[",j,"][",ai,bi,"]=",ztmp
    cdef int t_ch_r,t_ch_i,t_ch2_r,t_ch2_i,ni,lj,nis
    cdef mpc_t ch[1],ch2,ch21,ch22,tmpV,tmpV2
    if verbose>0:
        print "precAAA:",prec
    mpc_init2(tmpV2,prec)
    mpc_init2(tmpV,prec)
    mpc_init2(ch[0],prec)
    mpc_init2(ch2,prec)
    mpc_init2(ch21,prec)
    mpc_init2(ch22,prec)
    cdef MPComplexNumber tmpV1
    cdef int order = WM.finite_quadratic_module().order()
    cdef int* minus_ix
    minus_ix = <int*>sage_malloc(sizeof(int)*order)
    for i in range(order):
        minus_ix[i] = WM._neg_index(i)
    cdef int t
    tmpV1=CF(0)
    for l in range(Ml):
        for ai in range(Dl): #from 0 <= ai < Dl:
            for bi from 0 <= bi < Dl:
                lj=Ml*bi+l
                bii = bi + Ds
                mbii = minus_ix[bii]
                mbi = mbii - Ds
                if l==0 and Qvv[bi]<0:
                    if verbose>2:
                        print "Skipping l=",l," b=",bi, " lj=",lj
                    continue
                if bi==mbi and sym_type==-1:
                    raise ArithmeticError,"We have odd symmetry and D={0}".format(D)
                for j from 0 <= j < Ql:
                    if bi==mbi and sym_type==1:
                        mpc_set(ch[0],Cv[j][ai][bi],rnd)
                        if is_zero(Cv[j][ai][mbi])==0:  #<>0
                            mpc_set_ui(ch2,1,rnd)
                            mpc_div(ch2,ch2,Cv[j][ai][mbi],rnd)
                    else:
                        mpc_set(ch[0],Cv[j][ai][mbi],rnd)
                        mpc_mul_si(ch[0],ch[0],sym_type,rnd)
                        mpc_add(ch[0],ch[0],Cv[j][ai][bi],rnd)
                        mpc_set_si(ch21,0,rnd)
                        mpc_set_si(ch22,0,rnd)
                        if is_zero(Cv[j][ai][bi])==0: # i.e. is non-zero
                            mpc_set_ui(ch21,1,rnd)
                            mpc_div(ch21,ch21,Cv[j][ai][bi],rnd)
                        if is_zero(Cv[j][ai][mbi])==0:
                            mpc_set_ui(ch22,1,rnd)
                            mpc_div(ch22,ch22,Cv[j][ai][mbi],rnd)
                        if sym_type==1:
                            mpc_add(ch2,ch22,ch21,rnd)
                        else:
                            mpc_sub(ch2,ch22,ch21,rnd)
                    if is_zero(ch[0])==0:  ## ch<>0
                        mpc_set(tmp1[0],fak1[l][bi][j],rnd)
                        #if verbose>0:
                        #    mpc_set(tmpV1.value,tmp1[0],rnd)
                        #    print "fak1(",l,bi,j,")=",tmpV1
                        mpc_mul(tmp1[0],tmp1[0],ch[0],rnd)
                        #if verbose>0:
                        #    mpc_set(tmpV1.value,ch[0],rnd)
                        #    print "ch=",tmpV1
                        #    print "Cfak_plus(",j,")=",Cfak_plus[j]
                        mpc_mul(tmp1[0],tmp1[0],Cfak_plus._entries[j],rnd)
                    else:
                        mpc_set_si(tmp1[0],0,rnd)
                    if is_zero(ch2)==0:
                        mpc_set(tmp2[0],fak1[l][bi][j],rnd)
                        mpc_conj(tmp2[0],tmp2[0],rnd)
                        mpc_mul(tmp2[0],tmp2[0],ch2,rnd)
                        mpc_mul(tmp2[0],tmp2[0],Cfak_minus._entries[j],rnd)
                    else:
                        mpc_set_si(tmp2[0],0,rnd)
                    if is_zero(ch[0]) and is_zero(ch2):
                        continue
                    nis= Ml*ai
                    ## Do the case n=0 separately
                    if Qvv[ai]>=0:
                        mpc_set(tmpV,fak2[0][ai][j],rnd)
                        if verbose>0:
                            mpc_set(tmpV1.value,tmpV,rnd)
                            print "fak2(",0,ai,j,")=",tmpV1
                            mpc_set(tmpV1.value,tmp1[0],rnd)
                            print "tmp1[0](",j,n,l,")=",tmpV1
                        mpc_mul(tmpV2,tmp1[0],tmpV,rnd)
                        mpc_conj(tmpV,tmpV,rnd)
                        mpc_mul(tmpV,tmpV,tmp2[0],rnd)                        
                        mpc_add(tmpV1.value,tmpV,tmpV2,rnd)
                        if verbose>0:
                            print "tmpV1(",j,n,l,")=",tmpV1
                            print "V(",nis,lj,")=",V[nis,lj]
                            
                        mpc_add(V._matrix[nis][lj],V._matrix[nis][lj],tmpV1.value,rnd)                        
                    for n from 1 <= n < Ml:
                        ni=nis+n
                        mpc_set(tmpV,fak2[n][ai][j],rnd)
                        mpc_mul(tmpV2,tmp1[0],tmpV,rnd)
                        mpc_conj(tmpV,tmpV,rnd)
                        mpc_mul(tmpV,tmpV,tmp2[0],rnd)                        
                        mpc_add(tmpV1.value,tmpV,tmpV2,rnd)
                        mpc_add(V._matrix[ni][lj],V._matrix[ni][lj],tmpV1.value,rnd)
    if verbose>0:
        print "here1"
    if verbose>0:
        if V.nrows()>=9:
            print "V[9,9]=",V[9,9]
        if V.nrows()>=30:
            print "V[30,30]=",V[30,30]
    cdef mpfr_t twopiY
    mpfr_init2(twopiY,prec)
    mpfr_mul(twopiY,twopi.value,Y.value,rnd_re)
    ni = V.nrows()-1
    lj = V.ncols()-1
    for n from 0 <= n <= ni:
        for l from 0 <= l <= lj:
            mpc_div_ui(V._matrix[n][l],V._matrix[n][l],Qfaki,rnd)
    maxb=0
    maxw=0
    cdef MPComplexNumber an
    an = CF(0)
    for ai from 0 <= ai < Dl:
        for n from 0 <= n < Ml:
            #a=D[ai+Ds]
            mpfr_set(nr[0],Qvv._entries[ai],rnd_re)
            mpfr_add_si(nr[0],nr[0],n+Ms,rnd_re)
            mpc_set_fr(nri[0],nr[0],rnd)
            mpfr_swap(mpc_realref(nri[0]),mpc_imagref(nri[0]))
            ni=Ml*ai+n
            mpfr_mul(nrY2pi.value,nr[0],twopiY,rnd_re)
            if verbose>1:
                print "ni=",ni
                mpfr_set(tmp_r.value,nr[0],rnd_re)
                print "nr=",tmp_r
            if mpfr_sgn(nr[0])>=0:
                kbes=(-nrY2pi).exp()          
            else:
                if verbose>1:
                    print "kbes=0 (maybe skipping?)"
                kbes = RF(0) #continue
            mpc_sub_fr(V._matrix[ni][ni],V._matrix[ni][ni],kbes.value,rnd)
            mpc_set_si(summa.value,0,rnd) #=CF(0)
            

    # Clear variables
    mpc_clear(ch[0])
    mpc_clear(tmpV)
    mpc_clear(tmp1[0])
    mpc_clear(tmp2[0])
    mpc_clear(tmpV2)
    mpc_clear(ch2)
    mpc_clear(ch21)
    mpc_clear(ch22)
    mpfr_clear(nr[0])
    mpfr_clear(tmpr_t[0])
    mpfr_clear(twopiY)
    mpc_clear(nri[0])
    mpc_clear(tmpc_t[0])
    #sage_free(beta_int)
    #sage_free(mbeta_int)
    #sage_free(mm_real)
    #sage_free(mmi_cplx)
    sage_free(minus_ix)
    if Cv<>NULL:
        for j from 0 <= j < Ql: 
            if Cv[j]<>NULL:
                for n from 0 <= n <  Dl:
                    if Cv[j][n]<>NULL:
                        for l from 0<= l < rank:
                            #print "Cv[",j,n,l,"]=",printf("%p ",Cv[j][n][l])
                            #mpc_set(ztmp.value,Cv[j][n][l],rnd)
                            #print ztmp
                            mpc_clear(Cv[j][n][l])
                        sage_free(Cv[j][n])
                sage_free(Cv[j])
        sage_free(Cv)

    if fak1<>NULL:
        for n from 0<=n < Ml:
            if fak1[n]<>NULL:
                for ai from 0<=ai < Dl:
                    if fak1[n][ai]<>NULL:
                        for j from 0<=j<=Q-1:
                            mpc_clear(fak1[n][ai][j])
                        sage_free(fak1[n][ai])
                sage_free(fak1[n])            
        sage_free(fak1)
    if fak2<>NULL:
        for n from 0<=n < Ml:
            if fak2[n]<>NULL:
                for ai from 0<=ai < Dl:
                    if fak2[n][ai]<>NULL:
                        for j from 0<=j<=Q-1:
                            mpc_clear(fak2[n][ai][j])
                        sage_free(fak2[n][ai])
                sage_free(fak2[n])            
        sage_free(fak2)
    W=dict()
    #W['WR']=WR
    W['V']=V
    W['Ms']=Ms
    W['Mf']=Mf
    W['r']=Dl
    W['nc']=nc
    #W['RHS']=RHS
    W['sym_type']=sym_type
    return W



def is_int(q):
    r"""
    Find out if the rational number q is an integer.

    INPUT:
    -''q'' -- integer/rational/real

    OUTPUT:
    - logical -- True if q is an integer otherwise False

    EXAMPLES::

        sage: is_int(1)
        True
        sage: is_int(float(1.0))
        True
        sage: is_int(RR(1.0))   
        True
        sage: is_int(6/3) 
        True
        sage: is_int(6/4)
        False
        sage: is_int(1.5)    
        False
        sage: is_int(Gamma0(1))    
        False
        
    """
    if isinstance(q,Integer) or isinstance(q,int):
        return True
    if isinstance(q,Rational):
        n=q.denominator()
        if n==1:
            return True
    try:
        if floor(q)==pceil(q):
            return True
    except:
        pass
    return False

def pullback_pts_weil_rep(WR,Q,Y,weight,ds=None,df=None,dual_rep=True,deb=False):
    raise NotImplementedError,"This method is obsolete!"

def pullback_pts_weil_rep_ef(WR,Q,Y,weight,ds=None,df=None,dual_rep=1,verbose=0):
    r"""
    rho = representation of SL2Z, i.e. a map :SL2Z->Mat_r(C)


    -''WR'' -- Finite quadratic module or Discriminant form
    -''Q''  -- integer
    -''Y''  -- real
    -''weight''-- real
    -''ds'' -- integer
    -''df'' -- integer
    -''dual_rep''-- logical (default True)
    -''verbose''-- int (default 0)

    EXAMPLES::
    sage: DF=DiscriminantForm(1,False)
    sage: Y=mpmath.mp.mpf(0.5); Q=20; k=mpmath.mp.mpf(0.5)
    sage: pullback_pts_weil_rep(DF,Q,Y,k,False)
    
    """
    try:
        Y.ae; weight.ae
    except AttributeError:
        raise Exception,"Need Y and weight in mpmath format for pullback!"
    twopi=mpmath.mp.mpf(2)*mpmath.pi()
    twopii=mpmath.mp.mpc(0,twopi)
    mp_i=mpmath.mp.mpc(0,1)
    Xm=dict();  Xpb=dict() 
    Ypb=dict(); Cvec=dict()
    Cfak=dict()
    cdef int Qs,Qf,Ql,Qfaki
    Qs=1; Qf=Q
    Q = Qf - Qs + 1
    Qfaki = 4*Q
    fourQ=mpmath.mp.mpf(4*Q)
    cdef int j,ri,rj
    #if df==None:
    #    df=WR.rank()-1 # use full range
    #    ds=0
    ds = WR.D()[0];    df = WR.D()[-1]
    cdef int N,twoN
    #N = WR.rank()/2
    twoN=WR.weil_module().rank()
    N2=my_odd_part(twoN); N22=1-N2
    z=CyclotomicField(4).gen()        

    # Recall that Az=z* in F <=> JAJ(Jz)=Jz* in F
    # so we only need to do the pullback for half of the points
    # since Xm[1-j]=-X[j]
    #for j from Qs <= j <= Qf: #1-Q,Q):
    for j from 1 <= j <= Q: #1-Q,Q):
        Xm[j]=mpmath.mp.mpf(2*j-1)/fourQ
        [a,b,c,d]=pullback_to_psl2z_mat(float(Xm[j]),float(Y))
        ## for simplicity always choose the pullback matrix such that the inverse Tj has c>0
        ##

        ar=mpmath.mp.mpf(a); br=mpmath.mp.mpf(b)
        cr=mpmath.mp.mpf(c); dr=mpmath.mp.mpf(d)
        den=mpmath.power(cr*Xm[j]+dr,2)
        den=den+mpmath.power(cr*Y,2)
        y=Y/den
        x=(ar*cr*(Xm[j]*Xm[j]+Y*Y)+Xm[j]*(ar*dr+br*cr)+br*dr)/den        
        #Tj=SL2Z([a,b,c,d])
        #print "Tj=",Tj
        Xpb[j]=x*twopi
        Ypb[j]=y*twopi
        if verbose>0:
            print "Xm=",Xm[j]
            print "Xpb=",Xpb[j]
            print "Ypb=",Ypb[j]        
        Xm[j]=Xm[j]*twopi
        if c>=0:
            #ar=mpmath.mp.mpf(-a); br=mpmath.mp.mpf(-b)
            #cr=mpmath.mp.mpf(-c); dr=mpmath.mp.mpf(-d)
            Tj=SL2Z([-d,b,c,-a])
            # j_Tm(z)
            jfak=mpmath.power(mpmath.mp.mpc(cr*x-ar,cr*y),weight)
            # j_Tm*(z)
            jfak2=mpmath.power(mpmath.mp.mpc(-cr*x+ar,cr*y),weight)
            c2=my_odd_part(c)            
        else:
            Tj=SL2Z([d,-b,-c,a])
            jfak=mpmath.power(mpmath.mp.mpc(-cr*x+ar,-cr*y),weight)
            jfak2=mpmath.power(mpmath.mp.mpc(cr*x-ar,-cr*y),weight)
            c2=my_odd_part(-c)
        Nc=my_gcd(twoN,c)
        N2c=my_odd_part(Nc)
        if is_odd(c):
            arg=c2*N2
        else:
            arg=c2+N22
        arg=arg-1-N2/Nc
        if dual_rep==1:
            arg=-arg
        rr=WR(Tj)
        r=rr[0]
        fak=mpmath.mp.mpf(rr[1])*jfak
        fak2=mpmath.mp.mpf(rr[1])*jfak2
        if verbose>0:
            print "jfak=",jfak
            print "rr1=",rr[1]
            print "fak=",fak
            print "fak2=",fak2
        #return 
        Cvec[j]=mpmath.matrix(int(r.nrows()),int(r.ncols()))
        # j=1-j: J[a,b,c,d]J^-1=[-a,b,c,-d]        
        fstar=kronecker(-1,abs(c)/Nc)*z**arg
        #for ri from ds <= ri <= df-1:
        Cfak[j]=[fak]
        #if cr<0:
        Cfak[1-j]=[fak2*fstar,-1]
        #else:
        #    Cfak[1-j]=[fak2*fstar,1]
        #if j==1:
        #    print "Tj=",Tj
        #    print "Cfak[j]=",Cfak[j]
        #    print "Cfak[1-j]=",Cfak[1-j]
        #    print "x=",x
        #    print "y=",y
        #    print "jfak=",jfak
        #    print "jfak=",jfak2
        #    print "a,b,c,d=",a,b,c,d
        if dual_rep==1:
            for ri from ds <= ri <= df:
                for rj from 0 <= rj <= r.ncols()-1:                    
                    if r[ri,rj]<>0:
                        Cvec[j][ri,rj]=1/r[ri,rj] #mpmath.mp.mpc(r[ri,rj]).conjugate()
                    else:
                        Cvec[j][ri,rj]=0
        else:
            for ri from ds <= ri <= df:
                for rj from 0 <= rj <= r.ncols()-1:                    
                    Cvec[j][ri,rj]=r[ri,rj]
                    #print "Cvec0[",j,ri,rj,"]=",Cvec[j][ri,rj]
        #raise Exception
    return [Xm,Xpb,Ypb,Cvec,Cfak]




cdef pullback_pts_weil_rep_mpc(H,int Q,RealNumber Y,RealNumber weight,Vector_real_mpfr_dense Xm,Vector_complex_dense Zpb,Vector_complex_dense Cfak_plus,Vector_complex_dense Cfak_minus, mpc_t *** Cvec,int ds=-1,int df=-1):
    r"""
    rho = representation of SL2Z, i.e. a map :SL2Z->Mat_r(C)

    -''WR'' -- Finite quadratic module or Discriminant form
    -''Q''  -- integer
    -''Y''  -- real
    -''weight''-- real
    -''ds'' -- integer
    -''df'' -- integer
    -''verbose''-- int  (default 0)

    EXAMPLES::
    sage: DF=DiscriminantForm(1,False)
    sage: Y=mpmath.mp.mpf(0.5); Q=20; k=mpmath.mp.mpf(0.5)
    sage: pullback_pts_weil_rep(DF,Q,Y,k,False)
    
    """
    cdef RealField_class RF
    cdef MPComplexNumber fak2,jfak,jfak2,twopii,mp_i,fstar,ztmp,chi
    cdef RealNumber twopi,YY,r1,x,y,xtmp
    cdef int prec,Qs,Qf,fourQ,Ql,arg
    cdef int j,ri,rj
    cdef mpfr_t Ypt
    cdef int jj
    cdef int a,b,c,d,N
    cdef mpfr_t xx,yy,xpb,ypb
    cdef mpc_t mptmp1,mptmp2,z0
    cdef dict r_vals
    cdef int c2,N2,N22,Nc
    cdef int rjmax,dual_rep,aa,bb
    cdef tuple Tj
    cdef double xt,yt
    cdef int verbose = H._verbose
    cdef int sym_type = H.sym_type()
    r_vals=dict()
    #sig_on()
    WR = H.multiplier().weil_module()
    dual_rep = H.multiplier().is_dual()
    if hasattr(H,"_pullback_prec"):
        prec = H._pullback_prec ## Try out what happend with higher prec here
    else:
        prec = H._prec ## Try out what happend with higher prec here
    #prec = H._pullback_prec ## Try out what happend with higher prec here
    RF = RealField(prec)
    CF = MPComplexField(prec)
    mpc_init2(mptmp1,prec)
    mpc_init2(mptmp2,prec)
    mpc_init2(z0,prec)
    mpfr_init2(Ypt,prec)
    mpfr_init2(xx,prec)
    mpfr_init2(xpb,prec)
    mpfr_init2(ypb,prec)
    mpfr_init2(yy,prec)
    mpfr_set(Ypt,Y.value,rnd_re)
    x=RF(1); y=RF(1); r1=RF(0); YY = RF(Y)
    xtmp = RF(1); ztmp = CF(1)
    jfak = CF(1); jfak2=CF(1)
    chi = CF(1)
    twopi=RF(2)*RF.pi()
    twopii=CF(0,twopi)
    fstar=CF(1)
    mp_i=CF(0,1)


    #cdef Vector_real_mpfr_dense Xm
    #cdef Vector_complex_dense Zpb 
    #Xm=dict();  Xpb=dict(); Ypb=dict();
    #Xm = Vector_real_mpfr_dense(FreeModule(RF,Q+1),0)
    #Zpb= Vector_complex_dense(FreeModule(CF,Q+1),0)
    ## Cv has to be allocated outside this routine
    assert Cvec<>NULL
    #Cfak=dict()
    Qs=1; Qf=Q; Ql = Qf - Qs + 1; fourQ=4*Q
    #N = WR.rank()/2
    cdef int twoN=WR.weil_module().rank()
    if df <0:
        df=twoN-1 # use full range
        ds=0
    rjmax = len(WR.basis())-1

    mpfr_const_pi(xx,rnd_re)
    mpfr_div_si(xx,xx,2*twoN,rnd_re)
    mpc_mul_fr(z0,mp_i.value,xx,rnd)
    mpc_exp(z0,z0,rnd)
    N2=my_odd_part(twoN); N22=1-N2
    z=CF(0,1) #CyclotomicField(4).gen()        
    #MS = MatrixSpace(CF,int(len(WR.D)),int(len(WR.D)))
    # Recall that Az=z* in F <=> JAJ(Jz)=Jz* in F
    # so we only need to do the pullback for half of the points
    # since Xm[1-j]=-X[j]
    ## NOte: This works because we have a symmetric domain for SL(2,Z)
    for j from 0 <= j < Ql: #1-Q,Q):
        #jj = 2*j-1
        jj = 2*j+2*Qs-1
        mpfr_set_si(x.value,jj,rnd_re)
        mpfr_div_ui(x.value,x.value,fourQ,rnd_re)
        xt=mpfr_get_d(x.value,rnd_re)
        yt=mpfr_get_d(Ypt,rnd_re)
        pullback_to_psl2z_mat_c(&xt,&yt,&a,&b,&c,&d)
        if verbose>2:
            mpfr_set(xtmp.value,x.value,rnd_re)
            print "Xm[",j,"]=",xtmp
        mpfr_set(xpb,x.value,rnd_re)
        mpfr_set(ypb,Ypt,rnd_re)
        _apply_sl2z_map_mpfr(xpb,ypb,a,b,c,d)
        if verbose>2:
            mpfr_set(xtmp.value,xpb,rnd_re)
            print "Xpb[",j,"]=",xtmp            
            mpfr_set(xtmp.value,ypb,rnd_re)
            print "Xpb[",j,"]=",xtmp            
        if c>=0:
            Tj=(-d,b,c,-a)
            # j_Tm(z)
            # xx = cx-a+icy; yy = cy
            mpfr_mul_si(xx,xpb,c,rnd_re)
            mpfr_sub_si(xx,xx,a,rnd_re)
            mpfr_mul_si(yy,ypb,c,rnd_re)            
            mpc_set_fr_fr(jfak.value,xx,yy,rnd)
            mpc_pow_fr(jfak.value,jfak.value,weight.value,rnd)
            #jfak=CF(cr*x-ar,cr*y)**weight
            # j_Tm*(z)
            mpfr_mul_si(xx,xpb,-c,rnd_re)
            mpfr_add_si(xx,xx,a,rnd_re)
            mpfr_mul_si(yy,ypb,c,rnd_re)            
            mpc_set_fr_fr(jfak2.value,xx,yy,rnd)
            mpc_pow_fr(jfak2.value,jfak2.value,weight.value,rnd)
            #jfak2=CF(-cr*x+ar,cr*y)**weight
            c2=my_odd_part(c)            
        else:
            Tj=(d,-b,-c,a)
            mpfr_mul_si(xx,xpb,-c,rnd_re)
            mpfr_add_si(xx,xx,a,rnd_re)
            mpfr_mul_si(yy,ypb,-c,rnd_re)            
            mpc_set_fr_fr(jfak.value,xx,yy,rnd)
            mpc_pow_fr(jfak.value,jfak.value,weight.value,rnd)
            #jfak=CF(-cr*x+ar,-cr*y)**weight
            mpfr_mul_si(xx,xpb,c,rnd_re)
            mpfr_sub_si(xx,xx,a,rnd_re)
            mpfr_mul_si(yy,ypb,-c,rnd_re)            
            mpc_set_fr_fr(jfak2.value,xx,yy,rnd)
            mpc_pow_fr(jfak2.value,jfak2.value,weight.value,rnd)
            #jfak2=CF(cr*x-ar,-cr*y)**weight
            c2=my_odd_part(-c)

        mpfr_mul(xpb,xpb,twopi.value,rnd_re)
        mpfr_mul(ypb,ypb,twopi.value,rnd_re)
        mpfr_mul(Xm._entries[j],x.value,twopi.value,rnd_re)
        mpc_set_fr_fr(Zpb._entries[j],xpb,ypb,rnd)
        if verbose>2:
            print "Xm=",Xm[j]
            print "iZpb=",Zpb[j]
        Nc=my_gcd(twoN,c)
        N2c=my_odd_part(Nc)
        if my_is_odd(c):
            arg=c2*N2
        else:
            arg=c2+N22
        arg=arg-1-N2/Nc
        if dual_rep==1:
            arg=-arg
        if Tj not in r_vals.keys():
            rr = WR.matrix(SL2Z(Tj),numeric=-1,prec=prec)
            r = rr[0]; r1=RF(rr[1]); r_ind=rr[2];
            try:
                chi = CF(rr[3].complex_embedding(prec))
            except Exception as e:
                print "rr3=",rr[3]
                print "rr3.c_e()=",rr[3].complex_embedding(prec)
                raise e
            r_vals[Tj] = (r,r1,r_ind,chi) #(rr[0],rr[1],rr[2],CF(rr[3].complex_embedding(prec)))
        else:
            (r,r1,r_ind,chi)= r_vals[Tj]

        ## do not want to large dict...
        if len(r_vals.keys())>100:
            r_vals.pop(r_vals.keys()[0])
        #r=rr[0]
        #r1=rr[1]
        #r_ind=rr[2]
        #chi = rr[3]
        fak=r1*jfak
        fak2=r1*jfak2
        if verbose>2:
            print "arg=",arg
            print "jfak=",jfak
            print "rr1=",rr[1]
            print "fak=",fak
            print "fak2=",fak2
        #fstar=kronecker(-1,abs(c)/Nc)*z**arg
        fstar = CF(my_kronecker_minus_one_red(c,Nc))*z**arg
        #Cfak[j]=[fak]
        Cfak_plus[j]=CF(fak)
        fak2=fak2*fstar #.complex_embedding(prec)
        #Cfak[1-j]=[fak2,-1]
        Cfak_minus[j]=CF(fak2)

        for ri from ds <= ri <= df:
            for rj from 0 <= rj <= rjmax:                    
                #if rj == N or rj ==0:
                if r_ind[ri,rj]<>0:
                    aa = r[ri][rj]
                    mpc_pow_si(Cvec[j][ri-ds][rj],z0,aa,rnd)
                else:
                    mpc_set_si(Cvec[j][ri-ds][rj],0,rnd)
                #else:
                #    if r_ind[ri][rj]<>0:
                #        aa = r[ri][rj]
                #        mpc_pow_si(mptmp1,z0,aa,rnd)
                #    else: 
               #        mpc_set_ui(mptmp1,0,rnd)
                #    if r_ind[ri][rj]<>0:
                #        bb = r[ri][twoN-rj]
                #        mpc_pow_si(mptmp1,z0,bb,rnd)
                #    if sym_type==1:
                #        mpc_add(Cvec[j][ri-ds][rj],mptmp1,mptmp2,rnd)
                #    else:
                #        mpc_sub(Cvec[j][ri-ds][rj],mptmp1,mptmp2,rnd)
                mpc_mul(Cvec[j][ri-ds][rj],Cvec[j][ri-ds][rj],chi.value,rnd)
                if dual_rep==1:
                    mpc_conj(Cvec[j][ri-ds][rj],Cvec[j][ri-ds][rj],rnd)

        #raise Exception
    if verbose>0:
        print "num of maps=",len(r_vals.keys())
    mpfr_clear(xx)
    mpfr_clear(yy)
    mpfr_clear(xpb)
    mpfr_clear(ypb)
    mpfr_clear(Ypt)
    mpc_clear(mptmp1)
    mpc_clear(mptmp2)
    #sig_off()
    #return [Xm,Zpb,Cfak]


cpdef pullback_pts_vv_mpc(H,int Q,RealNumber Y,RealNumber weight,int verbose=0):
    r"""
    Calls the cdef function and returns the vectors of pullback points etc.
    """
    cdef int prec = Y.prec()
    cdef int n,l,st,Qs,Qf,Ql,Ds,Df,Dl,Qfaki
    CF = MPComplexField(prec)
    RF = RealField(prec)
    cdef Vector_real_mpfr_dense Xm
    cdef Vector_complex_dense Zpb,Cfak_plus,Cfak_minus
    cdef mpc_t ***Cv=NULL
    cdef int nc = H.ambient_rank()
    st=H.sym_type()
    if st<>0:
        Qs=1; Qf=Q; Ql=Qf-Qs+1; Qfaki=2*Q
    else:
        Qs=1-Q; Qf=Q; Ql=Qf-Qs+1; Qfaki=4*Q
        #Ds=0; Df=nc; Dl=Df-Ds+1
    Cv = <mpc_t***> sage_malloc ( sizeof(mpc_t**)*Ql)
    if Cv==NULL: raise MemoryError
    Xm = Vector_real_mpfr_dense(FreeModule(RF,Ql),0)
    Zpb= Vector_complex_dense(FreeModule(CF,Ql),0)
    Cfak_plus = Vector_complex_dense(FreeModule(CF,Ql),0)
    Cfak_minus = Vector_complex_dense(FreeModule(CF,Ql),0)
    Ds=H.D()[0]; Df=H.D()[-1]; Dl=Df-Ds+1
    for j from 0 <= j < Ql: 
        Cv[j] = NULL
        Cv[j]=<mpc_t**>sage_malloc(sizeof(mpc_t*)*Dl)
        if Cv[j]==NULL: raise MemoryError
        for n from 0 <= n <  Dl:
            Cv[j][n] = NULL
            Cv[j][n]=<mpc_t*>sage_malloc(sizeof(mpc_t)*nc) 
            if Cv[j][n]==NULL: raise MemoryError
            for l from 0<= l < nc:
                mpc_init2(Cv[j][n][l],prec)
    pullback_pts_weil_rep_mpc2(H,Q,Y,weight,Xm,Zpb,Cfak_plus,Cfak_minus,Cv,Ds,Df)
    cdef dict Cvec=dict()
    cdef MPComplexNumber ztmp
    ztmp = CF(0)
    MS = MatrixSpace(CF,Dl,nc)
    for j in range(Qs,Qf+1):
        Cvec[j]=Matrix_complex_dense(MS,0)
        for n in range(Dl):            
            for l in range(nc):
                mpc_set(ztmp.value,Cv[j][n][l],rnd)
                Cvec[j][n,l]=ztmp
    for j from 0 <= j < Ql:
        if Cv[j]<>NULL:
            for n in range(Dl):
                if Cv[j][n]<>NULL:
                    for l in range(nc):
                        if Cv[j][n][l]<>NULL:
                            mpc_clear(Cv[j][n][l])
                    sage_free(Cv[j][n])
            sage_free(Cv[j])
    sage_free(Cv)
    return Xm,Zpb,Cvec,Cfak_minus,Cfak_plus
    
cdef pullback_pts_weil_rep_mpc2(H,int Q,RealNumber Y,RealNumber weight,Vector_real_mpfr_dense Xm,Vector_complex_dense Zpb,Vector_complex_dense Cfak_plus,Vector_complex_dense Cfak_minus, mpc_t *** Cvec,int ds=-1,int df=-1):
    r"""
    rho = representation of SL2Z, i.e. a map :SL2Z->Mat_r(C)

    -''WR'' -- Finite quadratic module or Discriminant form
    -''Q''  -- integer
    -''Y''  -- real
    -''weight''-- real
    -''ds'' -- integer
    -''df'' -- integer
    -''verbose''-- int  (default 0)

    EXAMPLES::
    sage: DF=DiscriminantForm(1,False)
    sage: Y=mpmath.mp.mpf(0.5); Q=20; k=mpmath.mp.mpf(0.5)
    sage: pullback_pts_weil_rep(DF,Q,Y,k,False)
    
    """
    cdef RealField_class RF
    cdef MPComplexNumber fak2,jfak,jfak2,twopii,mp_i,fstar,ztmp
    cdef RealNumber twopi,YY,r1,x,y,xtmp
    cdef int prec,Qs,Qf,fourQ,Ql,arg
    cdef int j,ri,rj
    cdef mpfr_t Ypt
    cdef int jj
    cdef int a,b,c,d,N
    cdef mpfr_t xx,yy,xpb,ypb
    cdef mpc_t mptmp
    cdef dict r_vals
    cdef int c2,N2,N22,Nc
    cdef int rjmax,dual_rep
    cdef tuple Tj
    cdef int verbose = H._verbose
    cdef double xt,yt
    r_vals=dict()
    #sig_on()
    X = H.multiplier()
    WM = X.weil_module()
    cdef list D = X.D()
    dual_rep = H.multiplier().is_dual()
    if hasattr(H,"_pullback_prec"):
        prec = H._pullback_prec ## Try out what happend with higher prec here
    else:
        prec = H._prec ## Try out what happend with higher prec here
    RF = RealField(prec)
    CF = MPComplexField(prec)
    mpc_init2(mptmp,prec)
    mpfr_init2(Ypt,prec)
    mpfr_init2(xx,prec)
    mpfr_init2(xpb,prec)
    mpfr_init2(ypb,prec)
    mpfr_init2(yy,prec)
    mpfr_set(Ypt,Y.value,rnd_re)
    x=RF(1); y=RF(1); r1=RF(0); YY = RF(Y)
    xtmp = RF(1); ztmp = CF(1)
    jfak = CF(1); jfak2=CF(1)
    twopi=RF(2)*RF.pi()
    twopii=CF(0,twopi)
    fstar=CF(1)
    mp_i=CF(0,1)
    #cdef Vector_real_mpfr_dense Xm
    #cdef Vector_complex_dense Zpb 
    #Xm=dict();  Xpb=dict(); Ypb=dict();
    #Xm = Vector_real_mpfr_dense(FreeModule(RF,Q+1),0)
    #Zpb= Vector_complex_dense(FreeModule(CF,Q+1),0)
    ## Cv has to be allocated outside this routine
    Qs=1; Qf=Q; Ql = Qf - Qs + 1; fourQ=4*Q
    if verbose>1:        
        print "Xm.length()=",Xm.degree()
        print "Ql=",Ql
    assert Cvec<>NULL
    #Cfak=dict()

    #N = WR.rank()/2
    #cdef int twoN = WR.weil_module().rank()
    cdef int twoN = WM.rank()
    #ds = D[0]; df=D[-1]
    #if df <0:
    #    df=twoN-1 # use full range
    #    ds=0
    #rjmax = X.rank() #-1  #len(WM.basis())-1
    #rjmax = WM.matrix(SL2Z([1,0,0,1]),numeric=1,prec=53).nrows()
    rjmax = len(X.Qv)
    N2=my_odd_part(twoN); N22=1-N2
    z=CF(0,1) #CyclotomicField(4).gen()        
    #MS = MatrixSpace(CF,int(len(WR.D)),int(len(WR.D)))
    # Recall that Az=z* in F <=> JAJ(Jz)=Jz* in F
    # so we only need to do the pullback for half of the points
    # since Xm[1-j]=-X[j]
    for j from 0 <= j < Ql: #1-Q,Q):
        #jj = 2*j-1
        jj = 2*j+2*Qs-1
        mpfr_set_si(x.value,jj,rnd_re)
        mpfr_div_ui(x.value,x.value,fourQ,rnd_re)
        xt=mpfr_get_d(x.value,rnd_re)
        yt=mpfr_get_d(Ypt,rnd_re)
        pullback_to_psl2z_mat_c(&xt,&yt,&a,&b,&c,&d)
        if verbose>2:
            mpfr_set(xtmp.value,x.value,rnd_re)
            print "Xm[",j,"]=",xtmp
        mpfr_set(xpb,x.value,rnd_re)
        mpfr_set(ypb,Ypt,rnd_re)
        _apply_sl2z_map_mpfr(xpb,ypb,a,b,c,d)
        if verbose>2:
            mpfr_set(xtmp.value,xpb,rnd_re)
            print "Xpb[",j,"]=",xtmp            
            mpfr_set(xtmp.value,ypb,rnd_re)
            print "Xpb[",j,"]=",xtmp            
        if c>=0:
            Tj=(-d,b,c,-a)
            # j_Tm(z)
            # xx = cx-a+icy; yy = cy
            mpfr_mul_si(xx,xpb,c,rnd_re)
            mpfr_sub_si(xx,xx,a,rnd_re)
            mpfr_mul_si(yy,ypb,c,rnd_re)            
            mpc_set_fr_fr(jfak.value,xx,yy,rnd)
            mpc_pow_fr(jfak.value,jfak.value,weight.value,rnd)
            #jfak=CF(cr*x-ar,cr*y)**weight
            # j_Tm*(z)
            mpfr_mul_si(xx,xpb,-c,rnd_re)
            mpfr_add_si(xx,xx,a,rnd_re)
            mpfr_mul_si(yy,ypb,c,rnd_re)            
            mpc_set_fr_fr(jfak2.value,xx,yy,rnd)
            mpc_pow_fr(jfak2.value,jfak2.value,weight.value,rnd)
            #jfak2=CF(-cr*x+ar,cr*y)**weight
            c2=my_odd_part(c)            
        else:
            Tj=(d,-b,-c,a)
            mpfr_mul_si(xx,xpb,-c,rnd_re)
            mpfr_add_si(xx,xx,a,rnd_re)
            mpfr_mul_si(yy,ypb,-c,rnd_re)            
            mpc_set_fr_fr(jfak.value,xx,yy,rnd)
            mpc_pow_fr(jfak.value,jfak.value,weight.value,rnd)
            #jfak=CF(-cr*x+ar,-cr*y)**weight
            mpfr_mul_si(xx,xpb,c,rnd_re)
            mpfr_sub_si(xx,xx,a,rnd_re)
            mpfr_mul_si(yy,ypb,-c,rnd_re)            
            mpc_set_fr_fr(jfak2.value,xx,yy,rnd)
            mpc_pow_fr(jfak2.value,jfak2.value,weight.value,rnd)
            #jfak2=CF(cr*x-ar,-cr*y)**weight
            c2=my_odd_part(-c)

        mpfr_mul(xpb,xpb,twopi.value,rnd_re)
        mpfr_mul(ypb,ypb,twopi.value,rnd_re)
        mpfr_mul(Xm._entries[j],x.value,twopi.value,rnd_re)
        mpc_set_fr_fr(Zpb._entries[j],xpb,ypb,rnd)
        if verbose>2:
            print "Xm=",Xm[j]
            print "iZpb=",Zpb[j]
        Nc=my_gcd(twoN,c)
        N2c=my_odd_part(Nc)
        if my_is_odd(c):
            arg=c2*N2
        else:
            arg=c2+N22
        arg=arg-1-N2/Nc
        if dual_rep==1:
            arg=-arg
        if Tj not in r_vals.keys():
            r_vals[Tj] = WM.matrix(SL2Z(Tj),numeric=1,prec=prec)
        rr = r_vals[Tj]
        ## do not want to large dict...
        if len(r_vals.keys())>100:
            r_vals.pop(r_vals.keys()[0])
        r=rr[0]            
        r1=rr[1]
        fak=r1*jfak
        fak2=r1*jfak2
        if verbose>2:
            print "arg=",arg
            print "jfak=",jfak
            print "rr1=",rr[1]
            print "fak=",fak
            print "fak2=",fak2
        #fstar=kronecker(-1,abs(c)/Nc)*z**arg
        fstar = CF(my_kronecker_minus_one_red(c,Nc))*z**arg
        #Cfak[j]=[fak]
        Cfak_plus[j]=CF(fak)
        fak2=fak2*fstar #.complex_embedding(prec)
        #Cfak[1-j]=[fak2,-1]
        Cfak_minus[j]=CF(fak2)
        if dual_rep==1:
            for ri in range(ds,df+1): #from ds <= ri <= df:
                for rj in range(rjmax): #)from 0 <= rj <= rjmax:                    
                    if r[ri,rj]<>0:
                        mpc_set(mptmp,(<Matrix_complex_dense>r)._matrix[ri][rj],rnd)
                        mpc_pow_si(mptmp,mptmp,-1,rnd)
                        #mpc_set((<Matrix_complex_dense>Cvec[j])._matrix[ri][rj],mptmp,rnd)
                        mpc_set(Cvec[j][ri-ds][rj],mptmp,rnd)
                        #Cvec[j][ri,rj]=r[ri,rj]**-1 #(1/r[ri,rj]).complex_embedding(prec)
                    else:
                        #mpc_set_si((<Matrix_complex_dense>Cvec[j])._matrix[ri][rj],0,rnd)
                        mpc_set_si(Cvec[j][ri-ds][rj],0,rnd)
                        #Cvec[j][ri,rj]=0 #,Cvec[j]._rnd)
        else:
            for ri in range(ds,df+1): #from ds <= ri <= df:
                for rj in range(rjmax):# from 0 <= rj <= rjmax:                    
                    mpc_set(Cvec[j][ri-ds][rj],(<Matrix_complex_dense>r)._matrix[ri][rj],rnd)
        #raise Exception
    if verbose>0:
        print "num of maps=",len(r_vals.keys())
    mpfr_clear(xx)
    mpfr_clear(yy)
    mpfr_clear(xpb)
    mpfr_clear(ypb)
    mpfr_clear(Ypt)
    mpc_clear(mptmp)
    #sig_off()
    #return [Xm,Zpb,Cfak]


cpdef get_sym_type(m,weight):
    r"""
    Calculate the symmetry type (even/odd) for the combination
    of representation and weight.
    """
    t=weight-0.5
    if not is_int(t):
        raise ValueError, "Need half-integral value of weight! Got k=%s" %(weight)
    ti=Integer(float(t))
    if is_odd(ti):
        sym_type=-1
    else:
        sym_type=1
    if m.is_dual():
        sym_type=-sym_type
    return sym_type





    
def vv_harmonic_wmwf_phase2_2(M,PP,C,Ns,Is=None,prec=None,Yin=None,do_save=False):
    r"""
    Phase 2 for vector-valued harmonic weak Maass forms.
    """
    pass
    ## WR=M.WR;
    ## kappa=M.weight
    ## D=WR.D  
    ## Dsym=M.D # the symmetrized index set
    ## if(len(Dsym)<>len(C.keys())):
    ##     raise ValueError,"Got incompatible coefficient vector! indices=%s" % C.keys()
    ## import tempfile,os
    ## #we only use symmetrized values
    ## cdef int j,n,yj,l,bi,ai,Ms,Mf
    ## if(Is==None):
    ##     Is=[0,len(D)]
        
    ## N=WR.N
    ## t=-1
    ## sym_type=mpmath.mpf(M.sym_type)
    ## verbose=M.verbose
    ## if verbose>0:
    ##     print "In Phase2"
    ## #ndig=12
    ## if(prec == None):
    ##     prec=M.prec-2
    ## eps=mpmath.power(mpmath.mpf(10),mpmath.mpf(-prec))
    ## if verbose>0:
    ##     print "eps=",eps
    ##     print "Yin=",Yin
    ## betai=dict();mbeta=dict(); mbetai=dict(); mm=dict(); mmm=dict()
    ## mptwo=mpmath.mp.mpf(2); mpfour=mpmath.mp.mpf(4)
    ## twopi=mptwo*mpmath.pi(); twopii=mpmath.mp.mpc(0,twopi)
    ## fourpi=mptwo*twopi; mp0=mpmath.mpf(0)
    ## weight=mpmath.mp.mpf(kappa); weight_over_two=weight/mpmath.mp.mpf(2)
    ## K0=0; K1=0
    ## for (beta,m) in PP:
    ##     #if( (not beta in D) or (not 1-beta in D)):
    ##     if( (not beta in D) and (not 1-beta in D)):
    ##         raise Exception,"Need beta=%s in D=%s" %(beta,D)
    ##     betai[beta]=D.index(beta)
    ##     mbeta[beta]=1-beta
    ##     mbetai[beta]=D.index(1-beta)
    ##     mm[(beta,m)]=(m+WR.Qv[betai[beta]])
    ##     mmm[(beta,m)]=mpmath.mp.mpf(mm[(beta,m)])
    ##     if verbose>0:
    ##         print "beta,m=",beta,m
    ##         print "mm=",mm[(beta,m)]
    ##     if(mm[(beta,m)]>t):
    ##         t=mm
    ##     if(abs(mm[(beta,m)])>K0):
    ##         K0=abs(mm[(beta,m)])
    ##     if(abs(PP[(beta,m)])>K1):
    ##         K1=abs(PP[(beta,m)])
    ##     # One way to specify the principal part
    ##     # is to only specify half and let the rest be decided
    ##     # by the symmetry. If we specify the rest it must match
    ##     if(PP.has_key((mbeta[beta],m)) and PP.has_key((beta,m))):
    ##         test=abs(PP[(beta,m)]-sym_type*PP[(mbeta[beta],m)])
    ##         if(test>0 and not test.ae(mp0)):
    ##             raise ValueError,"The principal part has not correct symmetry: type=%s, PP=%s" %(sym_type,PP)
    ##     else:
    ##         pass
    ## abD=len(WR.D)
   
    ## if(Yin==None):
    ##     Y0=mpmath.mp.mpf(0.5)
    ## else:
    ##     Y0=mpmath.mp.mpf(Yin)
    ## kint=mpmath.mp.mpf(1-weight)
    ## sym_type=get_sym_type(WR,weight)
    ## NA=Ns[0]; NB=Ns[1]
    ## if(sym_type==1):
    ##     Ds=int(0); Df=int(WR.N) # 0,1,...,N
    ## elif(sym_type==-1):
    ##     Ds=int(1); Df=int(WR.N-1) # 1,2,...,N-1  (since -0=0 and -N=N)
    ## IA=int(max(Is[0],Ds)); IB=int(min(Is[1],Df))
    ## Ms=int(min(C[Ds].keys())); Mf=int(max(C[Ds].keys())); Ml=int(Mf-Ms+1)
    ## #Ms=int(min(C[D[Ds]].keys())); Mf=int(max(C[D[Ds]].keys())); Ml=int(Mf-Ms+1)
    ## #tmpfile="tmpC-N"+str(WR.N)+"-"+str(M._weight_rat)+"-"+str(PP)+".sobj"
    ## #tmpfile=tmpfile.replace(" ","")
    ## [f,tmpfile] = tempfile.mkstemp(prefix='tmpCoeffs',suffix='.txt')
    ## os.close(f)
    ## #fp=open(tmpfile,"append")
    ## if verbose>0:
    ##     print "Ms,Mf,Ml=",Ms,Mf,Ml
    ## K1=K1*2*N
    ## NAA=NA; IAA=IA
    ## numys=2
    ## # have    Y=mpmath.mp.mpf(Y_in)
    ## # have to find suitable Q for the given Y
    ## if verbose>0:
    ##     print "dps=",mpmath.mp.dps
    ## ctmp=dict(); ctmp_neg=dict()
    ## Cout=dict()
    ## for bi from IA <= bi <= IB:
    ##     Cout[bi]=dict()
    ## stw=str(weight)[0:5]
    ## Qadd=0; Yfak=mpmath.mpf(0.95); Yvold=[0.0,0.0]
    ## Xm=dict();Xpb=dict();Ypb=dict(); Cv=dict()
    ## Q=dict(); Qs=dict();Qf=dict(); QQ=dict()
    ## if(do_save):
    ##     pbfile="pbfile.sobj"
    ##     [f,tmpfilename] = tempfile.mkstemp(prefix=pbfile,suffix='.sobj')
    ##     os.close(f)
    ## for yloop in range(1000):
    ##     Yv=[Y0,Yfak*Y0]
    ##     for i from 0 <= i < numys:
    ##         Q[i]=int(M.get_M(Yv[i],K0,K1,prec)+Qadd)
    ##         Qs[i]=1-Q[i]; Qf[i]=Q[i]; QQ[i]=mpmath.mp.mpf(1)/mpmath.mp.mpf(2*Q[i])
    ##         if verbose>0:
    ##             print "Yv[",i,"]=[",mppr(Yv[0]),",",mppr(Yv[1]),"]"
    ##             print "Q(Y)[",i,"]=",Q[i],type(Q[i]),type(Qs[i])
    ##             #print "1/2Q[",i,"]=",mppr(QQ[i])
    ##     # Recall that the first Y-value is always the larger
    ##     if(verbose>1):
    ##         print "Yvold=",mppr(Yvold[0]),",",mppr(Yvold[1]),"]"
    ##     pbl=dict()
    ##     if(Yv[0].ae(Yvold[1])):
    ##         if(verbose>1):
    ##             print "do not evaluate for Yv=",Yv[0]
    ##         pbl[0]=[Xm[1],Xpb[1],Ypb[1],Cv[1]]
    ##         pbl[1]=pullback_pts_weil_rep(WR,Q[1],Yv[1],weight,Ds,Df)
    ##     else:
    ##         for i from 0 <= i < numys:
    ##             pbl[i]=pullback_pts_weil_rep(WR,Q[i],Yv[i],weight,Ds,Df)
    ##     if(do_save):
    ##         tmppb=pickle.dump( (Yv,Q,pbl),tmpfilename)
    ##     for i from 0 <= i < numys:
    ##         [Xm[i],Xpb[i],Ypb[i],Cv[i]]=pbl[i]

    ##     Yvold=Yv; Zipb=dict()
    ##     for yj from 0 <= yj <numys:
    ##         Zipb[yj]=dict()
    ##         for j from 1<= j <= Qf[yj]:
    ##             Zipb[yj][j]=mpmath.mp.mpc(-Ypb[yj][j],Xpb[yj][j])
                
    ##     gamma_fak=dict()
    ##     for yj from 0 <= yj <numys:
    ##         for bi from IA <= bi <= IB:                       
    ##             for l from Ms <= l <= Mf:
    ##                 lr=mpmath.mp.mpf(l+WR.Qv[bi])
    ##                 if(lr<0):
    ##                     lrtwo=lr*mptwo
    ##                     for j from 1 <= j <= Qf[yj]:
    ##                         gamma_fak[yj,bi,l,j]=mpmath.gammainc(kint,abs(lrtwo)*Ypb[yj][j])*mpmath.mp.exp(-lr*Ypb[yj][j])
    ##                         #gamma_fak[yj,bi,l,1-j]=gamma_fak[yj,bi,l,j].conjugate()
    ##                 else:
    ##                     for j from 1 <= j <= Qf[yj]:
    ##                         gamma_fak[yj,bi,l,j]=mpmath.mp.exp(-lr*Ypb[yj][j])
    ##                         #gamma_fak[yj,bi,l,1-j]=gamma_fak[yj,bi,l,j].conjugate()
    ##     if verbose>0:
    ##         print "Got pullback points!"
    ##         print "dps=",mpmath.mp.dps
    ##         print "NAA=",NAA
    ##     # If we want to do negative coefficients too we save time by computing simultaneously
    ##     do_neg=True
    ##     try:
    ##         for n from NAA <= n <= NB:
    ##             if verbose>0:
    ##                 print "n=",n
    ##             for ai from IAA <= ai <= IB:
    ##                 if verbose>0:
    ##                     print "ai=",ai
    ##                 mai=-ai % abD
    ##                 nr=mpmath.mp.mpf(n+WR.Qv[ai])
    ##                 nrtwo=mptwo*nr
    ##                 nri=mpmath.mp.mpc(0,nr)
    ##                 if(do_neg):
    ##                     nrm=mpmath.mp.mpf(-n+WR.Qv[ai])
    ##                     nrmi=mpmath.mp.mpc(0,nrm)
    ##                     nrmtwo=mptwo*nrm
    ##                 for yj from 0 <= yj <=numys-1:
    ##                     Y=Yv[yj]*twopi
    ##                     summa=mp0;
    ##                     summa_neg=mp0
    ##                     #print "IA,IB=",IA,IB
    ##                     fak=dict()
    ##                     for j from 1 <= j <= Qf[yj]:
    ##                         fak[j]=mpmath.mp.exp(-nri*Xm[yj][j])
    ##                         #fak[1-j]=fak[j].conjugate()
    ##                     if(do_neg):
    ##                         fak_neg=dict()
    ##                         for j from 1 <= j <= Qf[yj]:
    ##                             fak_neg[j]=mpmath.mp.exp(-nrmi*Xm[yj][j])
    ##                             #fak_neg[1-j]=fak_neg[j].conjugate()
    ##                     for bi from IA<= bi <=IB:
    ##                         mbi=-bi % abD
    ##                         #print "mbi=",mbi
    ##                         for l from Ms <= l <= Mf:
    ##                             if(C[bi][l]==0 or abs(C[bi][l]).ae(mp0)):
    ##                                 if verbose>0:
    ##                                     print "Skip coeff ",bi,l
    ##                                 continue
    ##                             lr=mpmath.mp.mpf(l+WR.Qv[bi])
    ##                             ilr=mpmath.mp.mpc(0,lr)
    ##                             Vtmp=mp0;
    ##                             Vtmp_neg=mp0
    ##                             for j from 1 <= j <= Qf[yj]:
    ##                                 if(mbi<>bi):
    ##                                     ch=Cv[yj][j][ai,bi]+sym_type*Cv[yj][j][ai,mbi]
    ##                                 else:
    ##                                     ch=Cv[yj][j][ai,bi]
    ##                                 if(ch==0 or ch.ae(mp0)):
    ##                                     continue
    ##                                 if(mbi<>bi):
    ##                                     ch2=Cv[yj][1-j][ai,bi]+sym_type*Cv[yj][1-j][ai,mbi]
    ##                                 else:
    ##                                     ch2=Cv[yj][1-j][ai,bi]
    ##                                 if(ch.ae(mp0) and ch2.ae(mp0)):
    ##                                     continue
    ##                                 tmp=mpmath.exp(ilr*Xpb[yj][j])*gamma_fak[yj,bi,l,j]
    ##                                 #print "ch[",n,ai,yj,bi,j,"=",ch
    ##                                 #print "ch2[",n,ai,yj,bi,j,"=",ch2
    ##                                 #tmp2=tmp.conjugate()
    ##                                 #(ch2*mpmath.exp(ilr*Xpb[yj][j]))*gamma_fak[yj,bi,l,j]
    ##                                 Vtmp=Vtmp+tmp*ch*fak[j]+ch2*(fak[j]*tmp).conjugate()
    ##                                 if(do_neg):
    ##                                     Vtmp_neg=Vtmp_neg+(ch*tmp*fak_neg[j]+ch2*(fak_neg[j]*tmp).conjugate())
    ##                             summa=summa+Vtmp*C[bi][l]
    ##                             if(do_neg):
    ##                                 summa_neg=summa_neg+Vtmp_neg*C[bi][l]
    ##                     if(verbose>1):
    ##                         print "summa(",yj,")=",summa
    ##                         if(do_neg):
    ##                             print "summa_neg(",yj,")=",summa_neg
    ##                     wsumma=mp0;
    ##                     wsumma_neg=mp0
    ##                     for (beta,m) in PP:
    ##                         app=mpmath.mp.mpf(PP[(beta,m)])
    ##                         lr=mpmath.mp.mpf(m+WR.Qv[betai[beta]])
    ##                         tmpsumma=mp0;
    ##                         tmpsumma_neg=mp0
    ##                         for j from 1 <= j <= Qf[yj]:
    ##                             if(betai[beta] <> mbetai[beta]):
    ##                                 ch=Cv[yj][j][ai,betai[beta]]+sym_type*Cv[yj][j][ai,mbetai[beta]]
    ##                             else:
    ##                                 ch=Cv[yj][j][ai,betai[beta]]
    ##                             if(betai[beta] <> mbetai[beta]):
    ##                                 ch2=Cv[yj][1-j][ai,betai[beta]]+sym_type*Cv[yj][1-j][ai,mbetai[beta]]
    ##                             else:
    ##                                 ch2=Cv[yj][1-j][ai,betai[beta]]
    ##                             if((ch.ae(mp0) and ch2.ae(mp0))):
    ##                                 continue
    ##                             # print "ch[",n,ai,yj,bi,j,"=",ch
    ##                             # print "ch2[",n,ai,yj,bi,j,"=",ch2
    ##                             tmp=mpmath.exp(lr*Zipb[yj][j])
    ##                             tmpsumma=tmpsumma+(ch*tmp*fak[j]+ch2*(tmp*fak[j]).conjugate())
    ##                             if(do_neg):
    ##                                 tmpsumma_neg=tmpsumma_neg+(ch*tmp*fak_neg[j]+ch2*(tmp*fak_neg[j]).conjugate())
    ##                         wsumma=wsumma+app*tmpsumma
    ##                         if(do_neg):
    ##                             wsumma_neg=wsumma_neg+app*tmpsumma_neg
    ##                     if(verbose>1):
    ##                         print "wsumma(",yj,")=",wsumma
    ##                         if(do_neg):
    ##                             print "wsumma_neg(",yj,")=",wsumma_neg
    ##                     sumtmp=(summa+wsumma)*QQ[yj]
    ##                     if(PP.has_key((D[ai],n))>0):
    ##                         sum_tmp=sum_tmp-PP[(D[ai],n)]*mpmath.mp.exp(-nr*Y)
    ##                     lhs=mpmath.mp.exp(nr*Y)
    ##                     if verbose>0:
    ##                         print "exp(2pinY)=",mppr(lhs)
    ##                     ctmp[yj]=sumtmp*lhs
    ##                     #return
    ##                     if(do_neg):
    ##                         sumtmp_neg=(summa_neg+wsumma_neg)*QQ[yj]
    ##                         if(PP.has_key((D[ai],-n))>0):
    ##                             sumtmp_neg=sumtmp_neg-PP[(D[ai],-n)]*mpmath.mp.exp(-nrm*Y)
    ##                         lhs=mpmath.gammainc(kint,abs(nrmtwo)*Y)*mpmath.mp.exp(-nrm*Y)
    ##                         if verbose>0:
    ##                             print "Gamma(1-k,4pinY)=",mppr(lhs)
    ##                         ctmp_neg[yj]=sumtmp_neg/lhs
    ##                 # end for yj
    ##                 if(verbose>-1):
    ##                     print "C1[",n,ai,"]=",ctmp[0].real,'+i*',mppr(ctmp[0].imag)
    ##                     print "C2[",n,ai,"]=",ctmp[1].real,'+i*',mppr(ctmp[1].imag)
    ##                     if(do_neg):
    ##                         print "C1[",-n,ai,"]=",ctmp_neg[0].real,'+i*',mppr(ctmp_neg[0].imag)
    ##                         print "C2[",-n,ai,"]=",ctmp_neg[1].real,'+i*',mppr(ctmp_neg[1].imag)
    ##                 if(do_neg):
    ##                     err_pos=abs(ctmp[1]-ctmp[0])
    ##                     err_neg=abs(ctmp_neg[1]-ctmp_neg[0])
    ##                     err=err_pos # max(err_pos,err_neg)
    ##                     if(verbose>-1):
    ##                         print "err_pos=",mppr(err_pos)
    ##                         print "err_neg=",mppr(err_neg)
    ##                 else:
    ##                     err=abs(ctmp[1]-ctmp[0])
    ##                     if(verbose>-1):
    ##                         print "err=",mppr(err)
    ##                 if verbose>0:                        
    ##                     if(C.keys().count(ai)>0):
    ##                         if(C[ai].keys().count(n)>0):
    ##                             print "Cin(",ai,n,")=",C[ai][n].real
    ##                             if(C[ai].keys().count(-n)>0):
    ##                                 print "Cin(",ai,-n,")=",C[ai][-n].real
    ##                 #sys.stdout.flush()
    ##                 if(err>eps):
    ##                     # Have to modify
    ##                     Y0=Yv[0]*Yfak
    ##                     if(verbose>-1):                        
    ##                         print " Need to decrease Y! new Y0=",Y0
    ##                     #Qadd=Qadd+10
    ##                     #Yv[0]=Y0; Yv[1]=Y0*mpmath.mp.mpf(0.95)
    ##                     NAA=n; IAA=ai
    ##                     raise StopIteration()
    ##                 else:
    ##                     Cout[ai][n]=ctmp[1]
    ##                     if(do_neg):
    ##                         Cout[ai][-n]=ctmp_neg[1]
    ##                     if verbose>0:                        
    ##                         print "OK! av=",(ctmp[1]+ctmp[0])/mpmath.mpf(2)
    ##                     #save(Cout,tmpfile)
    ##                     fp=open(tmpfile,"append")
    ##                     s="C["+str(ai)+"]["+str(n)+"]="+str(Cout[ai][n])+"\n"
    ##                     s=s+"C["+str(ai)+"]["+str(-n)+"]="+str(Cout[ai][-n])+"\n"
    ##                     fp.write(s)
    ##                     fp.close()
    ##                     # If we are in the range of the used C's we update
    ##                     #if(C.keys().count(ai)>0):
    ##                     #    if(C[ai].keys().count(n)>0):
    ##                     #        C[ai][n]=ctmp[1]                                
    ##                     #    if(do_neg and C[ai].keys().count(-n)>0):
    ##                     #        C[ai][-n]=ctmp_neg[1]                                
    ##                     #continue
    ##             # end for ai
    ##         #print "n at end=",n
    ##         # end for n
    ##         return Cout
    ##     except StopIteration:
    ##         if verbose>0:
    ##             print "Iteration stopped!"
    ##         continue
    ## raise StopIteration()

### Remember to uncomment this routine later!!!

def vv_harmonic_wmwf_phase2_2_ef(F,Ns,Is=None,prec=20,Yin=None,do_save=False,Qadd_in=None):
    r"""
    Phase 2 for vector-valued harmonic weak Maass forms.
    """
    pass
    ## C=F.coeffs
    ## PP=F.principal_part
    ## M=F.space
    ## WR=M.WR; kappa=M.weight
    ## D=WR.D  
    ## Dsym=M.D_as_int # the symmetrized index set
    ## if(len(Dsym)<>len(C.keys())):
    ##     raise ValueError,"Got incompatible coefficient vector! indices=%s" % C.keys()
    ## import tempfile,os
    ## #we only use symmetrized values
    ## cdef int j,n,yj,l,bi,ai,Ms,Mf
    ## if(Is==None):
    ##     Is=[0,len(D)]
        
    ## N=WR.N
    ## t=-1
    ## sym_type=mpmath.mpf(M.sym_type)
    ## verbose=M.verbose
    ## #ndig=12
    ## eps=mpmath.power(mpmath.mpf(10),mpmath.mpf(-prec))
    ## if verbose>0:
    ##     print "eps=",eps
    ##     print "Yin=",Yin
    ##     print "Is=",Is
    ##     print "Ns=",Ns
    ## betai=dict();mbeta=dict(); mbetai=dict(); mm=dict(); mmm=dict()
    ## mptwo=mpmath.mp.mpf(2); mpfour=mpmath.mp.mpf(4)
    ## twopi=mptwo*mpmath.pi(); twopii=mpmath.mp.mpc(0,twopi)
    ## fourpi=mptwo*twopi; mp0=mpmath.mpf(0)
    ## weight=mpmath.mp.mpf(kappa); weight_over_two=weight/mpmath.mp.mpf(2)
    ## K0=0; K1=0
    ## for (beta,m) in PP:
    ##     #if( (not beta in D) or (not 1-beta in D)):
    ##     if( (not beta in D) and (not 1-beta in D)):
    ##         raise Exception,"Need beta=%s in D=%s" %(beta,D)
    ##     betai[beta]=D.index(beta)
    ##     mbeta[beta]=1-beta
    ##     mbetai[beta]=D.index(1-beta)
    ##     mm[(beta,m)]=(m+WR.Qv[betai[beta]])
    ##     mmm[(beta,m)]=mpmath.mp.mpf(mm[(beta,m)])
    ##     if verbose>0:
    ##         print "beta,m=",beta,m
    ##         print "mm=",mm[(beta,m)]
    ##         #print "-beta=",minus_beta
    ##         #print "mm=",mm[(beta,m)]
    ##     if(mm[(beta,m)]>t):
    ##         t=mm
    ##     if(abs(mm[(beta,m)])>K0):
    ##         K0=abs(mm[(beta,m)])
    ##     if(abs(PP[(beta,m)])>K1):
    ##         K1=abs(PP[(beta,m)])
    ##     # One way to specify the principal part
    ##     # is to only specify half and let the rest be decided
    ##     # by the symmetry. If we specify the rest it must match
    ##     if(PP.has_key((mbeta[beta],m)) and PP.has_key((beta,m))):
    ##         test=abs(PP[(beta,m)]-sym_type*PP[(mbeta[beta],m)])
    ##         if(test>0 and not test.ae(mp0)):
    ##             raise ValueError,"The principal part has not correct symmetry: type=%s, PP=%s" %(sym_type,PP)
    ##     else:
    ##         pass
    ## abD=len(WR.D)
    ## if(Yin==None):
    ##     Y0=mpmath.mp.mpf(0.5)
    ## else:
    ##     Y0=mpmath.mp.mpf(Yin)
    ## kint=mpmath.mp.mpf(1-weight)
    ## sym_type=get_sym_type(WR,weight)
    ## NA=Ns[0]; NB=Ns[1]
    ## if(sym_type==1):
    ##     Ds=int(0); Df=int(WR.N) # 0,1,...,N
    ## elif(sym_type==-1):
    ##     Ds=int(1); Df=int(WR.N-1) # 1,2,...,N-1  (since -0=0 and -N=N)
    ## IA=int(max(Is[0],Ds)); IB=int(min(Is[1],Df))
    ## Ms=int(min(C[Ds].keys())); Mf=int(max(C[Ds].keys())); Ml=int(Mf-Ms+1)
    ## #Ms=int(min(C[D[Ds]].keys())); Mf=int(max(C[D[Ds]].keys())); Ml=int(Mf-Ms+1)
    ## tmpfile="tmpC-N"+str(WR.N)+"-"+str(M._weight_rat)+"-"+str(PP)+".sobj"
    ## [f,tmpfile] = tempfile.mkstemp(prefix='tmpCoeffs',suffix='.txt')
    ## os.close(f)
    ## if verbose>0:
    ##     print "IA,IB=",IA,IB
    ##     print "Ds,Df=",Ds,Df
    ##     print "Ms,Mf,Ml=",Ms,Mf,Ml
    ## K1=K1*2*N
    ## NAA=NA; IAA=IA
    ## numys=2
    ## # have    Y=mpmath.mp.mpf(Y_in)
    ## # have to find suitable Q for the given Y
    ## if verbose>0:
    ##     print "dps=",mpmath.mp.dps
    ## ctmp=dict(); ctmp_neg=dict()
    ## Cout=dict()
    ## for bi from IA <= bi <= IB:
    ##     Cout[bi]=dict()
    ## stw=str(weight)[0:5]
    ## if(Qadd_in<>None):
    ##     Qadd=Qadd_in
    ## else:
    ##     Qadd=0

    ## if(verbose>1):
    ##     print "Qadd=",Qadd
    ## Yfak=mpmath.mpf(0.95); Yvold=[0.0,0.0]
    ## Xm=dict();Xpb=dict();Ypb=dict(); Cv=dict(); Cfak=dict()
    ## Q=dict(); Qs=dict();Qf=dict(); QQ=dict()
    ## sys.stdout.flush()

    ## if(do_save):
    ##     pbfile="pbfile.sobj"
    ##     [f,tmpfilename] = tempfile.mkstemp(prefix=pbfile,suffix='.sobj')
    ##     os.close(f)
    ## for yloop in range(10):
    ##     Yv=[Y0,Yfak*Y0]
    ##     for i from 0 <= i < numys:
    ##         Q[i]=max(int(M.get_M(Yv[i],K0,K1,prec)),NAA+1)+Qadd
    ##         Qs[i]=1-Q[i]; Qf[i]=Q[i]; QQ[i]=mpmath.mp.mpf(1)/mpmath.mp.mpf(2*Q[i])
    ##         if verbose>0:
    ##             print "Yv[",i,"]=[",mppr(Yv[0]),",",mppr(Yv[1]),"]"
    ##             print "Q(Y)[",i,"]=",Q[i]
    ##             #print "1/2Q[",i,"]=",mppr(QQ[i])
    ##     # Recall that the first Y-value is always the larger
    ##     if(verbose>1):
    ##         print "Yvold=",mppr(Yvold[0]),",",mppr(Yvold[1]),"]"
    ##     pbl=dict()
    ##     if(Yv[0].ae(Yvold[1])):
    ##         if(verbose>1):
    ##             print "do not evaluate for Yv=",Yv[0]
    ##             sys.stdout.flush()
                    
    ##         pbl[0]=[Xm[1],Xpb[1],Ypb[1],Cv[1],Cfak[1]]
    ##         pbl[1]=pullback_pts_weil_rep_ef(WR,Q[1],Yv[1],weight,Ds,Df)
    ##     else:
    ##         for i from 0 <= i < numys:
    ##             pbl[i]=pullback_pts_weil_rep_ef(WR,Q[i],Yv[i],weight,Ds,Df)
    ##     #if(do_save):
    ##     #    tmppb=pickle.dump( (Yv,Q,pbl),tmpfilename)
    ##     for i from 0 <= i < numys:
    ##         [Xm[i],Xpb[i],Ypb[i],Cv[i],Cfak[i]]=pbl[i]

    ##     Yvold=Yv; Zipb=dict()
    ##     for yj from 0 <= yj <numys:
    ##         Zipb[yj]=dict()
    ##         for j from 1<= j <= Qf[yj]:
    ##             Zipb[yj][j]=mpmath.mp.mpc(-Ypb[yj][j],Xpb[yj][j])
                
    ##     gamma_fak=dict()
    ##     for yj from 0 <= yj <numys:
    ##         for bi in Dsym:                       
    ##             for l from Ms <= l <= Mf:
    ##                 lr=mpmath.mp.mpf(l+WR.Qv[bi])
    ##                 if(lr<0):
    ##                     lrtwo=lr*mptwo
    ##                     for j from 1 <= j <= Qf[yj]:
    ##                         gamma_fak[yj,bi,l,j]=mpmath.gammainc(kint,abs(lrtwo)*Ypb[yj][j])*mpmath.mp.exp(-lr*Ypb[yj][j])
    ##                         #gamma_fak[yj,bi,l,1-j]=gamma_fak[yj,bi,l,j].conjugate()
    ##                 else:
    ##                     for j from 1 <= j <= Qf[yj]:
    ##                         gamma_fak[yj,bi,l,j]=mpmath.mp.exp(-lr*Ypb[yj][j])
    ##                         #gamma_fak[yj,bi,l,1-j]=gamma_fak[yj,bi,l,j].conjugate()
    ##     if verbose>0:
    ##         print "Got pullback points!"
    ##         print "dps=",mpmath.mp.dps
    ##         print "NAA=",NAA
    ##         sys.stdout.flush()
                
    ##     # If we want to do negative coefficients too we save time by computing simultaneously
    ##     do_neg=True
    ##     try:
    ##         for n from NAA <= n <= NB:
    ##             if verbose>0:
    ##                 print "n=",n
    ##             for ai from IAA <= ai <= IB:
    ##                 if verbose>0:
    ##                     print "ai=",ai
    ##                 mai=-ai % abD
    ##                 nr=mpmath.mp.mpf(n+WR.Qv[ai])
    ##                 nrtwo=mptwo*nr
    ##                 nri=mpmath.mp.mpc(0,nr)
    ##                 if(do_neg):
    ##                     nrm=mpmath.mp.mpf(WR.Qv[ai]-n)
    ##                     nrmtwo=mptwo*nrm
    ##                     nrmi=mpmath.mp.mpc(0,nrm)
    ##                 for yj from 0 <= yj <=numys-1:
    ##                     Y=Yv[yj]*twopi
    ##                     summa=0
    ##                     summa_neg=0
    ##                     #print "IA,IB=",IA,IB
    ##                     fak=dict()
    ##                     for j from 1 <= j <= Qf[yj]:
    ##                         fak[j]=mpmath.mp.exp(-nri*Xm[yj][j])
    ##                         #fak[1-j]=fak[j].conjugate()
    ##                     if(do_neg):
    ##                         fak_neg=dict()
    ##                         for j from 1 <= j <= Qf[yj]:
    ##                             fak_neg[j]=mpmath.mp.exp(-nrmi*Xm[yj][j])
    ##                             #fak_neg[1-j]=fak_neg[j].conjugate()
    ##                     for bi in Dsym:
    ##                         mbi=-bi % abD
    ##                         #print "mbi=",mbi
    ##                         for l from Ms <= l <= Mf:
    ##                             if(C[bi][l]==0 or abs(C[bi][l]).ae(mp0)):
    ##                                 if verbose>0:
    ##                                     print "Skip coeff ",bi,l
    ##                                 continue
    ##                             lr=mpmath.mp.mpf(l+WR.Qv[bi])
    ##                             ilr=mpmath.mp.mpc(0,lr)
    ##                             Vtmp=0
    ##                             Vtmp_neg=0
    ##                             for j from 1 <= j <= Qf[yj]:
    ##                                 if(mbi==bi):
    ##                                     ch=Cv[yj][j][ai,bi]
    ##                                     if(Cv[yj][j][ai,mbi]<>0):
    ##                                         ch2=1/Cv[yj][j][ai,mbi]
    ##                                 else:
    ##                                     ch=(Cv[yj][j][ai,bi]+sym_type*Cv[yj][j][ai,mbi])
    ##                                     ch21=0; ch22=0
    ##                                     if(Cv[yj][j][ai,bi]<>0):
    ##                                         ch21=1/Cv[yj][j][ai,bi]
    ##                                     if(Cv[yj][j][ai,mbi]<>0):
    ##                                         ch22=1/Cv[yj][j][ai,mbi]
    ##                                     ch2=(ch22+sym_type*ch21)
    ##                                 if(ch==0 and ch2==0):
    ##                                     continue
    ##                                 tmp=mpmath.exp(ilr*Xpb[yj][j])*gamma_fak[yj,bi,l,j]
    ##                                 tmp1=tmp*fak[j]
    ##                                 ch=ch*Cfak[yj][j][0]; ch2=ch2*Cfak[yj][1-j][0]
    ##                                 Vtmp=Vtmp+ch*tmp1+ch2*tmp1.conjugate()
    ##                                 if(do_neg):
    ##                                     tmp2=tmp*mpmath.mp.exp(-nrmi*Xm[yj][j])
    ##                                     #tmp2=tmp*fak_neg[j]
    ##                                     Vtmp_neg=Vtmp_neg+ch*tmp2+ch2*tmp2.conjugate()
    ##                             summa=summa+Vtmp*C[bi][l]
    ##                             if(do_neg):
    ##                                 summa_neg=summa_neg+Vtmp_neg*C[bi][l]
    ##                     if(verbose>1):
    ##                         print "summa(",yj,")=",summa
    ##                         if(do_neg):
    ##                             print "summa_neg(",yj,")=",summa_neg
    ##                     wsumma=0
    ##                     wsumma_neg=0
    ##                     for (beta,m) in PP:
    ##                         app=mpmath.mp.mpf(PP[(beta,m)])
    ##                         lr=mpmath.mp.mpf(m+WR.Qv[betai[beta]])
    ##                         tmpsumma=0
    ##                         tmpsumma_neg=0
    ##                         bi=betai[beta]
    ##                         mbi=mbetai[beta]
    ##                         for j from 1 <= j <= Qf[yj]:
    ##                             ch=0; ch2=0
    ##                             if(bi == mbi):
    ##                                 ch=Cv[yj][j][ai,bi]
    ##                                 if(Cv[yj][j][ai,mbi]<>0): ch2=1/Cv[yj][j][ai,mbi]
    ##                             else:
    ##                                 ch=Cv[yj][j][ai,bi]+sym_type*Cv[yj][j][ai,mbi]
    ##                                 ch21=0; ch22=0
    ##                                 if(Cv[yj][j][ai,bi]<>0):  ch21=1/Cv[yj][j][ai,bi]
    ##                                 if(Cv[yj][j][ai,mbi]<>0): ch22=1/Cv[yj][j][ai,mbi]
    ##                                 ch2=ch22+sym_type*ch21
    ##                             if(ch==0 and ch2==0):
    ##                                 continue
    ##                             tmp=mpmath.exp(lr*Zipb[yj][j])
    ##                             tmp1=tmp*fak[j]
    ##                             ch=ch*Cfak[yj][j][0]; ch2=ch2*Cfak[yj][1-j][0]
    ##                             tmpsumma=tmpsumma+ch*tmp1+ch2*tmp1nn.conjugate()
    ##                             if(do_neg):
    ##                                 tmp2=mpmath.exp(lr*Zipb[yj][j]-nrmi*Xm[yj][j])
    ##                                 #tmp2=tmp*fak_neg[j]
    ##                                 tmpsumma_neg=tmpsumma_neg+ch*tmp2+ch2*tmp2.conjugate()
    ##                         wsumma=wsumma+app*tmpsumma
    ##                         if(do_neg):
    ##                             wsumma_neg=wsumma_neg+app*tmpsumma_neg
    ##                     if(verbose>1):
    ##                         print "wsumma(",yj,")=",wsumma
    ##                         if(do_neg):
    ##                             print "wsumma_neg(",yj,")=",wsumma_neg
    ##                     sumtmp=(summa+wsumma)*QQ[yj]
    ##                     if(PP.has_key((D[ai],n))>0):
    ##                         tmp=PP[(D[ai],n)]*mpmath.mp.exp(-nr*Y)
    ##                         if(verbose>1):                                
    ##                             print "subtracting:",tmp
    ##                         sum_tmp=sum_tmp-tmp
    ##                     lhs=mpmath.mp.exp(nr*Y)
    ##                     if verbose>0:
    ##                         print "exp(2pinY)=",mppr(lhs)
    ##                     ctmp[yj]=sumtmp*lhs
    ##                     #return 
    ##                     if(do_neg):
    ##                         sumtmp_neg=(summa_neg+wsumma_neg)*QQ[yj]
    ##                         if(verbose>1):
    ##                             print "summa_neg+wsumma_neg=",summa_neg+wsumma_neg
    ##                             print " / 2Q=",sumtmp_neg
    ##                         if(PP.has_key((D[ai],-n))>0):
    ##                             tmp=PP[(D[ai],-n)]*mpmath.mp.exp(-nrm*Y)
    ##                             if(verbose>1):                                
    ##                                 print "subtracting:",tmp
    ##                                 sumtmp_neg=sumtmp_neg-tmp
    ##                         lhs=mpmath.mp.gammainc(kint,abs(nrmtwo)*Y)*mpmath.mp.exp(-nrm*Y)
    ##                         if verbose>0:
    ##                             print "exp(-nrmY)=",mpmath.mp.exp(-nrm*Y)
    ##                             print "Gamma(1-k,4pinY)=",mppr(lhs)
    ##                         ctmp_neg[yj]=sumtmp_neg/lhs
    ##                 # end for yj
    ##                 if(verbose>-1):
    ##                     print "C1[",n,ai,"]=",ctmp[0].real,'+i*',mppr(ctmp[0].imag)
    ##                     print "C2[",n,ai,"]=",ctmp[1].real,'+i*',mppr(ctmp[1].imag)
    ##                     if(do_neg):
    ##                         print "C1[",-n,ai,"]=",ctmp_neg[0].real,'+i*',mppr(ctmp_neg[0].imag)
    ##                         print "C2[",-n,ai,"]=",ctmp_neg[1].real,'+i*',mppr(ctmp_neg[1].imag)
    ##                 if(do_neg):
    ##                     err_pos=abs(ctmp[1]-ctmp[0])
    ##                     err_neg=abs(ctmp_neg[1]-ctmp_neg[0])
    ##                     err=err_pos # max(err_pos,err_neg)
    ##                     if(verbose>-1):
    ##                         print "err_pos=",mppr(err_pos)
    ##                         print "err_neg=",mppr(err_neg)
    ##                 else:
    ##                     err=abs(ctmp[1]-ctmp[0])
    ##                     if(verbose>-1):
    ##                         print "err=",mppr(err)
    ##                 if verbose>0:                        
    ##                     if(C.keys().count(ai)>0):
    ##                         if(C[ai].keys().count(n)>0):
    ##                             print "Cin(",ai,n,")=",C[ai][n].real
    ##                             if(C[ai].keys().count(-n)>0):
    ##                                 print "Cin(",ai,-n,")=",C[ai][-n].real
    ##                 sys.stdout.flush()
    ##                 if(err>eps):
    ##                     # Have to modify
    ##                     Y0=Yv[0]*Yfak
    ##                     if(verbose>-1):                        
    ##                         print " Need to decrease Y! new Y0=",Y0
    ##                     #Qadd=Qadd+10
    ##                     #Yv[0]=Y0; Yv[1]=Y0*mpmath.mp.mpf(0.95)
    ##                     NAA=n; IAA=ai
    ##                     raise StopIteration()
    ##                 else:
    ##                     Cout[ai][n]=ctmp[1]
    ##                     if(do_neg):
    ##                         Cout[ai][-n]=ctmp_neg[1]
    ##                     if verbose>0:                        
    ##                         print "OK! av=",(ctmp[1]+ctmp[0])/mpmath.mpf(2)
    ##                     #save(Cout,tmpfile)
    ##                     fp=open(tmpfile,"append")
    ##                     s="C["+str(ai)+"]["+str(n)+"]="+str(Cout[ai][n])+"\n"
    ##                     s=s+"C["+str(ai)+"]["+str(-n)+"]="+str(Cout[ai][-n])+"\n"
    ##                     fp.write(s)
    ##                     fp.close()
    ##                     # If we are in the range of the used C's we update
    ##                     #if(C.keys().count(ai)>0):
    ##                     #    if(C[ai].keys().count(n)>0):
    ##                     #        C[ai][n]=ctmp[1]                                
    ##                     #    if(do_neg and C[ai].keys().count(-n)>0):
    ##                     #        C[ai][-n]=ctmp_neg[1]                                
    ##                     #continue
    ##             # end for ai
    ##         #print "n at end=",n
    ##         # end for n
    ##         return Cout
    ##     except StopIteration:
    ##         if verbose>0:
    ##             print "Iteration stopped!"
    ##         continue
    ## raise StopIteration()



def mppr(x):
    r"""
    Print fewer digits!

    """
    
    # We don't want to print 500 digits
    # So I truncate at two digits with a stupid simple rounding
    dpold=mpmath.mp.dps
    mpmath.mp.dps=5
    s=str(x)
    mpmath.mp.dps=dpold
    (s1,dot,s2)=s.partition(".")
    (s3,es,s4)=s2.partition("e")
    if len(s4)==0:  # No "e" format
        return s[0:15]
    try:
        ma=s3[0]
        if(len(s3)>1):
            if(s3[2]<5):
                ma=ma+s3[1]
            elif(s3[2]>=5):
                ma=ma+str(Integer(s3[1])+1)
        ss=s1+"."+ma+"e"+s4
    except:
        return s
    return ss



cpdef long my_odd_part(long a):
    return _odd_part(a)

cdef long _odd_part(long a):
    cdef mpz_t az[1],bz[1],two[1]
    cdef long b
    mpz_init(bz[0])
    mpz_init_set_si(az[0],a)
    mpz_init_set_ui(two[0],2)
    mpz_remove(bz[0],az[0],two[0])
    b = mpz_get_si(bz[0])
    mpz_clear(az[0])
    mpz_clear(bz[0])
    mpz_clear(two[0])
    return b

cpdef long my_gcd(long a,long b):
    return _gcd(a,b)

cdef long _gcd(long a,long b):
    cdef mpz_t az[1],bz[1],gcd[1]
    cdef long igcd
    mpz_init(gcd[0])
    mpz_init_set_si(bz[0],b)
    mpz_init_set_si(az[0],a)
    mpz_gcd (gcd[0],az[0],bz[0])
    igcd = mpz_get_ui(gcd[0])
    mpz_clear(az[0])
    mpz_clear(bz[0])
    mpz_clear(gcd[0])
    return igcd

cpdef int my_is_odd(long a):
    cdef mpz_t az[1]
    cdef int res
    mpz_init_set_si(az[0],a)
    res = mpz_odd_p (az[0])
    mpz_clear(az[0])
    return res

cpdef int my_kronecker_minus_one_red(long c,long Nc):
    r"""
    Nc must be positive and divide c otherwise the result will be wrong!
    We do not check this here.
    """
    cdef mpz_t cz[1]
    cdef int res
    mpz_init_set_si(cz[0],c)
    mpz_tdiv_q_ui(cz[0],cz[0],Nc)
    res = mpz_si_kronecker(-1,cz[0])
    mpz_clear(cz[0])
    return res
    

### Copied from sage/libs/mpmath/utils.pyx
cdef mpfr_from_mpfval(mpfr_t res, tuple x):
    """
    Set value of an MPFR number (in place) to that of a given mpmath mpf
    data tuple.
    """
    cdef int sign
    cdef cInteger man
    cdef long exp
    cdef long bc
    sign, man, exp, bc = x
    if man.__nonzero__():
        mpfr_set_z(res, man.value, MPFR_RNDZ)
        if sign:
            mpfr_neg(res, res, MPFR_RNDZ)
        mpfr_mul_2si(res, res, exp, MPFR_RNDZ)
        return
    from mpmath.libmp import finf, fninf
    if exp == 0:
        mpfr_set_ui(res, 0, MPFR_RNDZ)
    elif x == finf:
        mpfr_set_inf(res, 1)
    elif x == fninf:
        mpfr_set_inf(res, -1)
    else:
        mpfr_set_nan(res)

cdef mpfr_to_mpfval(mpfr_t value):
    """
    Given an MPFR value, return an mpmath mpf data tuple representing
    the same number.
    """
    if mpfr_nan_p(value):
        from mpmath.libmp import fnan
        return fnan
    if mpfr_inf_p(value):
        from mpmath.libmp import finf, fninf
        if mpfr_sgn(value) > 0:
            return finf
        else:
            return fninf
    if mpfr_sgn(value) == 0:
        from mpmath.libmp import fzero
        return fzero
    sign = 0
    cdef cInteger man = PY_NEW(Integer)
    exp = mpfr_get_z_exp(man.value, value)
    if mpz_sgn(man.value) < 0:
        mpz_neg(man.value, man.value)
        sign = 1
    cdef unsigned long trailing
    trailing = mpz_scan1(man.value, 0)
    if trailing:
        mpz_tdiv_q_2exp(man.value, man.value, trailing)
        exp += trailing
    bc = mpz_sizeinbase(man.value, 2)
    return (sign, man, int(exp), bc)


cdef int is_zero(mpc_t x):
    r"""
    returns 1 if x is zero and 0 if x is non-zero
    """
    cdef int t
    t = mpfr_zero_p(mpc_realref(x))
    if t==0:
        return 0
    t = mpfr_zero_p(mpc_imagref(x))
    if t==0:
        return 0
    return 1

##
#@cached_function
cpdef rn_from_D(m,D,verbose=0):
    r""" Find the pair(s) (r,n) s.t. +-D/4N= n +- q(r) for D in D 
         set r=-1 if D is not a square mod 4N
    
    INPUT:
    -''D'' -- integer or list of integers

    OUTPUT:
    -''t'' -- tuple (r,n) or list of tuples

    """
    cdef int N, l,sig = 1
    if m._dual==1: sig = -1
    N = int(m.weil_module().level())
    cdef double Df
    cdef list lout
    cdef tuple t
    cdef double* Qv
    cdef int r,n
    l = m.rank()
    if verbose>0:
        print "D=",D
    Qv = <double*>sage_malloc(sizeof(double)*l)
    for j in range(l):
        Qv[j]=float(m.Qv[j])
    if isinstance(D,list):
        lout=list()
        for DD in D: 
            Df = float(DD)/float(N)
            one_rn_from_D(N,sig,l,Df,Qv,&r,&n,verbose)
            if r>-1:
                lout.append((r,n))
        sage_free(Qv)
        return lout
    else:
        Df = float(D)/float(N)
        one_rn_from_D(N,sig,l,Df,Qv,&r,&n,verbose)
        sage_free(Qv)
        return r,n

cdef one_rn_from_D(int N,int sig,int l,double Df,double* Qv,int* r,int* n,int verbose=0):
    r""" Find the (r,n) s.t. +-D/4N= n +- q(r)
    return r=-1 if D is not a square mod 4N
    """            
    #Dv=QQ(D)/QQ(WR.level())
    #sig=1
    #if WR._is_dual_rep:
    #    sig=-1
    cdef double x
    cdef int i = 0
    n[0] = 0; r[0] = -1
    if verbose>0:
        print "Df=",Df
    for i in range(l):
        t = Df-Qv[i]
        if verbose>0:
            print "Df-Qv[{0}]={1}".format(i,t)
        if floor(t)==ceil(t):
            n[0]=sig*int(t)
            r[0] = i
            break




#def vv_harmonic_wmwf_phase2_tst1(M,PP,C,Ns,Is=None,prec=20,Yin=None):
#    try:
#        CC=vv_harmonic_wmwf_phase2_1(M,PP,C,Ns,Is,prec,Yin)
#        return CC
#    except KeyboardInterrupt:
#        pass
