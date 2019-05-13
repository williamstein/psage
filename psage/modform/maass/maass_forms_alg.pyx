# cython: profile=False
# cython: linetrace=False
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
Cython algorithms for Maass waveforms.
Used by routines in maass_forms.py

"""
from psage.rings.mp_cimports cimport *
from cysignals.memory cimport sig_free,sig_malloc
from cysignals.signals cimport sig_on,sig_off
from sage.all import save,load,bessel_K
#from sage.functions.gamma import gamma_inc as incomplete_gamma
import mpmath

from sage.rings.complex_mpc import MPComplexField
from sage.rings.real_mpfr import RealField
import sage.structure.element
cimport sage.structure.element
from sage.structure.element cimport Element, ModuleElement, RingElement
from sage.modular.cusps import Cusp
from sage.rings.infinity import infinity
from sage.rings.integer import Integer,is_Integer
from sage.rings.complex_double import CDF
from sage.all import log_b,RR
from numpy import array
from sage.all import exp,I,CC
from sage.all import find_root
## Since I want to import sinand cos from math.h
cimport numpy as cnp
import numpy as np
#cimport cython
import scipy.special
#from scipy.special import kv



import cython
cdef extern from "math.h":
    double fabs(double)
    double fmax(double,double)
    int ceil(double)
    double sqrt(double)
    double sin(double)
    double cos(double)
    double log(double)
    double erfc(double)
    double power(double,double)
    double M_PI

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


cdef int gcd( int a, int b ):
    cdef int c
    while ( a != 0 ):
        c = a; a = b%a;  b = c
    return b

cdef double complex cexpi(double x):
    return cexp(x*_I)
#ctypedef void (*foo_t)(int)
ctypedef complex (*complex_function_type)(complex)

DTYPE = np.float
ctypedef cnp.float_t DTYPE_t
CTYPE = np.complex128
ctypedef cnp.complex128_t CTYPE_t

cdef double complex _I = _Complex_I
cdef double complex _I2 = _I + _I

# some things that are not in the sage.libs.mpfr
cdef extern from "mpfr.h":
    int mpfr_mul_d (mpfr_t, mpfr_t, double, mpfr_rnd_t)

from sage.matrix.matrix_dense cimport *
from psage.rings.mpc_extras cimport *
from psage.modules.vector_complex_dense cimport Vector_complex_dense
from psage.modules.vector_real_mpfr_dense cimport Vector_real_mpfr_dense
from psage.matrix.matrix_complex_dense cimport Matrix_complex_dense
from sage.functions.all import ceil as pceil
from sage.matrix.all import MatrixSpace
from sage.modular.arithgroup.congroup_sl2z import SL2Z
from lpkbessel import besselk_dp
#from lpkbessel cimport besselk_dp_c
#from sage.modular.maass.all import MySubgroup,besselk_dp
from psage.modform.arithgroup.mysubgroup import MySubgroup
from pullback_algorithms cimport pullback_pts_mpc_new_c,pullback_pts_cplx_dp,pullback_pts_real_dp,pullback_pts_cplx_dp_sym
from lpkbessel cimport besselk_dp_c

from psage.modform.arithgroup.mysubgroups_alg import normalize_point_to_cusp_mpfr,pullback_to_Gamma0N_mpfr,apply_sl2z_map_mpfr,normalize_point_to_cusp_dp,apply_sl2z_map_dp
#from mysubgroups_alg cimport _apply_sl2z_map_mpfr

from pullback_algorithms import pullback_pts_dp,pullback_pts_mpc,pullback_pts_mpc_new

from maass_forms_parallel_alg cimport compute_V_cplx_dp_sym_par,SMAT_cplx_par_dp

# cpdef eval_maass_lp(F,double x,double y,int version = 1,int fi=0,int use_pb=1,int verbose=0):
#     r"""
#     Evaluate a Maass form
#     """
#     cdef float R = <float>F._R
#     cdef float Y = <float>F.group().minimal_height()
#     cdef float xx=<float>x
#     cdef float yy=<float>y
#     cdef int a,b,c,d
#     cdef int cj #G._cusps[cj]
#     cdef int ca,cb
#     RF=RealField(F.prec())
#     G=F.group()
#     # pullback
#     if use_pb == 1:
#         x1,y1,a,b,c,d =  G.pullback(x,y,version=version)
#         #print "pullback=",x1,y1
#     else:
#         x1 = x; y1 = y

#     v = G.closest_vertex(x1,y1)
#     cj= G._vertex_data[v]['cusp'] #representative[v]
#     a,b,c,d=G._vertex_data[v]['cusp_map']
#     if a<>1 or b<>0 or c<>0 or d<>1:
#         x2,y2 = apply_sl2z_map_mpfr(RF(x1),RF(y1),a,b,c,d)
#     else:
#         x2=x1;y2=y1
#         # And then normalize to the correct cusp
#     ca,cb = G._cusps[cj]
#     if cj<>0:
#         a,b,c,d=G._cusp_data[cj]['normalizer']
#         wi = RF(G._cusp_data[cj]['width'])
#         x2,y2 = apply_sl2z_map_mpfr(RF(x2),RF(y2),d,-b,-c,a)
#         x3 = x2/wi.sqrt()
#         y3 = y2/wi.sqrt()        
#         #x3,y3 = normalize_point_to_cusp_dp(G,(ca,cb),x2,y2,inv=1)
#     else:
#         x3 = x2; y3=y2
#     if verbose>0:
#         print "x3,y3=",x3,y3
#     #[x3,y3] = normalize_point_to_cusp_dp(G,(ca,cb),x2,y2,inv=1)
#     res=0
#     twopi=RF(2)*RF.pi()
#     if F._sym_type in [0,1]:
#         if F._sym_type==0:
#             fun=cos
#         elif F._sym_type==1:
#             fun=sin
#         arx=twopi*x3
#         ary=twopi*y3
#         for n in range(1,F._M0):
#             term=besselk_dp(R,ary*n,pref=1)*fun(arx*n)
#             res=res+F._coeffs[fi][cj][n]*term
#         res = res*sqrt(y3)
#     else:
#         arx=twopi*x3
#         ary=twopi*y3
#         for n in range(1,F._M0):
#             term=besselk_dp(R,ary*n,pref=1)
#             if verbose>0:
#                 print "term0 =",term
#                 print "term2 =",(cexpi(arx*n)*F._coeffs[fi][cj][n]+cexpi(-arx*n)*F._coeffs[fi][cj][-n])
#             res = res + term*(cexpi(arx*n)*F._coeffs[fi][cj][n]+cexpi(-arx*n)*F._coeffs[fi][cj][-n])
#         res = res*sqrt(y3)
#     ## we have trivial character here...
#     # Recall that we sum larger numbers (with the prefactor in there to avoid cancellation)
#     return res*exp(-RR.pi()*R*0.5)


# cpdef eval_maass_lp_vec(C,double R,int M0,int sym_type,double y,double x0,double x1,int nx):
#     r"""
#     Evaluate a Maass form on a set of nx equally distributed points z=x_i + iy
#     """
#     cdef double  xx
#     cdef double yy=<float>y
#     cdef double h
#     twopi=2.0*M_PI
#     h = float(x1-x0)/float(nx-1)
#     cdef double* kbvec, *coeffs
#     kbvec = <double*>sig_malloc(nx*sizeof(double))
#     kbvec[0] = 0
#     coeffs = <double*>sig_malloc(nx*sizeof(double))
#     ary=twopi*yy
#     for n in range(1,M0):
#         kbvec[n]=sqrt(y)*besselk_dp(R,ary*n)
#     cdef double tmp = 0
#     cdef list res
#     res = []
#     cdef double complex tmpc 
#     if sym_type  == 1:
#         for i in range(0,nx):
#             xx = x0+i*h
#             arx = twopi*xx
#             tmp = 0
#             for n in range(1,M0):
#                 tmp  += C[n].real_part()*kbvec[n]*cos(arx*n)
#             res.append(tmp)
#     elif  sym_type  == -1:
#         for i in range(nx):
#             xx = x0+i*h
#             arx = twopi*xx
#             tmp = 0
#             for n in range(1,M0):
#                 tmp  += C[n].real_part()*kbvec[n]*sin(arx*n)
#             res.append(tmp)
#     else:
#         for i in range(nx):
#             xx = x0+i*h
#             arx = CC(0,twopi*xx)
#             tmpc = 0
#             for n in range(1,M0):
#                 tmpc  += C[n]*kbvec[n]*exp(arx*n)
#             res.append(tmp)

#     ## we have trivial character here...
#     return res




cpdef whittaker_w_dp(double k,double R,double Y,int pref=0):
    rarg = mpmath.mp.mpc(0,R)
    res = mpmath.mp.whitw(k,rarg,Y)
    cdef double rres
    rres = RealField(53)(res.real)
    return <double>rres

@cython.boundscheck(False)
cdef int compute_V_cplx_wt_dp(double complex **V,double R,double Y,double weight,int** Mv,int** Qv,int nc, int cuspidal,int sym_type, int verbose,double *alphas, double *Xm,double ***Xpb,double ***Ypb, double complex ***Cvec):
    r"""
    Set up the matrix for the system of equations giving the Fourier coefficients of the Maass waveforms.
    INPUT:

    - ``R``   -- double (eigenvalue)
    - ``Y``   -- double (the height of the sampling horocycle)
    - ``Ms,Mf``  -- int (The coefficients we want to compute are C(n), Ms<=n<=Mf )
    - ``Qs,Qf``  -- int (The sampling points are X_m, Qs<=m<=Qf)
    - ``alphas`` -- [nc] double array (the shifts of the Fourier expansion at each cusp)
    - ``V``   -- [(Mf-Ms)*nc]^2 double complex matrix (allocated)
    - ``Xm``  -- [Qf-Qs] double array (allocated)
    - ``Xpb`` -- nc*nc*[Qf-Qs] double array (allocated)
    - ``Ypb`` -- nc*nc*[Qf-Qs] double array (allocated)
    - ``Cvec`` -- nc*nc*[Qf-Qs] double complex array (allocated)
    - `` cuspidal`` -- int (set to 1 if we compute cuspidal functions, otherwise zero)
    - `` sym_type`` -- int (set to 0/1 if we compute even/odd functions, otherwise -1)
    - ``verbose`` -- int (verbosity of output)

    """
    cdef int l,j,icusp,jcusp,n,ni,lj,Ml,Ql,s,Qs,Qf,Mf,Ms
    cdef double pi,Y2pi,nrY2pi,argm,argpb,twopi,two,kbes,Qfak,besarg,lr,nr
    cdef double complex ckbes,ctmpV,iargm,twopii,ctmp
    if not cuspidal in [0,1]:
        raise ValueError," parameter cuspidal must be 0 or 1"
    if weight==0.0:
        raise ValueError," Use this routine only for non-zero weight!"
    Ms=Mv[0][0]; Mf=Mv[0][1]
    Qs=Qv[0][0]; Qf=Qv[0][1]
    pi=M_PI #<double>RealField(53).pi() #3.141592653589793238
    two=<double>(2)
    Y2pi=Y*pi*two
    twopi=two*pi
    #cdef int Ql,Ml
    Ql=Qf-Qs+1
    Ml=Mf-Ms+1

    if sym_type in [0,1]:
        raise ValueError,"We haven't implemented symmetries for non-zero weights!"
    else:
        Qfak=<double>(Ql)
    if verbose>0:
        print "Q=",Qv[0][0],Qv[0][1]
        print "Ql=",Ql
        print "Qfak=",Qfak
    cdef double **nvec=NULL
    #print "Qv=",Qv[0],Qv[1]
    #print "Ql=",Ql
    nvec = <double**>sig_malloc(sizeof(double*)*Ml)
    if not nvec: raise MemoryError
    for n in range(Ml):
        nvec[n] = <double*>sig_malloc(sizeof(double)*nc)
    #cdef cnp.ndarray[CTYPE_t,ndim=4] ef1=np.zeros([Ql, nc,nc,Ml], dtype=CTYPE)
    cdef double complex ****ef1=NULL
    cdef double complex ***ef2_c=NULL
    cdef double ***ef2_r=NULL
    if sym_type not in [0,1]:
        ef2_c = <double complex***>sig_malloc(sizeof(double complex**)*Ml)
        if not ef2_c: raise MemoryError
        for n in range(Ml):
            ef2_c[n] = <double complex**>sig_malloc(sizeof(double complex*)*nc)
            for icusp in range(nc):
                ef2_c[n][icusp] = <double complex*>sig_malloc(sizeof(double complex)*Ql)
    else:
        ef2_r = <double***>sig_malloc(sizeof(double**)*Ml)
        if not ef2_r: raise MemoryError
        for n in range(Ml):
            ef2_r[n] = <double**>sig_malloc(sizeof(double*)*nc)
            for icusp in range(nc):
                ef2_r[n][icusp] = <double*>sig_malloc(sizeof(double)*Ql)

    ef1 = <double complex****>sig_malloc(sizeof(double complex***)*Ml)
    if ef1==NULL: raise MemoryError
    for n in range(Ml):
        ef1[n] = <double complex***>sig_malloc(sizeof(double complex**)*nc)
        if ef1[n]==NULL: raise MemoryError
        for jcusp in range(nc):
            ef1[n][jcusp] = <double complex**>sig_malloc(sizeof(double complex*)*nc)
            if ef1[n][jcusp]==NULL: raise MemoryError
            for icusp in range(nc):
                ef1[n][jcusp][icusp] = <double complex*>sig_malloc(sizeof(double complex)*Ql)
                if ef1[n][jcusp][icusp]==NULL: raise MemoryError
    for n in range(Ml):
        for jcusp in range(nc):
            nvec[n][jcusp]=<double>(n+Ms)+alphas[jcusp]
            for j in range(Ql):
                argm=nvec[n][jcusp]*Xm[j]
                ef2_c[n][jcusp][j]=cexpi(-argm)
    for n in range(Ml):
        for jcusp in range(nc):
            for icusp in range(nc):
                for j in range(Ql): #in range(Qs,Qf+1):
                    if Ypb[icusp][jcusp][j]==0: #not Xpb.has_key((icusp,jcusp,j):
                        continue
                    argpb=nvec[n][jcusp]*Xpb[icusp][jcusp][j]
                    ef1[n][icusp][jcusp][j]=cexpi(argpb)
                    #print "Cv=",Cv[icusp,jcusp,j],type(Cv[icusp,jcusp,j])
                    ctmp = Cvec[icusp][jcusp][j]
                    #ef1[j,icusp,jcusp,n]=ef1[j,icusp,jcusp,n]*<CTYPE_t>ctmp
                    ef1[n][icusp][jcusp][j]=ef1[n][icusp][jcusp][j]*ctmp

    cdef double besarg_old=0.0
    cdef double y,kbes_old=1.0
    #tmplist=[]
    cdef double ***kbesvec=NULL
    kbesvec=<double***>sig_malloc(sizeof(double**)*Ml)
    cdef double weight_sign,weight_half
    weight_half = (<double>0.5)*weight
    if kbesvec==NULL:
        raise MemoryError
    for l in range(Ml):
        kbesvec[l]=<double**>sig_malloc(sizeof(double*)*nc)
        if kbesvec[l]==NULL:
            raise MemoryError
        for jcusp in range(nc):
            kbesvec[l][jcusp]=<double*>sig_malloc(sizeof(double)*Ql)
            if kbesvec[l][jcusp]==NULL:
                raise MemoryError
        for jcusp in range(nc):
            lr=nvec[l][jcusp]*twopi
            if lr>0:
                weight_sign = weight_half
            else:
                weight_sign = -weight_half
            for j in range(Ql):
                # I am trying to make use of the fact that "many" pullback matrices
                # are invariant under the Fricke involution
                for icusp in range(nc):
                    if Ypb[icusp][jcusp][j]<>0:
                        besarg=two*abs(lr)*Ypb[icusp][jcusp][j]
                        if lr<>0.0:
                            kbes = whittaker_w_dp(weight_sign,R,besarg,pref=1)
                            kbesvec[l][icusp][j]=kbes
                        elif cuspidal==0:
                            kbesvec[l][icusp][j]=<double>1.0
                        else:
                            kbesvec[l][icusp][j]=<double>0.0

    for l in range(Ml):
        for jcusp in range(nc): 
            lr=nvec[l][jcusp]*twopi
            if lr==0.0 and cuspidal==1:
                continue
            lj=Ml*jcusp+l
            if lr>=0:
                weight_sign = weight_half
            else:
                weight_sign = -weight_half
            for icusp in range(nc):
                for j in range(Ql): 
                    if Ypb[icusp][jcusp][j]==0:
                        #if(not Ypb.has_key((icusp,jcusp,j))):
                        continue
                    ckbes=kbesvec[l][icusp][j]*ef1[l][icusp][jcusp][j]
                    #if do_real_norm and is_trivial:
                    #    pass
                    if verbose>2:
                        if abs(lr)>0:
                            besarg = two*fabs(lr)*Ypb[icusp][jcusp][j]
                            kbes=whittaker_w_dp(weight_sign,R,besarg,pref=1)
                            if abs(kbesvec[l][icusp][j]-kbes)>1E-12:
                                print "whitW(",weight_sign,R,besarg,")"
                                print "l,icusp,jcusp,j=",l,icusp,jcusp,j
                                print "whitWvec:=",kbesvec[l][icusp][j]
                                print "whitW=",kbes #/sqrt(fabs(nvec[l][jcusp]))
                    if verbose>3 and abs(l)<=2:
                        print "kbes:=",kbesvec[l][icusp][j]
                        print "kbes1=",whittaker_w_dp(weight_sign,R,besarg,pref=1) #/sqrt(fabs(nvec[l][jcusp]))
                        print "ef1=",ef1[l][icusp][jcusp][j]
                    for n in range(Ml):
                        if nvec[n][icusp]==0.0 and cuspidal==1:
                            continue
                        ni=Ml*icusp+n
                        ctmpV=ckbes*ef2_c[n][icusp][j]
                        V[ni][lj]=V[ni][lj]+ctmpV # <CTYPE_t> creal(ctmpV)+_I*cimag(ctmpV)
                        if ni==3 and lj==3 and verbose>2:
                            print "whitW(",weight_sign,R,besarg,")"
                            print "V[",ni,"][",lj,"][",j,"]=",V[ni][lj]
                            print "ctmp=",ctmpV
                            print "ckbes=",ckbes
                            print "kbesvec[",l,icusp,j,"=",kbesvec[l][icusp][j]
                            print "ef1[",l,icusp,jcusp,j,"]=",ef1[l][icusp][jcusp][j]
                            if sym_type==1 or sym_type==0:
                                print "ef2_r=",ef2_r[n][icusp][j]
                            else:
                                print "ef2_c=",ef2_c[n][icusp][j]
                        #    #print "ctmpV=",ctmpV

    #save(tmplist,"besarg.sobj")
    #if verbose>=0:
    #    print "V[12,12]0=",V[12][12]
    cdef double fak1
    for l in range(Ml):
        for jcusp in range(nc):            
            lj = Ml*jcusp+l
            if nvec[l][jcusp]==0:
                fak1=Qfak
            else:
                fak1 = Qfak*sqrt(fabs(nvec[l][jcusp]))
            for n in range(Ml):
                for icusp in range(nc):            
                    ni = Ml*icusp+n
                    V[ni][lj]=V[ni][lj]/fak1
    if verbose>2:
        print "V[3,3]1=",V[3][3]
        print "Y2pi=",Y2pi
    for n in range(Ml):
        #for icusp in range(nc):
        for icusp in range(nc):
            nr=fabs(nvec[n][icusp])
            if nvec[n][icusp]>0:                
                weight_sign = weight_half
            else:
                weight_sign = -weight_half
            if nr==0.0 and cuspidal==1:
                continue
            ni=Ml*icusp+n
            nrY2pi=nr*Y2pi
            if nr<>0.0:
                kbes=whittaker_w_dp(weight_sign,R,two*nrY2pi,pref=1)
                kbes = kbes/sqrt(nr)
            elif cuspidal==0:
                kbes=<double>1.0
            else:
                kbes=<double>0.0
                
            #if nrY2pi==0.0:
            #    print "arg=0!"
            if verbose>2 and ni==0:
                print "n,icusp=",n+Ms,icusp
                print "arg=",nrY2pi
                print "ckbes(",ni,")=whittaker_w_dp(",weight_sign,R,two*nrY2pi,")=",kbes
            V[ni][ni]=V[ni][ni] - kbes
    if verbose>2:
        print "V[3,3]2=",V[3][3]
    if kbesvec<>NULL:
        for l in range(Ml):
            if kbesvec[l]<>NULL:
                for icusp in range(nc):
                    if kbesvec[l][icusp]<>NULL:
                        sig_free(kbesvec[l][icusp])
                sig_free(kbesvec[l])
        sig_free(kbesvec)
    if ef1<>NULL:
        for n in range(Ml):
            if ef1[n]<>NULL:
                for jcusp in range(nc): 
                    if ef1[n][jcusp]<>NULL:
                        for icusp in range(nc):
                            if ef1[n][jcusp][icusp]<>NULL:
                                sig_free(ef1[n][jcusp][icusp])
                        sig_free(ef1[n][jcusp])
                sig_free(ef1[n])
        sig_free(ef1)
    if ef2_r<>NULL:
        for n in range(Ml):
            sig_free(ef2_r[n])
        sig_free(ef2_r)
    if ef2_c<>NULL:
        for n in range(Ml): 
            sig_free(ef2_c[n])
        sig_free(ef2_c)
    if nvec<>NULL:
        for n in range(Ml):
            if nvec[n]<>NULL:
                sig_free(nvec[n])
        sig_free(nvec)


cpdef my_kbes(double r,double x,double besprec = 1e-14,int pref=0):
    cdef double kbes
    if abs(r) < 500:
        besselk_dp_c(&kbes,r,x,besprec,pref=pref)
    else:
        ir = mpmath.mp.mpc(0.5,r)
        mpx = mpmath.mp.mpf(x)
        if pref == 1:
            arg = mpmath.mp.mpf(r)*mpmath.mp.pi*mpmath.mp.mpf(0.5)
            res = mpmath.mp.besselk(ir,mpx).real*mpmath.mp.exp(arg)
            kbes = <double> res
        else:
            kbes = <double> mpmath.mp.besselk(ir,mpx).real ## probably only going to be 0....
    return kbes
@cython.boundscheck(False)
cdef int compute_V_cplx_dp(double complex **V,double R,double Y,int** Mv,int** Qv,int nc, int cuspidal,int sym_type, int verbose,double *alphas, double *Xm,double ***Xpb,double ***Ypb, double complex ***Cvec,int is_exceptional=0):
    r"""
    Set up the matrix for the system of equations giving the Fourier coefficients of the Maass waveforms.
    INPUT:

    - ``R``   -- double (eigenvalue)
    - ``Y``   -- double (the height of the sampling horocycle)
    - ``Ms,Mf``  -- int (The coefficients we want to compute are C(n), Ms<=n<=Mf )
    - ``Qs,Qf``  -- int (The sampling points are X_m, Qs<=m<=Qf)
    - ``alphas`` -- [nc] double array (the shifts of the Fourier expansion at each cusp)
    - ``V``   -- [(Mf-Ms)*nc]^2 double complex matrix (allocated)
    - ``Xm``  -- [Qf-Qs] double array (allocated)
    - ``Xpb`` -- nc*nc*[Qf-Qs] double array (allocated)
    - ``Ypb`` -- nc*nc*[Qf-Qs] double array (allocated)
    - ``Cvec`` -- nc*nc*[Qf-Qs] double complex array (allocated)
    - `` cuspidal`` -- int (set to 1 if we compute cuspidal functions, otherwise zero)
    - `` sym_type`` -- int (set to 0/1 if we compute even/odd functions, otherwise -1)
    - ``verbose`` -- int (verbosity of output)

    """
    cdef int l,j,icusp,jcusp,n,ni,lj,Ml,Ql,s,Qs,Qf,Mf,Ms
    cdef double pi,Y2pi,nrY2pi,argm,argpb,twopi,two,kbes,Qfak,besarg,lr,nr,tmpr
    cdef double complex ckbes,ctmpV,iargm,twopii,ctmp
    if not cuspidal in [0,1]:
        raise ValueError," parameter cuspidal must be 0 or 1"
    Ms=Mv[0][0]; Mf=Mv[0][1]
    Qs=Qv[0][0]; Qf=Qv[0][1]
    pi=M_PI #<double>RealField(53).pi() #3.141592653589793238
    two=<double>(2)
    Y2pi=Y*pi*two
    twopi=two*pi
    #cdef int Ql,Ml
    Ql=Qf-Qs+1
    Ml=Mf-Ms+1
    cdef double besprec
    besprec=1.0E-14

    if sym_type in [0,1]:
        Qfak=<double>(Ql)/<double>(2)
    else:
        Qfak=<double>(Ql)
    cdef double **nvec=NULL
    #print "Qv=",Qv[0],Qv[1]
    #print "Ql=",Ql
    nvec = <double**>sig_malloc(sizeof(double*)*Ml)
    if not nvec: raise MemoryError
    for n in range(Ml):
        nvec[n] = <double*>sig_malloc(sizeof(double)*nc)
    #cdef cnp.ndarray[CTYPE_t,ndim=4] ef1=np.zeros([Ql, nc,nc,Ml], dtype=CTYPE)
    cdef double complex ****ef1=NULL
    cdef double complex ***ef2_c=NULL
    cdef double ***ef2_r=NULL
    if sym_type not in [0,1]:
        ef2_c = <double complex***>sig_malloc(sizeof(double complex**)*Ml)
        if not ef2_c: raise MemoryError
        for n in range(Ml):
            ef2_c[n] = <double complex**>sig_malloc(sizeof(double complex*)*nc)
            for icusp in range(nc):
                ef2_c[n][icusp] = <double complex*>sig_malloc(sizeof(double complex)*Ql)
    else:
        ef2_r = <double***>sig_malloc(sizeof(double**)*Ml)
        if not ef2_r: raise MemoryError
        for n in range(Ml):
            ef2_r[n] = <double**>sig_malloc(sizeof(double*)*nc)
            for icusp in range(nc):
                ef2_r[n][icusp] = <double*>sig_malloc(sizeof(double)*Ql)

    ef1 = <double complex****>sig_malloc(sizeof(double complex***)*Ml)
    if ef1==NULL: raise MemoryError
    for n in range(Ml):
        ef1[n] = <double complex***>sig_malloc(sizeof(double complex**)*nc)
        if ef1[n]==NULL: raise MemoryError
        for jcusp in range(nc):
            ef1[n][jcusp] = <double complex**>sig_malloc(sizeof(double complex*)*nc)
            if ef1[n][jcusp]==NULL: raise MemoryError
            for icusp in range(nc):
                ef1[n][jcusp][icusp] = <double complex*>sig_malloc(sizeof(double complex)*Ql)
                if ef1[n][jcusp][icusp]==NULL: raise MemoryError
    for n in range(Ml):
        for jcusp in range(nc):
            nvec[n][jcusp]=<double>(n+Ms)+alphas[jcusp]
            for j in range(Ql):
                argm=nvec[n][jcusp]*Xm[j]
                if sym_type==0:
                    ef2_r[n][jcusp][j]=cos(argm)
                elif sym_type==1:
                    ef2_r[n][jcusp][j]=sin(argm)
                else:
                    ef2_c[n][jcusp][j]=cexpi(-argm)
    for n in range(Ml):
        for jcusp in range(nc):
            for icusp in range(nc):
                for j in range(Ql): #in range(Qs,Qf+1):
                    if Ypb[icusp][jcusp][j]==0: #not Xpb.has_key((icusp,jcusp,j):
                        continue
                    argpb=nvec[n][jcusp]*Xpb[icusp][jcusp][j]
                    if sym_type==0:
                        ef1[n][icusp][jcusp][j]=2.0*cos(argpb)
                    elif sym_type==1:
                        ef1[n][icusp][jcusp][j]=sin(argpb)
                    else:
                        ef1[n][icusp][jcusp][j]=cexpi(argpb)
                    #print "Cv=",Cv[icusp,jcusp,j],type(Cv[icusp,jcusp,j])
                    ctmp = Cvec[icusp][jcusp][j]
                    #ef1[j,icusp,jcusp,n]=ef1[j,icusp,jcusp,n]*<CTYPE_t>ctmp
                    ef1[n][icusp][jcusp][j]=ef1[n][icusp][jcusp][j]*ctmp

    cdef double besarg_old=0.0
    cdef double y,kbes_old=1.0
    #tmplist=[]
    cdef double ***kbesvec=NULL
    kbesvec=<double***>sig_malloc(sizeof(double**)*Ml)
    if kbesvec==NULL:
        raise MemoryError
    for l in range(Ml):
        kbesvec[l]=<double**>sig_malloc(sizeof(double*)*nc)
        if kbesvec[l]==NULL:
            raise MemoryError
        for jcusp in range(nc):
            kbesvec[l][jcusp]=<double*>sig_malloc(sizeof(double)*Ql)
            if kbesvec[l][jcusp]==NULL:
                raise MemoryError
        for jcusp in range(nc):
            lr=nvec[l][jcusp]*twopi
            for j in range(Ql):
                # I am trying to make use of the fact that "many" pullback matrices
                # are invariant under the Fricke involution
                for icusp in range(nc):
                    if Ypb[icusp][jcusp][j]<>0:
                        besarg=abs(lr)*Ypb[icusp][jcusp][j]
                        if lr<>0.0: 
                            if is_exceptional == 0:
                                #besselk_dp_c(&tmpr,R,besarg,besprec,1)
                                tmpr = my_kbes(R,besarg,besprec,1)
                            else:
                                tmpr = scipy.special.kv(R,besarg)
                            #print "Ypb=",Ypb[icusp][jcusp][j],type(Ypb[icusp][jcusp][j])
                            #print "tmpr=",tmpr,type(tmpr)
                            kbesvec[l][icusp][j]=sqrt(Ypb[icusp][jcusp][j])*tmpr
                            #kbesvec[l][icusp][j]=whittaker_w_dp(weight_sign,R,besarg,pref=1)
                        else:
                            kbesvec[l][icusp][j]=<double>1.0

    for l in range(Ml):
        for jcusp in range(nc):
            lr=nvec[l][jcusp]*twopi
            if lr==0.0 and cuspidal==1:
                continue
            lj=Ml*jcusp+l
            for icusp in range(nc):
                for j in range(Ql):
                    if Ypb[icusp][jcusp][j]==0:
                        #if(not Ypb.has_key((icusp,jcusp,j))):
                        continue
                    ckbes=kbesvec[l][icusp][j]*ef1[l][icusp][jcusp][j]

                    #    pass
                    if verbose>2:
                        if abs(lr)>0:
                            if is_exceptional == 0:
                                #kbes=besselk_dp(R,abs(lr)*Ypb[icusp][jcusp][j],pref=1)
                                kbes=my_kbes(R,abs(lr)*Ypb[icusp][jcusp][j],pref=1)
                            else:
                                besarg = abs(lr)*Ypb[icusp][jcusp][j]
                                kbes = scipy.special.kv(R,besarg)
                            kbes = kbes*sqrt(Ypb[icusp][jcusp][j])
                            if abs(kbesvec[l][icusp][j]-kbes)>1E-12:
                                print "l,icusp,jcusp,j=",l,icusp,jcusp,j
                                print "kbes:=",kbesvec[l][icusp][j]
                                print "kbes1=",kbes
                    if verbose>3 and abs(l)<=2:
                        print "kbes:=",kbesvec[l][icusp][j]
                        #besselk_dp(R,abs(lr)*Ypb[icusp][jcusp][j],pref=1)
                        kbes = my_kbes(R,abs(lr)*Ypb[icusp][jcusp][j],pref=1)
                        print "kbes1=",sqrt(Ypb[icusp][jcusp][j])*kbes
                        print "ef1=",ef1[l][icusp][jcusp][j]
                    for n in range(Ml):
                        if n+Ms==0 and cuspidal==1:
                            continue
                        ni=Ml*icusp+n
                        if sym_type==1 or sym_type==0:
                            ctmpV=ckbes*ef2_r[n][icusp][j]
                        else:
                            ctmpV=ckbes*ef2_c[n][icusp][j]
                        V[ni][lj]=V[ni][lj]+ctmpV # <CTYPE_t> creal(ctmpV)+_I*cimag(ctmpV)
                        if ni==3 and lj==3 and verbose>2:
                            print "V[",ni,"][",lj,"][",j,"]=",V[ni][lj]
                            print "ctmp=",ctmpV
                            print "ckbes=",ckbes
                            print "kbesvec[",l,icusp,j,"=",kbesvec[l][icusp][j]
                            print "ef1[",l,icusp,jcusp,j,"]=",ef1[l][icusp][jcusp][j]
                            if sym_type==1 or sym_type==0:
                                print "ef2_r=",ef2_r[n][icusp][j]
                            else:
                                print "ef2_c=",ef2_c[n][icusp][j]
                        #    #print "ctmpV=",ctmpV

    #save(tmplist,"besarg.sobj")
    #if verbose>=0:
    #    print "V[12,12]0=",V[12][12]
    cdef double sqrtY
    sqrtY=sqrt(Y)
    for n in range(Ml*nc):
        for l in range(Ml*nc):
            V[n][l]=V[n][l]/Qfak
    if verbose>2:
        print "V[3,3]1=",V[3][3]
        print "Y2pi=",Y2pi
    for n in range(Ml):
        #for icusp in range(nc):
        for icusp in range(nc):
            nr=abs(nvec[n][icusp])
            if nr==0.0 and cuspidal==1:
                continue
            ni=Ml*icusp+n
            nrY2pi=nr*Y2pi
            if nr==0.0:
                kbes=<double>1.0
            else:
                if is_exceptional == 0:
                    kbes=sqrtY*my_kbes(R,nrY2pi,pref=1)
                    #besselk_dp(R,nrY2pi,pref=1)
                else:
                    kbes = sqrtY*scipy.special.kv(R,nrY2pi)
            #if nrY2pi==0.0:
            #    print "arg=0!"
            if verbose>2 and ni==0:
                print "n,icusp=",n+Ms,icusp
                print "arg=",nrY2pi
                print "ckbes(",ni,")=sqrt(Y)*kbessel(",R,nrY2pi,")=",kbes
            V[ni][ni]=V[ni][ni] - kbes
    if verbose>2:
        print "V[3,3]2=",V[3][3]
    if kbesvec<>NULL:
        for l in range(Ml):
            if kbesvec[l]<>NULL:
                for icusp in range(nc):
                    if kbesvec[l][icusp]<>NULL:
                        sig_free(kbesvec[l][icusp])
                sig_free(kbesvec[l])
        sig_free(kbesvec)
    if ef1<>NULL:
        for n in range(Ml):
            if ef1[n]<>NULL:
                for jcusp in range(nc):
                    if ef1[n][jcusp]<>NULL:
                        for icusp in range(nc):
                            if ef1[n][jcusp][icusp]<>NULL:
                                sig_free(ef1[n][jcusp][icusp])
                        sig_free(ef1[n][jcusp])
                sig_free(ef1[n])
        sig_free(ef1)
    if ef2_r<>NULL:
        for n in range(Ml):
            if ef2_r[n]<>NULL:
                for icusp in range(nc):
                    sig_free(ef2_r[n][icusp])
                sig_free(ef2_r[n])    
        sig_free(ef2_r)
    if ef2_c<>NULL:
        for n in range(Ml):
            sig_free(ef2_c[n])
        sig_free(ef2_c)
    if nvec<>NULL:
        for n in range(Ml):
            if nvec[n]<>NULL:
                sig_free(nvec[n])
        sig_free(nvec)





@cython.boundscheck(False)
@cython.cdivision(True)
cdef int compute_V_cplx_dp_sym_wt(double complex **V,
                           int N1,
                           double *Xm,
                           double ***Xpb,
                           double ***Ypb,
                           double complex ***Cvec,
                           double complex *cusp_evs,
                           double *alphas,
                           int **Mv,int **Qv,double *Qfak,
                           int *symmetric_cusps,
                           double R,double Y,
                           int nc, int ncols,
                           int cuspidal,
                           int verbose       ):


    r"""
    Set up the matrix for the system of equations giving the Fourier coefficients of the Maass waveforms.
    INPUT:

    - ``R``   -- double (eigenvalue)
    - ``Y``   -- double (the height of the sampling horocycle)
    - ``Ms,Mf``  -- int (The coefficients we want to compute are C(n), Ms<=n<=Mf )
    - ``Qs,Qf``  -- int (The sampling points are X_m, Qs<=m<=Qf)
    - ``alphas`` -- [nc] double array (the shifts of the Fourier expansion at each cusp)
    - ``V``   -- [(Mf-Ms)*nc]^2 double complex matrix (allocated)
    - ``Xm``  -- [Qf-Qs] double array (allocated)
    - ``Xpb`` -- nc*nc*[Qf-Qs] double array (allocated)
    - ``Ypb`` -- nc*nc*[Qf-Qs] double array (allocated)
    - ``Cvec`` -- nc*nc*[Qf-Qs] double complex array (allocated)
    - `` cuspidal`` -- int (set to 1 if we compute cuspidal functions, otherwise zero)
    - `` sym_type`` -- int (set to 0/1 if we compute even/odd functions, otherwise -1)
    - ``verbose`` -- int (verbosity of output)

    """
    raise NotImplementedError("Have not implemented symmetries and weight!")
    cdef int l,j,icusp,jcusp,n,ni,lj,Ml,Ql,s,Qs,Qf,Mf,Ms
    cdef double pi,sqrtY,Y2pi,nrY2pi,argm,argpb,twopi,two,kbes,besarg,lr,nr
    cdef double complex ckbes,ctmpV,iargm,twopii,ctmp
    if not cuspidal in [0,1]:
        raise ValueError," parameter cuspidal must be 0 or 1"
    if not cuspidal in [0,1]:
        raise ValueError," parameter cuspidal must be 0 or 1"
    if R < 0:
        ## In this case (corresponding to lambda in [0,1/4] we use the real parameter K-Bessel
        R = -R
        set_pref = -1
    else:
        set_pref = 1
    pi=M_PI 
    sqrtY=sqrt(Y)
    two=<double>(2)
    Y2pi=Y*pi*two
    twopi=two*pi
    Ml=0; Ql=0
    for i in range(nc):
        if Mv[i][2]>Ml:
            Ml=Mv[i][2]
        if Qv[i][2]>Ql:
            Ql=Qv[i][2]
        if verbose>0:
            print "Qv[",i,"]=(",Qv[i][0],",",Qv[i][1],",",Qv[i][2],")"
            print "Mv[",i,"]=(",Mv[i][0],",",Mv[i][1],",",Mv[i][2],")"
    if verbose>0:
        print "N1=",N1
        print "Ql=",Ql
    ## This is the effective offset at the
    cdef int* cusp_offsets=NULL
    cusp_offsets=<int*>sig_malloc(sizeof(int)*nc)
    if cusp_offsets==NULL: raise MemoryError
    for jcusp in range(nc):
        cusp_offsets[jcusp]=0
        for icusp in range(jcusp):
            if icusp==0 or cusp_evs[icusp]==0:
                cusp_offsets[jcusp]+=Mv[icusp][2]
        if verbose>0:
            print "cusp_offsets[",jcusp,"]=",cusp_offsets[jcusp]
    cdef int nc_sym=0
    for jcusp in range(nc):
        if verbose>0:
            print "cusp_evs[",jcusp,"]=",cusp_evs[jcusp]
        if jcusp==0 or cusp_evs[jcusp]<>0:
            nc_sym+=1
    cdef double **nvec=NULL
    nvec = <double**>sig_malloc(sizeof(double*)*nc)
    if not nvec: raise MemoryError
    for icusp in range(nc):
        nvec[icusp] = <double*>sig_malloc(sizeof(double)*Ml)
    cdef double complex ****ef1=NULL
    cdef double complex ***ef2_c=NULL
    ef2_c = <double complex***>sig_malloc(sizeof(double complex**)*nc)
    if not ef2_c: raise MemoryError
    for icusp in range(nc):
        ef2_c[icusp] = <double complex**>sig_malloc(sizeof(double complex*)*Mv[icusp][2])
        for n in range(Mv[icusp][2]):
            ef2_c[icusp][n] = <double complex*>sig_malloc(sizeof(double complex)*Qv[icusp][2])
    ef1 = <double complex****>sig_malloc(sizeof(double complex***)*nc)
    if ef1==NULL: raise MemoryError
    for icusp in range(nc):
        ef1[icusp] = <double complex***>sig_malloc(sizeof(double complex**)*nc)
        if ef1[icusp]==NULL: raise MemoryError
        for jcusp in range(nc):
            ef1[icusp][jcusp] = <double complex**>sig_malloc(sizeof(double complex*)*Mv[jcusp][2])
            if ef1[icusp][jcusp]==NULL: raise MemoryError
            for n in range(Mv[jcusp][2]):
                ef1[icusp][jcusp][n] = <double complex*>sig_malloc(sizeof(double complex)*Qv[jcusp][2])
                if ef1[icusp][jcusp][n]==NULL: raise MemoryError
    for jcusp in range(nc):
        for n in range(Ml):
            nvec[jcusp][n]=<double>(n+Mv[jcusp][0])+alphas[jcusp]
    cdef int twoQm1
    twoQm1= 2*Qv[0][1]-1
    for jcusp in range(nc):
        for n in range(Mv[jcusp][2]):
            for j in range(Qv[jcusp][2]):
                argm=nvec[jcusp][n]*Xm[j]
                if symmetric_cusps[jcusp]==0:
                    ef2_c[jcusp][n][j]=2.0*cos(argm)
                elif symmetric_cusps[jcusp]==1:
                    ef2_c[jcusp][n][j]=_I2*sin(-argm)
                else:
                    ef2_c[jcusp][n][j]=cexpi(-argm)

    cdef double argpb1
    for jcusp in range(nc):
        for icusp in range(nc):
            for n in range(Mv[jcusp][2]):
                for j in range(Qv[jcusp][2]): #in range(Qs,Qf+1):
                    if Ypb[icusp][jcusp][j]==0: #not Xpb.has_key((icusp,jcusp,j):
                        continue
                    argpb=nvec[jcusp][n]*Xpb[icusp][jcusp][j]
                    if symmetric_cusps[jcusp]==0:
                        ef1[icusp][jcusp][n][j]=2.0*cos(argpb)
                    elif symmetric_cusps[jcusp]==1:
                        ef1[icusp][jcusp][n][j]=_I2*sin(argpb)
                    else:
                        ef1[icusp][jcusp][n][j]=cexpi(argpb)
                    ctmp = Cvec[icusp][jcusp][j]
                    ef1[icusp][jcusp][n][j]=ef1[icusp][jcusp][n][j]*ctmp

    if verbose>1:
        print "here1121"
    cdef double besarg_old=0.0
    cdef double y,kbes_old=1.0
    cdef double ***kbesvec=NULL
    kbesvec=<double***>sig_malloc(sizeof(double**)*nc)
    if kbesvec==NULL:
        raise MemoryError
    for jcusp in range(nc):
        kbesvec[jcusp]=<double**>sig_malloc(sizeof(double*)*Ml)
        if kbesvec[jcusp]==NULL:
            raise MemoryError
        for l in range(Ml):
            kbesvec[jcusp][l]=<double*>sig_malloc(sizeof(double)*Ql) #Qv[jcusp][2])
            if kbesvec[jcusp][l]==NULL:
                raise MemoryError
    if verbose>0:
        print "here0"
        print "Ml=",Ml
        print "Ql=",Ql
    cdef double tmpr
    cdef double besprec
    besprec=1.0E-14
    ## Can we somehow use that "most" of Xpb[icusp][jcusp][2Q-j-1]=-Xpb[icusp][jcusp][j] ?
    ## uncomment the below lines to see this.
    # for jcusp in range(nc):
    #     for icusp in range(nc):
    #         for j from 0<=j<Qv[1][1]:
    #             print "diff[",jcusp,icusp,j,"]=",Xpb[icusp][jcusp][j]+Xpb[icusp][jcusp][2*Qv[jcusp][1]-j-1]

    for jcusp in range(nc):
        for icusp in range(nc):
            if icusp>0 and cusp_evs[icusp]<>0:
                continue
            for j in range(Qv[jcusp][2]):
                if Ypb[icusp][jcusp][j]==0:
                    continue
                for l in range(Mv[jcusp][2]):
                    lr=nvec[jcusp][l]*twopi
                    Mf = Mv[icusp][1]
                    besarg=fabs(lr)*Ypb[icusp][jcusp][j]
                    if lr<>0.0:
                        #besselk_dp_c(&tmpr,R,besarg,besprec,pref=set_pref)
                        tmpr = my_kbes(R,besarg,besprec,pref=set_pref)
                        kbesvec[icusp][l][j]=sqrt(Ypb[icusp][jcusp][j])*tmpr
                    else:
                        kbesvec[icusp][l][j]=<double>1.0

    cdef double complex cuspev

    for jcusp in range(nc):
        for l in range(Mv[jcusp][2]):
            lr=nvec[jcusp][l]*twopi
            lj=cusp_offsets[jcusp]+l
            if jcusp>0 and cusp_evs[jcusp]<>0:
                lj0=l; cuspev=cusp_evs[jcusp]
            else:
                lj0=lj;  cuspev=1.0
            if lr==0.0 and cuspidal==1:
                continue
            for j in range(Qv[jcusp][2]):
                for icusp in range(nc):
                    if icusp>0 and cusp_evs[icusp]<>0:
                        continue
                    if Ypb[icusp][jcusp][j]==0:
                        continue
                    ckbes=kbesvec[icusp][l][j]*ef1[icusp][jcusp][l][j]
                    for n in range(Mv[icusp][2]):
                        if nvec[icusp][n]==0 and cuspidal==1:
                            continue
                        #if n+Mv[icusp][0]==0 and cuspidal==1:
                        #    continue
                        ni=cusp_offsets[icusp]+n
                        ctmpV=ckbes*ef2_c[icusp][n][j]
                        if ni>N1 or lj0>N1:
                            raise ArithmeticError,"Index outside!"
                        V[ni][lj0]=V[ni][lj0]+ctmpV*cuspev
                        # if verbose>0 and lj0==324 and ni==46:
                        if verbose>2 and ni==3 and lj0==3:
                            print "l,jcusp,n,icusp=",l+Mv[jcusp][0],jcusp,n+Mv[icusp][0],icusp
                            print "Ypb(",j+Qv[jcusp][0],")=",Ypb[icusp][jcusp][j]
                            print "V[",ni,"][",lj0,"](",j,")=",V[ni][lj0]
                            print "ctmpV[",ni,"][",lj0,"]=",ctmpV
                            print "kbesvec[",icusp,l,j,"]=",kbesvec[icusp][l][j]
                            print "besarg=",fabs(lr)*Ypb[icusp][jcusp][j]
                            print "ef1=",ef1[icusp][jcusp][l][j]
                            print "ef2=",ef2_c[icusp][n][j]
                            print "V[",ni,"][",lj0,"]=",V[ni][lj0],"+",ctmpV,"*",cuspev

    if verbose>2:
    #    #for n from 30+cusp_offsets[3]<=n<40+cusp_offsets[3]:
    #    for n from 43<=n<49:
    #        print "V0[",0,n,"]=",V[0][n]
    #        #print "V0[",24,24,"]=",V[24][24]
        print "V0[",3,3,"]=",V[3][3]
        for icusp in range(nc):
            print "Qfak[{0}]={1}".format(icusp,Qfak[icusp])
    for icusp in range(nc):
        if icusp>0 and cusp_evs[icusp]<>0:
            continue
        for n in range(Mv[icusp][2]):
            ni=cusp_offsets[icusp]+n
            for jcusp in range(nc):
                if jcusp>0 and cusp_evs[jcusp]<>0:
                    continue
                for l in range(Mv[jcusp][2]):
                    lj=cusp_offsets[jcusp]+l
                    if lj>N1 or ni>N1: # some extra insurance...
                        print "ni=",cusp_offsets[icusp],"+",n,"=",ni
                        print "lj=",cusp_offsets[jcusp],"+",l,"=",lj
                        raise ArithmeticError,"Index outside!"
                    if verbose>2 and ni==3 and lj==3:
                        print "V[3][3]={0}/{1}".format(V[ni][lj],Qfak[jcusp])
                    V[ni][lj]=V[ni][lj]/<double complex>Qfak[jcusp]
                    if verbose>2 and ni==3 and lj==3:
                        print "V[3][3]/Qfak[{0}]={1}".format(jcusp,V[ni][lj])

    for icusp in range(nc):
        if icusp>0 and cusp_evs[icusp]<>0:
            continue
        for n in range(Mv[icusp][2]):
            nr=fabs(nvec[icusp][n])
            if nr==0.0 and cuspidal==1:
                continue
            ni=cusp_offsets[icusp]+n
            nrY2pi=nr*Y2pi
            if nr==0.0:
                kbes=<double>1.0
            else:
                #mpIR=mpmath.fp.mpc(0,R)
                #                kbes=float(mpmath.fp.besselk(mpIR,nrY2pi).real*exp(mpmath.fp.pi*R*0.5))
                #besselk_dp_c(&kbes,R,nrY2pi,besprec,pref=set_pref)
                kbes = my_kbes(R,nrY2pi,besprec,pref=set_pref)
                kbes=sqrtY*kbes # besselk_dp(R,nrY2pi,pref=1)
            if ni>N1:
                raise ArithmeticError,"Index outside!"
            if verbose>2 and ni==3:
                print "V[3][3]={0}-{1}={2}".format(V[ni][ni],kbes,V[ni][ni]-kbes)
            V[ni][ni]=V[ni][ni] - kbes

    if verbose>3:

        #        #for n from 30+cusp_offsets[3]<=n<40+cusp_offsets[3]:
        #        print "V3[",0,n,"]=",V[0][n]
        for n in range(min(100,Ml)):
            print "V3[",n,n,"]=",V[n][n]
    #if Qfak<>NULL:
    #    #for icusp in range(nc):
    #    #    if Qfak[icusp]<>NULL:
    #    #        sig_free(Qfak[icusp])
    #    sig_free(Qfak)
    if kbesvec<>NULL:
        for icusp in range(nc):
            if kbesvec[icusp]<>NULL:
                for l in range(Ml):
                    if kbesvec[icusp][l]<>NULL:
                        sig_free(kbesvec[icusp][l])
                sig_free(kbesvec[icusp])
        sig_free(kbesvec)
    #print "deal kbbes1"
    if ef1<>NULL:
        for jcusp in range(nc):
            if ef1[jcusp]<>NULL:
                for icusp in range(nc):
                    if ef1[jcusp][icusp]<>NULL:
                        for n in range(Mv[icusp][2]):
                            if ef1[jcusp][icusp][n]<>NULL:
                                sig_free(ef1[jcusp][icusp][n])
                        sig_free(ef1[jcusp][icusp])
                sig_free(ef1[jcusp])
        sig_free(ef1)
    if ef2_c<>NULL:
        for icusp in range(nc):
            if ef2_c[icusp]<>NULL:
                for n in range(Mv[icusp][2]):
                    if ef2_c[icusp][n]<>NULL:
                        sig_free(ef2_c[icusp][n])
                sig_free(ef2_c[icusp])
        sig_free(ef2_c)
    if nvec<>NULL:
        for icusp in range(nc):
            if nvec[icusp]<>NULL:
                sig_free(nvec[icusp])
        sig_free(nvec)
    if cusp_offsets<>NULL:
        sig_free(cusp_offsets)


@cython.boundscheck(False)
@cython.cdivision(True)
cdef int compute_V_cplx_dp_sym(double complex **V,
                           int N1,
                           double *Xm,
                           double ***Xpb,
                           double ***Ypb,
                           double complex ***Cvec,
                           double complex *cusp_evs,
                           double *alphas,
                           int **Mv,int **Qv,double *Qfak,
                           int *symmetric_cusps,
                           double R,double Y,
                           int nc, int ncols,
                           int cuspidal,
                           int verbose,
                           int is_exceptional=0) except -1:


    r"""
    Set up the matrix for the system of equations giving the Fourier coefficients of the Maass waveforms.
    INPUT:

    - ``R``   -- double (eigenvalue)
    - ``Y``   -- double (the height of the sampling horocycle)
    - ``Ms,Mf``  -- int (The coefficients we want to compute are C(n), Ms<=n<=Mf )
    - ``Qs,Qf``  -- int (The sampling points are X_m, Qs<=m<=Qf)
    - ``alphas`` -- [nc] double array (the shifts of the Fourier expansion at each cusp)
    - ``V``   -- [(Mf-Ms)*nc]^2 double complex matrix (allocated)
    - ``Xm``  -- [Qf-Qs] double array (allocated)
    - ``Xpb`` -- nc*nc*[Qf-Qs] double array (allocated)
    - ``Ypb`` -- nc*nc*[Qf-Qs] double array (allocated)
    - ``Cvec`` -- nc*nc*[Qf-Qs] double complex array (allocated)
    - `` cuspidal`` -- int (set to 1 if we compute cuspidal functions, otherwise zero)
    - `` sym_type`` -- int (set to 0/1 if we compute even/odd functions, otherwise -1)
    - ``verbose`` -- int (verbosity of output)

    """
    cdef int l,j,icusp,jcusp,n,ni,lj,Ml,Ql,s,Qs,Qf,Mf,Ms
    cdef double pi,sqrtY,Y2pi,nrY2pi,argm,argpb,twopi,two,kbes,besarg,lr,nr
    cdef double complex ckbes,ctmpV,iargm,twopii,ctmp
    if not cuspidal in [0,1]:
        raise ValueError," parameter cuspidal must be 0 or 1"
    if R < 0:
        ## In this case (corresponding to lambda in [0,1/4] we use the real parameter K-Bessel
        R = -R
        set_pref = -1
    else:
        set_pref = 1
    pi=M_PI 
    sqrtY=sqrt(Y)
    two=<double>(2)
    Y2pi=Y*pi*two
    twopi=two*pi
    Ml=0; Ql=0
    for i in range(nc):
        if Mv[i][2]>Ml:
            Ml=Mv[i][2]
        if Qv[i][2]>Ql:
            Ql=Qv[i][2]
        if verbose>0:
            print "Qv[",i,"]=(",Qv[i][0],",",Qv[i][1],",",Qv[i][2],")"
            print "Mv[",i,"]=(",Mv[i][0],",",Mv[i][1],",",Mv[i][2],")"
            print "Qfak[{0}]={1}".format(i,Qfak[i])
            print "symmetric_cusps[{0}]={1}".format(i,symmetric_cusps[i])
    if verbose>0:
        print "N1=",N1
        print "Ql=",Ql
    ## This is the effective offset at the
    cdef int* cusp_offsets=NULL
    cusp_offsets=<int*>sig_malloc(sizeof(int)*nc)
    if cusp_offsets==NULL: raise MemoryError
    for jcusp in range(nc):
        cusp_offsets[jcusp]=0
        for icusp in range(jcusp):
            if icusp==0 or cusp_evs[icusp]==0:
                cusp_offsets[jcusp]+=Mv[icusp][2]
        if verbose>0:
            print "cusp_offsets[",jcusp,"]=",cusp_offsets[jcusp]
    cdef int nc_sym=0
    for jcusp in range(nc):
        if verbose>0:
            print "cusp_evs[",jcusp,"]=",cusp_evs[jcusp]
        if jcusp==0 or cusp_evs[jcusp]<>0:
            nc_sym+=1
    cdef double **nvec=NULL
    nvec = <double**>sig_malloc(sizeof(double*)*nc)
    if not nvec: raise MemoryError
    for icusp in range(nc):
        nvec[icusp] = <double*>sig_malloc(sizeof(double)*Ml)
    cdef double complex ****ef1=NULL
    cdef double complex ***ef2_c=NULL
    ef2_c = <double complex***>sig_malloc(sizeof(double complex**)*nc)
    if not ef2_c: raise MemoryError
    for icusp in range(nc):
        ef2_c[icusp] = <double complex**>sig_malloc(sizeof(double complex*)*Mv[icusp][2])
        for n in range(Mv[icusp][2]):
            ef2_c[icusp][n] = <double complex*>sig_malloc(sizeof(double complex)*Qv[icusp][2])
    ef1 = <double complex****>sig_malloc(sizeof(double complex***)*nc)
    if ef1==NULL: raise MemoryError
    for icusp in range(nc):
        ef1[icusp] = <double complex***>sig_malloc(sizeof(double complex**)*nc)
        if ef1[icusp]==NULL: raise MemoryError
        for jcusp in range(nc):
            ef1[icusp][jcusp] = <double complex**>sig_malloc(sizeof(double complex*)*Mv[jcusp][2])
            if ef1[icusp][jcusp]==NULL: raise MemoryError
            for n in range(Mv[jcusp][2]):
                ef1[icusp][jcusp][n] = <double complex*>sig_malloc(sizeof(double complex)*Qv[jcusp][2])
                if ef1[icusp][jcusp][n]==NULL: raise MemoryError
    for jcusp in range(nc):
        for n in range(Ml):
            nvec[jcusp][n]=<double>(n+Mv[jcusp][0])+alphas[jcusp]
    cdef int twoQm1
    twoQm1= 2*Qv[0][1]-1
    for jcusp in range(nc):
        for n in range(Mv[jcusp][2]):
            for j in range(Qv[jcusp][2]):
                if j < Qv[0][1]:
                    argm=nvec[jcusp][n]*Xm[j]
                else:                    
                    argm=-nvec[jcusp][n]*Xm[2*Qv[0][1]-1-j]
                if symmetric_cusps[jcusp]==0:
                    ef2_c[jcusp][n][j]=cos(argm)
                elif symmetric_cusps[jcusp]==1:
                    ef2_c[jcusp][n][j]=_I*sin(-argm)
                else:
                    ef2_c[jcusp][n][j]=cexpi(-argm)

    cdef double argpb1
    for jcusp in range(nc):
        for icusp in range(nc):
            for n in range(Mv[jcusp][2]):
                for j in range(Qv[jcusp][2]): #in range(Qs,Qf+1):
                    if Ypb[icusp][jcusp][j]==0: #not Xpb.has_key((icusp,jcusp,j):
                        continue
                    argpb=nvec[jcusp][n]*Xpb[icusp][jcusp][j]
                    if symmetric_cusps[jcusp]==0:
                        ef1[icusp][jcusp][n][j]=cos(argpb)
                    elif symmetric_cusps[jcusp]==1:
                        ef1[icusp][jcusp][n][j]=_I*sin(argpb)
                    else:
                        ef1[icusp][jcusp][n][j]=cexpi(argpb)
                    ctmp = Cvec[icusp][jcusp][j]
                    ef1[icusp][jcusp][n][j]=ef1[icusp][jcusp][n][j]*ctmp

    if verbose>1:
        print "here1121"
    cdef double besarg_old=0.0
    cdef double y,kbes_old=1.0
    cdef double ***kbesvec=NULL
    kbesvec=<double***>sig_malloc(sizeof(double**)*nc)
    if kbesvec==NULL:
        raise MemoryError
    for jcusp in range(nc):
        kbesvec[jcusp]=<double**>sig_malloc(sizeof(double*)*Ml)
        if kbesvec[jcusp]==NULL:
            raise MemoryError
        for l in range(Ml):
            kbesvec[jcusp][l]=<double*>sig_malloc(sizeof(double)*Ql) #Qv[jcusp][2])
            if kbesvec[jcusp][l]==NULL:
                raise MemoryError
    if verbose>0:
        print "here0"
        print "Ml=",Ml
        print "Ql=",Ql
    cdef double tmpr
    cdef double besprec
    besprec=1.0E-14
    ## Can we somehow use that "most" of Xpb[icusp][jcusp][2Q-j-1]=-Xpb[icusp][jcusp][j] ?
    ## uncomment the below lines to see this.
    # for jcusp in range(nc):
    #     for icusp in range(nc):
    #         for j from 0<=j<Qv[1][1]:
    #             print "diff[",jcusp,icusp,j,"]=",Xpb[icusp][jcusp][j]+Xpb[icusp][jcusp][2*Qv[jcusp][1]-j-1]

    for jcusp in range(nc):
        for icusp in range(nc):
            if icusp>0 and cusp_evs[icusp]<>0:
                continue
            for j in range(Qv[jcusp][2]):
                if Ypb[icusp][jcusp][j]==0:
                    continue
                for l in range(Mv[jcusp][2]):
                    lr=nvec[jcusp][l]*twopi
                    Mf = Mv[icusp][1]
                    besarg=fabs(lr)*Ypb[icusp][jcusp][j]
                    if lr<>0.0:
                        #besselk_dp_c(&tmpr,R,besarg,besprec,pref=set_pref)
                        if is_exceptional==0:
                            tmpr = my_kbes(R,besarg,besprec,pref=set_pref)
                        else:
                            tmpr = scipy.special.kv(R,besarg)
                        kbesvec[icusp][l][j]=sqrt(Ypb[icusp][jcusp][j])*tmpr
                    else:
                        kbesvec[icusp][l][j]=<double>1.0

    cdef double complex cuspev

    for jcusp in range(nc):
        for l in range(Mv[jcusp][2]):
            lr=nvec[jcusp][l]*twopi
            lj=cusp_offsets[jcusp]+l
            if jcusp>0 and cusp_evs[jcusp]<>0:
                lj0=l; cuspev=cusp_evs[jcusp]
            else:
                lj0=lj;  cuspev=1.0
            if lr==0.0 and cuspidal==1:
                continue
            for j in range(Qv[jcusp][2]):
                for icusp in range(nc):
                    if icusp>0 and cusp_evs[icusp]<>0:
                        continue
                    if Ypb[icusp][jcusp][j]==0:
                        continue
                    ckbes=kbesvec[icusp][l][j]*ef1[icusp][jcusp][l][j]
                    for n in range(Mv[icusp][2]):
                        if nvec[icusp][n]==0 and cuspidal==1:
                            continue
                        #if n+Mv[icusp][0]==0 and cuspidal==1:
                        #    continue
                        ni=cusp_offsets[icusp]+n
                        ctmpV=ckbes*ef2_c[icusp][n][j]
                        if ni>N1 or lj0>N1:
                            raise ArithmeticError,"Index outside!"
                        V[ni][lj0]=V[ni][lj0]+ctmpV*cuspev
                        # if verbose>0 and lj0==324 and ni==46:
                        if verbose>2 and ni==3 and lj0==3:
                            print "l,jcusp,n,icusp=",l+Mv[jcusp][0],jcusp,n+Mv[icusp][0],icusp
                            print "Ypb(",j+Qv[jcusp][0],")=",Ypb[icusp][jcusp][j]
                            print "Cvec=",Cvec[icusp][jcusp][j]
                            print "V[",ni,"][",lj0,"](",j,")=",V[ni][lj0]
                            print "ctmpV[",ni,"][",lj0,"]=",ctmpV
                            print "kbesvec[",icusp,l,j,"]=",kbesvec[icusp][l][j]
                            print "besarg=",fabs(lr)*Ypb[icusp][jcusp][j]
                            print "ef1=",ef1[icusp][jcusp][l][j]
                            print "ef2=",ef2_c[icusp][n][j]
                            print "V[",ni,"][",lj0,"]=",V[ni][lj0],"+",ctmpV,"*",cuspev

    if verbose>2:
    #    #for n from 30+cusp_offsets[3]<=n<40+cusp_offsets[3]:
    #    for n from 43<=n<49:
    #        print "V0[",0,n,"]=",V[0][n]
    #        #print "V0[",24,24,"]=",V[24][24]
        print "V0[",3,3,"]=",V[3][3]
        for icusp in range(nc):
            print "Qfak[{0}]={1}".format(icusp,Qfak[icusp])
    for icusp in range(nc):
        if icusp>0 and cusp_evs[icusp]<>0:
            continue
        for n in range(Mv[icusp][2]):
            ni=cusp_offsets[icusp]+n
            for jcusp in range(nc):
                if jcusp>0 and cusp_evs[jcusp]<>0:
                    continue
                for l in range(Mv[jcusp][2]):
                    lj=cusp_offsets[jcusp]+l
                    if lj>N1 or ni>N1: # some extra insurance...
                        print "ni=",cusp_offsets[icusp],"+",n,"=",ni
                        print "lj=",cusp_offsets[jcusp],"+",l,"=",lj
                        raise ArithmeticError,"Index outside!"
                    if verbose>2 and ni==3 and lj==3:
                        print "V[3][3]={0}/{1}".format(V[ni][lj],Qfak[jcusp])
                    V[ni][lj]=V[ni][lj]/<double complex>Qfak[jcusp]
                    if verbose>2 and ni==3 and lj==3:
                        print "V[3][3]/Qfak[{0}]={1}".format(jcusp,V[ni][lj])

    for icusp in range(nc):
        if icusp>0 and cusp_evs[icusp]<>0:
            continue
        for n in range(Mv[icusp][2]):
            nr=fabs(nvec[icusp][n])
            if nr==0.0 and cuspidal==1:
                continue
            ni=cusp_offsets[icusp]+n
            nrY2pi=nr*Y2pi
            if nr==0.0:
                kbes=<double>1.0
            else:
                #mpIR=mpmath.fp.mpc(0,R)
                #                kbes=float(mpmath.fp.besselk(mpIR,nrY2pi).real*exp(mpmath.fp.pi*R*0.5))
                #besselk_dp_c(&kbes,R,nrY2pi,besprec,pref=set_pref)
                if is_exceptional==0:
                    kbes = my_kbes(R,nrY2pi,besprec,pref=set_pref)
                else:
                    kbes = scipy.special.kv(R,nrY2pi)
                kbes=sqrtY*kbes # besselk_dp(R,nrY2pi,pref=1)
                if verbose>2 and abs(nr)<= 3:
                   print "nr=",nr
                   print "nrY2pi=",nrY2pi
                   print "R=",R
                   print "kbes=",kbes
            if ni>N1:
                raise ArithmeticError,"Index outside!"
            if verbose>2 and ni==3:
                print "V[3][3]={0}-{1}={2}".format(V[ni][ni],kbes,V[ni][ni]-kbes)
            V[ni][ni]=V[ni][ni] - kbes

    if verbose>3:

        #        #for n from 30+cusp_offsets[3]<=n<40+cusp_offsets[3]:
        #        print "V3[",0,n,"]=",V[0][n]
        for n in range(min(100,Ml)):
            print "V3[",n,n,"]=",V[n][n]
    #if Qfak<>NULL:
    #    #for icusp in range(nc):
    #    #    if Qfak[icusp]<>NULL:
    #    #        sig_free(Qfak[icusp])
    #    sig_free(Qfak)
    if kbesvec<>NULL:
        for icusp in range(nc):
            if kbesvec[icusp]<>NULL:
                for l in range(Ml):
                    if kbesvec[icusp][l]<>NULL:
                        sig_free(kbesvec[icusp][l])
                sig_free(kbesvec[icusp])
        sig_free(kbesvec)
    #print "deal kbbes1"
    if ef1<>NULL:
        for jcusp in range(nc):
            if ef1[jcusp]<>NULL:
                for icusp in range(nc):
                    if ef1[jcusp][icusp]<>NULL:
                        for n in range(Mv[icusp][2]):
                            if ef1[jcusp][icusp][n]<>NULL:
                                sig_free(ef1[jcusp][icusp][n])
                        sig_free(ef1[jcusp][icusp])
                sig_free(ef1[jcusp])
        sig_free(ef1)
    if ef2_c<>NULL:
        for icusp in range(nc):
            if ef2_c[icusp]<>NULL:
                for n in range(Mv[icusp][2]):
                    if ef2_c[icusp][n]<>NULL:
                        sig_free(ef2_c[icusp][n])
                sig_free(ef2_c[icusp])
        sig_free(ef2_c)
    if nvec<>NULL:
        for icusp in range(nc):
            if nvec[icusp]<>NULL:
                sig_free(nvec[icusp])
        sig_free(nvec)
    if cusp_offsets<>NULL:
        sig_free(cusp_offsets)






cdef compute_V_real_dp(double **V,double R,double Y,int Ms,int Mf,int Qs,int Qf,int nc, double *alphas, double *Xm,double ***Xpb,double ***Ypb, double ***Cvec, int cuspidal=1,int sym_type=-1,int verbose=0):
    r"""
    Set up the matrix for the system of equations giving the Fourier coefficients of the Maass waveforms.s
    """
    cdef int l,j,icusp,jcusp,n,ni,lj,Ml,Ql,s
    cdef double sqrtY,Y2pi,nrY2pi,argm,argpb,twopi,two,kbes,Qfak
    cdef double ckbes,ctmpV,iargm,twopii,ctmp,lr,nr,besarg,pi
    cdef int set_pref
    if R < 0:
        ## In this case (corresponding to lambda in [0,1/4] we use the real parameter K-Bessel
        R = -R
        set_pref = -1
    else:
        set_pref = 1
    if not cuspidal in [0,1]:
        raise ValueError," parameter cuspidal must be 0 or 1"

    pi=M_PI #<double>RealField(53).pi() # 3.14159265358979323846264338328
    sqrtY=sqrt(Y)
    two=<double>(2)
    Y2pi=Y*pi*two
    twopi=two*pi
    #cdef int Ql,Ml
    Ql=Qf-Qs+1
    Ml=Mf-Ms+1
    if sym_type in [0,1]:
        Qfak=<double>(Ql)/<double>(2)
    else:
        raise ValueError,"Thius routine must be called with a completely symmetrised space so we get all parameters real!"
    #s=int(nc*Ml)
    #cdef nvec=np.arange(Ms,Mf+1,dtype=DTYPE)
    cdef double **nvec=NULL
    nvec = <double**>sig_malloc(sizeof(double*)*Ml)
    if nvec ==NULL: raise MemoryError
    for n in range(Ml):
        nvec[n] = NULL
        nvec[n] = <double*>sig_malloc(sizeof(double)*nc)
        if not nvec[n]: raise MemoryError
    cdef double **** ef1
    cdef double ***ef2=NULL
    # Here sym_type must be in 0,1
    if sym_type not in [0,1]:
        raise ValueError,"Need even or odd symmetry here!"
    ef2 = <double***>sig_malloc(sizeof(double**)*Ml)
    if not ef2: raise MemoryError
    for n in range(Ml):
        ef2[n] = <double**>sig_malloc(sizeof(double*)*nc)
        if not ef2[n]: raise MemoryError
        for jcusp in range(nc):
            ef2[n][jcusp] = <double*>sig_malloc(sizeof(double)*Ql)
            if not ef2[n][jcusp]: raise MemoryError
    #if verbose>0:
    #    print "Ml=",Ml
    #    print "Ql=",Ql
    ef1 = <double****>sig_malloc(sizeof(double***)*Ml)
    if ef1==NULL: raise MemoryError
    for n in range(Ml):
        ef1[n] = <double***>sig_malloc(sizeof(double**)*nc)
        if ef1[n]==NULL: raise MemoryError
        for jcusp in range(nc):
            ef1[n][jcusp] = <double**>sig_malloc(sizeof(double*)*nc)
            if ef1[n][jcusp]==NULL: raise MemoryError
            for icusp in range(nc):
                ef1[n][jcusp][icusp] = <double*>sig_malloc(sizeof(double)*Ql)
                if ef1[n][jcusp][icusp]==NULL: raise MemoryError

    for n in range(Ml): #in range(Ms,Mf+1):
        for jcusp in range(nc):
            nvec[n][jcusp]=<double>(n+Ms)+alphas[jcusp]
            for j in range(Ql): #in range(Qs,Qf+1):
                argm=nvec[n][jcusp]*Xm[j]
                if sym_type==0:
                    ef2[n][jcusp][j]=cos(argm)
                elif sym_type==1:
                    ef2[n][jcusp][j]=sin(argm)
    for n in range(Ml): #in range(Ms,Mf+1):
        #print "ef2[",n,j,"=",ef2_r[n][j]
        for jcusp in range(nc):
            for icusp in range(nc):
                for j in range(Ql): #in range(Qs,Qf+1):
                    if Ypb[icusp][jcusp][j]==0: #not Xpb.has_key((icusp,jcusp,j))
                        continue
                    argpb=nvec[n][jcusp]*Xpb[icusp][jcusp][j]
                    if sym_type==0:
                        ef1[n][icusp][jcusp][j]=cos(argpb)
                    elif sym_type==1:
                        ef1[n][icusp][jcusp][j]=sin(argpb)
                    #print "Cv=",Cv[icusp,jcusp,j],type(Cv[icusp,jcusp,j])
                    ctmp = Cvec[icusp][jcusp][j]
                    ef1[n][icusp][jcusp][j]=ef1[n][icusp][jcusp][j]*<DTYPE_t>ctmp
                    #print "ef1[",n-Ms,"][",icusp,"][",jcusp,"][",j-Qs,"]",ef1[j-Qs,icusp,jcusp,n-Ms]
    #if verbose>0:
    #                print "here!"
    for l in range(Ml):
        for n in range(Ml):
            V[n][l]=0.0
    for l in range(Ml):
        for jcusp in range(nc):
            lr=nvec[l][jcusp]*twopi
            if lr==0.0 and cuspidal==1:
                continue
            lj=Ml*jcusp+l
            for icusp in range(nc):
                for j in range(Ql):
                    if Ypb[icusp][jcusp][j]==0:
                        #if(not Ypb.has_key((icusp,jcusp,j))):
                        continue
                    #ypb=Ypb[icusp][jcusp][j]
                    besarg=abs(lr)*Ypb[icusp][jcusp][j]
                    if lr<>0.0:
                        kbes=sqrt(Ypb[icusp][jcusp][j])*my_kbes(R,besarg,pref=set_pref)
                    else:
                        kbes=<double>1.0
                    ckbes=kbes*ef1[l][icusp][jcusp][j]
                    #if verbose>2 and l==1:
                    #    print "kbes=",kbes
                    #print "ef1=",ef1[l][icusp][jcusp][j]
                    for n in range(Ml):
                        if nvec[n][icusp]==0 and cuspidal==1:
                            continue
                        ni=Ml*icusp+n
                        ctmpV=ckbes*ef2[n][jcusp][j]
                        V[ni][lj]=V[ni][lj]+ctmpV
                        # <CTYPE_t>creal(ctmpV)+_I*cimag(ctmpV)
                        if verbose>1 and lj==1 and (ni==1 or ni==12):
                            print "V[ni][lj]=",V[ni][lj]
                            print "Ypb[",icusp,jcusp,j,"]=",Ypb[icusp][jcusp][j]
                            print "ef1=",ef1[l][icusp][jcusp][j]

    for n in range(Ml*nc):
        for l in range(Ml*nc):
            V[n][l]=V[n][l]/Qfak
    for n in range(Ml):
        #for icusp in range(nc):
        for icusp in range(nc):
            nr=abs(nvec[n][icusp])
            if nr==0.0:
                if cuspidal==1:
                    continue
                else:
                    kbes=<double>1.0
            else:
                nrY2pi=nr*Y2pi
                kbes=sqrtY*my_kbes(R,nrY2pi,pref=set_pref)
            ni=Ml*icusp+n
            V[ni][ni]=V[ni][ni] - kbes
    if ef2<>NULL:
        for n in range(Ml):
            if ef2[n]<>NULL:
                for jcusp in range(nc):
                    if ef2[n][jcusp]<>NULL:
                        sig_free(ef2[n][jcusp])
                sig_free(ef2[n])
        sig_free(ef2)
    if ef1<>NULL:
        for n in range(Ml):
            if ef1[n] <> NULL:
                for jcusp in range(nc):
                    if ef1[n][jcusp]<>NULL:
                        for icusp in range(nc):
                            if ef1[n][jcusp][icusp]<>NULL:
                                sig_free(ef1[n][jcusp][icusp])
                        sig_free(ef1[n][jcusp])
                sig_free(ef1[n])
        sig_free(ef1)
    #my_dealloc_double_2(&nvec,Ml)
    #printf("nvec=%p \n",nvec)
    #print "hej nv=",nvec[0][0]
    if nvec<>NULL:
        for n in range(Ml):
            if nvec[n]<>NULL:
                sig_free(nvec[n])
        #print "End!"
        if nvec<>NULL:
            #printf("nvec=%p \n",nvec)
            sig_free(nvec)



cdef normalize_matrix_cplx_dp(double complex **V,int N,int comp_dim,int num_set,int *setc_list,double complex ** vals_list,int verbose=0):
    r"""
    W['V']
    an efficient matrix normalization....
    Assume that the right hand side is in the last column of the matrix
    """
    cdef int nrows,ncols
    cdef int i,j,do_cont,k,r,coffs,roffs
    if verbose>1:
        print "comp_dim=",comp_dim
        print "num_set=",num_set
        for i in range(num_set):
            print "setc_list=",setc_list[i]
            for j in range(comp_dim):
                print "vals_list=",vals_list[j][i]
        print "N=",N
    cdef double complex* tmp
    tmp=<double complex*>sig_malloc(sizeof(double complex)*comp_dim)
    roffs=0
    for r in range(N):
        do_cont=0
        for i in range(num_set):
            if setc_list[i]==r:
                do_cont=1
                break
        if do_cont==1:
            roffs=roffs+1
            continue
        coffs=0
        for j in range(comp_dim):
            tmp[j]=0
        for k in range(N):
            do_cont=0
            for i in range(num_set):
                if setc_list[i]==k:
                    for j in range(comp_dim):
                        tmp[j] = tmp[j] - vals_list[j][i]*V[r][k]
                    do_cont=1
                    break
            if do_cont==1:
                coffs=coffs+1
                continue
            V[r-roffs][k-coffs]=V[r][k]
        for j in range(comp_dim):
            V[r-roffs][N+j-coffs]=tmp[j]

    if tmp<>NULL:
        sig_free(tmp)

cdef normalize_matrix_real_dp(double **V,int N,int comp_dim,int num_set,int *setc_list,double ** vals_list,int verbose=0):
    r"""
    W['V']
    an efficient matrix normalization....
    Assume that the right hand side is in the last column of the matrix
    """
    cdef int nrows,ncols
    cdef int i,j,do_cont,k,r,coffs,roffs
    if verbose>0:
        print "comp_dim={0}".format(comp_dim)
        print "num_set=",num_set
        for i in range(num_set):
            print "setc_list=",setc_list[i]
            for j in range(comp_dim):
                print "vals_list=",vals_list[j][i]
        print "N=",N
    cdef double* tmp
    tmp=<double*>sig_malloc(sizeof(double)*comp_dim)
    roffs=0
    for r in range(N):
        do_cont=0
        for i in range(num_set):
            if setc_list[i]==r:
                do_cont=1
                break
        if do_cont==1:
            roffs=roffs+1
            continue
        coffs=0
        for j in range(comp_dim):
            tmp[j]=<double>0
        for k in range(N):
            do_cont=0
            for i in range(num_set):
                if setc_list[i]==k:
                    for j in range(comp_dim):
                        tmp[j] = tmp[j] - vals_list[j][i]*V[r][k]
                        #if r<=3:
                        #    print "tmpj=tmpj-",vals_list[j][i],"*",V[r][k]
                        #    print "tmp[",r,j,"]=",tmp[j]
                    do_cont=1
                    break
            if do_cont==1:
                coffs=coffs+1
                continue
            V[r-roffs][k-coffs]=V[r][k]
        for j in range(comp_dim):
            V[r-roffs][N+j-coffs]=tmp[j]
    if tmp<>NULL:
        sig_free(tmp)


cdef normalize_matrix_cplx_sym_dp(double complex **V,double complex **V1,int Ml,int N,int comp_dim,int num_set,int *setc_list,double complex ** vals_list,double complex * cusp_evs,int ncusps,int verbose=0):
    r"""
    W['V']
    an efficient matrix normalization....
    Assume that the right hand side is in the last column of the matrix
    V = (nc*M)x(nc*M)
    V1=(nc1*M')x(nc1*M') with nc1=nc-#[cusps which are set] and M'=M-#{coefficients which are set}
    """
    cdef int nrows,ncols
    cdef int i,j,do_cont,k,r,coffs,roffs
    cdef int ncolsV1
    cdef int icusp,jcusp
    ncolsV1=ncusps*Ml
    for i in range(1,ncusps):
        if cusp_evs[i]<>0:
            ncolsV1-=Ml
    ncolsV1=ncolsV1-num_set
    #for i from 0<=i<Ml:
    #    if
    if verbose>2:
        print "comp_dim=",comp_dim
        print "num_set=",num_set
        for i in range(ncusps):
            print "cusp_evs=",cusp_evs[i]
        for i in range(num_set):
            print "setc_list=",setc_list[i]
            for j in range(comp_dim):
                print "vals_list=",vals_list[j][i]
        print "N=",N
    cdef double complex* tmp
    tmp=<double complex*>sig_malloc(sizeof(double complex)*comp_dim)
    if not tmp: raise MemoryError
    roffs=0
    for icusp in range(ncusps):
        if icusp>0 and cusp_evs[icusp]<>0:
            roffs+=Ml
            continue
        for n in range(Ml):
            r=icusp*Ml+n
            #for r from 0<=r<N:
            do_cont=0
            for i in range(num_set):
                if setc_list[i]==r:
                    do_cont=1
                    break
            if do_cont==1:
                roffs=roffs+1
                continue
            for j in range(comp_dim):
                tmp[j]=0
            coffs=0
            for jcusp in range(ncusps):
                if jcusp>0 and cusp_evs[jcusp]<>0:
                    for l in range(Ml):
                        k = jcusp*Ml+l
                        V[r-roffs][k-coffs]+=V[r][k]*cusp_evs[jcusp]
                    coffs+=Ml
                for l in range(Ml):
                    k = jcusp*Ml+l
                    do_cont=0
                    for i in range(num_set):
                        if setc_list[i]==k:
                            for j in range(comp_dim):
                                tmp[j] = tmp[j] - vals_list[j][i]*V[r][k]
                            do_cont=1
                            break
                    if do_cont==1:
                        coffs=coffs+1
                        continue
                    V[r-roffs][k-coffs]=V[r][k]
                for j in range(comp_dim):
                    V[r-roffs][N+j-coffs]=tmp[j]

    if tmp<>NULL:
        sig_free(tmp)


cpdef get_coeff_fast_cplx_dp(S,double R,double Y,int M,int Q,dict Norm={},int gr=0,int norm_c=1,dict cusp_ev={},double eps=1e-12,int do_par=0,int ncpus=1):
        r"""
        Pick the correct method...
        """
        if cusp_ev == {}:
            cusp_ev = Norm.get('cusp_ev',{})
        #if cusp_ev=={} or not S.group().is_Gamma0() or S.weight()<>0:
#        if  not S.group().is_Gamma0() or S.weight()<>0:             
#            res = get_coeff_fast_cplx_dp_nosym(S,R,Y,M,Q,Norm,gr,norm_c,do_par=do_par,ncpus=ncpus)
#        else:
        res = get_coeff_fast_cplx_dp_sym(S,R,Y,M,Q,Norm,gr,norm_c,cusp_ev=cusp_ev,eps=1e-12,do_par=do_par,ncpus=ncpus)
        return res
            

cpdef get_coeff_fast_cplx_dp_sym(S,double R,double Y,int M,int Q,dict Norm={},int gr=0,int norm_c=1,dict cusp_ev={},double eps=1e-12,int do_par=0,int ncpus=1):
    r"""

    An efficient method to get coefficients in the double complex case.
    Trying to use as much symmetries as possible.

    """
    import mpmath
    cdef double complex **V=NULL
    cdef double complex **V1=NULL
    cdef double* Xm=NULL
    cdef double*** Xpb=NULL
    cdef double*** Ypb=NULL
    cdef double complex*** Cvec=NULL
    cdef double complex **C=NULL
    cdef int *setc_list=NULL
    cdef double complex *RHS=NULL
    cdef double complex **vals_list=NULL
    cdef double *alphas=NULL
    cdef int nc,Ql,Qs,Qf,Ml,Ms,Mf,j,k,l,N,sym_type,i,n,r
    cdef int num_rhs=0
    cdef double *Qfak,tmpr
    cdef double complex *cusp_evs=NULL
    cdef int **Mv=NULL
    cdef int **Qv=NULL
    cdef int verbose=S._verbose
    cdef int ncolsV1
    cdef int *symmetric_cusps=NULL
    cdef int N1=0
    cdef list SetCs
    cdef dict Vals
    cdef int comp_dim,num_set
    cdef int ncols,ncols1
    cdef int cuspidal=1
    cdef int q
    cdef double complex *sqch=NULL
    cdef int  is_exceptional = S._exceptional
    tmpr = <double>S._group.minimal_height()
    if Y<= 0 or Y >= tmpr:
        Y = 0.5 * tmpr
    if M<=0:
        M = get_M_for_maass_dp(R,Y,eps)

    if Q<M:
        Q=M+20
    #    Qfak=<double>(2*Q)
    sym_type = S._sym_type
    nc = int(S._group._ncusps)
    Qfak = <double *>sig_malloc(sizeof(double)*nc)
    Mv=<int**>sig_malloc(sizeof(int*)*nc)
    if not Mv: raise MemoryError
    Qv=<int**>sig_malloc(sizeof(int*)*nc)
    if not Qv: raise MemoryError
    for i in range(nc):
        Mv[i]=<int*>sig_malloc(sizeof(int)*3)
        if not Mv[i]: raise MemoryError
        Qv[i]=<int*>sig_malloc(sizeof(int)*3)
        if not Qv[i]: raise MemoryError
    symmetric_cusps=<int*> sig_malloc(sizeof(int)*nc)
    if not symmetric_cusps: raise MemoryError
    cusp_evs=<double complex*>sig_malloc(sizeof(double complex)*nc)
    if not cusp_evs: raise MemoryError
    N = 0; Ml=0; Ql=0
    cdef int* cusp_offsets=NULL
    cusp_offsets=<int*>sig_malloc(sizeof(int)*nc)
    if cusp_offsets==NULL: raise MemoryError
    cdef dict symmetries
    if verbose>0:
        print "INPUT: R={0}, Y={1}, M={2}, Q={3}, Norm={4}, gr={5}, norm_c={6}, cusp_ev={7}, eps={8}".format(<double>R,<double>Y,M,Q, Norm,gr,norm_c,cusp_ev,eps)
    set_Mv_Qv_symm(S,Mv,Qv,Qfak,symmetric_cusps,cusp_evs,cusp_offsets,&N,&Ml,&Ql,M,Q,verbose)
    
    Qs = 1-Q; Qf = Q
    if verbose>0:
        print "N=",N," Ml=",Ml," Ql=",Ql
    if cusp_ev=={}:
        cusp_ev = S.atkin_lehner_eigenvalues()
    for i in range(nc):
        cusp_evs[i]=CC(cusp_ev.get(S.group().cusps()[i],0))
        if i==0 or cusp_evs[i]==0:
            N1=N1+Mv[i][2]
        if verbose>0:
            print "cusp_ev[",i,"]=",cusp_evs[i]
    for jcusp in range(nc):
        cusp_offsets[jcusp]=0
        for icusp in range(jcusp):
            if icusp==0 or cusp_evs[icusp]==0:
                cusp_offsets[jcusp]+=Mv[icusp][2]
                if verbose>0:
                    print "cusp_offset[",jcusp,"]+=",Mv[icusp][2]
            #else:
            #    print "cusp_ev=",cusp_evs[icusp]
        if verbose>0:
            print "cusp_offset[",jcusp,"]=",cusp_offsets[jcusp]

    if verbose>0:
        print "N1=",N1
    if Norm=={}:
        Norm = S.set_norm(1)
    SetCs=Norm['SetCs'][0]
    Vals=Norm['Vals']
    comp_dim=Norm['comp_dim']
    num_set=0
    for r,n in SetCs:
        if r>0 and cusp_evs[r]<>0:
            continue
        num_set+=1
    V=<double complex**>sig_malloc(sizeof(double complex*)*N)
    if V==NULL: raise MemoryError
    V1=<double complex**>sig_malloc(sizeof(double complex*)*N1)
    if V1==NULL: raise MemoryError
    if num_set<comp_dim:
        ncols=N+comp_dim-num_set
        ncols1=N1+comp_dim-num_set
    else:
        ncols=N
        ncols1=N1
    if verbose>0:
        print "In get_coef_cplx_dp_sym R,Y,M,Q=",R,Y,M,Q
    if verbose>0:
        print "N=",N
        print "ncols=",ncols
        print "N1=",N1
        print "ncols1=",ncols1
        print "Vals=",Vals
        print "SetCs=",SetCs
        print "comp_dim=",comp_dim
        print "num_set=",num_set
        for i in range(nc):
            print "Ms,Mf,Ml[",i,"]=",Mv[i][0],Mv[i][1],Mv[i][2]
            print "Qs,Qf,Ql[",i,"]=",Qv[i][0],Qv[i][1],Qv[i][2]
    for j in range(N):
        V[j]=<double complex*>sig_malloc(sizeof(double complex)*(ncols))
        for k in range(ncols):
            V[j][k]=0
    for j in range(N1):
        V1[j]=<double complex*>sig_malloc(sizeof(double complex)*(ncols1))
        for k in range(ncols1):
            V1[j][k]=0
    alphas=<double*>sig_malloc(sizeof(double)*nc)
    if alphas==NULL: raise MemoryError
    #printf("alphas_init=%p \n",alphas)
    Xm=<double*>sig_malloc(sizeof(double)*Qv[0][1])
    if Xm==NULL: raise MemoryError
    Xpb = <double***> sig_malloc( sizeof(double** ) * nc )
    if Xpb==NULL: raise MemoryError
    Ypb = <double***> sig_malloc( sizeof(double** ) * nc )
    if Ypb==NULL: raise MemoryError
    for i in range(nc):
        Xpb[i]=NULL; Ypb[i]=NULL
        Xpb[i] = <double**>sig_malloc(sizeof(double*) * nc )
        Ypb[i] = <double**>sig_malloc(sizeof(double*) * nc )
        if Ypb[i]==NULL or Xpb[i]==NULL:
            raise MemoryError
        for j in range(nc):
            Xpb[i][j]=NULL; Ypb[i][j]=NULL
            Xpb[i][j] = <double*>sig_malloc(sizeof(double) * Ql )
            Ypb[i][j] = <double*>sig_malloc(sizeof(double) * Ql )
            if Ypb[i][j]==NULL or Xpb[i][j]==NULL:
                raise MemoryError
            for n in range(Ql):
                Xpb[i][j][n]=<double>0
                Ypb[i][j][n]=<double>0
    Cvec = <double complex***>sig_malloc(sizeof(double complex**) * nc )
    if Cvec==NULL: raise MemoryError
    for i in range(nc):
        Cvec[i] = <double complex**>sig_malloc(sizeof(double complex*) * nc )
        if Cvec[i]==NULL:
            raise MemoryError
        for j in range(nc):
            Cvec[i][j] = <double complex*>sig_malloc(sizeof(double complex) * Ql )
            if Cvec[i][j]==NULL:
                raise MemoryError
            for n in range(Ql):
                Cvec[i][j][n]=0
#    sig_on()
    cdef int t=0
    t = pullback_pts_cplx_dp_sym(S,Qv,Y,Xm,Xpb,Ypb,Cvec)
    if t==1:
        raise ArithmeticError,"Need smaller Y than {0}".format(Y)
#    sig_off()
    for j in range(nc):
        tmpr=float(S.alpha(j)[0])
        alphas[j]=<double>tmpr
        if verbose>0:
            print "alphas[",j,"]=",alphas[j],type(alphas[j])
#    sig_on()
    if do_par==1 and ncpus>=1:
        compute_V_cplx_dp_sym_par(V1,N1,Xm,Xpb,Ypb,Cvec,
                              cusp_evs,alphas,Mv,Qv,Qfak,
                              symmetric_cusps,
                              R,Y,nc,ncols,cuspidal,verbose,ncpus,
                                  is_exceptional=is_exceptional,
                                      is_trivial=0)                                      
    else:
        compute_V_cplx_dp_sym(V1,N1,Xm,Xpb,Ypb,Cvec,
                              cusp_evs,alphas,Mv,Qv,Qfak,
                              symmetric_cusps,
                              R,Y,nc,ncols,cuspidal,verbose,
                              is_exceptional=is_exceptional
                            )
#    sig_off()
    cdef Matrix_complex_dense VV
    #Vtmp = load("A.sobj")
    #for l from 0<=l<N1:
    #    for n from 0<=n<N1:
    #        V1[l][n]=<double complex>complex(Vtmp[l][n])
    if Qfak<>NULL:
        sig_free(Qfak)
    if alphas<>NULL:
        sig_free(alphas)
    if Xm<>NULL:
        sig_free(Xm)
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
    if sqch<>NULL:
        sig_free(sqch)
    if gr==1 or gr==4:
        CF=MPComplexField(53)
        MS=MatrixSpace(CF,N1,ncols1)
        VV=Matrix_complex_dense(MS,0)
        for n in range(N1):
            for l in range(ncols1):
                VV[n,l]=V1[n][l]
        #return VV
    C=<double complex**>sig_malloc(sizeof(double complex*)*(comp_dim))
    if C==NULL:
        raise MemoryError
    for i in range(comp_dim):
        C[i]=<double complex*>sig_malloc(sizeof(double complex)*(N1))
        if C[i]==NULL:
            raise MemoryError
    setc_list = <int*>sig_malloc(sizeof(int)*num_set)
    if setc_list==NULL:
        raise MemoryError
    vals_list = <double complex** >sig_malloc(sizeof(double complex*)*comp_dim)
    if vals_list==NULL:
        raise MemoryError
    for j in range(comp_dim):
        vals_list[j]=<double complex*>sig_malloc(sizeof(double complex)*num_set)
    i=0
    for r,n in SetCs:
        if cusp_evs[r]<>0 and r>0:
            continue
        l = cusp_offsets[r]+n-Mv[r][0]
        if l<0:
            continue
        setc_list[i]=l
        for j in range(comp_dim):
            vals_list[j][i]=<double complex>CC(Vals[j][(r,n)])
        i=i+1
    if num_rhs>0 and num_rhs<>comp_dim:
        raise ValueError,"Need same number of right hand sides (or just one) as the number of set coefficients!"
    if verbose>0:
        print "comp_dim=",comp_dim
    normalize_matrix_cplx_dp(V1,N1,comp_dim,num_set,setc_list,vals_list,verbose)
    if gr==2 or gr==5:
        CF=MPComplexField(53)
        MS=MatrixSpace(CF,N1,ncols1)
        VV=Matrix_complex_dense(MS,0)
        for n in range(N1):
            for l in range(ncols1):
                VV[n,l]=V1[n][l]
        #return VV
    sig_on()
    if verbose>0:
        print "comp_dim=",comp_dim
    if do_par<=1:
        SMAT_cplx_dp(V1,ncols1-num_set,comp_dim,num_set,C,vals_list,setc_list)
        if gr==5:
            ## check result...
            CF=MPComplexField(53)
            MS=MatrixSpace(CF,ncols1-num_set,ncols1-num_set)
            VVV=Matrix_complex_dense(MS,0)
            MSB=MatrixSpace(CF,ncols1-num_set,comp_dim)
            B=Matrix_complex_dense(MSB,0)
            MSC=MatrixSpace(CF,comp_dim,ncols1-num_set)
            C1=Matrix_complex_dense(MSC,0)
            for n in range(ncols1-num_set):
                for l in range(ncols1-num_set):
                    VVV[n,l]=VV[n,l]
                for l in range(comp_dim):
                    C1[l,n]=C[l][n]
                for l in range(comp_dim):
                    B[n,l]=VV[n,ncols1-num_set+l]
            print "return here"
            return VVV,C1,B
    else:
        if verbose>0:
            print "smat parallel!"
        SMAT_cplx_par_dp(V1,ncols1-num_set,comp_dim,num_set,C,vals_list,setc_list,ncpus)
    sig_off()
    if verbose>1:
        for k in range(ncols1-num_set):
            print "C[",k,"]=",C[0][k]
    if gr==3 or gr==4:
        CF=MPComplexField(53)
        MSC=MatrixSpace(CF,comp_dim,N1)
        C1=Matrix_complex_dense(MSC,0)
        for n in range(comp_dim):
            for l in range(N1):
                C1[n,l]=C[n][l]
    cdef dict res={}
    cdef int ki
    if gr==0:
        for j in range(comp_dim):
            res[j]=dict()
            for i in range(nc):
                if i>0 and cusp_evs[i]<>0:
                    res[j][i]=cusp_evs[i]
                    continue
                res[j][i]=dict()
                for k in range(Mv[i][2]):
                    ki=cusp_offsets[i]+k
                    res[j][i][k+Mv[i][0]]=C[j][ki]

    if setc_list<>NULL:
        sig_free(setc_list)
    if vals_list<>NULL:
        for j in range(comp_dim):
            if vals_list[j]<>NULL:
                sig_free(vals_list[j])
        sig_free(vals_list)
    if V<>NULL:
        for j in range(N):
            if V[j]<>NULL:
                sig_free(V[j])
        sig_free(V)
    if V1<>NULL:
        for j in range(N1):
            if V1[j]<>NULL:
                sig_free(V1[j])
        sig_free(V1)

    if C<>NULL:
        for i in range(comp_dim):
            if C[i]<>NULL:
                sig_free(C[i])
        sig_free(C)
    if Mv<>NULL:
        sig_free(Mv)
    if Qv<>NULL:
        sig_free(Qv)
    if symmetric_cusps<>NULL:
        sig_free(symmetric_cusps)
    if cusp_evs<>NULL:
        sig_free(cusp_evs)
    #if comp_dim==1:
    #    return res[0]
    if gr==1:
        return VV
    elif gr==2:
        return VV
    elif gr==3:
        return C1
    elif gr==4:
        return VV,C1
    return res

# cpdef get_coeff_fast_cplx_dp_nosym(S,double R,double Y,int M,int Q,dict Norm={},int gr=0,int norm_c=1,int do_par=0,int ncpus=1):
#     r"""
#     An efficient method to get coefficients in the double complex case.
#     """
#     import mpmath
#     cdef double complex **V=NULL
#     cdef double* Xm=NULL
#     cdef double*** Xpb=NULL
#     cdef double*** Ypb=NULL
#     cdef double complex*** Cvec=NULL
#     cdef double complex **C=NULL
#     cdef int *setc_list=NULL
#     cdef double complex *RHS=NULL
#     cdef double complex **vals_list=NULL
#     cdef int nc,Ql,Qs,Qf,Ml,Ms,Mf,j,k,l,N,sym_type,i,n,r
#     cdef int num_rhs=0
#     cdef double weight
#     cdef double *Qfak
#     cdef int verbose=S._verbose
#     cdef int  is_exceptional = S._exceptional
#     weight = <double>RealField(53)(S.weight())
#     if Q<M:
#         Q=M+20
#     sym_type = S._sym_type
#     if sym_type == 1 or sym_type == 0:
#         Ms=0;  Mf=M; Qs=1; Qf=Q
#     else:
#         Ms=-M;  Mf=M; Qs=1-Q; Qf=Q
#     Ml=Mf-Ms+1
#     Ql=Qf-Qs+1
#     for i in range(nc):
#         if sym_type == 0 or sym_type==1:
#             Qfak[i] =<double>(Q)/<double>(2)
#         else:
#             Qfak[i] =<double>(2*Q)    

#     #Qfak = <double*>sig_malloc(sizeof(double)*nc)
#     # if not Qfak: raise MemoryError    
#     # symmetric_cusps=<int*> sig_malloc(sizeof(int)*nc)
#     # if not symmetric_cusps: raise MemoryError
#     # cusp_evs=<double complex*>sig_malloc(sizeof(double complex)*nc)
#     # if not cusp_evs: raise MemoryError
#     # cusp_offsets=<int*>sig_malloc(sizeof(int)*nc)
#     # if cusp_offsets==NULL: raise MemoryError
#     # set_Mv_Qv_symm(S,Mv,Qv,Qfak,symmetric_cusps,cusp_evs,cusp_offsets,&N,&Ml,&Ql,M,Q,verbose)
#     nc = S._group._ncusps
#     N = nc*Ml
#     if Norm=={}:
#         Norm = S.set_norm(1)
#     cdef list SetCs
#     cdef dict Vals
#     SetCs=Norm['SetCs'][0]
#     Vals=Norm['Vals']
#     cdef int comp_dim,num_set
#     comp_dim=Norm['comp_dim']
#     num_set=len(SetCs)
#     V=<double complex**>sig_malloc(sizeof(double complex*)*N)
#     if V==NULL: raise MemoryError
#     cdef int ncols
#     if num_set<comp_dim:
#         ncols=N+comp_dim-num_set
#     else:
#         ncols=N
#     if verbose>1:
#         print "In get_coef_cplx_dp_nosym R=",R, "M=",M,"Y=",Y," Q=",Q,'do_par=',do_par,'ncpus=',ncpus
#     if verbose>2:
#         print "N=",N
#         print "Vals=",Vals
#         print "SetCs=",SetCs
#         print "comp_dim=",comp_dim
#         print "num_set=",num_set
#         print "ncols=",ncols
#     for j in range(N):
#         V[j]=<double complex*>sig_malloc(sizeof(double complex)*(ncols))
#         for k in range(ncols):
#             V[j][k]=0
#     Xm=<double*>sig_malloc(sizeof(double)*Ql)
#     if Xm==NULL: raise MemoryError
#     Xpb = <double***> sig_malloc( sizeof(double** ) * nc )
#     if Xpb==NULL: raise MemoryError
#     Ypb = <double***> sig_malloc( sizeof(double** ) * nc )
#     if Ypb==NULL: raise MemoryError
#     for i in range(nc):
#         Xpb[i]=NULL; Ypb[i]=NULL
#         Xpb[i] = <double**>sig_malloc(sizeof(double*) * nc )
#         Ypb[i] = <double**>sig_malloc(sizeof(double*) * nc )
#         if Ypb[i]==NULL or Xpb[i]==NULL:
#             raise MemoryError
#         for j in range(nc):
#             Xpb[i][j]=NULL; Ypb[i][j]=NULL
#             Xpb[i][j] = <double*>sig_malloc(sizeof(double) * Ql )
#             Ypb[i][j] = <double*>sig_malloc(sizeof(double) * Ql )
#             if Ypb[i][j]==NULL or Xpb[i][j]==NULL:
#                 raise MemoryError
#             for n in range(Ql):
#                 Xpb[i][j][n]=<double>0
#                 Ypb[i][j][n]=<double>0
#     Cvec = <double complex***>sig_malloc(sizeof(double complex**) * nc )
#     if Cvec==NULL: raise MemoryError
#     for i in range(nc):
#         Cvec[i] = <double complex**>sig_malloc(sizeof(double complex*) * nc )
#         if Cvec[i]==NULL:
#             raise MemoryError
#         for j in range(nc):
#             Cvec[i][j] = <double complex*>sig_malloc(sizeof(double complex) * Ql )
#             if Cvec[i][j]==NULL:
#                 raise MemoryError
#             for n in range(Ql):
#                 Cvec[i][j][n]=0
#     cdef int t = 0
#     t = pullback_pts_cplx_dp(S,Qs,Qf,Y,Xm,Xpb,Ypb,Cvec)
#     if t==1:
#         raise ArithmeticError,"Need smaller Y than {0}".format(Y)
#     cdef double * alphas=NULL
#     alphas=<double*>sig_malloc(sizeof(double)*nc)
#     for j in range(nc):
#         alphas[j]=<double>S.alpha(j)[0]
#     cdef int cuspidal=1
#     cdef int** Mv, **Qv
#     Mv=<int**>sig_malloc(sizeof(int*)*nc)
#     Qv=<int**>sig_malloc(sizeof(int*)*nc)
#     for i in range(nc):
#         Mv[i]=<int*>sig_malloc(sizeof(int)*3)
#         Qv[i]=<int*>sig_malloc(sizeof(int)*3)        
#         Mv[i][0]=Ms; Mv[i][1]=Mf; Mv[i][2]=Mf-Ms+1
#         Qv[i][0]=Qs; Qv[i][1]=Qf; Qv[i][2]=Qf-Qs+1
#     cdef int* symmetric_cusps
#     symmetric_cusps=<int*> sig_malloc(sizeof(int)*nc)
#     cdef double complex* cusp_evs
#     cusp_evs=<double complex*>sig_malloc(sizeof(double complex)*nc)
#     for i in range(nc):
#         cusp_evs[i]=0
#         symmetric_cusps[i]=-1
#     if do_par==0:
#         if weight==0.0:
#             compute_V_cplx_dp(V,R,Y,Mv,Qv,nc,cuspidal,sym_type,verbose,alphas,Xm,Xpb,Ypb,Cvec,is_exceptional=is_exceptional)
#         else:
#             compute_V_cplx_wt_dp(V,R,Y,weight,Mv,Qv,nc,cuspidal,sym_type,verbose,alphas,Xm,Xpb,Ypb,Cvec)
#     else:
#         if weight==0.0:
#             compute_V_cplx_dp_sym_par(V,N,Xm,Xpb,Ypb,Cvec,
#                               cusp_evs,alphas,Mv,Qv,Qfak,
#                               symmetric_cusps,
#                               R,Y,nc,ncols,cuspidal,verbose,ncpus,0)
# #compute_V_cplx_dp_sym_par(V,R,Y,Mv,Qv,nc,cuspidal,sym_type,verbose,alphas,Xm,Xpb,Ypb,Cvec)
#         else:
#             compute_V_cplx_wt_dp(V,R,Y,weight,Mv,Qv,nc,cuspidal,sym_type,verbose,alphas,Xm,Xpb,Ypb,Cvec)
            
#     ## Try to make coefficients real if possible.
#     cdef int q
#     cdef double complex *sqch=NULL
#     # if norm_c==1:
#     #     chi = S.multiplier()._character
#     #     q   = chi.conductor()
#     #     sqch=<double complex *>sig_malloc(sizeof(double complex)*q)
#     #     for i from 0<=i<q:
#     #         sqch[i]=csqrt(<double complex>CC(chi(i)))
#     #     for n from 0<=n< nc*Ml:
#     #         for l in range(Ml):
#     #             if gcd(l-Ms,q)==1:
#     #                 i =int(  (l-Ms)%q )
#     #                 V[n][l]=V[n][l]*sqch[i]

#     if alphas<>NULL:
#         sig_free(alphas)
#     if Ypb<>NULL:
#         for i in range(nc):
#             if Ypb[i]<>NULL:
#                 for j in range(nc):
#                     if Ypb[i][j]<>NULL:
#                         sig_free(Ypb[i][j])
#                 sig_free(Ypb[i])
#         sig_free(Ypb)
#     if Xpb<>NULL:
#         for i in range(nc):
#             if Xpb[i]<>NULL:
#                 for j in range(nc):
#                     if Xpb[i][j]<>NULL:
#                         sig_free(Xpb[i][j])
#                 sig_free(Xpb[i])
#         sig_free(Xpb)
#     if Cvec<>NULL:
#         for i in range(nc):
#             if Cvec[i]<>NULL:
#                 for j in range(nc):
#                     if Cvec[i][j]<>NULL:
#                         sig_free(Cvec[i][j])
#                 sig_free(Cvec[i])
#         sig_free(Cvec)
#     if sqch<>NULL:
#         sig_free(sqch)


#     #cdef cnp.ndarray[CTYPE_t,ndim=2] VV
#     cdef Matrix_complex_dense VV
#     if gr==1:
#         #VV=mpmath.fp.matrix(int(nc*Ml),int(nc*Ml),force_type=mpmath.fp.mpc)
#         CF=MPComplexField(53)
#         MS=MatrixSpace(CF,nc*Ml,nc*Ml)
#         VV=Matrix_complex_dense(MS,0)
#         #VV=np.zeros([nc*Ml, nc*Ml], dtype=CTYPE)
#         for n in range(nc*Ml):
#             for l in range(nc*Ml):
#                 VV[n,l]=V[n][l]
#         return VV
#     C=<double complex**>sig_malloc(sizeof(double complex*)*(comp_dim))
#     if C==NULL:
#         raise MemoryError
#     for i in range(comp_dim):
#         C[i]=<double complex*>sig_malloc(sizeof(double complex)*(nc*Ml))
#         if C[i]==NULL:
#             raise MemoryError
#     setc_list = <int*>sig_malloc(sizeof(int)*num_set)
#     if setc_list==NULL:
#         raise MemoryError
#     vals_list = <double complex** >sig_malloc(sizeof(double complex*)*comp_dim)
#     if vals_list==NULL:
#         raise MemoryError
#     for j in range(comp_dim):
#         vals_list[j]=<double complex*>sig_malloc(sizeof(double complex)*num_set)
#     i=0
#     for r,n in SetCs:
#         if r*Ml+n-Ms<0:
#             continue
#         setc_list[i]=r*Ml+n-Ms
#         for j in range(comp_dim):
#             vals_list[j][i]=<double complex>CC(Vals[j][(r,n)])
#         if verbose>1:
#             print "setc_list[",i,"]=",setc_list[i]
#             for j in range(comp_dim):
#                 print "vals_list[",i,j,"]=",vals_list[j][i]
#         i=i+1
#     if num_rhs>0 and num_rhs<>comp_dim:
#         raise ValueError,"Need same number of right hand sides (or just one) as the number of set coefficients!"
#     normalize_matrix_cplx_dp(V,N,comp_dim,num_set,setc_list,vals_list,verbose)
#     if gr==2:
#         VV=np.zeros([nc*Ml, nc*Ml], dtype=CTYPE)
#         #VV=mpmath.fp.matrix(int(N),int(N+comp_dim),force_type=mpmath.fp.mpc)
#         for n in range(N):
#             for l in range(ncols): #N+comp_dim:
#                 VV[n,l]=V[n][l]
#         return VV
#     SMAT_cplx_dp(V,ncols-num_set,comp_dim,num_set,C,vals_list,setc_list)
#     if verbose>1:
#         for i in range(comp_dim):
#             for k in range(ncols-num_set):
#                 print "C[{0}][{1}]={2}".format(i,k,C[i][k])

#     cdef dict res={}
#     for k in range(comp_dim):
#         res[k]=dict()
#         for i in range(nc):
#             res[k][i]=dict()
#             for j in range(Ml):
#                 res[k][i][j+Ms]=C[k][i*Ml+j]
#     if setc_list<>NULL:
#         sig_free(setc_list)
#     if vals_list<>NULL:
#         for j in range(comp_dim):
#             if vals_list[j]<>NULL:
#                 sig_free(vals_list[j])
#         sig_free(vals_list)
#     if V<>NULL:
#         for j in range(N):
#             if V[j]<>NULL:
#                 sig_free(V[j])
#         sig_free(V)
#     if C<>NULL:
#         for i in range(comp_dim):
#             if C[i]<>NULL:
#                 sig_free(C[i])
#         sig_free(C)
#     #if comp_dim==1:
#     #    return res[0]
#     return res


cpdef get_coeff_fast_real_dp(S,double R,double Y,int M,int Q,dict Norm={},int gr=0):
    r"""
    An efficient method to get coefficients in the double complex case.
    """
    cdef double **V=NULL
    cdef double* Xm=NULL
    cdef double*** Xpb=NULL
    cdef double*** Ypb=NULL
    cdef double *** Cvec=NULL
    cdef double **C=NULL
    cdef int *setc_list=NULL
    cdef double **vals_list=NULL
    cdef double * alphas=NULL
    cdef int nc,Ql,Qs,Qf,Ml,Ms,Mf,j,k,l,N,sym_type,i,n,r
    cdef double Qfak
    cdef int verbose=S._verbose
    ## Check that this "real" method is applicable here.
    if not S.multiplier().is_real():
        raise ValueError,"This method should only be called for real characters/multipliers!"
    sym_type = S._sym_type
    if sym_type not in [0,1]:
        raise ValueError,"This method should only be called for symmetrized (even/odd) functions!"
    if Q<M:
        Q=M+20
    Ms=0;  Mf=M; Qs=1; Qf=Q
    Qfak=<double>(Q)/<double>(2)
    Ml=Mf-Ms+1
    Ql=Qf-Qs+1
    nc = S.group().ncusps()
    N = nc*Ml
    if Norm=={}:
        Norm = S.set_norm(1)
    cdef list SetCs=[]
    cdef dict Vals={}
    SetCs=Norm['SetCs'][0]
    Vals=Norm['Vals']
    cdef int comp_dim,num_set
    comp_dim=Norm['comp_dim']
    num_set=len(SetCs)
    V=<double**>sig_malloc(sizeof(double*)*N)
    if V==NULL: raise MemoryError
    cdef int ncols
    if num_set<comp_dim:
        ncols=N+comp_dim-num_set
    else:
        ncols=N
    if verbose>1:
        print "In get_coef_real_dp R=",R
    if verbose>2:
        print "N=",N
        print "Vals=",Vals
        print "SetCs=",SetCs
        print "comp_dim=",comp_dim
        print "num_set=",num_set
        print "ncols=",ncols
    for j in range(N):
        V[j]=<double*>sig_malloc(sizeof(double)*(ncols))
        for k in range(ncols):
            V[j][k]=<double>0
    Xm=<double*>sig_malloc(sizeof(double)*Ql)
    if Xm==NULL: raise MemoryError
    #printf("Xm0: %p ", Xm)
    #my_alloc_double_2(&V,ncols,Ql)
    #my_alloc_double_1(&Xm,Ql)
    #printf("%p ", Xm)
    if Xm==NULL:
        print "Xm1=NULL"
    Xpb = <double***> sig_malloc( sizeof(double** ) * nc )
    if Xpb==NULL: raise MemoryError
    Ypb = <double***> sig_malloc( sizeof(double** ) * nc )
    if Ypb==NULL: raise MemoryError
    for i in range(nc):
        Xpb[i] = <double**>sig_malloc(sizeof(double*) * nc )
        Ypb[i] = <double**>sig_malloc(sizeof(double*) * nc )
        if Ypb[i]==NULL or Xpb[i]==NULL:
            raise MemoryError
        for j in range(nc):
            Xpb[i][j] = <double*>sig_malloc(sizeof(double) * Ql )
            Ypb[i][j] = <double*>sig_malloc(sizeof(double) * Ql )
            if Ypb[i][j]==NULL or Xpb[i][j]==NULL:
                raise MemoryError
            for n in range(Ql):
                Xpb[i][j][n]=<double>0
                Ypb[i][j][n]=<double>0
    Cvec = <double ***>sig_malloc(sizeof(double **) * nc )
    if Cvec==NULL: raise MemoryError
    for i in range(nc):
        Cvec[i] = <double **>sig_malloc(sizeof(double *) * nc )
        if Cvec[i]==NULL:
            raise MemoryError
        for j in range(nc):
            Cvec[i][j] = <double *>sig_malloc(sizeof(double) * Ql )
            if Cvec[i][j]==NULL:
                raise MemoryError
            for n in range(Ql):
                Cvec[i][j][n]=<double>0
    cdef int t=0
    t = pullback_pts_real_dp(S,Qs,Qf,Y,Xm,Xpb,Ypb,Cvec)
    if t==1:
        raise ArithmeticError,"Need smaller Y than {0}".format(Y)    
    alphas=<double*>sig_malloc(sizeof(double)*nc)
    for j in range(nc):
        alphas[j]=float(S.alpha(j)[0])
    compute_V_real_dp(V,R,Y,Ms,Mf,Qs,Qf,nc,alphas,Xm,Xpb,Ypb,Cvec,1,sym_type,verbose)
    cdef Matrix_complex_dense VV
    if gr==1:
        CF=MPComplexField(53)
        MS=MatrixSpace(CF,nc*Ml,nc*Ml)
        VV=Matrix_complex_dense(MS,0)
        for n in range(nc*Ml):
            for l in range(nc*Ml):
                VV[n,l]=V[n][l]
        return VV

    C=<double**>sig_malloc(sizeof(double*)*(comp_dim))
    if C==NULL:
        raise MemoryError
    for i in range(comp_dim):
        C[i]=<double*>sig_malloc(sizeof(double)*(nc*Ml))
        if C[i]==NULL:
            raise MemoryError
    setc_list = <int*>sig_malloc(sizeof(int)*num_set)
    vals_list = <double** >sig_malloc(sizeof(double*)*comp_dim)
    for j in range(comp_dim):
        vals_list[j]=<double*>sig_malloc(sizeof(double)*num_set)
    i=0
    for r,n in SetCs:
        if r*Ml+n-Ms<0:
            continue
        setc_list[i]=r*Ml+n-Ms
        for j in range(comp_dim):
            vals_list[j][i]=<double>CC(Vals[j][(r,n)])
        i=i+1
    cdef double *RHS=NULL
    cdef num_rhs=0
    if num_rhs>0 and num_rhs<>comp_dim:
        raise ValueError,"Need same number of right hand sides (or just one) as the number of set coefficients!"
    normalize_matrix_real_dp(V,N,comp_dim,num_set,setc_list,vals_list,verbose)
    if gr==2:
        CF=MPComplexField(53)
        MS=MatrixSpace(CF,N,ncols)
        VV=Matrix_complex_dense(MS,0)
        for n in range(N):
            for l in range(ncols):
                VV[n,l]=V[n][l]
        if V<>NULL:
            sig_free(V)
        if setc_list <>NULL:
            sig_free(setc_list)
        if vals_list<>NULL:
            for i in range(comp_dim):
                if vals_list[i]<>NULL:
                    sig_free(vals_list[i])
            sig_free(vals_list)
        return VV
    SMAT_real_dp(V,ncols-num_set,comp_dim,num_set,C,vals_list,setc_list)
    if setc_list <>NULL:
        sig_free(setc_list)
    if vals_list<>NULL:
        for i in range(comp_dim):
            if vals_list[i]<>NULL:
                sig_free(vals_list[i])
        sig_free(vals_list)
    my_dealloc_double_2(&V,ncols)
    my_dealloc_double_1(&Xm)
    #if V<>NULL:
    #    sig_free(V)

    cdef dict res={}
    for k in range(comp_dim):
        res[k]=dict()
        for i in range(nc):
            res[k][i]=dict()
            for j in range(Ml):
                #if j==0:
                #    print "C[",k,i,j+Ms,"]=",C[k][i*Ml+j]
                res[k][i][j+Ms]=C[k][i*Ml+j]
    if C<>NULL:
        for i in range(comp_dim):
            if C[i]<>NULL:
                sig_free(C[i])
        sig_free(C)
    if alphas<>NULL:
        sig_free(alphas)
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
    return res

## Routines to make allocation and deallocation easier


cdef int my_alloc_double_1(double **A,int d1):
    r"""
    Allocate a 1-dim array of doubles.
    """
    A[0]=<double*>sig_malloc(sizeof(double)*d1)
    if A==NULL or A[0]==NULL:
        raise MemoryError
    printf("Allocated: %p",A[0])
    return 1

cdef void my_dealloc_double_1(double **A):
    r"""
    Deallocate a 1-dim array of doubles.
    """
    if A[0]<>NULL:
        sig_free(A[0])


cdef my_alloc_double_2(double ***A,int d1,int d2):
    r"""
    Allocate a 2-dim array of doubles.
    """
    cdef int i
    A[0]=<double**>sig_malloc(sizeof(double*)*d1)
    if A==NULL or A[0]==NULL:
        raise MemoryError
    for i in range(d1):
        A[0][i]=<double*>sig_malloc(sizeof(double)*d2)

cdef my_dealloc_double_2(double ***A,int d1):
    r"""
    Allocate a 2-dim array of doubles.
    """
    cdef int i
    if A[0]<>NULL:
        for i in range(d1):
            if A[0][i]<>NULL:
                sig_free(A[0][i])
        sig_free(A[0])

cdef my_alloc_double_3(double ***A,int d1,int d2,int d3):
    r"""
    Allocate a 2-dim array of doubles.
    """
    cdef int i,j
    A=<double***>sig_malloc(sizeof(double**)*d1)
    if A==NULL:
        raise MemoryError
    for i in range(d1):
        A[i]=<double**>sig_malloc(sizeof(double*)*d2)
        if A[i]==NULL:
            raise MemoryError
        for j in range(d2):
            A[i][j]=<double*>sig_malloc(sizeof(double)*d3)
            if A[i][j]==NULL:
                raise MemoryError

cdef my_dealloc_double_3(double ***A,int d1,int d2):
    r"""
    Allocate a 2-dim array of doubles.
    """
    cdef int i,j
    if A<>NULL:
        for i in range(d1):
            if A[i]<>NULL:
                for j in range(d2):
                    if A[i][j]<>NULL:
                        sig_free(A[i][j])
                sig_free(A[i])
        sig_free(A)
### Start using MPC directly instead of mpmath....
### Cythonize....

cpdef setup_matrix_for_Maass_waveforms_np_cplx2(S,RealNumber R,RealNumber Yin,int M,int Q,int cuspidal=1,int tilde=1,do_hp=1,filter=None):
    r"""
    Filter (if set) is a list of the rows that I will compute.
    Ex.: filter=[0,1,1,0,,,,]
    """
    import mpmath
    from sage.all import bessel_K,ComplexField
    cdef int l,j,icusp,jcusp,n,ni,li,Ml,Ms,Mf,Qs,Qf,Ql,s,nc
    cdef mpc_t ckbes,ctmpV,iargm,iargpb,twopii
    cdef mpfr_t Qfak,sqrtY,Y2pi,nrY2pi,argm,argpb,twopi,two,kbes
    cdef mpfr_t besarg,lr,besfak,xtmp_t
    cdef Matrix_complex_dense VV
    cdef RealNumber Y,tmpx1,tmpx2,ypb,pi
    cdef MPComplexNumber tmpz1,tmpz2
    cdef ComplexNumber IR
    cdef mpfr_t *nvec
    cdef mpc_t *** ef1
    cdef mpc_t ** ef2
    #cdef dict ef1
    cdef double dbes,dbesarg,Rd
    cdef int prec
    RF = R.parent()
    prec=RF.prec()
    pi=RF.pi()
    CF = MPComplexField(prec)
    mpc_init2(ckbes,prec)
    mpc_init2(ctmpV,prec)
    mpc_init2(iargm,prec)
    mpc_init2(iargpb,prec)
    mpc_init2(twopii,prec)

    mpfr_init2(besarg,prec)
    mpfr_init2(besfak,prec)
    mpfr_init2(lr,prec)
    mpfr_init2(Qfak,prec)
    mpfr_init2(sqrtY,prec)
    mpfr_init2(Y2pi,prec)
    mpfr_init2(nrY2pi,prec)
    mpfr_init2(argm,prec)
    mpfr_init2(argpb,prec)
    mpfr_init2(two,prec)
    mpfr_init2(twopi,prec)
    mpfr_init2(kbes,prec)
    mpfr_init2(xtmp_t,prec)
    cdef int verbose = S._verbose
    #mp_ctx=mpmath.fp
    if do_hp==1 and prec>53:
        hprec=1
    else:
        hprec=0
        Rd=<double>R
    tmpx1=RF(0); tmpx2=RF(0);ypb=RF(0)
    tmpz1=CF(0); tmpz2=CF(0)
    CCF=ComplexField(prec)
    IR=CCF(0,R.real())
    if verbose>0:
        print "hprec=",hprec
        print "CF=",CF
        print "CCF=",CCF
        print "IR=",IR
        print "mpfr_prec=",mpfr_get_default_prec()
    if hprec:
        mpmath.mp.prec=prec
        mpfr_div_ui(besarg,pi.value,2,rnd_re)
        mpfr_mul(besarg,besarg,R.value,rnd_re)
        mpfr_exp(besfak,besarg,rnd_re)
        mpIR=mpmath.mp.mpc(0,R)
        if verbose>0:
            print "hpfak=",print_mpfr(besfak)
    else:
        mpfr_set_ui(besfak,1,rnd_re)
        #besfak=RF(exp(pi/2*R))
        if verbose>0:
            print "lpfak=",print_mpfr(besfak)
    Y = RF(Yin)
    mpfr_sqrt(sqrtY,Y.value,rnd_re)
    pi=RF.pi()
    mpfr_mul_ui(twopi,pi.value,2,rnd_re)
    mpfr_mul(Y2pi,Y.value,twopi,rnd_re)
    G = S.group()
    if verbose>0:
        print "S=",S
        print "G=",G
        print "Yin=",Yin
        print "Y=",Y
    character = S.character
    nc=int(G.ncusps())
    if(Q<M):
        Q=M+20
    mpfr_set_default_prec(prec) #else:
    #    Q=Q+10
    cdef int sym_type
    sym_type = <int>S._sym_type
    if sym_type in [0,1]:
        Ms=0;  Mf=M; Ml=Mf-Ms+1
        Qs=1; Qf=Q; Ql=Qf-Qs+1
        #mpfr_set_ui(Qfak,Q,rnd_re)
        #mpfr_div_ui(Qfak,Qfak,2,rnd_re)
        if verbose>0:
            print "using symmetry:",sym_type
    else:
        Ms=-M;  Mf=M; Ml=Mf-Ms+1
        Qs=1-Q; Qf=Q; Ql=Qf-Qs+1
    # Recall that we always use the same factor and change the functions sin/cos instead. 
    mpfr_set_ui(Qfak,Q,rnd_re)
    mpfr_mul_ui(Qfak,Qfak,2,rnd_re)
    if verbose>0:
        print "Qfak=",print_mpfr(Qfak)
        print "Qf=",Qf,"Qs=",Qs,"Ql=",Ql
        print "Mf=",Mf,"Ms=",Ms,"Ml=",Ml
    cdef mpfr_t* Xm
    cdef mpfr_t ***Xpb=NULL
    cdef mpfr_t ***Ypb=NULL
    cdef mpc_t ***Cvec=NULL
    Xm = <mpfr_t*> sig_malloc( sizeof(mpfr_t) * Ql )
    for n in range(Ql):
        mpfr_init2(Xm[n],prec)
    Xpb = <mpfr_t***> sig_malloc( sizeof(mpfr_t** ) * nc )
    if Xpb==NULL: raise MemoryError
    Ypb = <mpfr_t***> sig_malloc( sizeof(mpfr_t** ) * nc )
    if Ypb==NULL: raise MemoryError
    for i in range(nc):
        Xpb[i] = <mpfr_t**>sig_malloc(sizeof(mpfr_t*) * nc )
        Ypb[i] = <mpfr_t**>sig_malloc(sizeof(mpfr_t*) * nc )
        if Ypb[i]==NULL or Xpb[i]==NULL:
            raise MemoryError
        for j in range(nc):
            Xpb[i][j] = <mpfr_t*>sig_malloc(sizeof(mpfr_t) * Ql )
            Ypb[i][j] = <mpfr_t*>sig_malloc(sizeof(mpfr_t) * Ql )
            if Ypb[i][j]==NULL or Xpb[i][j]==NULL:
                raise MemoryError
            for n in range(Ql):
                mpfr_init_set_ui(Xpb[i][j][n],0,rnd_re)
                mpfr_init_set_ui(Ypb[i][j][n],0,rnd_re)
    Cvec = <mpc_t***>sig_malloc(sizeof(mpc_t**) * nc )
    if Cvec==NULL:
        raise MemoryError
    for i in range(nc):
        Cvec[i] = <mpc_t**>sig_malloc(sizeof(mpc_t*) * nc )
        if Cvec[i]==NULL:
            raise MemoryError
        for j in range(nc):
            Cvec[i][j] = <mpc_t*>sig_malloc(sizeof(mpc_t) * Ql )
            if Cvec[i][j]==NULL:
                raise MemoryError
            for n in range(Ql):
                mpc_init2(Cvec[i][j][n],prec)
    cdef RealNumber weight
    weight=RF(0)
    cdef int t=0
    t = pullback_pts_mpc_new_c(S,Qs,Qf,Y.value,Xm,Xpb,Ypb,Cvec)
    if t==1:
        raise ArithmeticError,"Need smaller Y than {0}".format(Y)
    #pb=pullback_pts_mpc_new(S,Qs,Qf,Y)
    #Xm=pb['xm']; Xpb=pb['xpb']; Ypb=pb['ypb']; Cv=pb['cvec']
    s=nc*Ml
    ef2 = <mpc_t**>sig_malloc(sizeof(mpc_t*) * Ml )
    for n in range(Ml):
        ef2[n ]= <mpc_t*>sig_malloc(sizeof(mpc_t) * Ql )
    # somehow sizeof(mpc_t ***) did not work...
    ef1 = <mpc_t***>sig_malloc(sizeof(mpc_t **) * Ml )
    for n in range(Ml):
        ef1[n]= <mpc_t**>sig_malloc(sizeof(mpc_t*) * nc*nc )
        for icusp in range(nc*nc):
            ef1[n][icusp]= <mpc_t*>sig_malloc(sizeof(mpc_t) * Ql  )
    nvec = <mpfr_t*>sig_malloc(sizeof(mpfr_t) * Ml)
    for n in range(Ml):
        mpfr_init2(nvec[n],prec)
        mpfr_set_si(nvec[n],n+Ms,rnd_re)
    MS=MatrixSpace(CF,s,s)
    VV=Matrix_complex_dense(MS,CF.zero(),True,True)
    if verbose>0:
        print "sym_type=",sym_type
    cdef int do_filter=0
    cdef int* filter_rows
    if filter<>None:
        do_filter=1
        filter_rows=<int*>sig_malloc(sizeof(int)*Ml)
        for n in range(Ml):
            filter_rows[n]=int(filter[n])
        if verbose>0:
            print "filter:",filter
    for n in range(Ml):
        for j in range(Ql):
            #tmpx1=Xm[j-Qs]
            mpfr_mul(argm,nvec[n],Xm[j],rnd_re) #tmpx1.value,rnd_re)
            # argm -> i*argm
            if sym_type==1:
                mpfr_neg(argm,argm,rnd_re)
                mpfr_sin(argm,argm,rnd_re)
                mpfr_add(argm,argm,argm,rnd_re) ## *2
                mpc_set_fr(iargm,argm,rnd)
            elif sym_type==0:
                mpfr_cos(argm,argm,rnd_re)
                mpfr_add(argm,argm,argm,rnd_re) ## *2
                mpc_set_fr(iargm,argm,rnd)
            else:
                mpfr_neg(argm,argm,rnd_re)
                mpfr_sin_cos(mpc_imagref(iargm),mpc_realref(iargm),argm,rnd_re) ## exp=sin+I*cos
            mpc_init2(ef2[n][j],prec)
            mpc_set(ef2[n][j],iargm,rnd)
            #mpc_set(tmpz1.value,ef2[n][j],rnd)
    for jcusp in range(nc):
        for icusp in range(nc):
            for j in range(Ql):
                # print "i,j=",icusp,jcusp
                if mpfr_zero_p(Ypb[icusp][jcusp][j])<>0:
                    continue
                #mpc_set(tmpz1.value,Cvec[icusp][jcusp][j],rnd)
                # mpfr_set(xtmp_t,Xpb[icusp][jcusp][j],rnd_re)
                for n in range(Ml):
                    mpfr_mul(argpb,nvec[n],Xpb[icusp][jcusp][j],rnd_re)
                    if sym_type==1:
                        mpfr_sin(argpb,argpb,rnd_re)
                        mpc_mul_fr(iargpb,tmpz1.value,argpb,rnd)
                    elif sym_type==0:
                        mpfr_cos(argpb,argpb,rnd_re)
                        mpc_mul_fr(iargpb,tmpz1.value,argpb,rnd)
                    else:
                        #mpfr_cos(mpc_realref(iargpb),argpb,rnd_re)
                        #mpfr_sin(mpc_imagref(iargpb),argpb,rnd_re)
                        mpfr_sin_cos(mpc_imagref(iargpb),mpc_realref(iargpb),argpb,rnd_re)
                    mpc_mul(iargpb,iargpb,Cvec[icusp][jcusp][j],rnd)
                    mpc_init2(ef1[n][icusp*nc+jcusp][j],prec)
                    mpc_set(ef1[n][icusp*nc+jcusp][j],iargpb,rnd)


    cdef double besprec
    besprec=1.0E-14

    for i in range(Ml):
        for l in range(Ml):
            if l+Ms==0 and cuspidal==1:
                continue
            mpfr_mul(lr,nvec[l],twopi,rnd_re)
            mpfr_abs(lr,lr,rnd_re)
            for j in range(Ql):
                for jcusp in range(nc):
                    lj=Ml*jcusp+l
                    for icusp in range(nc):
                        if mpfr_zero_p(Ypb[icusp][jcusp][j])<>0:
                            continue
                        mpfr_mul(besarg,lr,Ypb[icusp][jcusp][j],rnd_re)
                        if hprec==1:
                            mpfr_set(tmpx1.value,besarg,rnd_re)
                            tmpx2=RF(mpmath.mp.besselk(mpIR,tmpx1).real)
                            mpfr_sqrt(kbes,Ypb[icusp][jcusp][j],rnd_re)
                            mpfr_mul(kbes,kbes,besfak,rnd_re)
                            mpfr_mul(kbes,kbes,tmpx2.value,rnd_re)
                        else:
                            dbesarg = mpfr_get_d(besarg,rnd_re)
                            besselk_dp_c(&dbes,Rd,dbesarg,besprec,1)
                            mpfr_sqrt(kbes,Ypb[icusp][jcusp][j],rnd_re)
                            mpfr_mul_d(kbes,kbes,dbes,rnd_re)
                        mpc_mul_fr(ckbes,ef1[l][icusp*nc+jcusp][j],kbes,rnd)
                        for n in range(Ml):
                            if n+Ms==0 and cuspidal==1:
                                continue
                            #if filter<>None and (n > len(filter) or filter[n]==0):
                            #s    continue
                            ni=Ml*icusp+n
                            mpc_mul(ctmpV,ckbes,ef2[n][j],rnd)
                        mpc_add(VV._matrix[ni][lj],VV._matrix[ni][lj],ctmpV,rnd)
    for n in range(s):
        for l in range(s):
            mpc_div_fr(VV._matrix[n][l],VV._matrix[n][l],Qfak,rnd)
    if tilde==1:
        for n in range(Ms,Mf+1):
            if(n==0 and cuspidal==1):
                continue
            for icusp in range(nc):
                ni=Ml*icusp+n-Ms
                mpfr_mul(nrY2pi,nvec[n-Ms],Y2pi,rnd_re)
                mpfr_abs(nrY2pi,nrY2pi,rnd_re)
                if hprec==1:
                    mpfr_set(tmpx1.value,nrY2pi,rnd_re)
                    tmpx2=RF(mpmath.mp.besselk(mpIR,tmpx1).real)
                    #tmpx2=bessel_K(IR,tmpx1,prec=prec)
                    mpfr_mul(kbes,sqrtY,tmpx2.value,rnd_re)
                    mpfr_mul(kbes,kbes,besfak,rnd_re)
                else:
                    dbesarg = mpfr_get_d(nrY2pi,rnd_re)
                    #dbes = besselk_dp(R,dbesarg,pref=1)
                    besselk_dp_c(&dbes,Rd,dbesarg,besprec,1)
                    mpfr_mul_d(kbes,sqrtY,dbes,rnd_re)
                mpc_sub_fr(VV._matrix[ni][ni],VV._matrix[ni][ni],kbes,rnd)
    W=dict()
    #VV=mp_ctx.matrix(int(s),int(s),force_type=mp_ctx.mpc)
    #for n in range(s):
    #    for l in range(s):
    #        VV[n,l]=V[n,l]
    W['space']=S
    W['V']=VV
    W['Ms']=Ms
    W['Mf']=Mf
    W['Ml']=Ml
    W['nc']=nc
    W['Y']=Y
    W['R']=R
    W['Q']=Q
    mpc_clear(ckbes); mpc_clear(ctmpV)
    mpc_clear(iargm); mpc_clear(twopii)
    mpfr_clear(besarg); mpfr_clear(lr)
    mpfr_clear(Qfak); mpfr_clear(sqrtY)
    mpfr_clear(Y2pi); mpfr_clear(nrY2pi)
    mpfr_clear(argm); mpfr_clear(argpb)
    mpfr_clear(two); mpfr_clear(twopi)
    mpfr_clear(kbes)
    mpfr_clear(xtmp_t)
    for n in range(Ml):
        mpfr_clear(nvec[n])
        for icusp in range(nc):
            for jcusp in range(nc):
                for i in range(Ql):
                    # print n,icusp,i
                    if mpfr_zero_p(Ypb[icusp][jcusp][i]):
                        continue
                    mpc_clear(ef1[n][icusp*nc+jcusp][i])
                if ef1[n][icusp*nc+jcusp]<>NULL:
                    sig_free(ef1[n][icusp*nc+jcusp])
        for i in range(Ql):
            mpc_clear(ef2[n][i])
        if ef1[n]<>NULL:
            sig_free(ef1[n])
        if ef2[n]<>NULL:
            sig_free(ef2[n])
    if ef1<>NULL:
        sig_free(ef1)
    if ef2<>NULL:
        sig_free(ef2)
    if nvec<>NULL:
        sig_free(nvec)


    if Ypb<>NULL:
        for i in range(nc):
            if Ypb[i]<>NULL:
                for j in range(nc):
                    if Ypb[i][j]<>NULL:
                        for n in range(Ql):
                            if Ypb[i][j][n]<>NULL:
                                mpfr_clear(Ypb[i][j][n])
                        sig_free(Ypb[i][j])
                sig_free(Ypb[i])
        sig_free(Ypb)
    if Xpb<>NULL:
        for i in range(nc):
            if Xpb[i]<>NULL:
                for j in range(nc):
                    if Xpb[i][j]<>NULL:
                        for n in range(Ql):
                            if Xpb[i][j][n]<>NULL:
                                mpfr_clear(Xpb[i][j][n])
                        sig_free(Xpb[i][j])
                sig_free(Xpb[i])
        sig_free(Xpb)
    if Cvec<>NULL:
        for i in range(nc):
            if Cvec[i]<>NULL:
                for j in range(nc):
                    if Cvec[i][j]<>NULL:
                        for n in range(Ql):
                            mpc_clear(Cvec[i][j][n])
                        sig_free(Cvec[i][j])
                sig_free(Cvec[i])
        sig_free(Cvec)

    return W



cpdef setup_matrix_for_Maass_waveforms_np_cplx2_without_sym(S,RealNumber R,RealNumber Yin,int M,int Q,int cuspidal=1,int prec=53,int tilde=1):

    import mpmath
    cdef int l,j,icusp,jcusp,n,ni,li,Ml,Ms,Mf,Qs,Qf,Ql,s,nc
    cdef mpc_t ckbes,ctmpV,iargm,iargpb,twopii
    cdef mpfr_t Qfak,sqrtY,Y2pi,nrY2pi,argm,argpb,twopi,two,kbes
    cdef mpfr_t besarg,lr
    cdef Matrix_complex_dense VV
    cdef RealNumber Y,tmpx1,tmpx2,pi,tmp
    cdef MPComplexNumber tmpz1,tmpz2
    cdef mpfr_t *nvec
    cdef mpc_t *** ef1
    cdef mpc_t ** ef2
    #cdef dict ef1
    cdef double dbes,dbesarg
    mpc_init2(ckbes,prec)
    mpc_init2(ctmpV,prec)
    mpc_init2(iargm,prec)
    mpc_init2(iargpb,prec)
    mpc_init2(twopii,prec)

    mpfr_init2(besarg,prec)
    mpfr_init2(lr,prec)
    mpfr_init2(Qfak,prec)
    mpfr_init2(sqrtY,prec)
    mpfr_init2(Y2pi,prec)
    mpfr_init2(nrY2pi,prec)
    mpfr_init2(argm,prec)
    mpfr_init2(argpb,prec)
    mpfr_init2(two,prec)
    mpfr_init2(twopi,prec)
    mpfr_init2(kbes,prec)
    #mp_ctx=mpmath.fp
    RF = R.parent()
    CF = MPComplexField(RF.prec())
    tmpx1=RF(0); tmpx2=RF(0);
    tmp=RF(0)
    tmpz1=CF(0); tmpz2=CF(0)
    Y = RF(Yin)
    mpfr_sqrt(sqrtY,Y.value,rnd_re)
    pi=RF.pi()
    mpfr_mul_ui(twopi,pi.value,2,rnd_re)
    mpfr_mul(Y2pi,Y.value,twopi,rnd_re)
    G = S.group()
    if S._verbose>0:
        print "S=",S
        print "G=",G
        print "Yin=",Yin
        print "Y=",Y
    character = S.character
    nc=int(G.ncusps())
    if(Q<M):
        Q=M+20
    #else:
    #    Q=Q+10
    Ms=-M;  Mf=M; Ml=Mf-Ms+1
    Qs=1-Q; Qf=Q; Ql=Qf-Qs+1
    mpfr_set_ui(Qfak,Q,rnd_re)
    mpfr_mul_ui(Qfak,Qfak,2,rnd_re)
    # Qfak = 2*Q
    # Pullpack points

    pb=pullback_pts_mpc(S,Qs,Qf,Y,weight=S.weight())
    Xm=pb['xm']; Xpb=pb['xpb']; Ypb=pb['ypb']; Cv=pb['cvec']
    s=nc*Ml
    ef2 = <mpc_t**>sig_malloc(sizeof(mpc_t*) * Ml )
    for n in range(Ml):
        ef2[n ]= <mpc_t*>sig_malloc(sizeof(mpc_t) * Ql )
    # somehow sizeof(mpc_t ***) did not work...
    ef1 = <mpc_t***>sig_malloc(sizeof(mpc_t **) * Ml )
    #print "Ml,nc**2,Ql=",Ml,nc*nc,Ql
    for n in range(Ml):
        ef1[n]= <mpc_t**>sig_malloc(sizeof(mpc_t*) * nc*nc )
        for icusp in range(nc*nc):
            ef1[n][icusp]= <mpc_t*>sig_malloc(sizeof(mpc_t) * Ql  )
    nvec = <mpfr_t*>sig_malloc(sizeof(mpfr_t) * Ml)
    for n in range(Ml):
        mpfr_init2(nvec[n],prec)
        mpfr_set_si(nvec[n],n+Ms,rnd_re)
    MS=MatrixSpace(CF,s,s)
    VV=Matrix_complex_dense(MS,CF.zero(),True,True)
    for n in range(Ms,Mf+1):
        #nr=nvec[n-Ms]
        #print "nr(",n,")=",nr
        for j in range(Qs,Qf+1):
            tmpx1=RF(Xm[j])
            mpfr_mul(argm,nvec[n-Ms],tmpx1.value,rnd_re)
            mpfr_neg(argm,argm,rnd_re)
            #mpc_set_ui_fr(iargm,0,argm,rnd)
            mpc_set_fr(iargm,argm,rnd) # iargm = argm
            # argm -> i*argm
            mpfr_swap (mpc_realref (iargm), mpc_imagref (iargm))
            mpc_exp(iargm,iargm,rnd)
            mpc_init2(ef2[n-Ms][j-Qs],prec)
            mpc_set(ef2[n-Ms][j-Qs],iargm,rnd)
            #ef2[n,j]=mp_ctx.exp(mp_ctx.mpc(0,-argm))
            mpc_set(tmpz1.value,ef2[n-Ms][j-Qs],rnd)
            #print "ef2[",n,j,"=",tmpz1
            for jcusp in range(nc):
                for icusp in range(nc):
                    if(not Xpb.has_key((icusp,jcusp,j))):
                        continue
                    #argpb=(n)*Xpb[icusp,jcusp,j]
                    tmpx1=RF(Xpb[icusp,jcusp,j])
                    mpfr_mul(argpb,tmpx1.value,nvec[n-Ms],rnd_re)
                    mpc_set_fr(iargpb,argpb,rnd)
                    mpfr_swap (mpc_realref (iargpb), mpc_imagref (iargpb))
                    #mpc_swap(iargpb,rnd)
                    #argpb=nr*Xpb[icusp,jcusp,j]
                    #iargpb=i*argpb
                    #if(n==1):
                    #    print "Xpb[",icusp,jcusp,j,"]=",Xpb[icusp,jcusp,j]
                    #    print "argpb =",argpb,abs(argpb-Xpb[icusp,jcusp,j])
                    #    print "iargpb=",iargpb,abs(abs(iargpb)-abs(argpb))
                    mpc_exp(iargpb,iargpb,rnd)
                    #ef1[j,icusp,jcusp,n]=iargpb
                    #print n,icusp*nc+jcusp,j
                    #print "init:",n-Ms,icusp*nc+jcusp,j-Qs
                    try:
                        tmpz1=Cv[icusp,jcusp,j]
                    except:
                        tmpz1=CF(Cv[icusp,jcusp,j].real,Cv[icusp,jcusp,j].imag)
                    mpc_mul(iargpb,iargpb,tmpz1.value,rnd)
                    mpc_init2(ef1[n-Ms][icusp*nc+jcusp][j-Qs],prec)
                    mpc_set(ef1[n-Ms][icusp*nc+jcusp][j-Qs],iargpb,rnd)

    for l in range(Ms,Mf+1):
        if(l==0 and cuspidal==1):
            continue
        #lr=nvec[l-Ms]*twopi
        mpfr_mul(lr,nvec[l-Ms],twopi,rnd_re)
        for j in range(Qs,Qf+1):
            for jcusp in range(nc):
                lj=Ml*jcusp+l-Ms
                for icusp in range(nc):
                    if(not Ypb.has_key((icusp,jcusp,j))):
                        continue
                    #ypb=Ypb[icusp,jcusp,j]
                    tmpx1=RF(Ypb[icusp,jcusp,j])
                    if(tmpx1==0):
                        continue
                    mpfr_abs(besarg,lr,rnd_re)
                    mpfr_mul(besarg,besarg,tmpx1.value,rnd_re)
                    #besarg=abs(lr)*ypb
                    dbesarg = mpfr_get_d(besarg,rnd_re)
                    ## if lj==0 and icusp==0:
                    ##     mpfr_set(tmpx1.value,lr,rnd_re)
                    ##     print "lr=",tmpx1
                    ##     mpfr_set(tmpx1.value,besarg,rnd_re)
                    ##     print "besarg=",tmpx1
                    ##     print "dbesarg=",dbesarg
                    dbes = my_kbes(R,dbesarg,pref=1)
                    mpfr_sqrt(kbes,tmpx1.value,rnd_re)
                    mpfr_mul_d(kbes,kbes,dbes,rnd_re)
                    #kbes=mp_ctx.sqrt(ypb)*
                    #ckbes=kbes*ef1[j,icusp,jcusp,l]
                    #mpc_mul_fr(ckbes,ef1[j,icusp,jcusp,l],kbes,rnd)
                    mpc_mul_fr(ckbes,ef1[l-Ms][icusp*nc+jcusp][j-Qs],kbes,rnd)
                    for n in range(Ms,Mf+1):
                        if(n==0 and cuspidal==1):
                            continue
                        ni=Ml*icusp+n-Ms
                        mpc_mul(ctmpV,ckbes,ef2[n-Ms][j-Qs],rnd)
                        mpc_add(VV._matrix[ni][lj],VV._matrix[ni][lj],ctmpV,rnd)
                        #    print "V[16,16]=",V[16,16]
    for n in range(s):
        for l in range(s):
            mpc_div_fr(VV._matrix[n][l],VV._matrix[n][l],Qfak,rnd)
    #    print "V[16,16]=",V[16,16]
    #    print "sqrtY=",sqrtY
    if tilde==1:
        for n in range(Ms,Mf+1):
            if(n==0 and cuspidal==1):
                continue
            for icusp in range(nc):
                ni=Ml*icusp+n-Ms
                # nrY2pi=nr*Y2pi
                mpfr_mul(nrY2pi,nvec[n-Ms],Y2pi,rnd_re)
                mpfr_abs(nrY2pi,nrY2pi,rnd_re)
                dbesarg = mpfr_get_d(nrY2pi,rnd_re)
                dbes = my_kbes(R,dbesarg,pref=1)
                mpfr_mul_d(kbes,sqrtY,dbes,rnd_re)
                mpc_sub_fr(VV._matrix[ni][ni],VV._matrix[ni][ni],kbes,rnd)
            #V[ni,ni]=V[ni,ni]-ckbes
    W=dict()
    #print "V[16,16]=",V[16,16]
    #VV=mp_ctx.matrix(int(s),int(s),force_type=mp_ctx.mpc)
    #for n from 0<=n< s:
    #    for l from 0 <= l < s:
    #        VV[n,l]=V[n,l]
    W['space']=S
    W['V']=VV
    W['Ms']=Ms
    W['Mf']=Mf
    W['Ml']=Ml
    W['nc']=nc
    W['Y']=Y
    W['Q']=Q
    mpc_clear(ckbes); mpc_clear(ctmpV)
    mpc_clear(iargm); mpc_clear(twopii)
    mpfr_clear(besarg); mpfr_clear(lr)
    mpfr_clear(Qfak); mpfr_clear(sqrtY)
    mpfr_clear(Y2pi); mpfr_clear(nrY2pi)
    mpfr_clear(argm); mpfr_clear(argpb)
    mpfr_clear(two); mpfr_clear(twopi)
    mpfr_clear(kbes)

    for n in range(Ml):
        mpfr_clear(nvec[n])
        for icusp in range(nc):
            for jcusp in range(nc):
                for i in range(Ql):
                    # print n,icusp,i
                    if(not Xpb.has_key((icusp,jcusp,j))):
                        continue
                    mpc_clear(ef1[n][icusp*nc+jcusp][i])
                sig_free(ef1[n][icusp*nc+jcusp])
        for i in range(Ql):
            mpc_clear(ef2[n][i])
        sig_free(ef1[n])
        sig_free(ef2[n])
    sig_free(ef1)
    sig_free(ef2)
    sig_free(nvec)
    return W

@cython.boundscheck(True)

def phase_2(S,R,Cin,int NA,int NB,ndig=10,M0=None,Yin=None, cuspidal=True,sym_type=None,low_prec=False):
    r"""
    Computes more coefficients from the given initial set.

        INPUT:

            - ``S`` -- Space of Maas waveforms
            - ``R`` -- real (mpmath)
            - ``Cin`` -- complex vector (mpmath)
            - ``NA`` -- integer
            - ``NB`` -- integer
            - ``M0`` -- integer (default None) if set we only use this many coeffficients.
            - ``ndig`` -- integer (default 10) : want new coefficients with this many digits precision.
            - ``cuspidal`` -- (False) logical : Assume cusp forms
            - ``sym_type``  -- (None) integer

                -  0 -- even / sine
                -  1 -- odd  / cosine

            - ``low_prec`` -- (False) : Use low (double) precision K-Bessel

    OUTPUT:

    - ``Cout`` -- dictionary of Fourier coefficients

    EXAMPLE::


        sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
        sage: IR=mpmath.mpc(0,R)
        sage: M=MaassWaveForms(Gamma0(1))
        sage: C=Maassform_coeffs(M,R)
        sage: D=phase_2(SL2Z,R,C,2,10)



    """
    import mpmath
    mp_ctx=mpmath.mp
    pi=mpmath.mp.pi()
    mp0=mp_ctx.mpf(0)
    mp1=mp_ctx.mpf(1)
    mp2=mp_ctx.mpf(2)
    eps=mp_ctx.power(mp_ctx.mpf(10),mp_ctx.mpf(-ndig))
    twopi=mp2*pi
    G=S.group()
    nc=int(G.ncusps())
    IR=mp_ctx.mpc(0,R)
    if(M0<>None):
        M00=M0
        Mf=min ( M00, max(Cin[0].keys()))
    else:
        M00=infinity
        Mf= max(Cin[0].keys())
    if(sym_type<>None):
        Ms=1
    else:
        Ms=max ( -M00, min(Cin[0].keys()))

    Ml=Mf-Ms+1
    lvec=list()
    for l in range(Ms,Mf+1):
        lvec.append(mp_ctx.mpf(l))
    # starting value of Y
    if(Yin == None):
        Y0=twopi/mpmath.mp.mpf(NA)
        if(Y0>0.5):
            Y0=mp1/mp2
    else:
        Y0=Yin
    #
    NAa=NA
    ylp=0; Qp=0
    Cout=Cin # Start by assigning
    for yn in range(100):
        try:
            Yv=[Y0,Y0*mpmath.mpf(0.995)]

            Q=max(get_M_for_maass(R,Yv[1],eps)+5,Mf+10)+Qp
            Qs=1-Q; Qf=Q; Qfak=mp_ctx.mpf(2*Q)
            Xpb=dict(); Ypb=dict(); Cv=dict()
            pb=dict()
            pb[0]=pullback_pts_mpc(S,Qs,Qf,Yv[0])
            pb[1]=pullback_pts_mpc(S,Qs,Qf,Yv[1])
            Cnew=dict()
            print "Y[0]=",Yv,
            print "Ms,Mf=",Ms,Mf,"Qs,Qf=",Qs,Qf
            print "NA,NB=",NAa,NB
            kbdict=_setup_kbessel_dict(Ms,Mf,Qs,Qf,nc,IR,pb,lvec,mp_ctx)
            for n in range(NAa,NB+1):
                nr=mpmath.mp.mpf(n)
                for icusp in range(nc):
                    cvec=_get_tmpCs_twoY(n,icusp,Yv,IR,Ms,Mf,Qs,Qf,Qfak,nc,pb,Cin,kbdict,lvec,mp_ctx)
                    #                cvec=_get_tmpCs_twoY(n,icusp,Yv,IR,Ms,Mf,Qs,Qf,Qfak,nc,lvec,pb,Cin,kbdict,mp_ctx)
                    # print "c1,c2=",cvec
                    diff=abs(cvec[0]-cvec[1])
                    print "err[",n,"]=",diff
                    if(diff<eps):
                        Cnew[icusp,n]=cvec[1]
                        # # improve the coefficients we use
                        if(Cin.has_key((icusp,n))):
                            Cin[icusp,n]=cvec[1]
                    else:
                        NAa=n; ylp=ylp+1
                        if(ylp>2):
                            Qp=Qp+10; ylp=0
                        else:
                            Y0=Y0*mpmath.mpf(0.99)
                            #Qp=0

                        raise StopIteration()

        except StopIteration:
            continue
        return Cnew







def _setup_kbessel_dict(Ms,Mf,Qs,Qf,nc,IR,pb,lvec=None,mp_ctx=None):
    r"""
    Setup a dictionary of K-Bessel values.

    INPUT:

    -``Ms`` -- integer
    -``Mf`` -- integer
    -``Qs`` -- integer
    -``Qs`` -- integer
    -``nc`` -- integer
    -``IR`` -- complex (purely imaginary)
    -``pb`` -- dictionary
    -``lvec`` -- (optional) vector of mpmath real integers in range(Ms,Mf+1)
    -``mp_ctx`` -- mpmath context

    OUTPUT:
    -``kbdict`` dictionary of dimensions
        Ms:Mf,Qs:Qf,0:nc-1,0:nc-1,0:1

    EXAMPLES::


    Note that this function is not in the public name space.
    Therefore we need a workaround to call it from the command line.

        sage: Yv=[mpmath.mpf(0.5),mpmath.mpf(0.45)]
        sage: IR=mpmath.mpc(0,9.53369526135355755434423523592877032382125639510725198237579046413534)
        sage: M=20; Ms=-M; Mf=M; Q=30; Qs=1-Q; Qf=Q
        sage: pb=dict()
        sage: pb[0]=pullback_pts(Gamma0(1),Qs,Qf,Yv[0])
        sage: pb[1]=pullback_pts(Gamma0(1),Qs,Qf,Yv[1])
        sage: _temp = __import__('maass_forms_alg',globals(), locals(), ['_setup_kbessel_dict'],-1)
        sage: kbdict=_temp.__dict__['_setup_kbessel_dict'](Ms,Mf,Qs,Qf,1,IR,pb)
    """
    import mpmath
    if(mp_ctx==None):
        mp_ctx=mpmath.mp
    kbvec=dict()
    twopi=mp_ctx.mpf(2)*mp_ctx.pi()
    for l in range(Ms,Mf+1):
        if(lvec<>None):
            lr=lvec[l-Ms]
        else:
            lr=mp_ctx.mpf(l)
        for j in range(Qs,Qf+1):
            for icusp in range(nc):
                for jcusp in range(nc):
                    for yj in range(2):
                        ypb=pb[yj]['ypb'][icusp,jcusp,j]
                        ypbtwopi=ypb*twopi
                        besarg=abs(lr)*ypbtwopi
                        kbes=mp_ctx.sqrt(ypb)*(mp_ctx.besselk(IR,besarg).real)
                        kbvec[l,j,icusp,jcusp,yj]=kbes
    return kbvec


def _get_tmpCs_twoY(n,ic,Yv,IR,Ms,Mf,Qs,Qf,Qfak,nc,pb,Cin,kbdict,lvec=None,mp_ctx=None):
    r"""
    Compute two coeffcients given by the two Y-values.

    INPUT:

    -``n`` -- integer
    -``ic`` -- integer
    -``Yv`` -- list of two reals Y0,Y1
    -``IR`` -- complex  IR
    -``Ms`` -- integer
    -``Mf`` -- integer
    -``Qs`` -- integer
    -``Qf`` -- integer
    -``nc`` -- integer
    -``lvec`` -- list of mpmath reals in range(Ms,Mf)
    -``pb`` -- list of dictonaries of values related to the pullback
    -``Cin`` -- list of reals : Fourier coefficients
    -``kbdict-- dictionary of K-Besssel values
    -``mp_ctx`` -- mpmath context

    OUTPUT:

    --``cvec`` -- list of complex numbers : Fourier coefficients  c_1[n,ic] and c_2[n,ic]


    EXAMPLES::


    Note that this function is not in the public name space.
    Therefore we need a workaround to call it from the command line.
        sage: Yv=[mpmath.mpf(0.2),mpmath.mpf(0.17)]
        sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
        sage: IR=mpmath.mpc(0,R)
        sage: M=MaassWaveForms(Gamma0(1))
        sage: C=Maassform_coeffs(M,R)
        sage: M=max(C[0].keys()); Ms=-M; Mf=M
        sage: Q=max(get_M_for_maass(R,Yv[1],1E-12),M+10); Qs=1-Q; Qf=Q; nc=1; ic=0
        sage: Qfak=mpmath.mpf(2*Q)
        sage: pb=dict()
        sage: pb[0]=pullback_pts(Gamma0(1),Qs,Qf,Yv[0])
        sage: pb[1]=pullback_pts(Gamma0(1),Qs,Qf,Yv[1])
        sage: _temp = __import__('maass_forms_alg',globals(), locals(), ['_setup_kbessel_dict'],-1)
        sage: kbdict=_temp.__dict__['_setup_kbessel_dict'](Ms,Mf,Qs,Qf,nc,IR,pb)
        sage: n=10
        sage: c=_temp.__dict__['_get_tmpCs_twoY'](n,ic,Yv,IR,Ms,Mf,Qs,Qf,Qfak,nc,pb,C,kbdict)
        sage: C[0][n].real
        mpf('0.31053524297612709')
        sage: c[0]
        mpc(real='0.31053524291041912', imag='9.265227235203014e-17')
        sage: c[0]
        mpc(real='0.31053524291042567', imag='-3.8363741205534537e-17')
        sage:     sage: abs(c[0]-c[1])
        mpf('6.5516259713787926e-15')l


    """
    import mpmath
    if(mp_ctx==None):
        mp_ctx=mpmath.mp
    ctmp=dict(); tmpV=dict()
    nr=mp_ctx.mpf(n)
    mp0=mp_ctx.mpf(0)
    twopi=mp_ctx.mpf(2)*mp_ctx.pi()
    #print "keys=",Cin[0].keys()
    #print "Yv=",Yv
    for yj in range(2):
        Y=Yv[yj]
        summa=mp_ctx.mpf(0)
        for jc in range(nc):
            for l in range(Ms,Mf+1):
                if(Cin[jc][l]==0):
                    continue
                if(lvec<>None):
                    lr=lvec[l-Ms]
                else:
                    lr=mp_ctx.mpf(l)
                Vnlij=mp0
                for j in range(Qs,Qf+1):
                    if(not pb[yj]['ypb'].has_key((ic,jc,j))):
                        continue
                    ypb=pb[yj]['ypb'][ic,jc,j]
                    if(ypb==0):
                        continue
                    kbes=kbdict[l,j,ic,jc,yj]
                    xpb=pb[yj]['xpb'][ic,jc,j]
                    arg=lr*xpb-nr*pb[yj]['xm'][j]
                    tmp=mpmath.mp.exp(mpmath.mpc(0,arg))
                    Vnlij=Vnlij+kbes*tmp
                summa=summa+Vnlij*Cin[jc][l]
        tmpV[yj]=summa/Qfak
        besarg=nr*twopi*Y
        besY=mp_ctx.sqrt(Y)*(mp_ctx.besselk(IR,besarg).real)
        ctmp[yj]=tmpV[yj]/besY
        #print "summa[",yj,j,l,"]=",summa
        #print "tmpV[",yj,j,l,"]=",tmpV[yj]
        #print "C[",n,yj,"]=",ctmp[yj]
    return ctmp



def smallest_inf_norm_mpmath(V):
    r"""
    Computes the smallest of the supremum norms of the columns of a matrix.

    INPUT:

        - ``V`` -- matrix (real/complex)

    OUTPUT:

        - ``t`` -- minimum of supremum norms of the columns of V

    EXAMPLE::


        sage: A=mpmath.matrix([['0.1','0.3','1.0'],['7.1','5.5','4.8'],['3.2','4.4','5.6']])
        sage: smallest_inf_norm(A)
        mpf('5.5')


    """
    import mpmath
    minc=mpmath.mp.mpf(100)
    mi=0
    for j in range(V.cols):
        maxr=mpmath.mp.mpf(0)
        for k in range(V.rows):
            t=abs(V[k,j])
            if(t>maxr):
                maxr=t
        if(maxr<minc):
            minc=maxr
            mi=j
    return minc


def smallest_inf_norm(V):
    r"""
    Computes the smallest of the supremum norms of the columns of a matrix.

    INPUT:

        - ``V`` -- matrix (real/complex)

    OUTPUT:

        - ``t`` -- minimum of supremum norms of the columns of V

    EXAMPLE::


        sage: A=mpmath.matrix([['0.1','0.3','1.0'],['7.1','5.5','4.8'],['3.2','4.4','5.6']])
        sage: smallest_inf_norm(A)
        mpf('5.5')


    """
    RF = RealField(V[0,0].parent().prec())
    minc=RF(100)
    mi=0
    for j in range(V.ncols()):
        maxr=0
        for k in range(V.nrows()):
            t=abs(V[k,j])
            if t>maxr:
                maxr=t
        if maxr<minc:
            minc=maxr
            mi=j
    return minc

cpdef get_M_and_Y_v2(double R,double Y0,int M0,double eps,int mmax=1000,int verbose=0):
    r"""

    """
    R0 = max(R,1)
    kmax = besselk_dp(R,R0,pref=0)*eps
    try:
        for m in range(10,mmax):
            k = besselk_dp(R,m*Y0*2*M_PI,pref=0)
            if abs(k)<abs(kmax):
                raise StopIteration
        raise ArithmeticError,"Could not find suitable m0!"
    except StopIteration:
        pass
    if verbose>2:
        print "M0=",M0
    
    

        
cpdef get_M_and_Y(double R,double Y0,int M0,double eps,int mmax=10000,int verbose=0,int cuspidal=0):
    r""" Computes the  truncation point M>=M0 for Maass waveform with parameter R.
    using both a check that the diagonal term in the V-matrix is not too small
    and that the truncation gives the correct accuracy.

    CAVEAT: The estimate for the truncated series is done assuming that the asymptotic formula holds for 2piYM > R+12R^(1/3)
    which is not proven rigorously. Furthermore, the implied constant in this estimate is assumed to be 1 (which is of course only true asymptotically).

    INPUT:

        - ``R`` -- spectral parameter, real
        - ``Y`` -- height, real > 0
        - ``eps`` -- desired error


    OUTPUT:

        - ``M`` -- point for truncation

    EXAMPLES::


        sage: 
        sage: 


    """
    M = get_M_from_Y(R,Y0,M0,eps,maxm=mmax,verbose=verbose,cuspidal=cuspidal)
    return M,Y0
    # initial value of M
    MM=max(M0,ceil((12.0*R**0.3333+R)/(2*M_PI*Y0)))
    ## initial value of Y
    ## choosen s.t. we are to the right of the "hump" of the K-Bessel function
    Y = min(Y0, (log(2.0)-log(eps)-0.5*log(MM)+M_PI*R*0.5)/(2*M_PI*MM))
    if verbose>0:
        print "MM=",MM
        print "Y=",Y
    if M0>0 and Y>0: # check if we are ok to begin with
         if err_est_Maasswf(Y,M0,R,1)<eps:
             return M0,Y
         
    #if Y<0 or Y>Y0:
    #    raise ArithmeticError,"Could not get a good point"
    ## Now change M if necessary
    minm = ceil((12.0*R**0.3333+R)/(M_PI*Y*2.0))
    if verbose>0:
        print "minm=",minm
    ## Use low precision
    try:
        for m in range(minm+1,10000,3):
            if m>mmax and mmax>0:
                continue
            erest=err_est_Maasswf(Y,m,R,1)
            if verbose>2:
                print "erest({0},{1},{2},1)={3}".format(Y,m,R,erest)
            if erest<eps:
                YY = (log(2.0)-log(eps)-0.5*log(float(m))+M_PI*R*0.5)/(2*M_PI*float(m))
                YY = min(YY,Y0)
                erest=err_est_Maasswf(YY,m,R,1)
                #if verbose>2:                    
                #    print "erest({0},{1},{2},1)={3}".format(YY,m,R,erest)
                #if erest < eps:
                #    raise StopIteration()
    except StopIteration:
        Y = (log(2.0)-log(eps)-0.5*log(float(m))+M_PI*R*0.5)/(2*M_PI*float(m))                   
        Y = min(Y,Y0)
        return m,Y
    raise ArithmeticError," No good value for truncation was found!"

cpdef get_M_from_Y(double R,double Y0,int M0,double eps,int verbose=0,int maxm=1000,int cuspidal=0):
    r""" Computes the  truncation point M>=M0 for Maass waveform with parameter R.
    using both a check that the diagonal term in the V-matrix is not too small
    and that the truncation gives the correct accuracy.

    CAVEAT: The estimate for the truncated series is done assuming that the asymptotic formula holds for 2piYM > R+12R^(1/3)
    which is not proven rigorously. Furthermore, the implied constant in this estimate is assumed to be 1 (which is of course only true asymptotically).

    INPUT:

        - ``R`` -- spectral parameter, real
        - ``Y`` -- height, real > 0
        - ``eps`` -- desired error


    OUTPUT:

        - ``M`` -- point for truncation

    EXAMPLES::


        sage:

    """
    #maxm = 2000
    ## initial value of M
    M0 = ceil( 4.0/7.0/M_PI/Y0)
    if cuspidal==1:
        rhs = abs(eps)*sqrt(M_PI/8.0)
        beta = 2*M_PI*Y0
    else:
        rhs = abs(eps)*exp(-3*R)*sqrt(M_PI/8.0)
        beta = 7.0*M_PI*Y0/4.0
    minm=max(10,M0)
    if verbose>0:
        print "rhs={0}".format(rhs)
    try:
        for M in range(minm,minm+maxm):
            ## Want to find largest Y which satisfies both estimates for the given M0
            lhs = sqrt(float(M))*exp(-beta*M)
            if verbose>0:
                print "lhs({0})={1}".format(M,lhs)
            if lhs<rhs:
                raise StopIteration()
    except StopIteration:
        return M
    raise Exception," No good value for truncation was found!"

cpdef get_Y_from_M(double R,double Y0,int M0,double eps,int verbose=0,int maxny=1000,int cuspidal=0):
    r""" Computes the  truncation point M>=M0 for Maass waveform with parameter R.
    using both a check that the diagonal term in the V-matrix is not too small
    and that the truncation gives the correct accuracy.

    CAVEAT: The estimate for the truncated series is done assuming that the asymptotic formula holds for 2piYM > R+12R^(1/3)
    which is not proven rigorously. Furthermore, the implied constant in this estimate is assumed to be 1 (which is of course only true asymptotically).

    INPUT:

        - ``R`` -- spectral parameter, real
        - ``Y`` -- height, real > 0
        - ``eps`` -- desired error


    OUTPUT:

        - ``M`` -- point for truncation

    EXAMPLES::


        sage:

    """
    #maxm = 2000
    ## initial value of Y  s.t. 
    Y0 =  max(Y0,min(Y0,4.0/7.0/M_PI/float(M0)))
    if cuspidal==1:
        rhs = abs(eps)*sqrt(M_PI/8.0)/sqrt(float(M0))
        gamma = 2*M_PI
    else:
        rhs = abs(eps)*exp(-3*R)*sqrt(M_PI/8.0)/sqrt(float(M0))
        gamma = 7.0*M_PI/4.0
    Y = -RR(rhs).log()/float(M0)/gamma
    if verbose>0:
        print "rhs=",rhs
        print "Y0=",Y0
        print "Y=",Y
    return min(Y,Y0)

cpdef get_M_from_Y_old(double R,double Y0,int M0,double eps,int verbose=0,maxm=1000):
    r""" Computes the  truncation point M>=M0 for Maass waveform with parameter R.
    using both a check that the diagonal term in the V-matrix is not too small
    and that the truncation gives the correct accuracy.

    CAVEAT: The estimate for the truncated series is done assuming that the asymptotic formula holdsM0Y for 2piYM > R+12R^(1/3)
    which is not proven rigorously. Furthermore, the implied constant in this estimate is assumed to be 1 (which is of course only true asymptotically).

    INPUT:

        - ``R`` -- spectral parameter, real
        - ``Y`` -- height, real > 0
        - ``eps`` -- desired error


    OUTPUT:

        - ``M`` -- point for truncation

    EXAMPLES::


        sage:

    """
    #maxm = 2000
    ## initial value of M
    if M0<=0:
        M0 = 1
    minm = ceil( (log(2.0)-log(eps)-0.5*log(M0)+RR.pi()*R*0.5)/(2*RR.pi()*Y0))
    minm = max(M0,minm)
    minm= max(minm,ceil((12.0*R**0.3333+R)/RR.pi()/Y0/2.0))
    ## Use low precision
    twopiY0=Y0*RR.pi()*2
    if verbose>0:
        print "M0=",M0
        print "Y0=",Y0
        print "minm=",minm
    try:
        for M in range(minm,minm+maxm):
            ## Want to find largest Y which satisfies both estimates for the given M0
            erest = err_est_Maasswf(Y0,M,R,1)
            #t =besselk_dp(R,twopiY0*M,pref=1)
            t = (log(2.0)-log(eps)-0.5*log(M0)+M_PI*R*0.5)/(2*M_PI*Y0*M0)
            if verbose>0:
                print "erest=",erest
                print "t=",t
            #print "erest=",erest
            if erest<eps and t>1:
                raise StopIteration()
    except StopIteration:
        return M
    raise Exception," No good value for truncation was found!"


cpdef  get_Y_from_M_old(double R,double Y0,int M0,double eps,int verbose=0,int maxny=2000):
    r""" Computes the  truncation point M>=M0 for Maass waveform with parameter R.
    using both a check that the diagonal term in the V-matrix is not too small
    and that the truncation gives the correct accuracy.

    CAVEAT: The estimate for the truncated series is done assuming that the asymptotic formula holds for 2piYM > R+12R^(1/3)
    which is not proven rigorously. Furthermore, the implied constant in this estimate is assumed to be 1 (which is of course only true asymptotically).

    INPUT:

        - ``R`` -- spectral parameter, real
        - ``Y`` -- height, real > 0
        - ``eps`` -- desired error


    OUTPUT:

        - ``M`` -- point for truncation

    EXAMPLES::


        sage:

    """
    #maxny = 2000
    ## initial value of Y
    if Y0<=0:
        Y0 = (log(2.0)-log(eps)-0.5*log(M0)+M_PI*R*0.5)/(2*M_PI*M0)
    #if Y<0 or Y>Y0:
    #    raise ArithmeticError,"Could not get a good point"
    ## Now change M 
    #minm=ceil((12.0*R**0.3333+R)/RR.pi()/Y/2.0)
    ## Use low precision
    Y = Y0
    twopim0=M0*RR.pi()*2
    if verbose>0:
        print "M0=",M0
    try:
        for yi in range(maxny):
            Y = Y*0.995
            ## Want to find largest Y which satisfies both estimates for the given M0
            erest = err_est_Maasswf(Y,M0,R,1)
            #Y0 = (log(2.0)-log(eps)-0.5*log(M0)+M_PI*R*0.5)/(2*M_PI*M0)
            #t =besselk_dp(R,twopim0*Y,pref=1)
            #t = (log(2.0)-log(eps)-0.5*log(M0)+RR.pi()*R*0.5)/(2*RR.pi()*M0)
            #print "erest=",erest
            if verbose>0:
                print "Y=",Y
                print "errest=",erest
                print "Y0=",Y0
            if erest<eps and Y<Y0:
                raise StopIteration()
    except StopIteration:
        return Y
    raise ValueError," No good value for Y was found. Please increase M0!"



def get_M_for_maass(R,Y,eps):
    r""" Computes the  truncation point for Maass waveform with parameter R.

    INPUT:

        - ``R`` -- spectral parameter, real
        - ``Y`` -- height, real > 0
        - ``eps`` -- desired error


    OUTPUT:

        - ``M`` -- point for truncation

    EXAMPLES::ubuntuusers.de


        sage: get_M_for_maass(9.533,0.5,1E-16)
        12
        sage: get_M_for_maass(9.533,0.5,1E-25)
        18


    """
    ## Use low precision
    import mpmath
    dold=mpmath.mp.dps
    mpmath.mp.dps=int(mpmath.ceil(abs(mpmath.log10(eps))))+5
    twopi=2*mpmath.mp.pi()
    twopiy=twopi*mpmath.mp.mpf(Y)
    # an extra two for the accumulation of errors
    minv=(mpmath.mp.mpf(12)*mpmath.mp.power(R,0.3333)+R)/twopiy
    minm=mpmath.mp.ceil(minv)
    if mpmath.mp.isnan(minm):
        print "minm=",minm
        return -1
    minv = int(minm)
    #print "eps=",eps

    try:
        for m in range(minm+1,10000,3):
            erest=err_est_Maasswf_c(Y,m,R,1)
            #print "erest=",erest
            if erest<eps:
                raise StopIteration()
    except StopIteration:
        mpmath.mp.dps=dold
        return m
    mpmath.mp.dps=dold
    raise Exception," No good value for truncation was found!"


def err_est_Maasswf(Y,M,R=0,pref=1,proof=False):
    r"""
    Estimate the truncated series: 
      $|\Sum_{n\ge M} a(n) \sqrt{Y} K_{iR}(2\pi nY) e(nx)|$
    we use that K_{iR}(x) \sim sqrt(pi/2x)e^-x to estimate the sum with
     O_R(1) sqrt(pi/2)  |$\sum_{n\ge M} (2piny)^(-1/2)e^-(2piny)|
     which we then estimate by an integral and get: 
     sqrt(pi/2)  Gamma(3/2,2piMY) 
    and finally, for Y>1/(2piM) we have |Gamma(3/2,2piMY)|< 2x^(1/2)e^-x

    so that the final bound is 
    2O_R(1) M^(1/2) e^(-2piMY) 
    CAVEATS:
    - we assume the Ramanujan bound for the coefficients $a(n)$, i.e. $|a(n)|\le2$.
    - error estimate holds for 2piMY >> R
    - there is an implicit constant dependning on R which we neglect if proof=False
      (currently this is the only value supported)

    INPUT:

    - ``Y`` -- real > 0
    - ``M``  -- integer >0

    OUTPUT:

    - ``r``  -- error estimate

    EXAMPLES::


        sage: err_est_Maasswf(0.5,10)
        mpf('3.1837954148201557e-15')
        sage: err_est_Maasswf(0.85,10)
        mpf('5.3032651127544969e-25')

    """
    if proof==True:
        raise NotImplementedError
    import mpmath
    mpmath.mp.pres=103
    if R < 100: ## ordinary doubles should be sufficient
        if pref==1:
            f = exp(M_PI*R*0.5)
        else:
            f = 1
            #    r = f*mpmath.fp.erfc(arg)
        r = f*2*M**(0.5)*exp(-2*M_PI*Y*M)

    else:
        YY = mpmath.mp.mpf(Y)
        #arg=sqrt(2*mpmath.fp.pi*YY*M)
        #gamma_arg=2*mpmath.mp.pi*YY*M
        if pref==1:
            f = exp(mpmath.mp.pi*R/2.0)
            #*mpmath.mp.sqrt(2.0/YY)
        else:
            f = 1
            #    r = f*mpmath.fp.erfc(arg)
        r = f*2*M**(0.5)*exp(-2*M_PI*YY*M)
    #mpmath.mp.gammainc(0.75,gamma_arg)
    return r

cpdef err_est_Maasswf_c(double Y,int M,double R=0,int pref=1):
    if R < 100: ## ordinary doubles should be sufficient
        if pref==1:
            f = exp(M_PI*R*0.5)
        else:
            f = 1
            #    r = f*mpmath.fp.erfc(arg)
        r = f*2*M**(0.5)*exp(-2*M_PI*Y*M)

    else:
        YY = mpmath.mp.mpf(Y)
        #arg=sqrt(2*mpmath.fp.pi*YY*M)
        #gamma_arg=2*mpmath.mp.pi*YY*M
        if pref==1:
            f = exp(mpmath.mp.pi*R/2.0)
            #*mpmath.mp.sqrt(2.0/YY)
        else:
            f = 1
            #    r = f*mpmath.fp.erfc(arg)
        r = f*2*M**(0.5)*exp(-2*M_PI*YY*M)
    #mpmath.mp.gammainc(0.75,gamma_arg)
    return r


# cpdef get_Y_and_M_dp(S,double R,double eps,int verbose=0):
#     cdef int i,M0,MMAX=1000,m
#     cdef double Y
#     M0 = S.smallest_M0()
#     for m in range(M0,M0*MMAX):
#         Y=get_Y_for_M_dp(R,m,eps,S.group().minimal_height(),S.group().ncusps(),verbose=verbose)
#         if verbose>1:
#             print "M={0}, Y={1}".format(m,Y)
#         if Y<0:
#             continue
#         MM = get_M_for_maass_dp_c(R,Y,eps)
#         MM +=10
#         return Y,MM
#     raise ArithmeticError,"Could not find good Y and M for eps={0}".format(eps)

cpdef get_Y_for_M(double R,int M0,double eps,double Y0,int verbose=0,int cuspidal=0):
    r"""
    Computes a value of Y for a given M for which
    the estimated error is smaller than eps
    """
    Y = Y0
    for i in range(1,1000):
        Y = Y*0.999**i
        M=get_M_from_Y(R,Y,1,eps,cuspidal=cuspidal)
        if M > M0:
            return Y
    raise ArithmeticError,"Could not find a good Y!"
        
cpdef get_Y_for_M_dp_old(double R,int M,double eps,double Y,int nc,int verbose=0):
    r"""
    Computes a value of Y for a given M for which
    the estimated error is smaller than eps
    """
    #cdef double Y=S._group.minimal_height()
    #cdef int nc = S._group._ncusps
    cdef int yi,maxny=2000
    cdef double errest,errest_old,minv,twopim
    cdef int minm,m
    twopim=<double>2.0*M_PI*<double>M
    if R<0:
        raise ValueError,"Need non-negative R! Got R={0}".format(R)
    if M<=0:
        raise ValueError,"Need positive M! Got M={0}".format(M)
    Y = Y
    # Recall that the error estimate is only valid if 2piYM>>R
    # this is a good lower bound for this inequality
    minv=(12.0*R**0.3333+R)/twopim
    errest_old=1
    for yi in range(maxny):
        Y = Y*0.995
        errest=nc*err_est_Maasswf_dp_c(Y,M,R,1)
        if verbose>0:
            print "errest({0},{1})={2}".format(Y,M,errest)
        if errest<eps:
            break
        if errest>errest_old+eps:
            return -1  ### We need to get a larger Y
        if Y < minv:
            return -1
        
        errest_old=errest
    return Y

cpdef get_M_for_maass_dp(double R,double Y,double eps):
    if R<0:
        raise ValueError,"Need non-negative R! Got R={0}".format(R)
    return get_M_for_maass_dp_c(R,Y,eps)

@cython.cdivision(True)
cdef int get_M_for_maass_dp_c(double R,double Y,double eps):
    r""" Computes the  truncation point for Maass waveform with parameter R.

    INPUT:

        - ``R`` -- spectral parameter, real
        - ``Y`` -- height, real > 0
        - ``eps`` -- desired error


    OUTPUT:

        - ``M`` -- point for truncation

    EXAMPLES::ubuntuusers.de


        sage: get_M_for_maass(9.533,0.5,1E-16)
        12
        sage: get_M_for_maass(9.533,0.5,1E-25)
        18


    """

    cdef double errest,minv,twopiy
    cdef int minm,m
    twopiy=<double>2.0*M_PI*Y
    # an extra two for the accumulation of errors
    minv=(12.0*R**0.3333+R)/twopiy
    minm=ceil(minv)
    try:
        for m in range(minm+1,200000,3):
            errest=err_est_Maasswf_dp_c(Y,m,R,1)
            #print "erest=",erest
            if errest<eps:
                raise StopIteration()
    except StopIteration:
        return m
    raise ArithmeticError," No good value for truncation was found! For R={0}, Y={1} and eps={2}".format(R,Y,eps)

@cython.cdivision(True)
cdef double err_est_Maasswf_dp_c(double Y,int M,double R,int pref=1):
    r"""
    Estimate the truncated series: $|\Sum_{n\ge M} a(n) \sqrt{Y} K_{iR}(2\pi nY) e(nx)|$

    CAVEATS:
    - we assume the Ramanujan bound for the coefficients $a(n)$, i.e. $|a(n)|\le2$.
    - error estimate holds for 2piMY >> R
    - this only uses estimates on the Fourier series at infinity
    INPUT:

    - ``Y`` -- real > 0
    - ``M``  -- integer >0
    - ``R`` -- float >0
    - ``pref`` -- integer, set to 1 if we use the normalized K-Bessel function: e^(pi*R/2)*K_iR(x)
    OUTPUT:

    - ``r``  -- error estimate

    EXAMPLES::


        sage: err_est_Maasswf(0.5,10)
        mpf('3.1837954148201557e-15')
        sage: err_est_Maasswf(0.85,10)
        mpf('5.3032651127544969e-25')

    """
    cdef double r,arg=M_PI*Y
    arg=M_PI*Y*<double>(2*M)
    #cdef double incgg  # = Gamma(1/2,arg)
    if pref==1:
        r = exp(R*M_PI/2.0)/sqrt(2.0*Y)
    else:
        r = 1.0/sqrt(2.0*Y)
    return erfc(sqrt(arg))*r
    # mpmath.fp.gammainc(mpmath.fp.mpf(0.5),arg)*sqrt(pi)/sqrt(arg)


cpdef solve_system_for_Maass_waveforms_cp(W,N=None,deb=False):
    r"""
    Solve the linear system to obtain the Fourier coefficients of Maass forms

    INPUT:

    - ``W`` --   (system) dictionary

        - ``W['Ms']``  -- M start
        - ``W['Mf']``  -- M stop
        - ``W['nc']``  -- number of cusps
        - ``W['V']``   -- matrix of size ((Ms-Mf+1)*nc)**2
        - ``W['RHS']`` -- right hand side (for inhomogeneous system) matrix of size ((Ms-Mf+1)*nc)*(dim)

    - ``N`` -- normalisation (dictionary, output from the set_norm_for_maass function)

        - ``N['SetCs']``   -- Which coefficients are set
        - ``N['Vals'] ``   -- To which values are these coefficients set
        - ``N['comp_dim']``-- How large is the assumed dimension of the solution space

    - ``deb`` -- print debugging information (default False)

    OUTPUT:

    - ``C`` -- Fourier coefficients

    EXAMPLES::

        sage: S=MaassWaveForms(MySubgroup(Gamma0(1)))
        sage: mpmath.mp.dps=20
        sage: R=mpmath.mpf(9.533695261353557554344235235928770323821256395107251982375790464135348991298347781769255509975435366)
        sage: Y=mpmath.mpf(0.5)
        sage: W=setup_matrix_for_Maass_waveforms(S,R,Y,12,22)
        sage: N=S.set_norm_maass(1)
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
    import mpmath
    V=W['V']
    Ms=W['Ms']
    Mf=W['Mf']
    nc=W['nc']
    Ml=Mf-Ms+1
    M=W['space']
    if N==None:
        N=M.set_norm()
    if hasattr(V,'shape'):
        nrows,ncols=V.shape
    else:
        nrows=V.rows
        ncols=V.cols
    if ncols<>Ml*nc or nrows<>Ml*nc:
        raise Exception," Wrong dimension of input matrix!"
    SetCs=N['SetCs']
    Vals=N['Vals']
    comp_dim=N['comp_dim']
    if N['cuspidal']:
        for fn_j in range(comp_dim):
            for i in range(1,nc):
                if SetCs[fn_j].count(i*Ml)==0:
                    SetCs[fn_j].append(i*Ml)
                    Vals[fn_j][i*Ml]=0

    if Ms<0:
        use_sym=0
    else:
        use_sym=1
    if(use_sym==1 and SetCs.count(0)>0):
        num_set=len(N['SetCs'])-1
    else:
        num_set=len(N['SetCs'])
    t=V[0,0]
    if(isinstance(t,float)):
        mpmath_ctx=mpmath.fp
    else:
        mpmath_ctx=mpmath.mp
    if(W.has_key('RHS')):
        RHS=W['RHS']
    else:
        RHS=mpmath_ctx.matrix(int(Ml*nc-num_set),int(comp_dim))
    LHS=mpmath_ctx.matrix(int(Ml*nc-num_set),int(Ml*nc-num_set))
    roffs=0

    if(deb):
        print "num_set,use_sym=",num_set,use_sym
        print "SetCs,Vals=",SetCs,Vals
        print "V.rows,cols=",nrows,ncols
        print "LHS.rows,cols=",LHS.rows,LHS.cols
        print "RHS.rows,cols=",RHS.rows,RHS.cols
    for r in range(nrows):
        cr=r+Ms
        if(SetCs.count(r+Ms)>0):
            roffs=roffs+1
            continue
        for fn_j in range(comp_dim):
            RHS[r-roffs,fn_j]=mpmath_ctx.mpf(0)
            for cset in SetCs:
                v=Vals[fn_j][cset]
                if(mpmath_ctx==mpmath.mp):
                    tmp=mpmath_ctx.mpmathify(v)
                elif(isinstance(v,float)):
                    tmp=mpmath_ctx.mpf(v)
                else:
                    tmp=mpmath_ctx.mpc(v)
                #print "tmp=",tmp
                #print "V[",r,cset-Ms,"]=",V[r,cset-Ms]
                tmp=tmp*V[r,cset-Ms]
                RHS[r-roffs,fn_j]=RHS[r-roffs,fn_j]-tmp
        coffs=0
        for k in range(ncols):
            if(SetCs.count(k+Ms)>0):
                coffs=coffs+1
                continue
	    #print "r,k=",r,k
            #print "roffs,coffs=",roffs,coffs
            #print "r-roffs,k-coffs=",r-roffs,k-coffs
            LHS[r-roffs,k-coffs]=V[r,k]
            #print "LHS[",r,k,"]=",LHS[r-roffs,k-coffs]
    #print "RHS="
    #for j in range(RHS.rows):
    #    print j,RHS[j,0]#
    #return LHS
    try:
        A, p = mpmath_ctx.LU_decomp(LHS)
    except ZeroDivisionError:
        t1=smallest_inf_norm(LHS)
        print "n=",smallest_inf_norm(LHS)
        t2=mpmath_ctx.log10(smallest_inf_norm(LHS))
        t3=mpmath_ctx.ceil(-t2)
        isinf=False
        if(isinstance(t3,float)):
            isinf = (t3 == float(infinity))
        if(isinstance(t3,sage.libs.mpmath.ext_main.mpf)):
            isinf = ((t3.ae(mpmath.inf)) or t3==mpmath.inf)
        if(isinstance(t3,sage.rings.real_mpfr.RealLiteral)):
            isinf = t3.is_infinity()
        if(isinf):
            raise ValueError, " element in LHS is infinity! t3=%s" %t3
        #print "LHS="
        #for j in range(LHS.rows):
        #    print j,LHS.column(j)
        #print "t3=",t3
        t=int(t3)
        #t=int(mpmath_ctx.ceil(-mpmath_ctx.log10(smallest_inf_norm(LHS))))
        #raise ZeroDivisionError,"Need higher precision! Use > %s digits!" % t
        return t
    X=dict()
    for fn_j in range(comp_dim):
        X[fn_j] = dict() #mpmath.matrix(int(Ml),int(1))
        b = mpmath_ctx.L_solve(A, RHS.column(fn_j), p)
        TMP = mpmath_ctx.U_solve(A, b)
        roffs=0
        res = mpmath_ctx.norm(mpmath_ctx.residual(LHS, TMP, RHS.column(fn_j)))
        #print "res(",fn_j,")=",res
        for i in range(nc):
            X[fn_j][i]=dict()
        for n in range(Ml):
            if(SetCs.count(n+Ms)>0):
                roffs=roffs+1
                #print "X[",fn_j,",",n,",Vals[fn_j][n]
                X[fn_j][0][n+Ms]=Vals[fn_j][n+Ms]
                continue
            X[fn_j][0][n+Ms]=TMP[n-roffs,0]
            #print "X[",fn_j,",",n+Ms,"=",TMP[n-roffs,0]
        for i in range(1,nc):
            for n in range(Ml):
                if(SetCs.count(n+Ms+i*Ml)>0):
                    roffs=roffs+1
                    # print "X[",fn_j,",",n,",Vals[fn_j][n]
                    X[fn_j][i][n+Ms]=Vals[fn_j][n+Ms+i*Ml]
                    continue
                X[fn_j][i][n+Ms]=TMP[n+i*Ml-roffs,0]
    # return x
    return X

from sage.all import vector

# cpdef solve_system_for_Maass_waveforms_mpc2(W,N=None,gr=False,cn=False):
#     r"""
#     Solve the linear system to obtain the Fourier coefficients of Maass forms

#     INPUT:

#     - ``H`` -- Space of Maass waveforms
#     - ``W`` --   (system) dictionary
#         - ``W['Ms']``  -- M start
#         - ``W['Mf']``  -- M stop
#         - ``W['nc']``  -- number of cusps
#         - ``W['V']``   -- matrix of size ((Ms-Mf+1)*nc)**2
#         - ``W['RHS']`` -- right hand side (for inhomogeneous system) matrix of size ((Ms-Mf+1)*nc)*(dim)
#     - ``N`` -- normalisation (dictionary, output from the set_norm_for_maass function, default None)
#         - if N=None we assume that the solution is uniquely determined by the prinicpal part (in the right hand side)
#         - ``N['SetCs']``   -- Which coefficients are set
#         - ``N['Vals'] ``   -- To which values are these coefficients set
#         - ``N['comp_dim']``-- How large is the assumed dimension of the solution space
#         - ``N['num_set']`` -- Number of coefficients which are set



#     - ''gr''  -- only return the reduced matrix and right hand side. do not perform the solving .
#     - ''cn''  -- logical (default False) set to True to compute the max norm of V^-1
#     OUTPUT:

#     - ``C`` -- Fourier coefficients

#     EXAMPLES::

#         sage: G=MySubgroup(Gamma0(1))
#         sage: mpmath.mp.dps=20
#         sage: R=mpmath.mpf(9.533695261353557554344235235928770323821256395107251982375790464135348991298347781769255509975435366)
#         sage: Y=mpmath.mpf(0.5)
#         sage: W=setup_matrix_for_Maass_waveforms(G,R,Y,12,22)
#         sage: N=set_norm_maass(1)
#         sage: C=solve_system_for_Maass_waveforms(W,N)
#         sage: C[0][2]*C[0][3]-C[0][6]
#         mpc(real='-1.8055426724989656270259e-14', imag='1.6658248366482944572967e-19')

#     If M is too large and the precision is not high enough the matrix might be numerically singular

#         W=setup_matrix_for_Maass_waveforms(G,R,Y,20,40)
#         sage: C=solve_system_for_Maass_waveforms(W,N)
#         Traceback (most recent call last)
#         ...
#         ZeroDivisionError: Need higher precision! Use > 23 digits!

#     Increasing the precision helps

#         sage: mpmath.mp.dps=25
#         sage: R=mpmath.mpf(9.533695261353557554344235235928770323821256395107251982375790464135348991298347781769255509975435366)
#         sage: C=solve_system_for_Maass_waveforms(W,N)
#         sage: C[0][2]*C[0][3]-C[0][6]
#         mpc(real='3.780824715556976438911480324e-25', imag='2.114746048869188750991752872e-99')

#     """
#     cdef Matrix_complex_dense LHS,RHS,Q,R
#     cdef int nrows,ncols
#     cdef int i,j,r,n,k,roffs,coffs,fn_j,Ms,Mi,Ml,nc,verbose,comp_dim,num_set,rhs_j,num_rhs
#     cdef Vector_complex_dense y,v,b
#     V=W['V']
#     Ms=W['Ms']
#     Mf=W['Mf']
#     H=W['space']
#     #nc=W['nc']
#     Ml=Mf-Ms+1
#     get_reduced_matrix=gr
#     verbose = H._verbose
#     comp_norm=cn
#     nc=H.group().ncusps()
#     if V.ncols()<>Ml*nc or V.nrows()<>Ml*nc:
#         raise Exception," Wrong dimension of input matrix!"
#     if N==None:
#         N = H.set_norm(1)
#     SetCs=N['SetCs'][0]
#     Vals=N['Vals']
#     comp_dim=N['comp_dim']
#     num_set=len(SetCs[0])
#     cdef MPComplexField CF
#     CF=MPComplexField(H._prec)
#     MS = MatrixSpace(CF,int(Ml*nc-num_set),int(comp_dim))
#     RHS=Matrix_complex_dense(MS,0,True,True)
#     MS = MatrixSpace(CF,int(Ml*nc-num_set),int(Ml*nc-num_set))
#     LHS=Matrix_complex_dense(MS,0,True,True)
#     ncols=LHS.ncols()
#     nrows=LHS.nrows()
#     if(N['cuspidal']):
#         for i in range(1,nc):
#             if(SetCs.count((i,0))==0):
#                 SetCs.append((i,Ml))
#             for fn_j in range(comp_dim):
#                 Vals[fn_j][(i,Ml)]=0
#     if verbose>0:
#         print "SetCs=",SetCs
#     setc_list=list()
#     vals_list=dict()
#     for j in range(comp_dim):
#         vals_list[j]=dict()
#     for r,n in SetCs:
#         if r*Ml+n-Ms<0:
#             continue
#         setc_list.append(r*Ml+n-Ms)
#         for j in range(comp_dim):
#             vals_list[j][r*Ml+n-Ms]=Vals[j][(r,n)]

#     if verbose>0:
#         print "Ml=",Ml
#         print "num_set=",num_set
#         print "SetCs=",SetCs
#         print "Vals=",Vals
#         print "setc_list=",setc_list
#         print "vals_list=",vals_list
#         print "V.rows=",V.nrows()
#         print "V.cols=",V.ncols()
#         print "LHS.rows=",LHS.nrows()
#         print "LHS.cols=",LHS.ncols()
#         print "RHS.rows=",RHS.nrows()
#         print "RHS.cols=",RHS.ncols()
#         print "N=",N

#     num_rhs=0
#     if(W.has_key('RHS')):
#         num_rhs=W['RHS'].ncols()
#     if num_rhs>0 and num_rhs<>comp_dim:
#         raise ValueError,"Need same number of right hand sides (or just one) as the number of set coefficients!"
#     if V.nrows() <> nc*Ml:
#         raise ArithmeticError," Matrix does not have correct size!"
#     roffs=0
#     cdef MPComplexNumber tmp
#     tmp=CF(0)
#     for r in range(nrows):
#         #cr=r+Ms
#         if setc_list.count(r)>0:
#             roffs=roffs+1
#             continue
#         for fn_j in range(comp_dim):
#             if W.has_key('RHS'):
#                 if comp_dim==num_rhs:
#                     rhs_j=fn_j
#                 else:
#                     rhs_j=0
#                 RHS[r-roffs,fn_j]=-W['RHS'][r,rhs_j]
#             else:
#                 RHS[r-roffs,fn_j]=CF(0)
#             for cset in setc_list:
#                 tmp=CF(vals_list[fn_j][cset])
#                 tmp=tmp*V[r,cset]
#                 #RHS._matrix[r-roffs][fn_j]=RHS._matrix[r-roffs][fn_j]-tmp
#                 mpc_sub(RHS._matrix[r-roffs][fn_j],RHS._matrix[r-roffs][fn_j],tmp.value,rnd)
#         coffs=0
#         for k in range(ncols):
#             if setc_list.count(k)>0:
#                 coffs=coffs+1
#                 continue
#             mpc_set(LHS._matrix[r-roffs][k-coffs],(<Matrix_complex_dense>V)._matrix[r][k],rnd)    # for a in range(nc):
#     if get_reduced_matrix:
#         return [LHS,RHS]
#     cdef int maxit,done,dps,dps0
#     maxit=100;  i=0
#     done=0
#     dps0=CF.prec()
#     while (not done and i<=maxit):
#         try:
#             Q,R = LHS.qr_decomposition()
#             done=1
#         except ZeroDivisionError:
#             t=int(ceil(-log_b(smallest_inf_norm(LHS),10)))
#             dps=t+5*i; i=i+1
#             if verbose>0:
#                 print "raising number of digits to:",dps
#             LHS.set_prec(dps)
#             # raise ZeroDivisionError,"Need higher precision! Use > %s digits!" % t
#     if i>=maxit:
#         raise ZeroDivisionError,"Can not raise precision enough to solve system! Should need > %s digits! and %s digits was not enough!" % (t,dps)
#     if comp_norm:
#         max_norm=LHS.norm()
#         for j in range(LHS.rows):
#             #y=mpmath_ctx.matrix(LHS.rows,int(1)); y[j,0]=1
#             y = Vector_complex_dense(vector(CF,LHS.rows).parent(),0)
#             y[j]=1
#             TMP = RHS.solve(b) #pmath_ctx.U_solve(A, b)
#             tmpnorm=max(map(abs,TMP))
#             if tmpnorm>max_norm:
#                 max_norm=tmpnorm
#         if verbose>0:
#             print "max norm of V^-1=",max_norm
#     X=dict()
#     for fn_j in range(comp_dim):
#         X[fn_j] = dict() #mpmath.matrix(int(Ml),int(1))
#         v = RHS.column(fn_j)
#         if verbose>1:
#             print "len(B)=",len(v)
#         TMP = LHS.solve(v)
#         roffs=0
#         res = (LHS*TMP-v).norm()
#         if verbose>0:
#             print "res(",fn_j,")=",res
#         for i in range(nc):
#             X[fn_j][i]=dict()
#         for i in range(nc):
#             for n in range(Ml):
#                 if setc_list.count(n+i*Ml)>0:
#                     roffs=roffs+1
#                     X[fn_j][i][n+Ms]=vals_list[fn_j][n+i*Ml]
#                     continue
#                 X[fn_j][i][n+Ms]=TMP[n+i*Ml-roffs]
#     return X



# cpdef reduce_system_cplx(W):
#     r"""
#     We apply this method to reduce the system with respect to set coefficients etc.
#     """
#     cdef Matrix_complex_dense LHS,RHS,Q,R
#     cdef int nrows,ncols
#     cdef int i,j,r,n,k,roffs,coffs,fn_j,Ms,Mi,Ml,nc,verbose,comp_dim,num_set,rhs_j
#     cdef Vector_complex_dense y,v,b
#     V=W['V']
#     Ms=W['Ms']
#     Mf=W['Mf']
#     H=W['space']
#     #nc=W['nc']
#     Ml=Mf-Ms+1
#     #get_reduced_matrix=gr
#     verbose = H._verbose
#     #comp_norm=cn
#     nc=H.group().ncusps()
#     if V.ncols()<>Ml*nc or V.nrows()<>Ml*nc:
#         raise Exception," Wrong dimension of input matrix!"
#     if N==None:
#         N = H.set_norm(1)
#     SetCs=N['SetCs'][0]
#     Vals=N['Vals']
#     comp_dim=N['comp_dim']
#     num_set=len(SetCs[0])
#     CF=MPComplexField(H._prec)
#     MS = MatrixSpace(CF,int(Ml*nc-num_set),int(comp_dim))
#     RHS=Matrix_complex_dense(MS,0,True,True)
#     MS = MatrixSpace(CF,int(Ml*nc-num_set),int(Ml*nc-num_set))
#     LHS=Matrix_complex_dense(MS,0,True,True)
#     ncols=LHS.ncols()
#     nrows=LHS.nrows()
#     if(N['cuspidal']):
#         for i in range(1,nc):
#             if(SetCs.count((i,0))==0):
#                 SetCs.append((i,Ml))
#             for fn_j in range(comp_dim):
#                 Vals[fn_j][(i,Ml)]=0
#     if verbose>0:
#         print "SetCs=",SetCs
#     setc_list=list()
#     vals_list=dict()
#     for j in range(comp_dim):
#         vals_list[j]=dict()
#     for r,n in SetCs:
#         if r*Ml+n-Ms<0:
#             continue
#         setc_list.append(r*Ml+n-Ms)
#         for j in range(comp_dim):
#             vals_list[j][r*Ml+n-Ms]=Vals[j][(r,n)]

#     if verbose>0:
#         print "Ml=",Ml
#         print "num_set=",num_set
#         print "SetCs=",SetCs
#         print "Vals=",Vals
#         print "setc_list=",setc_list
#         print "vals_list=",vals_list
#         print "V.rows=",V.nrows()
#         print "V.cols=",V.ncols()
#         print "LHS.rows=",LHS.nrows()
#         print "LHS.cols=",LHS.ncols()
#         print "RHS.rows=",RHS.nrows()
#         print "RHS.cols=",RHS.ncols()
#         print "N=",N

#     num_rhs=0
#     if(W.has_key('RHS')):
#         num_rhs=W['RHS'].ncols()
#     if num_rhs>0 and num_rhs<>comp_dim:
#         raise ValueError,"Need same number of right hand sides (or just one) as the number of set coefficients!"
#     if V.nrows() <> nc*Ml:
#         raise ArithmeticError," Matrix does not have correct size!"
#     roffs=0
#     cdef MPComplexNumber tmp
#     tmp=CF(0)
#     for r in range(nrows):
#         #cr=r+Ms
#         if setc_list.count(r)>0:
#             roffs=roffs+1
#             continue
#         for fn_j in range(comp_dim):
#             if W.has_key('RHS'):
#                 if comp_dim==num_rhs:
#                     rhs_j=fn_j
#                 else:
#                     rhs_j=0
#                 RHS[r-roffs,fn_j]=-W['RHS'][r,rhs_j]
#             else:
#                 RHS[r-roffs,fn_j]=CF(0)
#             for cset in setc_list:
#                 tmp=CF(vals_list[fn_j][cset])
#                 tmp=tmp*V[r,cset]
#                 #RHS._matrix[r-roffs][fn_j]=RHS._matrix[r-roffs][fn_j]-tmp
#                 mpc_sub(RHS._matrix[r-roffs][fn_j],RHS._matrix[r-roffs][fn_j],tmp.value,rnd)
#         coffs=0
#         for k in range(ncols):
#             if setc_list.count(k)>0:
#                 coffs=coffs+1
#                 continue
#             mpc_set(LHS._matrix[r-roffs][k-coffs],(<Matrix_complex_dense>V)._matrix[r][k],rnd)    # for a in range(nc):






cpdef solve_using_Gauss_elem(LHS,RHS):
    cdef double complex **U=NULL
    cdef double complex **C=NULL
    cdef int n,m,i,j
    cdef int typ=0
    if hasattr(LHS,"nrows"):
        n=LHS.nrows()
        m=LHS.ncols()
    elif hasattr(LHS,"rows"):
        n=LHS.rows
        m=LHS.cols
        typ=1
    else:
        raise ValueError,"Incorrect matrix format in input! Got:{0}".format(type(LHS))
    assert m==n
    assert len(RHS)==n
    #assert RHS.ncols()==1
    U=<double complex**>sig_malloc(sizeof(double complex*)*n)
    if U==NULL:
        raise MemoryError
    if not isinstance(LHS[0,0],(complex,float)):
        for i in range(n):
            U[i]=<double complex*>sig_malloc(sizeof(double complex)*(n+1))
            if typ==0:
                for j in range(n):
                    U[i][j]=<double complex> CC(LHS[i,j].real(),LHS[i,j].imag())
                U[i][n]=<double complex> CC(RHS[i,0].real(),RHS[i,0].imag())
            else:
                for j in range(n):
                    U[i][j]=<double complex> CC(LHS[i,j].real,LHS[i,j].imag)
                    #print "V[",i,j,"]=",U[i][j],"=",LHS[i,j]
                U[i][n]=<double complex> CC(RHS[i,0].real,RHS[i,0].imag)
    else:
        for i in range(n):
            U[i]=<double complex*>sig_malloc(sizeof(double complex)*(n+1))
            for j in range(n):
                U[i][j]=<double complex> LHS[i,j]
            U[i][n]=<double complex> RHS[i]
    C=<double complex**>sig_malloc(sizeof(double complex*)*1)
    C[0]=<double complex*>sig_malloc(sizeof(double complex)*n)
    for i in range(n):
        print "V[",i,n,"]=",U[i][n]
    cdef double complex** values
    cdef int* cset
    cset=<int*>sig_malloc(sizeof(int)*2)
    cset[0]=0
    cset[0]=1

    values=<double complex**>sig_malloc(sizeof(double complex*)*1)
    values[0]=<double complex*>sig_malloc(sizeof(double complex*)*2)
    values[0][0]=0
    values[0][1]=1
    SMAT_cplx_dp(U,n,1,1,C,values,cset)
    #SMAT_cplx(V,ncols-num_set,comp_dim,num_set,C,vals_list,setc_list)
    cdef dict res={}
    #cdef double complex tmp
    for i in range(n):
        #tmp = CC(C[i])
        res[i]=C[0][i]
    if U<>NULL:
        sig_free(U)
    if C<>NULL:
        sig_free(C)
    return res



@cython.cdivision(True)
cdef SMAT_cplx_dp(double complex** U,int N,int num_rhs,int num_set,double complex **C,double complex** values,int* setc):
    r"""
    Use Gauss elimination to solve a linear system AX=B
    U = (A|B) is a N x (N+num_rhs) double complex matrix
    setc and values should be allocated of length num_set
    """
    cdef int m,maxi,j,k,i
    cdef double complex TT,temp2
    #cdef double complex **tmpu
    cdef double temp
    cdef int *piv
    cdef int *used
    piv=<int*>sig_malloc(sizeof(int)*N)
    used=<int*>sig_malloc(sizeof(int)*N)
    if C==NULL:
        C=<double complex**>sig_malloc(sizeof(double complex*)*num_rhs)
        for j in range(num_rhs):
            C[j]=<double complex*>sig_malloc(sizeof(double complex*)*N)
    for j in range(N):
        piv[j]=0
        used[j]=0
    for m in range(N):
        temp=0.0
        maxi=0
        #Locate maximum
        for j in range(N):
            if used[j]<>0:
                continue
            #print "U[",j,m,"]=",U[j][m],abs(U[j][m])
            if cabs(U[j][m]) <= temp:
                continue
            maxi=j
            temp=cabs(U[j][m])
            #print "temp=",temp
        piv[m]=maxi
        #print "piv[",m,"]=",maxi
        #print "norm=",temp
        #return
        used[maxi]=1
        temp2=U[maxi][m]
        if cabs(temp2)==0.0:
            print 'ERROR: pivot(',m,') == 0, system bad!!!'
            raise ArithmeticError
        for j in range(m+1,N+num_rhs): # do j=M+1,N+1
            U[maxi][j]=U[maxi][j]/temp2
            #! eliminate from all other rows
        for j in range(maxi):
            TT=U[j][m]
            for k in range(m+1,N+num_rhs): #K=M+1,N+1
                U[j][k]=U[j][k]-U[maxi][k]*TT
        for j in range(maxi+1,N): #DO J=Maxi+1,N
            TT=U[j][m]
            for k in range(m+1,N+num_rhs): #do K=M+1,N+1
                U[j][k]=U[j][k]-U[maxi][k]*TT
      #!! now remember we have pivot for x_j at PIV(j) not j
    cdef int do_cont,m_offs
    #print "N+num_rhs=",N+num_rhs
    for i in range(num_rhs):
        m_offs=0
        for m in range(N+num_set): #DO M=1,N
            do_cont=0
            for j in range(num_set):
                if setc[j]==m:
                    do_cont=1
                    break
            if do_cont==1:
                m_offs=m_offs+1
                C[i][m]=values[i][j]
                #continue
            else:
                C[i][m]=U[piv[m-m_offs]][N+i]
            #print "i,m,N+i=",i,m,N+i
            #print "C[",i,m,"]=",C[i][m]
            #print "next"
            #print "C[",i,m,"]=",C[i][m]
            #print "=U[",piv[m-m_offs],"][",N+i,"]"
        #print "C0[",m,"]=",C[m]
    if piv<>NULL:
        sig_free(piv)
    if used<>NULL:
        sig_free(used)
    # if tmpu<>NULL:
    #     for j in range(N):
    #         if tmpu[j]<>NULL:
    #             sig_free(tmpu[j])
    #     sig_free(tmpu)

@cython.cdivision(True)
cdef SMAT_real_dp(double ** U,int N,int num_rhs,int num_set,double **C,double ** values,int* setc):
    r"""
    Use Gauss elimination to solve a linear system AX=B
    U = (A|B) is a N x (N+num_rhs) double matrix
    setc and values should be allocated of length num_set
    """
    if C==NULL or values==NULL:
        raise ValueError,"Need an allocated vector C as input!"
    cdef int m,maxi,j,k
    cdef double TT,temp2
    #cdef double **tmpu
    cdef double temp
    cdef int *piv
    cdef int *used
    piv=<int*>sig_malloc(sizeof(int)*N)
    used=<int*>sig_malloc(sizeof(int)*N)
    for j in range(N):
        piv[j]=0
        used[j]=0
    for m in range(N): #DO m=1,N
        temp=0.0
        maxi=0
        #Locate maximum
        for j in range(N):
            if used[j]<>0:
                continue
            #print "U[",j,m,"]=",U[j][m],abs(U[j][m])
            if fabs(U[j][m]) <= temp:
                continue
            maxi=j
            temp=fabs(U[j][m])
            #print "temp=",temp
        piv[m]=maxi
        used[maxi]=1
        temp2=U[maxi][m]
        if temp2==0.0:
            print 'ERROR: pivot(',m,') == 0, system bad!!!'
            raise ArithmeticError
        for j in range(m+1,N+num_rhs): # do j=M+1,N+1
            U[maxi][j]=U[maxi][j]/temp2
            #! eliminate from all other rows
        for j in range(maxi-1):
            TT=U[j][m]
            for k in range(m+1,N+num_rhs): #K=M+1,N+1
               U[j][k]=U[j][k]-U[maxi][k]*TT
        for j in range(maxi+1,N): #DO J=Maxi+1,N
            TT=U[j][m]
            for k in range(m+1,N+num_rhs): #do K=M+1,N+1
               U[j][k]=U[j][k]-U[maxi][k]*TT
    #!! now remember we have pivot for x_j at PIV(j) not j
    cdef int do_cont,m_offs
    m_offs=0
    for i  in range(num_rhs):
        m_offs=0
        for m in range(N): #DO M=1,N
            do_cont=0
            for j in range(num_set):
                if setc[j]==m:
                    do_cont=1
                    break
            if do_cont==1:
                m_offs=m_offs+1
                C[i][m]=values[i][j]
                continue
            #print "i,m,N+i=",i,m,N+i
            #print "next"
            C[i][m]=U[piv[m-m_offs]][N+i]
            #print "C[",i,m,"]=",C[i][m]


    if piv<>NULL:
        sig_free(piv)
    if used<>NULL:
        sig_free(used)
    # if tmpu<>NULL:
    #     for j in range(N):
    #         if tmpu[j]<>NULL:
    #             sig_free(tmpu[j])
    #     sig_free(tmpu)












def mppr(x):
    r"""
    Print fewer digits!

    """

    # We don't want to print 500 digits
    # So I truncate at two digits with a stupid simple rounding
    import mpmath
    cdef int dpold=mpmath.mp.dps
    mpmath.mp.dps=5
    cdef str s,ss,s1,s2,es,dot,s3,s4,ma
    s=str(x)
    mpmath.mp.dps=dpold
    (s1,dot,s2)=s.partition(".")
    (s3,es,s4)=s2.partition("e")
    if(len(s4)==0):  # No "e" format
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

cpdef Ttest_mul(x,y,test=1):
    cdef double z=0.0
    cdef RealNumber w
    if test==1:
        _test_mul1(<double>x,<double>y,z)
    else:
        w = RealField(53)(0)
        _test_mul2(<RealNumber>x,<RealNumber>y,w)


cdef void _test_mul1(double x,double y,double z):
    z= x*y

cdef void _test_mul2(RealNumber x,RealNumber y,RealNumber z):
    z= x*y

cdef test_mul22(RealNumber x,RealNumber y,RealNumber z):
    z= x*y

cpdef split_interval(H,double R1,double R2):
    # First we find the next zero
    # First split into intervals having at most one zero
    ivs=list()
    cdef double Y00,k,rnew,rold,r1,r2,Y0,Y1,r11,t,pi
    pi=M_PI
    cdef int i,verbose
    verbose=H._verbose
    rnew=R1; rold=R1
    Y00=0.995*sqrt(3.0)/<double>(2 *H._group._level)
    while rnew < R2:
        rnew=min(R2,H.next_eigenvalue(rold))
        if abs(rold-rnew)==0.0:
            if verbose>0:
                print "ivs=",ivs
            break
        iv=(rold,rnew)
        ivs.append(iv)
        rold=rnew
        # We now need to split these intervals into pieces with at most one zero of the K-Bessel function
        new_ivs=list()
        for (r1,r2) in ivs:
            if verbose>0:
                print "r1,r2=",r1,r2
            Y0=Y00; r11=r1
            i=0
            while(r11 < r2 and i<1000):
                t=next_kbessel_zero(r11,r2,Y0*pi);i=i+1
                if verbose>0:
                    print "r11,r2,Y0*pi,t=",r11,r2,Y0*pi,t
                iv=(r11,t,Y0); new_ivs.append(iv)
                # must find Y0 s.t. |besselk(it,Y0)| is large enough
                Y1=Y0
                k=my_kbes(t,Y1)
                j=0
                while j<1000 and abs(k)<1e-3:
                    Y1=Y1*0.999;j=j+1
                    k=my_kbes(t,Y1)
                Y0=Y1
                r11=t+1E-08
        return new_ivs

def next_kbessel_zero(double r1,double r2,double y,int verbose=0):
    import mpmath
    base=mpmath.fp
    cdef double h,t1,t2,r0,kd0,kd,xtol,rtol
    cdef int i
    h=(r2-r1)/500.0
    t1=-1.0; t2=-1.0
    r0=r1
    my_kbes_diff_real_dp(r0,y,&kd0)
    xtol=1E-6
    rtol=1E-6

    if verbose>0:
        print "r0,y=",r0,y,kd0
    while(t1<r1 and r0<r2):
        # Let us first find a R-value for which the derivative changed sign
        my_kbes_diff_real_dp(r0,y,&kd)
        i=0
        while(kd*kd0>0 and i<500 and r0<r2):
            i=i+1
            r0=r0+h
            my_kbes_diff_real_dp(r0,y,&kd)
            #t1=base.findroot(lambda x :  besselk_dp(x,y),r0)
            t1=find_root(lambda x :  my_kbes(x,y),r0-h,r0+h,xtol=xtol,rtol=rtol)
            r0=r0+h
            if verbose>0:
                print "r0,y,t1=",r0,y,t1
        if(r0>=r2 or t1>=r2):
            t1=r2
        r0=r1
        #kd0=my_kbes_diff_r(r0,y,base)
        my_kbes_diff_real_dp(r0,y,&kd0)
        while(t2<r1 and r0<r2):
            my_kbes_diff_real_dp(r0,y,&kd)
            #kd=my_kbes_diff_r(r0,y,base)
            i=0
            while(kd*kd0>0 and i<500 and r0<r2):
                i=i+1
                r0=r0+h
                my_kbes_diff_real_dp(r0,y,&kd)
                #kd=my_kbes_diff_r(r0,y,base)
                #t2=base.findroot(lambda x :
                t2=find_root(lambda x :  my_kbes(x,2.0*y),r0-h,r0+h,xtol=xtol,rtol=rtol)
            #t2=base.findroot(lambda x :  base.besselk(base.mpc(0,x),base.mpf(2*y),verbose=True).real,r0)
            r0=r0+h
        if r0>=r2 or t2>=r2:
            t2=r2
            #print "zero(besselk,y1,y2)(",r1,r2,")=",t1,t2
        t=min(min(max(r1,t1),max(r1,t2)),r2)
        return t


cdef my_kbes_diff_real_dp(double r,double x,double *diff):
    cdef double h,f1,f2
    h=1e-8
    #print "r,x,h=",r,x,h
    f1 = my_kbes(r,x+h)
    f2 = my_kbes(r,x-h)
    diff[0]=0.5*(f2-f1)/h

### Algorithms for detecting eigenvalues by signchanges
cpdef get_coeff_and_signs_fast_real_dp(S,double R,double Y,int M,int Q,double Y2=0.0,dict Norm={},int gr=0,int num_tests=5):
    r"""
    An efficient method to get coefficients in the double real case.
    """
    cdef double **V=NULL

    cdef double* Xm=NULL
    cdef double*** Xpb=NULL
    cdef double*** Ypb=NULL
    cdef double *** Cvec=NULL
    cdef int nc,Ql,Qs,Qf,Ml,Ms,Mf,j,k,l,N,sym_type,i,n,r
    cdef double Qfak
    cdef int verbose=S._verbose
    ## Check that this "real" method is applicable here.
    if not S.multiplier().is_real():
        raise ValueError,"This method should only be called for real characters/multipliers!"
    sym_type=S.sym_type()
    if sym_type not in [0,1]:
        raise ValueError,"This method should only be called for symmetrized (even/odd) functions!"
    if Q<M:
        Q=M+20
    sym_type = S._sym_type
    Ms=0;  Mf=M; Qs=1; Qf=Q
    Qfak=<double>(Q)/<double>(2)
    Ml=Mf-Ms+1
    Ql=Qf-Qs+1
    nc = S.group().ncusps()
    N = nc*Ml
    if Norm=={}:
        Norm = S.set_norm(1)
    cdef list SetCs
    cdef dict Vals
    SetCs=Norm['SetCs'][0]
    Vals=Norm['Vals']
    cdef int comp_dim,num_set
    comp_dim=Norm['comp_dim']
    num_set=len(SetCs)
    V=<double**>sig_malloc(sizeof(double*)*N)
    if V==NULL: raise MemoryError
    cdef int ncols
    if num_set<comp_dim:
        ncols=N+comp_dim-num_set
    else:
        ncols=N
    if verbose>2:
        print "In get_coef_real_dp R=",R
    if verbose>2:
        print "N=",N
        print "Vals=",Vals
        print "SetCs=",SetCs
        print "comp_dim=",comp_dim
        print "num_set=",num_set
        print "ncols=",ncols
    for j in range(N):
        V[j]=<double*>sig_malloc(sizeof(double)*(ncols))
        for k in range(ncols):
            V[j][k]=<double>0
    Xm=<double*>sig_malloc(sizeof(double)*Ql)
    if Xm==NULL: raise MemoryError
    Xpb = <double***> sig_malloc( sizeof(double** ) * nc )
    if Xpb==NULL: raise MemoryError
    Ypb = <double***> sig_malloc( sizeof(double** ) * nc )
    if Ypb==NULL: raise MemoryError
    for i in range(nc):
        Xpb[i] = <double**>sig_malloc(sizeof(double*) * nc )
        Ypb[i] = <double**>sig_malloc(sizeof(double*) * nc )
        if Ypb[i]==NULL or Xpb[i]==NULL:
            raise MemoryError
        for j in range(nc):
            Xpb[i][j] = <double*>sig_malloc(sizeof(double) * Ql )
            Ypb[i][j] = <double*>sig_malloc(sizeof(double) * Ql )
            if Ypb[i][j]==NULL or Xpb[i][j]==NULL:
                raise MemoryError
            for n in range(Ql):
                Xpb[i][j][n]=<double>0
                Ypb[i][j][n]=<double>0
    Cvec = <double ***>sig_malloc(sizeof(double **) * nc )
    if Cvec==NULL: raise MemoryError
    for i in range(nc):
        Cvec[i] = <double **>sig_malloc(sizeof(double *) * nc )
        if Cvec[i]==NULL:
            raise MemoryError
        for j in range(nc):
            Cvec[i][j] = <double *>sig_malloc(sizeof(double) * Ql )
            if Cvec[i][j]==NULL:
                raise MemoryError
            for n in range(Ql):
                Cvec[i][j][n]=<double>0
    cdef int t=0
    t = pullback_pts_real_dp(S,Qs,Qf,Y,Xm,Xpb,Ypb,Cvec)
    if t==1:
        raise ArithmeticError,"Need smaller Y than {0}".format(Y)
    cdef double * alphas=NULL
    alphas=<double*>sig_malloc(sizeof(double)*nc)
    for j in range(nc):
        alphas[j]=float(S.alpha(j)[0])
    compute_V_real_dp(V,R,Y,Ms,Mf,Qs,Qf,nc,alphas,Xm,Xpb,Ypb,Cvec,1,sym_type,verbose)
    cdef cnp.ndarray[DTYPE_t,ndim=2] VV
    if gr==1:
        VV=np.zeros([nc*Ml, nc*Ml], dtype=DTYPE)
        #VV=mpmath.fp.matrix(int(nc*Ml),int(nc*Ml),force_type=mpmath.fp.mpc)
        for n in range(nc*Ml):
            for l in range(nc*Ml):
               VV[n,l]=V[n][l]
        return VV
    #cdef double *C=NULL
    #C=<double*>sig_malloc(sizeof(double)*(nc*Ml))
    #if C==NULL:
    #    raise MemoryError
    cdef double **C=NULL
    C=<double**>sig_malloc(sizeof(double*)*(comp_dim))
    if C==NULL:
        raise MemoryError
    for i in range(comp_dim):
        C[i]=<double*>sig_malloc(sizeof(double)*(nc*Ml))
        if C[i]==NULL:
            raise MemoryError


    cdef int *setc_list=NULL
    cdef double **vals_list=NULL
    setc_list = <int*>sig_malloc(sizeof(int)*num_set)
    vals_list = <double** >sig_malloc(sizeof(double*)*comp_dim)
    for j in range(comp_dim):
        vals_list[j]=<double*>sig_malloc(sizeof(double)*num_set)
    i=0
    for r,n in SetCs:
        if r*Ml+n-Ms<0:
            continue
        setc_list[i]=r*Ml+n-Ms
        for j in range(comp_dim):
            vals_list[j][i]=<double>CC(Vals[j][(r,n)])
        i=i+1
    cdef double *RHS=NULL
    cdef num_rhs=0
    if num_rhs>0 and num_rhs<>comp_dim:
        raise ValueError,"Need same number of right hand sides (or just one) as the number of set coefficients!"
    normalize_matrix_real_dp(V,N,comp_dim,num_set,setc_list,vals_list,verbose)
    if gr==2:
        VV=np.zeros([nc*Ml, nc*Ml], dtype=DTYPE)
        for n in range(N):
           for l in range(N): #+comp_dim:
                VV[n,l]=V[n][l]
        return VV
    SMAT_real_dp(V,ncols-num_set,comp_dim,num_set,C,vals_list,setc_list)
    ## We now compute part of V(Y2,R)*C
    if Y2==0:
        Y2=0.94*Y

    compute_V_real_dp(V,R,Y2,Ms,Mf,Qs,Qf,nc,alphas,Xm,Xpb,Ypb,Cvec,1,sym_type,verbose)
    #normalize_matrix_real_dp(V,N,comp_dim,num_set,setc_list,vals_list,verbose)
    signs={}


    cdef int Mmax=max(Ml,num_tests)
    cdef dict res={}
    cdef double tmp
    res[-1]=dict()
    for i in range(Mmax):
        # check if i in setc_list
        in_list=0
        for j in range(num_set):
            if setc_list[j]==i:
                in_list=1
                break
        if in_list==0:
            tmp=0
            #print "i:",i
            for j in range(N-num_set):
                #in_list=0
                #for k from 0<=k< num_set:
                #    if setc_list[j]==i:
                #        in_list=1
                #        break
                #if in_list==0:
                #if i==2:
                #    print "C[",j,"]=",C[j]
                #    print "V[",i,",",j,"]=",V[i][j]
                tmp=tmp+C[0][j]*V[i][j]
            res[-1][i]=tmp
    for k in range(comp_dim):
        res[k]=dict()
        for i in range(nc):
            res[k][i]=dict()
            for j in range(Ml):
                res[k][i][j+Ms]=C[k][i*nc+j]

    if V<>NULL:
        sig_free(V)
    if C<>NULL:
        sig_free(C)
    if alphas<>NULL:
        sig_free(alphas)
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
    return res


cpdef get_coeff_and_signs_fast_cplx_dp(S,double R,double Y,int M,int Q,double Y2=0,dict Norm={},int gr=0,int num_tests=5,int ncpus=1):
    r"""
    An efficient method to get coefficients in the double complex case.
    """
    import mpmath
    cdef double complex **V=NULL
    cdef double* Xm=NULL
    cdef double*** Xpb=NULL
    cdef double*** Ypb=NULL
    cdef double complex*** Cvec=NULL
    cdef int nc,Ql,Qs,Qf,Ml,Ms,Mf,j,k,l,N,sym_type,i,n,r
    #    cdef double Qfak
    cdef double complex tmp
    cdef int verbose=S._verbose
    cdef int is_exceptional = S._exceptional
    if Q<M:
        Q=M+20
    sym_type = S._sym_type
    if sym_type in [0,1]:
        Ms=0;  Mf=M; Qs=1; Qf=Q
        #Qfak=<double>(Q)/<double>(2)
    else:
        Ms=-M;  Mf=M; Qs=1-Q; Qf=Q
        #Qfak=<double>(2*Q)
    Ml=Mf-Ms+1
    Ql=Qf-Qs+1
    nc = S.group().ncusps()
    N = nc*Ml
    if Norm=={}:
        Norm = S.set_norm(1)
    cdef list SetCs
    cdef dict Vals
    SetCs=Norm['SetCs'][0]
    Vals=Norm['Vals']
    cdef int comp_dim,num_set
    comp_dim=Norm['comp_dim']
    num_set=len(SetCs)
    V=<double complex**>sig_malloc(sizeof(double complex*)*N)
    if V==NULL: raise MemoryError
    cdef int ncols
    if num_set<comp_dim:
        ncols=N+comp_dim-num_set
    else:
        ncols=N
    if verbose>2:
        print "In get_coef_cplx_dp R=",R
    if verbose>2:
        print "N=",N
        print "Vals=",Vals
        print "SetCs=",SetCs
        print "comp_dim=",comp_dim
        print "num_set=",num_set
        print "ncols=",ncols
    for j in range(N):
        V[j]=<double complex*>sig_malloc(sizeof(double complex)*(ncols))
        for k in range(ncols):
            V[j][k]=0
    Xm=<double*>sig_malloc(sizeof(double)*Ql)
    if Xm==NULL: raise MemoryError
    Xpb = <double***> sig_malloc( sizeof(double** ) * nc )
    if Xpb==NULL: raise MemoryError
    Ypb = <double***> sig_malloc( sizeof(double** ) * nc )
    if Ypb==NULL: raise MemoryError
    for i in range(nc):
        Xpb[i]=NULL; Ypb[i]=NULL
        Xpb[i] = <double**>sig_malloc(sizeof(double*) * nc )
        Ypb[i] = <double**>sig_malloc(sizeof(double*) * nc )
        if Ypb[i]==NULL or Xpb[i]==NULL:
            raise MemoryError
        for j in range(nc):
            Xpb[i][j]=NULL; Ypb[i][j]=NULL
            Xpb[i][j] = <double*>sig_malloc(sizeof(double) * Ql )
            Ypb[i][j] = <double*>sig_malloc(sizeof(double) * Ql )
            if Ypb[i][j]==NULL or Xpb[i][j]==NULL:
                raise MemoryError
            for n in range(Ql):
                Xpb[i][j][n]=<double>0
                Ypb[i][j][n]=<double>0
    Cvec = <double complex***>sig_malloc(sizeof(double complex**) * nc )
    if Cvec==NULL: raise MemoryError
    for i in range(nc):
        Cvec[i] = <double complex**>sig_malloc(sizeof(double complex*) * nc )
        if Cvec[i]==NULL:
            raise MemoryError
        for j in range(nc):
            Cvec[i][j] = <double complex*>sig_malloc(sizeof(double complex) * Ql )
            if Cvec[i][j]==NULL:
                raise MemoryError
            for n in range(Ql):
                Cvec[i][j][n]=0
    cdef int t=0
    t = pullback_pts_cplx_dp(S,Qs,Qf,Y,Xm,Xpb,Ypb,Cvec)
    if t==1:
        raise ArithmeticError,"Need smaller Y than {0}".format(Y)
    cdef double * alphas=NULL
    alphas=<double*>sig_malloc(sizeof(double)*nc)
    for j in range(nc):
        alphas[j]=<double>S.alpha(j)[0]
    cdef int cuspidal=1
    cdef int **Mv,**Qv
    Mv=<int**>sig_malloc(sizeof(int*)*nc)
    Qv=<int**>sig_malloc(sizeof(int*)*nc)
    for i in range(nc):
        Mv[i]=<int*>sig_malloc(sizeof(int)*3)
        Qv[i]=<int*>sig_malloc(sizeof(int)*3)        
        Mv[i][0]=Ms; Mv[i][1]=Mf; Mv[i][2]=Mf-Ms+1
        Qv[i][0]=Qs; Qv[i][1]=Qf; Qv[i][2]=Qf-Qs+1
    compute_V_cplx_dp(V,R,Y,Mv,Qv,nc,cuspidal,sym_type,verbose,alphas,Xm,Xpb,Ypb,Cvec,is_exceptional=is_exceptional)
    cdef cnp.ndarray[CTYPE_t,ndim=2] VV
    if gr==1:
        #VV=mpmath.fp.matrix(int(nc*Ml),int(nc*Ml),force_type=mpmath.fp.mpc)
        VV=np.zeros([nc*Ml, nc*Ml], dtype=CTYPE)
        for n in range(nc*Ml):
            for l in range(nc*Ml):
                VV[n,l]=V[n][l]
        return VV
    cdef double complex **C=NULL
    C=<double complex**>sig_malloc(sizeof(double complex*)*(comp_dim))
    for i in range(comp_dim):
        C[i]=<double complex*>sig_malloc(sizeof(double complex)*(nc*Ml))
    if C==NULL:
        raise MemoryError

    cdef int *setc_list=NULL
    cdef double complex **vals_list=NULL
    setc_list = <int*>sig_malloc(sizeof(int)*num_set)
    vals_list = <double complex** >sig_malloc(sizeof(double complex*)*comp_dim)
    for j in range(comp_dim):
        vals_list[j]=<double complex*>sig_malloc(sizeof(double complex)*num_set)
    i=0
    for r,n in SetCs:
        if r*Ml+n-Ms<0:
            continue
        setc_list[i]=r*Ml+n-Ms
        for j in range(comp_dim):
            vals_list[j][i]=<double complex>CC(Vals[j][(r,n)])
        i=i+1

    cdef double complex *RHS=NULL
    cdef num_rhs=0
    if num_rhs>0 and num_rhs<>comp_dim:
        raise ValueError,"Need same number of right hand sides (or just one) as the number of set coefficients!"
    normalize_matrix_cplx_dp(V,N,comp_dim,num_set,setc_list,vals_list,verbose)
    if gr==2:
        VV=np.zeros([nc*Ml, nc*Ml], dtype=CTYPE)
        #VV=mpmath.fp.matrix(int(N),int(N+comp_dim),force_type=mpmath.fp.mpc)
        for n in range(N):
            for l in range(N+comp_dim):
                VV[n,l]=V[n][l]
        return VV
    SMAT_cplx_dp(V,ncols-num_set,comp_dim,num_set,C,vals_list,setc_list)
    if verbose>2:
        for k in range(Ml):
            print "C[",k,"]=",C[0][k]
    cdef dict res={}
    if Y2==0:
        Y2=0.94*Y

#    compute_V_cplx_dp(V,R,Y,Mv,Qv,nc,cuspidal,sym_type,verbose,alphas,Xm,Xpb,Ypb,Cvec)
    cdef int Mmax=max(Ml,num_tests)
    res[-1]=dict()
    for i in range(Mmax):
        # check if i in setc_list
        in_list=0
        for j in range(num_set):
            if setc_list[j]==i:
                in_list=1
                break
        if in_list==0:
            tmp=0
            for j in range(N-num_set):
                if i==2:
                    print "C[",j,"]=",C[0][j]
                    print "V[",i,",",j,"]=",V[i][j]
                tmp=tmp+C[0][j]*V[i][j]
            res[-1][i]=tmp

    for j in range(comp_dim):
        res[j]=dict()
        for i in range(nc):
            res[j][i]=dict()
            for k in range(Ml):
                res[j][i][k+Ms]=C[j][i*Ml+k]


    if V<>NULL:
        sig_free(V)
    if C<>NULL:
        sig_free(C)
    if alphas<>NULL:
        sig_free(alphas)
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
    #if comp_dim==1:
    #    return res[0]
    return res

cdef void set_Mv_Qv_symm(S,int **Mv,int **Qv,double *Qfak,int * symmetric_cusps,double complex* cusp_evs,int *cusp_offsets,int *N,int *Ml,int *Ql,int M,int Q,int verbose=0):
    cdef int nc=S._group.ncusps()
    cdef int sym_type = S._sym_type
    N[0]=0
    Ml[0]=0
    Ql[0]=0
    for i in range(nc):
        symmetric_cusps[i]=0
        # print "sym_type=",sym_type
        symmetries=S.even_odd_symmetries()
        if sym_type<>-1 and symmetries[i][0]<>0:
            if symmetries[i][1]==1:
                symmetric_cusps[i]=int(sym_type)
            elif sym_type==0 and symmetries[i][1]==-1:
                symmetric_cusps[i]=1
            elif sym_type==1 and symmetries[i][1]==-1:
                symmetric_cusps[i]=0
            Mv[i][0]=0; Mv[i][1]=M;  Mv[i][2]=M+1
            Qv[i][0]=1; Qv[i][1]=Q;  Qv[i][2]=Q
            #Qv[i][0]=1-Q; Qv[i][1]=Q;  Qv[i][2]=2*Q
            N[0]=N[0]+M+1
        else:
            #print "S._group._symmetrizable_cusp[",i,"]=",S._group._symmetrizable_cusp[i]
            symmetric_cusps[i]=-1
            #sym_type=-1
            Mv[i][0]=-M; Mv[i][1]=M;  Mv[i][2]=2*M+1
            Qv[i][0]=1-Q; Qv[i][1]=Q;  Qv[i][2]=2*Q
            N[0]=N[0]+2*M+1
        if Mv[i][2]>Ml[0]:
            Ml[0]=Mv[i][2]
        if Qv[i][2]>Ql[0]:
            Ql[0]=Qv[i][2]
        if verbose>0:
            #print "symmetries=",symmetries
            print "Mv[",i,"]=",Mv[i][0],Mv[i][1],Mv[i][2]
            print "Qv[",i,"]=",Qv[i][0],Qv[i][1],Qv[i][2]
            print "sym_cusps_evs[",i,"]=",symmetric_cusps[i]
    for j in range(nc):
        if symmetric_cusps[j]<0:
            #    #if Qv[j][2]==2*Q:
            Qfak[j]=<double>(2*Q) # = 2Q # factor for 1-Q<=j<=Q
        else: #elif Qv[j][2]==Q:
            #Qfak[j]=<double>(Qv[j][1]) # = Q/2# factor for 1<=j<=Q
            Qfak[j]=<double>(Q)/2.0 # We use the same factor but different cos/sin/exp functions instead
        if verbose>0:
            print "Qfak[{0}]={1}".format(j,Qfak[j])
        #else:
        #    raise ArithmeticError,"Got unknown Qv!"
    for i in range(nc):
        cusp_offsets[i]=0
        for j in range(i):
            if j==0 or cusp_evs[j]==0:
                cusp_offsets[i]+=Mv[j][2]

    if verbose>0:
        print "N,Ml,Ql=",N[0],Ml[0],Ql[0]




# cdef void set_Mv_Qv_symm_real(S,int **Mv,int **Qv,double *Qfak,int * symmetric_cusps,double * cusp_evs,int *cusp_offsets,int *N,int *Ml,int *Ql,int M,int Q,int verbose=0):
#     r"""
#     To be used in the completely real case. 
#     """
#     cdef int nc=S._group.ncusps()
#     cdef int sym_type = S._sym_type
#     N[0]=0
#     Ml[0]=0
#     Ql[0]=0
#     for i in range(nc):
#         symmetric_cusps[i]=0
#         # print "sym_type=",sym_type
#         symmetries=S.even_odd_symmetries()
#         if sym_type<>-1 and symmetries[i][0]<>0:
#             if symmetries[i][1]==1:
#                 symmetric_cusps[i]=int(sym_type)
#             elif sym_type==0 and symmetries[i][1]==-1:
#                 symmetric_cusps[i]=1
#             elif sym_type==1 and symmetries[i][1]==-1:
#                 symmetric_cusps[i]=0
#             Mv[i][0]=0; Mv[i][1]=M;  Mv[i][2]=M+1
#             #Qv[i][0]=1; Qv[i][1]=Q;  Qv[i][2]=Q
#             Qv[i][0]=1-Q; Qv[i][1]=Q;  Qv[i][2]=2*Q
#             N[0]=N[0]+M+1
#         else:
#             #print "S._group._symmetrizable_cusp[",i,"]=",S._group._symmetrizable_cusp[i]
#             symmetric_cusps[i]=-1
#             #sym_type=-1
#             Mv[i][0]=-M; Mv[i][1]=M;  Mv[i][2]=2*M+1
#             Qv[i][0]=1-Q; Qv[i][1]=Q;  Qv[i][2]=2*Q
#             N[0]=N[0]+2*M+1
#         if Mv[i][2]>Ml[0]:
#             Ml[0]=Mv[i][2]
#         if Qv[i][2]>Ql[0]:
#             Ql[0]=Qv[i][2]
#         if verbose>0:
#             #print "symmetries=",symmetries
#             print "Mv[",i,"]=",Mv[i][0],Mv[i][1],Mv[i][2]
#             print "Qv[",i,"]=",Qv[i][0],Qv[i][1],Qv[i][2]
#             print "sym_cusps_evs[",i,"]=",symmetric_cusps[i]
#     for j in range(nc):
#         if symmetric_cusps[j]<0:
#             #if Qv[j][2]==2*Q:
#             Qfak[j]=<double>(Qv[j][2]) # = 2Q # factor for 1-Q<=j<=Q
#         else: #elif Qv[j][2]==Q:
#             Qfak[j]=<double>(Qv[j][1]) # = Q/2# factor for 1<=j<=Q
#         if verbose>0:
#             print "Qfak[{0}]={1}".format(j,Qfak[j])
#         #else:
#         #    raise ArithmeticError,"Got unknown Qv!"
#     for i in range(nc):
#         cusp_offsets[i]=0
#         for j in range(i):
#             if j==0 or cusp_evs[j]==0:
#                 cusp_offsets[i]+=Mv[j][2]

#     if verbose>0:
#         print "N,Ml,Ql=",N[0],Ml[0],Ql[0]


        
# ########################## TEST
# cpdef get_coeff_fast_cplx_dp_sym2(S,double R,double Y,int M,int Q,dict Norm={},int gr=0,int norm_c=1,dict cusp_ev={},double eps=1e-12,int do_par=0,int ncpus=1):
#     r"""

#     An efficient method to get coefficients in the double complex case.
#     Trying to use as much symmetries as possible.

#     """
#     import mpmath
#     cdef double complex **V=NULL
#     cdef double complex **V1=NULL
#     cdef double* Xm=NULL
#     cdef double*** Xpb=NULL
#     cdef double*** Ypb=NULL
#     cdef double complex*** Cvec=NULL
#     cdef double complex **C=NULL
#     cdef int *setc_list=NULL
#     cdef double complex *RHS=NULL
#     cdef double complex **vals_list=NULL
#     cdef double *alphas=NULL
#     cdef int nc,Ql,Qs,Qf,Ml,Ms,Mf,j,k,l,N,sym_type,i,n,r
#     cdef int num_rhs=0
#     cdef double *Qfak,tmpr
#     cdef double complex *cusp_evs=NULL
#     cdef int **Mv=NULL
#     cdef int **Qv=NULL
#     cdef int verbose=S._verbose
#     cdef int ncolsV1
#     cdef int *symmetric_cusps=NULL
#     cdef int N1=0
#     cdef list SetCs
#     cdef dict Vals
#     cdef int comp_dim,num_set
#     cdef int ncols,ncols1
#     cdef int cuspidal=1
#     cdef int q
#     cdef double complex *sqch=NULL
#     tmpr = <double>S._group.minimal_height()
#     if Y<= 0 or Y >= tmpr:
#         Y = 0.5 * tmpr
#     if M<=0:
#         M = get_M_for_maass_dp(R,Y,eps)

#     if Q<M:
#         Q=M+20
#     #    Qfak=<double>(2*Q)
#     sym_type = S._sym_type
#     nc = int(S._group._ncusps)
#     Qfak = <double *>sig_malloc(sizeof(double)*nc)
#     Mv=<int**>sig_malloc(sizeof(int*)*nc)
#     if not Mv: raise MemoryError
#     Qv=<int**>sig_malloc(sizeof(int*)*nc)
#     if not Qv: raise MemoryError
#     for i in range(nc):
#         Mv[i]=<int*>sig_malloc(sizeof(int)*3)
#         if not Mv[i]: raise MemoryError
#         Qv[i]=<int*>sig_malloc(sizeof(int)*3)
#         if not Qv[i]: raise MemoryError
#     symmetric_cusps=<int*> sig_malloc(sizeof(int)*nc)
#     if not symmetric_cusps: raise MemoryError
#     cusp_evs=<double complex*>sig_malloc(sizeof(double complex)*nc)
#     if not cusp_evs: raise MemoryError
#     N = 0; Ml=0; Ql=0
#     cdef int* cusp_offsets=NULL
#     cusp_offsets=<int*>sig_malloc(sizeof(int)*nc)
#     if cusp_offsets==NULL: raise MemoryError
#     cdef dict symmetries
#     if verbose>0:
#         print "INPUT: R={0}, Y={1}, M={2}, Q={3}, Norm={4}, gr={5}, norm_c={6}, cusp_ev={7}, eps={8}".format(<double>R,<double>Y,M,Q, Norm,gr,norm_c,cusp_ev,eps)
#     set_Mv_Qv_symm(S,Mv,Qv,Qfak,symmetric_cusps,cusp_evs,cusp_offsets,&N,&Ml,&Ql,M,Q,verbose)
#     Qs = 1-Q; Qf = Q
#     if verbose>0:
#         print "N=",N," Ml=",Ml," Ql=",Ql
#     if cusp_ev=={}:
#         cusp_ev = S.atkin_lehner_eigenvalues()
#     for i in range(nc):
#         Qfak[i]=<double>(2*Q)
#         cusp_evs[i]=CC(cusp_ev.get(S.group().cusps()[i],0))
#         if i==0 or cusp_evs[i]==0:
#             N1=N1+Mv[i][2]
#         if verbose>0:
#             print "cusp_ev[",i,"]=",cusp_evs[i]
#     for jcusp in range(nc):
#         cusp_offsets[jcusp]=0
#         for icusp in range(jcusp):
#             if icusp==0 or cusp_evs[icusp]==0:
#                 cusp_offsets[jcusp]+=Mv[icusp][2]
#                 if verbose>0:
#                     print "cusp_offset[",jcusp,"]+=",Mv[icusp][2]
#             #else:
#             #    print "cusp_ev=",cusp_evs[icusp]
#         if verbose>0:
#             print "cusp_offset[",jcusp,"]=",cusp_offsets[jcusp]

#     if verbose>0:
#         print "N1=",N1
#     if Norm=={}:
#         Norm = S.set_norm(1)
#     SetCs=Norm['SetCs'][0]
#     Vals=Norm['Vals']
#     comp_dim=Norm['comp_dim']
#     num_set=0
#     for r,n in SetCs:
#         if r>0 and cusp_evs[r]<>0:
#             continue
#         num_set+=1
#     V=<double complex**>sig_malloc(sizeof(double complex*)*N)
#     if V==NULL: raise MemoryError
#     V1=<double complex**>sig_malloc(sizeof(double complex*)*N1)
#     if V1==NULL: raise MemoryError
#     if num_set<comp_dim:
#         ncols=N+comp_dim-num_set
#         ncols1=N1+comp_dim-num_set
#     else:
#         ncols=N
#         ncols1=N1
#     if verbose>0:
#         print "In get_coef_cplx_dp_sym R,Y,M,Q=",R,Y,M,Q
#     if verbose>0:
#         print "N=",N
#         print "ncols=",ncols
#         print "N1=",N1
#         print "ncols1=",ncols1
#         print "Vals=",Vals
#         print "SetCs=",SetCs
#         print "comp_dim=",comp_dim
#         print "num_set=",num_set
#         for i in range(nc):
#             print "Ms,Mf,Ml[",i,"]=",Mv[i][0],Mv[i][1],Mv[i][2]
#             print "Qs,Qf,Ql[",i,"]=",Qv[i][0],Qv[i][1],Qv[i][2]
#     for j in range(N):
#         V[j]=<double complex*>sig_malloc(sizeof(double complex)*(ncols))
#         for k in range(ncols):
#             V[j][k]=0
#     for j in range(N1):
#         V1[j]=<double complex*>sig_malloc(sizeof(double complex)*(ncols1))
#         for k in range(ncols1):
#             V1[j][k]=0
#     alphas=<double*>sig_malloc(sizeof(double)*nc)
#     if alphas==NULL: raise MemoryError
#     #printf("alphas_init=%p \n",alphas)
#     Xm=<double*>sig_malloc(sizeof(double)*Ql)
#     if Xm==NULL: raise MemoryError
#     Xpb = <double***> sig_malloc( sizeof(double** ) * nc )
#     if Xpb==NULL: raise MemoryError
#     Ypb = <double***> sig_malloc( sizeof(double** ) * nc )
#     if Ypb==NULL: raise MemoryError
#     for i in range(nc):
#         Xpb[i]=NULL; Ypb[i]=NULL
#         Xpb[i] = <double**>sig_malloc(sizeof(double*) * nc )
#         Ypb[i] = <double**>sig_malloc(sizeof(double*) * nc )
#         if Ypb[i]==NULL or Xpb[i]==NULL:
#             raise MemoryError
#         for j in range(nc):
#             Xpb[i][j]=NULL; Ypb[i][j]=NULL
#             Xpb[i][j] = <double*>sig_malloc(sizeof(double) * Ql )
#             Ypb[i][j] = <double*>sig_malloc(sizeof(double) * Ql )
#             if Ypb[i][j]==NULL or Xpb[i][j]==NULL:
#                 raise MemoryError
#             for n in range(Ql):
#                 Xpb[i][j][n]=<double>0
#                 Ypb[i][j][n]=<double>0
#     Cvec = <double complex***>sig_malloc(sizeof(double complex**) * nc )
#     if Cvec==NULL: raise MemoryError
#     for i in range(nc):
#         Cvec[i] = <double complex**>sig_malloc(sizeof(double complex*) * nc )
#         if Cvec[i]==NULL:
#             raise MemoryError
#         for j in range(nc):
#             Cvec[i][j] = <double complex*>sig_malloc(sizeof(double complex) * Ql )
#             if Cvec[i][j]==NULL:
#                 raise MemoryError
#             for n in range(Ql):
#                 Cvec[i][j][n]=0
# #    sig_on()
#     cdef int t=0
#     t = pullback_pts_cplx_dp(S,Qs,Qf,Y,Xm,Xpb,Ypb,Cvec)
#     if t==1:
#         raise ArithmeticError,"Need smaller Y than {0}".format(Y)
# #    sig_off()
#     for j in range(nc):
#         tmpr=float(S.alpha(j)[0])
#         alphas[j]=<double>tmpr
#         if verbose>0:
#             print "alphas[",j,"]=",alphas[j],type(alphas[j])
# #    sig_on()
#     if do_par==1 and ncpus>=1:
#         compute_V_cplx_dp_sym_par(V1,N1,Xm,Xpb,Ypb,Cvec,
#                               cusp_evs,alphas,Mv,Qv,Qfak,
#                               symmetric_cusps,
#                               R,Y,nc,ncols,cuspidal,verbose,ncpus)
#     else:
#         compute_V_cplx_dp_sym2(V1,N1,Xm,Xpb,Ypb,Cvec,
#                               cusp_evs,alphas,Mv,Qv,Qfak,
#                               symmetric_cusps,
#                               R,Y,nc,ncols,cuspidal,verbose)
# #    sig_off()
#     cdef Matrix_complex_dense VV
#     #Vtmp = load("A.sobj")
#     #for l from 0<=l<N1:
#     #    for n from 0<=n<N1:
#     #        V1[l][n]=<double complex>complex(Vtmp[l][n])
#     if Qfak<>NULL:
#         sig_free(Qfak)
#     if alphas<>NULL:
#         sig_free(alphas)
#     if Xm<>NULL:
#         sig_free(Xm)
#     if Ypb<>NULL:
#         for i in range(nc):
#             if Ypb[i]<>NULL:
#                 for j in range(nc):
#                     if Ypb[i][j]<>NULL:
#                         sig_free(Ypb[i][j])
#                 sig_free(Ypb[i])
#         sig_free(Ypb)
#     if Xpb<>NULL:
#         for i in range(nc):
#             if Xpb[i]<>NULL:
#                 for j in range(nc):
#                     if Xpb[i][j]<>NULL:
#                         sig_free(Xpb[i][j])
#                 sig_free(Xpb[i])
#         sig_free(Xpb)
#     if Cvec<>NULL:
#         for i in range(nc):
#             if Cvec[i]<>NULL:
#                 for j in range(nc):
#                     if Cvec[i][j]<>NULL:
#                         sig_free(Cvec[i][j])
#                 sig_free(Cvec[i])
#         sig_free(Cvec)
#     if sqch<>NULL:
#         sig_free(sqch)
#     if gr==1:
#         CF=MPComplexField(53)
#         MS=MatrixSpace(CF,N1,ncols1)
#         VV=Matrix_complex_dense(MS,0)
#         for n in range(N1):
#             for l in range(ncols1):
#                 VV[n,l]=V1[n][l]
#         return VV
#     C=<double complex**>sig_malloc(sizeof(double complex*)*(comp_dim))
#     if C==NULL:
#         raise MemoryError
#     for i in range(comp_dim):
#         C[i]=<double complex*>sig_malloc(sizeof(double complex)*(N1))
#         if C[i]==NULL:
#             raise MemoryError
#     setc_list = <int*>sig_malloc(sizeof(int)*num_set)
#     if setc_list==NULL:
#         raise MemoryError
#     vals_list = <double complex** >sig_malloc(sizeof(double complex*)*comp_dim)
#     if vals_list==NULL:
#         raise MemoryError
#     for j in range(comp_dim):
#         vals_list[j]=<double complex*>sig_malloc(sizeof(double complex)*num_set)
#     i=0
#     for r,n in SetCs:
#         if cusp_evs[r]<>0 and r>0:
#             continue
#         l = cusp_offsets[r]+n-Mv[r][0]
#         if l<0:
#             continue
#         setc_list[i]=l
#         for j in range(comp_dim):
#             vals_list[j][i]=<double complex>CC(Vals[j][(r,n)])
#         i=i+1
#     if num_rhs>0 and num_rhs<>comp_dim:
#         raise ValueError,"Need same number of right hand sides (or just one) as the number of set coefficients!"
#     if verbose>0:
#         print "comp_dim=",comp_dim
#     normalize_matrix_cplx_dp(V1,N1,comp_dim,num_set,setc_list,vals_list,verbose)
#     if gr==2:
#         CF=MPComplexField(53)
#         MS=MatrixSpace(CF,N1,ncols1)
#         VV=Matrix_complex_dense(MS,0)
#         for n in range(N1):
#             for l in range(ncols1):
#                 VV[n,l]=V1[n][l]
#         return VV
#     sig_on()
#     if verbose>0:
#         print "comp_dim=",comp_dim
#     if do_par<=1:
#         SMAT_cplx_dp(V1,ncols1-num_set,comp_dim,num_set,C,vals_list,setc_list)
#     else:
#         if verbose>0:
#             print "smat parallel!"
#         SMAT_cplx_par_dp(V1,ncols1-num_set,comp_dim,num_set,C,vals_list,setc_list,ncpus)
#     sig_off()
#     if verbose>1:
#         for k in range(ncols1-num_set):
#             print "C[",k,"]=",C[0][k]
#     cdef dict res={}
#     cdef int ki
#     for j in range(comp_dim):
#         res[j]=dict()
#         for i in range(nc):
#             if i>0 and cusp_evs[i]<>0:
#                 res[j][i]=cusp_evs[i]
#                 continue
#             res[j][i]=dict()
#             for k in range(Mv[i][2]):
#                 ki=cusp_offsets[i]+k
#                 res[j][i][k+Mv[i][0]]=C[j][ki]

#     if setc_list<>NULL:
#         sig_free(setc_list)
#     if vals_list<>NULL:
#         for j in range(comp_dim):
#             if vals_list[j]<>NULL:
#                 sig_free(vals_list[j])
#         sig_free(vals_list)
#     if V<>NULL:
#         for j in range(N):
#             if V[j]<>NULL:
#                 sig_free(V[j])
#         sig_free(V)
#     if V1<>NULL:
#         for j in range(N1):
#             if V1[j]<>NULL:
#                 sig_free(V1[j])
#         sig_free(V1)

#     if C<>NULL:
#         for i in range(comp_dim):
#             if C[i]<>NULL:
#                 sig_free(C[i])
#         sig_free(C)
#     if Mv<>NULL:
#         sig_free(Mv)
#     if Qv<>NULL:
#         sig_free(Qv)
#     if symmetric_cusps<>NULL:
#         sig_free(symmetric_cusps)
#     if cusp_evs<>NULL:
#         sig_free(cusp_evs)
#     #if comp_dim==1:
#     #    return res[0]
#     return res
        

# cdef int compute_V_cplx_dp_sym2(double complex **V,
#                            int N1,
#                            double *Xm,
#                            double ***Xpb,
#                            double ***Ypb,
#                            double complex ***Cvec,
#                            double complex *cusp_evs,
#                            double *alphas,
#                            int **Mv,int **Qv,double *Qfak,
#                            int *symmetric_cusps,
#                            double R,double Y,
#                            int nc, int ncols,
#                            int cuspidal,
#                            int verbose):


#     r"""
#     Set up the matrix for the system of equations giving the Fourier coefficients of the Maass waveforms.
#     INPUT:

#     - ``R``   -- double (eigenvalue)
#     - ``Y``   -- double (the height of the sampling horocycle)
#     - ``Ms,Mf``  -- int (The coefficients we want to compute are C(n), Ms<=n<=Mf )
#     - ``Qs,Qf``  -- int (The sampling points are X_m, Qs<=m<=Qf)
#     - ``alphas`` -- [nc] double array (the shifts of the Fourier expansion at each cusp)
#     - ``V``   -- [(Mf-Ms)*nc]^2 double complex matrix (allocated)
#     - ``Xm``  -- [Qf-Qs] double array (allocated)
#     - ``Xpb`` -- nc*nc*[Qf-Qs] double array (allocated)
#     - ``Ypb`` -- nc*nc*[Qf-Qs] double array (allocated)
#     - ``Cvec`` -- nc*nc*[Qf-Qs] double complex array (allocated)
#     - `` cuspidal`` -- int (set to 1 if we compute cuspidal functions, otherwise zero)
#     - `` sym_type`` -- int (set to 0/1 if we compute even/odd functions, otherwise -1)
#     - ``verbose`` -- int (verbosity of output)

#     """
#     cdef int l,j,icusp,jcusp,n,ni,lj,Ml,Ql,s,Qs,Qf,Mf,Ms
#     cdef double pi,sqrtY,Y2pi,nrY2pi,argm,argpb,twopi,two,kbes,besarg,lr,nr
#     cdef double complex ckbes,ctmpV,iargm,twopii,ctmp
#     if not cuspidal in [0,1]:
#         raise ValueError," parameter cuspidal must be 0 or 1"
#     if R < 0:
#         ## In this case (corresponding to lambda in [0,1/4] we use the real parameter K-Bessel
#         R = -R
#         set_pref = -1
#     else:
#         set_pref = 1
#     pi=M_PI 
#     sqrtY=sqrt(Y)
#     two=<double>(2)
#     Y2pi=Y*pi*two
#     twopi=two*pi
#     Ml=0; Ql=0
#     for i in range(nc):
#         if Mv[i][2]>Ml:
#             Ml=Mv[i][2]
#         if Qv[i][2]>Ql:
#             Ql=Qv[i][2]
#         if verbose>0:
#             print "Qv[",i,"]=(",Qv[i][0],",",Qv[i][1],",",Qv[i][2],")"
#             print "Mv[",i,"]=(",Mv[i][0],",",Mv[i][1],",",Mv[i][2],")"
#             print "Qfak[{0}]={1}".format(i,Qfak[i])
#             print "symmetric_cusps[{0}]={1}".format(i,symmetric_cusps[i])
#             print "in sym2"
#     if verbose>0:
        
#         print "N1=",N1
#         print "Ql=",Ql
#     ## This is the effective offset at the
#     cdef int* cusp_offsets=NULL
#     cusp_offsets=<int*>sig_malloc(sizeof(int)*nc)
#     if cusp_offsets==NULL: raise MemoryError
#     for jcusp in range(nc):
#         cusp_offsets[jcusp]=0
#         for icusp in range(jcusp):
#             if icusp==0 or cusp_evs[icusp]==0:
#                 cusp_offsets[jcusp]+=Mv[icusp][2]
#         if verbose>0:
#             print "cusp_offsets[",jcusp,"]=",cusp_offsets[jcusp]
#     cdef int nc_sym=0
#     for jcusp in range(nc):
#         if verbose>0:
#             print "cusp_evs[",jcusp,"]=",cusp_evs[jcusp]
#         if jcusp==0 or cusp_evs[jcusp]<>0:
#             nc_sym+=1
#     cdef double **nvec=NULL
#     nvec = <double**>sig_malloc(sizeof(double*)*nc)
#     if not nvec: raise MemoryError
#     for icusp in range(nc):
#         nvec[icusp] = <double*>sig_malloc(sizeof(double)*Ml)
#     cdef double complex ****ef1=NULL
#     cdef double complex ***ef2_c=NULL
#     ef2_c = <double complex***>sig_malloc(sizeof(double complex**)*nc)
#     if not ef2_c: raise MemoryError
#     for icusp in range(nc):
#         ef2_c[icusp] = <double complex**>sig_malloc(sizeof(double complex*)*Mv[icusp][2])
#         for n in range(Mv[icusp][2]):
#             ef2_c[icusp][n] = <double complex*>sig_malloc(sizeof(double complex)*Qv[icusp][2])
#     ef1 = <double complex****>sig_malloc(sizeof(double complex***)*nc)
#     if ef1==NULL: raise MemoryError
#     for icusp in range(nc):
#         ef1[icusp] = <double complex***>sig_malloc(sizeof(double complex**)*nc)
#         if ef1[icusp]==NULL: raise MemoryError
#         for jcusp in range(nc):
#             ef1[icusp][jcusp] = <double complex**>sig_malloc(sizeof(double complex*)*Mv[jcusp][2])
#             if ef1[icusp][jcusp]==NULL: raise MemoryError
#             for n in range(Mv[jcusp][2]):
#                 ef1[icusp][jcusp][n] = <double complex*>sig_malloc(sizeof(double complex)*Qv[jcusp][2])
#                 if ef1[icusp][jcusp][n]==NULL: raise MemoryError
#     for jcusp in range(nc):
#         for n in range(Ml):
#             nvec[jcusp][n]=<double>(n+Mv[jcusp][0])+alphas[jcusp]
#     cdef int twoQm1
#     twoQm1= 2*Qv[0][1]-1
#     for jcusp in range(nc):
#         for n in range(Mv[jcusp][2]):
#             for j in range(Qv[jcusp][2]):
#                 argm=nvec[jcusp][n]*Xm[j]
#                 if symmetric_cusps[jcusp]==0:
#                     ef2_c[jcusp][n][j]=2.0*cos(argm)
#                 elif symmetric_cusps[jcusp]==1:
#                     ef2_c[jcusp][n][j]=_I2*sin(-argm)
#                 else:
#                     ef2_c[jcusp][n][j]=cexpi(-argm)

#     cdef double argpb1
#     for jcusp in range(nc):
#         for icusp in range(nc):
#             for n in range(Mv[jcusp][2]):
#                 for j in range(Qv[jcusp][2]): #in range(Qs,Qf+1):
#                     if Ypb[icusp][jcusp][j]==0: #not Xpb.has_key((icusp,jcusp,j):
#                         continue
#                     argpb=nvec[jcusp][n]*Xpb[icusp][jcusp][j]
#                     if symmetric_cusps[jcusp]==0:
#                         ef1[icusp][jcusp][n][j]=2.0*cos(argpb)
#                     elif symmetric_cusps[jcusp]==1:
#                         ef1[icusp][jcusp][n][j]=_I2*sin(argpb)
#                     else:
#                         ef1[icusp][jcusp][n][j]=cexpi(argpb)
#                     ctmp = Cvec[icusp][jcusp][j]
#                     ef1[icusp][jcusp][n][j]=ef1[icusp][jcusp][n][j]*ctmp

#     if verbose>1:
#         print "here1121"
#     cdef double besarg_old=0.0
#     cdef double y,kbes_old=1.0
#     cdef double ***kbesvec=NULL
#     kbesvec=<double***>sig_malloc(sizeof(double**)*nc)
#     if kbesvec==NULL:
#         raise MemoryError
#     for jcusp in range(nc):
#         kbesvec[jcusp]=<double**>sig_malloc(sizeof(double*)*Ml)
#         if kbesvec[jcusp]==NULL:
#             raise MemoryError
#         for l in range(Ml):
#             kbesvec[jcusp][l]=<double*>sig_malloc(sizeof(double)*Ql) #Qv[jcusp][2])
#             if kbesvec[jcusp][l]==NULL:
#                 raise MemoryError
#     if verbose>0:
#         print "here0"
#         print "Ml=",Ml
#         print "Ql=",Ql
#     cdef double tmpr
#     cdef double besprec
#     besprec=1.0E-14
#     ## Can we somehow use that "most" of Xpb[icusp][jcusp][2Q-j-1]=-Xpb[icusp][jcusp][j] ?
#     ## uncomment the below lines to see this.
#     # for jcusp in range(nc):
#     #     for icusp in range(nc):
#     #         for j from 0<=j<Qv[1][1]:
#     #             print "diff[",jcusp,icusp,j,"]=",Xpb[icusp][jcusp][j]+Xpb[icusp][jcusp][2*Qv[jcusp][1]-j-1]

#     for jcusp in range(nc):
#         for icusp in range(nc):
#             if icusp>0 and cusp_evs[icusp]<>0:
#                 continue
#             for j in range(Qv[jcusp][2]):
#                 if Ypb[icusp][jcusp][j]==0:
#                     continue
#                 for l in range(Mv[jcusp][2]):
#                     lr=nvec[jcusp][l]*twopi
#                     Mf = Mv[icusp][1]
#                     besarg=fabs(lr)*Ypb[icusp][jcusp][j]
#                     if lr<>0.0:
#                         besselk_dp_c(&tmpr,R,besarg,besprec,pref=set_pref)
#                         kbesvec[icusp][l][j]=sqrt(Ypb[icusp][jcusp][j])*tmpr
#                     else:
#                         kbesvec[icusp][l][j]=<double>1.0

#     cdef double complex cuspev

#     for jcusp in range(nc):
#         for l in range(Mv[jcusp][2]):
#             lr=nvec[jcusp][l]*twopi
#             lj=cusp_offsets[jcusp]+l
#             if jcusp>0 and cusp_evs[jcusp]<>0:
#                 lj0=l; cuspev=cusp_evs[jcusp]
#             else:
#                 lj0=lj;  cuspev=1.0
#             if lr==0.0 and cuspidal==1:
#                 continue
#             for j in range(Qv[jcusp][2]):
#                 for icusp in range(nc):
#                     if icusp>0 and cusp_evs[icusp]<>0:
#                         continue
#                     if Ypb[icusp][jcusp][j]==0:
#                         continue
#                     ckbes=kbesvec[icusp][l][j]*ef1[icusp][jcusp][l][j]
#                     for n in range(Mv[icusp][2]):
#                         if nvec[icusp][n]==0 and cuspidal==1:
#                             continue
#                         #if n+Mv[icusp][0]==0 and cuspidal==1:
#                         #    continue
#                         ni=cusp_offsets[icusp]+n
#                         ctmpV=ckbes*ef2_c[icusp][n][j]
#                         if ni>N1 or lj0>N1:
#                             raise ArithmeticError,"Index outside!"
#                         V[ni][lj0]=V[ni][lj0]+ctmpV*cuspev
#                         # if verbose>0 and lj0==324 and ni==46:
#                         if verbose>2 and ni==3 and lj0==3:
#                             print "l,jcusp,n,icusp=",l+Mv[jcusp][0],jcusp,n+Mv[icusp][0],icusp
#                             print "Ypb(",j+Qv[jcusp][0],")=",Ypb[icusp][jcusp][j]
#                             print "Cvec=",Cvec[icusp][jcusp][j]
#                             print "V[",ni,"][",lj0,"](",j,")=",V[ni][lj0]
#                             print "ctmpV[",ni,"][",lj0,"]=",ctmpV
#                             print "kbesvec[",icusp,l,j,"]=",kbesvec[icusp][l][j]
#                             print "besarg=",fabs(lr)*Ypb[icusp][jcusp][j]
#                             print "ef1=",ef1[icusp][jcusp][l][j]
#                             print "ef2=",ef2_c[icusp][n][j]
#                             print "V[",ni,"][",lj0,"]=",V[ni][lj0],"+",ctmpV,"*",cuspev

#     if verbose>2:
#     #    #for n from 30+cusp_offsets[3]<=n<40+cusp_offsets[3]:
#     #    for n from 43<=n<49:
#     #        print "V0[",0,n,"]=",V[0][n]
#     #        #print "V0[",24,24,"]=",V[24][24]
#         print "V0[",3,3,"]=",V[3][3]
#         for icusp in range(nc):
#             print "Qfak[{0}]={1}".format(icusp,Qfak[icusp])
#     for icusp in range(nc):
#         if icusp>0 and cusp_evs[icusp]<>0:
#             continue
#         for n in range(Mv[icusp][2]):
#             ni=cusp_offsets[icusp]+n
#             for jcusp in range(nc):
#                 if jcusp>0 and cusp_evs[jcusp]<>0:
#                     continue
#                 for l in range(Mv[jcusp][2]):
#                     lj=cusp_offsets[jcusp]+l
#                     if lj>N1 or ni>N1: # some extra insurance...
#                         print "ni=",cusp_offsets[icusp],"+",n,"=",ni
#                         print "lj=",cusp_offsets[jcusp],"+",l,"=",lj
#                         raise ArithmeticError,"Index outside!"
#                     if verbose>2 and ni==3 and lj==3:
#                         print "V[3][3]={0}/{1}".format(V[ni][lj],Qfak[jcusp])
#                     V[ni][lj]=V[ni][lj]/<double complex>Qfak[jcusp]
#                     if verbose>2 and ni==3 and lj==3:
#                         print "V[3][3]/Qfak[{0}]={1}".format(jcusp,V[ni][lj])

#     for icusp in range(nc):
#         if icusp>0 and cusp_evs[icusp]<>0:
#             continue
#         for n in range(Mv[icusp][2]):
#             nr=fabs(nvec[icusp][n])
#             if nr==0.0 and cuspidal==1:
#                 continue
#             ni=cusp_offsets[icusp]+n
#             nrY2pi=nr*Y2pi
#             if nr==0.0:
#                 kbes=<double>1.0
#             else:
#                 #mpIR=mpmath.fp.mpc(0,R)
#                 #                kbes=float(mpmath.fp.besselk(mpIR,nrY2pi).real*exp(mpmath.fp.pi*R*0.5))
#                 besselk_dp_c(&kbes,R,nrY2pi,besprec,pref=set_pref)
#                 kbes=sqrtY*kbes # besselk_dp(R,nrY2pi,pref=1)
#             if ni>N1:
#                 raise ArithmeticError,"Index outside!"
#             if verbose>2 and ni==3:
#                 print "V[3][3]={0}-{1}={2}".format(V[ni][ni],kbes,V[ni][ni]-kbes)
#             V[ni][ni]=V[ni][ni] - kbes

#     if verbose>3:

#         #        #for n from 30+cusp_offsets[3]<=n<40+cusp_offsets[3]:
#         #        print "V3[",0,n,"]=",V[0][n]
#         for n in range(min(100,Ml)):
#             print "V3[",n,n,"]=",V[n][n]
#     #if Qfak<>NULL:
#     #    #for icusp in range(nc):
#     #    #    if Qfak[icusp]<>NULL:
#     #    #        sig_free(Qfak[icusp])
#     #    sig_free(Qfak)
#     if kbesvec<>NULL:
#         for icusp in range(nc):
#             if kbesvec[icusp]<>NULL:
#                 for l in range(Ml):
#                     if kbesvec[icusp][l]<>NULL:
#                         sig_free(kbesvec[icusp][l])
#                 sig_free(kbesvec[icusp])
#         sig_free(kbesvec)
#     #print "deal kbbes1"
#     if ef1<>NULL:
#         for jcusp in range(nc):
#             if ef1[jcusp]<>NULL:
#                 for icusp in range(nc):
#                     if ef1[jcusp][icusp]<>NULL:
#                         for n in range(Mv[icusp][2]):
#                             if ef1[jcusp][icusp][n]<>NULL:
#                                 sig_free(ef1[jcusp][icusp][n])
#                         sig_free(ef1[jcusp][icusp])
#                 sig_free(ef1[jcusp])
#         sig_free(ef1)
#     if ef2_c<>NULL:
#         for icusp in range(nc):
#             if ef2_c[icusp]<>NULL:
#                 for n in range(Mv[icusp][2]):
#                     if ef2_c[icusp][n]<>NULL:
#                         sig_free(ef2_c[icusp][n])
#                 sig_free(ef2_c[icusp])
#         sig_free(ef2_c)
#     if nvec<>NULL:
#         for icusp in range(nc):
#             if nvec[icusp]<>NULL:
#                 sig_free(nvec[icusp])
#         sig_free(nvec)
#     if cusp_offsets<>NULL:
#         sig_free(cusp_offsets)
