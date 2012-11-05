# cython: profile=False
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
Cython algorithms for Eisenstein series.
Used by routines in maass_forms.py

"""
include 'sage/ext/stdsage.pxi'
include "sage/ext/cdefs.pxi"
include 'sage/ext/interrupt.pxi'
#include "sage/ext/gmp.pxi"
include "sage/rings/mpc.pxi"
from sage.all import save,incomplete_gamma,load,bessel_K,vector
import mpmath    
from sage.libs.mpfr cimport *

cdef mpc_rnd_t rnd
cdef mpfr_rnd_t rnd_re
rnd = MPC_RNDNN
rnd_re = GMP_RNDN
from sage.rings.complex_mpc cimport MPComplexNumber
from sage.rings.complex_mpc import MPComplexField
from sage.rings.real_mpfr cimport RealNumber,RealField_class
from sage.rings.real_mpfr import RealField
from sage.rings.complex_number cimport ComplexNumber
from sage.rings.complex_field import ComplexField
import sage.structure.element
cimport sage.structure.element
from sage.structure.element cimport Element, ModuleElement, RingElement
from sage.modular.cusps import Cusp
#xfrom sage.modular.cusps cimport Cusp
from sage.rings.infinity import infinity
from sage.rings.integer import Integer,is_Integer
from sage.rings.integer_ring cimport Integer
from sage.rings.complex_double import CDF
from sage.all import log_b
from numpy import array
from sage.all import exp,I,CC,find_root
from sage.matrix.all import MatrixSpace
## Since I want to import sinand cos from math.h 
#from sage.all import sin as ssin
#from sage.all import cos as scos
#from libc.math cimport sin as dsin
#from libc.math cimport cos as dcos
import numpy as np
cimport numpy as cnp
#cimport cython

from maass_forms_alg cimport SMAT_cplx_dp,set_Mv_Qv_symm
from maass_forms_alg import get_M_for_maass_dp
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


cdef double complex cexpi(double x):
    return cexp(x*_I)

cdef double complex _I = _Complex_I


   
from sage.matrix.matrix_dense cimport *
from psage.rings.mpc_extras cimport *
from psage.modules.vector_complex_dense cimport Vector_complex_dense
from psage.modules.vector_real_mpfr_dense cimport Vector_real_mpfr_dense
from psage.matrix.matrix_complex_dense cimport Matrix_complex_dense
from mysubgroups_alg import normalize_point_to_cusp_mpfr,pullback_to_Gamma0N_mpfr,apply_sl2z_map_mpfr,normalize_point_to_cusp_dp,apply_sl2z_map_dp
from pullback_algorithms import pullback_pts_dp,pullback_pts_mpc,pullback_pts_mpc_new
from pullback_algorithms cimport pullback_pts_cplx_dp




cpdef eisenstein_series(S,double sigma,double R,double Y,int M,int Q,int gr=0,int use_fak=0,double eps=1e-12):
    r"""
    Compute the Eisenstein series for the Hecke triangle group G_q
    """
    import mpmath
    cdef double complex **V=NULL
    cdef double* Xm=NULL
    cdef double*** Xpb=NULL
    cdef double*** Ypb=NULL
    cdef double complex*** Cvec=NULL
    cdef double complex **C=NULL
    cdef int *setc_list=NULL
    cdef double complex **RHS=NULL
    cdef double complex **vals_list=NULL
    cdef double *alphas=NULL
    cdef int nc,Ql,Qs,Qf,Ml,Ms,Mf,j,k,l,N,sym_type,i,n,r
    cdef int num_rhs=0
    cdef double *Qfak,tmpr
    cdef double complex *cusp_evs=NULL
    cdef int **Mv=NULL
    cdef int **Qv=NULL
    cdef int verbose=S._verbose
    cdef int *symmetric_cusps=NULL
    cdef int N1=0
    cdef list SetCs
    cdef dict Vals
    cdef int comp_dim,num_set
    cdef int ncols,ncols1
    cdef int cuspidal=1
    cdef int q
    cdef double complex *sqch=NULL
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
    Qfak = <double *>sage_malloc(sizeof(double)*nc)
    Mv=<int**>sage_malloc(sizeof(int*)*nc)
    if not Mv: raise MemoryError
    Qv=<int**>sage_malloc(sizeof(int*)*nc)    
    if not Qv: raise MemoryError
    for i from 0<=i<nc:
        Mv[i]=<int*>sage_malloc(sizeof(int)*3)
        if not Mv[i]: raise MemoryError
        Qv[i]=<int*>sage_malloc(sizeof(int)*3)
        if not Qv[i]: raise MemoryError
    symmetric_cusps=<int*> sage_malloc(sizeof(int)*nc)
    if not symmetric_cusps: raise MemoryError
    cusp_evs=<double complex*>sage_malloc(sizeof(double complex)*nc)
    if not cusp_evs: raise MemoryError
    N = 0; Ml=0; Ql=0
    cdef int* cusp_offsets=NULL
    cusp_offsets=<int*>sage_malloc(sizeof(int)*nc)
    if cusp_offsets==NULL: raise MemoryError
    cdef dict symmetries
    set_Mv_Qv_symm(S,Mv,Qv,Qfak,symmetric_cusps,cusp_evs,cusp_offsets,&N,&Ml,&Ql,M,Q,verbose)
    Qs = 1-Q; Qf = Q
    if verbose>0:
        print "N=",N," Ml=",Ml," Ql=",Ql
    comp_dim=nc
    num_set=0
    for jcusp in range(nc):
        cusp_offsets[jcusp]=0
        for icusp in range(jcusp):
            if icusp==0 or cusp_evs[icusp]==0:
                cusp_offsets[jcusp]+=Mv[icusp][2]
                if verbose>0:
                    print "cusp_offset[",jcusp,"]+=",Mv[icusp][2]
        if verbose>0:
            print "cusp_offset[",jcusp,"]=",cusp_offsets[jcusp]
    ncols = N + comp_dim
    V=<double complex**>sage_malloc(sizeof(double complex*)*N)
    if V==NULL: raise MemoryError
    for j in range(N):
        V[j]=<double complex*>sage_malloc(sizeof(double complex)*(ncols))
        for k in range(ncols): #from 0<=k<ncols:
            V[j][k]=<double complex>0

    RHS=<double complex**>sage_malloc(sizeof(double complex*)*N)
    if RHS==NULL: raise MemoryError
    for n in range(N):
        RHS[n]=NULL
        RHS[n]=<double complex*>sage_malloc(sizeof(double complex)*nc)
        if RHS[n]==NULL: raise MemoryError
    if verbose>0:
        print "In get_coef_cplx_dp_sym R,Y,M,Q=",R,Y,M,Q
    if verbose>1:
        print "N=",N
        print "ncols=",ncols
        print "comp_dim=",comp_dim
        print "num_set=",num_set
        for i from 0<=i<nc:
            print "Ms,Mf,Ml[",i,"]=",Mv[i][0],Mv[i][1],Mv[i][2]
            print "Qs,Qf,Ql[",i,"]=",Qv[i][0],Qv[i][1],Qv[i][2]
    alphas=<double*>sage_malloc(sizeof(double)*nc)
    if alphas==NULL: raise MemoryError
    Xm=<double*>sage_malloc(sizeof(double)*Ql)
    if Xm==NULL: raise MemoryError
    Xpb = <double***> sage_malloc( sizeof(double** ) * nc )
    if Xpb==NULL: raise MemoryError
    Ypb = <double***> sage_malloc( sizeof(double** ) * nc )
    if Ypb==NULL: raise MemoryError
    for i from 0<=i<nc:
        Xpb[i]=NULL; Ypb[i]=NULL
        Xpb[i] = <double**>sage_malloc(sizeof(double*) * nc )
        Ypb[i] = <double**>sage_malloc(sizeof(double*) * nc )
        if Ypb[i]==NULL or Xpb[i]==NULL:
            raise MemoryError
        for j from 0<=j<nc:
            Xpb[i][j]=NULL; Ypb[i][j]=NULL
            Xpb[i][j] = <double*>sage_malloc(sizeof(double) * Ql )
            Ypb[i][j] = <double*>sage_malloc(sizeof(double) * Ql )
            if Ypb[i][j]==NULL or Xpb[i][j]==NULL:
                raise MemoryError
            for n from 0<=n<Ql:
                Xpb[i][j][n]=<double>0
                Ypb[i][j][n]=<double>0
    Cvec = <double complex***>sage_malloc(sizeof(double complex**) * nc )
    if Cvec==NULL: raise MemoryError
    for i in range(nc):
        Cvec[i] = <double complex**>sage_malloc(sizeof(double complex*) * nc )
        if Cvec[i]==NULL:
            raise MemoryError
        for j in range(nc):
            Cvec[i][j] = <double complex*>sage_malloc(sizeof(double complex) * Ql )
            if Cvec[i][j]==NULL:
                raise MemoryError
            for n in range(Ql):
                Cvec[i][j][n]=<double complex>0
    sig_on()
    pullback_pts_cplx_dp(S,Qs,Qf,Y,Xm,Xpb,Ypb,Cvec)
    sig_off()
    for j in range(nc):
        tmpr=float(S.alpha(j)[0])
        alphas[j]=<double>tmpr
        if verbose>0:
            print "alphas[",j,"]=",alphas[j],type(alphas[j])
    sig_on()
    compute_V_cplx_eis_dp_sym(V,RHS,N,Xm,Xpb,Ypb,Cvec,
                              cusp_evs,alphas,Mv,Qv,Qfak,
                              symmetric_cusps,sigma,
                              R,Y,nc,ncols,cuspidal,
                              verbose,use_fak,0)   
    sig_off()
    cdef Matrix_complex_dense VV,RRHS
    if Qfak<>NULL:
        sage_free(Qfak)
    if alphas<>NULL:
        sage_free(alphas)
    if Xm<>NULL:
        sage_free(Xm)
    if Ypb<>NULL:
        for i from 0<=i<nc:
            if Ypb[i]<>NULL:
                for j from 0<=j<nc:
                    if Ypb[i][j]<>NULL:
                        sage_free(Ypb[i][j])
                sage_free(Ypb[i])
        sage_free(Ypb)
    if Xpb<>NULL:
        for i from 0<=i<nc:
            if Xpb[i]<>NULL:
                for j from 0<=j<nc:
                    if Xpb[i][j]<>NULL:
                        sage_free(Xpb[i][j])
                sage_free(Xpb[i])
        sage_free(Xpb)
    if Cvec<>NULL:
        for i from 0<=i<nc:
            if Cvec[i]<>NULL:
                for j from 0<=j<nc:
                    if Cvec[i][j]<>NULL:
                        sage_free(Cvec[i][j])
                sage_free(Cvec[i])
        sage_free(Cvec)
    if sqch<>NULL:
        sage_free(sqch)
    if gr==1:
        CF=MPComplexField(53)
        MS=MatrixSpace(CF,N,N)
        VV=Matrix_complex_dense(MS,0)
        MS=MatrixSpace(CF,N,nc)
        RRHS=Matrix_complex_dense(MS,0)
        for n in range(N):
            for l in range(N):
                VV[n,l]=V[n][l]
            for l in range(nc):
                RRHS[n,l]=RHS[n][l]
        return VV,RRHS
    C=<double complex**>sage_malloc(sizeof(double complex*)*(comp_dim))
    if C==NULL:
        raise MemoryError
    for i in range(comp_dim):
        C[i]=<double complex*>sage_malloc(sizeof(double complex)*(N))
        if C[i]==NULL:
            raise MemoryError
    #setc_list = <int*>sage_malloc(sizeof(int)*num_set)
    #if setc_list==NULL:
    #    raise MemoryError
    #vals_list = <double complex** >sage_malloc(sizeof(double complex*)*comp_dim)
    #if vals_list==NULL:
    #    raise MemoryError
    #for j in range(comp_dim):
    #    vals_list[j]=<double complex*>sage_malloc(sizeof(double complex)*num_set) 
    i=0
    # for r,n in SetCs:
    #     if cusp_evs[r]<>0 and r>0:
    #         continue
    #     l = cusp_offsets[r]+n-Mv[r][0]
    #     if l<0:
    #         continue
    #     setc_list[i]=l
    #     for j from 0<=j<comp_dim:
    #         vals_list[j][i]=<double complex>CC(Vals[j][(r,n)])
    #     i=i+1
    if num_rhs>0 and num_rhs<>comp_dim:
        raise ValueError,"Need same number of right hand sides (or just one) as the number of set coefficients!"
    #normalize_matrix_cplx_dp(V1,N,comp_dim,num_set,setc_list,vals_list,verbose)
    for n in range(N):
        for jcusp in range(comp_dim):            
            V[n][N+jcusp]=-RHS[n][jcusp]
    if gr==2:
        CF=MPComplexField(53)
        MS=MatrixSpace(CF,N,ncols)
        VV=Matrix_complex_dense(MS,0)
        for n from 0<=n< N:
            for l from 0 <= l < ncols:
                VV[n,l]=V[n][l]
        return VV
    #sig_on()
    SMAT_cplx_dp(V,N,comp_dim,0,C,vals_list,setc_list)
    #sig_off()
    if verbose>1:
        for k in range(ncols):
            print "C[",k,"]=",C[0][k]
    cdef dict res={}
    cdef int ki
    for j from 0<=j<comp_dim:
        res[j]=dict()
        for i from 0<=i<nc:
            if i>0 and cusp_evs[i]<>0:
                res[j][i]=cusp_evs[i]
                continue
            res[j][i]=dict()
            for k from 0<=k<Mv[i][2]:
                ki=cusp_offsets[i]+k
                res[j][i][k+Mv[i][0]]=C[j][ki]

    #if setc_list<>NULL:
    #    sage_free(setc_list)
    #if vals_list<>NULL:
    #    for j from 0<=j<comp_dim:
    #        if vals_list[j]<>NULL:
    #            sage_free(vals_list[j])
    #    sage_free(vals_list)
    if V<>NULL:
        for j from 0<=j<N:
            if V[j]<>NULL:
                sage_free(V[j])
        sage_free(V)
    if C<>NULL:
        for i from 0<=i<comp_dim:
            if C[i]<>NULL:
                sage_free(C[i])
        sage_free(C)
    if RHS<>NULL:
        for n in range(N):
            if RHS[n]<>NULL:
                sage_free(RHS[n])
        sage_free(RHS)
    if Mv<>NULL:
        sage_free(Mv)
    if Qv<>NULL:
        sage_free(Qv)                
    if symmetric_cusps<>NULL:
        sage_free(symmetric_cusps)
    if cusp_evs<>NULL:
        sage_free(cusp_evs)
    if comp_dim==1:
        return res[0]
    return res


@cython.boundscheck(False)
@cython.cdivision(True)
cdef compute_V_cplx_eis_dp_sym(double complex **V,
                               double complex **RHS,
                           int N,
                           double *Xm,
                           double ***Xpb,
                           double ***Ypb,
                           double complex ***Cvec,
                           double complex *cusp_evs,
                           double *alphas,
                           int **Mv,int **Qv,double *Qfak,
                           int *symmetric_cusps,
                           double sigma,
                           double R,double Y,
                           int nc, int ncols,
                           int cuspidal,
                           int verbose, int use_fak=0,
                           int is_trivial=0):


    r"""
    Set up the matrix for the system of equations giving the Fourier coefficients of a non-holomorphic Eisenstein series with
    spectral parameter s=sigma+iR
    ( we set up a system for all Eisensteins series wrt all cusps)
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
    #mpmath.mp.dps=25
    cdef int l,j,icusp,jcusp,n,ni,lj,Ml,Ql,Qs,Qf,Mf,Ms,prec
    cdef double pi,sqrtY,Y2pi,nrY2pi,argm,argpb,twopi,two,besarg,lr,nr
    cdef double complex ckbes,ctmpV,iargm,twopii,kbes,ctmp
    cdef double RlnY,x,y
    cdef double complex Ys
    cdef mpc_t ztmp,stmp
    cdef mpfr_t xtmp,ytmp
    if not cuspidal in [0,1]:
        raise ValueError," parameter cuspidal must be 0 or 1"
    pi=M_PI #<double>RealField(53).pi() #3.141592653589793238
    sqrtY=sqrt(Y)
    RlnY=R*log(Y)
    prec = 103
    mpc_init2(ztmp,prec);    mpc_init2(stmp,prec)
    mpfr_init2(xtmp,prec);   mpfr_init2(ytmp,prec)
    mpfr_set_d(xtmp,sigma,rnd_re)
    mpfr_set_d(ytmp,R,rnd_re)
    mpc_set_fr_fr(stmp,xtmp,ytmp,rnd)
    mpfr_set_d(ytmp,Y,rnd_re)
    mpc_set_fr(ztmp,ytmp,rnd)    
    mpc_pow(ztmp,ztmp,stmp,rnd)
    mpfr_set(xtmp,mpc_realref(ztmp),rnd_re)
    mpfr_set(ytmp,mpc_imagref(ztmp),rnd_re)
    x = mpfr_get_d(xtmp,rnd_re)
    y = mpfr_get_d(ytmp,rnd_re)
    mpc_clear(ztmp); mpc_clear(stmp)
    mpfr_clear(xtmp); mpfr_clear(ytmp)
    Ys = x + _I*y
    if verbose>0:
        print "Ys=",Ys
        print "Y**s=",Y**<double complex>(sigma+_I*R)
    mpy=mpmath.mp.mpf(Y)
    mpsqY=mpmath.mp.sqrt(mpy)
    mpone_minus_s=mpmath.mp.mpc(1.0-sigma,-R)
    mp_y_1ms=mpmath.mp.power(mpy,mpone_minus_s)
    two=<double>(2)
    Y2pi=Y*pi*two
    twopi=two*pi
    Ml=0; Ql=0
    for i from 0<=i<nc:
        if Mv[i][2]>Ml:
            Ml=Mv[i][2]
        if Qv[i][2]>Ql:
            Ql=Qv[i][2]
        if verbose>0:
            print "Qv[",i,"]=(",Qv[i][0],",",Qv[i][1],",",Qv[i][2],")"
            print "Mv[",i,"]=(",Mv[i][0],",",Mv[i][1],",",Mv[i][2],")"
    if verbose>0:
        print "N=",N
        print "Ql=",Ql
    ## This is the effective offset at the 
    cdef int* cusp_offsets=NULL
    cusp_offsets=<int*>sage_malloc(sizeof(int)*nc)
    if cusp_offsets==NULL: raise MemoryError
    for jcusp from 0 <= jcusp < nc:
        cusp_offsets[jcusp]=0
        for icusp from 0 <= icusp < jcusp:
            if icusp==0 or cusp_evs[icusp]==0:
                cusp_offsets[jcusp]+=Mv[icusp][2]
        if verbose>0:
            print "cusp_offsets[",jcusp,"]=",cusp_offsets[jcusp]
    cdef double **nvec=NULL
    nvec = <double**>sage_malloc(sizeof(double*)*nc)
    if not nvec: raise MemoryError
    for icusp in range(nc):
        nvec[icusp] = <double*>sage_malloc(sizeof(double)*Ml)
    cdef double complex ****ef1=NULL
    cdef double complex ***ef2_c=NULL
    ef2_c = <double complex***>sage_malloc(sizeof(double complex**)*nc)
    if not ef2_c: raise MemoryError
    for icusp  in range(nc):
        ef2_c[icusp] = <double complex**>sage_malloc(sizeof(double complex*)*Mv[icusp][2])
        for n  in range(Mv[icusp][2]):
            ef2_c[icusp][n] = <double complex*>sage_malloc(sizeof(double complex)*Qv[icusp][2])
    ef1 = <double complex****>sage_malloc(sizeof(double complex***)*nc)
    if ef1==NULL: raise MemoryError        
    for icusp  in range(nc):
        ef1[icusp] = <double complex***>sage_malloc(sizeof(double complex**)*nc) 
        if ef1[icusp]==NULL: raise MemoryError        
        for jcusp in range(nc):
            ef1[icusp][jcusp] = <double complex**>sage_malloc(sizeof(double complex*)*Mv[jcusp][2])
            if ef1[icusp][jcusp]==NULL: raise MemoryError        
            for n from 0<=n<Mv[jcusp][2]:
                ef1[icusp][jcusp][n] = <double complex*>sage_malloc(sizeof(double complex)*Qv[jcusp][2])
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
                    ef2_c[jcusp][n][j]=cos(argm)
                elif symmetric_cusps[jcusp]==1:
                    ef2_c[jcusp][n][j]=_I*sin(-argm)
                else:
                    ef2_c[jcusp][n][j]=cexpi(-argm)
    cdef double argpb1
    for jcusp in range(nc):
        for icusp in range(nc):
            for n from 0<=n<Mv[jcusp][2]:
                for j from 0<=j<Qv[jcusp][2]: #in range(Qs,Qf+1):
                    if Ypb[icusp][jcusp][j]==0:
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
    #cdef double y #,kbes_old=1.0
    cdef double complex ***kbesvec=NULL
    kbesvec=<double complex***>sage_malloc(sizeof(double complex**)*nc)
    if kbesvec==NULL:
        raise MemoryError
    for jcusp in range(nc):
        kbesvec[jcusp]=<double complex**>sage_malloc(sizeof(double complex*)*Ml)
        if kbesvec[jcusp]==NULL:
                raise MemoryError
        for l in range(Ml):
            kbesvec[jcusp][l]=<double complex*>sage_malloc(sizeof(double complex)*Ql)
            if kbesvec[jcusp][l]==NULL:
                raise MemoryError
    if verbose>0:
        print "here0"
        print "Ml=",Ml
        print "Ql=",Ql
    cdef double besprec
    besprec=1.0E-14
    cdef double complex s,one_minus_s,s_minus_half
    cdef double fak
    s = CC(sigma,R)
    s_minus_half=CC(sigma-0.5,R)
    one_minus_s=CC(1.0-sigma,-R) #-s
    if use_fak==1:
        fak = mpmath.mp.exp(R*mpmath.mp.pi*0.5)
    else:
        fak = 1
    for jcusp in range(nc):
        for icusp in range(nc):
            for j in range(Qv[jcusp][2]):            
                ypb=Ypb[icusp][jcusp][j]
                if ypb==0:
                    continue
                for l in range(Mv[jcusp][2]):
                    lr=nvec[jcusp][l]*twopi
                    besarg=fabs(lr)*ypb
                    if lr<>0.0:
                        try:
                            kbes=mpmath.mp.besselk(s_minus_half,besarg)
                        except:
                            print "s=",s,type(s)
                            print "besarg=",besarg,type(besarg)
                            raise ArithmeticError
                        kbesvec[icusp][l][j]=sqrt(ypb)*kbes*fak
                    else:
                        if abs(s-0.5)>1E-14:
                            kbesvec[icusp][l][j]=ypb**one_minus_s
                        else:
                            kbesvec[icusp][l][j]=sqrt(ypb)*log(ypb)
    #cdef double complex cuspev
    for jcusp in range(nc):
        for l in range(Mv[jcusp][2]):
            lr=nvec[jcusp][l]*twopi
            lj=cusp_offsets[jcusp]+l
            for j in range(Qv[jcusp][2]):
                for icusp in range(nc):
                    if Ypb[icusp][jcusp][j]==0: 
                        continue
                    ckbes=kbesvec[icusp][l][j]*ef1[icusp][jcusp][l][j]
                    for n in range(Mv[icusp][2]):
                        ni=cusp_offsets[icusp]+n
                        ctmpV=ckbes*ef2_c[icusp][n][j]
                        if ni>N or lj>N:
                            raise ArithmeticError,"Index outside!"
                        V[ni][lj]=V[ni][lj]+ctmpV
    if verbose>3:
        print "V0[",3,3,"]=",V[3][3]
    for icusp in range(nc): 
        for n in range(Mv[icusp][2]):
            ni=cusp_offsets[icusp]+n
            for jcusp in range(nc):
                for l in range(Mv[jcusp][2]):
                    lj=cusp_offsets[jcusp]+l
                    if lj>N or ni>N: # some extra insurance...
                        print "ni=",cusp_offsets[icusp],"+",n,"=",ni
                        print "lj=",cusp_offsets[jcusp],"+",l,"=",lj
                        raise ArithmeticError,"Index outside!"                    
                    V[ni][lj]=V[ni][lj]/Qfak[jcusp]
    cdef double complex summa,term
    #summa = CF(0); term = CF(0)
    for icusp in range(nc):
        for n in range(Mv[icusp][2]):
            nr=fabs(nvec[icusp][n])
            ni=cusp_offsets[icusp]+n
            if nr==0.0:
                 if abs(s-0.5)>1E-14:
                     #kbes=Y**one_minus_s
                     kbes=mp_y_1ms 
                 else:
                     kbes=sqrtY*log(Y)
            else:
                besarg=nr*Y2pi
                try:
                    #kbes=bessel_K(s,besarg)
                    kbes=mpmath.fp.besselk(s_minus_half,besarg)
                except:
                    print "s=",s,type(s)
                    print "besarg=",besarg,type(besarg)
                    raise ArithmeticError
                kbes=sqrtY*kbes*fak # besselk_dp(R,nrY2pi,pref=1)
            if ni>N:
                raise ArithmeticError,"Index outside!"
            V[ni][ni]=V[ni][ni] - kbes
            ## Also get the RHS:
            for jcusp in range(nc):
                if jcusp<>icusp:
                    RHS[ni][jcusp]=0
                else:
                    summa=0.0
                    for j in range(Qv[jcusp][2]):            
                        ypb=Ypb[icusp][jcusp][j]
                        if abs(s-0.5)>1E-14:
                            #term = ypb**s
                            term = mpmath.mp.power(mpmath.mp.mpf(ypb),s)
                        else:
                            term = sqrt(ypb)
                        term = term*ef2_c[jcusp][n][j]
                        summa = summa + term
                    RHS[ni][jcusp]=summa/Qfak[jcusp]
                if nr==0 and jcusp==icusp:
                    ## Subtracting of....
                    if abs(s-0.5)>1E-14:
                        #kbes=mpmath.mp.power(Y,s)
                        #kbes = (Y**sigma)
                        #kbes = kbes * (cos(RlnY)+_I*sin(RlnY))
                        kbes = Ys
                    else:
                        kbes=mpmath.mp.sqrt(Y)
                    RHS[ni][jcusp]=RHS[ni][jcusp]-kbes
                if verbose>2:
                    print 'RHS[',ni,jcusp,']=',RHS[ni][jcusp]
    if verbose>2:
        print "V3[",3,3,"]=",V[3][3]
    if kbesvec<>NULL:
        for icusp from 0<=icusp<nc:
            if kbesvec[icusp]<>NULL:
                for l from 0<=l<Ml:
                    if kbesvec[icusp][l]<>NULL:
                        sage_free(kbesvec[icusp][l])
                sage_free(kbesvec[icusp])
        sage_free(kbesvec)
    #print "deal kbbes1"
    if ef1<>NULL:
        for jcusp from 0<=jcusp<nc:
            if ef1[jcusp]<>NULL:
                for icusp from 0<=icusp<nc:
                    if ef1[jcusp][icusp]<>NULL:
                        for n from 0<=n<Mv[icusp][2]:
                            if ef1[jcusp][icusp][n]<>NULL: 
                                sage_free(ef1[jcusp][icusp][n])
                        sage_free(ef1[jcusp][icusp])
                sage_free(ef1[jcusp])                                
        sage_free(ef1)
    if ef2_c<>NULL:
        for icusp from 0<=icusp<nc:
            if ef2_c[icusp]<>NULL:
                for n from 0<=n<Mv[icusp][2]:
                    if ef2_c[icusp][n]<>NULL:
                        sage_free(ef2_c[icusp][n])
                sage_free(ef2_c[icusp])
        sage_free(ef2_c)
    if nvec<>NULL:
        for icusp from 0<=icusp<nc:
            if nvec[icusp]<>NULL:
                sage_free(nvec[icusp])
        sage_free(nvec)
    if cusp_offsets<>NULL:
        sage_free(cusp_offsets)





### Algorithms for testing the functional equation of the Selberg zeta function


from pullback_algorithms import pullback_pts_mp
from sage.functions.special import Bessel

cpdef Eisenstein_series_one_cusp(S,RealNumber sigma,RealNumber R,RealNumber Y,int M,int verbose=0,int gr=0,int tprec=0):
    cdef int Q,Qs,Qf,Ql,Ms,Mf,Ml,n,j,l,prec
    cdef RealNumber llambda
    if tprec==0:
        prec = sigma.prec()
    else:
        prec = tprec
    # Make sure that the precision is the same for all parameters
    RF = RealField(prec)
    CF = MPComplexField(prec)
    sigma = RF(sigma); R = RF(R); Y = RF(Y)
    Q = M + 10
    Qs = 1-Q; Qf = Q; Ql=2*Q
    Ms = -M; Mf = M; Ml=2*M+1
    ## Recall that we have scaled the fundamental domain by 1/lambda
    cdef double Ymax
    llambda = RF(2)*(RF.pi()/RF(S._group._q)).cos()
    Ymax=<double>S._group.minimal_height()/<double>float(llambda)
    if Y>Ymax:
        raise ValueError,"Need Y<{0}. We have Y={1}".format(Ymax,Y)
    pb = pullback_pts_mp(S,Qs,Qf,Y)
    Xm = pb['xm']
    Xpb = pb['xpb']
    Ypb = pb['ypb']
    MS = MatrixSpace(CF,Ml,Ml)
    cdef Matrix_complex_dense V
    V = Matrix_complex_dense(MS,0)
    mpmath.mp.prec=CF.prec()
    cdef RealNumber pi,twopi,mpsqY,one_minus_sigma
    cdef MPComplexNumber cmpY
    pi = RF.pi(); twopi = RF(2)*pi
    mpY = RF(Y); cmpY = CF(Y) 
    mpsqY = mpY.sqrt();  s = CF(sigma,R) 
    one_minus_sigma=RF(1-sigma)
    cdef ComplexNumber s2
    CCF = ComplexField(prec)
    s2 = CCF(sigma,R)-CCF(1,0)/CCF(2,0)
    s2_mp = mpmath.mp.mpc(sigma-0.5,R)
    #MyKBessel = Bessel(s2,typ='K',prec=prec)    
    one_minus_s=CF(1,0)-s #mpmath.mp.mpc(1,0)-s
    cdef RealNumber xpb,ypb
    cdef MPComplexNumber argpb,argm,tmpc,kbes,summa,zero
    xpb = RF(0); ypb = RF(0); kbes=CF(0)
    argpb=CF(0); tmpc=CF(0); argm=CF(0); summa=CF(0); zero=CF(0)
    cdef mpfr_t nn,rz
    mpfr_init2(nn,prec);    mpfr_init2(rz,prec)
    mpfr_set_si(rz,0,rnd_re)
    cdef mpc_t **ef2=NULL
    ef2=<mpc_t**>sage_malloc(sizeof(mpc_t*)*Ml)
    if ef2==NULL: raise MemoryError
    for n in range(Ml):
        ef2[n]=<mpc_t*>sage_malloc(sizeof(mpc_t)*Ql)
        if ef2[n]==NULL: raise MemoryError
        mpfr_set_si(nn,n+Ms,rnd_re)
        mpfr_neg(nn,nn,rnd_re)
        for j in range(Ql):
            mpc_set_fr_fr(argm.value,rz,nn,rnd)
            mpc_mul_fr(argm.value,argm.value,(<RealNumber>Xm[j]).value,rnd)
            mpc_exp(argm.value,argm.value,rnd)
            mpc_init2(ef2[n][j],prec)
            mpc_set(ef2[n][j],argm.value,rnd)
        #argm = CF(0,-nn*Xm[j])
    
    cdef RealNumber logypb,ccr,cci,RY,RS,besarg
    logypb=RF(0);ccr=RF(0); cci=RF(0); RY=RF(0);RS=RF(0)
    besarg = RF(0)
    for l in range(Ml):
        ll=RF(l+Ms)
        for j in range(Ql):
            ypb = Ypb[0][0][j]
            xpb = Xpb[0][0][j]
            if l+Ms==0:
                logypb=ypb.log()
                RY=R*logypb
                RS=one_minus_sigma*logypb
                RS=RS.exp()
                ccr=(RY).cos()
                cci=(RY).sin()
                kbes = RS*CF(ccr,-cci)
                #kbes = ypb**one_minus_s
            else:
                if l+Ms>0:
                    besarg = ypb*ll*twopi
                else:
                    besarg = -ypb*ll*twopi
                bes = mpmath.mp.besselk(s2_mp,besarg)
                kbes = CF(bes.real,bes.imag)*ypb.sqrt()
            argpb = CF(0,ll*xpb)
            kbes = kbes*argpb.exp()
            for n in range(Ml):
                mpc_mul(tmpc.value,kbes.value,ef2[n][j],rnd)
                mpc_add(V._matrix[n][l],V._matrix[n][l],tmpc.value,rnd)
    for n in range(Ml):
        for l in range(Ml):
            mpc_div_ui(V._matrix[n][l],V._matrix[n][l],2*Q,rnd)
    for n in range(Ml):
        if n+Ms==0:
            logypb=mpY.log()
            RY=R*logypb
            RS=one_minus_sigma*logypb
            RS=RS.exp()
            ccr=(RY).cos()
            cci=(RY).sin()
            kbes = RS*CF(ccr,-cci)
            #kbes = cmpY**one_minus_s #mpmath.mp.power(mpY,one_minus_s)
        else:
            mpfr_set_si(nn,n+Ms,rnd_re)
            mpfr_mul(nn,nn,twopi.value,rnd_re)
            mpfr_abs(nn,nn,rnd_re)
            mpfr_mul(besarg.value,Y.value,nn,rnd_re)
            bes = mpmath.mp.besselk(s2_mp,besarg)
            kbes = mpsqY*CF(bes.real,bes.imag)
                
        #bess = CF(bes.real,bes.imag)
        mpc_sub(V._matrix[n][n],V._matrix[n][n],kbes.value,rnd)
    #MS = MatrixSpace(CF,Ml,1)
    cdef Vector_complex_dense RHS
    VS = vector(CF,Ml).parent()
    RHS = Vector_complex_dense(VS,0)    
    for j in range(Ql):
        ypb = Ypb[0][0][j]
        logypb=ypb.log()
        RY=R*logypb
        RS=sigma*logypb
        RS=RS.exp()
        ccr=(RY).cos()
        cci=(RY).sin()
        kbes = RS*CF(ccr,cci)
        #if abs(kbes-ypb**s)>0.01:
        #    print "kbes=",kbes
        #    print "ypb**s=",ypb**s
        #    print "ypb=",ypb
        #    raise ArithmeticError
        #kbes = ypb**s
        for n in range(Ml):    
            mpc_mul(tmpc.value,kbes.value,ef2[n][j],rnd)
            mpc_add(RHS._entries[n],RHS._entries[n],tmpc.value,rnd)
    for n in range(Ml):    
        mpc_div_ui(RHS._entries[n],RHS._entries[n],2*Q,rnd)
        if n+Ms==0:
            logypb=mpY.log()
            RY=R*logypb
            RS=sigma*logypb
            RS=RS.exp()
            ccr=(RY).cos()
            cci=(RY).sin()
            kbes = RS*CF(ccr,cci)
            #print "err0=",abs(kbes-mpY**s)
            #kbes = cmpY**s #mpmath.mp.power(mpY,s)
            mpc_sub(RHS._entries[n],RHS._entries[n],kbes.value,rnd)
        mpc_neg(RHS._entries[n],RHS._entries[n],rnd)
    mpfr_clear(rz); mpfr_clear(nn)
    if ef2<>NULL:
        for n in range(Ml):
            if ef2[n]<>NULL:
                for j in range(Ql):
                    mpc_clear(ef2[n][j])
                sage_free(ef2[n])
        sage_free(ef2)
    #B = RHS.column(0)
    if gr==1:
        return V,RHS
    else:
        X1 = V.solve(RHS)
        X = dict()
        for n in range(Ml):
            X[n+Ms]=X1[n]
        return X

#from sage.libs.mpmath.utils import mpfr_to_mpfval,mpfr_from_mpfval
#from sage.libs.mpmath.utils import mpfr_to_mpfval,mpfr_from_mpfval
#from sage.libs.mpmath.ext_main import make_mpc,make_mpf

cpdef mpmath_to_mpc(z,prec):
    cdef MPComplexNumber resz
    cdef RealNumber resx
    if hasattr(z, "_mpf_"):
        resx = RealField(prec)()
        mpfr_from_mpfval(resx.value, z._mpf_)
        return resx
    elif hasattr(z, "_mpc_"):
        resz = MPComplexField(prec)(0)
        re, im = z._mpc_
        
        mpfr_from_mpfval(resz.value.re, re)
        mpfr_from_mpfval(resz.value.im, im)
        return resz
    else:
        raise TypeError("cannot convert %r to Sage", z)


cdef mpfr_t_to_mpmath(mpfr_t x):
    return mpmath.mp.make_mpf(mpfr_to_mpfval(x))
    
cdef mpc_t_to_mpmath(mpc_t z):
    return mpmath.mp.make_mpc(mpfr_to_mpfval(z.re),mpfr_to_mpfval(z.im))

cpdef mpc_to_mpmath(MPComplexNumber z,prec):
    return mpc_t_to_mpmath(z.value)

cpdef mpfr_to_mpmath(RealNumber x,prec):
    return mpfr_t_to_mpmath(x.value)

from sage.rings.integer_ring cimport Integer
cdef mpfr_from_mpfval(mpfr_t res, tuple x):
    """
    Set value of an MPFR number (in place) to that of a given mpmath mpf
    data tuple.
    """
    cdef int sign
    cdef Integer man
    cdef long exp
    cdef long bc
    sign, man, exp, bc = x
    if man.__nonzero__():
        mpfr_set_z(res, man.value, GMP_RNDZ)
        if sign:
            mpfr_neg(res, res, GMP_RNDZ)
        mpfr_mul_2si(res, res, exp, GMP_RNDZ)
        return
    from mpmath.libmp import finf, fninf
    if exp == 0:
        mpfr_set_ui(res, 0, GMP_RNDZ)
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
    cdef Integer man = PY_NEW(Integer)
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


#@cython.cdivision(True) 
cpdef ComplexNumber besselk_mpc_rec(RealNumber sigma, RealNumber R,RealNumber  x,double eps=1e-14,int pref=0):   # |r| <<1000, x >> 0 !
    r"""
    Modified K-Bessel function in double precision using the backwards Miller-recursion algorithm. 

    INPUT
     - ''R'' -- RealNumber
     - ''sigma'' -- RealNumber
     - ''x'' -- RealNumber
     - ''prec'' -- (default 1E-12) double
     - ''pref'' -- (default 1) int
                =1 => computes K_iR(x)*exp(pi*R/2)
                =0 => computes K_iR(x)

    OUTPUT:
     - exp(Pi*R/2)*K_{i*R}(x)  -- double
     
    EXAMPLES::



    """
    cdef ComplexNumber p,s,llambda,rr
    cdef int prec
    prec = R.prec()
    CF = ComplexField(prec)
    RF = RealField(prec)
    s = CF(sigma,R); p = CF(0)
    llambda = s*(RF(1)-s)
    rr = (llambda-RF(1)/RF(4)).sqrt()
    if x<0:
        raise ValueError,"X<0"
    #! To make a specific test is time-consuming
    cdef RealNumber q,err,one,ef1,two,pi,mr
    one=RF(1); two=RF(2)
    pi=RF.pi()
    p=llambda #RF(1)/RF(4)+R*R
    q=RF(2)*(x-one)
    print "p=",p
    print "rr=",rr
    err=one
    cdef int NMAX=5000
    cdef int n_start=40 #128
    ef1=(two*x/pi).log()
    mr=RF(0)
    cdef int n=n_start
    cdef RealNumber nr,mr_p1,mp_m1
    cdef ComplexNumber cone,d,tmp,tmp2,ef,efarg,den,k,y,t
    tmp=CF(0); tmp2=CF(0); ef=CF(0); nr=RF(0)
    t=CF(0); cone=CF(1); mr_p1=RF(1); mr_m1=RF(0) 
    k=CF(1); y = CF(0) 
    d = CF(1)
    mr_m1=CF(0); efarg = CF(0); den=CF(0);
    cdef int nn
    for nn in range(1,NMAX+1):
        err=abs(t-k)
        if err < eps:
            if n>n_start+40:
                break
        n=n+20       
        t=k 
        y=cone #               !/*arbitrary*/1; 
        k=cone
        d=cone
        tmp=two*x-rr*pi
        nr=RF(n)
        if tmp>1300.0:
            return RF(0)
        ef = ((ef1 + tmp)/(two*nr)).exp()
        mr_p1=RF(n+1)   # m+1
        mr=RF(n)        # m 
        for m from n >=m>=1:
            mr_m1=RF(m-1) # m-1
            den = (q+(mr_p1)*(two-y))
            y=(mr_m1+p/mr)/den
            k=ef*(d+y*k)
            d=d*ef
            mr_p1=mr
            mr=mr_m1

        if k==0:
            st='the K-bessel routine failed (k large) for x,R=%s,%s, value=%s' 
            raise ValueError,st%(x,R,k)
        k=k**-1
    if  abs(k) > 1E30: #something big... 
        st='the K-bessel routine failed (k large) for x,R=%s,%s, value=%s' 
        raise ValueError,st%(x,R,k)
    if(nn>=NMAX):
        st='the K-bessel routine failed (too many iterations) for x,R=%s,%s, value=%s' 
        raise ValueError,st%(x,R,k)
    #!k= exp( pi*r/2) * K_ir( x) !
    if(pref==1):
        return k
    else:
        return k*exp(-pi*rr/two)
