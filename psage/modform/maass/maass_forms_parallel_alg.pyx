# cython: profile=True
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

Parallel routines for computing Maass forms.

"""
include "stdsage.pxi"  
include "cdefs.pxi"

cdef extern from "stdio.h":
    cdef extern void printf(char *fmt,...) nogil
    
from psage.rings.mpfr_nogil cimport *
import cython
from cython.parallel cimport parallel, prange
from lpkbessel cimport besselk_dp_c
cimport openmp

cdef extern from "math.h" nogil:
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

cdef extern from "complex.h" nogil:
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
    return cexp(x*_Complex_I)







@cython.boundscheck(False)
@cython.cdivision(True)
cdef int compute_V_cplx_dp_sym_par(double complex **V,
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
                           int ncpus=1,
                           int is_trivial=0):


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
    - ``ncpus`` -- the number of cpus (cores)
    """
    cdef int l,i,j,icusp,jcusp,n,ni,lj,Ml,Ql,s,Qs,Qf,Mf,Ms
    cdef double pi,sqrtY,Y2pi,nrY2pi,argm,argpb,twopi,two,kbes,besarg,lr,nr
    cdef double complex ckbes,ctmpV,iargm,twopii,ctmp
    if not cuspidal in [0,1]:
        raise ValueError," parameter cuspidal must be 0 or 1"
    cdef int **endpts
    cdef int* cusp_offsets=NULL
    cdef double **nvec=NULL
    cdef double complex ****ef1=NULL
    cdef double complex ***ef2_c=NULL
    cdef double ***kbesvec=NULL
    ## Set constants
    pi=M_PI 
    sqrtY=sqrt(Y)
    two=<double>(2)
    Y2pi=Y*pi*two
    twopi=two*pi
    Ml=0; Ql=0
    if verbose>0:
        printf("in compute_V_cplx_dp_sym_par \n")
    for i in range(nc):
        if Mv[i][2]>Ml:
            Ml=Mv[i][2]
        if Qv[i][2]>Ql:
            Ql=Qv[i][2]
        if verbose>0:
            printf("symmetric cusp[%d]=%d\n",i,symmetric_cusps[i])
            printf("Qv[%d](%d,%d,%d)\n",i,Qv[i][0],Qv[i][1],Qv[i][2])
            printf("Mv[%d](%d,%d,%d)\n",i,Mv[i][0],Mv[i][1],Mv[i][2])
    if verbose>0:
        printf("Ql=%d \n",Ql)
    ## Allocate arrays 
    cusp_offsets=<int*>sage_malloc(sizeof(int)*nc)
    if cusp_offsets==NULL: raise MemoryError
    nvec = <double**>sage_malloc(sizeof(double*)*nc)
    if not nvec: raise MemoryError
    for icusp in range(nc):
        nvec[icusp] = <double*>sage_malloc(sizeof(double)*Ml)
    ef2_c = <double complex***>sage_malloc(sizeof(double complex**)*nc)
    if not ef2_c: raise MemoryError
    for icusp in range(nc):
        ef2_c[icusp] = <double complex**>sage_malloc(sizeof(double complex*)*Mv[icusp][2])
        for n in range(Mv[icusp][2]):
            ef2_c[icusp][n] = <double complex*>sage_malloc(sizeof(double complex)*Qv[icusp][2])
    ef1 = <double complex****>sage_malloc(sizeof(double complex***)*nc)
    if ef1==NULL: raise MemoryError
    for icusp in range(nc):
        ef1[icusp] = <double complex***>sage_malloc(sizeof(double complex**)*nc)
        if ef1[icusp]==NULL: raise MemoryError
        for jcusp in range(nc):
            ef1[icusp][jcusp] = <double complex**>sage_malloc(sizeof(double complex*)*Mv[jcusp][2])
            if ef1[icusp][jcusp]==NULL: raise MemoryError
            for n in range(Mv[jcusp][2]):
                ef1[icusp][jcusp][n] = <double complex*>sage_malloc(sizeof(double complex)*Qv[jcusp][2])
                if ef1[icusp][jcusp][n]==NULL: raise MemoryError
    endpts = <int**>sage_malloc(ncpus*sizeof(int*))
    if endpts==NULL:
        raise MemoryError
    for i in range(ncpus):
        endpts[i]=<int*>sage_malloc(2*sizeof(int))

    ## Assigning values to arrays
    for jcusp  in range(nc):
        cusp_offsets[jcusp]=0
        for icusp in range(jcusp):
            if icusp==0 or cusp_evs[icusp]==0:
                cusp_offsets[jcusp]+=Mv[icusp][2]
        if verbose>0:
            printf("cusp_offsets[%d]=%d \n",jcusp,cusp_offsets[jcusp])
    cdef int nc_sym=0
    for jcusp  in range(nc):
        if verbose>0:
            printf("cusp_evs[%d]=%d \n",jcusp,cusp_evs[jcusp])
        if jcusp==0 or cusp_evs[jcusp]<>0:
            nc_sym+=1
    for jcusp in range(nc):
        for n in range(Ml):
            nvec[jcusp][n]=<double>(n+Mv[jcusp][0])+alphas[jcusp]
    cdef int twoQm1
    twoQm1= 2*Qv[0][1]-1
    for jcusp in range(nc):
        #printf("symmetric cusps[%d]=%d \n",jcusp,symmetric_cusps[jcusp])        
        for n in range(Mv[jcusp][2]):
            for j in range(Qv[jcusp][2]):
                argm=nvec[jcusp][n]*Xm[j]
                if symmetric_cusps[jcusp]==0:
                    ef2_c[jcusp][n][j]=cos(argm)
                elif symmetric_cusps[jcusp]==1:
                    ef2_c[jcusp][n][j]=_Complex_I*sin(-argm)
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
                        ef1[icusp][jcusp][n][j]=_Complex_I*sin(argpb)
                    else:
                        ef1[icusp][jcusp][n][j]=cexpi(argpb)
                    ctmp = Cvec[icusp][jcusp][j]
                    ef1[icusp][jcusp][n][j]=ef1[icusp][jcusp][n][j]*ctmp

    if verbose>1:
        printf("here1121")
    cdef double besarg_old=0.0
    cdef double y,kbes_old=1.0
    kbesvec=<double***>sage_malloc(sizeof(double**)*nc)
    if kbesvec==NULL:
        raise MemoryError
    for jcusp  in range(nc):
        kbesvec[jcusp]=<double**>sage_malloc(sizeof(double*)*Ml)
        if kbesvec[jcusp]==NULL:
            raise MemoryError
        for l in range(Ml):
            kbesvec[jcusp][l]=<double*>sage_malloc(sizeof(double)*Ql) #Qv[jcusp][2])
            if kbesvec[jcusp][l]==NULL:
                raise MemoryError
    if verbose>0:
        printf("here0 \n")
        printf("Ml=%d \n",Ml)
        printf("Ql=%d \n",Ql)
    cdef double tmpr
    cdef double besprec
    besprec=1.0E-14
    ## Can we somehow use that "most" of Xpb[icusp][jcusp][2Q-j-1]=-Xpb[icusp][jcusp][j] ?
    ## uncomment the below lines to see this.
    # for jcusp from 0<=jcusp<nc:
    #     for icusp from 0<=icusp<nc:
    #         for j from 0<=j<Qv[1][1]:
    #             print "diff[",jcusp,icusp,j,"]=",Xpb[icusp][jcusp][j]+Xpb[icusp][jcusp][2*Qv[jcusp][1]-j-1]

    # for jcusp  in range(nc):
    #     for icusp  in range(nc):
    #         if icusp>0 and cusp_evs[icusp]<>0:
    #             continue
    #         for j in range(Qv[jcusp][2]):
    #             if Ypb[icusp][jcusp][j]==0:
    #                 continue
    #             for l in range(Mv[jcusp][2]):
    #                 lr=nvec[jcusp][l]*twopi
    #                 Mf = Mv[icusp][1]
    #                 besarg=fabs(lr)*Ypb[icusp][jcusp][j]
    #                 if lr<>0.0:
    #                     besselk_dp_c(&tmpr,R,besarg,besprec,1)
    #                     kbesvec[icusp][l][j]=sqrt(Ypb[icusp][jcusp][j])*tmpr
    #                 else:
    #                     kbesvec[icusp][l][j]=<double>1.0
    ## compute the number of chunks:
    cdef int chunksize = (Ml+1)/ncpus
    endpts[0][0]=0
    endpts[0][1]=chunksize
    for i in range(1,ncpus):
        endpts[i][0]=endpts[i-1][1]
        if i==ncpus-1:
            endpts[i][1]=Ml
        else:
            endpts[i][1]=endpts[i][0]+chunksize
        if endpts[i][1]>Ml:
            endpts[i][1]=Ml
    cdef double **kbesvec2
    kbesvec2 = <double**> sage_malloc(Ml*sizeof(double*))
    for icusp in range(nc):
        kbesvec2[icusp] = <double*> sage_malloc(Ml*sizeof(double*))
        for l in range(Ml):
            #printf("nvec(%d,%d)=%e \n",icusp,l,nvec[icusp][l])              
            if nvec[icusp][l]==0:
                if cuspidal==1:
                    kbesvec2[icusp][l]=0.0
                else:
                    kbesvec2[icusp][l]=1.0
            else:
                nrY2pi=fabs(nvec[icusp][l])*Y2pi
                besselk_dp_c(&kbes,R,nrY2pi,besprec,1)
                #printf("K(%d,%d)=%e \n",icusp,l,kbes)  
                kbes=sqrtY*kbes 
                kbesvec2[icusp][l]=kbes
                #printf("R=%e, arg=%e \n",R,nrY2pi)                  
    if verbose>0:
        for i in range(ncpus):
            printf("ranges[%d]=[%d:%d] \n",i,endpts[i][0],endpts[i][1])
    cdef int maxthreads = openmp.omp_get_max_threads()
    #printf("max threads: %d, ncpus: %d \n",maxthreads,ncpus)
#    openmp.omp_set_num_threads(ncpus)
    for i in prange(ncpus, nogil=True,num_threads=ncpus):
        #printf("i=%d \n",i)
        setV(V,ef1,ef2_c,kbesvec,kbesvec2,nvec,Mv,Qv,cusp_offsets,cusp_evs,
             Ypb,Xpb,Qfak,R,Y,sqrtY,nc,cuspidal,endpts[i][0],endpts[i][1],verbose)
    ###    
    ### Deallocate arrays
    ###
    if endpts<>NULL:
        for i in range(ncpus):
            sage_free(endpts[i])
        sage_free(endpts)
    if kbesvec<>NULL:
        for icusp in range(nc):
            if kbesvec[icusp]<>NULL:
                for l in range(Ml):
                    if kbesvec[icusp][l]<>NULL:
                        sage_free(kbesvec[icusp][l])
                sage_free(kbesvec[icusp])
        sage_free(kbesvec)
    if kbesvec2<>NULL:
        for icusp in range(nc):
            if kbesvec2[icusp]<>NULL:
                sage_free(kbesvec2[icusp])
        sage_free(kbesvec2)
    #print "deal kbbes1"
    if ef1<>NULL:
        for jcusp in range(nc):
            if ef1[jcusp]<>NULL:
                for icusp in range(nc):
                    if ef1[jcusp][icusp]<>NULL:
                        for n in range(Mv[icusp][2]):
                            if ef1[jcusp][icusp][n]<>NULL:
                                sage_free(ef1[jcusp][icusp][n])
                        sage_free(ef1[jcusp][icusp])
                sage_free(ef1[jcusp])
        sage_free(ef1)
    if ef2_c<>NULL:
        for icusp in range(nc):
            if ef2_c[icusp]<>NULL:
                for n in range(Mv[icusp][2]):
                    if ef2_c[icusp][n]<>NULL:
                        sage_free(ef2_c[icusp][n])
                sage_free(ef2_c[icusp])
        sage_free(ef2_c)
    if nvec<>NULL:
        for icusp in range(nc):
            if nvec[icusp]<>NULL:
                sage_free(nvec[icusp])
        sage_free(nvec)
    if cusp_offsets<>NULL:
        sage_free(cusp_offsets)

@cython.boundscheck(False)
@cython.cdivision(True)
cdef int compute_V_cplx_dp_par(double complex **V,
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
                           int ncpus=1,
                           int is_trivial=0):


    r"""
    Set up the matrix for the system of equations giving the Fourier coefficients of the Maass waveforms (assume no symmetry).
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
    - ``ncpus`` -- the number of cpus (cores)
    """
    cdef int l,i,j,icusp,jcusp,n,ni,lj,Ml,Ql,s,Qs,Qf,Mf,Ms
    cdef double pi,sqrtY,Y2pi,nrY2pi,argm,argpb,twopi,two,kbes,besarg,lr,nr
    cdef double complex ckbes,ctmpV,iargm,twopii,ctmp
    if not cuspidal in [0,1]:
        raise ValueError," parameter cuspidal must be 0 or 1"
    cdef int **endpts
    cdef int* cusp_offsets=NULL
    cdef double **nvec=NULL
    cdef double complex ****ef1=NULL
    cdef double complex ***ef2_c=NULL
    cdef double ***kbesvec=NULL
    ## Set constants
    pi=M_PI 
    sqrtY=sqrt(Y)
    two=<double>(2)
    Y2pi=Y*pi*two
    twopi=two*pi
    Ml=0; Ql=0
    if verbose>0:
        printf("in compute_V_cplx_dp_par \n")
    for i in range(nc):
        if Mv[i][2]>Ml:
            Ml=Mv[i][2]
        if Qv[i][2]>Ql:
            Ql=Qv[i][2]
        if verbose>0:
            printf("Qv[%d](%d,%d,%d)\n",i,Qv[i][0],Qv[i][1],Qv[i][2])
            printf("Mv[%d](%d,%d,%d)\n",i,Mv[i][0],Mv[i][1],Mv[i][2])
    if verbose>0:
        printf("N1=%d \n",N1)
        printf("Ql=%d \n",Ql)
    ## Allocate arrays 
    cusp_offsets=<int*>sage_malloc(sizeof(int)*nc)
    if cusp_offsets==NULL: raise MemoryError
    nvec = <double**>sage_malloc(sizeof(double*)*nc)
    if not nvec: raise MemoryError
    for icusp in range(nc):
        nvec[icusp] = <double*>sage_malloc(sizeof(double)*Ml)
    ef2_c = <double complex***>sage_malloc(sizeof(double complex**)*nc)
    if not ef2_c: raise MemoryError
    for icusp in range(nc):
        ef2_c[icusp] = <double complex**>sage_malloc(sizeof(double complex*)*Mv[icusp][2])
        for n in range(Mv[icusp][2]):
            ef2_c[icusp][n] = <double complex*>sage_malloc(sizeof(double complex)*Qv[icusp][2])
    ef1 = <double complex****>sage_malloc(sizeof(double complex***)*nc)
    if ef1==NULL: raise MemoryError
    for icusp in range(nc):
        ef1[icusp] = <double complex***>sage_malloc(sizeof(double complex**)*nc)
        if ef1[icusp]==NULL: raise MemoryError
        for jcusp in range(nc):
            ef1[icusp][jcusp] = <double complex**>sage_malloc(sizeof(double complex*)*Mv[jcusp][2])
            if ef1[icusp][jcusp]==NULL: raise MemoryError
            for n in range(Mv[jcusp][2]):
                ef1[icusp][jcusp][n] = <double complex*>sage_malloc(sizeof(double complex)*Qv[jcusp][2])
                if ef1[icusp][jcusp][n]==NULL: raise MemoryError
    endpts = <int**>sage_malloc(ncpus*sizeof(int*))
    if endpts==NULL:
        raise MemoryError
    for i in range(ncpus):
        endpts[i]=<int*>sage_malloc(2*sizeof(int))

    ## Assigning values to arrays
    for jcusp  in range(nc):
        cusp_offsets[jcusp]=0
        for icusp in range(jcusp):
            if icusp==0 or cusp_evs[icusp]==0:
                cusp_offsets[jcusp]+=Mv[icusp][2]
        if verbose>0:
            printf("cusp_offsets[%d]=%d \n",jcusp,cusp_offsets[jcusp])
    cdef int nc_sym=0
    for jcusp  in range(nc):
        if verbose>0:
            printf("cusp_evs[%d]=%d \n",jcusp,cusp_evs[jcusp])
        if jcusp==0 or cusp_evs[jcusp]<>0:
            nc_sym+=1
    for jcusp in range(nc):
        for n in range(Ml):
            nvec[jcusp][n]=<double>(n+Mv[jcusp][0])+alphas[jcusp]
    cdef int twoQm1
    twoQm1= 2*Qv[0][1]-1
    for jcusp in range(nc):
        #printf("symmetric cusps[%d]=%d \n",jcusp,symmetric_cusps[jcusp])        
        for n in range(Mv[jcusp][2]):
            for j in range(Qv[jcusp][2]):
                argm=nvec[jcusp][n]*Xm[j]
                if symmetric_cusps[jcusp]==0:
                    ef2_c[jcusp][n][j]=cos(argm)
                elif symmetric_cusps[jcusp]==1:
                    ef2_c[jcusp][n][j]=_Complex_I*sin(-argm)
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
                        ef1[icusp][jcusp][n][j]=_Complex_I*sin(argpb)
                    else:
                        ef1[icusp][jcusp][n][j]=cexpi(argpb)
                    ctmp = Cvec[icusp][jcusp][j]
                    ef1[icusp][jcusp][n][j]=ef1[icusp][jcusp][n][j]*ctmp

    if verbose>1:
        printf("here1121")
    cdef double besarg_old=0.0
    cdef double y,kbes_old=1.0
    kbesvec=<double***>sage_malloc(sizeof(double**)*nc)
    if kbesvec==NULL:
        raise MemoryError
    for jcusp  in range(nc):
        kbesvec[jcusp]=<double**>sage_malloc(sizeof(double*)*Ml)
        if kbesvec[jcusp]==NULL:
            raise MemoryError
        for l in range(Ml):
            kbesvec[jcusp][l]=<double*>sage_malloc(sizeof(double)*Ql) #Qv[jcusp][2])
            if kbesvec[jcusp][l]==NULL:
                raise MemoryError
    if verbose>0:
        printf("here0 \n")
        printf("Ml=%d \n",Ml)
        printf("Ql=%d \n",Ql)
    cdef double tmpr
    cdef double besprec
    besprec=1.0E-14
    ## Can we somehow use that "most" of Xpb[icusp][jcusp][2Q-j-1]=-Xpb[icusp][jcusp][j] ?
    ## uncomment the below lines to see this.
    # for jcusp from 0<=jcusp<nc:
    #     for icusp from 0<=icusp<nc:
    #         for j from 0<=j<Qv[1][1]:
    #             print "diff[",jcusp,icusp,j,"]=",Xpb[icusp][jcusp][j]+Xpb[icusp][jcusp][2*Qv[jcusp][1]-j-1]

    for jcusp  in range(nc):
        for icusp  in range(nc):
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
                        besselk_dp_c(&tmpr,R,besarg,besprec,1)
                        kbesvec[icusp][l][j]=sqrt(Ypb[icusp][jcusp][j])*tmpr
                    else:
                        kbesvec[icusp][l][j]=<double>1.0
    ## compute the number of chunks:
    cdef int chunksize = (Ml+1)/ncpus
    endpts[0][0]=0
    endpts[0][1]=chunksize
    for i in range(1,ncpus):
        endpts[i][0]=endpts[i-1][1]
        if i==ncpus-1:
            endpts[i][1]=Ml
        else:
            endpts[i][1]=endpts[i][0]+chunksize
        if endpts[i][1]>Ml:
            endpts[i][1]=Ml
    cdef double **kbesvec2
    kbesvec2 = <double**> sage_malloc(Ml*sizeof(double*))
    if verbose>0:
        for i in range(ncpus):
            printf("ranges[%d]=[%d:%d] \n",i,endpts[i][0],endpts[i][1])
    cdef int maxthreads = openmp.omp_get_max_threads()
    #printf("max threads: %d, ncpus: %d \n",maxthreads,ncpus)
#    openmp.omp_set_num_threads(ncpus)
    for i in prange(ncpus, nogil=True,num_threads=ncpus):
        #printf("i=%d \n",i)
        setV(V,ef1,ef2_c,kbesvec,kbesvec2,nvec,Mv,Qv,cusp_offsets,cusp_evs,
             Ypb,Xpb,Qfak,R,Y,sqrtY,nc,cuspidal,endpts[i][0],endpts[i][1],verbose)
    ###    
    ### Deallocate arrays
    ###
    if endpts<>NULL:
        for i in range(ncpus):
            sage_free(endpts[i])
        sage_free(endpts)
    if kbesvec<>NULL:
        for icusp in range(nc):
            if kbesvec[icusp]<>NULL:
                for l in range(Ml):
                    if kbesvec[icusp][l]<>NULL:
                        sage_free(kbesvec[icusp][l])
                sage_free(kbesvec[icusp])
        sage_free(kbesvec)
    if kbesvec2<>NULL:
        for icusp in range(nc):
            if kbesvec2[icusp]<>NULL:
                sage_free(kbesvec2[icusp])
        sage_free(kbesvec2)
    #print "deal kbbes1"
    if ef1<>NULL:
        for jcusp in range(nc):
            if ef1[jcusp]<>NULL:
                for icusp in range(nc):
                    if ef1[jcusp][icusp]<>NULL:
                        for n in range(Mv[icusp][2]):
                            if ef1[jcusp][icusp][n]<>NULL:
                                sage_free(ef1[jcusp][icusp][n])
                        sage_free(ef1[jcusp][icusp])
                sage_free(ef1[jcusp])
        sage_free(ef1)
    if ef2_c<>NULL:
        for icusp in range(nc):
            if ef2_c[icusp]<>NULL:
                for n in range(Mv[icusp][2]):
                    if ef2_c[icusp][n]<>NULL:
                        sage_free(ef2_c[icusp][n])
                sage_free(ef2_c[icusp])
        sage_free(ef2_c)
    if nvec<>NULL:
        for icusp in range(nc):
            if nvec[icusp]<>NULL:
                sage_free(nvec[icusp])
        sage_free(nvec)
    if cusp_offsets<>NULL:
        sage_free(cusp_offsets)
        
@cython.cdivision(True)
@cython.boundscheck(False)
cdef int setV(double complex **V,
              double complex ****ef1,
              double complex ***ef2_c,
              double ***kbesvec,
              double  **kbesvec2,
              double **nvec,
              int **Mv,
              int **Qv,
              int *cusp_offsets,
              double complex *cusp_evs,
              double ***Ypb,
              double ***Xpb,
              double *Qfak,
              double R,
              double Y,
              double sqrtY,
              int nc,
              int cuspidal,
              int la, int lb,
              int verbose=0
              ) nogil:

#               int ***CSvec, mpfr_t **** besv, mpfr_t *** Ypb, mpfr_t ****ef1cosv, mpfr_t ****ef1sinv, mpfr_t ***ef2cosv, mpfr_t ***ef2sinv, int nc, int Ql, int Ml, int l, mpfr_prec_t prec)n ogil:
    r"""
    Set a chunk of rows of the matrix
    """
    cdef int icusp,jcusp,lj,l,lj0,j,n,ni,N1
    cdef double lr,kbes,nr,twopi,besarg,tmpr,Y2pi
    cdef double besprec=1.0E-14
    cdef double complex ctmpV,cuspev,ckbes
    twopi=<double>6.2831853071795864769252867666
    cdef int id=cython.parallel.threadid()
    if verbose>0:
        printf("%d : setV in ranges:%d -- %d \n",id,la,lb)
    ## We also set the K-Bessel function here.
    Y2pi=Y*twopi
    for jcusp  in range(nc):
        for icusp  in range(nc):
            if icusp>0 and cusp_evs[icusp]<>0:
                continue
            for j in range(Qv[jcusp][2]):
                if Ypb[icusp][jcusp][j]==0:
                    continue
                for l in range(la,lb):
                    lr=nvec[jcusp][l]*twopi
                    Mf = Mv[icusp][1]
                    besarg=fabs(lr)*Ypb[icusp][jcusp][j]
                    if lr<>0.0:
                        besselk_dp_c(&tmpr,R,besarg,besprec,1)
                        kbesvec[icusp][l][j]=sqrt(Ypb[icusp][jcusp][j])*tmpr
                    else:
                        kbesvec[icusp][l][j]=<double>1.0
        for l in range(la,lb):
            if nvec[jcusp][l]==0:
                if cuspidal==1:
                    kbesvec2[jcusp][l]=0.0
                else:
                    kbesvec2[jcusp][l]=1.0
            else:
                nrY2pi=fabs(nvec[jcusp][l])*Y2pi
                besselk_dp_c(&kbes,R,nrY2pi,besprec,1)
                #printf("K(%d,%d)=%e \n",icusp,l,kbes)  
                #kbes=sqrtY*kbes 
                kbesvec2[icusp][l]=sqrtY*kbes
    for l in range(la,lb):
        for jcusp in range(nc):
            lr=nvec[jcusp][l]*twopi
            lj=cusp_offsets[jcusp]+l
            if jcusp>0 and cusp_evs[jcusp]<>0:
                lj0=l; cuspev=cusp_evs[jcusp]
            else:
                lj0=lj;  cuspev=1.0
            if lr==0.0 and cuspidal==1:
                continue
            for j in range(Qv[jcusp][2]):
                for icusp in range(0,nc):
                    if icusp>0 and cusp_evs[icusp]<>0:
                        continue
                    if Ypb[icusp][jcusp][j]==0:
                        continue
                    ckbes=kbesvec[icusp][l][j]*ef1[icusp][jcusp][l][j]
                    for n in range(0,Mv[icusp][2]):
                        if nvec[icusp][n]==0 and cuspidal==1:
                            continue
                        ni=cusp_offsets[icusp]+n
                        ctmpV=ckbes*ef2_c[icusp][n][j]
                        V[ni][lj0]=V[ni][lj0]+ctmpV*cuspev
    for icusp in range(nc):
        if icusp>0 and cusp_evs[icusp]<>0:
            continue
        for n in range(0,Mv[icusp][2]):
            ni=cusp_offsets[icusp]+n
            for jcusp in range(0,nc):
                if jcusp>0 and cusp_evs[jcusp]<>0:
                    continue
                for l in range(la,lb):
                    lj=cusp_offsets[jcusp]+l
                    V[ni][lj]=V[ni][lj]/<double complex>Qfak[jcusp]
    cdef int nthreads
    if verbose>0:
        nthreads = openmp.omp_get_num_threads()
        printf("%d : here1 in ranges:%d -- %d nthreads=%d\n",id,la,lb,nthreads)        
    
    for icusp in range(nc):
        if icusp>0 and cusp_evs[icusp]<>0:
            continue
        for n in range(la,lb):
            if nvec[icusp][n]==0.0 and cuspidal==1:
                continue
            ni=cusp_offsets[icusp]+n
            if verbose>2:
                printf("%d : K(%d,%d)=%e \n",id,icusp,n,kbesvec2[icusp][n])
            V[ni][ni]=V[ni][ni] - kbesvec2[icusp][n]            
    if verbose>0:
        printf("%d : here2 in ranges:%d -- %d nthreads=%d\n",id,la,lb,nthreads)
        
    return 0

cdef int setV_nosym(double complex **V,
              double complex ****ef1,
              double complex ***ef2_c,
              double ***kbesvec,
              double  **kbesvec2,
              double **nvec,
              int **Mv,
              int **Qv,
              int *cusp_offsets,
              double complex *cusp_evs,
              double ***Ypb,
              double ***Xpb,
              double *Qfak,
              double R,
              double Y,
               double sqrtY,                    
              int nc,
              int cuspidal,
              int la, int lb,
              int verbose=0
              ) nogil:

#               int ***CSvec, mpfr_t **** besv, mpfr_t *** Ypb, mpfr_t ****ef1cosv, mpfr_t ****ef1sinv, mpfr_t ***ef2cosv, mpfr_t ***ef2sinv, int nc, int Ql, int Ml, int l, mpfr_prec_t prec) nogil:
    r"""
    Set a chunk of rows of the matrix
    """
    cdef int icusp,jcusp,lj,l,lj0,j,n,ni,N1
    cdef double lr,kbes,nr,twopi
    cdef double complex ctmpV,cuspev,ckbes
    twopi=<double>6.2831853071795864769252867666
    cdef int id=cython.parallel.threadid()
    if verbose>0:
        printf("%d : setV in ranges:%d -- %d \n",id,la,lb)
    for jcusp in range(nc):
        for l in range(la,lb):
            lr=nvec[jcusp][l]*twopi
            lj=cusp_offsets[jcusp]+l
            if jcusp>0 and cusp_evs[jcusp]<>0:
                lj0=l; cuspev=cusp_evs[jcusp]
            else:
                lj0=lj;  cuspev=1.0
            if lr==0.0 and cuspidal==1:
                continue
            for j in range(Qv[jcusp][2]):
                for icusp in range(0,nc):
                    if icusp>0 and cusp_evs[icusp]<>0:
                        continue
                    if Ypb[icusp][jcusp][j]==0:
                        continue
                    ckbes=kbesvec[icusp][l][j]*ef1[icusp][jcusp][l][j]
                    for n in range(0,Mv[icusp][2]):
                        if nvec[icusp][n]==0 and cuspidal==1:
                            continue
                        ni=cusp_offsets[icusp]+n
                        ctmpV=ckbes*ef2_c[icusp][n][j]
                        V[ni][lj0]=V[ni][lj0]+ctmpV*cuspev
    for icusp in range(nc):
        if icusp>0 and cusp_evs[icusp]<>0:
            continue
        for n in range(0,Mv[icusp][2]):
            ni=cusp_offsets[icusp]+n
            for jcusp in range(0,nc):
                if jcusp>0 and cusp_evs[jcusp]<>0:
                    continue
                for l in range(la,lb):
                    lj=cusp_offsets[jcusp]+l
                    V[ni][lj]=V[ni][lj]/<double complex>Qfak[jcusp]
    cdef int nthreads = openmp.omp_get_num_threads()
    if verbose>0:
        printf("%d : here1 in ranges:%d -- %d nthreads=%d\n",id,la,lb,nthreads)        
    
    for icusp in range(nc):
        if icusp>0 and cusp_evs[icusp]<>0:
            continue
        for n in range(la,lb):
            if nvec[icusp][n]==0.0 and cuspidal==1:
                continue
            ni=cusp_offsets[icusp]+n
            if verbose>0:
                printf("%d : K(%d,%d)=%e \n",id,icusp,n,kbesvec2[icusp][n])
            V[ni][ni]=V[ni][ni] - kbesvec2[icusp][n]            
    if verbose>0:
        printf("%d : here2 in ranges:%d -- %d nthreads=%d\n",id,la,lb,nthreads)
        
    return 0


@cython.cdivision(True)
cdef int SMAT_cplx_par_dp(double complex** U,int N,int num_rhs,int num_set,double complex **C,double complex** values,int* setc,int ncpus):
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
    piv=<int*>sage_malloc(sizeof(int)*N)
    used=<int*>sage_malloc(sizeof(int)*N)
    if C==NULL:
        C=<double complex**>sage_malloc(sizeof(double complex*)*num_rhs)
        for j in range(num_rhs):
            C[j]=<double complex*>sage_malloc(sizeof(double complex*)*N)
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
        used[maxi]=1
        temp2=U[maxi][m]
        if cabs(temp2)==0.0:
            printf('ERROR: pivot(%d) == 0, system bad!!!',m)
            return 1
        #for j in range(m+1,N+num_rhs): # do j=M+1,N+1
        for j in prange(m+1,N+num_rhs,nogil=True,num_threads=ncpus): # do j=M+1,N+1        
            U[maxi][j]=U[maxi][j]/temp2
            #! eliminate from all other rows
        #for j in range(maxi):
        for j in prange(maxi,nogil=True,num_threads=ncpus): # do j=M+1,N+1
            TT=U[j][m]
            for k in range(m+1,N+num_rhs): #K=M+1,N+1
                U[j][k]=U[j][k]-U[maxi][k]*TT
        #for j in range(maxi+1,N): #DO J=Maxi+1,N
        for j in prange(maxi+1,N,nogil=True,num_threads=ncpus): #DO J=Maxi+1,N                
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
        sage_free(piv)
    if used<>NULL:
        sage_free(used)
    return 0
