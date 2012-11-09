# cython: profile=True
# -*- coding: utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2010 Fredrik Strömberg <fredrik314@gmail.com>
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
Class containing the linear system we use to solve for Fourier coefficients.
AUTHOR:

 - Fredrik Stroemberg


"""
#from sage.rings.complex_mpc import MPComplexField
#from psage.matrix import *
cdef mpc_rnd_t rnd
cdef mpfr_rnd_t rnd_re
rnd = MPC_RNDNN
rnd_re = GMP_RNDN
include "sage/ext/interrupt.pxi"  # ctrl-c interrupt block support
include "sage/ext/stdsage.pxi"  # ctrl-c interrupt block support
include "sage/ext/cdefs.pxi"
include "sage/ext/gmp.pxi"
from sage.rings.complex_mpc cimport MPComplexNumber
from psage.matrix.matrix_complex_dense cimport Matrix_complex_dense
from sage.all import copy,MatrixSpace
from psage.modform.maass.vv_harmonic_weak_maass_forms_alg import pullback_pts_vv_mpc
from sage.rings.real_mpfr import RealField
from sage.rings.complex_mpc import MPComplexField
#cdef class LinearSystem(object):
import cython
cdef class LinearSystem(object):
    r"""
    Linear system.
    """
    def __cinit__(self,X,int prec=53,int verbose=0,**kwds):
       print "cinit!"
       self._verbose = verbose
       self._M0=0
       self._prec = prec
       self._weight = Rational(X.weight())

       
    def __init__(self,X,int prec=53,int verbose=0,**kwds):
        print "init"
        self._multiplier = X
        self._CF = MPComplexField(prec)
        self._RF = RealField(prec)


        #self._matrix = None #Matrix_complex_dense(M,0)

    def __dealloc__(self):
        pass
        
    def multiplier(self):
        return self._multiplier

    cpdef setup_system_vv(self,int M0,RealNumber Y):
        r"""
        Set up a symmetrized system for the Weil representation for holomorphic functions
        SLOW version. Used for testing.
        
        """
        assert hasattr(self._multiplier,"D")
        #if prec>0:
        #    self._CF = MPComplexField(prec)
        #    self._RF = RealField(prec)
        self._M0 = M0
        #RF = self._RF
        s = len(self.multiplier().D())*(M0+1)

        MS = MatrixSpace(self._CF,s,s)
        self._matrix = Matrix_complex_dense(MS,0)
        D = self._multiplier.D()
        ### Get Pullback points
        Q = M0+15
        Ml = M0+1
        X = self._multiplier;
        weight=self._RF(self._weight)
        Ds = D[0]; Df=D[-1]
        st = X._sym_type
        Xm,Zpb,Cvec,Cfak_m,Cfak_p=pullback_pts_vv_mpc(X,Q,Y,weight,self._verbose)     
        
        ### Now assign to the matrix
        for l in range(s):
            for k in range(s):
                self._matrix[l,k]=0            
        for l in range(M0+1):
            for b in D:
                mbi = X.neg_index(b)
            if self._verbose>1:
                print "b=",b
                print "-b=",mbi
            for l in range(M0+1):
                ll = self._RF(l)+self._RF(X.Qv[b])
                lj = (b-Ds)*Ml+l
                for a in D:
                    mai = X.neg_index(a)
                    if self._verbose>1:
                        print "a=",a
                        print "-a=",mai
                    for j in range(1-Q,Q+1):
                        #f1 = CF(0,ll*Xpb[j]).exp()                    
                        lx = ll*Zpb[j].real()
                        kk = self._RF(-ll*Zpb[j].imag()).exp()
                        for n in range(M0+1):
                            nn = self._RF(n)+self._RF(X.Qv[a])
                            f2 = self._CF(0,lx-nn*Xm[j]).exp()             
                            if b<>mbi:
                                ch = Cvec[j][a,b]+st*Cvec[j][a,mbi]
                            elif st==1:
                                ch = Cvec[j][a,b]
                            else:
                                ch = 0
                            if ch<>0:
                                if j<=0:
                                    ch = ch*Cfak_m[1-j]
                                else:
                                    ch = ch*Cfak_p[j]
                                tmp = kk*f2*ch
                                ni = (a-Ds)*Ml+n
                                self._matrix[ni,lj]=self._matrix[ni,lj]+tmp
        twoQ=self._RF(2*Q)
        if self._verbose>0:
            print "V1[0,0]=",self._matrix[0,0]
        for n in range(s):
            for l in range(s):
                self._matrix[n,l]=self._matrix[n,l]/twoQ
        if self._verbose>0:
            print "V2[0,0]=",self._matrix[0,0]
        twopi=self._RF(2)*self._RF.pi()
        for a in D:
            #mai = X.neg_index(a)
            for n in range(M0+1):
                nn=self._RF(n)+self._RF(X.Qv[a])
                nn = nn*twopi
                ni = (a-Ds)*Ml+n
                tmp = self._RF(-nn*Y).exp()
                self._matrix[ni,ni] = self._matrix[ni,ni]-tmp            

    #def setup_system_sv(
        

                
    def simple_solve(self,cuspidal=1):
        ## Check which components we have to set to zero
        X = self._multiplier; weight=X._weight
        D = X.D()
        zero_set=[]; one_set=[]
        for a in D:
            if X.Qv[a]==0 and cuspidal==1:
                zero_set.append((self._M0+1)*a)
            else:
                one_set.append((self._M0+1)*a)
        zero_set.reverse()
        A = copy(self._matrix)
        if self._verbose>0:
            print "zero_set=",zero_set
            print "one_set=",one_set
        if one_set==[]:
            raise ValueError
        for j in zero_set:
            A.delete_row(j)
        A.delete_row(one_set[0])
        B = -A.column(one_set[0])
        A.delete_column(one_set[0])
        for j in zero_set:
            A.delete_column(j)        
        AA = A.inverse()
        if self._verbose>0:
            print "rows:",A.nrows()
            print "cols:",A.ncols()
            print "A.prec()=",A.prec()
            print "B.prec()=",B.prec()
            print "AA.prec()=",AA.prec()
        
        C = copy(AA*B)
        return C

                            
### Efficient solving routine
@cython.cdivision(True)
cdef SMAT_mpc(mpc_t** U,int N,int num_rhs,int num_set,mpc_t** C,mpc_t** values,int* setc):
    r"""
    Use Gauss elimination to solve a linear system AX=B
    U = (A|B) is a N x (N+num_rhs) double complex matrix 
    setc and values should be allocated of length num_set
    """
    cdef int m,maxi,j,k,i
    cdef mpc_t ctemp
    cdef int prec
    prec = mpc_get_prec(U[0][0])
    #mpc_init2(TT,prec);
    mpc_init2(ctemp,prec)    
    #cdef double complex **tmpu
    cdef mpfr_t rtemp,tabs
    mpfr_init2(rtemp,prec); mpfr_init2(tabs,prec)
    cdef int *piv
    cdef int *used
    piv=<int*>sage_malloc(sizeof(int)*N)
    used=<int*>sage_malloc(sizeof(int)*N)
    if C==NULL:
        C=<mpc_t**>sage_malloc(sizeof(mpc_t*)*num_rhs)
        for j from 0<=j<num_rhs:
            C[j]=<mpc_t*>sage_malloc(sizeof(mpc_t*)*N)
    for j in range(N):
        piv[j]=0
        used[j]=0
    for m in range(N):
        mpfr_set_ui(rtemp,0,rnd_re)
        maxi=0
        #Locate maximum
        for j in range(N):
            if used[j]<>0:
                continue
            #print "U[",j,m,"]=",U[j][m],abs(U[j][m])
            mpc_abs(tabs,U[j][m],rnd_re) 
            #if cabs(U[j][m]) <= temp:
            if mpfr_cmp(tabs,rtemp)<=0:
                continue
            maxi=j
            mpfr_set(rtemp,tabs,rnd_re)
            #temp=cabs(U[j][m])
            #print "temp=",temp
        piv[m]=maxi
        #print "piv[",m,"]=",maxi
        #print "norm=",temp
        #return
        used[maxi]=1
        mpc_set(ctemp,U[maxi][m],rnd)
        mpc_abs(rtemp,ctemp,rnd_re)
        if mpfr_zero_p(rtemp)<>0: #==0.0:
            print 'ERROR: pivot(',m,') == 0, system bad!!!'
            raise ArithmeticError
        for j in range(m+1,N+num_rhs): # do j=M+1,N+1
            mpc_div(U[maxi][j],U[maxi][j],ctemp,rnd)
            #! eliminate from all other rows
        for j in range(maxi):
            #TT=U[j][m]
            for k in range(m+1,N+num_rhs): #K=M+1,N+1
                mpc_mul(ctemp,U[maxi][k],U[j][m],rnd)
                mpc_sub(U[j][k],U[j][k],ctemp,rnd)
                #U[j][k]=U[j][k]-U[maxi][k]*TT
        for j in range(maxi+1,N): #DO J=Maxi+1,N
            #TT=U[j][m]
            for k in range(m+1,N+num_rhs): #do K=M+1,N+1
                mpc_mul(ctemp,U[maxi][k],U[j][m],rnd)
                mpc_sub(U[j][k],U[j][k],ctemp,rnd)
                #U[j][k]=U[j][k]-U[maxi][k]*TT
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
                mpc_set(C[i][m],values[i][j],rnd)
                #continue
            else:
                mpc_set(C[i][m],U[piv[m-m_offs]][N+i],rnd)
        #print "C0[",m,"]=",C[m]
    if piv<>NULL:
        sage_free(piv)
    if used<>NULL:
        sage_free(used)
    mpc_clear(ctemp)
    mpfr_clear(tabs); mpfr_clear(rtemp)

cpdef test_lin_solve(Matrix_complex_dense A,Matrix_complex_dense RHS,dict setC):
    cdef mpc_t** U=NULL,**C=NULL,**values=NULL
    cdef int* setc
    cdef int m,n,nr,nc
    cdef int nrow,ncol,nrhs,nset,prec
    cdef MPComplexNumber tmpc    
    prec = A.prec()
    CF = MPComplexField(prec)
    tmpc = CF(1)
    nrow = A.nrows()
    ncol = A.ncols()
    nrhs = RHS.ncols()
    nset = len(setC[0].keys())
    if RHS.nrows()<>nrow:
        raise ArithmeticError,"Incompatible RHS and LHS!"
    setc = <int*>sage_malloc(sizeof(int)*nset)    
    U = <mpc_t **>sage_malloc(sizeof(mpc_t*)*nrow)
    C = <mpc_t **>sage_malloc(sizeof(mpc_t*)*nrow)
    values = <mpc_t **>sage_malloc(sizeof(mpc_t*)*nrhs)
    for i in range(nrow):
        C[i]=<mpc_t*>sage_malloc(sizeof(mpc_t)*(nset))
        for j in range(nset):
            mpc_init2(C[i][j],prec)
        U[i]=<mpc_t*>sage_malloc(sizeof(mpc_t)*(ncol+nrhs))
        for j in range(ncol):
            mpc_init2(U[i][j],prec)
            mpc_set(U[i][j],A._matrix[i][j],rnd)
        for j in range(ncol,ncol+nrhs):
            mpc_init2(U[i][j],prec)
            mpc_set(U[i][j],RHS._matrix[i][j],rnd)

    for i in range(nrhs):
        values[i] = <mpc_t *>sage_malloc(sizeof(mpc_t)*nset)
        for j in range(nset):
            c,n = setC[i].keys()[j]
            setc[j]=n                
            tmpc = CF(setC[i]((c,n)))
            mpc_set(values[i][j],tmpc.value,rnd)
    #    SMAT_mpc(U,N,num_rhs,num_set,C,values,setc)
    if U<>NULL:
        for i in range(nrow):
            if U[i]<>NULL:
                for j in range(ncol+nrhs):
                    mpc_clear(U[i][j])
                sage_free(U[i])
        sage_free(U)
        

#cdef SMAT_mpc(mpc_t** U,int N,int num_rhs,int num_set,mpc_t** C,mpc_t** values,int* setc):
