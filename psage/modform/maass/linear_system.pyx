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

include "sage/ext/interrupt.pxi"  # ctrl-c interrupt block support
include "sage/ext/stdsage.pxi"  # ctrl-c interrupt block support
include "sage/ext/cdefs.pxi"
include "sage/ext/gmp.pxi"
from psage.matrix.matrix_complex_dense cimport Matrix_complex_dense
from sage.all import copy,MatrixSpace
from psage.modform.maass.vv_harmonic_weak_maass_forms_alg import pullback_pts_vv_mpc
from sage.rings.real_mpfr import RealField
from sage.rings.complex_mpc import MPComplexField
#cdef class LinearSystem(object):
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

                            
