# cython: profile=True
#clib mpc gmp
#include "/home/fredrik/install/sage/devel/sage-mpc_test/sage/ext/interrupt.pxi"
#import cython



include '../modform/maass/common_defs.pxd'


#include "sage/rings/mpc.pxi"
from psage.rings.mpfr_nogil cimport *
from sage.all import save,incomplete_gamma,load
import mpmath    
import cython
from libc.stdlib cimport abort, malloc, free
from sage.libs.flint.fmpz_mat cimport *
from sage.libs.flint.fmpz cimport *

from cython.parallel cimport parallel, prange

#from sage.libs.mpfr cimport *
import sys

from sage.rings.complex_mpc import MPComplexField
from sage.rings.real_mpfr import RealField

from sage.rings.complex_field import ComplexField
from sage.rings.integer import Integer,is_Integer
from sage.all import exp,I,CC,vector
from sage.functions.transcendental import zeta
from sage.all import find_root
from psage.rings.mpc_extras cimport _mpc_mul,_mpc_mul_fr,_mpc_add,_mpc_set
## Since I want to import sinand cos from math.h 
#from sage.all import sin as ssin
#from sage.all import cos as scos
#from libc.math cimport sin as dsin
#from libc.math cimport cos as dcos
import numpy as np
cimport numpy as cnp
#cimport cython


import cython
cdef extern from "math.h":
    int iabs(int)
    double fabs(double)
    double fmax(double,double)
    int ceil(double) 
    int floor(double)
    double sqrt(double)
    double sin(double)
    double cos(double)
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



@cython.cdivision(True)
cdef int gcd( int a, int b ) nogil:
  cdef int c
  while a <>0:
      c = a; a = b % a;  b = c
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

    
# some things that are not in the sage.libs.mpfr
#cdef extern from "mpfr.h":
#    int mpfr_mul_d (mpfr_t, mpfr_t, double, mp_rnd_t) 
   
from sage.matrix.matrix_dense cimport *
#from psage.rings.mpc_extras cimport *
from psage.modules.vector_complex_dense cimport Vector_complex_dense
from psage.modules.vector_real_mpfr_dense cimport Vector_real_mpfr_dense
from psage.matrix.matrix_complex_dense cimport Matrix_complex_dense
from sage.functions.all import ceil as pceil
from sage.matrix.all import MatrixSpace

from sage.all import  Parent,RR,ZZ,QQ,is_even,is_odd
from sage.rings.complex_mpc import _mpfr_rounding_modes,_mpc_rounding_modes
from sage.all import zeta
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense

from sage.rings.complex_mpc cimport MPComplexNumber

cpdef pochammer_over_fak(MPComplexNumber z,int n):
    cdef MPComplexNumber res
    cdef mpc_t resv
    cdef int prec = z.parent().prec()
    mpc_init2(resv,prec)
    res = MPComplexField(prec)(0)
    
    _pochammer_over_fak(&resv,z.value,n,rnd,rnd_re)
    mpc_set(res.value,resv,rnd)
    mpc_clear(resv)
    return res

cdef inline void _pochammer_over_fak(mpc_t *res, mpc_t z, int n,mpc_rnd_t rnd,mpfr_rnd_t rnd_re) nogil:
    r""" Compute the Pochammer symbol (s)_{n} / n!   """
    cdef int j
    cdef mpc_t x,zz
    cdef mpc_t t[2]
    cdef int prec= mpc_get_prec(z) +50 
    mpc_init2(t[0],prec)
    mpc_init2(t[1],prec)
    mpc_init2(x,prec)
    mpc_init2(zz,prec)
    mpc_set_ui(res[0],1,rnd)
    mpc_set(zz,z,rnd)
    for j in range(n):
        #mpc_add_ui(z.re,z.re,1,rnd_re)
        mpc_add_ui(x,zz,j,rnd)
        mpc_div_ui(x,x,j+1,rnd)
        _mpc_mul(res,res[0],x,t,rnd,rnd_re)
        #mpc_mul(res[0],res[0],x,rnd)
    mpc_clear(t[0])
    mpc_clear(t[1])
    mpc_clear(x)
    mpc_clear(zz)




cpdef trop_approximation(Matrix_complex_dense A, Matrix_integer_dense B, int M0, int q,int dim, int Nmax, MPComplexNumber s, RealNumber llambda, int verbose=0,int approx_type=0,int eps=1, dict alphas={},dict rhos={}):
    r"""
    Compute approximation of transfer operator.
    """
    cdef mpc_t z
    cdef int n = dim*(M0+1)
    F = A.parent().base_ring()
    cdef mpc_rnd_t rnd = _mpc_rounding_modes.index(F.rounding_mode())
    cdef mpfr_rnd_t rnd_re = _mpfr_rounding_modes.index(F.rounding_mode_real())
    cdef mpfr_t llv
    cdef mpfr_t *alphas_t,*rhos_t
    cdef RealNumber x
    cdef int prec = F.prec()
    mpc_init2(z,prec)
    mpfr_init2(llv,prec)
    mpc_set(z,s.value,rnd)
    mpfr_set(llv,llambda.value,rnd_re)
    if verbose>2:
        for i in range(dim):
            for j in range(dim):
                print "B[{0}][{1}]={2}".format(i,j,B[i][j])
    #Nij =  <mpc_t *> sage_malloc(sizeof(*int)*(dim))
    #print "in trop_approx 1!"
    if approx_type==0:
        setup_approximation(A._matrix, M0,  B._matrix, dim,  q, Nmax, z,llv,rnd,rnd_re,verbose)
    elif approx_type==3:
        RF = RealField(prec)

        assert len(alphas)==dim
        assert len(rhos)==dim
        alphas_t=<mpfr_t*>sage_malloc(sizeof(mpfr_t)*dim)
        rhos_t=<mpfr_t*>sage_malloc(sizeof(mpfr_t)*dim)
        for i in range(dim):
            mpfr_init2(alphas_t[i],prec)
            mpfr_init2(rhos_t[i],prec)
            x = RF(alphas[i])
            mpfr_set(alphas_t[i],x.value,rnd_re)
            x = RF(rhos[i])
            mpfr_set(rhos_t[i],x.value,rnd_re)
        setup_approximation_sym(A._matrix, M0,  B._matrix, dim,  q,eps, Nmax, z,llv,rnd,rnd_re,alphas_t,rhos_t,verbose)
    mpc_clear(z)
    mpfr_clear(llv)



cdef setup_approximation(mpc_t** A, int M,  fmpz_mat_t Nij, int dim, int q, int Nmax, mpc_t s, mpfr_t llambda, mpc_rnd_t rnd,mpfr_rnd_t rnd_re,int verbose=0):
    r""" Setup the Matrix approximation of the transfer operator.  """
    cdef int i,j,k,l,n,nn
    cdef mpc_t twos,tmp,tmp1,clambda,AA,poc
    cdef mpz_t tmp_mpz 
    cdef int prec = mpc_get_prec(s)
    cdef mpc_t t[4]
    cdef mpc_t *Z=NULL,*lpow=NULL,*twozn=NULL
    cdef mpc_t **npow=NULL,***lsum=NULL
    cdef mpfr_t **B
    cdef int **abN
    cdef MPComplexNumber z,twoz,tmpz
    cdef mpc_t mpc_zv

    sig_on()
    #print "in setup_approx!"
    mpc_init2(mpc_zv,prec)
    mpc_init2(AA,prec)
    mpc_init2(poc,prec)
    mpc_init2(tmp,prec); mpc_init2(tmp1,prec)
    mpc_init2(twos,prec)
    mpc_init2(clambda,prec)
    mpc_init2(t[0],prec); mpc_init2(t[1],prec)
    mpc_init2(t[2],prec); mpc_init2(t[3],prec)
    mpc_set_fr(clambda,llambda,rnd)
    mpc_set(twos,s,rnd)
    mpc_mul_ui(twos,twos,2,rnd)
    mpz_init(tmp_mpz)
    #mpc_mul_ui(twoz,2,rnd)
    CF = MPComplexField(prec)
    tmpz=CF(0)
    z = CF(0)
    twoz = MPComplexField(prec+20)(0)
    mpc_set(twoz.value,twos,rnd)
    Z  =  <mpc_t *> sage_malloc(sizeof(mpc_t)*(2*M+1))
    if Z==NULL: raise MemoryError
    twozn  =  <mpc_t *> sage_malloc(sizeof(mpc_t)*(2*M+1))
    if twozn==NULL: raise MemoryError
    # Try to evaluate all powers outside the double loop
    lsum=<mpc_t ***> sage_malloc(sizeof(mpc_t**)*(2*M+1))
    if lsum==NULL: raise MemoryError
    if q > 3:
        lpow = <mpc_t *> sage_malloc(sizeof(mpc_t)*(2*M+1))
        if lpow==NULL: raise MemoryError
        npow = <mpc_t **> sage_malloc(sizeof(mpc_t*)*(Nmax-1))
        if npow==NULL: raise MemoryError
        for i in range(Nmax):
            npow[i] = NULL
            npow[i] = <mpc_t*> sage_malloc(sizeof(mpc_t)*(2*M+1))
            if npow[i]==NULL: raise MemoryError
            for n in range(2*M+1):
                mpc_init2(npow[i][n],prec)
                mpc_set_ui(npow[i][n],0,rnd)
                
    abN = <int **>sage_malloc(sizeof(int*)*(dim))
    if abN==NULL: raise MemoryError
    for i in range(dim):
        abN[i]==NULL
        abN[i] = <int *> sage_malloc(sizeof(int)*(dim))
        if abN[i]==NULL: raise MemoryError
        for j in range(dim):
            fmpz_get_mpz(tmp_mpz,fmpz_mat_entry(Nij,i,j))
            abN[i][j]=mpz_get_si(tmp_mpz)
            if verbose>0:
                print "N[{0}][{1}]={2}".format(i,j,abN[i][j])
            abN[i][j]=abs(abN[i][j])

    for n in range(2*M+1):
        lsum[n]=<mpc_t **> sage_malloc(sizeof(mpc_t*)*(dim))
        mpc_init2(twozn[n],prec)
        mpc_add_ui(twozn[n],twos,n,rnd)
        #_mpc_set(&twozn[i],tmpprec)
#        _mpc_set(&twoz.value,twozn[n],rnd_re)
        mpc_set(twoz.value,twozn[n],rnd)
        z = twoz.zeta()
        mpc_init2(Z[n],prec)
        mpc_set(Z[n],z.value,rnd)
        mpc_neg(tmp,twozn[n],rnd)        
        if q>3:
            mpc_init2(lpow[n],prec)
            mpc_pow(lpow[n],clambda,tmp,rnd)
            # lpow[n] = lambda**(-2s-n)
            # npow[n][i]=(i+2)**(-2s-n) 
            for i from 2<= i <= Nmax: 
                mpfr_set_ui(tmp1.re,i,rnd_re)
                mpfr_set_ui(tmp1.im,0,rnd_re)
                mpc_pow(npow[i-2][n],tmp1,tmp,rnd)
        for i from 0<= i < dim:
            lsum[n][i]=<mpc_t *> sage_malloc(sizeof(mpc_t)*(dim))
            for j from 0<= j < dim:
                mpc_init2(lsum[n][i][j],prec)
                mpc_set_ui(lsum[n][i][j],0,rnd)
                for l from 1 <= l < abN[i][j]:
                    mpfr_set_ui(tmp1.re,l,rnd_re)
                    mpfr_set_ui(tmp1.im,0,rnd_re)
                    mpc_pow(tmp1,tmp1,tmp,rnd)
                    _mpc_add(&lsum[n][i][j],lsum[n][i][j],tmp1,rnd_re)
               
    #print "N=",abN[0][0]
    cdef int ni,kj
    for n in range(M+1): 
        for k in range(M+1): 
            _pochammer_over_fak(&poc,twozn[k],n,rnd,rnd_re)
            #poc=_pochammer(twoz,n)/RF.factorial(n)
            mpfr_add_ui(tmp.re,tmp.re,n,rnd_re) # tmp = 2s+n+k
              #    #poc=poc*lambda**(-n-k-2s)
            if q<>3:
                _mpc_mul(&poc,poc,lpow[n+k],t,rnd,rnd_re)
            if (n + k) % 2 == 1:
                sg = -1
            else:
                sg = 1
            for i in range(dim):
                for j in range(dim):
                    if abN[i][j]==0:
                        mpc_set_ui(AA,0,rnd)
                    else:
                        if 2*j < dim -2 or 2*j>dim:
                            # never happens for q = 3
                            nn = <int>abN[i][j]
                            if nn==0:
                                mpc_set_ui(AA,0,rnd)
                            elif nn==1:
                                mpc_set_ui(AA,1,rnd)
                            else:
                                _mpc_set(&AA,npow[nn-2][n+k],rnd_re)
                        else:
                            mpc_set(AA,Z[k+n],rnd)
                            _mpc_sub(&AA,AA,lsum[k+n][i][j],rnd_re)
                        fmpz_get_mpz(tmp_mpz,fmpz_mat_entry(Nij,i,j))                        
                        if mpz_sgn(tmp_mpz)>0 and sg==-1:
                            #AA=-AA
                            mpfr_neg(AA.re,AA.re,rnd_re)
                            mpfr_neg(AA.im,AA.im,rnd_re)
                    _mpc_mul(&tmp1,poc,AA,t,rnd,rnd_re)
                    ni = i*(M+1)+n
                    kj = j*(M+1)+k
                    mpc_set(A[ni][kj],tmp1,rnd)
                    if ni==0 and kj==7 and verbose>0:
                        mpc_set(tmpz.value,AA,rnd)
                        print "AA=",tmpz
                        print "sg=",sg
                        print "B=",<int>abN[i][j]
                        mpc_set(tmpz.value,poc,rnd)
                        print "poc=",tmpz
                        mpc_set(tmpz.value,A[ni][kj],rnd)
                        print "A=",tmpz
    mpc_clear(t[0]);mpc_clear(t[1])
    mpc_clear(t[2]);mpc_clear(t[3])
    mpc_clear(twos); mpc_clear(clambda)
    mpc_clear(tmp); mpc_clear(tmp1)
    mpc_clear(AA);  mpc_clear(poc)

    if Z<>NULL:
        for n from 0 <= n < 2*M+1:
            mpc_clear(Z[n])
        sage_free(Z)
    if lsum<>NULL:
        for n from 0 <= n < 2*M+1:
            if lsum[n]<>NULL:
                for i from 0<= i < dim:
                    if lsum[n][i]<>NULL:
                        for j from 0<= j < dim:
                            mpc_clear(lsum[n][i][j])
                        sage_free(lsum[n][i])
                sage_free(lsum[n])
        sage_free(lsum)
    if abN<>NULL:
        for i in range(dim):
            if abN[i]<>NULL:
                sage_free(abN[i])
        sage_free(abN)    
    if twozn<>NULL:
        for n in range(2*M+1):
            mpc_clear(twozn[n])
        sage_free(twozn)
    if q>3:
        if lpow<>NULL:
            for n in range(2*M+1):
                mpc_clear(lpow[n])
            sage_free(lpow)
        if npow<>NULL:
            for i in range(Nmax):
                if npow[i]<>NULL:
                    for n in range(2*M+1):
                        mpc_clear(npow[i][n])
                    sage_free(npow[i])
            sage_free(npow)
    sig_off()
    #return




cdef setup_approximation_sym(mpc_t** A, int M,  fmpz_mat_t Nij, int dim, int q, int eps, int Nmax, mpc_t s, mpfr_t llambda, mpc_rnd_t rnd,mpfr_rnd_t rnd_re,mpfr_t* alphas,mpfr_t* rhos,verbose=0):
    r"""    Setup the Matrix approximation of the transfer operator (version 3).  """
    cdef int i,j,k,l,n,nn
    cdef mpc_t twos,tmp,tmp1,clambda,poc, *AA
    cdef int prec = mpc_get_prec(s)
    cdef mpc_t t[4]
    cdef mpc_t ****Z=NULL,*lpow=NULL,*twozn=NULL
    cdef mpc_t*** dummy
    cdef mpfr_t **ajpow
    cdef Vector_real_mpfr_dense B
    cdef int **abN
    cdef int sym_dim = dim / 2
    cdef MPComplexNumber twoz,tmpz,z
    cdef ComplexNumber twozz
    cdef mpz_t binc
    cdef mpfr_t fak1,fak2,fak
    cdef RealNumber tmpx
    cdef mpz_t tmp_mpz
    sig_on()
    RF = RealField(prec)
    CF = MPComplexField(prec)
    CFF = ComplexField(prec)
    B = Vector_real_mpfr_dense(vector(RF,2).parent(),0)
    twozz = CFF(0)
    tmpx = RF(0)
    mpmath.mp.prec=prec
    mpz_init(tmp_mpz)
    z = CF(0)
    assert abs(eps)==1
    if verbose>0:
        print "eps=",eps
        print "prec=",prec
    mpfr_init2(fak1,prec);mpfr_init2(fak2,prec);mpfr_init2(fak,prec)
    #B[0]=RF(0); B[1]=RF(0)
    mpz_init(binc)
    #print "in setup_approx!"
    AA =  <mpc_t *> sage_malloc(sizeof(mpc_t)*(2))
    mpc_init2(AA[0],prec);mpc_init2(AA[1],prec)
    mpc_init2(poc,prec)
    mpc_init2(tmp,prec); mpc_init2(tmp1,prec)
    mpc_init2(twos,prec)
    mpc_init2(clambda,prec)
    mpc_init2(t[0],prec); mpc_init2(t[1],prec)
    mpc_init2(t[2],prec); mpc_init2(t[3],prec)
    mpc_set_fr(clambda,llambda,rnd)
    mpc_set(twos,s,rnd)
    mpc_mul_ui(twos,twos,2,rnd)
    tmpz=CF(0)
    #z = CF(0)
    twoz = MPComplexField(prec+20)(0)
    mpc_set(twoz.value,twos,rnd)
    ## Recall that sizeof(mpc_t ***) does not work
    Z  =  <mpc_t ****> sage_malloc(sizeof(dummy)*(sym_dim))
    if Z==NULL: raise MemoryError
    for i in range(sym_dim):
        Z[i] = <mpc_t ***> sage_malloc(sizeof(mpc_t**)*(sym_dim))
        if Z[i]==NULL: raise MemoryError
        for j in range(sym_dim):
            Z[i][j]  =  <mpc_t **> sage_malloc(sizeof(mpc_t*)*(2*M+1))
            if Z[i][j]==NULL: raise MemoryError
            for n in range(2*M+1):
                Z[i][j][n]  =  <mpc_t *> sage_malloc(sizeof(mpc_t)*(2))
                if Z[i][j][n]==NULL: raise MemoryError
                mpc_init2(Z[i][j][n][0],prec)
                mpc_init2(Z[i][j][n][1],prec)
    twozn  =  <mpc_t *> sage_malloc(sizeof(mpc_t)*(2*M+1))
    if twozn==NULL: raise MemoryError
    # Try to evaluate all powers outside the double loop
    ajpow=<mpfr_t **> sage_malloc(sizeof(mpfr_t*)*(sym_dim))
    if ajpow==NULL: raise MemoryError
    for i in range(sym_dim):        
        ajpow[i]=<mpfr_t *> sage_malloc(sizeof(mpfr_t)*(M+1))
        if ajpow[i]==NULL: raise MemoryError
        for n in range(M+1):
            mpfr_init2(ajpow[i][n],prec)
            mpfr_pow_si(ajpow[i][n],alphas[i],n,rnd_re)

    if q > 3:
        lpow = <mpc_t *> sage_malloc(sizeof(mpc_t)*(2*M+1))
        if lpow==NULL: raise MemoryError
    cdef RealNumber xarg
    cdef mpc_t zarg,minus_twoz
    mpc_init2(zarg,prec)
    mpc_init2(minus_twoz,prec)    
    xarg = RF(0)
    for n in range(2*M+1):
        #lsum[n]=<mpc_t **> sage_malloc(sizeof(mpc_t*)*(dim))
        mpc_init2(twozn[n],prec)
        mpc_add_ui(twozn[n],twos,n,rnd)
        #_mpc_set(&twoz.value,twozn[n],rnd_re)
        mpc_set(twoz.value,twozn[n],rnd)
        mpc_neg(minus_twoz,twoz.value,rnd)
        twozz=CFF(twoz.real(),twoz.imag())
        for i in range(sym_dim):
            for j in range(sym_dim):
                fmpz_get_mpz(tmp_mpz,fmpz_mat_entry(Nij,i,j))
                mpfr_set_z(B._entries[0],tmp_mpz,rnd_re)
                fmpz_get_mpz(tmp_mpz,fmpz_mat_entry(Nij,i,dim-1-j))
                #mpfr_set_z(B._entries[1],Nij[i][dim-1-j],rnd_re)
                mpfr_set_z(B._entries[1],tmp_mpz,rnd_re)
                #B[0] = RF(self._Nij[i,j])
                #B[1] = RF(self._Nij[i,dim-1-j]) 
                mpfr_div(xarg.value,alphas[i],llambda,rnd_re)
                mpfr_add(xarg.value,xarg.value,B._entries[0],rnd_re)
                if j < sym_dim - 1:
                    #z = xarg**-twoz
                    mpc_set_fr(zarg,xarg.value,rnd)
                    mpc_pow(Z[i][j][n][0],zarg,minus_twoz,rnd)
                    #z = mpmath.mp.zeta(twoz,alphas[i]/llambda+B[0])
                    # xarg = alphas[i]/llambda+B[0])
                else:
                    mpz = mpmath.mp.zeta(twozz,xarg)
                    z = CF(mpz.real,mpz.imag)
                    #Z[i][j][l][0]=CF(z.real,z.imag)                
                    mpc_set(Z[i][j][n][0],z.value,rnd)
                mpfr_div(xarg.value,alphas[i],llambda,rnd_re)
                mpfr_neg(xarg.value,xarg.value,rnd_re)
                mpfr_sub(xarg.value,xarg.value,B._entries[1],rnd_re)
                if j < sym_dim - 1:
                    mpc_set_fr(zarg,xarg.value,rnd)
                    mpc_pow(Z[i][j][n][1],zarg,minus_twoz,rnd)
                    #z = xarg**-twoz
                else:
                    #z = mpmath.mp.zeta(zarg,-alphas[i]/llambda-B[1])
                    mpz = mpmath.mp.zeta(twozz,xarg)
                    z = CF(mpz.real,mpz.imag)
                    #Z[i][j][l][1]=CF(z.real,z.imag)
                    mpc_set(Z[i][j][n][1],z.value,rnd)
                #mpc_set(z.value,Z[i][j][n][0],rnd)
                #print i,j,n,0,CC(z.real(),z.imag())
                #mpc_set(z.value,Z[i][j][n][1],rnd)
                #print i,j,n,1,CC(z.real(),z.imag())
        #z = twoz.zeta()
        #mpc_init2(Z[n],prec)        
        mpc_neg(tmp,twozn[n],rnd)        
        if q>3:
            mpc_init2(lpow[n],prec)
            mpc_pow(lpow[n],clambda,tmp,rnd)

    if verbose>0:
        print "Here1"
    #print "N=",abN[0][0]
    cdef int ni,kj,ii
    for n in range(M+1): 
        for k in range(M+1): 
            if (n + k) % 2 == 1:
                sg = -1
            else:
                sg = 1
            for i in range(sym_dim):
                for j in range(sym_dim):
                    ni = i*(M+1)+n
                    kj = j*(M+1)+k
                    fmpz_get_mpz(tmp_mpz,fmpz_mat_entry(Nij,i,j))
                    mpfr_set_z(B._entries[0],tmp_mpz,rnd_re)
                    fmpz_get_mpz(tmp_mpz,fmpz_mat_entry(Nij,i,dim-1-j))
                    mpfr_set_z(B._entries[1],tmp_mpz,rnd_re)
                    #B[1] = RF(Nij[i][dim-1-j])
                    for ii in range(2):
                        mpc_set_ui(AA[ii],0,rnd)
                        if mpfr_sgn(B._entries[ii])<>0:
                            for l in range(k+1):
                                _pochammer_over_fak(&poc,twozn[l],n,
                                                rnd,rnd_re)
                                mpz_bin_uiui(binc,k,l)
                                mpfr_mul_z(poc.re,poc.re,binc,rnd_re)
                                mpfr_mul_z(poc.im,poc.im,binc,rnd_re)
                                mpc_mul(poc,poc,Z[i][j][l+n][ii],rnd)
                                #if ni==0 and kj==32:
                                #mpc_set(z.value,Z[i][j][l+n][ii],rnd)
                                #    print "z[",i,j,l+n,ii,"]=",z
                                if q<>3:
                                    _mpc_mul(&poc,poc,lpow[n+l],t,
                                             rnd,rnd_re)
                                mpc_mul_fr(poc,poc,ajpow[j][k-l],rnd)
                                #if ni==0 and kj==32:
                                #    mpc_set(z.value,poc,rnd)
                                #    print "tmp[",i,j,l+n,ii,"]=",z

                                mpc_add(AA[ii],AA[ii],poc,rnd)
                            if mpfr_sgn(B._entries[ii])>0 and sg==-1:
                                #AA=-AA
                                mpc_neg(AA[ii],AA[ii],rnd)
                    #if ni==0 and kj==32:
                    #    mpc_set(z.value,AA[0],rnd)
                    #    print "AA[0]=",z
                    #    mpc_set(z.value,AA[1],rnd)
                    #    print "AA[1]=",z
                    if (k % 2) == 1:
                        mpc_neg(AA[1],AA[1],rnd)
                    if eps == -1:
                        mpc_neg(AA[1],AA[1],rnd)
                    mpc_add(A[ni][kj],AA[0],AA[1],rnd)
                    mpfr_pow_si(fak1,rhos[j],-k,rnd_re)
                    mpfr_pow_si(fak2,rhos[i],n,rnd_re)
                    mpfr_mul(fak,fak1,fak2,rnd_re)
                    #if ni==0 and kj==32:
                    #    mpc_set(z.value,A[ni][kj],rnd)
                    #    print "A[0,32]=",z
                    #    mpfr_set(tmpx.value,fak,rnd_re)
                    #    print "fak=",tmpx
                    mpc_mul_fr(A[ni][kj],A[ni][kj],fak,rnd)
                    
    if verbose>0:
        print "About to clear"
    mpc_clear(t[0]);mpc_clear(t[1])
    mpc_clear(t[2]);mpc_clear(t[3])
    mpc_clear(twos); mpc_clear(clambda)
    mpc_clear(tmp); mpc_clear(tmp1)
    mpc_clear(poc)
    mpc_clear(zarg); mpc_clear(minus_twoz)
    mpfr_clear(fak1);mpfr_clear(fak2);mpfr_clear(fak)
    if AA<>NULL:
        mpc_clear(AA[0]);  mpc_clear(AA[1])
    if Z<>NULL:
        for i in range(sym_dim):
            if Z[i]<>NULL:
                for j in range(sym_dim):      
                    if Z[i][j]<>NULL:
                        for n from 0 <= n < 2*M+1:
                            mpc_clear(Z[i][j][n][0])
                            mpc_clear(Z[i][j][n][1])
                            sage_free(Z[i][j][n])
                        sage_free(Z[i][j])
                sage_free(Z[i])
        sage_free(Z)
    # if lsum<>NULL:
    #     for n from 0 <= n < 2*M+1:
    #         if lsum[n]<>NULL:
    #             for i from 0<= i < dim:
    #                 if lsum[n][i]<>NULL:
    #                     for j from 0<= j < dim:
    #                         mpc_clear(lsum[n][i][j])
    #                     sage_free(lsum[n][i])
    #             sage_free(lsum[n])
    #     sage_free(lsum)
    
    if ajpow<>NULL:
        for i in range(sym_dim):
            if ajpow[i]<>NULL:
                for n in range(M+1):
                    mpfr_clear(ajpow[i][n])
                sage_free(ajpow[i])
        sage_free(ajpow)    
    if twozn<>NULL:
        for n in range(2*M+1):
            mpc_clear(twozn[n])
        sage_free(twozn)
    if q>3:
        if lpow<>NULL:
            for n in range(2*M+1):
                mpc_clear(lpow[n])
            sage_free(lpow)     
    sig_off()
    #return
                    #print "A[",i*(M+1)+n,",",j*(M+1)+k,"]=",print_mpc(w)
                    #A[i*(M+1)+n,j*(M+1)+k]=w

                        

# def FE_Psi(s,q,prec=0):
#     r"""  Compute the factor Psi in the functional equation Z(1-s)=Psi(s)*Z(s). """
#     if hasattr(s,"prec"):
#         prec = s.prec()
#     elif prec>0:
#         prec = prec
#     else:
#         prec = mpmath.mp.prec
#     mpmath.mp.dps = prec
#     ss=s-0.5; sigma=ss.real(); T=ss.imag()
#     fak1=mpmath.mp.gamma(1.5-s)/mpmath.mp.gamma(s+0.5)
#     mpi = mpmath.mp.mpc(0,1)
#     mp1 = mpmath.mp.mpf(1)
#     mppi = mpmath.mp.pi
#     twopi = mppi*mpmath.mp.mpf(2)
#     twopii = mppi*mpi*mpmath.mp.mpf(2)
#     A=mpmath.mp.mpf(q-2)/mpmath.mp.mpf(q)*mppi
#     f1 = lambda y: mpi*mpmath.mp.tan(mpi*mppi*y)
#     IH1=mpi*mpmath.mp.quad(f1,[0,T])
#     f2 = lambda x: mpmath.mp.mpc(x,T)*mpmath.mp.tan(mppi*mpmath.mp.mpc(x,T))
#     IH2=mpmath.mp.quad(f2,[0,sigma])
#     H1=-A*(IH1+IH2)
#     f3 = lambda y : mpmath.mp.cos(mppi*mpi*y)**-1
#     IE1=mpi*mpmath.mp.quad(f3,[0,T])
#     f4 = lambda x: mpmath.mp.cos(mppi*mpmath.mp.mpc(x,T))**-1
#     IE2=  mpmath.mp.quad(f4,[0,sigma])
#     m=q
#     E1=mppi*(IE1+IE2)/mpmath.mp.mpf(2)
#     for k in range(1,m): #from 1 to m-1 do:
#         km = mpmath.mp.mpf(k)/mpmath.mp.mpf(m)
#         g1 = lambda t: mpmath.mp.exp(twopi*km*t)  / (mp1+mpmath.mp.exp(twopi*t ))+mpmath.mp.exp(-twopi*km*t)/(mp1+mpmath.mp.exp(-twopi*t))
#         IE11 = mpmath.mp.quad(g1,[0,T])
#         #IE11:=I*int( exp(-2*PI*I*k/m*(  I*t))/(1+exp(-2*PI*I*(  I*t )))+
# 	    #     exp(2*PI*I*k/m*(  I*t))/(1+exp(2*PI*I*(  I*t))),t=0..T):
#         g2 = lambda x: mpmath.mp.exp(-twopii*km*mpmath.mp.mpc(x,T))  / (mp1+mpmath.mp.exp(-twopii*mpmath.mp.mpc(x,T) ))+mpmath.mp.exp(twopii*km*mpmath.mp.mpc(x,T))/(mp1+mpmath.mp.exp(twopii*mpmath.mp.mpc(x,T)))

#         IE12 = mpmath.mp.quad(g2,[0,sigma])
#         #IE12:=  int( exp(-2*PI*I*k/m*(x+I*T))/(1+exp(-2*PI*I*(x+I*T )))+
#         #exp(2*PI*I*k/m*(x+I*T))/(1+exp(2*PI*I*(x+I*T))),x=0..sigma): 
#         E1=E1+mppi*(IE11+IE12)/mpmath.mp.mpf(m)/mpmath.mp.sin(mppi*km)
#     P1=(1-2*s)*mpmath.mp.ln(2.0)
#     P=fak1*mpmath.mp.exp(H1+E1+P1)
#     return ComplexField(prec)(P.real,P.imag)


cpdef Gauss_transfer_operator_mpc(RealNumber s,RealNumber t, int r1=0,int r2=0,int k1=0,int k2=0,int verbose=0,ret="dict"):
    r"""
    Compute Gauss transfer operator matrix operator approximation A[r,k] with r1<=r<=r2 and k1<=k<=k2.
    """
    cdef int prec = s.parent().prec()
    cdef MPComplexNumber tmpc
    cdef RealNumber alpha,rho
    cdef mpc_t** A=NULL    
    cdef int r,k
    alpha = RealField(prec)(1)
    rho = RealField(prec)(1.5)
    if verbose>0:
        print "alpha=",alpha
        print "rho=",rho
    if k2==0 and r2==0 and r1>0:
        k2 = r1; r2=r1; k1=0
        r1=0
    A = <mpc_t**> sage_malloc(sizeof(mpc_t*)*(r2-r1))
    for r in range(0,r2-r1):
        A[r] = <mpc_t*> sage_malloc(sizeof(mpc_t)*(k2-k1))
        for k in range(0,k2-k1):
            mpc_init2(A[r][k],prec)
    setup_Gauss_c(A,s.value,t.value,alpha.value,rho.value,r1,r2,k1,k2,verbose)
    cdef dict resd
    if ret=="dict":
        resd = {}
        tmpc = MPComplexField(prec)(1)
        for r in range(0,r2-r1):
            for k in range(0,k2-k1):
                mpc_set(tmpc.value,A[r][k],rnd)
                resd[(r,k)]=1*tmpc
        return resd
    MS = MatrixSpace(MPComplexField(prec),r2-r1,k2-k1)
    resm = Matrix_complex_dense(MS,0)
    for r in range(0,r2-r1):
        for k in range(0,k2-k1):
            mpc_set(resm._matrix[r][k],A[r][k],rnd)
    return resm
    
@cython.cdivision(True)
cdef setup_Gauss_c(mpc_t** A, mpfr_t s, mpfr_t t, mpfr_t alpha, mpfr_t rho,int r1,int r2,int k1,int k2,int verbose=0):
    r"""    Setup the Matrix approximation of the Gauss transfer operator  """
    cdef int i,j,k,l,n,nn
    cdef mpc_t twos,tmp,tmp1,poc, AA
    cdef int prec = mpfr_get_prec(s)
    cdef mpc_t tt[2]
    cdef mpc_t *Z=NULL,*twozn=NULL
    cdef mpfr_t *ajpow
    cdef MPComplexNumber twoz,tmpz,z
    cdef ComplexNumber twozz
    cdef mpz_t binc
    cdef mpfr_t fak1,fak2,fak,rtemp
    cdef RealNumber tmpx
    cdef int lmin,lmax
    lmin = r1+k1; lmax = r2+k2
    sig_on()
    RF = RealField(prec)
    CF = MPComplexField(prec)
    CFF = ComplexField(prec)
    twoz = CF(0); tmpz = CF(0); z = CF(0)
    twozz = CFF(0); tmpx = RF(0)
    mpmath.mp.prec=prec
    if verbose>0:
        print "prec=",prec
        mpfr_set(tmpx.value,alpha,rnd_re)
        print "alpha=",tmpx
        mpfr_set(tmpx.value,rho,rnd_re)
        print "rho=",tmpx
        print "r1,r2=",r1,r2
        print "k1,k2=",k1,k2
    mpfr_init2(fak1,prec);mpfr_init2(fak2,prec);mpfr_init2(fak,prec)
    mpfr_init2(rtemp,prec)
    mpz_init(binc)
    mpc_init2(AA,prec)
    mpc_init2(poc,prec)
    mpc_init2(tmp,prec); mpc_init2(tmp1,prec)
    mpc_init2(twos,prec)
    mpc_init2(tt[0],prec); mpc_init2(tt[1],prec)
    #mpc_init2(t[2],prec); mpc_init2(t[3],prec)
    mpc_set_fr_fr(twos,s,t,rnd)
    mpc_mul_ui(twos,twos,2,rnd)
    tmpz=CF(0)
    #z = CF(0)
    twoz = MPComplexField(prec+20)(0)
    mpc_set(twoz.value,twos,rnd)
    ## Recall that sizeof(mpc_t ***) does not work
    # ajpow[n]=(-alpha)**n
    ajpow=<mpfr_t *> sage_malloc(sizeof(mpfr_t)*(k2))
    if ajpow==NULL: raise MemoryError
    for n in range(k2):
        mpfr_init2(ajpow[n],prec)
        mpfr_neg(ajpow[n],alpha,rnd_re)
        mpfr_pow_si(ajpow[n],ajpow[n],n,rnd_re)

    cdef RealNumber xarg
    cdef mpc_t zarg,minus_twoz
    mpc_init2(zarg,prec)
    mpc_init2(minus_twoz,prec)    
    xarg = RF(0)
    mpfr_add_ui(xarg.value,alpha,1,rnd_re)
    Z  =  <mpc_t *> sage_malloc(sizeof(mpc_t)*(lmax-lmin+2))
    if Z==NULL: raise MemoryError
    twozn  =  <mpc_t *> sage_malloc(sizeof(mpc_t)*(lmax-lmin+2))
    if twozn==NULL: raise MemoryError
    if mpfr_cmp_ui(xarg.value,2)==0:
        if verbose>1:
            print "alpha=1 n0=1"            
        for n in range(0,lmax-lmin+1):
            mpc_init2(Z[n],prec)        
            mpc_init2(twozn[n],prec)
            mpc_add_ui(twozn[n],twos,n+lmin,rnd)
            #_mpc_set(&twoz.value,twozn[n],rnd_re)
            mpc_set(twoz.value,twozn[n],rnd)
            mpc_neg(minus_twoz,twoz.value,rnd)
            #twozz=CFF(twoz.real(),twoz.imag())
            # mpz = mpmath.mp.zeta(twozz,xarg)        
            if verbose==-7:
                mpc_set_ui(z.value,1,rnd)
            else:
                z = twoz.zeta()
                
            #z = CF(mpz.real,mpz.imag)
            _mpc_set(&Z[n],z.value,rnd_re)
            mpc_sub_ui(Z[n],Z[n],1,rnd)
    else:
        for n in range(0,lmax-lmin+1):
            mpc_init2(Z[n],prec)        
            mpc_init2(twozn[n],prec)
            mpc_add_ui(twozn[n],twos,n+lmin,rnd)
            #_mpc_set(&twoz.value,twozn[n],rnd_re)
            mpc_set(twoz.value,twozn[n],rnd)
            mpc_neg(minus_twoz,twoz.value,rnd)
            twozz=CFF(twoz.real(),twoz.imag())
            mpz = mpmath.mp.zeta(twozz,xarg)        
            z = CF(mpz.real,mpz.imag)
            _mpc_set(&Z[n],z.value,rnd_re)
    cdef int tenpercent,chunks
    if verbose==-5:
        return    
    if verbose>0:
        #sys.stdout.write('%d / %d\r' % (i, total))
        sys.stdout.write('Computed Z\n')
        sys.stdout.flush()
        #print "computed Z"
        tenpercent = int ( float(r2-r1)/float(10))
        chunks = 0
    cdef int ni,kj,ii
    #for n in range(r1,r2+1):
    cdef mpc_t **pochammer_vec=NULL
    pochammer_vec = <mpc_t**>sage_malloc(sizeof(mpc_t*)*(k2-k1))
    for k in range(k2-k1):
        pochammer_vec[k] = <mpc_t*>sage_malloc(sizeof(mpc_t)*(r2-r1))
    for k in prange(k2-k1,nogil=True):
        for n in range(r2-r1):
            mpc_init2(pochammer_vec[k][n],prec)            
            _pochammer_over_fak(&pochammer_vec[k][n],twozn[k],n+r1,rnd,rnd_re)
    for n in prange(r2-r1,nogil=True): 
        #for k in range(k1,k2+1):
        for k in range(k2-k1): #
            mpc_set_ui(AA,0,rnd)
            for l in range(k+1):
                #_pochammer_over_fak(&poc,twozn[l],n+r1,rnd,rnd_re)
                #mpc_set_ui(poc,1,rnd)
                _mpc_set(&poc,pochammer_vec[l][n],rnd_re)
                mpz_bin_uiui(binc,k+k1,l)
                mpfr_mul_z(poc.re,poc.re,binc,rnd_re)
                mpfr_mul_z(poc.im,poc.im,binc,rnd_re)
                _mpc_mul(&poc,poc,Z[l+n],tt,rnd,rnd_re)
                _mpc_mul_fr(&poc,poc,ajpow[k-l],rnd,rnd_re)
                _mpc_add(&AA,AA,poc,rnd_re)

            mpfr_pow_si(rtemp,rho,n-k+r1-k1,rnd_re)
            if verbose>1:
                #mpc_set(z.value,AA,rnd)
                printf("A[%d][%d]=%d+i%d",n,k,mpfr_get_d(AA.re,rnd_re),mpfr_get_d(AA.im,rnd_re))
                #    print "A[{0}][{1}]={2}".format(n,k,z)
                #mpfr_set(tmpx.value,rtemp,rnd_re)
                #printf("A[%d][%d]=%d",n,k,z)
            #    print "r^{0}[{1}][{2}]={3}".format(n-k+r1-k1,n,k,tmpx)
                
            _mpc_mul_fr(&AA,AA,rtemp,rnd,rnd_re)
            if (n % 2) == 1:
                mpc_neg(AA,AA,rnd)
            #if verbose>1:
            #   mpc_set(z.value,AA,rnd)
            #    print "A[{0}][{1}]*r^(n-k)={2}".format(n,k,z)
            _mpc_set(&A[n][k],AA,rnd_re)
            #if ni==0 and kj==32:
            #    mpc_set(z.value,A[ni][kj],rnd)
            #    print "A[0,32]=",z
            #    mpfr_set(tmpx.value,fak,rnd_re)
            #    print "fak=",tmpx
        if verbose>0:
            if n % tenpercent == 0:
                #chunks+=1
                printf("Computing row: %d / %d\r",n, r2-r1)
                #sys.stdout.write('Computing row: %d / %d\r' % (n, r2-r1))
                #sys.stdout.flush()
    if verbose>0:
        printf("\n")
        #sys.stdout.write('\n')
        #sys.stdout.flush()
        #print "About to clear"
    mpc_clear(tt[0]);mpc_clear(tt[1])
    mpc_clear(twos); 
    mpc_clear(tmp); mpc_clear(tmp1)
    mpc_clear(poc)
    mpc_clear(zarg); mpc_clear(minus_twoz)
    mpfr_clear(fak1);mpfr_clear(fak2);mpfr_clear(fak)
    mpfr_clear(rtemp)
    mpc_clear(AA)
    if Z<>NULL:
        for n in range(r2-r1+1):
            mpc_clear(Z[n])
        sage_free(Z)
    if ajpow<>NULL:
        for k in range(k2):
            mpfr_clear(ajpow[k])
        sage_free(ajpow)    
    if twozn<>NULL:
        for n in range(0,lmax-lmin):
            mpc_clear(twozn[n])
        sage_free(twozn)
    if pochammer_vec <> NULL:
        for k in range(k2-k1):
            if pochammer_vec[k]<>NULL:
                for n in range(r2-r1):
                    mpc_clear(pochammer_vec[k][n])
                sage_free(pochammer_vec[k])
        sage_free(pochammer_vec)
    sig_off()

       
cpdef my_floor(double x):
    r"""  Return n where n<x<=n+1 and x>0 or  n<=x<n+1 and x<=0.  """
    cdef int nc,nf
    nc = ceil(x)
    nf = floor(x)
    if nc == nf:
        if x<=0:
            return nf
        else:
            return nc
    return nf

    
