# cython: profile=False
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
Algorithms for use together with MySubgroup -- an extension to the standard implementation of subgroups of the modular group.

AUTHOR:

 - Fredrik Stroemberg


 CONVENTION:
 internal functions, beginning with underscore might modify input variables, while other functions does (should) not.


"""


include "sage/ext/interrupt.pxi"  # ctrl-c interrupt block support
include "sage/ext/stdsage.pxi"  # ctrl-c interrupt block support
include "sage/ext/cdefs.pxi"
include "sage/ext/gmp.pxi"
#include "sage/rings/mpc.pxi"

## For multiprecision support
from sage.libs.mpfr cimport *
cdef mpc_rnd_t rnd
cdef mpfr_rnd_t rnd_re
rnd = MPC_RNDNN
rnd_re = GMP_RNDN
from sage.rings.complex_mpc cimport * #MPComplexNumber
from sage.rings.complex_mpc import MPComplexField
from sage.rings.real_mpfr cimport RealNumber,RealField_class
from sage.rings.real_mpfr import RealField
from sage.modular.cusps import Cusp
from sage.rings.infinity import infinity
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.matrix.matrix_integer_2x2 cimport Matrix_integer_2x2
from sage.matrix.matrix_integer_2x2 import Matrix_integer_2x2
from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.integer_ring cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.rational_field import QQ
from sage.all import RR,ZZ,SL2Z,matrix
from copy import deepcopy
from sage.combinat.permutation import Permutation_class
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.functions.all import ceil as pceil
#from sage.rings.rational.Rational import floor as qq_floor
import cython
cdef extern from "math.h":
    double fabs(double)
    double fmax(double,double)
    int ceil(double)
    int floor(double)
    double M_LN10
    double log(double)


cdef extern from "stdio.h":
    int snprintf (char *s, size_t maxlen, char* format, int value)
    #int sprintf (char *s, char* format, int value)
    int sprintf (char *s, char* format, ...)

#cdef extern class sage.matrix.matrix_integer_dense.Matrix_integer_dense as Matrix_integer_dense_class#:
#    pass
#cdef extern class sage.matrix.matrix_integer_2x2.Matrix_integer_2x2 as Matrix_integer_2x2_class:
#    pass



    
# We import these rather than the implementation since it is limited to n<12 
from sage.combinat.permutation_cython cimport *

### A primitive (hopefully efficient) class for permutations
### Note that for small N the SymmetricGroup is more efficient.

from sage.structure.sage_object cimport SageObject
  
import mpmath

### Optimized class for SL2(Z) elements
### represented as integer pointers with 4 elements
cdef class SL2Z_elt(object):
    r"""
    Class for efficient computations with elements of SL(2,Z)
    """
    def __cinit__(self,int a=1,int b=0,int c=0,int d=1):
        assert a*d-b*c==1
        self.ent=NULL
        if self.ent==NULL:
            self.ent=<int*>sage_malloc(sizeof(int)*4)
            if self.ent==NULL:
                raise MemoryError
        self.ent[0]=a
        self.ent[1]=b
        self.ent[2]=c
        self.ent[3]=d
        self._str=''

    cpdef a(self):
        return self.ent[0]
    cpdef b(self):
        return self.ent[1]
    cpdef c(self):
        return self.ent[2]
    cpdef d(self):
        return self.ent[3]

    cpdef is_in_Gamma0N(self,int N):
        if self.ent[2] % N ==0:
            return 1
        return 0
        
    def __init__(self,int a=1,int b=0,int c=0,int d=1):
        pass

    def __dealloc__(self):
        if self.ent<>NULL:
            sage_free(self.ent)
            self.ent=NULL

    def __del__(self):
        self.__dealloc__()
        


    def __repr__(self):
        return "[{0}, {1}, {2}, {3}]".format(self.ent[0],self.ent[1],self.ent[2],self.ent[3])
        

    def __hash__(self):
        if not self._str:
            self._str="{0} {1} {2} {3}".format(self.ent[0],self.ent[1],self.ent[2],self.ent[3])
        return self._str.__hash__()

    
    def __reduce__(self):
        #return (SL2Z_elt,(self.a(),self.b(),self.c(),self.d()))
        return (SL2Z_elt,(self.ent[0],self.ent[1],self.ent[2],self.ent[3]))
    
    
    def __getitem__(self,i):
        if i>=0 and i<4:
            return self.ent[i]
        raise IndexError
    
    def __list__(self):
        return [self.ent[0],self.ent[1],self.ent[2],self.ent[3]]

    def __richcmp__(self,other,op):
        if op==2 or op==3:
            if type(self)<>type(other):
                res=False
            else:
                res = bool(self._eq(other))
        else:
            raise NotImplementedError,"No ordering of SL2Z elements is implemented!"
        if op==3:
            res=not res
        return res

    cpdef _eq(self,SL2Z_elt other):
        if self.ent[0] <> other.ent[0]:
            res=0
        elif self.ent[1]<>other.ent[1]:
            res=0
        elif self.ent[2]<>other.ent[2]:
            res=0
        elif self.ent[3]<>other.ent[3]:
            res=0
        else:
            res=1        
        return res
    
    cpdef acton(self,z):
        den = self.ent[2]*z+self.ent[3]
        if den==0.0:
            return infinity
        else:
            return (self.ent[0]*z+self.ent[1])/den


    def SL2Z(self):
        return SL2Z([self.ent[0],self.ent[1],self.ent[2],self.ent[3]])

    def matrix(self,mtype=0):
        r"""
        Return a matrix representation of self.
        mtype=0 => return matrix of type Matrix_integer_2x2
        mtype=1 => return matrix of type matrix(Integer,2,2,[a,b,c,d])
        """
        if mtype==0:
            MS = MatrixSpace_ZZ_2x2()
            return Matrix_integer_2x2(MS,[self.ent[0],self.ent[1],self.ent[2],self.ent[3]],True,True)
        MS = MatrixSpace(ZZ,2,2)
        return MS([self.ent[0],self.ent[1],self.ent[2],self.ent[3]])
    # cdef _acton_d(self,double complex z):
    #     den = self.ent[2]*z+self.ent[3]
    #     if den==0:
    #         return infinity
    #     else:
    #         return (self.ent[0]*z+self.ent[1])/den
        
    def __pow__(self,k,dummy):
        return self._pow(int(k))

    cpdef _pow(self,int k):
        cdef int i
        cdef SL2Z_elt res
        if k==0:
            return SL2Z_elt(1,0,0,1)
        elif k==1:
            return self
        elif k==-1:
           return self.inverse() 
        else:
            if k>1:
                res=self
                for i from 1<=i<k:
                    res=self._mul(res)
            elif k<-1:
                res= self.inverse()
                return res**(-k)
            return res
        
    cpdef inverse(self):
        return SL2Z_elt(self.ent[3],-self.ent[1],-self.ent[2],self.ent[0])
                    
    def __mul__(self,other):
        if isinstance(other,SL2Z_elt):
            return self._mul(other)
        #elif hasattr(other,'BKZ') and other.nrows()==2 and other.ncols()==0:
        elif isinstance(other,Matrix_integer_dense) and other.nrows()==2 and other.ncols()==2 and other.det()==1:
            return self._mul_mat_id(<Matrix_integer_dense?>other)
        elif isinstance(other,Matrix_integer_2x2) and other.det()==1:
        #elif hasattr(other,'nrows') and other.nrows()==2 and other.ncols()==2:
            return self._mul_mat_i_2x2(<Matrix_integer_2x2?>other)
        raise ValueError,"Can not multiply {0} with {1}".format(type(self),type(other))
        
    
    cpdef _mul(self,SL2Z_elt other,int inv=0):
        r"""
        inv=1 => A^-1*B
        inv=2 => A*B^-1
        """
        cdef SL2Z_elt res
        cdef int* resent
        resent=<int*>sage_malloc(sizeof(int)*4)
        self._mul_c(other.ent,resent,inv)
        res=SL2Z_elt(resent[0],resent[1],resent[2],resent[3])
        if resent<>NULL:
            sage_free(resent)
        return res

    cpdef _mul_mat_id(self,Matrix_integer_dense other,int inv=0):
        r"""
        inv=1 => A^-1*B
        inv=2 => A*B^-1
        """
        cdef SL2Z_elt res
        cdef int* resent=NULL
        resent=<int*>sage_malloc(sizeof(int)*4)
        if resent==NULL: raise MemoryError
        self._mul_c_mpz(other._entries,resent,inv)
        res=SL2Z_elt(resent[0],resent[1],resent[2],resent[3])
        if resent<>NULL:
            sage_free(resent)
        return res



    cpdef _mul_mat_i_2x2(self,Matrix_integer_2x2 other,int inv=0):
        r"""
        inv=1 => A^-1*B
        inv=2 => A*B^-1
        """
        cdef SL2Z_elt res
        cdef int* resent=NULL
        resent=<int*>sage_malloc(sizeof(int)*4)
        if resent==NULL: raise MemoryError
        self._mul_c_mpz(other._entries,resent,inv)
        res=SL2Z_elt(resent[0],resent[1],resent[2],resent[3])
        if resent<>NULL:
            sage_free(resent)
        return res

    cdef _mul_c(self,int* other,int* res,int inv=0):
        if inv==0:
            res[0]=self.ent[0]*other[0]+self.ent[1]*other[2]
            res[1]=self.ent[0]*other[1]+self.ent[1]*other[3]
            res[2]=self.ent[2]*other[0]+self.ent[3]*other[2]
            res[3]=self.ent[2]*other[1]+self.ent[3]*other[3]
        elif inv==2:
            res[0]=self.ent[0]*other[3]-self.ent[1]*other[2]
            res[1]=-self.ent[0]*other[1]+self.ent[1]*other[0]
            res[2]=self.ent[2]*other[3]-self.ent[3]*other[2]
            res[3]=-self.ent[2]*other[1]+self.ent[3]*other[0]        
        elif inv==1:
            res[0]=self.ent[3]*other[0]-self.ent[1]*other[2]
            res[1]=self.ent[3]*other[1]-self.ent[1]*other[3]
            res[2]=-self.ent[2]*other[0]+self.ent[0]*other[2]
            res[3]=-self.ent[2]*other[1]+self.ent[0]*other[3]

    cdef _mul_c_mpz(self,mpz_t* other,int* res,int inv=0):
       cdef int a,b,c,d
       a = mpz_get_si(other[0]); b = mpz_get_si(other[1]);
       c = mpz_get_si(other[2]); d = mpz_get_si(other[3]);
       if inv==0:
           res[0]=self.ent[0]*a+self.ent[1]*c
           res[1]=self.ent[0]*b+self.ent[1]*d
           res[2]=self.ent[2]*a+self.ent[3]*c
           res[3]=self.ent[2]*b+self.ent[3]*d
       elif inv==2:
           res[0]=self.ent[0]*d-self.ent[1]*c
           res[1]=-self.ent[0]*b+self.ent[1]*a
           res[2]=self.ent[2]*d-self.ent[3]*c
           res[3]=-self.ent[2]*b+self.ent[3]*a
       elif inv==1:
           res[0]=self.ent[3]*a-self.ent[1]*c
           res[1]=self.ent[3]*b-self.ent[1]*d
           res[2]=-self.ent[2]*a+self.ent[0]*c
           res[3]=-self.ent[2]*b+self.ent[0]*d
            
## Factoring of matrix in SL2(Z) in S and T using continued fractions.

cpdef factor_matrix_in_sl2z(A,B=None,C=None,D=None,int verbose=0):
    r"""
    Factor a matrix from SL2Z in S and T.
    INPUT:

    - A -- 2x2 integer matrix with determinant 1 (Note: using SL2Z elements are much slower than using integer matrices)

    OUTPUT:
    - l -- list of the form l=[z,n,[a_1,...,a_n]] if A=z T^n ST^{a_1}\cdots ST^{a_n}.

    EXAMPLES::

    sage: A=SL2Z([-28,-5,-67,-12])
    sage: factor_matrix_in_sl2z(A)
    [1, 0, [-2, 2, -2, -6, 0]]
    

    """
    cdef int a,b,c,d
    if D<>None: # If we provide 4 arguments they should be integers. isinstance(A,tuple):
        if A*D-B*C<>1:
            raise ValueError,"Matrix does not have determinant 1!"
        #test = abs(4*C**2+D**2)**
        return fast_sl2z_factor(A,B,C,D)
    else:
        if isinstance(A,SL2Z_elt):
            a=A[0]; b=A[1]; c=A[2]; d=A[3]
        elif isinstance(A,(list,tuple)):
            a,b,c,d=A
        else:
            a=A[0,0]; b=A[0,1]; c=A[1,0]; d=A[1,1]
        if a*d-b*c<>1:
            raise ValueError,"Matrix does not have determinant 1!"
        return fast_sl2z_factor(a,b,c,d)


cpdef factor_matrix_in_sl2z_ncf(A,B=None,C=None,D=None,int check=1,int verbose=0):
    r"""
    Factor a matrix from SL2Z in S and T.
    INPUT:

    - A -- 2x2 integer matrix with determinant 1 (Note: using SL2Z elements are much slower than using integer matrices)

    OUTPUT:
    - l -- list of the form l=[z,n,[a_1,...,a_n]] if A=z T^n ST^{a_1}\cdots ST^{a_n}.

    EXAMPLES::

    sage: A=SL2Z([-28,-5,-67,-12])
    sage: factor_matrix_in_sl2z(A)
    [1, 0, [-2, 2, -2, -6, 0]]
    

    """
    cdef int a,b,c,d
    if D<>None: # If we provide 4 arguments they should be integers. isinstance(A,tuple):
        a,b,c,d=A,B,C,D
    elif isinstance(A,(SL2Z_elt,list,tuple)):
        a,b,c,d=A
    else:
        a=A[0,0]; b=A[0,1]; c=A[1,0]; d=A[1,1]
    if a*d-b*c<>1:
        raise ValueError,"Matrix does not have determinant 1!"
        #test = abs(4*C**2+D**2)**
    if c==0:
        if a==1:
            return [0,b,[]] # T^b
        else:
            return [-1,b,[]] # -T^b
    elif d==0:
        if c==1:
            return [1,a,[0]]
        else:
            return [-1,-a,[0]]
    if verbose>0:
        print "a/c=",Rational(a)/Rational(c)
    nc = nearest_integer_continued_fraction(Rational(a)/Rational(c),verbose=verbose-1)
    if verbose>0:
        print "nc=",nc
    trans = nc[0]
    cdef SL2Z_elt AA,BB
    AA = ncf_to_SL2Z_element(nc)
    BB = AA.inverse()*SL2Z_elt(a,b,c,d) # = +-ST^n= SL2Z([0,-1,1,n])
    #if B.b()<>0:
    if verbose>0:
        print "A=",AA
        print "B=",BB
    # Get sign
    if AA[1]==a:
        sgn = 1
    elif AA[1]==-a:
        sgn = -1
    else:
        raise ArithmeticError,"Could not factor A=[{0},{1},{2},{3}]. Got cf={4}".format(a,b,c,d,nc) 
    if check==1:
        AA = AA*SL2Z_elt(0,-sgn,sgn,sgn*BB.d())
        if not AA==SL2Z_elt(a,b,c,d):
            AA = AA*SL2Z_elt(0,-sgn,sgn,sgn*BB.d())
            raise ArithmeticError,"Could not factor A=[{0},{1},{2},{3}]. Got B={4} and cf={5}".format(a,b,c,d,AA,nc)   
    nc.append(sgn*BB.d())
    nc.remove(trans)
    return [sgn,trans,nc]


    
cdef fast_sl2z_factor(int a,int b,int c,int d,int verbose=0):
    r"""
    Factor a matrix in S and T.
    INPUT:
    - a,b,c,d -- integers with ad-bc = 1, representing A in SL2Z
    OUTPUT:
    - pref,ntrans,l -- tuple with
        - pref == integer, 1 or -1
        - ntrans -- integer
        - l -- list of the form l=[ntrans,[a_1,...,a_n]] if A=z pref*T^ntrans ST^{a_1}\cdots ST^{a_n}.

    EXAMPLES::

    sage: A=SL2Z([-28,-5,-67,-12])
    sage: a,b,c,d=A
    sage: fast_sl2z_factor(a,b,c,d)
    [1, 0, [-2, 2, -2, -6, 0]]

    """
    cdef double x,y
    x = 0.0; y = 2.0
    _apply_sl2z_map_dp(&x,&y,d,-b,-c,a)    
    cdef int aa,bb,cc,dd
    cdef int mapping
    aa=1; bb=0; cc=0; dd=1
    # Now z = A(2i) and we want to pullback z
    l=list()
    cdef int maxc,n,pref,ntrans
    cdef char ch
    ## I think (have not proven) that this number of steps suffices in all cases...
    ## TODO: look up actual bounds
    if d<>0:
        maxc = (abs(c)+abs(d)+1)
    else:
        maxc = (abs(c)+abs(a)+1)
    ntrans=0; pref=1
    for i from 0 <= i <= maxc:
        _apply_one_pb_map(&x,&y,&mapping,&n,&aa,&bb,&cc,&dd)
        # Note that we will never use the same mapping twice in a row here 
        if mapping==1: # then we had a single translation, which might appear at the beginning
            ntrans = n #l.insert(0,('T',n))
        elif mapping==2:
            l.insert(0,n)
        else:
            break
    if i==maxc:
            raise ArithmeticError," Pullback failed! need to increse number of iterations. Used {0}".format(maxc)
    if aa<>a or dd<>d or cc<>c or bb<>b:
        # check -A
        if aa==-a and dd==-d and cc==-c and bb==-b:
            pref=-1
        else:
            raise ArithmeticError," Could not pullback! A={0}, AA={1}".format((a,b,c,d),(aa,bb,cc,dd))
    return [pref,ntrans,l]


cpdef ncf_to_SL2Z_element(l):
    r"""
    EXAMPLE:
    sage: cf=nearest_integer_continued_fraction(pi.n(100),nmax=10);cf
        [3, -7, 16, 294, 3, 4, 5, 15, -3, -2, 2]
        sage: A=ncf_to_matrix_in_SL2Z(cf); A
        [ -411557987 -1068966896]
        [ -131002976  -340262731]
        sage: factor_matrix_in_sl2z_in_S_and_T(A)
        [[3, -7, 16, 294, 3, 4, 5, 15, -2, 2, 3], -1]        
    """
    cdef int j
    cdef SL2Z_elt A
    A=SL2Z_elt(1,l[0],0,1)
    for j in range(1,len(l)):
        A=A*SL2Z_elt(0,-1,1,l[j])
    return A


cpdef ncf_to_SL2Z_matrix(l):
    r"""
    RETURNS matrix
    """

    S=Matrix_integer_2x2([0,-1,1,0])
    T=Matrix_integer_2x2([1,1,0,1])
    A=Matrix_integer_2x2([1,l[0],0,1])
    for j in range(1,len(l)):
        A=A*S*Matrix_integer_2x2([1,l[j],0,1]) #T**l[j]
    return A


cpdef ncf_to_SL2Z_list(list l):
    r"""
    EXAMPLE:

    """
    S=[0,-1,1,0];T=[1,1,0,1]
    A=[1,l[0],0,1]
    cdef int j,ll
    ll = len(l)
    for j from 1<=j<ll: #in range(1,len(l)):
        A = mul_list_maps(A,S)
        A = mul_list_maps_c(A,[1,l[j],0,1])
        #A=A*S*T**l[j]
    return A

cpdef mul_list_maps(list l1,list l2,int inv=0):
    return mul_list_maps_c(l1,l2,inv)

cdef list mul_list_maps_c(list l1,list l2,int inv=0):
    a,b,c,d=l1
    x,y,z,w=l2
    if inv ==1:
        tmp=a
        a=d
        d=tmp
        b=-b; c=-c
    if inv==2:
        tmp=x
        x=w
        w=tmp
        y=-y; z=-z
    aa=a*x+b*z
    bb=a*y+b*w
    cc=c*x+d*z
    dd=c*y+d*w
    return [aa,bb,cc,dd]


cdef _apply_one_pb_map(double *x,double *y,int *mapping,int *n,int* a,int* b,int* c,int* d):
    r"""
    Do one (or rather two) steps in the pullback algorithm
    """
    
    cdef double absval,minus_one,half
    #print "x,y=",x[0],y[0]
    #print "before: a,b,c,d=",a[0],b[0],c[0],d[0]
    mapping[0]=0
    if fabs(x[0]) > 0.5:
        mapping[0] = 1
        half = <double>0.5
        if x[0] >0:
            dx = x[0]-half
            n[0] = -ceil(dx)
        else:
            dx = -x[0]-half
            n[0] = ceil(dx)
        x[0] = x[0]+<double>n[0]
        #print "n=",n[0]
        a[0] = a[0] + n[0]*c[0]
        b[0] = b[0] + n[0]*d[0]
    else:
        n[0] = 0
    absval=x[0]*x[0]+y[0]*y[0]
    if(absval<1.0-1E-15):
        mapping[0] = 2  # meaning that we applied S*T^n
        minus_one = <double>-1.0
        x[0]=minus_one*x[0]/absval
        y[0]=y[0]/absval
        aa=a[0]
        bb=b[0]
        a[0]=-c[0]
        b[0]=-d[0]
        c[0]=aa
        d[0]=bb

    

    #print "after: a,b,c,d=",a[0],b[0],c[0],d[0]              
    #print "ch=",ch

cpdef pullback_to_psl2z_dble_py(double x,double y):
    r""" Pullback to the fundamental domain of SL2Z
         interfacing a fast (typed) Cython version
    INPUT:

        - ``x``  -- double
        - ``x``  -- double

    OUTPUT:

    - ``[u,v,[a,b,c,d]`` --
        u,v -- double
        [a,b,c,d] 4-tuple of integers giving the matrix A
                  such that A(x+iy)=u+iv is inside the fundamental domain of SL2Z

    EXAMPLES::

        sage: pullback_to_psl2z_dble(0.1,0.1)
        [6.9388939039072274e-16, 4.9999999999999991, [5, -1, 1, 0]]


    """
    return pullback_to_psl2z_dble(x,y)


cpdef pullback_to_psl2z_mpfr(RealNumber x,RealNumber y):
    cdef RealNumber xx,yy
    xx = RealNumber(x.parent(),x); yy = RealNumber(y.parent(),y)
    _pullback_to_psl2z_mpfr(xx.value,yy.value)

cdef void _pullback_to_psl2z_mpfr(mpfr_t x,mpfr_t y):
    cdef int a,b,c,d
    cdef double xtmp,ytmp
    xtmp=mpfr_get_d(x,rnd_re)
    ytmp=mpfr_get_d(y,rnd_re)
    pullback_to_psl2z_mat_c(&xtmp,&ytmp,&a,&b,&c,&d)
    _apply_sl2z_map_mpfr(x,y,a,b,c,d)
    
@cython.cdivision(True) 
cdef tuple pullback_to_psl2z_dble(double x,double y):
    r""" Pullback to the fundamental domain of SL2Z
         Fast (typed) Cython version
    INPUT:

        - ``x``  -- double
        - ``x``  -- double

    OUTPUT:

    - ``[u,v,[a,b,c,d]`` --
        u,v -- double
        [a,b,c,d] 4-tuple of integers giving the matrix A
                  such that A(x+iy)=u+iv is inside the fundamental domain of SL2Z

    EXAMPLES::

        sage: pullback_to_psl2z_dble(0.1,0.1)
        [6.9388939039072274e-16, 4.9999999999999991, [5, -1, 1, 0]]


    """
    cdef int a,b,c,d
    [a,b,c,d]=pullback_to_psl2z_mat(x,y)
    cdef double ar=<double>a
    cdef double br=<double>b
    cdef double cr=<double>c
    cdef double dr=<double>d
    cdef double den=(cr*x+dr)**2+(cr*y)**2
    cdef double yout=y/den
    cdef double xout=(ar*cr*(x*x+y*y)+x*(ar*dr+br*cr)+br*dr)/den
    cdef tuple t=(a,b,c,d)
    cdef tuple tt=(xout,yout,t)
    return tt
    #return [xout,yout,[a,b,c,d]]
    

cpdef tuple pullback_to_psl2z_mat(double x,double y):
    cdef int a,b,c,d
    pullback_to_psl2z_mat_c(&x,&y,&a,&b,&c,&d)
    return a,b,c,d

cpdef pullback_to_psl2z_mat2(double x,double y,int a,int b,int c,int d):
    pullback_to_psl2z_mat_c(&x,&y,&a,&b,&c,&d)
    #return a,b,c,d

# unsafe division by zero...
@cython.cdivision(True) 
cdef void pullback_to_psl2z_mat_c(double *x,double *y,int *a,int *b, int *c, int *d,int verbose=0):
    r""" Pullback to the fundamental domain of SL2Z
         Fast (typed) Cython version
    INPUT:

        - ``x``  -- double*
        - ``x``  -- double*

    OUTPUT:

    - ``[a,b,c,d]`` -- 4-tuple of integers giving the matrix A
                      such that A(x+iy) is inside the fundamental domain of SL2Z

    EXAMPLES::

        sage: pullback_to_psl2z_mat(0.1,0.1)
        [5, -1, 1, 0]

        
    """
    cdef int imax=10000
    cdef int done=0
    cdef int aa,bb,i
    a[0]=1; b[0]=0; c[0]=0; d[0]=1
    cdef int res[4]
    cdef double half=<double> 0.5
    cdef double minus_one=<double> -1.0
    cdef int numT
    cdef double absval,dx
    for i from 0<=i<=imax:
        if(fabs(x[0])>half):
            dx=max(x[0]-half,-x[0]-half)
            numT=ceil(dx)
            if(x[0]>0):
                numT=-numT
            x[0]=x[0]+<double>numT
            a[0]=a[0]+numT*c[0]
            b[0]=b[0]+numT*d[0]
        else:
            # We might have to flip
            absval=x[0]*x[0]+y[0]*y[0]
            if(absval<1.0-1E-15):
                x[0]=minus_one*x[0]/absval
                y[0]=y[0]/absval
                aa=a[0]
                bb=b[0]
                a[0]=-c[0]
                b[0]=-d[0]
                c[0]=aa
                d[0]=bb
            else:
                break
    if verbose>0:
        print "(mod)xpb,ypb=",x[0],y[0]
    #cdef tuple t=(a,b,c,d)

cpdef pullback_to_hecke_triangle_mat_dp(double x,double y,double lambdaq):
    #,double *a,double *b,double *c,double *d):
    cdef double a,b,c,d
    pullback_to_hecke_triangle_mat_c(&x,&y,lambdaq,&a,&b,&c,&d)
    return a,b,c,d

cpdef pullback_to_hecke_triangle_mat_mpfr(RealNumber x,RealNumber y,RealNumber lambdaq):
    #,double *a,double *b,double *c,double *d):
    cdef RealNumber a,b,c,d
    cdef int prec = x.prec()
    RF = RealField(prec)
    a = RF(1); b = RF(0); c = RF(0); d = RF(1)
    pullback_to_hecke_triangle_mat_c_mpfr(x.value,y.value,
                                          lambdaq.value,
                                          a.value,b.value,
                                          c.value,d.value)
    return a,b,c,d



cdef void pullback_to_hecke_triangle_mat_c(double *x,double *y,double lambdaq,double *a,double *b,double *c,double *d):
    r""" Pullback to the fundamental domain of Hecke triangle group G_q=<S,T_q>, T_q:z -> z+lambda
         Fast (typed) Cython version
    INPUT:s

        - ``x``  -- double
        - ``x``  -- double

    OUTPUT:

    - ``[a,b,c,d]`` -- 4-tuple of integers giving the matrix A
                      such that A(x+iy) is inside the fundamental domain of SL2Z

    EXAMPLES::

        sage: pullback_to_psl2z_mat(0.1,0.1)
        [5, -1, 1, 0]

    
    """
    cdef int imax=10000
    cdef int i,done=0
    cdef double aa,bb
    a[0]=1; b[0]=0; c[0]=0; d[0]=1
    cdef double half=<double> 0.5*lambdaq
    cdef double minus_one=<double> -1.0
    cdef int numT
    cdef double absval,dx
    for i from 0<=i<=imax:
        if(fabs(x[0])>half):
            dx=max(x[0]-half,-x[0]-half)
            numT=ceil(dx)
            if(x[0]>0):
                numT=-numT
            x[0]=x[0]+<double>numT*lambdaq
            a[0]=a[0]+numT*lambdaq*c[0]
            b[0]=b[0]+numT*lambdaq*d[0]
        else:
            # We might have to flip
            absval=x[0]*x[0]+y[0]*y[0]
            if(absval<1.0-1E-15):
                x[0]=minus_one*x[0]/absval
                y[0]=y[0]/absval
                aa=a[0]
                bb=b[0]
                a[0]=-c[0]
                b[0]=-d[0]
                c[0]=aa
                d[0]=bb
            else:
                break


cdef void pullback_to_hecke_triangle_mat_c_mpfr(mpfr_t x,mpfr_t y,mpfr_t lambdaq,mpfr_t a,mpfr_t b,mpfr_t c,mpfr_t d):
    r""" Pullback to the fundamental domain of Hecke triangle group G_q=<S,T_q>, T_q:z -> z+lambda
         Fast (typed) Cython version
    INPUT:s

        - ``x``  -- double
        - ``x``  -- double

    OUTPUT:

    - ``[a,b,c,d]`` -- 4-tuple of integers giving the matrix A
                      such that A(x+iy) is inside the fundamental domain of SL2Z

    EXAMPLES::

        sage: pullback_to_psl2z_mat(0.1,0.1)
        [5, -1, 1, 0]

    
    """
    cdef int imax=10000
    cdef int i,done=0
    cdef mpfr_t aa,bb,absval,dx,tmpx,tmpy,lambda_half
    prec = mpfr_get_prec(x)
    mpfr_init2(aa,prec);mpfr_init2(bb,prec)
    mpfr_init2(absval,prec);mpfr_init2(dx,prec)
    mpfr_init2(tmpx,prec);mpfr_init2(tmpy,prec)
    mpfr_init2(lambda_half,prec)
    mpfr_set_ui(a,1,rnd_re)
    mpfr_set_ui(b,0,rnd_re)
    mpfr_set_ui(c,0,rnd_re)
    mpfr_set_ui(d,1,rnd_re)
    mpfr_div_ui(lambda_half,lambdaq,2,rnd_re) 
    #cdef double minus_one=<double> -1.0
    cdef int numT
    for i in range(imax+1):
        #if fabs(x)>half:
        mpfr_abs(tmpx,x,rnd_re)
        if mpfr_cmp(tmpx,lambda_half)>0:
            mpfr_sub(tmpx,x,lambda_half,rnd_re)
            mpfr_neg(tmpy,x,rnd_re)
            mpfr_sub(tmpy,tmpy,lambda_half,rnd_re)
            if mpfr_cmp(tmpx,tmpy)>0:
                mpfr_set(dx,tmpx,rnd_re)
            else:
                mpfr_set(dx,tmpy,rnd_re)
            #dx=max(x-half,-x-half)
            numT=ceil(mpfr_get_d(dx,rnd_re))
            if mpfr_sgn(x)>0:
                numT=-numT
            mpfr_mul_si(tmpx,lambdaq,numT,rnd_re)
            mpfr_add(x,x,tmpx,rnd_re)
            #x=x+<double>numT*lambdaq
            mpfr_mul(tmpy,tmpx,c,rnd_re)
            mpfr_add(a,a,tmpy,rnd_re)
            #a=a+numT*lambdaq*c
            mpfr_mul(tmpy,tmpx,d,rnd_re)
            mpfr_add(b,b,tmpy,rnd_re)
            #b=b+numT*lambdaq*d
        else:
            # We might have to flip
            mpfr_mul(tmpx,x,x,rnd_re)
            mpfr_mul(tmpy,y,y,rnd_re)
            mpfr_add(absval,tmpx,tmpy,rnd_re)
            #absval=x*x+y*y
            if mpfr_get_d(absval,rnd_re)<1.0-1E-15:
                mpfr_neg(x,x,rnd_re)
                mpfr_div(x,x,absval,rnd_re)
                #x=minus_one*x/absval
                mpfr_div(y,y,absval,rnd_re)
                #y=y/absval
                mpfr_set(aa,a,rnd_re)
                mpfr_set(bb,b,rnd_re)
                mpfr_neg(a,c,rnd_re)
                mpfr_neg(b,d,rnd_re)
                mpfr_set(c,aa,rnd_re)
                mpfr_set(d,bb,rnd_re)
            else:
                break
    mpfr_clear(aa);mpfr_clear(bb)
    mpfr_clear(dx);mpfr_clear(absval)
    mpfr_clear(tmpx);mpfr_clear(tmpy)
    mpfr_clear(lambda_half)
##
## Faster algorithms for Gamma0(N)
##
@cython.cdivision(True) 
cpdef pullback_to_Gamma0N_mpfr(G,RealNumber x,RealNumber y):
    r"""
    Pullback of x and y to fund. dom. of G=Gamma0(N)
    
    """
    cdef int*** reps
    cdef int nreps,N,j
    cdef int a,b,c,d
    cdef RealNumber xx,yy
    N=G.level()
    nreps=G.index()
    if G._coset_reps_v0==None:
        G._coset_reps_v0 = G._get_coset_reps_from_perms(G.permS,G.permR)
    reps= <int ***> sage_malloc(sizeof(int**) * nreps)
    for j from 0 <=j<nreps:
        reps[j]=<int **> sage_malloc(sizeof(int*) * 2)
        reps[j][0]=<int *> sage_malloc(sizeof(int) * 2)
        reps[j][1]=<int *> sage_malloc(sizeof(int) * 2)
        reps[j][0][0]=G._coset_reps_v0[j][0]
        reps[j][0][1]=G._coset_reps_v0[j][1]
        reps[j][1][0]=G._coset_reps_v0[j][2]
        reps[j][1][1]=G._coset_reps_v0[j][3]

    xx=RealNumber(x.parent(),x)
    yy=RealNumber(x.parent(),y)
    _pullback_to_Gamma0N_mpfr(reps ,nreps, N,xx.value,yy.value)
    a=reps[0][0][0];b=reps[0][0][1];c=reps[0][1][0];d=reps[0][1][1]
    sage_free(reps)
    return xx,yy,a,b,c,d


@cython.cdivision(True) 
cdef void pullback_to_Gamma0N_mpfr_c(G,mpfr_t xout,mpfr_t yout, mpfr_t xin,mpfr_t yin,int *a,int *b,int *c,int *d):
    r"""
    Pullback of x and y to fund. dom. of G=Gamma0(N)
    
    """
    cdef int*** reps
    cdef int nreps,N,j
    #cdef int a,b,c,d
    #cdef RealNumber xx,yy
    N=G.level()
    nreps=G.index()
    if G._coset_reps_v0==None:
        G._coset_reps_v0 = G._get_coset_reps_from_perms(G.permS,G.permR)
    reps= <int ***> sage_malloc(sizeof(int**) * nreps)
    for j in range(nreps):
        reps[j]=<int **> sage_malloc(sizeof(int*) * 2)
        reps[j][0]=<int *> sage_malloc(sizeof(int) * 2)
        reps[j][1]=<int *> sage_malloc(sizeof(int) * 2)
        reps[j][0][0]=G._coset_reps_v0[j][0]
        reps[j][0][1]=G._coset_reps_v0[j][1]
        reps[j][1][0]=G._coset_reps_v0[j][2]
        reps[j][1][1]=G._coset_reps_v0[j][3]

    mpfr_set(xout,xin,rnd_re)
    mpfr_set(yout,yin,rnd_re)
    #xx=RealNumber(x.parent(),x)
    #yy=RealNumber(x.parent(),y)
    _pullback_to_Gamma0N_mpfr(reps ,nreps, N,xout,yout)
    a[0]=reps[0][0][0]
    b[0]=reps[0][0][1]
    c[0]=reps[0][1][0]
    d[0]=reps[0][1][1]
    sage_free(reps)
    #print "here: a,b,c,d=",a[0],b[0],c[0],d[0]
    #return xx,yy,a,b,c,d


cdef void _pullback_to_Gamma0N_mpfr(int*** reps ,int nreps, int N,mpfr_t x,mpfr_t y):
    r""" Mpfr version of pullback alg for Gamma_0(N).
    (the reason for this restriction is the more complicated membership tests otherwise involved)

    INPUT:
    
     - ``reps`` -- list of coset-representatives of Gamma_0(N)
     - ``nreps`` -- number of coset-representatives of Gamma_0(N)
     - ``N`` -- level
     - ``x`` -- real
     - ``y`` -- real > 0

     ON RETURN:
     
     - x,y are overwritten with pull-backed points
     - reps[0] is overwritten with the pullback map
     

    EXAMPLES::

    

    """
    cdef int a,b,c,d
    cdef double x1,y1
    x1=mpfr_get_d(x,rnd_re)
    y1=mpfr_get_d(y,rnd_re)
    pullback_to_psl2z_mat_c(&x1,&y1,&a,&b,&c,&d)
    #A=SL2Z([a,b,c,d])
    cdef int a1,b1,c1,d1
    cdef int** V
    V= <int **> sage_malloc(sizeof(int*) * 2)
    V[0]= <int *> sage_malloc(sizeof(int) * 2)
    V[1]= <int *> sage_malloc(sizeof(int) * 2)
    # V is now a 2x2 matrix over int 
    cdef int j
    for j from 0<=j<=nreps:
        V=reps[j]
        c1 = V[1][0]*a+V[1][1]*c
        if c1 % N == 0:
            a1 = V[0][0]*a+V[0][1]*c
            b1 = V[0][0]*b+V[0][1]*d
            d1 = V[1][0]*b+V[1][1]*d
            a=a1; b=b1; c=c1; d=d1
            #_apply_sl2z_map_mpc(xout,yout,xin,yin,a,b,c,d)
            _apply_sl2z_map_mpfr(x,y,a,b,c,d)
            reps[0][0][0]=a; reps[0][0][1]=b; reps[0][1][0]=c; reps[0][1][1]=d; 
            return 
    if V[0]<>NULL:
        sage_free(V[0])
    if V[1]<>NULL:
        sage_free(V[1])    
    if V<>NULL:
        sage_free(V)
    raise Exception,"Did not find pullback! A=[%s %s, %s %s]" %(a,b,c,d)


cpdef tuple pullback_to_Gamma0N_dp(G,double x,double y,int verbose=0):
    r"""
    Pullback of x and y to fund. dom. of G=Gamma0(N)
    
    """
    cdef int*** reps=NULL
    cdef int nreps,N,j
    cdef int a,b,c,d
    cdef double xx,yy
    N=G.level()
    nreps=G.index()
    if G._coset_reps_v0==None:
        G._coset_reps_v0 = G._get_coset_reps_from_perms(G.permS,G.permR)
    reps= <int ***> sage_malloc(sizeof(int**) * nreps)
    for j from 0 <=j<nreps:
        reps[j]=NULL
        reps[j]=<int **> sage_malloc(sizeof(int*) * 2)
        if reps[j]==NULL:
            raise MemoryError
        reps[j][0]=NULL;reps[j][1]=NULL
        reps[j][0]=<int *> sage_malloc(sizeof(int) * 2)
        if reps[j][0]==NULL:
            raise MemoryError
        reps[j][1]=<int *> sage_malloc(sizeof(int) * 2)
        if reps[j][1]==NULL:
            raise MemoryError
        reps[j][0][0]=G._coset_reps_v0[j][0]
        reps[j][0][1]=G._coset_reps_v0[j][1]
        reps[j][1][0]=G._coset_reps_v0[j][2]
        reps[j][1][1]=G._coset_reps_v0[j][3]
        #print "reps[",j,",00=",reps[j][0][0]
        #print "reps[",j,",01=",reps[j][0][1]
        #print "reps[",j,",10=",reps[j][1][0]
        #print "reps[",j,",11=",reps[j][1][1]        
    xx=x
    yy=y
    _pullback_to_Gamma0N_dp(reps ,nreps, N,&xx,&yy,&a,&b,&c,&d,verbose)
    #a=reps[0][0][0];b=reps[0][0][1];c=reps[0][1][0];d=reps[0][1][1]
    if reps<>NULL:
        for j from 0 <=j<nreps:
            if reps[j]<>NULL:
                if reps[j][0]<>NULL:
                    sage_free(reps[j][0])
                if reps[j][1]<>NULL:
                    sage_free(reps[j][1])
                sage_free(reps[j])
        sage_free(reps)
    return xx,yy,a,b,c,d


@cython.cdivision(True) ## N should never be 0
cdef void _pullback_to_Gamma0N_dp(int*** reps ,int nreps, int N,double *x,double *y,int *a,int *b,int *c,int *d,int verbose):
    r""" Mpfr version of pullback alg for Gamma_0(N).
    (the reason for this restriction is the more complicated membership tests otherwise involved)

    INPUT:
    
     - ``reps`` -- list of coset-representatives of Gamma_0(N)
     - ``nreps`` -- number of coset-representatives of Gamma_0(N)
     - ``N`` -- level
     - ``x`` -- real
     - ``y`` -- real > 0

     ON RETURN:
     
     - x,y are overwritten with pull-backed points
     - reps[0] is overwritten with the pullback map
     

    EXAMPLES::

    """
    #cdef int a,b,c,d
    if verbose>0:
        print "x_in,y_in=",x[0],y[0]
    pullback_to_psl2z_mat_c(x,y,a,b,c,d,verbose)
    if verbose>0:
        print "x_pb,y_pb(mod)=",x[0],y[0]
        print "map:",a[0],b[0],c[0],d[0]
    cdef int a1,b1,c1,d1
    cdef int** V=NULL
    V= <int **> sage_malloc(sizeof(int*) * 2)
    if V==NULL:
        raise MemoryError
    V[0]=NULL; V[1]=NULL
    V[0]= <int *> sage_malloc(sizeof(int) * 2)
    if V[0]==NULL:
        raise MemoryError
    V[1]= <int *> sage_malloc(sizeof(int) * 2)
    if V[1]==NULL:
        raise MemoryError
    # V is now a 2x2 matrix over int 
    cdef int j,done=0
    for j from 0<=j<nreps:
        #print "reps[",j,",00=",reps[j][0][0]
        #print "reps[",j,",01=",reps[j][0][1]
        #print "reps[",j,",10=",reps[j][1][0]
        #print "reps[",j,",11=",reps[j][1][1]        
        V[0][0]=reps[j][0][0]
        V[0][1]=reps[j][0][1]
        V[1][0]=reps[j][1][0]
        V[1][1]=reps[j][1][1]
        c1 = V[1][0]*a[0]+V[1][1]*c[0]
        # Check if A in Gamma_0(N)*V_j^-1
        # <=> A*V_j  in Gamma_0(N)
        # and we then apply VjA to the point
        if c1 % N == 0:
            a1 = V[0][0]*a[0]+V[0][1]*c[0]
            b1 = V[0][0]*b[0]+V[0][1]*d[0]
            c1 = V[1][0]*a[0]+V[1][1]*c[0]
            d1 = V[1][0]*b[0]+V[1][1]*d[0]
            a[0]=a1; b[0]=b1; c[0]=c1; d[0]=d1
            #_apply_sl2z_map_mpc(xout,yout,xin,yin,a,b,c,d)
            if verbose>0:
                print "Coset rep nr. ",j,"=",a[0],b[0],c[0],d[0]
                #_apply_sl2z_map_dp(x,y,a,b,c,d)
            _apply_sl2z_map_dp(x,y,V[0][0],V[0][1],V[1][0],V[1][1])
            #reps[0][0][0]=a; reps[0][0][1]=b; reps[0][1][0]=c; reps[0][1][1]=d; 
            done=1
            break #return 
    if verbose>0:
        print "pb_to_Gamma:",x[0],y[0]
    if V[0]<>NULL:
        sage_free(V[0])
    if V[1]<>NULL:
        sage_free(V[1])
    if V<>NULL:
        sage_free(V)
    if done==0:
        raise Exception,"Did not find pullback! A=[%s %s, %s %s]" %(a[0],b[0],c[0],d[0])

cpdef apply_sl2z_map(x,y,A):
    cdef int a,b,c,d
    cdef RealNumber xx,yy
    [a,b,c,d]=A
    #if isinstance(x,RealNumber):
    xx = RealNumber(x.parent(),x); yy = RealNumber(y.parent(),y)
    _apply_sl2z_map_mpfr(xx.value,yy.value,a,b,c,d)
    return xx,yy

cpdef apply_sl2z_map_dp(double x,double y,int a, int b,int c,int d):
    #cdef int a,b,c,d
    cdef double xx,yy
    #[a,b,c,d]=A
    xx = <double>x; yy=<double>y
    _apply_sl2z_map_dp(&xx,&yy,a,b,c,d)
    return xx,yy

@cython.cdivision(True) ## den should never be 0
cdef void _apply_sl2z_map_dp(double *x,double *y,int a, int b, int c,int d):
    cdef double ar,br,cr,dr
    cdef double den,tmp1,tmp2
    ar=<double>a
    br=<double>b
    cr=<double>c
    dr=<double>d
    den=(cr*x[0]+dr)**2+(cr*y[0])**2
    tmp1 = ar*cr*(x[0]*x[0]+y[0]*y[0])+(a*d+b*c)*x[0]+br*dr
    x[0] = tmp1/den
    y[0] = y[0]/den


cpdef apply_sl2r_map_dp(double xin,double yin, double a,double b,double c,double d):
    cdef double xout,yout
    xout=xin; yout=yin
    _apply_sl2r_map_dp(&xout,&yout,a,b,c,d)
    return xout,yout
    #print "xout,yout=",xout,yout

@cython.cdivision(True) ## den should never be 0
cdef void _apply_sl2r_map_dp(double *x,double *y,double a,double b,double c,double d):
    cdef double den,tmp1,tmp2
    den=(c*x[0]+d)**2+(c*y[0])**2
    tmp1 = a*c*(x[0]*x[0]+y[0]*y[0])+(a*d+b*c)*x[0]+b*d
    x[0] = tmp1/den
    y[0] = y[0]/den    
    #print "x,y=",x[0],y[0]
    #return x,y


cpdef apply_sl2z_map_mpfr(RealNumber x,RealNumber y,int a,int b,int c,int d):
    cdef RealNumber xx,yy
    xx = RealNumber(x.parent(),x); yy = RealNumber(y.parent(),y)
    _apply_sl2z_map_mpfr(xx.value,yy.value,a,b,c,d)
    return xx,yy

cdef void _apply_sl2z_map_mpfr(mpfr_t x,mpfr_t y,int a,int b,int c,int d):
    cdef RealNumber xx,yy
    RF = RealField(mpfr_get_prec(x))
    xx=RF(0);yy=RF(0)
    mpfr_set(xx.value,x,rnd_re);    mpfr_set(yy.value,y,rnd_re)
    #print "apply_map:prec=",mpfr_get_prec(x)
    #print "In: x,y,a,b,c,d=",xx,yy,a,b,c,d
    cdef mpfr_t den,tmp1,tmp2
    cdef prec=mpfr_get_prec(x)
    #cdef RealNumber xtmp
    #xtmp = RealField(prec)(1)
    mpfr_init2(den,prec)
    mpfr_init2(tmp1,prec)
    mpfr_init2(tmp2,prec)
    mpfr_mul_si(den,x,c,rnd_re)
    mpfr_add_si(den,den,d,rnd_re)
    mpfr_mul(den,den,den,rnd_re)
    mpfr_mul_si(tmp1,y,c,rnd_re)
    mpfr_mul(tmp1,tmp1,tmp1,rnd_re)
    mpfr_add(den,den,tmp1,rnd_re)
    #den=(cr*x+dr)**2+(cr*y)**2
    mpfr_mul(tmp1,y,y,rnd_re)
    mpfr_mul(tmp2,x,x,rnd_re)
    mpfr_add(tmp1,tmp1,tmp2,rnd_re)
    cdef int itmp
    itmp=a*c
    mpfr_mul_si(tmp1,tmp1,itmp,rnd_re)
    itmp=a*d+b*c
    mpfr_mul_si(tmp2,x,itmp,rnd_re)
    itmp=b*d
    mpfr_add_si(tmp2,tmp2,itmp,rnd_re)
    mpfr_add(tmp2,tmp2,tmp1,rnd_re)
    mpfr_div(x,tmp2,den,rnd_re)
    mpfr_div(y,y,den,rnd_re)
    mpfr_set(xx.value,x,rnd_re);    mpfr_set(yy.value,y,rnd_re)
    #print "Out: x,y=",xx,yy
    #yy=y/den
    mpfr_clear(den)
    mpfr_clear(tmp1)
    mpfr_clear(tmp2)


cpdef apply_gl2z_map_mpfr(RealNumber x,RealNumber y,int a,int b,int c,int d):
    cdef RealNumber xx,yy
    xx = RealNumber(x.parent(),x)
    yy = RealNumber(y.parent(),y)
    _apply_gl2z_map_mpfr(xx.value,yy.value,a,b,c,d)
    return xx,yy
    
cdef void _apply_gl2z_map_mpfr(mpfr_t x,mpfr_t y,int a,int b,int c,int d):
    #    cdef mpfr_t dtmp1,tmp2
    #cdef prec=mpfr_get_prec(x)
    _apply_sl2z_map_mpfr(x,y,a,b,c,d)
    #mpfr_set(det,a*d-b*c,rnd_re)
    mpfr_mul_si(y,y,a*d-b*c,rnd_re)


cpdef normalize_point_to_cusp_mpfr(G ,int ca,int cb,RealNumber x,RealNumber y,int inv=0):    
    r"""
    Compute the normalized point with respect to the cusp cu
    
    """
    cdef int wi=G._cusp_data[(ca,cb)]['width']
    if cb==0 and wi==1: # Inf is the first cusp
        return [x,y]
    #
    cdef RealNumber xx,yy
    cdef int prec=x.parent().prec()
    xx = RealNumber(x.parent(),x); yy = RealNumber(y.parent(),y)
    if inv==1:
        [d,b,c,a]=G.cusp_normalizer((ca,cb))
        b=-b; c=-c
    else:
        [a,b,c,d]=G.cusp_normalizer((ca,cb))
    _normalize_point_to_cusp_mpfr(xx.value,yy.value,a,b,c,d,wi,inv)
    return xx,yy

    #return [u,v]
cdef void normalize_point_to_cusp_mpfr_c(mpfr_t xout, mpfr_t yout, int * N ,int wi, mpfr_t xin, mpfr_t yin, int inv):    
    r"""
    Compute the normalized point with respect to the cusp cu
    
    """
    #xx = RealNumber(x.parent(),x); yy = RealNumber(y.parent(),y)
    cdef int a,b,c,d
    if inv==1:
        d = N[0]; b = N[1]; c = N[2]; a = N[3]
        #[d,b,c,a]=G.cusp_normalizer((ca,cb))
        b=-b; c=-c
    else:
        a = N[0]; b = N[1]; c = N[2]; d = N[3]
        #[a,b,c,d]=G.cusp_normalizer((ca,cb))
    mpfr_set(xout,xin,rnd_re)
    mpfr_set(yout,yin,rnd_re)
    #print "wi=",wi
    #print "xin,yin=",mpfr_get_d(xin,rnd_re),mpfr_get_d(yin,rnd_re)
    #print "xout,yout0=",mpfr_get_d(xout,rnd_re),mpfr_get_d(yout,rnd_re)
    _normalize_point_to_cusp_mpfr(xout,yout,a,b,c,d,wi,inv)
    #print "xout1,yout1=",mpfr_get_d(xout,rnd_re),mpfr_get_d(yout,rnd_re)

cdef void _normalize_point_to_cusp_mpfr(mpfr_t x,mpfr_t y,int a,int b,int c,int d, int wi, int inv=0):    
    r"""
    Compute the normalized point with respect to the cusp cu
    """
    if inv<>1 and wi<>1:
        mpfr_mul_ui(x,x,wi,rnd_re)
        mpfr_mul_ui(y,y,wi,rnd_re)
        #x=x*wi; y=y*wi
    #assert a*d-b*c==1
    #print "z=",mpfr_get_d(x,rnd_re),mpfr_get_d(y,rnd_re)
    #print "N=",a,b,c,d
    _apply_sl2z_map_mpfr(x,y,a,b,c,d)
    #print "Nz=",mpfr_get_d(x,rnd_re),mpfr_get_d(y,rnd_re)
    if inv and wi<>1:
        mpfr_div_ui(x,x,wi,rnd_re)
        mpfr_div_ui(y,y,wi,rnd_re)
        #u=u/wi
        #v=v/wi

cpdef normalize_point_to_cusp_dp(G,cu,double x,double y,int inv=0):    
    r"""
    Compute the normalized point with respect to the cusp cu
    """
    if(cu==Cusp(infinity) and G.cusp_width(cu)==1): # Inf is the first cusp
        return [x,y]
    cdef int wi=G.cusp_width(cu)
    cdef double xx,yy
    xx = x; yy=y
    #xx = RealNumber(x.parent(),x); yy = RealNumber(y.parent(),y)
    if inv==1:
        [d,b,c,a]=G.cusp_normalizer(cu)
        b=-b; c=-c
    else:
        [a,b,c,d]=G.cusp_normalizer(cu)
    _normalize_point_to_cusp_dp(&xx,&yy,a,b,c,d,wi,inv)
    return xx,yy
    #return [u,v]


cdef void _normalize_point_to_cusp_dp(double *x,double *y,int a,int b,int c,int d, int wi, int inv=0):    
    r"""
    Compute the normalized point with respect to the cusp cu
    If inv=0 we use sigma_j
    If inv=0 we use sigma_j^-1
    where sigma_j = rho_j A_j
    """
    cdef double wir = <double>wi
    #print "x,y0=",x[0],y[0]
    if inv<>1 and wi<>1:
        x[0]=x[0]*wir
        y[0]=y[0]*wir
    #assert a*d-b*c==1
    #print "x,y1=",x[0],y[0]
    _apply_sl2z_map_dp(x,y,a,b,c,d)
    if inv and wi<>1:
        x[0]=x[0]/wir
        y[0]=y[0]/wir



cdef void _normalize_point_to_cusp_real_dp(double *x,double *y,int a,int b,int c,int d, double wi, int inv=0):    
    r"""
    Compute the normalized point with respect to the cusp cu
    where the width is a real number
    """
    cdef double wir = <double>wi
    if inv<>1 and wi<>1:
        x[0]=x[0]*wir
        y[0]=y[0]*wir
    #assert a*d-b*c==1
    _apply_sl2z_map_dp(x,y,a,b,c,d)
    if inv and wi<>1:
        x[0]=x[0]/wir
        y[0]=y[0]/wir  



cpdef normalize_point_to_cusp_mpmath(G,x,y,cu, inv=0,mp_ctx=mpmath.mp):    
    r"""
    Compute the normalized point with respect to the cusp cu
    """
    #print "cu=",cu
    if(cu==Cusp(infinity) and G.cusp_width(cu)==1): # Inf is the first cusp
        return [x,y]
    cdef int wi=G.cusp_width(cu)
    #print "wi=",wi
    #cdef double xx,yy
    #xx = x; yy=y
    #cdef int prec=x.parent().prec()
    #xx = RealNumber(x.parent(),x); yy = RealNumber(y.parent(),y)
    #print "wi=",wi
    wir = mp_ctx.mpf(wi)
    #print "wir=",wir
    if inv==1:
        [d,b,c,a]=G.cusp_normalizer(cu)
        b=-b; c=-c
    else:
        [a,b,c,d]=G.cusp_normalizer(cu)
    if inv<>1 and wi<>1:
        x=x*wir
        y=y*wir
    z = mp_ctx.mpc(x,y)
    tmp = (a*z+b)/(c*z+d)
    x=tmp.real
    y=tmp.imag
    if inv and wi<>1:
        x=x/wir
        y=y/wir
    return x,y


cpdef closest_vertex(in_vertex_maps, wids,int nv_in,double x,double y,int verbose=0): 
    #cdef double** vertex_maps=NULL
    cdef int** vertex_maps =NULL
    cdef double* widths=NULL
    cdef double xx,yy,a,b,c,d
    cdef int vmax
    cdef int nv=int(nv_in)
    if y<=0:
        raise ArithmeticError,"Can not have y<=0! Got y={0}".format(y)
    vertex_maps=<int**>sage_malloc(sizeof(int*)*nv)
    if vertex_maps==NULL: raise MemoryError
    for i from 0<=i<nv:
        vertex_maps[i]=<int*>sage_malloc(sizeof(int)*4)
        #if vertex_maps[i]==NULL: raise MemoryError
        #a,b,c,d=in_vertex_maps[i]
        vertex_maps[i][0]=<SL2Z_elt>in_vertex_maps[i].a()
        vertex_maps[i][1]=<SL2Z_elt>in_vertex_maps[i].b()
        vertex_maps[i][2]=<SL2Z_elt>in_vertex_maps[i].c()
        vertex_maps[i][3]=<SL2Z_elt>in_vertex_maps[i].d()
    widths=<double*>sage_malloc(sizeof(double)*nv)
    for i from 0<=i<nv:
        widths[i]=<double>wids[i]
    xx=<double>x
    yy=<double>y
    vmax = closest_vertex_dp_c(nv,vertex_maps,widths,&xx,&yy,int(verbose))     
    #print "vmax=",vmax
    if vmax==-1:
        raise ArithmeticError," Did not find closest vertex to z={0}+i{1}".format(xx,yy)
    if vertex_maps<>NULL:
        for i from 0<=i<nv:
            if vertex_maps[i]<>NULL:
                sage_free(vertex_maps[i])
        sage_free(vertex_maps)
    if widths<>NULL:
        sage_free(widths)
    return vmax

cdef int closest_vertex_dp_c(int nv,int **vertex_maps,double *widths,double* x,double* y,int verbose=0):
    cdef double y2,a,b,c,d,den,w,ymax
    cdef int i
    vmax=-1; ymax=-1
    for i from 0<=i <nv:
        w=<double>widths[i]
        a=<double>(vertex_maps[i][0])
        b=<double>(vertex_maps[i][1])
        c=<double>(vertex_maps[i][2])
        d=<double>(vertex_maps[i][3])
        den=(c*x[0]+d)**2 +(c*y[0])**2
        y2=y[0]/den/w # We only want the y coordinate
        #print "a,b,c,d,w=",a,b,c,d,w
        if verbose>0:
            print "i=",i
            print "vertex_map:",a,b,c,d
            print "width=",w
            print "y2=",y2
        if y2>ymax:
            ymax=y2
            vmax=i
    #print "vmax=",vmax
    return vmax

## def are_transitive_permutations(E,R):
##     r""" Check that E,R are transitive permutations, i.e. that <E,R>=S_N

##     INPUT:

##          - ``E`` -- permutation on N letters 
##          - ``R`` -- permutation on N letters

##              - E and R can be in any of the following formats:

##                  - list [a1,a2,...,aN]
##                  - member of Permutations(N)
##                  - member of SymmetricGroup(N)

##      OUTPUT:

##      - bool  


##      EXAMPLES::

##          sage: E=Permutations(4)([1,2,4,3]); E.to_cycles()
##          [(1,), (2,), (3, 4)]
##          sage: R=Permutations(4)([2,1,3,4]); R.to_cycles()
##          [(1, 2), (3,), (4,)]
##          sage: are_transitive_permutations(E,R)
##          False
##          sage: R=Permutations(4)([2,3,1,4]); R.to_cycles()
##          [(1, 2, 3), (4,)]
##          sage: are_transitive_permutations(E,R)
##          True
##          sage: ES=SymmetricGroup(4)([1,2,4,3]);ES
##          (3,4)
##          sage: ER=SymmetricGroup(4)([2,3,1,4]);ER
##          (1,2,3)
##          sage: are_transitive_permutations(ES,RS)
##          True

##     """
##     gotten=list()    
##     t0=isinstance(E,list)
##     if(t0):
##         El=E; Rl=R
##     else:
##         t1=isinstance(E,Permutation_class) # constructed from Permutations(N)
##         if(t1):
##             El=list(E)
##             Rl=list(R)
##         else:
##             t2=isinstance(E,PermutationGroupElement) # constructed from SymmetricGroup(N)
##             if(t2):
##                 El=E.list()
##                 Rl=R.list()           
##             else:
##                 raise TypeError, "Indata need to be Permuatations! Got E=%s of type=%s" %(E,type(E))
##     N=len(El)
##     gotten.append(Rl[0])
##     gotten0=deepcopy(gotten)
##     for j in range(N):
##         for x in gotten0:
##             y=El[x-1]
##             if(gotten.count(y)==0):
##                 gotten.append(y)
##             yy=Rl[y-1]
##             if(gotten.count(yy)==0):
##                 gotten.append(yy)
##             yyy=Rl[yy-1]
##             if(gotten.count(yyy)==0):
##                 gotten.append(yyy)
##         if(gotten0==gotten):
##             if(len(gotten)<>N):
##                 return False
##             else:
##                 return True
##         gotten0=deepcopy(gotten)
##     return False


cpdef nearest_integer_continued_fraction(x,int nmax=0,int verbose=0):
    r""" Nearest integer continued fraction of x
    where x is a rational number, and at n digits

    EXAMPLES::

        sage: nearest_integer_continued_fraction(0.5)
        [1, 2]
        nearest_integer_continued_fraction(pi.n(100),nmax=10)
        [3, -7, 16, 294, 3, 4, 5, 15, -3, -2, 2]

    """
    if isinstance(x,float):
        return nearest_integer_continued_fraction_dble(<double>x,nmax)
    elif isinstance(x,RealNumber):
        return nearest_integer_continued_fraction_real(<RealNumber>x,nmax)
    if nmax <= 0:
        nmax=100000    
    cdef int jj=0 
    cdef list cf=list()
    cdef int n
    cdef Rational x0,x1,zero
    n=nearest_integer(x)
    cf.append(n)
    x1=Rational(x-n)
    zero=Rational(0)
    if verbose>0:
        print "x=",x1
        print "nmax=",nmax
    for jj in range(nmax): #while jj<nmax and x1<>0:
        x0 = -x1**-1
        n=nearest_integer(x0)
        if verbose>0:
            print "x0=",x0
            print "n=",n
        x1=x0 - n
        if verbose>0:
            print "x1=",x1
        cf.append(n)
        if x1==0:
            break        
        #jj=jj+1
        
    return cf


@cython.cdivision(True) 
cpdef nearest_integer_continued_fraction_dble(double x,int nmax=0):
    r""" Nearest integer continued fraction of x
    where x is a rational number, and at n digits

    EXAMPLES::

        sage: nearest_integer_continued_fraction(0.5)
        [1, 2]
        nearest_integer_continued_fraction(pi.n(100),nmax=10)
        [3, -7, 16, 294, 3, 4, 5, 15, -3, -2, 2]

    """
    if nmax <= 0:
        nmax=100  # For non-rational numbrs  we don't want so many convergents
    cdef int jj=0 
    cdef list cf=[]
    cdef int n
    cdef double minus_one = <double>-1
    cdef double x0,x1
    cdef double eps = 2.0**-45
    n=nearest_integer_dble(x)
    cf.append(n)
    x1=<double>x-n
    while jj<nmax and abs(x1)>eps:
        x0 = minus_one/x1
        n=nearest_integer_dble(x0)
        x1=x0-<double>n
        cf.append(n)
        jj=jj+1 
    return cf

cpdef nearest_integer_continued_fraction_real(RealNumber x,int nmax=0):
    r""" Nearest integer continued fraction of x
    where x is a rational number, and at n digits

    EXAMPLES::

        sage: nearest_integer_continued_fraction(0.5)
        [1, 2]
        nearest_integer_continued_fraction(pi.n(100),nmax=10)
        [3, -7, 16, 294, 3, 4, 5, 15, -3, -2, 2]

    """
    if nmax <= 0:
        nmax=100  # For non-rational numbrs  we don't want so many convergents
    cdef int jj=0 
    cdef list cf=[]
    cdef Integer n
    cdef RealNumber x0,x1,eps
    RF = x.parent()
    eps = RF(2.0)**(8-RF.prec())
    n=nearest_integer_real(x)
    cf.append(n)
    x1=RF(x-n)
    x0 = RF(0)
    while jj<nmax and abs(x1)>eps:
        x0 = -x1**-1
        n=nearest_integer_real(x0)
        x1=x0-n
        cf.append(n)
        jj=jj+1 
    return cf

cpdef nearest_integer_real(RealNumber x):
    return (x+x.parent()(0.5)).floor()

cpdef nearest_integer(Rational x):
    r""" Returns the nearest integer to x: [x]
    using the convention that 
    [1/2]=0 and [-1/2]=0


    EXAMPLES::

        sage: nearest_integer(0)
        0
        sage: nearest_integer(0.5)
        1
        sage: nearest_integer(-0.5)
        0
        sage: nearest_integer(-0.500001)
        -1
    """
    return (x+Rational(1)/Rational(2)).floor()

cpdef nearest_integer_dble(double x):
    r""" Returns the nearest integer to x: [x]
    using the convention that 
    [1/2]=0 and [-1/2]=0


    EXAMPLES::

        sage: nearest_integer(0)
        0
        sage: nearest_integer(0.5)
        1
        sage: nearest_integer(-0.5)
        0
        sage: nearest_integer(-0.500001)
        -1
    """
    return floor(x+<double>0.5)

