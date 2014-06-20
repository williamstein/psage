# cython: profile=False
# -*- coding=utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2013 Fredrik Strömberg <stroemberg@mathematik.tu-darmstadt.de>
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
Cython class for for working with H^n, i.e. n copies of upper half-plane.


"""
include 'sage/ext/stdsage.pxi'
#include "sage/ext/cdefs.pxi"
include 'sage/ext/interrupt.pxi'
include "../../rings/double_prec_math.pxi"


cdef extern from "complex.h":
    ### Only in new versions of gcc....
    ### cdef complex CMPLX(double,double)
    cdef double complex _Complex_I
    
cdef double complex CMPLX(double x,double y):
    return x+y*_Complex_I
    

from sage.all import real,imag,Integer,Rational,ComplexField,Infinity
from sage.rings.complex_number import is_ComplexNumber
from sage.rings.real_mpfr import is_RealNumber
from sage.groups.perm_gps.permgroup_element import is_PermutationGroupElement
from sage.rings.number_field.number_field_ideal import is_NumberFieldFractionalIdeal as is_NFIdeal
from sage.rings.number_field.number_field_element import is_NumberFieldElement 


cdef class  Hn(object):
    r"""
    Class of elements in products of complex upper  half-planes
    with additional ring structure given by:
    - component-wise multiplication and division
    - multiplication and division by elements of a number field of the same degree as the dimension.

    TODO: Inherit from Ring or something? (for speed probably NO!)

    """
    def __cinit__(self,x,y=[],degree=0,prec=53,verbose=0):
        r"""
        Init self from a list or an algebraic point.
        Currently we only work with double (53 bits) precision.
        """
        self._verbose = verbose
        u = []
        if isinstance(x,list) and isinstance(y,list) and y<>[] and x<>[]:
            degree = len(x)
            assert len(y) == degree
        elif not isinstance(x,list):
            y=[]
            if hasattr(x,'complex_embeddings'):
                v = x.vector()  ## split the real and imaginary parts
                v0 = v[0].complex_embeddings()
                v1 = v[1].complex_embeddings()
                sqD = abs(x.parent().power_basis()[1].complex_embedding()) ## Recall that sqD is purely imaginary
                degree = max(len(v0),len(v1))
                for i in range(degree):
                    u.append(v0[i])
                    if v1[i]>0:
                        y.append(v1[i]*sqD)
                    else:
                        y.append(-v1[i]*sqD)
            elif isinstance(x,Hn):
                u = x.real_list(); y = x.imag_list()
                self._prec = x.prec()
                degree = x.degree()
            elif hasattr(x,"Hn"):
                u = x.Hn().real_list(); y = x.Hn().imag_list()
                self._prec = x.Hn().prec()
                degree = x.Hn().degree()
            else:
                raise TypeError,"Can not construct point in H^n from %s!" % x
            if verbose>0:
                print "u=",u
                print "v=",y
        else:
            ## WE now have to see if we are called with one of the three possibilities:
            ## v=[xlist,ylist]=[[x[0],x[1],...,],[y[0],y[1],...]]
            y = []
            if isinstance(x[0],list):
                if verbose>0:
                    print "x[0] is list"
                if(len(x))<>2:
                    raise ValueError
                u = x[0]; y = x[1]
                degree = len(x)
            ## v=zlist=[z[0],z[1],...,]
            ## Fast cases first
            elif isinstance(x[0],complex):
                if verbose>0:
                    print "x[0] is complex"
                degree = len(x)
                for i in range(degree):
                    u.append(x[i].real)
                    y.append(x[i].imag)
            elif is_ComplexNumber(x[0]):
                if verbose>0:
                    print "x[0] is ComplexNumber"
                degree = len(x)
                for i in range(degree):
                    u.append(x[i].real())
                    y.append(x[i].imag())
            elif isinstance(x[0],float) or is_RealNumber(x[0]):
                degree = len(x)
                for i in range(degree):
                    u.append(x[i])
                    y.append(0)
            elif isinstance(x,list):
                degree = len(x)
                try: 
                    for i in range(degree):
                        u.append(real(x[i]))
                        y.append(imag(x[i]))
                    self._degree = len(x)
                except:
                    raise TypeError,"Can not construct point in H^n from %s!" % x
            else:
                raise TypeError,"Can not construct point in H^n from %s!" % x
        if u<>[]:
            self._xlist = u
        else:
            self._xlist = x
        self._ylist = y
        self._prec = 53  ## Nothing else implemented for now
        if verbose>0:
            print "Set x=",self._xlist
            print "Set y=",self._ylist
            print "set degree=",degree
        self._degree = degree
        assert len(self._xlist)==self._degree and len(self._ylist)==self._degree
        ## norms and indicators
        self._imag_norm = 0.0
        self._real_norm = 0.0
        self._norm = 0.0
        self._imag_norm_set = 0
        self._real_norm_set = 0
        self._norm_set = 0
        self.c_new(self._xlist,self._ylist)
        if verbose>0:
            print "c_new successful!"
        

    def __init__(self,x,y=[],degree=0,prec=53,verbose=0):
        pass

    def __reduce__(self):
        return (Hn,(self._xlist,self._ylist,self._degree,self._prec,self._verbose))
        
    cdef c_new(self,list x,list y):
        self._x = NULL; self._y=NULL
        self._x = <double*>sage_malloc(sizeof(double)*self._degree)
        if self._x==NULL:
            raise MemoryError
        self._y = <double*>sage_malloc(sizeof(double)*self._degree)
        if self._y==NULL:
            raise MemoryError
        cdef int i
        for i in range(self._degree):
            self._x[i] = <double>x[i]
            self._y[i] = <double>y[i]
            if self._verbose>1:
                print "x[{0}]={1}".format(i,x[i])
                print "y[{0}]={1}".format(i,y[i])
            if y[i]<0:
                raise ValueError,"Not in H^n^*! y[{0}]={1}".format(i,y[i])
        if self._verbose>0:
            print "allocated x and y!"


    def __richcmp__(self, right, int op):
        res=1
        if op <>2 and op <>3:
            raise NotImplementedError,"Ordering of points in H^n is not implemented!"
        if type(self)<>type(right):
            res=0
        elif right.degree() <> self.degree():
            res=0
        else:
            for i in range(self.degree()):
                if self.x(i)<>right.x(i) or self.y(i)<>right.y(i):
                    res = 0
                    break
            
        if op==3: # != # op=2 is ==
            if res==0:
                return True
            else:
                return False
        if res==0:
            return False
        else:
            return True

    def __copy__(self):
        x = self.real_list()
        y = self.imag_list()
        return Hn(x,y,self._degree,self._prec,self._verbose)
        
    def __repr__(self):
        return str(list(self))

    def __getitem__(self,i):        
        if isinstance(i,(int,Integer)) and i>=0 and i<self._degree:
            return self.z(i)
        else:
            raise IndexError
        
    def __dealloc__(self):
        self._free()
            
    cdef _free(self):
        if self._x<>NULL:
            sage_free(self._x)
        if self._y<>NULL:
            sage_free(self._y)
        self._degree = 0        

    def __add__(self,other):
        r"""
        Add two points in the upper half-plane produce another point
        """
        if not is_Hn(other):
            other = Hn(other)
        x = self.real_list(); y = self.imag_list()
        assert self.degree()==other.degree()
        for i in range(self.degree()):
            x[i]+=other.x(i)
            y[i]+=other.y(i)
        return Hn(x,y)

    def __sub__(self,other):
        r"""
        Substract two points in the upper half-plane may sometimes produce another point.
        """
        if not is_Hn(other):
            other = Hn(other)
        x = self.real_list(); y = self.imag_list()
        assert self.degree()==other.degree()
        for i in range(self.degree()):
            x[i]-=other.x(i)
            y[i]-=other.y(i)
        return Hn(x,y)
                
    def __imag__(self):
        return self.imag()

    def __real__(self):
        return self.real()
    
    def __len__(self):
        return self.degree()


    def addto(self,other):
        #x = self.real_list(); y = self.imag_list()
        if hasattr(other,'complex_embeddings'):
            c = other.complex_embeddings(self._prec)
            assert len(c)==self.degree()
            for i in range(self.degree()):
                self._x[i]+=c[i].real()
                self._y[i]+=c[i].imag()
        elif hasattr(other,'_xlist'):
            assert self.degree()==other.degree()
            for i in range(self.degree()):
                self._x[i]+=other.x(i)
                self._y[i]+=other.y(i)
        elif is_ComplexNumber(other):
            for i in range(self.degree()):
                self._x[i]+=other
        else:
            raise NotImplementedError,"Can not add Hn and %s of type %s!" %(other,type(other))
        for i in range(self.degree()):
            self._xlist[i]=self._x[i]
            self._ylist[i]=self._y[i]
            
        #return Hn(res)

    def prec(self):
        return self._prec

    cpdef imag_norm(self):
        cdef int i
        if self._imag_norm_set==0:
            self._imag_norm=1.0
            for i in range(self._degree):
                self._imag_norm = self._imag_norm*self._y[i]
            self._imag_norm_set=1
        return self._imag_norm


    cpdef real_norm(self):        
        cdef int i
        if self._real_norm_set==0:
            self._real_norm=1.0
            for i in range(self._degree):
                self._real_norm = self._real_norm*self._x[i]
            self._real_norm_set=1
        return self._real_norm

    cpdef imag_list(self):
        return self._ylist


    cpdef real_list(self):
        return self._xlist

    cpdef norm(self):
        cdef int i
        if self._norm_set==0:
            self._norm = CMPLX(1.0,0.0) #<double complex>1.0
            for i in range(self._degree):
                self._norm = self._norm*(self._x[i]+_I*self._y[i])
            self._norm_set=1
        return self._norm 

    cpdef vector_norm(self):
        r"""
        Return the Euclidean norm of self as a vector in C^n
        """
        cdef double res = 1.0
        cdef int i
        for i in range(self._degree):
            res = res+(cabs(self._x[i]+_I*self._y[i]))**2
        return sqrt(res)

    cpdef trace(self):
        cdef double complex res
        res = 1.0
        cdef int i
        for i in range(self._degree):
            res = res+self._x[i]+_I*self._y[i]
        return res

    cpdef real_trace(self):
        cdef double res = 1.0
        cdef int i
        for i in range(self._degree):
            res = res+self._x[i]
        return res

    cpdef imag_trace(self):
        cdef double res = 1.0
        cdef int i
        for i in range(self._degree):
            res = res+self._y[i]
        return res
    

    cpdef degree(self):
        return self._degree
    
    
    cpdef x(self,int i):
        if i>=0 and i<self._degree:
            return self._x[i]
        else:
            raise IndexError

    cpdef y(self,int i):
        if i>=0 and i<self._degree:
            return self._y[i]
        else:
            raise IndexError

    cpdef z(self,int i):
        if i>=0 and i<self._degree:
            return self._x[i]+_I*self._y[i]
        else:
            raise IndexError
        
    cpdef imag(self):        
        return Hn(self._ylist)

    cpdef real(self):
        return Hn(self._xlist)

    def __list__(self):
        res = []
        for i in range(self.degree()):
            res.append(self.z(i))
        return res

####  HERE
    
    def permute(self,s):
        r"""
        The permutation group acts by permuting the components.
        """
        assert is_PermutationGroupElement(s)
        znew = [0 for i in range(self._degree)]
        for i in range(self._degree):
            si = s.list()[i]-1
            znew[si]=self.z(i)
        return Hn(znew)




    
    def acton_by(self,A):
        r"""
        act on self by the matrix A in SL(2,K)
        """
        assert A.det()==1
        res=[]
        prec = self._prec
        if isinstance(A[0,0],Rational):
            a = [A[0,0] for i in range(self._degree)]
            b = [A[0,1] for i in range(self._degree)]
            c = [A[1,0] for i in range(self._degree)]
            d = [A[1,1] for i in range(self._degree)]
        elif is_NumberFieldElement(A[0,0]):
            a = A[0,0].complex_embeddings(prec)
            b = A[0,1].complex_embeddings(prec)
            c = A[1,0].complex_embeddings(prec)
            d = A[1,1].complex_embeddings(prec)
        for i in range(self._degree):
            den = self.z(i)
            den = den*c[i]+d[i]
            if den<>0:
                w = (self.z(i)*a[i]+b[i])/den
            else:
                w = Infinity
            res.append(w)
        return Hn(res)





    
    def parent(self):
        return None

        
        
    cpdef is_in_upper_half_plane(self):
        r"""
        Returns true if all components are in the upper half-plane.
        """
        cdef double y
        for i in range(self._degree):
            y = self._y[i]
            if y<0:
                return False
        return True
    
    #def log(self):
    #    res = map(log,list(self))
    #    return Hn(res)
##    HERE 2
    def hyp_dist(self,w,dtype=1):
        r"""
        If self and w are in the upper half-plane then
        return d(self.z(i),z.z(i))
        dtype = 0 => dist = dist1+dist2
        dtype = 1 => dist = max(dist1,dist2)
        TODO: FASTER VERSION...
        """
        if dtype == 1:
            maxd = 0
        if hasattr(w,'_xlist') and w._degree == self._degree:
            if w.z == self.z:
                return 0
            distances=[]
            for i in range(self.degre()):
                ab1 =abs(self.z(i)-w.z(i))
                ab2 =abs(self.z(i)-ComplexField(self._prec)(w.x(i),-w.y(i)))
                l = log((ab1+ab2)/(ab2-ab1))
                distances.append(l)

            if dtype==0:
                return sum(distances)
            else:
                return max(distances)

    def reflect(self):
        r"""
        z -> -conjugate(z)
        """
        for i in range(self.degree()):
            self._x[i]=-self._x[i]
            self._xlist[i]=-self._xlist[i]

    cpdef diff_abs_norm(self,z):
        r"""
        Return the euclidean norm of the difference of self and z as vectors in C^n: ||self-z||
        """
        cdef double norm = 0.0
        if hasattr(z,'complex_embeddings'):
            c = z.complex_embeddings(self._prec)
            for i in range(self.degree()):
                norm += abs(self.z(i)-c[i])**2
        elif is_Hn(z):
            for i in range(self.degree()):
                norm += abs(self.z(i)-z.z(i))**2
        elif not isinstance(z,list):  ## Assume z is a complex number...
            try:
                for i in range(self.degree()):
                    norm += abs(self.z(i)-z)**2
            except:
                raise NotImplementedError,"Can not substract point in Hn and %s of type %s!" %(z,type(z))
        else:
            raise NotImplementedError,"Can not substract point in Hn and %s of type %s!" %(z,type(z))
        return sqrt(norm)




    ### Below this point some algorithms might be inapproprate leftovers... 
    ### 
    cpdef addto_re(self,x):
        r"""
        Add the real x to self
        """
        if hasattr(x,'complex_embeddings'):
            c = x.complex_embeddings(self._prec)
            for i in range(self.degree()):
                self._x[i]+=c[i]
                self._xlist[i]+=c[i]
        elif is_Hn(x):
            for i in range(self.degree()):
                self._x[i]+=x._x[i]
                self._xlist[i]+=x._xlist[i]
        elif is_RealNumber(x):
            for i in range(self.degree()):
                self._x[i]+=x
                self._xlist[i]+=x
        else:
            raise NotImplementedError,"Can not add Hn and %s of type %s!" %(x,type(x))
        
    def __div__(self,other):
        if hasattr(other,'_is_Hn') and self._degree == other._degree:            
            res = [self.z(i)/other.z(i) for i in xrange(self._degree) ]
        elif hasattr(other,'complex_embeddings') and other.parent().degree()==self._degree:
            w = other.complex_embeddings(self._prec)
            res = [self.z(i)/w[i] for i in range(self.degree()) ]
        elif isinstance(other,list) and len(other)==self.degree():
                res = [self.z(i)/other[i] for i in range(self.degree()) ]
        else:
            raise NotImplementedError,"Can not divide Hn by %s of type %s!" %(other,type(other))
        return Hn(res)


    def __rmul__(self,other):
        if isinstance(other,type(self)):
            assert self._degree == other._degree
            res = [self.z(i)*other.z(i) for i in xrange(self._degree) ]
        elif hasattr(other,'complex_embeddings'):
            w = other.complex_embeddings(self._prec)
            if len(w) == self._degree:
                res = [w[i]*self.z(i) for i in xrange(self._degree) ]
        elif isinstance(other,list) and len(other)==self._degree:
                res = [self.z(i)*other[i] for i in xrange(self._degree) ]
        else:
            res = [self.z(i)*other for i in xrange(self._degree) ]
        return Hn(res)

    def __lmul__(self,other):
        if isinstance(other,type(self)):
            assert self.degree() == other.degree()
            res = [self.z(i)*other.z(i) for i in range(self.degree()) ]
        elif hasattr(other,'complex_embeddings'):
            w = other.complex_embeddings(self._prec)
            assert len(w) == self.degree()
            res = [w[i]*self.z(i) for i in range(self.degree()) ]
        elif isinstance(other,list) and len(other)==self.degree():
            res = [self.z(i)*other[i] for i in range(self.degree()) ]
        else:
            res = [self.z(i)*other for i in range(self.degree()) ]
        return Hn(res)


    def __mul__(self,other):
        if isinstance(other,type(self)):
            assert self.degree() == other._degree
            res = [self.z(i)*other.z(i) for i in range(self.degree()) ]
        elif hasattr(other,'complex_embeddings'):
            w = other.complex_embeddings(self.prec())
            assert len(w) == self.degree()
            res = [w[i]*self.z(i) for i in range(self.degree()) ]
        elif isinstance(other,list) and len(other)==self.degree():
            res = [self.z(i)*other[i] for i in range(self.degree()) ]
        else:
            res = [self.z(i)*other for i in range(self.degree()) ]
        return Hn(res)
 
        

cpdef is_Hn(z):
    if hasattr(z,"real_list"):
        return 1
    return 0
