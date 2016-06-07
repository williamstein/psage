"""
Vectors with mpfr entries. 

AUTHOR:
    -- Fredrik Stroemberg (2011)

NOTE: based on vector_complex_dense.pyx which in turn was copied from sage/modules/vector_rational_dense.pyx
and


EXAMPLES:


TESTS:

"""

###############################################################################
#   Sage: System for Algebra and Geometry Experimentation    
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###############################################################################

include 'cysignals/signals.pxi'
include '../ext/stdsage.pxi'


# set rounding to be nearest integer
# TODO: make t possible to change rounding 
cdef mpfr_rnd_t rnd
rnd = GMP_RNDN

from sage.structure.element cimport Element, ModuleElement, RingElement, Vector
from sage.all import FreeModule

from sage.rings.integer cimport Integer

#from sage.rings.rational cimport Rational
#from sage.rings.complex_mpc cimport MPComplexNumber,MPComplexField
#from sage.rings.complex_mpc cimport MPComplexNumber
#from sage.rings.complex_mpc cimport MPComplexField_class
from sage.modules.free_module_element cimport FreeModuleElement
#from psage.matrix.matrix_complex_dense cimport Matrix_complex_dense
from sage.rings.real_mpfr cimport RealNumber,RealField_class

#from free_module_element import vector

cdef class Vector_real_mpfr_dense(FreeModuleElement):
    cdef bint is_dense_c(self):
        return 1
    cdef bint is_sparse_c(self):
        return 0

    cdef _new_c(self,v=None):
        cdef Vector_real_mpfr_dense y
        #print "in new_c! v=",v
        #y = PY_NEW(Vector_complex_dense)
        y = Vector_real_mpfr_dense.__new__(Vector_real_mpfr_dense,self._parent,v,False,False)
        return y

    def __copy__(self):
        cdef Vector_real_mpfr_dense y
        y = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            mpfr_init2(y._entries[i],self._prec)
            mpfr_set(y._entries[i], self._entries[i],rnd)
        return y

    cdef _init(self, Py_ssize_t degree, parent):
        self._degree = degree
        self._parent = parent
        self._base_ring=parent._base
        self._prec = self._base_ring.__prec
        self._entries = <mpfr_t *> sage_malloc(sizeof(mpfr_t) * degree)
        #print "_inite entries=",<int>self._entries
        #print "_inite pprec=",self._prec
        if self._entries == NULL:
            raise MemoryError
        for i from 0 <= i < self._degree:
            mpfr_init2(self._entries[i],self._prec)
        
    def __cinit__(self, parent=None, x=None, coerce=True,copy=True):
        self._entries = NULL
        self._is_mutable = 1
        #print "in __cinit__"
        if not parent is None:
            self._init(parent.degree(), parent)
        #print "in __cinit__ entries=",printf("%p",self._entries)


        
    def __init__(self, parent, x, coerce=True, copy=True):
        cdef Py_ssize_t i
        cdef RealNumber z
        # we have to do this to avoid a garbage collection error in dealloc
        #print "in init parent=",parent
        #print "in init x.parent=",x.parent

        #print "in init x=",x
        #
        #for i from 0 <= i < self._degree:
        #    mpfr_init2(self._entries[i],self._prec)
        #    print "init ",i
        if isinstance(x, (list, tuple)) or (hasattr(x,'parent') and x.parent() == parent) :
            if len(x) != self._degree:
                raise TypeError, "entries must be a list of length %s"%self._degree
            for i from 0 <= i < self._degree:
                z = RealNumber(self._base_ring,x[i])
                mpfr_set(self._entries[i], z.value,rnd)
            return
        elif x != 0:
            print "type(x)=",type(x)
            print "z=",RealNumber(self._base_ring,x[0])
            raise TypeError, "can't initialize vector from nonzero non-list"
        else:
            for i from 0 <= i < self._degree:
                mpfr_set_ui(self._entries[i], 0,rnd)
                
    def __dealloc__(self):
        cdef Py_ssize_t i
        if self._entries<>NULL:
            sig_on()
            for i from 0 <= i < self._degree:
                #print "clearing mpfr's entry %s"%i

                mpfr_clear(self._entries[i])
            sig_off()
            #print "clearing python entries"
            sage_free(self._entries)

    cpdef base_ring(self):
        return self._base_ring

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        EXAMPLES:
            sage: v = vector(QQ, [0,0,0,0])
            sage: v == 0
            True
            sage: v == 1
            False
            sage: v == v
            True
            sage: w = vector(QQ, [-1,3/2,0,0])
            sage: w < v
            True
            sage: w > v
            False
        """
        cdef Py_ssize_t i
        cdef int c
        for i from 0 <= i < left.degree():
            c = mpfr_cmp(left._entries[i], (<Vector_real_mpfr_dense>right)._entries[i])
            if c < 0:
                return -1
            elif c > 0:
                return 1
        return 0

    def __len__(self):
        return self._degree

    def __setitem__(self, Py_ssize_t i, x):
        if not self._is_mutable:
            raise ValueError, "vector is immutable; please change a copy instead (use copy())"
        cdef RealNumber z
        if i < 0 or i >= self._degree:
            raise IndexError
        else:
            z = RealNumber(self._base_ring,x)
            mpfr_set(self._entries[i], z.value,rnd)
            
    def __getitem__(self, Py_ssize_t i):
        """
        Return the ith entry of self.

        EXAMPLES:
            sage: v = vector([1/2,2/3,3/4]); v
            (1/2, 2/3, 3/4)
            sage: v[0]
            1/2
            sage: v[2]
            3/4
            sage: v[-2]
            2/3
            sage: v[5]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
            sage: v[-5]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        cdef RealNumber z
        #z = PY_NEW(MPComplexNumber)
        z = RealNumber(self._base_ring,0)
        if i < 0:
            i += self._degree
        if i < 0 or i >= self._degree:
            raise IndexError, 'index out of range'
        else:
            mpfr_set(z.value, self._entries[i],rnd)
            return z

    def __reduce__(self):
        return (unpickle_v1, (self._parent, self.list(), self._degree, self._is_mutable))

    cpdef ModuleElement _add_(self, ModuleElement right):
        cdef Vector_real_mpfr_dense z, r
        r = right
        #print "in add!"
        z = self._new_c()
        cdef Py_ssize_t i
        #prec=self.parent().base_ring().prec()
        for i from 0 <= i < self._degree:
            mpfr_init2(z._entries[i],self._prec)
            mpfr_add(z._entries[i], self._entries[i], r._entries[i],rnd)
        return z
        

    cpdef ModuleElement _sub_(self, ModuleElement right):
        cdef Vector_real_mpfr_dense z, r
        r = right
        #print "in sub!"
        z = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            mpfr_init2(z._entries[i],self._prec)
            mpfr_sub(z._entries[i], self._entries[i], r._entries[i],rnd)
        return z
        
    cpdef Element _dot_product_(self, Vector right):
        """
        Dot product of dense vectors over mpfr.
        
        EXAMPLES:
            sage: v = vector(QQ, [1,2,-3]); w = vector(QQ,[4,3,2])
            sage: v*w
            4
            sage: w*v
            4
        """
        cdef Vector_real_mpfr_dense r = right
        cdef RealNumber z
        z = RealNumber(self._base_ring,0)
        cdef mpfr_t t
        mpfr_init2(t,self._prec)
        mpfr_set_si(z.value, 0,rnd)
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            mpfr_mul(t, self._entries[i], r._entries[i],rnd)
            mpfr_add(z.value, z.value, t,rnd)
        mpfr_clear(t)
        return z
    
    def scalar_product(self,right):
        return self._scalar_product_(right)

    cdef RealNumber _scalar_product_(self, Vector right):
        """
        Hermitian scalar product of dense vectors over mpfr.

        Note: conjugation on the second component.
        EXAMPLES:

        """
        cdef Vector_real_mpfr_dense r = right
        cdef RealNumber z
        z = RealNumber(self._base_ring,0)
        cdef mpfr_t t
        mpfr_init2(t,self._prec)
        mpfr_set_si(z.value, 0,rnd)
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            #mpfr_conj(x,r._entries[i],rnd)
            mpfr_mul(t, self._entries[i], r._entries[i],rnd)
            mpfr_add(z.value, z.value, t,rnd)
        mpfr_clear(t)
        return z
        

    cpdef Vector _pairwise_product_(self, Vector right):
        """
        EXAMPLES:
            sage: v = vector(QQ, [1,2,-3]); w = vector(QQ,[4,3,2])
            sage: v.pairwise_product(w)
            (4, 6, -6)
        """
        print "in pwp"
        cdef Vector_real_mpfr_dense z, r
        r = right
        z = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            mpfr_init2(z._entries[i],self._prec)
            mpfr_mul(z._entries[i], self._entries[i], r._entries[i],rnd)
        return z

        
    def __mul__(self, ModuleElement right):
         #cdef Vector_real_mpfr_dense z
         cdef RealNumber a
         if isinstance(right,Vector_real_mpfr_dense):
             return self._dot_product_(right)
         if isinstance(right,RealNumber):
             return self._lmul_(right)
         raise NotImplemented,"Real matrix times vector is not implemented!"
         #if isinstance(right,Matrix_complex_dense):
         #   return Matrix_complex_dense._vector_times_matrix_(right,self)

        
    cpdef ModuleElement _rmul_(self, RingElement left):
        cdef Vector_real_mpfr_dense z
        cdef Py_ssize_t i
        cdef RealNumber a
        z = self._new_c()
        # we can convert almost anything to MPComplexNumber
        if not isinstance(left,RealNumber):
            a = RealNumber(self._base_ring,left)
            for i from 0 <= i < self._degree:
                mpfr_init2(z._entries[i],self._prec)
                mpfr_mul(z._entries[i], self._entries[i], a.value,rnd)
        else:
            for i from 0 <= i < self._degree:
                mpfr_init2(z._entries[i],self._prec)
                mpfr_mul(z._entries[i], self._entries[i], (<RealNumber>(left)).value,rnd)
        return z


    cpdef ModuleElement _lmul_(self, RingElement right):
        cdef Vector_real_mpfr_dense z
        cdef RealNumber a
        # we can convert almost anything to MPComplexNumber
        if not isinstance(right,RealNumber):
             a = RealNumber(self._base_ring,right)
        else:
            a = right
        z = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            mpfr_init2(z._entries[i],self._prec)
            mpfr_mul(z._entries[i], self._entries[i], a.value,rnd)
        return z

    cpdef ModuleElement _neg_(self):
        cdef Vector_real_mpfr_dense z
        z = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            mpfr_init2(z._entries[i],self._prec)
            mpfr_neg(z._entries[i], self._entries[i],rnd)
        return z


    cpdef RealNumber norm(self,int ntype=2):
        r"""
        The norm of self.
        """
        cdef RealNumber res
        cdef mpfr_t x,s
        mpfr_init2(x,self._prec)
        mpfr_init2(s,self._prec)
        res = RealNumber(self._base_ring._base,0)
        if ntype==2:
            mpfr_set_si(s, 0, rnd)
            for i from 0 <= i < self._degree:
                mpfr_abs(x,self._entries[i],rnd)
                mpfr_sqr(x,x,rnd)
                mpfr_add(s,s,x,rnd)
            mpfr_sqrt(s,s,rnd)
            mpfr_set(res.value,s,rnd)
            mpfr_clear(x); mpfr_clear(s)
            return res
        else:
            mpfr_clear(x); mpfr_clear(s)
            raise NotImplementedError,"Only 2-norm is currently implemented"

    def set_prec(self,int prec):
        r"""
        Change the defualt precision for self.
        """        
        cdef Py_ssize_t i,j
        cdef mpfr_t z
        from sage.rings.complex_mpc import MPComplexField
        from sage.matrix.all import MatrixSpace
        mpfr_init2(z,prec)
        self._prec = prec
        self._base_ring=MPComplexField(prec)
        self._parent = FreeModule(self._base_ring,self._degree)
        #self._parent = FreeModule(self._base_ring,self._nrows,self._ncols)
        for i from 0 <= i < self._degree:
                mpfr_set(z,self._entries[i],rnd)
                mpfr_set_prec(self._entries[i],prec)
                mpfr_set(self._entries[i],z,rnd)
        #RF= self._base_ring._base
        #self._eps = RF(2)**RF(5-self._prec)
        #mpfr_init2(self._mpeps,self._prec)
        #mpfr_set(self._mpeps,self._eps.value,self._rnd_re)
        mpfr_clear(z)
        

    def n(self, *args, **kwargs):
        """
        Returns a numerical approximation of self by calling the n()
        method on all of its entries.

        EXAMPLES:
            sage: v = vector(QQ, [1,2,3])
            sage: v.n()
            (1.00000000000000, 2.00000000000000, 3.00000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 53 bits of precision
            sage: v.n(prec=75)
            (1.000000000000000000000, 2.000000000000000000000, 3.000000000000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 75 bits of precision
        """
        return Vector( [e.n(*args, **kwargs) for e in self] )


def make_FreeModuleElement_complex_dense_v1(parent, entries, degree):
    # If you think you want to change this function, don't.
    # Instead make a new version with a name like
    #    make_FreeModuleElement_generic_dense_v1
    # and changed the reduce method below.
    cdef Vector_real_mpfr_dense v
    v = PY_NEW(Vector_real_mpfr_dense)
    v._init(degree, parent)
    prec=parent.prec()
    cdef RealNumber z
    base_ring=parent.base_ring()
    for i from 0 <= i < degree:
        z = RealNumber(base_ring,entries[i])
        mpfr_init2(v._entries[i],prec)
        mpfr_set(v._entries[i], z.value,rnd)
    return v

def unpickle_v1(parent, entries, degree, is_mutable):
    cdef Vector_real_mpfr_dense v
    v = PY_NEW(Vector_real_mpfr_dense)
    v._init(degree, parent)
    prec=parent.prec()
    cdef RealNumber z
    base_ring=parent.base_ring()
    for i from 0 <= i < degree:
        z = RealNumber(base_ring,entries[i])
        mpfr_init2(v._entries[i],prec)
        mpfr_set(v._entries[i], z.value,rnd)
    v._is_mutable = is_mutable
    return v
