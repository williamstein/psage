# encoding: utf-8
# cython: profile=True
# filename: matrix_complex_dense.pyx

r"""
Dense matrices over the complex field.

EXAMPLES
"""
include "stdsage.pxi"
include "cysignals/signals.pxi"

from psage.rings.mp_cimports cimport *


# set rounding to be nearest integer
# TODO: make t possible to change rounding


from sage.rings.real_mpfr import RealField as RFF
from sage.matrix.matrix cimport Matrix
from sage.matrix.matrix_dense cimport Matrix_dense
from sage.rings.complex_mpc import MPComplexField

from psage.modules.vector_complex_dense cimport Vector_complex_dense
from psage.matrix.linalg_complex_dense cimport _eigenvalues,qr_decomp,_norm,_hessenberg_reduction,init_QR,QR_set
from psage.rings.mpc_extras cimport *
from psage.matrix.linalg_complex_dense cimport qr_decomp,_reconstruct_matrix,solve_upper_triangular
from psage.modules.vector_complex_dense cimport Vector_complex_dense

from sage.structure.element cimport ModuleElement, RingElement, Element, Vector

from sage.all import MatrixSpace
from  sage.modules.free_module_element import vector
from sage.all import FreeModule
from sage.functions.other import ceil
from sage.rings.ring import is_Ring
from sage.rings.rational_field import QQ
from sage.rings.complex_double import CDF
from sage.arith.all import gcd,valuation
from sage.matrix.matrix import is_Matrix
from sage.structure.element import is_Vector
from sage.matrix.matrix2 import cmp_pivots, decomp_seq
from sage.matrix.matrix0 import Matrix as Matrix_base

from sage.misc.misc import verbose, get_verbose
from sage.misc.all import prod

## #########################################################

    
#########################################################

cdef class Matrix_complex_dense(Matrix_dense):  

    ########################################################################
    # LEVEL 1 functionality
    # x * __cinit__  
    # x * __dealloc   
    # x * __init__      
    # x * set_unsafe
    # x * get_unsafe
    # x * cdef _pickle
    # x * cdef _unpickle
    ########################################################################

    def __cinit__(self, parent,entries=0, coerce=True, copy=True,verbose=0):
        """
        Create and allocate memory for the matrix.
        
        Unlike over matrix_integer_dense, mpc_init() is called (as there is no mpc_init_set function).
        
        INPUT:
        
        
        -  ``parent, entries, coerce, copy`` - as for
           __init__.
        
        
        EXAMPLES::
        
            sage: from sage.matrix.matrix_rational_dense import Matrix_rational_dense
            sage: a = Matrix_rational_dense.__new__(Matrix_rational_dense, Mat(ZZ,3), 0,0,0)
            sage: type(a)
            <type 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>
        
        .. warning::

           This is for internal use only, or if you really know what
           you're doing.
        """
        # This is called before __init__
        self._entries_are_allocated = 0
        self._double_matrix_is_set = 0
        self._verbose=int(verbose)
        if self._verbose>0:
            print "in cinit"
        if self._verbose>1:
            print "entries=",entries,parent
        cdef double x
#        if not isinstance(parent,sage.matrix.matrix_space.MatrixSpace):
        if not hasattr(parent,'base_ring') or not hasattr(parent,'nrows'):
            if verbose>0:
                print "no matrix space!"
            raise ValueError,"Need MatrixSpace as parent!" 

        #Matrix_dense.__init__(self, parent) #,None,copy,coerce) #,entries=entries) #coerce=coerce,copy=copy)
        if self._verbose>0:
            print "before base ring"
        self._parent = parent
        if self._verbose>0:
            print "before base ring"
        self._base_ring = parent.base_ring()
        if self._verbose>0:
            print "after base ring"
        self._nrows = parent.nrows()
        self._ncols = parent.ncols()
        if self._verbose>0:
            print "matrix dense inited!"
        cdef Py_ssize_t i, k
        self._entries = <mpc_t *> sage_malloc(sizeof(mpc_t)*(self._nrows * self._ncols))
        if self._verbose>0:
            print "entries alloc!"
        if self._entries == NULL:
            raise MemoryError, "Out of memory allocating entries for a matrix of size {0} x {1}".format(self._nrows,self._ncols)
        sig_on()        
        self._matrix =  <mpc_t **> sage_malloc(sizeof(mpc_t*) * self._nrows)
        sig_off()
        if self._verbose>0:
            print "matrix alloc!"
        if self._matrix == NULL:
            sage_free(self._entries)
            self._entries = NULL
            raise MemoryError, "Out of memory allocating a matrix of size {0} x {1}".format(self._nrows,self._ncols)
        # store pointers to the starts of the rows

        k = 0
        for i from 0 <= i < self._nrows:
            self._matrix[i] = self._entries + k
            k = k + self._ncols
        #cdef int prec
        if self._verbose>0:
            print "here0: base = ",self._parent._base
        #prec=self.parent().base_ring().prec()
        self._prec= self._base_ring.prec()
        self._base_ring = MPComplexField(self._prec)

        ## TODO: acces rounding modes of parents
        mpfr_init2(self._mpeps,self._prec)
        self._rnd=MPC_RNDNN ## self._base_ring.__rnd
        #print "here! rnd=",self._rnd
        self._rnd_re=MPFR_RNDN #self._base_ring.__real_field.__rnd
        #print "here! rnd_re=",self._rnd_re
        self._rnd_im=MPFR_RNDN # self._base_ring.__imag_field.__rnd
        #print "here! rnd_im=",self._rnd_im
        ## we also set something like an "effective" machine epsilon
        self._base_for_str_rep=32
        self._truncate=0
        RF=self._base_ring._base
        ## This is an "efficient" epsilon for self.
        ## i.e. we allow for sloppy implementations...
        ##
        self._eps = RF(2)**RF(1-self._prec)
        mpfr_set(self._mpeps,self._eps.value,self._rnd_re)
        if self._nrows<>self._ncols:
            self._is_square = 0
        else:
            self._is_square = 1

        self._transformation_to_hessenberg=NULL
        self._error_qr=RF(0)

        for i from 0 <= i < self._nrows * self._ncols:
            mpc_init2(self._entries[i],self._prec)
        self._entries_are_allocated = 1
        if self._verbose>0:
            print "allocated {0}".format(self._nrows * self._ncols)

    def  __dealloc__(self):
        mpfr_clear(self._mpeps)
        if self._entries == NULL:
            return 
        cdef Py_ssize_t i
        if self._entries_are_allocated == 1:
            for i from 0 <= i < self._nrows * self._ncols:
                if self._verbose>1:
                    print "dealloc i=",i
                mpc_clear(self._entries[i])
        sage_free(self._entries)
        if self._matrix<> NULL:
            sage_free(self._matrix)

    
    def __init__(self, parent, entries=0, coerce=True, copy=True,verbose=0):
    
        cdef Py_ssize_t i,j
        cdef MPComplexNumber z
        if self._verbose>0:
            print "in init"
        self._norms=dict()
        self._double_matrix=None
        self._double_matrix_is_set=0
        if isinstance(entries, (list, tuple)):
            if len(entries) != self._nrows * self._ncols:
                raise TypeError, "entries has the wrong length"
            sig_on()
            if coerce:
                #print  "len=",self._nrows * self._ncols
                for i from 0 <= i < self._nrows * self._ncols:
                    # TODO: Should use an unsafe un-bounds-checked array access here.
                    z = self._base_ring(entries[i])
                    mpc_set(self._entries[i], z.value, self._rnd)
            else: 
                for i from 0 <= i < self._nrows * self._ncols:
                    # TODO: Should use an unsafe un-bounds-checked array access here.
                    z = entries[i]
                    mpc_set(self._entries[i], z.value, self._self._rnd)
            sig_off()
        else: 
            # is it a scalar?
            try:
                # Try to coerce entries to a scalar 
                if self._verbose>0:
                    print "entries=",entries,type(entries)
                    print "base_ring=",self.base_ring()
                z = self._base_ring(entries)
                is_list = False
            except TypeError as er:
                # print "er=",er
                #print "did not work!"
                raise TypeError, "entries must be coercible to a list or complex number. got:{0} of type:{1}!".format(entries,type(entries))
            
            if not z.is_zero():
                if self._nrows != self._ncols:
                    raise TypeError, "nonzero scalar matrix must be square"
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < self._ncols:
                        mpc_set(self._entries[i*self._ncols+j], z.value,self._rnd)
            else:
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < self._ncols:
                        mpc_set(self._entries[i*self._ncols+j], z.value,self._rnd)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        cdef MPComplexNumber y
        y = value
        mpc_set(self._matrix[i][j], y.value,self._rnd)


    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        r"""
        TODO: Make faster if possible
        """
        cdef MPComplexNumber z
        #x = PY_NEW(MPComplexNumber)
        #x = MPComplexNumber(self.parent().base_ring(),0)
        z = self._base_ring(0)
        #print "get_unsafe ",i,j
        mpc_set(z.value, self._matrix[i][j],self._rnd)
        return z

    cpdef int base_for_str_rep(self):
        return self._base_for_str_rep

    cpdef int set_base_for_str_rep(self,int base):
        self._base_for_str_rep=base
        return self._base_for_str_rep

    def zero_matrix(self):
        r"""
        Return the zero matrix in this space.
        """
        assert self._is_square
        cdef Matrix_complex_dense res
        res=Matrix_complex_dense.__new__(Matrix_complex_dense,self._parent,None,None,None)
        cdef int n,i,j,m
        n = self._nrows
        m = self._ncols
        for i from 0 <= i < n:
            for j from 0 <= j < m:
                mpc_set_ui(res._matrix[i][j],0,self._rnd)
        return res
    def identity_matrix(self):
        r"""
        Return the identity matrix in this space.
        """
        assert self._is_square
        cdef Matrix_complex_dense res
        res=Matrix_complex_dense.__new__(Matrix_complex_dense,self._parent,None,None,None)
        cdef int n,i,j,m
        n = self._nrows
        m = self._ncols
        if n != m:
            raise ValueError,"Only square matrices can have identity!"
        for i from 0 <= i < n:
            for j from 0 <= j < n:
                mpc_set_ui(res._matrix[i][j],0,self._rnd)
            mpc_set_ui(res._matrix[i][i],1,self._rnd)
        return res


    cpdef set_zero_elements(self,double tol=0):
        r"""
        Set all entries self with absolute value < self.eps() to zero 
        TODO: come up with a better name for this function... 
        """
        cdef int i,j
        cdef mpfr_t tmp,mptol
        mpfr_init2(tmp,self._prec)
        mpfr_init2(mptol,self._prec)
        if tol>0:
            mpfr_set_d(mptol,tol,self._rnd_re)
        else:
            mpfr_set(mptol,self._eps.value,self._rnd_re)
        for i from 0<=i<self._nrows:
            for j from 0<=j<self._ncols:
                mpc_abs(tmp,self._matrix[i][j],self._rnd_re)
                if mpfr_cmp(tmp,mptol)<=0:
                    mpc_set_si(self._matrix[i][j],0,self._rnd)
        # Since we may have changed pivots We now have to recompute some necessary things...
        self.clear_cache()
        mpfr_clear(tmp)
        mpfr_clear(mptol)
        #piv=self._cache['pivots']=None

    cpdef prec(self):
        return self._prec
#    cpdef eps(self):
#        return self._eps
    
    cpdef int numerical_rank(self,double tol=0):
        r"""
        Find the numerical rank of self up to given tolerance.
        We find it through computing the Q,R decoposition and observe that 
        the rank of R is the same as the ran of self. 

        """
        cdef int i,j,rank,row_is_nonzero
        Q,R=self.qr_decomposition()
        rank = 0
        if tol==0:
            tol=self._eps
        for i in range(0,R.nrows()):
            if abs(R[i,i])>tol:
                rank+=1
        return rank
        #cdef mpfr_t tmp
        #cdef mpfr_t mptol
        #mpfr_init2(tmp,self._prec)
        #mpfr_init2(mptol,self._prec)
        #rank=0
        #if tol>0:
        #    mpfr_set_d(mptol,tol,self._rnd_re)
        #else:
        #    mpfr_set(mptol,self._eps.value,self._rnd_re)
        #for i from 0<=i<self._nrows:
        #    row_is_nonzero=0
        #    for j from i<=j<self._nrows:
        #        mpc_abs(tmp,R._matrix[i][j],self._rnd_re)
        #        if mpfr_cmp(tmp,mptol)>0:
        #            row_is_nonzero=1
        #            break
        #    if row_is_nonzero:
        #        rank+=1
        #return rank
    
    cpdef Vector_complex_dense column(self,int n):
        r""" return column nr. n of self.
        """
        cdef Vector_complex_dense res
        if n>self._ncols or n<0:
            raise IndexError,"Index of column outside range! n={0}, ncols={1}".format(n,self._ncols)
        res = Vector_complex_dense(FreeModule(self._base_ring,self._nrows),0)
        #print "res=",res
        self._column(res._entries,n)
        #print "res1=",res
        return res

    cdef void _column(self,mpc_t  *colv,int n):
        r"""return column of self as pointer of mpc_t """
        cdef int i
        for i from 0 <= i <= self._nrows-1:
            #mpc_init2(colv[i],self._prec)
            mpc_set(colv[i],self._matrix[i][n],self._rnd)

    cpdef Vector_complex_dense row(self,int n):
        r""" return row nr. n of self.
        """
        cdef Vector_complex_dense res
        cdef int i
        if n>self._nrows or n<0:
            raise IndexError,"Index of row outside range! n={0}, nrows={1}".format(n,self._nrows)
        res = Vector_complex_dense(FreeModule(self._base_ring,self._ncols),0)
        for i from 0 <= i <= self._ncols-1:
            mpc_set(res._entries[i],self._matrix[n][i],self._rnd)
        return res

    cpdef delete_row(self,int n,int clear=1):
        r"""
        Delete row nr. n of self.
        Note: I do not deallocate anything unless clear = 1
        """
        cdef int i
        cdef mpc_t *entries_tmp
        cdef int nn,r,k,nnr,nn1
        sig_on()
        if clear == 0:
            for i from n <= i<= self._nrows-2:
                self._matrix[i]=self._matrix[i+1]
                #print i+1,"=>",i
            for i from 0 <= i <= self._ncols-1:
                nn=(self._nrows-1)*self._ncols+i
                mpc_clear(self._entries[nn])
            self._nrows = self._nrows -1
        else:
            self._nrows = self._nrows -1
            entries_tmp = <mpc_t *> sage_malloc(sizeof(mpc_t)*(self._nrows * self._ncols))
            for r from 0 <= r <= self._nrows-1:
                if r < n:
                    nnr=r *self._ncols
                else:
                    nnr=(r+1) *self._ncols
                for k from 0<=k <=self._ncols-1:
                    nn =r *self._ncols + k
                    nn1=nnr + k
                    #print " setting",nn,' to ',nn1
                    mpc_init2(entries_tmp[nn],self._prec)
                    mpc_set(entries_tmp[nn],self._entries[nn1],self._rnd)
            for r from 0 <= r <= self._ncols*(self._nrows+1)-1:
                mpc_clear(self._entries[r])
            sage_free(self._entries)
            self._entries = <mpc_t *> sage_malloc(sizeof(mpc_t)*(self._nrows * self._ncols))
            for r from 0 <= r <= self._nrows*self._ncols-1:
                mpc_init2(self._entries[r],self._prec)
                mpc_set(self._entries[r],entries_tmp[r],self._rnd)
                #print "set entries[",r,"]:",print_mpc(self._entries[r])
            for i from 0 <= i <= self._nrows * self._ncols-1:
                mpc_clear(entries_tmp[i])
            sage_free(entries_tmp)
            ### delete matrix and add the pointers again
            sage_free(self._matrix)
            self._matrix =  <mpc_t **> sage_malloc(sizeof(mpc_t*) * self._nrows)
            if self._matrix == NULL:
                sage_free(self._entries)
                self._entries = NULL
                raise MemoryError, "out of memory allocating a matrix"
            k=0
            for i from 0 <= i < self._nrows:
                self._matrix[i] = self._entries + k
                #for j from 0<= j < self._ncols:
                #    print "set matrix[",i,j,"]:",print_mpc(self._matrix[i][j])
                k = k + self._ncols
        sig_off()
        if self._nrows<>self._ncols:
            self._is_square = 0
        else:
            self._is_square = 1
        new_parent = MatrixSpace(self._parent.base_ring(),self._nrows,self._ncols)
        self._parent = new_parent


    cpdef delete_column(self,int n):
        r"""
        Delete column nr. n of self.
        Note: I do not deallocate anything unless clear = 1
        """
        cdef int i
        cdef mpc_t *entries_tmp
        cdef int nn,r,k

        self._ncols = self._ncols -1
        sig_on()
        entries_tmp = <mpc_t *> sage_malloc(sizeof(mpc_t)*(self._nrows * self._ncols))
            
        for r from 0 <= r <= self._nrows-1:
            for k from 0<=k <=self._ncols-1:
                nn =r *(self._ncols) + k
                if k <= n-1:
                    nn1=nn+r
                else:
                    nn1=nn+r+1
                mpc_init2(entries_tmp[nn],self._prec)
                mpc_set(entries_tmp[nn],self._entries[nn1],self._rnd)
    
        for i from 0 <= i <= self._nrows * self._ncols:
            mpc_clear(self._entries[i])
        sage_free(self._entries)
        self._entries = <mpc_t *> sage_malloc(sizeof(mpc_t)*(self._nrows * self._ncols))
        for r from 0 <= r <= self._nrows*self._ncols-1:
            mpc_init2(self._entries[r],self._prec)
            mpc_set(self._entries[r],entries_tmp[r],self._rnd)                
        for i from 0 <= i <= self._nrows*self._ncols-1:
            mpc_clear(entries_tmp[i])
        sage_free(entries_tmp)
            ### delete matrix and add the pointers again
        sage_free(self._matrix)
        self._matrix =  <mpc_t **> sage_malloc(sizeof(mpc_t*) * self._nrows)
        if self._matrix == NULL:
            sage_free(self._entries)
            self._entries = NULL
            raise MemoryError, "out of memory allocating a matrix"
        sig_off()
        k=0
        for i from 0 <= i < self._nrows:
            self._matrix[i] = self._entries + k
            k = k + self._ncols
        if self._nrows<>self._ncols:
            self._is_square = 0
        else:
            self._is_square = 1
        ## Also change parent of sel
        new_parent = MatrixSpace(self._parent.base_ring(),self._nrows,self._ncols)
        self._parent = new_parent

            
    def eps(self):
        r"""
        """
        #cdef RealNumber x
        #print "eps=",self._eps
        #x=RealNumber(self._base_ring._base,1)
        #mpfr_set(x.value,self._mpeps,self._rnd_re)
        #print "mpeps=",x
        return self._eps
    #cdef _add_ui_unsafe_assuming_int(self, Py_ssize_t i, Py_ssize_t j, unsigned long int n):
    #    # doesn't check immutability
    #    # doesn't do bounds checks.
    #    # assumes that self[i,j] is an integer.
    #    mpc_add_ui(mpc_numref(self._matrix[i][j]), mpc_numref(self._matrix[i][j]), n)
    #
    #cdef _sub_ui_unsafe_assuming_int(self, Py_ssize_t i, Py_ssize_t j, unsigned long int n):
    #    # doesn't check immutability
    #    # doesn't do bounds checks.
    #    # assumes that self[i,j] is an integer.
    #    mpc_sub_ui(mpc_numref(self._matrix[i][j]), mpc_numref(self._matrix[i][j]), n)

    def _pickle(self):
        return self._pickle_version0(), 0
    
    def _unpickle(self, data, int version=0):
        if version == 0:
            self._unpickle_version0(data)
        else:
            raise RuntimeError, "unknown matrix version (=%s)"%version
        
    cdef _pickle_version0(self):
        return self._export_as_string()

    cpdef _export_as_string(self, int base=0, int truncate=0):
        """
        Return space separated string of the entries in this matrix, in the
        given base. This is optimized for speed.
        
        INPUT: base -an integer = 32; (default: 16)
        
        EXAMPLES::
        
            sage: m = matrix(QQ,2,3,[1,2/3,-3/4,1,-2/3,-45/17])
            sage: m._export_as_string(10)
            '1 2/3 -3/4 1 -2/3 -45/17'
            sage: m._export_as_string(16)
            '1 2/3 -3/4 1 -2/3 -2d/11'

            TODO: Use C functions for pickling instead of pyton
            
        """
        cdef Py_ssize_t i, j, len_so_far, m, n
        cdef mpfr_t x,y
        cdef char *a
        cdef char *s
        cdef char *t
        cdef char *tmp
        #print "exporting as string!"
        cdef int reqdigits
        cdef int p,k
        cdef int mpfr_rnd
        cdef mp_exp_t exponent
        #cdef char* data
        # since there are some exceptions for other bases
        # i.e. 7 nd 49 I simply assume either 10 or power of 2 as basis
        cdef list allowed_bases=[2,4,6,8,10,16,32]
            # since for base>62 we do not know the exact length of the strings...
        if base==0:            
            base=self._base_for_str_rep

        if not base in allowed_bases:
            raise ValueError,"Please use either 10 or a power of 2 <=32 for base!"
        if truncate<>self._truncate:
            self._truncate=truncate
        # we want to store the base used for pickling
        mpfr_init2(x,self._prec)
        mpfr_init2(y,self._prec)

        sig_on()
        if self._nrows <> 0 and self._ncols <> 0:
            # how long is the entire string going to be?
            n = 2*self._nrows*self._ncols # = number of real entries
            reqdigits=self._reqdigits(base,truncate)
            #print "required digits is:",reqdigits
            #print "trunc=",self._truncate
            #print "base=",base
            # with exponent stored we also need an extra size of
            emin2 = mpfr_get_emin()
            emax2 = mpfr_get_emax()
            emax2 = max(emax2,emin2)
            emax = ceil(emax2*self._base_ring._base.log2()/self._base_ring(base).log())
            #print "max exponent in base=",emax
            n = n*(reqdigits+1+len(str(emax)))
            s = <char*> sage_malloc(n * sizeof(char))
            #print "alloc len=",n * sizeof(char)
            t = s
        
            exponent = 1
            
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < self._ncols:
                    # the length of a complex number is fixed, depending only on the precision
                    # use builtin functions for making strings...
                    mpc_real(x,self._matrix[i][j],self._rnd_re)
                    mpc_imag(y,self._matrix[i][j],self._rnd_im)
                    mpfr_get_str(t, &exponent, base, reqdigits,x,self._rnd_re)
                    if(exponent>emax):
                        strcpy(t,'@inf@')
                        t=t+strlen('@inf@')
                    elif(exponent<-emax):
                        strcpy(t,'-@inf@')
                        t=t+strlen('-@inf@')
                    else:
                        m=strlen(t)
                        t = t+m
                        # save exponent too
                        #print "exponent.re=",exponent
                        expstr='@'+str(exponent)
                        strcpy(t,expstr)
                        m = strlen(expstr)
                        t=t+m
                    t[0] = <char>32
                    t[1] = <char>0
                    t = t + 1
                    mpfr_get_str(t, &exponent, base, reqdigits,y,self._rnd_im)
                    if(exponent>emax):
                        strcpy(t,'@inf@')
                        t=t+strlen('@inf@')
                    elif(exponent<-emax):
                        strcpy(t,'-@inf@')
                        t=t+strlen('-@inf@')
                    else:
                        m=strlen(t)
                        t = t+m
                        #print "exponent.re=",exponent
                        #raise ValueError,"Can only store exponents smaller than the maximum emax= %s in base %s" % (emax,base)
                        expstr='@'+str(exponent)
                        strcpy(t,expstr)
                        m = strlen(expstr)
                        t=t+m
                    t[0] = <char>32
                    t[1] = <char>0
                    t = t + 1                    
                    #print "t=",str(t)
            #print "data=",data
            #print "len=",strlen(s)
            #print "base=",str(base)+' '
            data=str(base)+' '+str(s)[:-1]
            #print "data=",data
            #data[i:-1] = str(s)[:-1]
            #print "data=",data
        else:
            data=''
            #print "len(data)=",len(data)
            #print "s=",str(s)
        sig_off()
        mpfr_clear(x); mpfr_clear(y)
        sage_free(s)
        return data


    cpdef _reqdigits(self,int base=16,int truncate=0):
        r"""
        returns the required number of digits to represent
        a real number in base with precision self._prec
        """
        cdef int p,reqdigits
        p=self._prec
        if base==2 or base==4 or base==8 or base==16 or base==32:
            p-=1
            k = valuation(base,2)
            reqdigits = 1 + ceil(<double>p/<double>k)
        if base == 10 and truncate<>0:
            # see real_mpfr.pyx for explanation of this...
            reqdigits = 1+ceil( <double>(p-1) * 0.3010299956)
        elif base==10:
            reqdigits = 1+ceil( (<double>p) * 0.3010299956)
        if reqdigits <= 1:
            reqdigits = 2

        return reqdigits
    
    cdef _unpickle_version0(self, data):
        cdef Py_ssize_t i, n,xi,yi,l
        cdef int base
        cdef mpfr_t x,y
        cdef RealNumber tmpr
        tmpr=RealNumber(self._base_ring._base,0)
        base=self._base_for_str_rep
        data = data.split()
        mpfr_init2(x,self._prec);mpfr_init2(y,self._prec)
        n = self._nrows * self._ncols
        reqdigits=self._reqdigits(base,self._truncate)
        #print "trunc=",self._truncate
        #print "base=",base,"=",self._base_for_str_rep
        #print "data=",data
        base = int(data[0])
        if len(data) != 2*n+1:
            raise RuntimeError, "invalid pickle data %s<>%s." %(len(data),self._nrows*self._ncols)
        l=1
        for i from 0 <= i < n:
            #print "data=",data[l]
            # we need to add the decimal point
            if data[l][0]=='-':
                s = '-.'+data[l][1:]
            else:
                s = '.'+data[l]
            #print "s=",str(s)
            xi=mpfr_set_str(x, s, base, self._rnd_re)
            mpfr_set(tmpr.value,x,self._rnd_re)
            #print "x=",tmpr
            if data[l+1][0]=='-':
                s = '-.'+data[l+1][1:]
            else:
                s = '.'+data[l+1]
            yi=mpfr_set_str(y, s, base, self._rnd_im)
            mpfr_set(tmpr.value,y,self._rnd_re)
            #print "y=",tmpr            
            if yi==-1 or xi==-1:
                mpfr_clear(x); mpfr_clear(y)
                raise RuntimeError, "invalid pickle data s=%s" %str(s)
            mpc_set_fr_fr(self._entries[i], x,y,self._rnd)
            l=l+2
        mpfr_clear(x); mpfr_clear(y)
        
    def __richcmp__(Matrix self, right, int op):
        
        
        t1 = type(right)==type(self)
        if t1:
            t2 = self._ncols == (<Matrix_complex_dense>right)._ncols
            t3 = self._nrows == (<Matrix_complex_dense>right)._nrows
        else:
            t2 = 1; t3=1
        if not t1 or not t2 or not t3:
            if op == 2:
                return 0
            if op == 3:
                return 1
            else:
                raise NotImplementedError,"Can not compare Matrix_complex_dense with {0}!".format(type(right))
            #return 0
        else:
            return self._richcmp_(right, op)
        
    # cpdef _eq_(self, right):

    #         return  self._eq_c(<Matrix_complex_dense>right)
    
    
    
    def __hash__(self):
        return self._hash()



   ########################################################################
    # LEVEL 2 functionality
    # x * cdef _add_         
    # x * cdef _mul_
    # x * cdef _vector_times_matrix_
    # x * cdef _cmp_c_impl
    # x * __neg__
    #   * __invert__
    # x * __copy__
    # x * _multiply_classical
    #   * _list -- list of underlying elements (need not be a copy)
    #   * _dict -- sparse dictionary of underlying elements (need not be a copy)
    ########################################################################

    cpdef _lmul_(self, Element right):
        """
        EXAMPLES::
        
            sage: a = matrix(QQ,2,range(6))
            sage: (3/4) * a
            [   0  3/4  3/2]
            [ 9/4    3 15/4]
        """
        cdef Py_ssize_t i
        cdef MPComplexNumber _x
        _x = self.base_ring()(right)
        cdef Matrix_complex_dense M
        M = Matrix_complex_dense.__new__(Matrix_complex_dense, self._parent, None, None, None)
        for i from 0 <= i < self._nrows * self._ncols:
            mpc_mul(M._entries[i], self._entries[i], _x.value,self._rnd)
        return M

    
    cpdef _rmul_(self, Element left):
        """
        EXAMPLES::
        
            sage: a = matrix(QQ,2,range(6))
            sage: (3/4) * a
            [   0  3/4  3/2]
            [ 9/4    3 15/4]
        """
        cdef Py_ssize_t i
        cdef MPComplexNumber _x
        _x = self.base_ring()(left)
        cdef Matrix_complex_dense M
        M = Matrix_complex_dense.__new__(Matrix_complex_dense, self._parent, None, None, None)
        for i from 0 <= i < self._nrows * self._ncols:
            mpc_mul(M._entries[i], self._entries[i], _x.value,self._rnd)
        return M
    
    cpdef ModuleElement _add_(self, _right):
        """
        Add two dense matrices over MPC.
        
        EXAMPLES::
        
            sage: a = MatrixSpace(QQ,3)(range(9))
            sage: b = MatrixSpace(QQ,3)([1/n for n in range(1,10)])
            sage: a+b
            [   1  3/2  7/3]
            [13/4 21/5 31/6]
            [43/7 57/8 73/9]
            sage: b.swap_rows(1,2)
            sage: #a+b
        """
        cdef Py_ssize_t i, j
        cdef Matrix_complex_dense M
        cdef Matrix_complex_dense right = _right
        M = Matrix_complex_dense.__new__(Matrix_complex_dense, self._parent, None, None, None)

        cdef mpc_t *M_row
        cdef mpc_t *self_row
        cdef mpc_t *right_row
        sig_on()
        for i from 0 <= i < self._nrows:
            M_row = M._matrix[i]
            self_row = self._matrix[i]
            right_row = (<Matrix_complex_dense>right)._matrix[i]
            for j from 0 <= j < self._ncols:
                mpc_add(M_row[0], self_row[0], right_row[0],self._rnd)
                M_row = M_row + 1
                self_row = self_row + 1
                right_row = right_row + 1
        sig_off()
        return M
        
    cpdef ModuleElement _sub_(self, _right):
        """
        Subtract two dense matrices over MPC.
        
        EXAMPLES::
        
            sage: a = MatrixSpace(QQ,3)(range(9))
            sage: b = MatrixSpace(QQ,3)([1/n for n in range(1,10)])
            sage: a-b
            [  -1  1/2  5/3]
            [11/4 19/5 29/6]
            [41/7 55/8 71/9]
        """
        cdef Py_ssize_t i, j
        cdef Matrix_complex_dense M
        cdef Matrix_complex_dense right = _right
        cdef MPComplexNumber z
        M = Matrix_complex_dense.__new__(Matrix_complex_dense, self._parent, None, None, None)
        cdef mpc_t *M_row
        cdef mpc_t *self_row
        cdef mpc_t *right_row
        CF=self._base_ring
        if not isinstance(right,Matrix_complex_dense):
            ## then we have to coerce its elements
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < self._ncols:
                    z=CF(right[i][j])
                    mpc_sub(M._matrix[i][j], self._matrix[i][j], z.value,self._rnd)
        else:
            sig_on()
            for i from 0 <= i < self._nrows:
                M_row = M._matrix[i]
                self_row = self._matrix[i]
                right_row = (<Matrix_complex_dense>right)._matrix[i]
                for j from 0 <= j < self._ncols:
                    mpc_sub(M_row[0], self_row[0], right_row[0],self._rnd)
                    M_row = M_row + 1
                    self_row = self_row + 1
                    right_row = right_row + 1
            sig_off()

        return M

    cdef int _cmp_c_impl(self, Element right) except -2:
        cdef mpc_t *a
        cdef mpc_t *b
        cdef Py_ssize_t i, j
        cdef int k
        for i from 0 <= i < self._nrows:
            a = self._matrix[i]
            b = (<Matrix_complex_dense>right)._matrix[i]
            for j from 0 <= j < self._ncols:
                k = mpc_cmp(a[j], b[j])
                if k:
                    if k < 0:
                        return -1
                    else:
                        return 1
        return 0

    # cdef int _eq_c(self, Matrix_complex_dense right) except -2:
    #     cdef mpc_t *a, *b
    #     cdef Py_ssize_t i, j
    #     cdef int k
    #     for i from 0 <= i < self._nrows:
    #         a = self._matrix[i]
    #         b = right._matrix[i]
    #         for j from 0 <= j < self._ncols:
    #             k = mpc_cmp(a[j], b[j])
    #             if k:
    #                 if k < 0:
    #                     return -1
    #                 else:
    #                     return 1
    #     return 0
    cdef Vector _vector_times_matrix_(self, Vector v):
        """
        Returns the vector times matrix product.
        
        INPUT:
        
        
        -  ``v`` - a free module element.
        
        
        OUTPUT: The vector times matrix product v\*A.
        
        EXAMPLES::
        
            sage: B = matrix(QQ,2, [1,2,3,4])
            sage: V = QQ^2
            sage: w = V([-1,5/2])
            sage: w*B
            (13/2, 8)
        """
        cdef Vector_complex_dense w, ans
        cdef Py_ssize_t i, j
        cdef mpc_t x, y
        
        M = self._row_ambient_module()
        w = <Vector_complex_dense> v
        ans = M.zero_vector()

        mpc_init2(x,self._prec)
        mpc_init2(y,self._prec)
        for i from 0 <= i < self._ncols: 
            mpc_set_si(x, 0,self._rnd)
            for j from 0 <= j < self._nrows:
                mpc_mul(y, w._entries[j], self._matrix[j][i],self._rnd)
                mpc_add(x, x, y,self._rnd)
            mpc_set(ans._entries[i], x,self._rnd)
        mpc_clear(x)
        mpc_clear(y)
        return ans


    cdef Vector _matrix_times_vector_(self, Vector v):
        """
        Returns the vector times matrix product.
        
        INPUT:
        
        
        -  ``v`` - a free module element.
        
        
        OUTPUT: The vector times matrix product v\*A.
        
        EXAMPLES::
        

        """
        cdef Vector_complex_dense w, ans
        cdef Py_ssize_t i, j
        cdef mpc_t x, y
        ans = Vector_complex_dense(self._column_ambient_module(),0)
        if isinstance(v,Vector_complex_dense):
            return self._matrix_times_vector_complex_dense(v)
        
        w =  Vector_complex_dense(v.parent(),v)
        mpc_init2(x,self._prec)
        mpc_init2(y,self._prec)
        for i from 0 <= i < self._ncols: 
            mpc_set_si(x, 0,self._rnd)
            for j from 0 <= j < self._nrows:
                mpc_mul(y, w._entries[j], self._matrix[i][j],self._rnd)
                mpc_add(x, x, y,self._rnd)
            mpc_set(ans._entries[i], x,self._rnd)
        mpc_clear(x)
        mpc_clear(y)
        return ans


    cdef Vector _matrix_times_vector_complex_dense(self, Vector_complex_dense v):
        """
        Returns the vector times matrix product.
        
        INPUT:
        
        
        -  ``v`` - a free module element.
        
        
        OUTPUT: The vector times matrix product v\*A.
        
        EXAMPLES::
        

        """
        cdef Vector_complex_dense w, ans
        cdef Py_ssize_t i, j
        cdef mpc_t x, y
        cdef mpc_t* entries
        #print "here v=",v
        ans = Vector_complex_dense(self._column_ambient_module(),0)
        mpc_init2(x,self._prec)
        mpc_init2(y,self._prec)
        for i from 0 <= i < self._ncols: 
            mpc_set_si(x, 0,self._rnd)
            for j from 0 <= j < self._nrows:
                mpc_mul(y, v._entries[j], self._matrix[i][j],self._rnd)
                mpc_add(x, x, y,self._rnd)
            mpc_set(ans._entries[i], x,self._rnd)
        mpc_clear(x)
        mpc_clear(y)
        return ans

    def __mul__(self, ModuleElement right):
         #cdef Vector_complex_dense z
         cdef MPComplexNumber a
         if isinstance(right,Matrix_complex_dense):
             return Matrix_complex_dense._multiply_classical(self,right)
         if is_Vector(right):
             if right.degree() <>self.ncols():
                 raise ValueError,"Dimensions do not match! %s<>%s" %(right.degree(),self.ncols())
             return Matrix_complex_dense._matrix_times_vector_(self,right)
         if is_Matrix(right): ## if any other kind of matrix we have to coerce
             return Matrix_complex_dense._multiply_classical(self,right)
         #return right._vector_times_matrix(self)
         if not isinstance(right,MPComplexNumber):
             z = self._base_ring(right)
         else:
             z = right
         return self._lmul_(z)

    def __neg__(self):
        """
        Negate a matrix over MPComplexFIeld.
        
        EXAMPLES::
        
            sage: a = MatrixSpace(QQ,3)([1/n for n in range(1,10)])
            sage: -a
            [  -1 -1/2 -1/3]
            [-1/4 -1/5 -1/6]
            [-1/7 -1/8 -1/9]
        """
        cdef Py_ssize_t i, j
        cdef Matrix_complex_dense M
        M = Matrix_complex_dense.__new__(Matrix_complex_dense, self._parent, None, None, None)

        cdef mpc_t *M_row
        cdef mpc_t *self_row
        sig_on()
        for i from 0 <= i < self._nrows:
            M_row = M._matrix[i]
            self_row = self._matrix[i]
            for j from 0 <= j < self._ncols:
                mpc_neg(M_row[0], self_row[0],self._rnd)
                M_row = M_row + 1
                self_row = self_row + 1
        sig_off()
        return M
        
    def __copy__(self):
        """
        Copy self .
        
        """
        cdef Py_ssize_t i
        cdef Matrix_complex_dense M
        M = Matrix_complex_dense(self.parent(),0)
        #M = Matrix_complex_dense.__new__(Matrix_complex_dense, self._parent, None, None, None)

        cdef mpc_t *M_row
        cdef mpc_t *self_row
        sig_on()
        for i from 0 <= i < self._nrows*self._ncols:
            mpc_set(M._entries[i],self._entries[i],self._rnd)
        sig_off()
        if self.subdivisions is not None:
            M.subdivide(*self.get_subdivisions())
        return M

    cpdef conjugate_transpose(self):
        res = Matrix_complex_dense(self.parent(),0)
        self._conjugate_transpose(res)
        return res
    
    cdef _conjugate_transpose(self,Matrix_complex_dense res):
        #cdef Matrix_complex_dense res
        cdef int j,k
        cdef mpc_t tmp
        #res = Matrix_complex_dense.__new__(Matrix_complex_dense, self._parent, None, None, None)
        mpc_init2(tmp,self._prec)
        for j from 0<=j<=self._nrows-1:
            for k from 0<=k< j:
                mpc_set(tmp,self._matrix[j][k],self._rnd)
                mpc_conj(tmp,tmp,self._rnd)
                mpc_set(res._matrix[k][j],tmp,self._rnd)
                mpc_set(tmp,self._matrix[k][j],self._rnd)
                mpc_conj(tmp,tmp,self._rnd)
                mpc_set(res._matrix[j][k],tmp,self._rnd)

            mpc_set(tmp,self._matrix[j][j],self._rnd)
            mpc_conj(tmp,tmp,self._rnd)
            mpc_set(res._matrix[j][j],tmp,self._rnd)
                
        mpc_clear(tmp)
        return res

    cpdef transpose(self):
        MS=MatrixSpace(self.base_ring(),self.ncols(),self.nrows())
        res = Matrix_complex_dense(MS,0)
        self._transpose(res)
        return res
    
    cdef _transpose(self,Matrix_complex_dense res):
        #cdef Matrix_complex_dense res
        cdef int j,k
        cdef mpc_t tmp
        mpc_init2(tmp,self._prec)
        #res = Matrix_complex_dense.__new__(Matrix_complex_dense, self._parent, None, None, None)
        #print "transponerar!"
        for j from 0<=j<=self._nrows-1:
            for k from 0<=k<self._ncols-1:
                #mpc_set(tmp,self._matrix[k][j],self._rnd)
                #print "tmp[",k,j,"]=",print_mpc(tmp)
                mpc_set(res._matrix[k][j],self._matrix[j][k],self._rnd)
#                mpc_set(res._matrix[j][k],self._matrix[k][j],self._rnd)
#            mpc_set(res._matrix[j][j],self._matrix[j][j],self._rnd)
        mpc_clear(tmp)
        return res



    def _multiply_classical(self, Matrix_complex_dense right):
        """
        Use the standard `O(n^3)` matrix multiplication algorithm.
        This is never used by default, since it is slow.
        
        EXAMPLES::
        
            sage: n = 3
            sage: a = matrix(QQ,n,range(n^2))/3
            sage: b = matrix(QQ,n,range(1, n^2 + 1))/5
            sage: a._multiply_classical(b)
            [ 6/5  7/5  8/5]
            [18/5 22/5 26/5]
            [   6 37/5 44/5]
        """
        #print "In multiply classical!!"
        if self._ncols != right._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of right."

        cdef Py_ssize_t i, j, k, l, nr, nc, snc
        cdef mpc_t *v
        cdef object parent
        
        nr = self._nrows
        nc = right._ncols
        snc = self._ncols

        if snc == nc:
            # right acts on the space of self
            parent = self.parent()
        elif nr == nc:
            # self acts on the space of right
            parent = right.parent()
        else:
            parent = self.matrix_space(nr, nc)

        cdef Matrix_complex_dense M, _right
        _right = right
        
        M = Matrix_complex_dense.__new__(Matrix_complex_dense, parent, None, None, None)
        #Matrix.__init__(M, parent)    
        # M has memory allocated *and* entries are initialized
        #cdef mpc_t *entries
        #entries = M._entries

        cdef mpc_t s, z
        mpc_init2(s,self._prec)
        mpc_init2(z,self._prec)
        sig_on()
        l = 0
        for i from 0 <= i < nr:
            v = self._matrix[i]
            for j from 0 <= j < nc:
                mpc_set_si(s,0,self._rnd)   # set s = 0
                for k from 0 <= k < snc:
                    mpc_mul(z, v[k], _right._matrix[k][j],self._rnd)
                    mpc_add(s, s, z,self._rnd)
                mpc_set(M._entries[l], s,self._rnd)
                l += 1
        sig_off()
        mpc_clear(s)
        mpc_clear(z)
        return M
    
    # cdef _mul_(self, Matrix right):
    # cdef int _cmp_c_impl(self, Matrix right) except -2:
    # def __invert__(self):
    # def _list(self):
    # def _dict(self):    
    
        
  
   
#### Level 4 functions. Linear algebra routines


    cpdef det(self):
        r"""
        Computing the determinant using QR-decomposition of self.
        """
        if not self.is_square():
            return self.determinant() ### Use generic method
            #raise NotImplementedError,"We only compute determinants of square matrices!"
        if self.is_upper_triangular():
            det = 1
            for n in range(self.nrows()):
                det = det*self[n,n]
            return det

        l=self.eigenvalues()
        det = 1
        for x in l:
            det=det*x
        return det


    cpdef RealNumber error_qr(self):
        r""" Return an estimate for the error in the qr-decomp. of self.
        Reference: Gentleman, 'Error Analysis of QR decomposition by Givens Transformations',
        Linear Algebra and its applications, Vol 10, Issue 3, 1975

        TODO: Check that my implementation actually satisfies the conditions for this bound and make it sharper. 
        """
        if self._error_qr==0:
            RF = self._base_ring._base
            eta = max(self._eps,RF(2)**RF(-self._prec))
            nn = RF(2*self._nrows-2)
            errest = self.norm()*eta*nn**(1+eta)**(nn-RF(1))
            self._error_qr=errest
        return self._error_qr

    #def qr_decomposition(self,overwrite=0,check=0,num_threads=1,schedule=0):
    #    return self._qr_decomposition(overwrite,check,num_threads,schedule)
    
    cpdef tuple qr_decomposition(self,int overwrite=0, int check=0,int num_threads=1,int schedule=0):
        r"""
        Computes the QR-dcomposition of self using Givens rotations.
        INPUT:

         - overwrite -- set to 1 if you don't mind self._matrix being overwritten.
                        othrwise 0.
        NOTE: num_threads>1 is not implemented.
        """
        cdef int n,m,i,j,ii
        cdef mpc_t** Q
        cdef mpc_t** R
        cdef mpc_t** A
        cdef Matrix_complex_dense QQ,RR
        cdef mpc_t s,sc,t
        cdef mpfr_t x,y,c
        cdef mpc_t tt[3]
        if num_threads>1:
            raise NotImplementedError,"Parallel linear algebra is not implemented!"
        mpc_init2(s,self._prec); mpc_init2(t,self._prec)
        mpc_init2(sc,self._prec)        
        mpfr_init2(x,self._prec);mpfr_init2(y,self._prec)
        mpfr_init2(c,self._prec)
        m = self._nrows
        n = self._ncols
        if overwrite:
            qr_decomp(self._matrix,n,n,self._prec,self._mpeps,self._rnd,self._rnd_re)
        else:
            A = <mpc_t **> sage_malloc(sizeof(mpc_t*) * n)
            if A==NULL:
                raise MemoryError            
            for i from 0 <= i <  n:
                A[i]=<mpc_t *> sage_malloc(sizeof(mpc_t) * n)
                for j from 0 <= j < n:
                    mpc_init2(A[i][j],self._prec+100)
                    mpc_set(A[i][j],self._matrix[i][j],self._rnd)
            if num_threads==1:
                qr_decomp(A,n,n,self._prec+100,self._mpeps,self._rnd,self._rnd_re)
            else:
                raise NotImplementedError,"Parallel linear algebra is not implemented!"
                #qr_decomp_par(A,n,n,self._prec+100,self._mpeps,num_threads,schedule,self._rnd,self._rnd_re)
            #if i<>1:
            #    raise ArithmeticError,"Could not compute QR-decomposition of self!"
            # reconstruct Q,R from self or A
            # RR=Matrix_complex_dense.__new__(Matrix_complex_dense,self._parent,None,None,None)
            RR=Matrix_complex_dense(self._parent,None,None)
            QQ=Matrix_complex_dense(self._parent,None,None)
        #QQ=Matrix_complex_dense.__new__(Matrix_complex_dense,self._parent,None,None,None)
        for i from 0 <= i < n:
            for j from i <= j < n:
                if overwrite:
                    mpc_set(RR._matrix[i][j],self._matrix[i][j],self._rnd)
                else:
                    mpc_set(RR._matrix[i][j],A[i][j],self._rnd)
            for j from 0 <= j < i:
                mpc_set_ui(RR._matrix[i][j],0,self._rnd)
        if overwrite:
            _reconstruct_matrix(QQ._matrix,self._matrix, m, n, 1, self._prec, self._rnd, self._rnd_re)
        else:
            _reconstruct_matrix(QQ._matrix,A, m, n, 1, self._prec, self._rnd, self._rnd_re)
        if not overwrite:
            for i from 0 <= i < n:
                for j from 0 <= j < n:
                    mpc_clear(A[i][j])
            sage_free(A)
        ## Set the errors in RR and QQ as estimated by error_qr
        cdef RealNumber errest=self.error_qr()
        cdef int prec = ceil(errest.log2())
        if prec<0 and -prec>self._prec:
            QQ.set_prec(-prec)
            RR.set_prec(-prec)
        mpc_clear(s); mpc_clear(t); mpc_clear(sc)   
        mpfr_clear(x);mpfr_clear(y); mpfr_clear(c)
        return QQ,RR
    

    cpdef hessenberg(self,int return_transformation=0,int check=0,int num_threads=1,int schedule=0):
        r"""
        Reduce self to Hessenberg form and optionally store
        the transformations.
        This method should primarily be used when we need the
        transformations.
        NOTE: num_threads>1 is not implemented
        """
        cdef int i,j,n,m
        cdef QR_set q = init_QR(self._prec)
        cdef Matrix_complex_dense Q
        cdef mpc_t t,s,sc
        cdef mpc_t tt[3]
        cdef mpfr_t x,c
        m = self._nrows
        n = self._ncols
        if num_threads>1:
            raise NotImplementedError,"Parallel linear algebra is not implemented!"
        if num_threads==1:
            _hessenberg_reduction(self._matrix, n, q, self._prec, self._rnd, self._rnd_re, return_transformation)
        else:
            raise NotImplementedError,"Parallel linear algebra is not implemented!"
            #_hessenberg_reduction_par(self._matrix, n, q, self._prec, num_threads,schedule,self._rnd, self._rnd_re, return_transformation)       
        if return_transformation:
            Q = Matrix_complex_dense.__new__(Matrix_complex_dense,self._parent,None,None,None)
            _reconstruct_matrix(Q._matrix,self._matrix, m, n, 2, self._prec, self._rnd, self._rnd_re)
            ## mpc_init2(tt[0],self._prec); mpc_init2(tt[1],self._prec)
            return Q


    cpdef _balance(self):
        r""" Balance self. """
        cdef int n,i,ii,j,last
        cdef mpfr_t x,c,r,g,f,t,s,q
        cdef mpc_t z
        cdef RealNumber tmpr
        #tmpr=self._base_ring._base(1)
        mpc_init2(z,self._prec)
        mpfr_init2(c,self._prec)
        mpfr_init2(r,self._prec)
        mpfr_init2(f,self._prec)
        mpfr_init2(g,self._prec)
        mpfr_init2(t,self._prec)
        mpfr_init2(s,self._prec)
        mpfr_init2(x,self._prec)
        mpfr_init2(q,self._prec)        
        mpfr_set_d(q,0.95,self._rnd_re)
        assert self._is_square
        n = self._nrows
         # we shouldn't need more tries than this
        for ii in range(n*n):
            last=1
            for i in range(n):
                mpfr_set_ui(c,0,self._rnd_re)
                mpfr_set_ui(r,0,self._rnd_re)
                for j in range(n):
                    mpc_abs(x,self._matrix[j][i],self._rnd_re)
                    mpfr_add(c,c,x,self._rnd_re)
                    mpc_abs(x,self._matrix[i][j],self._rnd_re)
                    mpfr_add(r,r,x,self._rnd_re)
                    #c=sum(abs(A[j,i]) for j in xrange(n))
                    #r=sum(abs(A[i,j]) for j in xrange(n))
                    # if c<>0 and r<>0
                if  mpfr_zero_p(c)==0 and mpfr_zero_p(r)==0:
                    mpfr_div_ui(g,r,2,self._rnd_re)
                    mpfr_set_ui(f,1,self._rnd_re) #f=ctx.one()
                    mpfr_add(s,c,r,self._rnd_re)
                    while mpfr_cmp(c,g)<0 :
                        mpfr_mul_ui(f,f,2,self._rnd_re)
                        mpfr_mul_ui(c,c,2,self._rnd_re)
                    mpfr_mul_ui(g,r,2,self._rnd_re)
                    while mpfr_cmp(c,g)>0 :
                        mpfr_div_ui(f,f,2,self._rnd_re)
                        mpfr_div_ui(c,c,2,self._rnd_re)
                    mpfr_add(t,c,r,self._rnd_re)
                    mpfr_div(t,t,f,self._rnd_re)
                    mpfr_div(t,t,q,self._rnd_re)
                    if mpfr_cmp(t,s)<0:
                        # if((c+r)/f<0.95*s):
                        last=0
                        mpfr_ui_div(g,1,f,self._rnd_re) #=ctx.one()/f
                        for j in range(n):
                            #print "i,j=",i,j
                            #mpfr_set(tmpr.value,g,self._rnd_re)
                            #print "g=",tmpr
                            mpc_mul_fr(self._matrix[i][j],self._matrix[i][j],g,self._rnd)
                            #mpc_set(self._matrix[i][j],z,self._rnd)
                            #mpfr_set(tmpr.value,f,self._rnd_re)
                            #print "f=",tmpr
                            mpc_mul_fr(self._matrix[i][j],self._matrix[j][i],f,self._rnd)
                            #mpc_set(self._matrix[j][i],z,self._rnd)
                        #A[j,i]=A[j,i]*f
            if(last==0):
                continue
            else:
                break
        mpc_clear(z)
        mpfr_clear(c);mpfr_clear(r);mpfr_clear(f);mpfr_clear(g);mpfr_clear(t)
        mpfr_clear(s);mpfr_clear(x);mpfr_clear(q)
        

        
    ## ### test various properties of self
    cpdef int is_square(self):
        r"""

        """
        return self._is_square



    cpdef int is_hessenberg(self,double maxerr=0,int show_err=0):
        r"""
        Is self of Hessenberg form up to prec?
        """
        cdef int t
        cdef mpfr_t y
        t=_is_hessenberg(self._entries,self._nrows,self._prec,self._mpeps,self._rnd,self._rnd_re,y,show_err)
        if t==0 and show_err:
            maxerr=mpfr_get_d(y,self._rnd_re)
            return 0
        else:
            return t
   

    cpdef int is_upper_triangular(self,int return_err=0):
        r"""
        Returns 1 if self is upper-triangular, else 0
        """
        cdef mpfr_t maxerr
        cdef RealNumber err
        cdef int t
        mpfr_init2(maxerr,self._prec)
        assert self._is_square
        t = _is_upper_triangular(self._matrix,self._nrows,self._prec,self._mpeps, self._rnd, self._rnd_re,maxerr,return_err)
        if return_err:
            err=self._base_ring._base(0)
            mpfr_set(err.value,maxerr,self._rnd_re)
            mpfr_clear(maxerr)
            return (t,err)
        mpfr_clear(maxerr)
        return t
        
    
    



    cpdef int _pivot_element(self,int k,int r):
        r"""
        Returns 1 if self is upper-triangular, else 0
        """
        cdef int kpiv,n,m
        n = <int>self._nrows
        m = <int>self._ncols
        #kpiv=-1
        kpiv= _pivot_element(k,r,self._matrix,&n,&m,&self._prec,self._rnd_re)
        return kpiv





    cpdef int is_orthogonal(self):
        r"""
        """
        cdef mpc_t z,w,s
        cdef mpfr_t x,one
        cdef int i,j,kk
        #assert self._is_square
        if self._is_square==0:
            return 0
        mpc_init2(s,self._prec);   mpc_init2(z,self._prec)
        mpc_init2(w,self._prec);   mpfr_init2(x,self._prec)
        mpfr_init2(one,self._prec)
        mpfr_set_ui(one,1,self._rnd_re)
        #n = (self*self.transpose().conjugate()-self._parent.one()).norm()
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                mpc_set_ui(s,0,self._rnd)
                for kk from 0 <= kk < self._ncols:
                   mpc_set(z,self._matrix[j][kk],self._rnd)
                   mpc_mul(z,z,self._matrix[i][kk],self._rnd)
                   mpc_add(s,s,z,self._rnd)
                if i<>j:
                    mpc_abs(x,s,self._rnd_re)
                    if mpfr_cmp(x,self._mpeps)>0:
                        mpc_clear(s);mpc_clear(z);mpc_clear(w);mpfr_clear(x);mpfr_clear(one)
                        return 0
                if i==j:
                    mpc_sub_fr(s,s,one,self._rnd)
                    mpc_abs(x,s,self._rnd_re)
                    if mpfr_cmp(x,self._mpeps)>0:
                        mpc_clear(s);mpc_clear(z);mpc_clear(w);mpfr_clear(x);mpfr_clear(one)
                        return 0
        mpc_clear(s);mpc_clear(z);mpc_clear(w);mpfr_clear(x);mpfr_clear(one)
        return 1

    cpdef int is_unitary(self,int return_err=0):
        r"""
        """
        cdef mpfr_t maxerr
        cdef RealNumber err
        cdef int t
        assert self._is_square
        mpfr_init2(maxerr,self._prec)
        t = _is_unitary(self._matrix,self._nrows,self._prec,self._mpeps, self._rnd, self._rnd_re,maxerr,return_err)        
        if return_err:
            err = self._base_ring._base(0)
            mpfr_set(err.value,maxerr,self._rnd_re)
            mpfr_clear(maxerr)
            return (t,err)
        mpfr_clear(maxerr)
        return t
        
    
    cpdef int is_hermitian(self):
        r"""
        Is self[i,j]=self[j,i].conjugate() or not?
        """
        assert self._is_square
        return _is_hermitian(self._entries,self._nrows,self._prec,self._mpeps, self._rnd, self._rnd_re)


    cpdef RealNumber norm(self,int ntype=2):
        r"""
        Compute the norm of self.
        ntype = 2 => Frobenius (or 2-) norm
        ntype = 0 => max norm
        """
        cdef RealNumber res
        cdef mpfr_t s
        if not self._norms:
            self._norms=dict()
        if self._norms.has_key(ntype):
            return self._norms[ntype]
        res = RealNumber(self._base_ring._base,0)
        mpfr_init2(s,self._prec)
        self._norm(s,ntype)
        mpfr_set(res.value,s,self._rnd_re)
        mpfr_clear(s)
        self._norms[ntype]=res
        return res

    cdef void _norm(self,mpfr_t norm,int ntype):
        r"""
        """
        _norm(norm,self._matrix,self._nrows,self._ncols, self._prec,self._rnd_re,ntype)
        
        #return s       

   
        


    cpdef Vector_complex_dense solve(self,Vector_complex_dense b,int overwrite=0,int num_threads=1,int schedule=0):
        r"""
        Self should be n x n. 
        Solve self*X=b using QR-decomposition.
        If overwrite = 1 then we overwrite A in the process off solving.
        Note: num_threads>1 is not implemented.
        """
        #cdef mpc_t** Q
        #cdef mpc_t** R
        #cdef mpc_t** A
        cdef Vector_complex_dense res
        assert self._is_square
        cdef int n = self._nrows
        assert len(b)==n
        if num_threads>1: raise NotImplementedError,"Parallel linear algebra is not implemented!"
        cdef mpc_t s,sc
        cdef mpc_t t[3]
        cdef mpc_t w[2]
        cdef mpfr_t c,x
        cdef mpc_t **A
        cdef mpc_t *v
        
        res = Vector_complex_dense.__new__(Vector_complex_dense,b._parent,None,False,False)
        #cdef int k
        #print "prec=",self._prec
        #print "b=",b
        mpc_init2(s,self._prec); mpc_init2(sc,self._prec)
        mpc_init2(t[0],self._prec); mpc_init2(t[1],self._prec); mpc_init2(t[2],self._prec)
        mpfr_init2(c,self._prec); mpfr_init2(x,self._prec) 
        mpc_init2(w[0],self._prec); mpc_init2(w[1],self._prec)
        #if not overwrite:
        #
        v = <mpc_t *> sage_malloc(sizeof(mpc_t) * n)
        for i from 0 <= i <  n:
            
            mpc_init2(v[i],self._prec)
            mpc_set(v[i],b._entries[i],self._rnd)
        if overwrite == 1:
            if num_threads==1:
                qr_decomp(self._matrix,n,n,self._prec,self._mpeps,self._rnd,self._rnd_re)
            else:
                raise NotImplementedError,"Parallel linear algebra isnot yet implemented!"
                #with nogil:
                #    qr_decomp_par(self._matrix,n,n,self._prec,self._mpeps,num_threads,schedule,self._rnd,self._rnd_re)
        else:
            A = <mpc_t **> sage_malloc(sizeof(mpc_t*) * n)
            for i from 0 <= i <  n:
                A[i]=<mpc_t *> sage_malloc(sizeof(mpc_t) * n)
                for j from 0 <= j < n:
                    mpc_init2(A[i][j],self._prec)
                    mpc_set(A[i][j],self._matrix[i][j],self._rnd)
            if num_threads==1:
                qr_decomp(A,n,n,self._prec,self._mpeps,self._rnd,self._rnd_re)
            else:
                raise NotImplementedError,"Parallel linear algebra isnot yet implemented!"
                #with nogil:
                #    qr_decomp_par(A,n,n,self._prec,self._mpeps,num_threads,schedule,self._rnd,self._rnd_re)
        # compute b' = Q.conjugate().transpose()*b 
        #
        #        for j in xrange(n-2,-1,-1):
        #            for i in xrange(n-1,j,-1):
        for j from 0 <= j < n-1:
            for i from j+1 <= i <n:
                if overwrite==1:
                    mpc_set(t[0],self._matrix[i][j],self._rnd)
                    mpc_set_ui(self._matrix[i][j],0,self._rnd)
                else:
                    mpc_set(t[0],A[i][j],self._rnd)
                    mpc_set_ui(A[i][j],0,self._rnd)
                if mpfr_inf_p(t[0].re) or mpfr_inf_p(t[0].im):
                    mpfr_set_ui(c,0,self._rnd_re)
                    mpc_set_ui(s,1,self._rnd)
                    mpc_set_ui(sc,1,self._rnd)
                else:
                    mpc_abs(x,t[0],self._rnd_re)
                    mpfr_mul(x,x,x,self._rnd_re)
                    mpfr_add_ui(x,x,1,self._rnd_re)
                    mpfr_sqrt(x,x,self._rnd_re)
                    mpfr_ui_div(c,1,x,self._rnd_re)
                    mpc_mul_fr(s,t[0],c,self._rnd)
                    mpc_conj(sc,s,self._rnd)
                # multiuply from left with the same rotation
                # as in the reduction
                ## k = j 
       
                _mpc_mul_fr(&t[0], v[j], c,self._rnd, self._rnd_re)
                _mpc_mul(&t[1], s, v[i],w,self._rnd,self._rnd_re)
                _mpc_add(&t[2], t[1], t[0],self._rnd_re)
                _mpc_mul(&t[0], sc, v[j],w,self._rnd,self._rnd_re)
                _mpc_mul_fr(&t[1], v[i], c,self._rnd, self._rnd_re)
                _mpc_sub(&v[i], t[1], t[0],self._rnd_re)
                _mpc_set(&v[j],t[2],self._rnd_re)


        #print "Q'b="
        #for j from 0 <= j < n:
        #    print "V[",j,"]=",print_mpc(v[j])
        #print "R=",self=
        # and then solve Rx=b'
        if overwrite:
            solve_upper_triangular(self._matrix, v,n, self._prec, self._rnd, self._rnd_re)
        else:
            solve_upper_triangular(A, v,n, self._prec, self._rnd, self._rnd_re)
        #print "b=",b
        #res = b
        #res = vector(self._base_ring,n)
        #for i from 0 <= i < n:
        #    mpc_set(res[i],b._entries[i],self._rnd)
        #
        for i from 0 <= i < n:
            mpc_init2(res._entries[i],self._prec)
            mpc_set(res._entries[i],v[i],self._rnd)
        for i from 0 <= i < n:
            mpc_clear(v[i])
        sage_free(v)
        mpc_clear(s); mpc_clear(sc)
        mpc_clear(t[0]); mpc_clear(t[1]); mpc_clear(t[2])
        mpfr_clear(c); mpfr_clear(x)
        mpc_clear(w[0]);mpc_clear(w[1])
        if overwrite <> 1:
            for i from 0 <= i < n:
                for j from 0 <= j < n:
                    mpc_clear(A[i][j])
            sage_free(A)
        return res
  


    cpdef Matrix_complex_dense inverse(self,int overwrite=0):
        r"""
        Inverse of self.
        """
        assert self._is_square
        cdef Matrix_complex_dense res,B
        #B = self.parent().one()
        B=Matrix_complex_dense.__new__(Matrix_complex_dense,self._parent,None,None,None)
        cdef int n,i,j
        n = self._nrows
        ## print "B0=",B
        for i from 0 <= i < n:
            for j from 0 <= j < n:
                mpc_set_ui(B._matrix[i][j],0,self._rnd)
            for j from i+1 <= j < n:
                mpc_set_ui(B._matrix[i][j],0,self._rnd)
            mpc_set_ui(B._matrix[i][i],1,self._rnd)
        ## print "B1=",B
        res = self.mat_solve(B,overwrite)
        #print "B0=",B
        return res

    cpdef Matrix_complex_dense mat_solve(self,Matrix_complex_dense B,int overwrite=0):
        r"""
        Self should be n x n. 
        Solve self*X=B using QR-decomposition.
        If overwrite = 1 then we overwrite A in the process of solving.
        """
        assert self._is_square
        cdef Matrix_complex_dense res
        res = Matrix_complex_dense.__new__(Matrix_complex_dense,self._parent,None,None,None)
        cdef mpc_t **A
        cdef mpc_t **BB
        cdef int n = self._nrows
        cdef mpc_t s,sc
        cdef mpc_t t[3]
        cdef mpc_t w[2]
        cdef mpfr_t c,x
        assert B._nrows == n
        #print "B=",B
        cdef int k,j,i
        #print "prec=",self._prec
        #print "b=",b
        mpc_init2(s,self._prec); mpc_init2(sc,self._prec)
        mpc_init2(t[0],self._prec); mpc_init2(t[1],self._prec); mpc_init2(t[2],self._prec)
        mpfr_init2(c,self._prec); mpfr_init2(x,self._prec) 
        mpc_init2(w[0],self._prec); mpc_init2(w[1],self._prec)
        if overwrite==1:
            qr_decomp(self._matrix,n,n,self._prec,self._mpeps,self._rnd,self._rnd_re)
        else:
            A = <mpc_t **> sage_malloc(sizeof(mpc_t*) * n)
            for i from 0 <= i <  n:
                A[i]=<mpc_t *> sage_malloc(sizeof(mpc_t) * n)
                for j from 0 <= j < n:
                    mpc_init2(A[i][j],self._prec)
                    mpc_set(A[i][j],self._matrix[i][j],self._rnd)
            qr_decomp(A,n,n,self._prec,self._mpeps,self._rnd,self._rnd_re)
        for j from 0 <= j < B._nrows:
            for k from 0 <= k < B._ncols:
                mpc_set(res._matrix[j][k],B._matrix[j][k],self._rnd)
        # compute b' = Q.conjugate().transpose()*b 
        v =<mpc_t *> sage_malloc(sizeof(mpc_t) * n)
        for j from 0 <= j < n:
            mpc_init2(v[j],self._prec)

        for j from 0 <= j < n-1:
            for i from j+1 <= i <n:
                if overwrite==1:
                    mpc_set(t[0],self._matrix[i][j],self._rnd)
                    mpc_set_ui(self._matrix[i][j],0,self._rnd)
                else:
                    mpc_set(t[0],A[i][j],self._rnd)
                    mpc_set_ui(A[i][j],0,self._rnd)
                if mpfr_inf_p(t[0].re) or mpfr_inf_p(t[0].im):
                    mpfr_set_ui(c,0,self._rnd_re)
                    mpc_set_ui(s,1,self._rnd)
                    mpc_set_ui(sc,1,self._rnd)
                else:
                    mpc_abs(x,t[0],self._rnd_re)
                    mpfr_mul(x,x,x,self._rnd_re)
                    mpfr_add_ui(x,x,1,self._rnd_re)
                    mpfr_sqrt(x,x,self._rnd_re)
                    mpfr_ui_div(c,1,x,self._rnd_re)
                    mpc_mul_fr(s,t[0],c,self._rnd)
                    mpc_conj(sc,s,self._rnd)
                        # multiuply from left with the same rotation
                        # as in the reduction
                        ## k = j 
                for k from 0 <= k < B._ncols:
                    #for l from 0 <= l < n:
                    #    mpc_set(v[l],B._matrix[l][k],self._rnd)
                    #    print "v[",l,"]=",print_mpc(v[l])
                    #print "c=",print_mpfr(c)
                    #print "s=",print_mpc(s)
                    #print "sc=",print_mpc(sc)
                    mpc_set(v[i],res._matrix[i][k],self._rnd)
                    mpc_set(v[j],res._matrix[j][k],self._rnd)
                    #print "v[",i,"][",k,"]0=",print_mpc(v[i])
                    #print "v[",j,"][",k,"]0=",print_mpc(v[j])
                    _mpc_mul_fr(&t[0], v[j], c,self._rnd, self._rnd_re)
                    _mpc_mul(&t[1], s, v[i],w,self._rnd,self._rnd_re)
                    _mpc_add(&t[2], t[1], t[0],self._rnd_re)
                    _mpc_mul(&t[0], sc, v[j],w,self._rnd,self._rnd_re)
                    _mpc_mul_fr(&t[1], v[i], c,self._rnd, self._rnd_re)
                    _mpc_sub(&v[i], t[1], t[0],self._rnd_re)
                    _mpc_set(&v[j],t[2],self._rnd_re)
                    #print "v[",i,"][",k,"]1=",print_mpc(v[i])
                    #print "v[",j,"][",k,"]1=",print_mpc(v[j])
                    mpc_set(res._matrix[i][k],v[i],self._rnd)
                    mpc_set(res._matrix[j][k],v[j],self._rnd)
        for k from 0 <= k < B._ncols:
            #print "Q'b [k:",k,"]"
            for j from 0 <= j < n:
                mpc_set(v[j],res._matrix[j][k],self._rnd)
            #print "V[",j,",",k,"]=",print_mpc(v[j])
            if overwrite:
                solve_upper_triangular(self._matrix, v,n, self._prec, self._rnd, self._rnd_re)
            else:
                solve_upper_triangular(A, v,n, self._prec, self._rnd, self._rnd_re)
            for j from 0 <= j < n:
                mpc_set(res._matrix[j][k],v[j],self._rnd)
                #print "res[",j,k,"]=",print_mpc(res._matrix[j][k])
            # print "b=",b
            # res = b
            # res = vector(self._base_ring,n)
            # for i from 0 <= i < n:
            #    mpc_set(res[i],b._entries[i],self._rnd)
            #
        for j from 0 <= j < n:
            mpc_clear(v[j])
        sage_free(v)
        mpc_clear(s); mpc_clear(sc); mpc_clear(t[0]); mpc_clear(t[1]); mpc_clear(t[2])
        mpfr_clear(c); mpfr_clear(x)
        mpc_clear(w[0]);mpc_clear(w[1])
        if overwrite <> 1:
            for i from 0 <= i < n:
                for j from 0 <= j < n:
                    mpc_clear(A[i][j])
            sage_free(A)
        return res

    cpdef list eigenvalues(self,int check=0,int sorted=0,int overwrite=0,int num_threads=1,int schedule=0):
        r"""
        Compute the eigenvalues of self.
        sorted = 0  => return unsorted list
        sorted = 1  => return sorted list according to abs-value
        sorted = 2  => return sorted list according to first real, then imaginary value
        num_threads = 1 => use one thread (i.e. no parallelization)
                 = n>1 -- use n threads. 
                 =-1   -- use automatic number of threads
        schedule = 0 => dynamic scheduling
                 = 1 => static scheduling
        NOTE: num_threads >1 is not implemented
        """
        from sage.all import deepcopy
        if num_threads>1: raise NotImplementedError,"Parallel linear algebra is not implemented!"
        cdef mpc_t* evs
        cdef list res
        cdef int i
        cdef MPComplexNumber z
        cdef mpfr_t tmp
        cdef mpc_t** A
        assert self._is_square
        assert schedule == 0 or schedule == 1
        mpfr_init2(tmp,self._prec)
        #from linalg_complex_dense cimport Eigenvalues_of_M
        z = self._base_ring(1)
        evs = <mpc_t *> sage_malloc(sizeof(mpc_t) * self._nrows)
        for i in xrange(self._nrows):
            mpc_init2(evs[i],self._prec)
        if self.is_immutable() or overwrite==0:
            A = <mpc_t **> sage_malloc(sizeof(mpc_t*) * self._nrows)
            for i from 0 <= i <  self._nrows:
                A[i]=<mpc_t *> sage_malloc(sizeof(mpc_t) * self._nrows)
                for j from 0 <= j < self._nrows:
                    mpc_init2(A[i][j],self._prec)
                    mpc_set(A[i][j],self._matrix[i][j],self._rnd)
            if num_threads==1:
                _eigenvalues(evs,A,self._nrows,self._prec,self._rnd,self._rnd_re,self._verbose)
            else:
                raise NotImplementedError,"Parallel linear algebra is not implemented!"
                #_eigenvalues_par(evs,A,self._nrows,self._prec,num_threads,schedule,self._verbose,self._rnd,self._rnd_re)
            #for i in range(self._nrows):
            #    print "evs[",i,"]=",print_mpc(evs[i])
            sage_free(A)
        else:
            if num_threads==1:
                _eigenvalues(evs,self._matrix,self._nrows,self._prec,self._rnd,self._rnd_re,self._verbose)
            else:
                raise NotImplementedError,"Parallel linear algebra is not implemented!"
                #_eigenvalues_par(evs,self._matrix,self._nrows,self._prec,num_threads,schedule,self._verbose,self._rnd,self._rnd_re)
                #_eigenvalues_par(evs,A,self._nrows,self._prec,self._rnd,self._rnd_re,nthreads)
        z=self._base_ring(1)
        res  = list()
        for i from 0 <= i < self._nrows:
            z = self._base_ring(0)
            mpc_abs(tmp,evs[i],self._rnd_re)
            if mpfr_cmp(tmp,self._eps.value)>0:
                mpc_set(z.value,evs[i],self._rnd)
            if self._verbose>0:
                mpc_set(z.value,evs[i],self._rnd)
                print "z[",i,"]=",z,type(z)
            #res.append(deepcopy(z))
            res.append(z)
            
        #for i from 0 <= i < self._nrows:
        #    print "RES1[",i,"]=",res[i],type(res[i])
        if sorted==0:
            res.sort(cmp=my_abscmp)
        elif sorted==1:
            res.sort(cmp=my_lexcmp_re)
        elif sorted==-1:
            res.sort(cmp=my_lexcmp_im)
        if evs<>NULL:
            for i from 0<=i<self._nrows:
                if evs[i]<>NULL:
                    mpc_clear(evs[i])
            sage_free(evs)
        mpfr_clear(tmp)
        return res

    cpdef list eigenvectors(self,int check=1,int overwrite=0,int sorted=0,int verbose=0,int depth=0,double old_tol=0.0,double old_eps=0.0,int num_threads=1,int schedule=0):
        r"""
        Ineffective (I guess...) method of computing eigenvectors one by one

        INPUT::
        
        
         - 'sorted' -- integer (default 0). Set to 0 = > eigenvalues are sorted according to abs-value
                       sorted = 1 = > eigenvalues are sorted according to first real, then imaginary value
                       sorted = -1 = > eigenvalues are sorted according to first imaginary, then real value

         - 'check' -- integer (default 0). Set to 1 => Increase precision until machine epsilon precision is reached.
        
        """
        assert self._is_square
        from sage.all import deepcopy
        #cdef mpc_t* evs
        cdef list res=[],vecs=[]
        cdef int i,j,evi,jj
        cdef mpfr_t tmp
        cdef MPComplexNumber z,z2,ev
        cdef Vector_complex_dense rhs,eigenvec
        cdef RealNumber tol,new_tol,eps,new_eps
        cdef int prec00,prec0 = self._prec
        cdef double tolfak
        cdef dict swaps
        cdef tuple pivots
        swaps={}
        if num_threads>1: raise NotImplementedError,"Parallel linear algebra is not implemented!"
        ## This is the error we estimate from the eigenvalues
        tolfak=10.0
        if depth>1 and check==0:
            ## We don't care about loosing some precision in each step...
            ## TODO: figure out exactly when, where and why we lose precision...
            tolfak=tolfak*100.0
        # Get eigenvalues first
        new_tol=self._base_ring._base(tolfak)*self.error_qr()
        if old_tol>0 and old_tol<1.0:
            tol = self._base_ring._base(old_tol)
        else:
            tol = new_tol
        ## Note that if A is known to self.eps() then we can never expect more than this precision
        ## in the eigenvalue or eigenvector 
        new_eps = self.nrows()*self.eps()
        if old_eps>0 and old_eps<1.0:
            eps = self._base_ring._base(old_eps)
        else:
            eps = new_eps    
        evs=self.eigenvalues(sorted=sorted,num_threads=num_threads,schedule=schedule)
        #_eigenvalues(evs,self._matrix,self._nrows,self._prec,self._rnd,self._rnd_re)
        if verbose>0:
            print "tol=",tol
            print "eps=",eps
            print "ev0=",evs
        z = self._base_ring(0)
        z2 = self._base_ring(0)
        ## We set numerically zero eigenvalues to exactly zero        
        cdef int numz=0
        for evi from 0<=evi<self._nrows:
            if abs(evs[evi])<tol:
                evs[evi]=self._base_ring(0)
                numz+=1
        if verbose>0:
            print "ev1=",evs
        ## The point with computing this is that we know that
        ## This number of diagonal elements should be zero in R below
        cdef int dim,do_cont=0,nume=0
        for evi from 0<=evi<self._nrows:
            z=self._base_ring(evs[evi])
            if verbose>0:
                print "z=",z
            if verbose>1:
                print "have res:",res
            ## If we already have this eigenvalue (with multiplicity) in the list we skip it now.
            do_cont=0
            nume=0
            if len(res)>0:
                for z2,l in res:
                    nume+=len(l)
                    if z2==z and len(l)==evs.count(z):
                        if verbose>0:
                            print "Already have this eigenvalue with corr. mult.!"
                        do_cont=1
                        break
            if nume>=self.ncols():
                do_cont=1
            if do_cont==1:
                continue
            for i from 0<=i<self._nrows:
                mpc_sub(self._matrix[i][i],self._matrix[i][i],z.value,rnd)

                # then solve (A-lambda*I)=0 to find eigenvectors
            #print "self-lamba*E=",self
            Q,R=self.qr_decomposition()
            ## The eigenvectors span the kernel of R
            if check>=2:
                test = (Q*R-self).norm()
                if verbose>1:
                    print "||QR-A||=",test                
            if verbose>1:
                print "A-lambda*I=",self
                print "det(A-lambdaI)=",self.determinant()
                print "Q=",Q
                print "R=",R
                
            if verbose>0:
                print "prec(self)=",self._prec
                print "prec(Q)=",Q.prec()
                print "prec(R)=",R.prec()
            assert R.is_upper_triangular()==1
            R.set_zero_elements(tol) # Set all numerically zero elements to zero
            R.echelonize()
            #print "R=",R
            #print "ev=",eigenvec
            dim=self._ncols-R.rank() # Number of parameters
            if verbose>1:
                print "echelon(R)=",R
            if verbose>0:                
                print "dim=",dim
            if dim==0 and abs(z)>0:
                # If this is the case then A-z*I appears to have full rank (which it should not)
                # To redeem this we try to increase precision.
                if depth<10:
                    ## Substitute back and retry:
                    for i from 0<=i<self._nrows:
                        mpc_add(self._matrix[i][i],self._matrix[i][i],z.value,rnd)
                    prec00=self._prec
                    self.set_prec(self._prec+10*self._nrows)
                    if verbose>0:
                        print "Trying to get eigenvectors with higher precision!"
                    vecs = self.eigenvectors(check=check,overwrite=overwrite,sorted=sorted,verbose=verbose,depth=depth+1,old_tol=tol,old_eps=new_eps)
                    res2=[]
                    self.set_prec(prec00)
                    if verbose>0:
                        print "reset prec0 to:{0}".format(self._prec)
                    for i in range(len(vecs)):
                        z2 = self._base_ring(vecs[i][0])
                        vecs2=[]
                        for j in range(len(vecs[i][1])):
                            v = vecs[i][1][j]                
                            v.set_prec(prec0)
                        vecs2.append(v)
                        res2.append((z2,vecs2))
                    return res2
                
                raise ArithmeticError,"Could not compute eigenvector to sufficient precision! A={0} dim={1}, ev={2}".format(self,dim,z)
            if verbose>0:
                print "dim=",dim
            vecs=[]
            if dim==0:  ### We have a unique solution [0,...,0]
                eigenvec=Vector_complex_dense(self.row(0).parent(),0)
                vecs.append(eigenvec)
            #cdef Matrix_complex_dense P
            #P = Matrix_complex_dense(self.parent(),0)
                ### We might need to do a permutation of the variables
            pivots = R.pivots()
            #perms = {}
            swaps={}
            for j in range(R.ncols()):
                swaps[j]=j
            for j in range(R.rank()):
                if pivots[j]<>j:
                    ## Need to swap columns (and record what we did swap)
                    jj = swaps[j]
                    swaps[j]=swaps[pivots[j]]
                    swaps[pivots[j]]=jj
                    R.swap_columns(j,swaps[j])
                    if verbose>0:
                        print "Swap: {0} and {1}".format(j,pivots[j])
            
            for i from 0<=i<dim:

                eigenvec=Vector_complex_dense(self.row(0).parent(),0)
                for j from 0<=j<dim:
                    if j<>i:
                        eigenvec[self._ncols-j-1]=self._base_ring(0)
                    else:
                        eigenvec[self._ncols-j-1]=self._base_ring(1)
                # Now solve the dependent variables
                for j from 0<=j<R.rank():
                    #print "Rank-j-1=",R.rank()-j-1
                    if verbose>1:
                        print "R[",R.rank()-j-1,"][",self._ncols-i-1,"]=",R[R.rank()-j-1][self._ncols-i-1]
                    eigenvec[R.rank()-j-1]=-R[R.rank()-j-1][self._ncols-i-1]
                if verbose>0:
                    if len(swaps)>0:
                        print "swaps={0}".format(swaps)
                    print "eigenvector=",eigenvec
                ## See if we need to swap back
                ## TODO: use a vector or pointers to do the swap and inverse more efficient
                for j in range(R.ncols()):
                    if swaps[j]==j:
                        continue
                    for jj in range(j,R.ncols()):
                        if swaps[jj]==j:
                            if verbose>1:
                                print "swap[{0}]={1}".format(jj,j)
                                print "swapping: {0} and {1}".format(eigenvec[j],eigenvec[jj])
                            z2 = eigenvec[j]
                            eigenvec[j]=eigenvec[jj]
                            eigenvec[jj]=z2
                            #ii = swaps[jj]
                            #swaps[jj]=swaps[j]
                            #swaps[j]=j
                            break
                if verbose>0:
                    if len(swaps)>0:
                        print "after swap: swaps={0}".format(swaps)
                        #print "after swap: eigenvector=",eigenvec
                vecs.append(eigenvec)  
            #mpc_set(z.value,evs[evi],rnd)
            if verbose>0:
                for v in vecs:
                    #if verbose>1:
                    #    print "v=",v
                    #print "|(A-lambda*I)v|=",self*v-
                    print "||QR-A||=",(Q*R-self).norm()
            for i from 0<=i<self._nrows:
                mpc_add(self._matrix[i][i],self._matrix[i][i],z.value,rnd)
            if check>=1:
                for v in vecs:
                    #if verbose>1:
                    #    print "type(self*v)=",type(self*v),(self*v).list()
                    #    print "v*z=",type(v*z),(v*z).list()
                    test=self*v-v*z
                    if verbose>0:
                        print "|self*z-v*z|=",test.norm()
                        if verbose>1:
                            print "v=",v
                            print "z=",z
                            print "self*v-v*z=",test
                            print "R*v=",R*v
                            print "(QR)*v=",(Q*R)*v
                            print "self=",self
                        print "tol*||v||=",tol*v.norm()
                    #if test.norm()<=self._ncols*eps:
                    if test.norm()<=tol*v.norm():
                        continue
                    if depth<10:
                        ## Substitute back and retry:
                        #for i from 0<=i<self._nrows:
                        #    mpc_add(self._matrix[i][i],self._matrix[i][i],z.value,rnd)
                        prec00=self._prec
                        self.set_prec(self._prec+10*self._nrows)
                        if verbose>0:
                            print "Trying to get eigenvectors with higher precision!  new prec:{0}".format(self._prec)
                        vecs = self.eigenvectors(check=check,overwrite=overwrite,sorted=sorted,verbose=verbose,depth=depth+1,old_tol=tol,old_eps=new_eps)
                        res2=[]
                        if verbose>0:
                            print "reset prec1 to:{0}".format(prec00)
                        self.set_prec(prec00)
                        for i in range(len(vecs)):
                            z2 = self._base_ring(vecs[i][0])
                            vecs2=[]
                            for j in range(len(vecs[i][1])):
                                v = vecs[i][1][j]                
                                v.set_prec(prec0)
                            vecs2.append(v)
                            res2.append((z2,vecs2))
                        return res2
                
                        #return vecs
                    raise ArithmeticError,"Could not compute eigenvector to sufficient precision! A={0} dim={1}, ev={2}".format(self,dim,z)
                #if self._verbose>0:
                #    print "test=",test
                #    print "test.norm()=",test.norm()
                #                    raise ArithmeticError,"Could not compute eigenvector to sufficient precision!"
            #eigenvec=self.solve(rhs,overwrite=0)
            ## We should return vectors with the same precision as we started with
            if verbose>0:
                print "reset prec2 to:{0}".format(self._prec)
            vecs2=[]
            for i in range(len(vecs)):
                v = vecs[i]                
                v.set_prec(prec0)
                vecs2.append(v)
                z2 = self._base_ring(z)
            res.append((z2,vecs2))
            # Reset self
        if check>=1:  ## We also check that we have the correct number of eigenvalues (with mult)
            i = 0
            for z,l in res:
                i+=len(l)
            if i<>self.ncols():
                if verbose>1:
                    print "res={0}".format(res)
                raise ArithmeticError,"We got wrong number of eigenvalues! Got:{0} and wanted:{1}".format(i,self.ncols())
            
            
        if verbose>0:
            print "reset rec of self to:",prec0
        if self._prec>prec0:
            self.set_prec(prec0)
        return res
                       
    cpdef list eigenvalues2(self):
        cdef mpc_t* evs
        cdef list res
        cdef int i 
        cdef MPComplexNumber z
        from sage.all import deepcopy
        raise NotImplementedError
        res  = list()
        #from linalg_complex_dense cimport Eigenvalues_of_M
        z = self._base_ring(1)
        evs = <mpc_t *> sage_malloc(sizeof(mpc_t) * self._nrows)
        for i in xrange(self._nrows):
            mpc_init2(evs[i],self._prec)
        #_eigenvalues_(evs,self._matrix,self._nrows,self._prec,self._mpeps,self._rnd,self._rnd_re)
        for i from 0 <= i < self._nrows:
            z=self._base_ring(1)
            mpc_set(z.value,evs[i],rnd)
            res.append(deepcopy(z))
        res.sort(cmp=my_abscmp)
        return res

    def set_prec(self,int prec):
        r"""
        Change the defualt precision for self.
        """        
        cdef int i,j
        cdef mpc_t z
        from sage.rings.complex_mpc import MPComplexField
        from sage.matrix.all import MatrixSpace
        mpc_init2(z,prec)
        self._prec = prec
        self._base_ring=MPComplexField(prec)
        self._parent = MatrixSpace(self._base_ring,self._nrows,self._ncols)
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                mpc_set(z,self._matrix[i][j],self._rnd)
                mpc_set_prec(self._matrix[i][j],prec)
                mpc_set(self._matrix[i][j],z,self._rnd)
        RF= self._base_ring._base
        self._error_qr=RF(0)
        self._eps = RF(2)**RF(1-self._prec)
        mpfr_init2(self._mpeps,self._prec)
        mpfr_set(self._mpeps,self._eps.value,self._rnd_re)
        mpc_clear(z)

    def to_double(self):
        r"""
        Return a Matrix_complex_double_dense approximation to self.
        """
        if self._double_matrix_is_set==0:            
            dbl_entries=[]
            for i in range(self._nrows):
                for j in range(self._ncols):
                    z = self.get_unsafe(i,j)
                    dbl_entries.append(CDF(z.real(),z.imag()))
            MS = MatrixSpace(CDF,self._nrows,self._ncols)
            self._double_matrix = Matrix(MS,dbl_entries)
        return self._double_matrix
                    
        
    def to_numpy(self):
        return self.to_double().numpy()
        
    
####  Helper functions
        
#from sage.rings.complex_mpc cimport MPComplexField_class

cpdef test_matrix_met(int n=20):
    from sage.rings.complex_mpc import MPComplexField
    from sage.matrix.matrix_space import MatrixSpace
    F=MPComplexField(103)
    
    A=MatrixSpace(F,n).random_element()
    for i in xrange(100):
        B=A.qr_decomp()

cdef _my_sign(mpc_t alpha, mpc_t z, prec,mpc_rnd_t rnd_cplx,mpfr_rnd_t rnd_re):
    r""" Sign of a complex number: sign(z)=exp(iArg(z))
    """
    #cdef mpc_t alpha
    cdef mpfr_t arg,y
    #mpc_init2(alpha,prec)
    mpfr_init2(arg,prec)
    mpfr_init2(y,prec)
    mpc_arg(arg,z,rnd_re)
    mpfr_set_ui(y,0,rnd_re)
    mpc_set_fr_fr(alpha,y,arg,rnd_cplx)
    mpc_exp(alpha,alpha,rnd_cplx)
    mpfr_clear(arg)
    mpfr_clear(y)

    
cpdef my_lexcmp_re(MPComplexNumber a,MPComplexNumber b,int reverse=0):
    r"""
    Compare a and b for complex a and b.
    Using lexicographical ordering starting with real part
    """
    cdef int t1,t2
    cdef int c = mpc_cmp(a.value,b.value)
    cdef int cx = MPC_INEX_RE(c)
    cdef int cy = MPC_INEX_IM(c)
    if cx<0:
        return -1
    elif cx>0:
        return 1
    elif cy>0:
        return 1
    elif cy<0:
        return -1
    return 0

cpdef my_lexcmp_im(MPComplexNumber a,MPComplexNumber b,int reverse=0):
    r"""
    Compare a and b for complex a and b.
    Using lexicographical ordering starting with imaginary part
    """
    cdef int t1,t2
    cdef int c = mpc_cmp(a.value,b.value)
    cdef int cx = MPC_INEX_RE(c)
    cdef int cy = MPC_INEX_IM(c)
    if cy<0:
        return -1
    elif cy>0:
        return 1
    elif cx>0:
        return 1
    elif cx<0:
        return -1
    return 0

cpdef my_abscmp(a,b,reverse=0):
    r"""
    Compare |a| and |b| for complex a and b.
    """
    if reverse==1:
        sig=-1
    else:
        sig=1
    aa=abs(a); ab=abs(b)
    if(aa<ab):
        return -1*sig 
    if(aa>ab):
        return 1*sig
    return 0

cpdef my_abscmp2(a,b,reverse=0):
    r"""
    Compare |a[0]| and |b[0]| for tuples a and b with complex a[0] and b[0].
    """
    if reverse==1:
        sig=-1
    else:
        sig=1
    aa=abs(a[0]); ab=abs(b[0])
    if(aa<ab):
        return -1*sig 
    if(aa>ab):
        return 1*sig
    return 0





cdef _norm_vector(mpfr_t* norm, mpc_t* v,int n,mpfr_rnd_t rnd_re):
    cdef int i
    cdef mpfr_t x
    i = <int> mpc_get_prec(v[0])
    #print "prec=",i
    mpfr_init2(x,i)
    mpfr_set_ui(norm[0],0,rnd_re)
    for i from 0 <= i < n:
        #print "v[",i,"]=",print_mpc(v[i])
        mpc_abs(x,v[i],rnd_re)
        mpfr_mul(x,x,x,rnd_re)
        #print "|x[",i,"]|**2=",print_mpfr(x)
        mpfr_add(norm[0],norm[0],x,rnd_re)
    mpfr_clear(x)
    mpfr_sqrt(norm[0],norm[0],rnd_re)

        
cdef _normalize_vector(mpfr_t *norm, mpc_t* v, int n, mpc_rnd_t rnd,mpfr_rnd_t rnd_re):
    cdef mpfr_t x,s
    cdef int i 
    _norm_vector(norm,v,n,rnd_re)
    for i from 0 <= i < n:
        mpc_div_fr(v[i],v[i],norm[0],rnd)


cdef int _is_hermitian(mpc_t* A,int n,int prec, mpfr_t eps, mpc_rnd_t rnd, mpfr_rnd_t rnd_re):
    r"""
    Is self[i,j]=self[j,i].conjugate() or not?
    """
    cdef mpc_t z,w
    cdef mpfr_t x
    cdef int i,j
    mpc_init2(z,prec)
    mpc_init2(w,prec)
    mpfr_init2(x,prec)
    mpc_clear(z); mpc_clear(w); mpfr_clear(x)
    for i from 0 <= i < n:
        for j from 0 <= j < n:
            mpc_set(z, A[i*n+j],rnd)
            mpc_conj(w,z,rnd)
            mpc_sub(z,w,A[j*n+i],rnd)
            mpc_abs(x,z,rnd_re)
            if mpfr_cmp(x,eps)>0:
                mpc_clear(z); mpc_clear(w); mpfr_clear(x)
                return 0
    mpc_clear(z); mpc_clear(w); mpfr_clear(x)
    return 1

cdef int _is_unitary(mpc_t** A,int n,int prec, mpfr_t eps, mpc_rnd_t rnd, mpfr_rnd_t rnd_re,mpfr_t maxerr,int return_err=0):
    r"""
    """
    cdef mpc_t z,w,s
    cdef mpfr_t x,one
    cdef int i,j,kk
    cdef MPComplexNumber tmpc
    from sage.rings.complex_mpc import MPComplexField
    tmpc = MPComplexField(prec)(0)
    mpc_init2(s,prec);  mpc_init2(z,prec)
    mpc_init2(w,prec);  mpfr_init2(x,prec)
    mpfr_init2(one,prec)
    if return_err:
        mpfr_set_ui(maxerr,0,rnd_re)
    mpfr_set_ui(one,1,rnd_re)
    for i from 0 <= i < n:
        for j from 0 <= j < n:
            mpc_set_ui(s,0,rnd)
            for kk from 0 <= kk < n:
                mpc_conj(z,A[j][kk],rnd)
                mpc_mul(z,z,A[i][kk],rnd)
                mpc_add(s,s,z,rnd)
                #mpc_set(tmpc.value,s,self._rnd)
            #print "QQt[",i,",",j,"]=",tmpc
            if i<>j:
                mpc_abs(x,s,rnd_re)
                if return_err:
                    if mpfr_cmp(x,maxerr)>0:
                        mpfr_set(maxerr,x,rnd_re)
                if mpfr_cmp(x,eps)>0:
                    mpc_set(tmpc.value,s,rnd)
                    print "|Qt[",i,",",j,"]|=",abs(tmpc)
                    mpc_clear(s); mpc_clear(z); mpc_clear(w); mpfr_clear(x);mpfr_clear(one)
                    return 0
            if i==j:
                mpc_sub_fr(s,s,one,rnd)
                mpc_abs(x,s,rnd_re)
                if return_err:
                    if mpfr_cmp(x,maxerr)>0:
                        mpfr_set(maxerr,x,rnd_re)
                if mpfr_cmp(x,eps)>0:
                    mpc_set(tmpc.value,s,rnd)
                    print "Qt[",i,",",j,"]=",tmpc
                    mpc_clear(s); mpc_clear(z); mpc_clear(w); mpfr_clear(x);mpfr_clear(one)
                    return 0
    mpc_clear(s); mpc_clear(z); mpc_clear(w); mpfr_clear(x);mpfr_clear(one)
    return 1



cdef int _is_upper_triangular(mpc_t** A,int n,int prec, mpfr_t eps, mpc_rnd_t rnd, mpfr_rnd_t rnd_re,mpfr_t maxerr,int return_err=0):
    r"""
    """
    cdef mpc_t z
    cdef mpfr_t x
    cdef int i,j
    cdef MPComplexNumber tmpc
    from sage.rings.complex_mpc import MPComplexField
    tmpc = MPComplexField(prec)(0)
    mpc_init2(z,prec); mpfr_init2(x,prec)
    if return_err:
        mpfr_set_ui(maxerr,0,rnd_re)
    for i from 1<=i < n:
        for j from 0 <= j < i:
            mpc_set(z, A[i][j],rnd)
            mpc_abs(x,z,rnd_re)
            if return_err:
                if mpfr_cmp(x,maxerr)>0:
                    mpfr_set(maxerr,x,rnd_re)
            if mpfr_cmp(x,eps)>0:
                mpc_set(tmpc.value,z,rnd)
                #print "|R[",i,",",j,"]|=",tmpc
                mpc_clear(z); mpfr_clear(x)
                return 0

    mpc_clear(z); mpfr_clear(x)
    return 1



cdef int _is_hessenberg(mpc_t* A,int n,int prec, mpfr_t eps, mpc_rnd_t rnd, mpfr_rnd_t rnd_re,mpfr_t maxerr,int return_err=0):
    r"""
    """
    cdef mpc_t z
    cdef mpfr_t x
    cdef int i,j
    #cdef MPComplexNumber tmpc
    #from sage.rings.complex_mpc import MPComplexField
    #tmpc = MPComplexField(prec)(0)
    mpc_init2(z,prec)
    mpfr_init2(x,prec)
    if return_err:
        mpfr_set_ui(maxerr,0,rnd_re)
    for j from 1<=j < n:
        for i from j+2 <= i < n:
            mpc_set(z, A[i*n+j],rnd)
            mpc_abs(x,z,rnd_re)
            if return_err:
                if mpfr_cmp(x,maxerr)>0:
                    mpfr_set(maxerr,x,rnd_re)
            if mpfr_cmp(x,eps)>0:
                #mpc_set(tmpc.value,z,rnd)
                #print "|R[",i,",",j,"]|=",tmpc
                mpc_clear(z); mpfr_clear(x)
                return 0
    mpc_clear(z); mpfr_clear(x)
    return 1
    

                     
#    cpdef MPComplexNumber _wilkinson_shift(self,MPComplexNumber a,MPComplexNumber b,MPComplexNumber c,MPComplexNumber d):




cdef int _pivot_element(int k, int r, mpc_t** A,int  *nrows,int *ncols,int *prec,mpfr_rnd_t rnd_re):
    r"""
    Find the pivot row in column k and row >= r
    I.e. the row number i>=k with maximum abs. value self[i][k]
    """
    # TODO: use double precsion for speeed ?
    cdef mpfr_t xmax,xtmp
    cdef int i
    mpfr_init2(xmax,prec[0])
    mpfr_init2(xtmp,prec[0])
    mpfr_set_ui(xmax,0,rnd_re)
    kpiv=-1
    for i from r<=i <nrows[0]:
        mpc_abs(xtmp,A[i][k],rnd_re)
        if mpfr_cmp(xtmp,xmax) > 0 :
            mpfr_set(xmax,xtmp,rnd_re)
            kpiv=i
    mpfr_clear(xmax);mpfr_clear(xtmp)
    return kpiv



cdef int _pivot_element2(int k, int r, mpc_t* A,int  nrows,int ncols,int prec,mpfr_rnd_t rnd_re):
    r"""
    Find the pivot row in column k and row >= r
    I.e. the row number i>=k with maximum abs. value self[i][k]
    """
    # TODO: use double precsion for speeed ?
    cdef mpfr_t xmax,xtmp
    cdef int i
    #cdef RealNumber tmpr
    #tmpr=RealField(prec)(0)
    mpfr_init2(xmax,prec)
    mpfr_init2(xtmp,prec)
    mpfr_set_ui(xmax,0,rnd_re)
    kpiv=-1
    for i from r<=i <nrows:
        mpc_abs(xtmp,A[i*nrows+k],rnd_re)
        #mpfr_set(tmpr.value,xtmp,rnd_re)
        #print "|A[",i,",",i*nrows+k,"]|=",tmpr
        if mpfr_cmp(xtmp,xmax) > 0 :
            ##print "kpiv=",i
            mpfr_set(xmax,xtmp,rnd_re)
            kpiv=i
            #mpfr_set(tmpr.value,xmax,rnd_re)
            #print "|xmax|=",tmpr
            # permuting rows k and kpiv
    #print "kpiv=",kpiv
    mpfr_clear(xmax);mpfr_clear(xtmp)
    return kpiv




cdef _multiply_classical(mpc_t* AB, mpc_t* A,mpc_t* B,int m, int n,int l,int prec,mpc_rnd_t rnd):
        """
        Multiply matrices A,B given that they
        are represented as m*n and n*l  vectors of mpc_t pointers
        with the same precision.
        Checks must be performed *before* calling this routine!
        INPUT:
        -- AB - m x l matrix: <mpc_t *> sage_malloc(sizeof(mpc_t) * m*l)
        -- A - m x n matrix: <mpc_t *> sage_malloc(sizeof(mpc_t) * m*n)
        -- B - n x l matrix: <mpc_t *> sage_malloc(sizeof(mpc_t) * n*l)
        """
        #print "In multiply classical!!"
        cdef int i,j,k,ii
        cdef mpc_t s, z        
        #AB = <mpc_t *> sage_malloc(sizeof(mpc_t) * m*l)
        mpc_init2(z,prec)
        for i from 0 <= i < m:
            for j from 0 <= j < l:
                ii=i*l+j
                mpc_init2(AB[ii],prec)
                mpc_set_ui(AB[ii],0,rnd)
                for k from 0 <= k < n:
                    mpc_mul(z,A[i*n+k],B[k*l+j],rnd)
                    mpc_add(AB[ii],AB[ii],z,rnd)
        mpc_clear(z)

cdef Gershgorin_disks(mpc_t* B, int n, int prec,mpfr_rnd_t rnd_re):
        r"""
        Compute the G.-disks forcontain  self.
        Each eigenvalue of self is contained in one of these disks.
        If they are disjoint then they each contain precisely
        one eigenvalue.
        """
        cdef int i
        cdef mpfr_t r,x
        cdef RealNumber rr
        from sage.all import deepcopy
        mpfr_init2(r,prec)
        mpfr_init2(x,prec)
        rr = RFF(prec)(1)
        res=dict()
        for i from 0 <= i < n:
            mpfr_set_ui(r,0,rnd_re)
            for j from 0 <= j < n:
                if j<>i:
                    mpc_abs(x,B[i*n+j],rnd_re)
                    mpfr_add(r,r,x,rnd_re)
            mpfr_set(rr.value,r,rnd_re)
            res[i]=deepcopy(rr)
        mpfr_clear(r);mpfr_clear(x)
        return res

            
 
##### to be removed



cdef int _is_unitary2(mpc_t* A,int n,int prec, mpfr_t eps, mpc_rnd_t rnd, mpfr_rnd_t rnd_re,mpfr_t maxerr,int return_err=0):
    r"""
    """
    cdef mpc_t z,w,s
    cdef mpfr_t x,one
    cdef int i,j,kk
    cdef MPComplexNumber tmpc
    from sage.rings.complex_mpc import MPComplexField
    tmpc = MPComplexField(prec)(0)
    mpc_init2(s,prec);  mpc_init2(z,prec)
    mpc_init2(w,prec);  mpfr_init2(x,prec); mpfr_init2(one,prec)
    if return_err:
        mpfr_set_ui(maxerr,0,rnd_re)
    mpfr_set_ui(one,1,rnd_re)
    for i from 0 <= i < n:
        for j from 0 <= j < n:
            mpc_set_ui(s,0,rnd)
            for kk from 0 <= kk < n:
                mpc_conj(z,A[j*n+kk],rnd)
                mpc_mul(z,z,A[i*n+kk],rnd)
                mpc_add(s,s,z,rnd)
                #mpc_set(tmpc.value,s,self._rnd)
            #print "QQt[",i,",",j,"]=",tmpc
            if i<>j:
                mpc_abs(x,s,rnd_re)
                if return_err:
                    if mpfr_cmp(x,maxerr)>0:
                        mpfr_set(maxerr,x,rnd_re)
                if mpfr_cmp(x,eps)>0:
                    mpc_set(tmpc.value,s,rnd)
                    print "Qt[",i,",",j,"]=",tmpc
                    mpc_clear(s);mpc_clear(z);mpc_clear(w)
                    mpfr_clear(x); mpfr_clear(one)
                    return 0
            if i==j:
                mpc_sub_fr(s,s,one,rnd)
                mpc_abs(x,s,rnd_re)
                if return_err:
                    if mpfr_cmp(x,maxerr)>0:
                        mpfr_set(maxerr,x,rnd_re)
                if mpfr_cmp(x,eps)>0:
                    mpc_set(tmpc.value,s,rnd)
                    print "Qt[",i,",",j,"]=",tmpc
                    mpc_clear(s);mpc_clear(z);mpc_clear(w)
                    mpfr_clear(x); mpfr_clear(one)
                    return 0
    mpc_clear(s);mpc_clear(z);mpc_clear(w)
    mpfr_clear(x); mpfr_clear(one)
    return 1


cdef int _is_upper_triangular2(mpc_t* A,int n,int prec, mpfr_t eps, mpc_rnd_t rnd, mpfr_rnd_t rnd_re,mpfr_t maxerr,int return_err=0):
    r"""
    """
    cdef mpc_t z
    cdef mpfr_t x
    cdef int i,j
    cdef MPComplexNumber tmpc
    from sage.rings.complex_mpc import MPComplexField
    tmpc = MPComplexField(prec)(0)
    mpc_init2(z,prec)
    mpfr_init2(x,prec)
    if return_err:
        mpfr_set_ui(maxerr,0,rnd_re)
    for i from 1<=i < n:
        for j from 0 <= j < i:
            mpc_set(z, A[i*n+j],rnd)
            mpc_abs(x,z,rnd_re)
            if return_err:
                if mpfr_cmp(x,maxerr)>0:
                    mpfr_set(maxerr,x,rnd_re)
            if mpfr_cmp(x,eps)>0:
                mpc_set(tmpc.value,z,rnd)
                print "|R[",i,",",j,"]|=",tmpc
                mpfr_clear(x);mpc_clear(z)
                return 0

    mpfr_clear(x);mpc_clear(z)        
    return 1

def random_matrix_eigenvalues(F,n):
    r"""
    Give a random matrix together with its eigenvalues.
    """
    from sage.all import MatrixSpace,copy
    l=list()
    M = MatrixSpace(F,n)
    #U = Matrix(ring=F,n)
    D = copy(M.zero())
    for i in xrange(n):
        x = F.random_element()
        l.append(x)
        D[i,i]=x
    # now we need a unitary matrix:
    # make a random choice of vectors and use Gram-Schmidt to orthogonolize
    U = random_unitary_matrix(F,n)
    UT = U.transpose().conjugate()
    A = U*D*UT
    l.sort(cmp=my_abscmp)
    return   A,U,l



def random_unitary_matrix(F,n):
    r"""
    Give a random matrix together with its eigenvalues.
    """
    from sage.all import VectorSpace
    from sage.matrix.matrix_space import MatrixSpace 
    M = MatrixSpace(F,n)    
    U=Matrix_complex_dense.__new__(Matrix_complex_dense,M,None,None,None)

    # make a random choice of vectors and use Gram-Schmidt to orthogonolize
    V = VectorSpace(F,n)
    v = dict(); u=dict()
    for i in xrange(n):
        v[i] = V.random_element()
        nv= v[i].norm()
        v[i]=v[i]/nv
    # check that the v's are indeed linear independent
    VV= V.vector_space_span(v.values())
    assert VV.dimension()==V.dimension()
    for i in xrange(n):
        u[i] = v[i]
        for j in xrange(i):
            u[i]=u[i]-u[j]*u[i].scalar_product(u[j])
        u[i]=u[i]/u[i].norm()
    for i in xrange(n):
        for j in xrange(n):
            U[i,j]=u[i][j]
    return U


def RandomComplexMatrix(sz,prec=53,**kwds):
    r"""
    Compute a random matrix of type Matrix_complex_dense
    """
    CF = MPComplexField(prec)
    MS = MatrixSpace(CF,sz,sz)
    A = Matrix_complex_dense(MS,MS.random_element(**kwds).list())
    return A

def test_eigenvalues(num_test=5,sz=10,prec=102,verbose=0):
    r"""
    Test the computation of eigenvalues.
    """
    eps = 2.0**-prec
    if verbose>0:
        print "eps=",eps
    for j in range(num_test):
        A = RandomComplexMatrix(sz,prec)
        ev = A.eigenvalues(check=1)
        tr = sum(ev)
        tr1 = A.trace()
        er = abs(tr-tr1)
        if verbose>0:
            print "er=",er
        assert er < eps

def test_eigenvalues(num_test=5,sz=10,prec=102,verbose=0):
    r"""
    Test the computation of eigenvalues.
    """
    CF = MPComplexField(prec)
    MS = MatrixSpace(CF,sz,sz)
    eps = sz*2.0**(5-prec)
    if verbose>0:
        print "eps=",eps
    for j in range(num_test):
        A = RandomComplexMatrix(sz,prec)
        ev = A.eigenvalues(check=1)
        tr = sum(ev)
        tr1 = A.trace()
        er = abs(tr-tr1)
        if verbose>0:
            print "er=",er
        assert er < eps
