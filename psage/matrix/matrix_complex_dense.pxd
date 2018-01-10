#include "sage/ext/stdsage.pxi"
include "sage/ext/cdefs.pxi"
#include "sage/rings/mpc.pxi"
#include "sage/ext/gmp.pxi"
from psage.rings.mpfr_nogil cimport *
from psage.modules.vector_complex_dense cimport Vector_complex_dense
from sage.matrix.matrix_dense cimport Matrix_dense
from sage.structure.element cimport Vector,Element
from sage.rings.complex_mpc cimport MPComplexField_class,MPComplexNumber
from sage.rings.real_mpfr cimport RealNumber
from sage.matrix.matrix_complex_double_dense cimport Matrix_complex_double_dense

cdef class Matrix_complex_dense(Matrix_dense):

    #cdef mpc_t tmp
    cpdef _add_(self, other)
    cpdef _sub_(self, other)
    
    cdef mpc_t *_entries
    cdef int _entries_are_allocated
    cdef mpc_t ** _matrix
    cdef object __pivots
    cdef Vector _matrix_times_vector_(matrix_left, Vector v_right)
    cdef Vector _matrix_times_vector_complex_dense(matrix_left, Vector_complex_dense v_right)
    cdef Vector _vector_times_matrix_(matrix_right, Vector v_left)
    cdef int _prec
    cdef int _verbose
    cdef RealNumber _eps
    cdef RealNumber _error_qr
    cdef mpfr_t _mpeps
    #cpdef RealNumber eps(self)
    cdef int _is_square
    cdef mpc_rnd_t _rnd      # complex rounding mode
    cdef mpfr_rnd_t _rnd_re # real base class rounding mode
    cdef mpfr_rnd_t _rnd_im # real base class rounding mode
    cdef int _base_for_str_rep
    cpdef int base_for_str_rep(self)
    cpdef int set_base_for_str_rep(self,int base)
    cdef int _truncate
    cdef int _double_matrix_is_set
    cdef Matrix_complex_double_dense _double_matrix
    #    cpdef Matrix_complex_dense zero_matrix()
    #cpdef Matrix_complex_dense identity_matrix()
    #cdef MPComplexField_class _base_ring
    #cdef object _base_ring
    cdef int _cmp_c_impl(self, Element right) except -2

    cdef mpc_t* _transformation_to_hessenberg
    #cpdef Matrix_complex_dense transformation_to_hessenberg(self)
    
    # cdef int _rescale(self, mpz_t a) except -1
    #cdef pickle(self)
    cdef _pickle_version0(self)
    cdef _unpickle_version0(self, data)    
    cpdef _export_as_string(self, int base=?,int truncate=?)
    cpdef _reqdigits(self,int base=?,int truncate=?)

    cpdef int _pivot_element(self,int k,int r)
    cpdef RealNumber error_qr(self)
    cpdef tuple qr_decomposition(self, int overwrite=?, int check=?,int num_threads=?,int schedule=?)
    cpdef hessenberg(self,int return_transformation=?, int check=?,int num_threads=?,int schedule=?)

    cpdef Vector_complex_dense row(self,int n)
    cpdef Vector_complex_dense column(self,int n)
    cdef void _column(self,mpc_t *colv,int n)
    cpdef delete_column(self,int n)
    cpdef delete_row(self,int n,int clear=?)

    cpdef  set_zero_elements(self,double tol=?)
    cpdef int numerical_rank(self,double tol=?)
    cpdef _balance(self)
    cpdef int is_hessenberg(self,double maxerr=?,int show_err=?)

    cpdef int is_upper_triangular(self,int return_err=?)
    cpdef int is_square(self)
    cpdef int is_orthogonal(self)
    cpdef int is_unitary(self,int return_err=?)
    cpdef int is_hermitian(self)
    cpdef RealNumber norm(self,int ntype=*)
    cdef void _norm(self,mpfr_t,int ntype)
    cpdef dict _norms
    cpdef list eigenvalues(self,int check=?,int sorted=?,int overwrite=?,int num_threads=?,int schedule=?)
    cpdef list eigenvectors(self,int check=?,int overwrite=?,int sorted=?,int verbose=?,int depth=?,double old_tol=?,double old_eps=?,int num_threads=?,int schedule=?)
    cpdef list eigenvalues2(self)
    cpdef conjugate_transpose(self)
    cdef _conjugate_transpose(self, Matrix_complex_dense res)    
    cpdef transpose(self)
    cdef _transpose(self, Matrix_complex_dense res)     
    cpdef det(self)
    cpdef prec(self)
  
    cpdef Vector_complex_dense solve(self,Vector_complex_dense b,int overwrite=?,int num_threads=?,int schedule=?)
    cpdef Matrix_complex_dense mat_solve(self,Matrix_complex_dense B,int overwrite=*)
    cpdef Matrix_complex_dense inverse(self,int overwrite=?)

#    cpdef tuple qr_decomp2(self,int check=?)
  
    #cdef RotateLeft(self,mpc_t s,mpcl_t skk, mpfr_t c,int i,int k,mpc_t t[3])
        
    #cdef RotateRight(self,mpc_t s,mpc_t skk, mpfr_t c,int i,int k,mpc_t t[3])



    #cdef _add_ui_unsafe_assuming_int(self, Py_ssize_t i, Py_ssize_t j, unsigned long int n)
    #cdef _sub_ui_unsafe_assuming_int(self, Py_ssize_t i, Py_ssize_t j, unsigned long int n)
    
#cdef class MatrixWindow:
#    cdef Matrix_complex_dense _matrix
#    cdef int _row, _col, _nrows, _ncols

################################################################
# fast conversion to pari on the stack
################################################################
## ctypedef long* GEN
## cdef GEN pari_GEN(Matrix_complex_dense B)
cpdef test_matrix_met(int n=?)


