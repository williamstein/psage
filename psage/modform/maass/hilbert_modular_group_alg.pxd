include "sage/ext/interrupt.pxi" 
include "sage/ext/stdsage.pxi"  
include "sage/ext/cdefs.pxi"
include "sage/rings/mpc.pxi"

cdef class Hn_dble(object):
cdef class Hn(object):
    cdef double *_x
    cdef double *_y
    cdef int _degree
    cdef list _xlist
    cdef list _ylist
    cdef int _prec
    cdef int _verbose
    cdef double complex _norm
    cdef double _imag_norm
    cdef double _real_norm
    cdef int _norm_set
    cdef int _imag_norm_set
    cdef int _real_norm_set
    
    cdef _free(self)
    cdef c_new(self,list x,list y)
    cpdef x(self,int i)
    cpdef y(self,int i)
    cpdef z(self,int i)
    cpdef real_list(self)
    cpdef imag_list(self)
    cpdef imag(self)
    cpdef real(self)
    cpdef imag_norm(self)
    cpdef real_norm(self)
    cpdef norm(self)
    cpdef vector_norm(self)
    cpdef imag_trace(self)
    cpdef real_trace(self)
    cpdef trace(self)
    cpdef degree(self)
    cpdef diff_abs_norm(self,z)
    cpdef addto_re(self,x)
    cdef double *_x
    cdef double *_y
    cdef int _degree
    cdef list _xlist
    cdef list _ylist
    cdef int _prec
    cdef int _verbose
    cdef double complex _norm
    cdef double _imag_norm
    cdef double _real_norm
    cdef int _norm_set
    cdef int _imag_norm_set
    cdef int _real_norm_set
    
    cdef _free(self)
    cdef c_new(self,list x,list y)
    cpdef x(self,int i)
    cpdef y(self,int i)
    cpdef z(self,int i)
    cpdef real_list(self)
    cpdef imag_list(self)
    cpdef imag(self)
    cpdef real(self)
    cpdef imag_norm(self)
    cpdef real_norm(self)
    cpdef norm(self)
    cpdef vector_norm(self)
    cpdef imag_trace(self)
    cpdef real_trace(self)
    cpdef trace(self)
    cpdef degree(self)
    cpdef diff_abs_norm(self,z)
    cpdef addto_re(self,x)
    cpdef is_in_upper_half_plane(self)

cdef class Hn(object):
    cdef mpfr_t *_x
    cdef mpfr_t *_y
    cdef mpc_t *_z
    cdef int _degree
    cdef list _xlist
    cdef list _ylist
    cdef int _prec
    cdef int _verbose
    cdef mpc_t _norm
    cdef mpfr_t _imag_norm
    cdef mpfr_t _real_norm
    cdef int _norm_set
    cdef int _imag_norm_set
    cdef int _real_norm_set
    
    cdef _free(self)
    cdef c_new(self,list x,list y)
    cpdef x(self,int i)
    cpdef y(self,int i)
    cpdef z(self,int i)
    cpdef real_list(self)
    cpdef imag_list(self)
    cpdef imag(self)
    cpdef real(self)
    cpdef imag_norm(self)
    cpdef real_norm(self)
    cpdef norm(self)
    cpdef vector_norm(self)
    cpdef imag_trace(self)
    cpdef real_trace(self)
    cpdef trace(self)
    cpdef degree(self)
    cpdef diff_abs_norm(self,z)
    cpdef addto_re(self,x)
    cpdef addto_im(self,x)
