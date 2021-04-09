from sage.libs.mpfr cimport *
#include "sage/rings/mpc.pxi"
#include "sage/ext/gmp.pxi"

from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense

cdef class GL2Z_elt(object):
    cdef int* ent
    cdef int det
    cdef str _str
    cpdef a(self)
    cpdef b(self)
    cpdef c(self)
    cpdef d(self)
    cpdef determinant(self)
    cpdef is_in_Gamma0N(self,int N)
    cpdef acton(self,z)
    cpdef _mul(self,GL2Z_elt other,int inv=?)
    cpdef _mul_mat_id(self,Matrix_integer_dense other,int inv=?)
    #cpdef _mul_mat_i_2x2(self,Matrix_integer_2x2 other,int inv=?)
    cdef _mul_c(self,int* other,int *res,int inv=?)
    cdef _mul_c_mpz(self,mpz_t* other,int *res,int inv=?)
    cpdef inverse(self)
    cpdef _eq(self,GL2Z_elt other)
    cpdef _pow(self,int k)
    cpdef _conjugate(self,GL2Z_elt other)

cdef class SL2Z_elt(GL2Z_elt):
    pass
    
from sage.structure.sage_object cimport *

cdef void _apply_sl2z_map_mpfr(mpfr_t x,mpfr_t y,int a,int b,int c,int d)
cdef void _apply_gl2z_map_mpfr(mpfr_t x,mpfr_t y,int a,int b,int c,int d)
cdef void _apply_sl2z_map_dp(double *x,double *y,int a,int b,int c,int d)

cdef void _pullback_to_psl2z_mpfr(mpfr_t x,mpfr_t y)

cdef void _pullback_to_Gamma0N_mpfr(int*** reps ,int nreps, int N,mpfr_t x,mpfr_t y)

cdef void _pullback_to_Gamma0N_dp(int*** reps ,int nreps, int N,double *x,double *y,int *a,int *b,int *c,int *d,int verbose=?)

#cdef void _pullback_to_Gamma0N_dp(int*** reps ,int nreps, int N,double *x,double *y,int *a,int *b,int *c,int *d,int verbose)
cdef int closest_vertex_dp_c(int nv,int **vertex_maps,double *widths,double* x,double* y,int verbose=?)



cdef void pullback_to_hecke_triangle_mat_c(double *x,double *y,double lambdaq,double *a,double *b,double *c,double *d)

cdef void pullback_to_psl2z_mat_c(double *x,double *y,int *a,int *b, int *c, int *d,int verbose=?)

cpdef pullback_to_psl2z_mat2(double x,double y,int a,int b,int c,int d)

cpdef tuple pullback_to_psl2z_mat(double x,double y)

cdef void normalize_point_to_cusp_mpfr_c(mpfr_t xout, mpfr_t yout, int * N ,int wi, mpfr_t xin, mpfr_t yin, int inv)

cdef void _normalize_point_to_cusp_mpfr(mpfr_t x,mpfr_t y,int a,int b,int c,int d, int wi, int inv=?)


cdef void pullback_to_Gamma0N_mpfr_c(G,mpfr_t xout,mpfr_t yout, mpfr_t xin,mpfr_t yin,int *a,int *b,int *c,int *d)



cdef void _normalize_point_to_cusp_dp(double *x,double *y,int a,int b,int c,int d, int wi, int inv=?)    

cdef void _normalize_point_to_cusp_real_dp(double *x,double *y,int a,int b,int c,int d, double wi, int inv=?)    

cdef void pullback_to_hecke_triangle_mat_c_mpfr(mpfr_t x,mpfr_t y,mpfr_t lambdaq,mpfr_t a,mpfr_t b,mpfr_t c,mpfr_t d)

cdef tuple fast_sl2z_factor(int a,int b,int c,int d,int verbose=?)
