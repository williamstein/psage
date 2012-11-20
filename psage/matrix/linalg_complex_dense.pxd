include "sage/ext/stdsage.pxi"
include "sage/ext/cdefs.pxi"

#include "sage/rings/mpc.pxi"

#from psage.modules.vector_complex_dense cimport *

#cdef void givens(mpc_t *s, mpc_t *skk, mpfr_t *c, mpc_t *r, mpc_t *f, mpc_t *g, mpc_t *t, mpc_rnd_t rnd,mpfr_rnd_t rnd_re)

#cdef void RotateLeft_(mpc_t s, mpc_t skk, mpfr_t c,mpc_t** A, int i, int k, int min_j,int max_j, mpc_t *t,mpc_rnd_t rnd, mpfr_rnd_t rnd_re)

#cdef void RotateRight_(mpc_t s, mpc_t skk, mpfr_t c,mpc_t** A, int i, int k, int min_j, int max_j, mpc_t *t,mpc_rnd_t rnd, mpfr_rnd_t rnd_re)

from psage.rings.mpc_extras cimport *

ctypedef struct QR_set:
    mpc_t s
    mpc_t skk
    mpfr_t c
    mpc_t s_t
    mpc_t skk_t
    mpfr_t c_t
    mpc_t r
    mpc_t t[5]
    mpfr_t x[4]
    mpc_t mu1
    mpc_t mu2

cdef QR_set init_QR(prec)

#cdef print_mat(mpc_t **A,int m,int n)

cdef void qr_decomp(mpc_t** A,int m, int n, int prec, mpfr_t eps, mpc_rnd_t rnd, mpfr_rnd_t rnd_re)

cdef void _eigenvalues(mpc_t* res, mpc_t** A,int nrows,int prec,mpc_rnd_t rnd,mpfr_rnd_t rnd_re)

cdef _reconstruct_matrix(mpc_t** Q, mpc_t** A, int n, int m, int k, int prec, mpc_rnd_t rnd, mpfr_rnd_t rnd_re)

cdef void _hessenberg_reduction(mpc_t** A, int nrows, QR_set  q,int prec, mpc_rnd_t rnd,mpfr_rnd_t rnd_re, int rt =?)

cdef void _norm(mpfr_t norm,mpc_t** A, int nrows,int ncols, int prec,mpfr_rnd_t rnd_re,int ntype)

cdef void solve_upper_triangular(mpc_t** R,mpc_t* b,int n, int prec, mpc_rnd_t rnd, mpfr_rnd_t rnd_re)

## cdef void _wilkinson_shift(mpc_t* mu,mpc_t* a,mpc_t* b,mpc_t* c,mpc_t* d,int prec,mpc_rnd_t rnd,mpfr_rnd_t rnd_re, mpc_t zz[4], mpfr_t xx[4])


## cdef inline mpc_zero_p(mpc_t z):
##     return mpfr_zero_p(z.re) and mpfr_zero_p(z.im)




## cdef inline void _mpc_div(mpc_t * z, mpc_t a, mpc_t b,  mpc_t t[2], mpfr_rnd_t rnd_re)

## cdef inline void _mpc_mul(mpc_t* z, mpc_t a, mpc_t b, mpc_t *t, mpc_rnd_t rnd, mpfr_rnd_t rnd_re)

## cdef inline void _mpc_mul_fr(mpc_t* z, mpc_t a, mpfr_t b, mpc_rnd_t rnd, mpfr_rnd_t rnd_re)

## cdef inline void _mpc_add(mpc_t *res, mpc_t a, mpc_t b,mpfr_rnd_t rnd_re)

## cdef inline void _mpc_sub(mpc_t *res, mpc_t a, mpc_t b,mpfr_rnd_t rnd_re)

## cdef inline void _mpc_set(mpc_t *res, mpc_t z, mpfr_rnd_t rnd_re)

## cdef print_mpc(mpc_t z)

## cdef print_mpfr(mpfr_t z)
