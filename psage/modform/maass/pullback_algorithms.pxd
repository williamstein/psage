include 'sage/ext/stdsage.pxi'
include "sage/ext/cdefs.pxi"
include 'sage/ext/interrupt.pxi'
#include "sage/ext/gmp.pxi"
include "sage/rings/mpc.pxi"
from sage.libs.mpfr cimport *

from psage.modules.vector_real_mpfr_dense cimport Vector_real_mpfr_dense
from sage.rings.real_mpfr cimport RealNumber

cdef int pullback_pts_mpc_new_c(S,int Qs,int Qf,mpfr_t Y, mpfr_t* Xm,mpfr_t*** Xpb,mpfr_t*** ypb,mpc_t ***Cvec)
cdef int pullback_pts_mpc_new_c_sym(S,int Qs,int Qf,mpfr_t Y, mpfr_t* Xm,mpfr_t*** Xpb,mpfr_t*** ypb,mpfr_t ****RCvec,int*** CSvec)


cdef int pullback_pts_cplx_dp(S,int Qs,int Qf,double Y,double *Xm,double *** Xpb,double*** Ypb,double complex ***Cvec)

cdef int pullback_pts_real_dp(S,int Qs,int Qf,double Y,double *Xm,double *** Xpb,double*** Ypb,double ***Cvec)

cdef int pullback_pts_cplx_dp_sym(S,int **Qv,double Y,double *Xm,double *** Xpb,double*** Ypb,double complex ***Cvec)
