include 'sage/ext/stdsage.pxi'
include "sage/ext/cdefs.pxi"
include 'sage/ext/interrupt.pxi'
#include "sage/ext/gmp.pxi"
include "sage/rings/mpc.pxi"
from sage.libs.mpfr cimport *

from psage.modules.vector_real_mpfr_dense cimport Vector_real_mpfr_dense
from sage.rings.real_mpfr cimport RealNumber

cdef pullback_pts_mpc_new_c(S,int Qs,int Qf,mpfr_t Y, mpfr_t* Xm,mpfr_t*** Xpb,mpfr_t*** ypb,mpc_t ***Cvec)
cdef pullback_pts_mpc_new_c_sym(S,int Qs,int Qf,mpfr_t Y, mpfr_t* Xm,mpfr_t*** Xpb,mpfr_t*** ypb,mpfr_t ****RCvec,int*** CSvec)


cdef void pullback_pts_cplx_dp(S,int Qs,int Qf,double Y,double *Xm,double *** Xpb,double*** Ypb,double complex ***Cvec,double weight=*,int holo=*)

cdef void pullback_pts_real_dp(S,int Qs,int Qf,double Y,double *Xm,double *** Xpb,double*** Ypb,double ***Cvec,double weight=*,int holo=*)

