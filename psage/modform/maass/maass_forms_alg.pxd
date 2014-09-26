cdef void set_Mv_Qv_symm(S,int **Mv,int **Qv,double *Qfak,int * symmetric_cusps,double complex *cusp_evs,int *cusp_offsets,int *N,int *Ml,int *Ql,int M,int Q,int verbose=?)

cdef int get_M_for_maass_dp_c(double R,double Y,double eps)
cdef double err_est_Maasswf_dp_c(double Y,int M,double R,int pref=?)
#cpdef get_Y_for_M_dp(S,double R,int M,double eps,int verbose=?)

cdef void set_Mv_Qv_symm(S,int **Mv,int **Qv,double *Qfak,int * symmetric_cusps,double complex* cusp_evs,int *cusp_offsets,int *N,int *Ml,int *Ql,int M,int Q,int verbose=?)
cdef void set_Mv_Qv_symm_real(S,int **Mv,int **Qv,double *Qfak,int * symmetric_cusps,double *cusp_evs,int *cusp_offsets,int *N,int *Ml,int *Ql,int M,int Q,int verbose=?)

cpdef get_M_for_maass_dp(double R,double Y,double eps)

#cpdef get_Y_for_M_dp(S,double R,int M,double eps,int verbose=?)

cdef SMAT_cplx_dp(double complex** U,int N,int num_rhs,int num_set,double complex **C,double complex** values,int* setc)


#cpdef get_Y_and_M_dp(S,double R,double eps,int verbose=?)
