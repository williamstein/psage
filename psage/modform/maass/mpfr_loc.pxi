cdef extern from "mpfr.h" nogil:
    ctypedef int mpfr_prec_t
    ctypedef unsigned int mpfr_uprec_t
    ctypedef int mpfr_sign_t
    ctypedef short mpfr_exp_t
    ctypedef unsigned short mpfr_uexp
    
    struct __mpfr_struct:
        mpfr_prec_t  _mpfr_prec
        mpfr_sign_t  _mpfr_sign
        mpfr_exp_t   _mpfr_exp
        mp_limb_t   *_mpfr_d

    ctypedef enum mpfr_rnd_t:
        MPFR_RNDN=0,
        MPFR_RNDZ,
        MPFR_RNDU,
        MPFR_RNDD,
        MPFR_RNDA,
        MPFR_RNDF,
        MPFR_RNDNA=-1
        
    ctypedef __mpfr_struct mpfr_t[1]
    ctypedef __mpfr_struct *mpfr_ptr
    
    ctypedef __mpfr_struct *mpfr_srcptr
    
    void mpfr_init2 (mpfr_ptr, mpfr_prec_t)
    void mpfr_init (mpfr_ptr)
    void mpfr_clear (mpfr_ptr)
    int mpfr_const_pi (mpfr_ptr, mpfr_rnd_t)
    int mpfr_mul_ui (mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t)
    int mpfr_sqrt (mpfr_ptr, mpfr_srcptr, mpfr_rnd_t)
    int mpfr_add (mpfr_ptr, mpfr_srcptr,mpfr_srcptr, mpfr_rnd_t)
    int mpfr_get_prec(mpfr_srcptr)
    int mpfr_add_ui (mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t)
    void mpfr_set4 (mpfr_ptr, mpfr_srcptr, mpfr_rnd_t, int)
    void mpfr_set (mpfr_ptr, mpfr_srcptr, mpfr_rnd_t)
    unsigned long mpfr_get_ui (mpfr_srcptr, mpfr_rnd_t)
    int mpfr_set_ui (mpfr_ptr, unsigned long, mpfr_rnd_t)
    int mpfr_cos (mpfr_ptr, mpfr_srcptr,mpfr_rnd_t)
    int mpfr_sin (mpfr_ptr, mpfr_srcptr,mpfr_rnd_t)
    int mpfr_exp (mpfr_ptr, mpfr_srcptr,mpfr_rnd_t)
    double mpfr_get_d (mpfr_srcptr, mpfr_rnd_t)
    int mpfr_const_pi (mpfr_ptr, mpfr_rnd_t)
    int mpfr_pow (mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t)
    int mpfr_pow_si (mpfr_ptr, mpfr_srcptr, long int, mpfr_rnd_t)
    int mpfr_pow_ui (mpfr_ptr, mpfr_srcptr, unsigned long int, mpfr_rnd_t)
    int mpfr_log (mpfr_ptr, mpfr_srcptr,mpfr_rnd_t)
    int mpfr_set_si (mpfr_ptr, long, mpfr_rnd_t)
    int mpfr_cmp  (mpfr_srcptr, mpfr_srcptr)
    int mpfr_cmp3 (mpfr_srcptr, mpfr_srcptr, int)
    int mpfr_cmp_d (mpfr_srcptr, double)
    int mpfr_cmp_ld (mpfr_srcptr, long double)
    int mpfr_cmpabs (mpfr_srcptr, mpfr_srcptr)
    int mpfr_cmp_ui (mpfr_srcptr, unsigned long)
    int mpfr_cmp_si (mpfr_srcptr, long)
    int mpfr_mul (mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t)
    int mpfr_abs (mpfr_ptr, mpfr_srcptr, mpfr_rnd_t)
    int mpfr_neg (mpfr_ptr, mpfr_srcptr, mpfr_rnd_t)
    int mpfr_div_ui (mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t)
    int mpfr_set_si_2exp (mpfr_ptr, long, mpfr_exp_t, mpfr_rnd_t)
    void mpfr_set_prec (mpfr_ptr, mpfr_prec_t)
    int mpfr_div (mpfr_ptr, mpfr_srcptr,mpfr_srcptr, mpfr_rnd_t)
    int mpfr_sqrt (mpfr_ptr, mpfr_srcptr,mpfr_rnd_t)
    int mpfr_erfc (mpfr_ptr, mpfr_srcptr,mpfr_rnd_t)
    int mpfr_const_euler (mpfr_ptr, mpfr_rnd_t)
    int mpfr_set_d(mpfr_ptr, double, mpfr_rnd_t)
    int mpfr_sub_si (mpfr_ptr, mpfr_srcptr, long int, mpfr_rnd_t)
    int mpfr_mul_si (mpfr_ptr, mpfr_srcptr, long int, mpfr_rnd_t)
    int mpfr_zero_p (mpfr_srcptr)
    int mpfr_div_si (mpfr_ptr, mpfr_srcptr, long int, mpfr_rnd_t)
    int mpfr_mul_d (mpfr_ptr, mpfr_srcptr,double, mpfr_rnd_t)
    float mpfr_get_flt (mpfr_srcptr, mpfr_rnd_t)
    int mpfr_zero_p (mpfr_srcptr)
    int mpfr_sub (mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t)
    long mpfr_get_si (mpfr_srcptr, mpfr_rnd_t)

cdef extern from "mpc.h" nogil:
    int  mpc_set_fr (mpc_ptr, mpfr_srcptr, mpc_rnd_t)
    int  mpc_set_si_si (mpc_ptr, long int, long int, mpc_rnd_t)
    int  mpc_mul_fr (mpc_ptr, mpc_srcptr, mpfr_srcptr, mpc_rnd_t)
    int  mpc_add (mpc_ptr, mpc_srcptr, mpc_srcptr, mpc_rnd_t)
    void mpc_init2 (mpc_ptr, mpfr_prec_t)
