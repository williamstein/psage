cdef extern from "mpfr.h" nogil:
    ctypedef long mpfr_prec_t
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
        MPFR_RNDZ=1,
        MPFR_RNDU=2,
        MPFR_RNDD=3,
        MPFR_RND_MAX=4,
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
    #int mpfr_zero_p (mpfr_srcptr)
    int mpfr_sub (mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t)
    long mpfr_get_si (mpfr_srcptr, mpfr_rnd_t)

cdef extern from "mpc.h" nogil:
    #ctypedef void* mpc_t
    ctypedef struct __mpc_struct:
        mpfr_t re
        mpfr_t im
    ctypedef __mpc_struct* mpc_t
    ctypedef mpc_t mpc_ptr
    ctypedef mpc_t mpc_srcptr

    ctypedef enum mpc_rnd_t:
        MPC_RNDNN = 0
        MPC_RNDZN = 1
        MPC_RNDUN = 2
        MPC_RNDDN = 3
        MPC_RNDNZ = 16
        MPC_RNDZZ = 17
        MPC_RNDUZ = 18
        MPC_RNDDZ = 19
        MPC_RNDNU = 32
        MPC_RNDZU = 33
        MPC_RNDUU = 34
        MPC_RNDDU = 35
        MPC_RNDND = 48
        MPC_RNDZD = 49
        MPC_RNDUD = 50
        MPC_RNDDD = 51

    # Memory management
    void mpc_init (mpc_t)
    void mpc_init2 (mpc_t, mpfr_prec_t)
    void mpc_clear (mpc_t)

    # Precision accessors
    mpfr_prec_t mpc_get_prec (mpc_t)
    void mpc_set_prec (mpc_t, mpfr_prec_t)

    # Set real part to given value and imaginary part to +0
    int  mpc_set_ui (mpc_t, unsigned long int, mpc_rnd_t)
    int  mpc_set_si (mpc_t, long int, mpc_rnd_t)
    int  mpc_set_z (mpc_t, mpz_t, mpc_rnd_t)
    int  mpc_set_d (mpc_t, double, mpc_rnd_t)
    int  mpc_set_fr (mpc_t, mpfr_t, mpc_rnd_t)
    # Set value
    int  mpc_set (mpc_t, mpc_t, mpc_rnd_t)
    int  mpc_set_ui_ui (mpc_t, unsigned long int, unsigned long int, mpc_rnd_t)
    int  mpc_set_si_si (mpc_t, long int, long int, mpc_rnd_t)
    int  mpc_set_d_d (mpc_t, double, double, mpc_rnd_t)
    int  mpc_set_fr_fr (mpc_t, mpfr_t, mpfr_t, mpc_rnd_t)
    void mpc_set_nan(mpc_t)
    void mpc_swap(mpc_t, mpc_t)

    # Comparisons
    int  mpc_cmp (mpc_t, mpc_t)
    int  mpc_cmp_si_si (mpc_t, long int, long int)
    int  mpc_cmp_si (mpc_t, long int)

    # Projection
    int mpc_real (mpfr_t rop, mpc_t op, mpfr_rnd_t rnd)
    int mpc_imag (mpfr_t rop, mpc_t op, mpfr_rnd_t rnd)
    mpfr_t mpc_realref (mpc_t op)
    mpfr_t mpc_imagref (mpc_t op)
    int mpc_arg (mpfr_t rop, mpc_t op, mpfr_rnd_t rnd)
    int mpc_proj (mpc_t rop, mpc_t op, mpc_rnd_t rnd)

    # Arithmetic
    int  mpc_neg (mpc_t, mpc_t, mpc_rnd_t)
    int  mpc_conj (mpc_t, mpc_t, mpc_rnd_t)

    int  mpc_add (mpc_t, mpc_t, mpc_t, mpc_rnd_t)
    int  mpc_sub (mpc_t, mpc_t, mpc_t, mpc_rnd_t)
    int  mpc_mul (mpc_t, mpc_t, mpc_t, mpc_rnd_t)
    int  mpc_div (mpc_t, mpc_t, mpc_t, mpc_rnd_t)
    int  mpc_sqr (mpc_t, mpc_t, mpc_rnd_t)

    int  mpc_add_fr (mpc_t, mpc_t, mpfr_t, mpc_rnd_t)
    int  mpc_add_ui (mpc_t, mpc_t, unsigned long int, mpc_rnd_t)
    int  mpc_add_si (mpc_t, mpc_t, unsigned long, mpc_rnd_t)
    int  mpc_sub_fr (mpc_t, mpc_t, mpfr_t, mpc_rnd_t)
    int  mpc_fr_sub (mpc_t, mpfr_t, mpc_t, mpc_rnd_t)
    int  mpc_sub_ui (mpc_t, mpc_t, unsigned long int, mpc_rnd_t)
    int  mpc_ui_sub (mpc_t, unsigned long int, mpc_t, mpc_rnd_t)
    int  mpc_ui_ui_sub (mpc_t, unsigned long int, unsigned long int, mpc_t, mpc_rnd_t)
    int  mpc_mul_fr (mpc_t, mpc_t, mpfr_t, mpc_rnd_t)
    int  mpc_mul_ui (mpc_t, mpc_t, unsigned long int, mpc_rnd_t)
    int  mpc_mul_si (mpc_t, mpc_t, long int, mpc_rnd_t)
    int  mpc_mul_i  (mpc_t, mpc_t, int, mpc_rnd_t)
    int  mpc_ui_div (mpc_t, unsigned long int, mpc_t, mpc_rnd_t)
    int  mpc_div_fr (mpc_t, mpc_t, mpfr_t, mpc_rnd_t)
    int  mpc_fr_div (mpc_t, mpfr_t, mpc_t, mpc_rnd_t)
    int  mpc_div_ui (mpc_t, mpc_t, unsigned long int, mpc_rnd_t)
    int  mpc_ui_div (mpc_t, unsigned long int, mpc_t, mpc_rnd_t)
    int  mpc_div_2ui (mpc_t, mpc_t, unsigned long int, mpc_rnd_t)
    int  mpc_div_2si (mpc_t, mpc_t, long int, mpc_rnd_t)
    int  mpc_mul_2ui (mpc_t, mpc_t, unsigned long int, mpc_rnd_t)
    int  mpc_mul_2si (mpc_t, mpc_t, long int, mpc_rnd_t)
    #
    int  mpc_abs (mpfr_t, mpc_t, mp_rnd_t)
    int  mpc_norm (mpfr_t, mpc_t, mp_rnd_t)


    # Power functions and logarithm
    int  mpc_sqrt (mpc_t, mpc_t, mpc_rnd_t)
    int mpc_exp (mpc_t, mpc_t, mpc_rnd_t)
    int mpc_log (mpc_t, mpc_t, mpc_rnd_t)
    int mpc_pow (mpc_t rop, mpc_t op1, mpc_t op2, mpc_rnd_t rnd)
    int mpc_pow_si (mpc_t rop, mpc_t op1, long op2, mpc_rnd_t rnd)
    int mpc_pow_ui (mpc_t rop, mpc_t op1, unsigned long op2, mpc_rnd_t rnd)
    int mpc_pow_z (mpc_t rop, mpc_t op1, mpz_t, mpc_rnd_t rnd)
    int mpc_pow_d (mpc_t rop, mpc_t op1, double, mpc_rnd_t rnd)
    int mpc_pow_fr (mpc_t rop, mpc_t op1, mpfr_t op2, mpc_rnd_t rnd)

    # Trigonometric functions
    void mpc_sin (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_cos (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_tan (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_sinh (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_cosh (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_tanh (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_asin (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_acos (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_atan (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_asinh (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_acosh (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_atanh (mpc_t, mpc_t, mpc_rnd_t)

    # Random Function
    int  mpc_urandom (mpc_t, gmp_randstate_t)

    # utility
    int MPC_INEX_RE(int)
    int MPC_INEX_IM(int)
   
    #int  mpc_set_fr (mpc_ptr, mpfr_srcptr, mpc_rnd_t)
    #int  mpc_set_si_si (mpc_ptr, long int, long int, mpc_rnd_t)
    #int  mpc_mul_fr (mpc_ptr, mpc_srcptr, mpfr_srcptr, mpc_rnd_t)
    #int  mpc_add (mpc_ptr, mpc_srcptr, mpc_srcptr, mpc_rnd_t)
    #void mpc_init2 (mpc_ptr, mpfr_prec_t)
    #void mpc_clear (mpc_ptr)
