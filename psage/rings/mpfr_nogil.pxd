
from sage.libs.gmp.all cimport mp_exp_t, mp_size_t, gmp_randstate_t, mpz_t, mpq_t

## First stuff from sage/libs/mpfr.pxd
cdef extern from "mpfr.h" nogil:
    ctypedef struct __mpfr_struct:
        pass
        
    #ctypedef __mpfr_struct mpfr_t[1]
    ctypedef __mpfr_struct* mpfr_t
    ctypedef mpfr_t mpfr_ptr 
    ctypedef mpfr_t mpfr_srcptr
    ctypedef enum mpfr_rnd_t:  ## Somehow this does not make mpfr_rnd_t into a valid type
        MPFR_RNDN #= 0
        MPFR_RNDZ #= 1
        MPFR_RNDU #= 2
        MPFR_RNDD #= 3
        MPFR_RND_MAX = 4
        MPFR_RNDNA   = -1
    
    ctypedef mpfr_rnd_t mp_rnd_t
    ctypedef long mpfr_prec_t
#    cdef mp_rnd_t  GMP_RNDN = MP_RNDN, GMP_RNDZ = 1, GMP_RNDU = 2, GMP_RNDD = 3, GMP_RND_MAX = 4, GMP_RNDNA = -1
    int MPFR_PREC_MIN, MPFR_PREC_MAX

    # Initialization Functions
    void mpfr_init2 (mpfr_t x, mpfr_prec_t prec) 
    void mpfr_clear (mpfr_t x) 
    void mpfr_init (mpfr_t x) 
    void mpfr_set_default_prec (mpfr_prec_t prec) 
    mpfr_prec_t mpfr_get_default_prec () 
    void mpfr_set_prec (mpfr_t x, mpfr_prec_t prec) 
    mpfr_prec_t mpfr_get_prec (mpfr_t x) 

    # Assignment Functions
    #int mpfr_set (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd)
    int mpfr_set (mpfr_ptr, mpfr_srcptr, mpfr_rnd_t) 
    int mpfr_set_ui (mpfr_t rop, unsigned long int op, mpfr_rnd_t rnd) 
    int mpfr_set_si (mpfr_t rop, long int op, mpfr_rnd_t rnd) 
    # int mpfr_set_uj (mpfr_t rop, uintmax_t op, mpfr_rnd_t rnd) 
    # int mpfr_set_sj (mpfr_t rop, intmax_t op, mpfr_rnd_t rnd) 
    int mpfr_set_d (mpfr_t rop, double op, mpfr_rnd_t rnd) 
    int mpfr_set_ld (mpfr_t rop, long double op, mpfr_rnd_t rnd) 
    # int mpfr_set_decimal64 (mpfr_t rop, _Decimal64 op, mpfr_rnd_t rnd) 
    int mpfr_set_z (mpfr_t rop, mpz_t op, mpfr_rnd_t rnd) 
    int mpfr_set_q (mpfr_t rop, mpq_t op, mpfr_rnd_t rnd) 
    # int mpfr_set_f (mpfr_t rop, mpf_t op, mpfr_rnd_t rnd) 
    int mpfr_set_ui_2exp (mpfr_t rop, unsigned long int op, mp_exp_t e, mpfr_rnd_t rnd) 
    int mpfr_set_si_2exp (mpfr_t rop, long int op, mp_exp_t e, mpfr_rnd_t rnd) 
    # int mpfr_set_uj_2exp (mpfr_t rop, uintmax_t op, intmax_t e, mpfr_rnd_t rnd) 
    # int mpfr_set_sj_2exp (mpfr_t rop, intmax_t op, intmax_t e, mpfr_rnd_t rnd) 
    int mpfr_set_str (mpfr_t rop,  char *s, int base, mpfr_rnd_t rnd) 
    int mpfr_strtofr (mpfr_t rop, char *nptr, char **endptr, int base, mpfr_rnd_t rnd) 
    void mpfr_set_inf (mpfr_t x, int sign) 
    void mpfr_set_nan (mpfr_t x) 
    void mpfr_swap (mpfr_t x, mpfr_t y) 

    # Combined Initialization and Assignment Functions
    int mpfr_init_set (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_init_set_ui (mpfr_t rop, unsigned long int op, mpfr_rnd_t rnd) 
    int mpfr_init_set_si (mpfr_t rop, signed long int op, mpfr_rnd_t rnd) 
    int mpfr_init_set_d (mpfr_t rop, double op, mpfr_rnd_t rnd) 
    int mpfr_init_set_ld (mpfr_t rop, long double op, mpfr_rnd_t rnd) 
    int mpfr_init_set_z (mpfr_t rop, mpz_t op, mpfr_rnd_t rnd) 
    int mpfr_init_set_q (mpfr_t rop, mpq_t op, mpfr_rnd_t rnd) 
    # int mpfr_init_set_f (mpfr_t rop, mpf_t op, mpfr_rnd_t rnd) 
    int mpfr_init_set_str (mpfr_t x, char *s, int base, mpfr_rnd_t rnd) 

    # Conversion Functions
    double mpfr_get_d (mpfr_t op, mpfr_rnd_t rnd) 
    long double mpfr_get_ld (mpfr_t op, mpfr_rnd_t rnd) 
    # _Decimal64 mpfr_get_decimal64 (mpfr_t op, mpfr_rnd_t rnd) 
    double mpfr_get_d_2exp (long *exp, mpfr_t op, mpfr_rnd_t rnd) 
    long double mpfr_get_ld_2exp (long *exp, mpfr_t op, mpfr_rnd_t rnd) 
    long mpfr_get_si (mpfr_t op, mpfr_rnd_t rnd) 
    unsigned long mpfr_get_ui (mpfr_t op, mpfr_rnd_t rnd) 
    # intmax_t mpfr_get_sj (mpfr_t op, mpfr_rnd_t rnd) 
    # uintmax_t mpfr_get_uj (mpfr_t op, mpfr_rnd_t rnd) 
    mp_exp_t mpfr_get_z_exp (mpz_t rop, mpfr_t op) 
    void mpfr_get_z (mpz_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    # int mpfr_get_f (mpf_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    char * mpfr_get_str (char *str, mp_exp_t *expptr, int b, size_t n, mpfr_t op, mpfr_rnd_t rnd) 
    void mpfr_free_str (char *str) 
    bint mpfr_fits_ulong_p (mpfr_t op, mpfr_rnd_t rnd) 
    bint mpfr_fits_slong_p (mpfr_t op, mpfr_rnd_t rnd) 
    bint mpfr_fits_uint_p (mpfr_t op, mpfr_rnd_t rnd) 
    bint mpfr_fits_sint_p (mpfr_t op, mpfr_rnd_t rnd) 
    bint mpfr_fits_ushort_p (mpfr_t op, mpfr_rnd_t rnd) 
    bint mpfr_fits_sshort_p (mpfr_t op, mpfr_rnd_t rnd) 
    bint mpfr_fits_intmax_p (mpfr_t op, mpfr_rnd_t rnd) 
    bint mpfr_fits_uintmax_p (mpfr_t op, mpfr_rnd_t rnd) 

    # Basic Arithmetic Functions
    int mpfr_add (mpfr_t rop, mpfr_t op1, mpfr_t op2, mpfr_rnd_t rnd) 
    int mpfr_add_ui (mpfr_t rop, mpfr_t op1, unsigned long int op2, mpfr_rnd_t rnd) 
    int mpfr_add_si (mpfr_t rop, mpfr_t op1, long int op2, mpfr_rnd_t rnd) 
    int mpfr_add_z (mpfr_t rop, mpfr_t op1, mpz_t op2, mpfr_rnd_t rnd) 
    int mpfr_add_q (mpfr_t rop, mpfr_t op1, mpq_t op2, mpfr_rnd_t rnd) 
    int mpfr_sub (mpfr_t rop, mpfr_t op1, mpfr_t op2, mpfr_rnd_t rnd) 
    int mpfr_ui_sub (mpfr_t rop, unsigned long int op1, mpfr_t op2, mpfr_rnd_t rnd) 
    int mpfr_sub_ui (mpfr_t rop, mpfr_t op1, unsigned long int op2, mpfr_rnd_t rnd) 
    int mpfr_si_sub (mpfr_t rop, long int op1, mpfr_t op2, mpfr_rnd_t rnd) 
    int mpfr_sub_si (mpfr_t rop, mpfr_t op1, long int op2, mpfr_rnd_t rnd) 
    int mpfr_sub_z (mpfr_t rop, mpfr_t op1, mpz_t op2, mpfr_rnd_t rnd) 
    int mpfr_sub_q (mpfr_t rop, mpfr_t op1, mpq_t op2, mpfr_rnd_t rnd) 
    int mpfr_mul (mpfr_t rop, mpfr_t op1, mpfr_t op2, mpfr_rnd_t rnd) 
    int mpfr_mul_ui (mpfr_t rop, mpfr_t op1, unsigned long int op2, mpfr_rnd_t rnd) 
    int mpfr_mul_si (mpfr_t rop, mpfr_t op1, long int op2, mpfr_rnd_t rnd) 
    int mpfr_mul_z (mpfr_t rop, mpfr_t op1, mpz_t op2, mpfr_rnd_t rnd) 
    int mpfr_mul_q (mpfr_t rop, mpfr_t op1, mpq_t op2, mpfr_rnd_t rnd)
    int mpfr_mul_d (mpfr_ptr, mpfr_srcptr,double, mpfr_rnd_t)   ## This was not in original
    int mpfr_sqr (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_div (mpfr_t rop, mpfr_t op1, mpfr_t op2, mpfr_rnd_t rnd) 
    int mpfr_ui_div (mpfr_t rop, unsigned long int op1, mpfr_t op2, mpfr_rnd_t rnd) 
    int mpfr_div_ui (mpfr_t rop, mpfr_t op1, unsigned long int op2, mpfr_rnd_t rnd) 
    int mpfr_si_div (mpfr_t rop, long int op1, mpfr_t op2, mpfr_rnd_t rnd) 
    int mpfr_div_si (mpfr_t rop, mpfr_t op1, long int op2, mpfr_rnd_t rnd) 
    int mpfr_div_z (mpfr_t rop, mpfr_t op1, mpz_t op2, mpfr_rnd_t rnd) 
    int mpfr_div_q (mpfr_t rop, mpfr_t op1, mpq_t op2, mpfr_rnd_t rnd) 
    int mpfr_sqrt (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_sqrt_ui (mpfr_t rop, unsigned long int op, mpfr_rnd_t rnd) 
    int mpfr_cbrt (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_root (mpfr_t rop, mpfr_t op, unsigned long int k, mpfr_rnd_t rnd) 
    int mpfr_pow (mpfr_t rop, mpfr_t op1, mpfr_t op2, mpfr_rnd_t rnd) 
    int mpfr_pow_ui (mpfr_t rop, mpfr_t op1, unsigned long int op2, mpfr_rnd_t rnd) 
    int mpfr_pow_si (mpfr_t rop, mpfr_t op1, long int op2, mpfr_rnd_t rnd) 
    int mpfr_pow_z (mpfr_t rop, mpfr_t op1, mpz_t op2, mpfr_rnd_t rnd) 
    int mpfr_ui_pow_ui (mpfr_t rop, unsigned long int op1, unsigned long int op2, mpfr_rnd_t rnd) 
    int mpfr_ui_pow (mpfr_t rop, unsigned long int op1, mpfr_t op2, mpfr_rnd_t rnd) 
    int mpfr_neg (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_abs (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_dim (mpfr_t rop, mpfr_t op1, mpfr_t op2, mpfr_rnd_t rnd) 
    int mpfr_mul_2ui (mpfr_t rop, mpfr_t op1, unsigned long int op2, mpfr_rnd_t rnd) 
    int mpfr_mul_2si (mpfr_t rop, mpfr_t op1, long int op2, mpfr_rnd_t rnd) 
    int mpfr_div_2ui (mpfr_t rop, mpfr_t op1, unsigned long int op2, mpfr_rnd_t rnd) 
    int mpfr_div_2si (mpfr_t rop, mpfr_t op1, long int op2, mpfr_rnd_t rnd) 

    # Comparison Functions
    int mpfr_cmp (mpfr_t op1, mpfr_t op2) 
    int mpfr_cmp_ui (mpfr_t op1, unsigned long int op2) 
    int mpfr_cmp_si (mpfr_t op1, signed long int op2) 
    int mpfr_cmp_d (mpfr_t op1, double op2) 
    int mpfr_cmp_ld (mpfr_t op1, long double op2) 
    int mpfr_cmp_z (mpfr_t op1, mpz_t op2) 
    int mpfr_cmp_q (mpfr_t op1, mpq_t op2) 
    # int mpfr_cmp_f (mpfr_t op1, mpf_t op2) 
    int mpfr_cmp_ui_2exp (mpfr_t op1, unsigned long int op2, mp_exp_t e) 
    int mpfr_cmp_si_2exp (mpfr_t op1, long int op2, mp_exp_t e) 
    int mpfr_cmpabs (mpfr_t op1, mpfr_t op2) 
    bint mpfr_nan_p (mpfr_t op) 
    bint mpfr_inf_p (mpfr_t op) 
    bint mpfr_number_p (mpfr_t op) 
    bint mpfr_zero_p (mpfr_t op) 
    int mpfr_sgn (mpfr_t op) 
    bint mpfr_greater_p (mpfr_t op1, mpfr_t op2) 
    bint mpfr_greaterequal_p (mpfr_t op1, mpfr_t op2) 
    bint mpfr_less_p (mpfr_t op1, mpfr_t op2) 
    bint mpfr_lessequal_p (mpfr_t op1, mpfr_t op2) 
    bint mpfr_lessgreater_p (mpfr_t op1, mpfr_t op2) 
    bint mpfr_equal_p (mpfr_t op1, mpfr_t op2) 
    bint mpfr_unordered_p (mpfr_t op1, mpfr_t op2) 

    # Special Functions
    int mpfr_log (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_log2 (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_log10 (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_exp (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_exp2 (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_exp10 (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_cos (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_sin (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_tan (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_sec (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_csc (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_cot (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_sin_cos (mpfr_t sop, mpfr_t cop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_acos (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_asin (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_atan (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_atan2 (mpfr_t rop, mpfr_t y, mpfr_t x, mpfr_rnd_t rnd) 
    int mpfr_cosh (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_sinh (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_tanh (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_sech (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_csch (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_coth (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_acosh (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_asinh (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_atanh (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_fac_ui (mpfr_t rop, unsigned long int op, mpfr_rnd_t rnd) 
    int mpfr_log1p (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_expm1 (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_eint (mpfr_t y, mpfr_t x, mpfr_rnd_t rnd) 
    int mpfr_gamma (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_lngamma (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_lgamma (mpfr_t rop, int *signp, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_zeta (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_zeta_ui (mpfr_t rop, unsigned long op, mpfr_rnd_t rnd) 
    int mpfr_erf (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_erfc (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_j0 (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_j1 (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_jn (mpfr_t rop, long n, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_y0 (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_y1 (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_yn (mpfr_t rop, long n, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_fma (mpfr_t rop, mpfr_t op1, mpfr_t op2, mpfr_t op3, mpfr_rnd_t rnd) 
    int mpfr_fms (mpfr_t rop, mpfr_t op1, mpfr_t op2, mpfr_t op3, mpfr_rnd_t rnd) 
    int mpfr_agm (mpfr_t rop, mpfr_t op1, mpfr_t op2, mpfr_rnd_t rnd) 
    int mpfr_hypot (mpfr_t rop, mpfr_t x, mpfr_t y, mpfr_rnd_t rnd) 
    int mpfr_const_log2 (mpfr_t rop, mpfr_rnd_t rnd) 
    int mpfr_const_pi (mpfr_t rop, mpfr_rnd_t rnd) 
    int mpfr_const_euler (mpfr_t rop, mpfr_rnd_t rnd) 
    int mpfr_const_catalan (mpfr_t rop, mpfr_rnd_t rnd) 
    void mpfr_free_cache () 
    int mpfr_sum (mpfr_t rop, mpfr_ptr tab[], unsigned long n, mpfr_rnd_t rnd) 

    # Input and Output Functions
    # size_t mpfr_out_str (file *stream, int base, size_t n, mpfr_t op, mpfr_rnd_t rnd) 
    # size_t mpfr_inp_str (mpfr_t rop, file *stream, int base, mpfr_rnd_t rnd) 

    # Integer Related Functions
    int mpfr_rint (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_ceil (mpfr_t rop, mpfr_t op) 
    int mpfr_floor (mpfr_t rop, mpfr_t op) 
    int mpfr_round (mpfr_t rop, mpfr_t op) 
    int mpfr_trunc (mpfr_t rop, mpfr_t op) 
    int mpfr_rint_ceil (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_rint_floor (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_rint_round (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_rint_trunc (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_frac (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd) 
    int mpfr_remainder (mpfr_t r, mpfr_t x, mpfr_t y, mpfr_rnd_t rnd) 
    int mpfr_remquo (mpfr_t r, long* q, mpfr_t x, mpfr_t y, mpfr_rnd_t rnd) 
    bint mpfr_integer_p (mpfr_t op) 

    # Miscellaneous Functions
    void mpfr_nexttoward (mpfr_t x, mpfr_t y) 
    void mpfr_nextabove (mpfr_t x) 
    void mpfr_nextbelow (mpfr_t x) 
    int mpfr_min (mpfr_t rop, mpfr_t op1, mpfr_t op2, mpfr_rnd_t rnd) 
    int mpfr_max (mpfr_t rop, mpfr_t op1, mpfr_t op2, mpfr_rnd_t rnd) 
    int mpfr_urandomb (mpfr_t rop, gmp_randstate_t state) 
    void mpfr_random (mpfr_t rop) 
    void mpfr_random2 (mpfr_t rop, mp_size_t size, mp_exp_t exp) 
    mp_exp_t mpfr_get_exp (mpfr_t x) 
    int mpfr_set_exp (mpfr_t x, mp_exp_t e) 
    int mpfr_signbit (mpfr_t op) 
    int mpfr_setsign (mpfr_t rop, mpfr_t op, int s, mpfr_rnd_t rnd) 
    int mpfr_copysign (mpfr_t rop, mpfr_t op1, mpfr_t op2, mpfr_rnd_t rnd) 
    char * mpfr_get_version () 
    long MPFR_VERSION_NUM (major, minor, patchlevel) 
    char * mpfr_get_patches () 

    # Rounding Mode Related Functions
    void mpfr_set_default_rounding_mode (mpfr_rnd_t rnd) 
    mpfr_rnd_t mpfr_get_default_rounding_mode () 
    int mpfr_prec_round (mpfr_t x, mpfr_prec_t prec, mpfr_rnd_t rnd) 
    int mpfr_round_prec (mpfr_t x, mpfr_rnd_t rnd, mpfr_prec_t prec) 
    char * mpfr_print_rnd_mode (mpfr_rnd_t rnd) 

    # Exception Related Functions
    mp_exp_t mpfr_get_emin () 
    mp_exp_t mpfr_get_emax () 
    int mpfr_set_emin (mp_exp_t exp) 
    int mpfr_set_emax (mp_exp_t exp) 
    mp_exp_t mpfr_get_emin_min () 
    mp_exp_t mpfr_get_emin_max () 
    mp_exp_t mpfr_get_emax_min () 
    mp_exp_t mpfr_get_emax_max () 
    int mpfr_check_range (mpfr_t x, int t, mpfr_rnd_t rnd) 
    int mpfr_subnormalize (mpfr_t x, int t, mpfr_rnd_t rnd) 
    void mpfr_clear_underflow () 
    void mpfr_clear_overflow () 
    void mpfr_clear_nanflag () 
    void mpfr_clear_inexflag () 
    void mpfr_clear_erangeflag () 
    void mpfr_set_underflow () 
    void mpfr_set_overflow () 
    void mpfr_set_nanflag () 
    void mpfr_set_inexflag () 
    void mpfr_set_erangeflag () 
    void mpfr_clear_flags () 
    bint mpfr_underflow_p () 
    bint mpfr_overflow_p () 
    bint mpfr_nanflag_p () 
    bint mpfr_inexflag_p () 
    bint mpfr_erangeflag_p () 
    bint mpfr_divby0_p()

    
    # Advanced Functions
    MPFR_DECL_INIT (name, prec) 
    void mpfr_inits (mpfr_t x, ...) 
    void mpfr_inits2 (mpfr_prec_t prec, mpfr_t x, ...) 
    void mpfr_clears (mpfr_t x, ...) 

    # Compatibility With MPF
    void mpfr_set_prec_raw (mpfr_t x, mpfr_prec_t prec) 
    int mpfr_eq (mpfr_t op1, mpfr_t op2, unsigned long int op3) 
    void mpfr_reldiff (mpfr_t rop, mpfr_t op1, mpfr_t op2, mpfr_rnd_t rnd) 
    int mpfr_mul_2exp (mpfr_t rop, mpfr_t op1, unsigned long int op2, mpfr_rnd_t rnd) 
    int mpfr_div_2exp (mpfr_t rop, mpfr_t op1, unsigned long int op2, mpfr_rnd_t rnd) 

    # Custom Interface
    size_t mpfr_custom_get_size (mpfr_prec_t prec) 
    void mpfr_custom_init (void *significand, mpfr_prec_t prec) 
    void mpfr_custom_init_set (mpfr_t x, int kind, mp_exp_t exp, mpfr_prec_t prec, void *significand) 
    int mpfr_custom_get_kind (mpfr_t x) 
    void * mpfr_custom_get_mantissa (mpfr_t x) 
    mp_exp_t mpfr_custom_get_exp (mpfr_t x) 
    void mpfr_custom_move (mpfr_t x, void *new_position) 
    # Internals
    int mpfr_can_round (mpfr_t b, mp_exp_t err, mpfr_rnd_t rnd1, mpfr_rnd_t rnd2, mpfr_prec_t prec) 
    double mpfr_get_d1 (mpfr_t op) 

#     ctypedef long mpfr_prec_t
#     #ctypedef unsigned int mpfr_uprec_t
#     #ctypedef int mpfr_sign_t
#     #ctypedef short mpfr_exp_t
#     #ctypedef unsigned short mpfr_uexp
    
#     ctypedef struct __mpfr_struct:
#         pass

#     ctypedef __mpfr_struct mpfr_t[1]
#     ctypedef __mpfr_struct *mpfr_ptr
    
#     ctypedef __mpfr_struct *mpfr_srcptr

#     #        mpfr_prec_t  _mp_prec
#     #        mpfr_sign_t  _mpfr_sign
#     #        mpfr_exp_t   _mpfr_exp
#     #        mp_limb_t   *_mpfr_d

#     ctypedef enum mpfr_rnd_t:
# #        pass
#         GMP_RNDN = 0
#         GMP_RNDZ = 1
#         GMP_RNDU = 2
#         GMP_RNDD = 3
#         GMP_RND_MAX = 4
#         GMP_RNDNA = -1
# #        MPFR_RNDN=0,
# #        MPFR_RNDZ=1,
# #        MPFR_RNDU=2,
# #        MPFR_RNDD=3,
# #        MPFR_RND_MAX=4,
# #        MPFR_RNDA,
# #        MPFR_RNDF,
# #        MPFR_RNDNA=-1
        
#     ctypedef mpfr_rnd_t mpfr_rnd_t
#     void mpfr_init2 (mpfr_ptr, mpfr_prec_t)
#     void mpfr_init (mpfr_ptr)
#     void mpfr_clear (mpfr_ptr)
#     int mpfr_const_pi (mpfr_t, mpfr_rnd_t)
#     int mpfr_mul_ui (mpfr_t, mpfr_t, unsigned long int, mpfr_rnd_t)
#     int mpfr_sqrt (mpfr_ptr, mpfr_srcptr, mpfr_rnd_t)
#     int mpfr_add (mpfr_ptr, mpfr_srcptr,mpfr_srcptr, mpfr_rnd_t)
#     int mpfr_get_prec(mpfr_srcptr)
#     int mpfr_add_ui (mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t)
#     void mpfr_set4 (mpfr_ptr, mpfr_srcptr, mpfr_rnd_t, int)
#     void mpfr_set (mpfr_ptr, mpfr_srcptr, mpfr_rnd_t)
#     unsigned long mpfr_get_ui (mpfr_srcptr, mpfr_rnd_t)
#     int mpfr_set_ui (mpfr_ptr, unsigned long, mpfr_rnd_t)
#     int mpfr_cos (mpfr_ptr, mpfr_srcptr,mpfr_rnd_t)
#     int mpfr_sin (mpfr_ptr, mpfr_srcptr,mpfr_rnd_t)
#     int mpfr_exp (mpfr_ptr, mpfr_srcptr,mpfr_rnd_t)
#     double mpfr_get_d (mpfr_srcptr, mpfr_rnd_t)
#     #int mpfr_const_pi (mpfr_ptr, mpfr_rnd_t)
#     int mpfr_pow (mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t)
#     int mpfr_pow_si (mpfr_ptr, mpfr_srcptr, long int, mpfr_rnd_t)
#     int mpfr_pow_ui (mpfr_ptr, mpfr_srcptr, unsigned long int, mpfr_rnd_t)
#     int mpfr_log (mpfr_ptr, mpfr_srcptr,mpfr_rnd_t)
#     int mpfr_set_si (mpfr_ptr, long, mpfr_rnd_t)
#     int mpfr_cmp  (mpfr_srcptr, mpfr_srcptr)
#     int mpfr_cmp3 (mpfr_srcptr, mpfr_srcptr, int)
#     int mpfr_cmp_d (mpfr_srcptr, double)
#     int mpfr_cmp_ld (mpfr_srcptr, long double)
#     int mpfr_cmpabs (mpfr_srcptr, mpfr_srcptr)
#     int mpfr_cmp_ui (mpfr_srcptr, unsigned long)
#     int mpfr_cmp_si (mpfr_srcptr, long)
#     int mpfr_mul (mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t)
#     int mpfr_abs (mpfr_ptr, mpfr_srcptr, mpfr_rnd_t)
#     int mpfr_neg (mpfr_ptr, mpfr_srcptr, mpfr_rnd_t)
#     int mpfr_div_ui (mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t)
#     int mpfr_set_si_2exp (mpfr_ptr, long, mpfr_exp_t, mpfr_rnd_t)
#     void mpfr_set_prec (mpfr_ptr, mpfr_prec_t)
#     int mpfr_div (mpfr_ptr, mpfr_srcptr,mpfr_srcptr, mpfr_rnd_t)
#     int mpfr_sqrt (mpfr_ptr, mpfr_srcptr,mpfr_rnd_t)
#     int mpfr_erfc (mpfr_ptr, mpfr_srcptr,mpfr_rnd_t)
#     int mpfr_const_euler (mpfr_ptr, mpfr_rnd_t)
#     int mpfr_set_d(mpfr_ptr, double, mpfr_rnd_t)
#     int mpfr_sub_si (mpfr_ptr, mpfr_srcptr, long int, mpfr_rnd_t)
#     int mpfr_mul_si (mpfr_ptr, mpfr_srcptr, long int, mpfr_rnd_t)
#     int mpfr_zero_p (mpfr_srcptr)
#     int mpfr_div_si (mpfr_ptr, mpfr_srcptr, long int, mpfr_rnd_t)
#     int mpfr_mul_d (mpfr_ptr, mpfr_srcptr,double, mpfr_rnd_t)
#     float mpfr_get_flt (mpfr_srcptr, mpfr_rnd_t)
#     #int mpfr_zero_p (mpfr_srcptr)
#     int mpfr_sub (mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t)
#     long mpfr_get_si (mpfr_srcptr, mpfr_rnd_t)


## and then mpc.pxi
cdef extern from "mpc.h" nogil:
    ctypedef void* mpc_t
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
    int  mpc_abs (mpfr_t, mpc_t, mpfr_rnd_t)
    int  mpc_norm (mpfr_t, mpc_t, mpfr_rnd_t)


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

    


 ### MPZ
cdef extern from "gmp.h" nogil:
    void mpz_bin_uiui (mpz_t rop, unsigned long int n, unsigned long int k) 


cdef extern from "stdio.h" nogil:
    cdef extern void printf(char *fmt,...) nogil
    
