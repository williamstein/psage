r"""
Extra functions for mpc numbers.

Note: When the builtin mpc_mul etc. catch up to mine you can delete these routines...

"""

from sage.rings.complex_mpc import _mpfr_rounding_modes,_mpc_rounding_modes
from sage.rings.complex_mpc cimport MPComplexNumber
from sage.rings.real_mpfr cimport RealNumber
from sage.rings.real_mpfr import RealField



cdef inline void mpc_z_to_pol(mpc_t * res, mpc_t z,mpfr_rnd_t rnd_re):
    #assert_nsame(res[0].re, z.re)
    mpfr_hypot(res[0].re, z.re, z.im, rnd_re)
    mpfr_atan2(res[0].im, z.im, z.re, rnd_re)



## This is about 20% faster than standard mpc_mul
cdef inline void _mpc_mul(mpc_t* z, mpc_t a, mpc_t b, mpc_t *t, mpc_rnd_t rnd, mpfr_rnd_t rnd_re):
    # __inline__ void compl_mul(complex_t * z, complex_t a, complex_t b)
    #assert_nsame(z->re, a.re);
    #assert_nsame(z->re, b.re);
    if mpfr_zero_p(a.re):
        if mpfr_zero_p(a.im):
            mpc_set_si(z[0], 0, rnd)
            return
        if mpfr_zero_p(b[0].im):
            mpfr_set_si(z[0].re, 0, rnd_re)
        else:
            mpfr_mul(z[0].re, a.im, b[0].im, rnd_re)
            mpfr_neg(z[0].re, z[0].re, rnd_re)            
        if mpfr_zero_p(b[0].re):
            mpfr_set_si(z[0].im, 0, rnd_re)
        else:
            mpfr_mul(z[0].im, a.im, b[0].re, rnd_re)
        return
    if mpfr_zero_p(a.im):
        if mpfr_zero_p(b[0].re):
            mpfr_set_si(t[0].re, 0, rnd_re)
        else:
            mpfr_mul(t[0].re, a.re, b[0].re, rnd_re)
        if mpfr_zero_p(b[0].im):
            mpfr_set_si(t[0].im, 0, rnd_re)
        else:
            mpfr_mul(t[0].im, a.re, b[0].im, rnd_re)
        mpfr_set(z[0].re,t[0].re,rnd_re)
        mpfr_set(z[0].im,t[0].im,rnd_re)
        return
    if mpfr_zero_p(b[0].re):
        if mpfr_zero_p(b[0].im):
            mpc_set_si(z[0], 0, rnd)
            return
        mpfr_mul(t[0].re, a.im, b[0].im, rnd_re)
        mpfr_neg(t[0].re, t[0].re, rnd_re)
        mpfr_mul(z[0].im, a.re, b[0].im, rnd_re)
        mpfr_set(z[0].re,t[0].re,rnd_re)
        return
    if mpfr_zero_p(b[0].im):
        mpfr_mul(z[0].re, a.re, b[0].re, rnd_re)
        mpfr_mul(z[0].im, a.im, b[0].re, rnd_re)
        return
    mpfr_mul(t[0].re, a.re, b[0].re, rnd_re)
#/* z[0].re = a.re*b[0].re */
    mpfr_mul(t[0].im, a.im, b[0].im, rnd_re)
#/* z[0].im = a.im*b[0].im */
    mpfr_sub(t[0].re, t[0].re, t[0].im, rnd_re)
#/* z->re = a.re*b[0].re - a.im*b[0].im */

    mpfr_mul(t[0].im, a.im, b[0].re, rnd_re)
#/* z->im = a.im*b[0].re */
    mpfr_fma(z[0].im, a.re, b[0].im, t[0].im, rnd_re)	#/* z->im = a.re*b[0].im + a.im*b[0].re */
    mpfr_set(z[0].re,t[0].re,rnd_re)




cdef inline void _mpc_mul_fr(mpc_t* z, mpc_t a, mpfr_t b, mpc_rnd_t rnd, mpfr_rnd_t rnd_re):
    cdef int kk
    mpfr_mul(z[0].re, a.re, b, rnd_re)
    mpfr_mul(z[0].im, a.im, b, rnd_re)
    return


cdef inline void _mpc_div(mpc_t * z, mpc_t a, mpc_t b, mpc_t t[2], mpfr_rnd_t rnd_re):
    r"""

    
    
    z = a/b
    This works if z is a pointer to a but not b.
    """
    #   # t is a temporary variable
    #	assert_nsame(z->re, a.re);
    #	assert_nsame(z->re, b.re);
    # if 1
    #	m_assert(!mpfr_zero_p(b.re) || !mpfr_zero_p(b.im))
    if mpfr_zero_p(a.re):
        if mpfr_zero_p(a[0].im):
            mpfr_set_si(z[0].re, 0, rnd_re)
            mpfr_set_si(z[0].im, 0, rnd_re)
            return
        if mpfr_zero_p(b[0].re):
            mpfr_set_si(z[0].im, 0, rnd_re)
            mpfr_div(z[0].re, a[0].im, b[0].im, rnd_re)
            return
            
        if mpfr_zero_p(b[0].im):
            mpfr_set_si(z[0].re, 0, rnd_re)
            mpfr_div(z[0].im, a[0].im, b[0].re, rnd_re)
            return
        mpfr_mul(t[0].re, a[0].im, b[0].im, rnd_re)
        mpfr_mul(t[0].im, a[0].im, b[0].re, rnd_re)
        
        mpfr_mul(t[1].re, b[0].re, b[0].re, rnd_re)
        mpfr_fma(t[1].re, b[0].im, b[0].im, t[0].re, rnd_re)
        mpfr_div(z[0].re, t[0].re, t[1].re, rnd_re)
        mpfr_div(z[0].im, t[0].im, t[1].re, rnd_re)
        return
    

    if mpfr_zero_p(a[0].im):
        if mpfr_zero_p(b[0].re):
            mpfr_set_si(z[0].re, 0, rnd_re)
            mpfr_div(z[0].im, a[0].re, b[0].im, rnd_re)
            mpfr_neg(z[0].im, z[0].im, rnd_re)
            return
        if mpfr_zero_p(b[0].im):
            mpfr_set_si(z[0].im, 0, rnd_re)
            mpfr_div(z[0].re, a[0].re, b[0].re, rnd_re)
            return
        #print "here: a=",print_mpc(a)
        #print "here: b=",print_mpc(b)
        mpfr_mul(t[0].re, a[0].re, b[0].re, rnd_re)
        mpfr_mul(t[0].im, a[0].re, b[0].im, rnd_re)
        mpfr_neg(t[0].im, t[0].im, rnd_re)
        
        mpfr_mul(t[1].re, b[0].re, b[0].re, rnd_re)
        mpfr_fma(t[1].re, b[0].im, b[0].im, t[1].re, rnd_re)
        mpfr_div(z[0].re, t[0].re, t[1].re, rnd_re)
        mpfr_div(z[0].im, t[0].im, t[1].re, rnd_re)
        return
    if mpfr_zero_p(b[0].re):
        mpfr_div(t[0].re, a[0].im, b[0].im, rnd_re)
        
        mpfr_div(t[0].im, a[0].re, b[0].im, rnd_re)
        mpfr_neg(t[0].im, t[0].im, rnd_re)
        mpfr_set(z[0].re,t[0].re,rnd_re)
        mpfr_set(z[0].im,t[0].im,rnd_re)
        return
    if mpfr_zero_p(b[0].im):
        _mpc_div_fr(&z[0], a, b[0].re,rnd_re)
        return
    # endif
    mpfr_mul(t[0].im, a[0].im, b[0].re, rnd_re)
#/* t[0].im = a.im*b.re */
    mpfr_mul(t[0].re, a[0].re, b[0].im, rnd_re)
#/* t[0].re = a.re*b.im */
    mpfr_sub(t[0].im, t[0].im, t[0].re, rnd_re)
#/* z.im = a.im*b.re - a.re*b.im */
    
    mpfr_mul(t[0].re, a[0].re, b[0].re, rnd_re)
#/* t[0].re = a.re*b.re */
    mpfr_fma(t[0].re, a[0].im, b[0].im, t[0].re, rnd_re)
#/* z.im = a.re*b.re + a.im*b.im */
    mpfr_mul(t[1].re, b[0].re, b[0].re, rnd_re)
    mpfr_fma(t[1].re, b[0].im, b[0].im, t[1].re, rnd_re)
    
    assert(not mpfr_zero_p(t[0].re))
    mpfr_div(z[0].re, t[0].re, t[1].re, rnd_re)
    mpfr_div(z[0].im, t[0].im, t[1].re, rnd_re)




cdef inline void _mpc_add(mpc_t *res, mpc_t a, mpc_t b,mpfr_rnd_t rnd_re):
    mpfr_add(res[0].re, a[0].re, b[0].re, rnd_re)
    mpfr_add(res[0].im, a[0].im, b[0].im, rnd_re)


cdef inline void _mpc_sub(mpc_t *res, mpc_t a, mpc_t b,mpfr_rnd_t rnd_re):
    mpfr_sub(res[0].re, a[0].re, b[0].re, rnd_re)
    mpfr_sub(res[0].im, a[0].im, b[0].im, rnd_re)


cdef inline void _mpc_div_fr(mpc_t * r, mpc_t a, mpfr_t d, mpfr_rnd_t rnd):
    mpfr_div(r[0].re, a[0].re, d, rnd)
    mpfr_div(r[0].im, a[0].im, d, rnd)


cdef inline void _mpc_div_ui(mpc_t *res,mpc_t z, unsigned int i, mpfr_rnd_t rnd_re):
    mpfr_div_ui(res[0].re,z.re,i,rnd_re)
    mpfr_div_ui(res[0].im,z.im,i,rnd_re)



cdef inline void _mpc_set(mpc_t *res, mpc_t z, mpfr_rnd_t rnd_re):
    mpfr_set(res[0].re, z[0].re,rnd_re)
    mpfr_set(res[0].im, z[0].im,rnd_re)

cdef inline void _mpc_conj(mpc_t *res, mpc_t z, mpfr_rnd_t rnd_re):
    mpfr_set(res[0].re, z[0].re,rnd_re)
    mpfr_set(res[0].im, z[0].im,rnd_re)
    mpfr_neg(res[0].im, res[0].im,rnd_re)

# Note: mpc_abs and mpc_sqrt are more or less optimal.
#       we don't need to make any new versions of these


cdef inline void _pochammer(mpc_t *res, mpc_t z, int n,mpc_rnd_t rnd,mpfr_rnd_t rnd_re):
    r"""
    Pochammer symbol
    """
    cdef int j
    cdef int prec= mpc_get_prec(z)
    cdef mpc_t x
    cdef mpc_t t[2]
    mpc_init2(t[0],prec+10)
    mpc_init2(t[1],prec+10)
    mpc_init2(x,prec+10)
    mpc_set_ui(res[0],1,rnd)
    for j from 0 <= j < n:
        #mpc_add_ui(z.re,z.re,1,rnd_re)
        mpc_add_ui(x,z,j,rnd)
        _mpc_mul(res,res[0],x,t,rnd,rnd_re)
    mpc_clear(t[0])
    mpc_clear(t[1])
    mpc_clear(x)

def mpc_pochammer(MPComplexNumber z,n):
    assert isinstance(z,MPComplexNumber)
    cdef mpc_t x
    cdef MPComplexNumber res
    #assert ZZ(n)==n
    F = z.parent()
    mpc_init2(x,F.prec())
    res = F.one()
    rnd_re = _mpfr_rounding_modes.index(F.rounding_mode_real())
    rnd = _mpc_rounding_modes.index(F.rounding_mode())
    _pochammer(&x,z.value,n,rnd,rnd_re)
    mpc_set(res.value,x,rnd)
    mpc_clear(x)
    return res



cdef print_mpc(mpc_t z):
    from sage.rings.complex_mpc import MPComplexField
    cdef int prec
    cdef MPComplexNumber zz
    prec = mpc_get_prec(z)
    zz = MPComplexField(prec)(1)
    mpc_set(zz.value,z,MPC_RNDNN)
    return str(zz)

cdef print_mpfr(mpfr_t x):
    cdef int prec
    cdef RealNumber xx
    prec = mpfr_get_prec(x)
    xx = RealField(prec)(1)
    mpfr_set(xx.value,x,GMP_RNDN)
    return str(xx)


cdef void mpc_sign(mpc_t *res, mpc_t *z,mpfr_rnd_t rnd_re):
    r"""
    The complex sign function. Returns z/abs(z)=exp(i*Arg(z))
    """
    mpc_abs(res[0].re, z[0],rnd_re)
    if mpfr_zero_p(res[0].re):
        mpfr_set_ui(res[0].re,1,rnd_re)
        mpfr_set_ui(res[0].im,0,rnd_re)
        #assert not 
    else:
        mpfr_div(res[0].im, z[0].im, res[0].re, rnd_re)
        mpfr_div(res[0].re, z[0].re, res[0].re, rnd_re)


