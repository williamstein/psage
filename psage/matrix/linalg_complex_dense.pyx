r"""
Helper subroutines for linear algebra for dense multi-preision complex 
matrices (using mpc and mpfr).
Part of the original C-code is due to Markus Fraczek (2010).

Algorithms: 
- QR-decomposition of A
- Eigenvalues of A
- Eigenvalues and vectors of A
The main idea is to use QR-reduction computed by Givens rotations for complex matrices.

"""

###
###  efficient subroutines for complex dense linear algebra
### 


#mpfr_set_default_prec(103)

include "sage/ext/gmp.pxi"
include "sage/ext/random.pxi"
#include "sage/rings/mpc.pxi"

from sage.rings.complex_mpc cimport MPComplexNumber
from sage.rings.real_mpfr cimport RealNumber,RealField_class
from sage.rings.real_mpfr import RealField
from sage.rings.complex_mpc import MPComplexField
from sage.rings.complex_mpc import _mpfr_rounding_modes,_mpc_rounding_modes
from psage.rings.mpc_extras cimport *


cdef void givens(mpc_t *s, mpc_t *skk, mpfr_t *c, mpc_t *r, mpc_t *f, mpc_t *g, mpc_t *t, mpc_rnd_t rnd,mpfr_rnd_t rnd_re):
    r"""
    Return parameters for Givens rotations rotating the vector (f , g)
    to (r , 0). I.e. [[ c, s],[ -\bar{s}, c]] (f,g)^t = (r,0)^t 

    See: Bindel, Demmel and Kahan, 'On Computing Givens Rotations Reliably and Efficiently', ACM Transactions n Mathematical Software, Vol. 28, No. 2, 2002, p. 206-238.
    The below is Algorithm 1 of this paper.
    TODO: Implement a 'safe' (slower) version, i.e. the read the rest of the paper and implement all over/underflow checks.

    """
    if mpc_zero_p(g[0]):
        mpfr_set_si(c[0], 1, rnd_re)
        mpc_set_si(s[0], 0, rnd)
        _mpc_set(r, f[0],rnd_re)
    elif mpc_zero_p(f[0]):
        mpfr_set_si(c[0], 0, rnd_re)
        _mpc_conj(r, g[0],rnd_re)
        #mpfr_set(r[0].re, g[0].re,rnd_re)
        #mpfr_set(r[0].im, g[0].im,rnd_re)
        #mpfr_neg(r[0].im, r[0].im,rnd_re)
        mpc_sign(s, r,rnd_re)
        mpc_abs(r[0].re, g[0],rnd_re)
        mpfr_set_si(r[0].im, 0, rnd_re)
    else:
        #print "f=",print_mpc(f[0])
        #print "g=",print_mpc(g[0])
        _mpc_conj(r, g[0],rnd_re)
        #mpfr_set(r[0].re, g[0].re,rnd_re)
        #mpfr_set(r[0].im, g[0].im,rnd_re)
        #mpfr_neg(r[0].im, r[0].im,rnd_re)
        _mpc_mul(s, f[0], r[0],t,rnd,rnd_re)
        #print "f*conj(g)=",print_mpc(s[0])
        mpc_abs(r[0].re, f[0],rnd_re)
        mpc_abs(r[0].im, g[0],rnd_re)
        mpfr_hypot(c[0], r[0].re, r[0].im, rnd_re)
        #print "|f|^2+|g|^2=",print_mpfr(c[0])
        _mpc_div_fr(s, s[0], r[0].re,rnd_re)
        #print "f*conj(g)/|f|=",print_mpc(s[0])
        _mpc_div_fr(s, s[0], c[0],rnd_re)
        #print "f*conj(g)/|f|/hyp(|f|,|g|)=",print_mpc(s[0])
        mpfr_div(c[0], r[0].re, c[0], rnd_re)
        _mpc_div_fr(r, f[0], c[0],rnd_re)
    _mpc_conj(skk, s[0],rnd_re)


cdef void RotateLeft_(mpc_t s, mpc_t skk, mpfr_t c,mpc_t** A, int i, int k, int min_j,int max_j, mpc_t *t,mpc_rnd_t rnd,mpfr_rnd_t rnd_re):
    r"""
    Apply the Givens rotation [[ c, s ], [-\bar{s},c]] from the left. Acts on rows i and k from columns min_j to max_j
    """
    cdef int j
    for j from  min_j <= j < max_j:
        ## the below might seem compliated but is
        ## necessary to minimize number of refs
        _mpc_mul_fr(t, A[i][j], c,rnd, rnd_re)
        _mpc_mul(t+1, s, A[k][j],t+2,rnd,rnd_re)
        _mpc_add(t+1, t[0], t[1],rnd_re)
        _mpc_mul(t, skk, A[i][j],t+2,rnd,rnd_re)
        _mpc_set(&A[i][j], t[1],rnd_re)
        _mpc_mul_fr(t+1, A[k][j], c,rnd, rnd_re)
        _mpc_sub(&A[k][j], t[1], t[0],rnd_re)


cdef void RotateRight_(mpc_t s, mpc_t skk, mpfr_t c,mpc_t** A, int i, int k, int min_j, int max_j, mpc_t *t,mpc_rnd_t rnd, mpfr_rnd_t rnd_re):
    r"""
    Apply the inverse Givens rotation [[ c, -s ], [\bar{s},c]] from the right. Acts on columns i and k from rows min_j to max_j
    """
    cdef int j
    cdef int prec = mpc_get_prec(s)
    for j from min_j<= j < max_j:
        _mpc_mul_fr(t, A[j][i], c,rnd,rnd_re)
        _mpc_mul(t+1, skk, A[j][k],t+2,rnd,rnd_re)
        _mpc_add(t+1, t[0], t[1],rnd_re)       
        _mpc_mul(t, s, A[j][i],t+2,rnd,rnd_re)
        _mpc_set(&A[j][i], t[1],rnd_re)
        _mpc_mul_fr(t+1, A[j][k], c,rnd,rnd_re)
        _mpc_sub(&A[j][k], t[1], t[0],rnd_re)

# The QR_set struct contains various parameters useful in computing
# QR-decomposition.


cdef QR_set init_QR(prec):
    r"""
    Initialize all variables in a QR_set.
    """
    cdef QR_set res
    mpc_init2(res.s,prec) # = mpc_init()
    mpc_init2(res.skk,prec) 
    mpfr_init2(res.c,prec)
    mpc_init2(res.r,prec)
    # t's and x's are used as temporary variables
    mpc_init2(res.t[0],prec)
    mpc_init2( res.t[1],prec)
    mpc_init2(res.t[2],prec)
    mpc_init2(res.t[3],prec)
    mpc_init2(res.t[4],prec)
    mpfr_init2(res.x[0],prec)
    mpfr_init2(res.x[1],prec)
    mpfr_init2(res.x[2],prec)
    mpfr_init2(res.x[3],prec)
    mpc_init2(res.mu1,prec)
    mpc_init2(res.mu2,prec)
    mpc_init2(res.s_t,prec)
    mpc_init2(res.skk_t,prec)
    mpfr_init2(res.c_t,prec)
    return res


cdef void clear_QR(QR_set * q):
    r"""
    Clear all variables in a QR_set.
    """
    mpc_clear(q.s)
    mpc_clear(q.skk)
    mpfr_clear(q.c)
    mpc_clear(q.r)
    mpc_clear(q.t[0])
    mpc_clear(q.t[1])
    mpc_clear(q.t[2])
    mpc_clear(q.t[3])
    mpc_clear(q.t[4])
    mpfr_clear(q.x[0])
    mpfr_clear(q.x[1])
    mpfr_clear(q.x[2])
    mpfr_clear(q.x[3])    
    mpc_clear(q.mu1)
    mpc_clear(q.mu2)
    mpc_clear(q.s_t)
    mpc_clear(q.skk_t)
    mpfr_clear(q.c_t)


cdef void qr_decomp(mpc_t** A,int m, int n, int prec, mpfr_t eps, mpc_rnd_t rnd, mpfr_rnd_t rnd_re):
    r"""
    Compute QR-decomposition of an m x n matrix using complex givens rotations.
    Overwrites A with R and keeps the information about Q
    in the zero-set entries.
    """
    cdef QR_set q = init_QR(prec)
    cdef int i,j,k
    cdef mpc_t t
    cdef mpfr_t x
    cdef MPComplexNumber z
    z=MPComplexField(prec)(1)
    assert  m >= n
    mpc_init2(t,prec)
    mpfr_init2(x,prec)
    for j from 0 <= j < n:
        for i from j+1 <= i < m:
            givens(&q.s, &q.skk, &q.c, &q.r,&A[j][j], &A[i][j],q.t,rnd,rnd_re)
            _mpc_div(&t,A[i][j],A[j][j],q.t,rnd_re)
            _mpc_conj(&t,t,rnd_re)
            if j< n-1:
                RotateLeft_( q.s, q.skk, q.c, A, j, i, j+1, n, q.t,rnd,rnd_re)
            _mpc_set(&A[i][j], t,rnd_re)
            _mpc_set(&A[j][j], q.r,rnd_re)
            mpc_set(z.value,q.s,MPC_RNDNN)
            # print "q.s=",z
            # mpc_set(z.value,q.skk,MPC_RNDNN)
            # print "q.skk=",z
            # mpc_set(z.value,q.t[0],MPC_RNDNN)
            # print "q.t[0]=",z
            # mpc_set(z.value,q.t[1],MPC_RNDNN)
            # print "q.t[1]=",z
            # mpc_set(z.value,q.t[2],MPC_RNDNN)
            # print "q.t[2]=",z
            # mpc_set(z.value,q.t[3],MPC_RNDNN)
            # print "q.t[3]=",z
            # mpc_set(z.value,q.t[4],MPC_RNDNN)
            # print "q.t[4]=",z

            # mpc_set(z.value,q.r,MPC_RNDNN)
            # print "q.r=",z
            
    # Now A consists of R on the n x n upper-triangular part
    # and the A[i][j]'s are t_ij corresponding to each step of the





cdef _swap_rows(mpc_t** A, int *m, int *n):
    cdef mpc_t* tmprow = A[m[0]]
    A[n[0]] = A[m[0]]
    A[m[0]] = tmprow
    
cdef _swap_columns(mpc_t** A, int *nrows, int i, int m):
    cdef int k
    for k from 0<= k < nrows[0]:
        mpc_swap(A[k][i],A[k][m])
          

cdef void _qr_step(mpc_t** A, int *nrows, int *start, int *end, QR_set * q, mpc_t mu,mpc_rnd_t rnd,mpfr_rnd_t rnd_re):
    r"""
    Compute one QR-step:
     Q * R = A + mu * Id
     A     = R * Q +mu * Id
     
    """
    cdef int j
    assert start[0] >= 0
    assert end[0] <= nrows[0]
    for j from start[0] <= j < end[0]:
        #_mpc_sub(&A[j][j], A[j][j], mu,rnd_re)
        mpc_sub(A[j][j], A[j][j], mu,rnd)

    for j from start[0] <= j < end[0]:
        if j < end[0] - 1:
            #/*if (mpc_zero(A[j + 1][j]))   continue; */
            givens(&q.s, &q.skk, &q.c, &q.r, &A[j][j], &A[j + 1][j],q.t,rnd,rnd_re)
            mpc_set_si(A[j + 1][j], 0, rnd)
            _mpc_set(&A[j][j], q.r,rnd_re)
            RotateLeft_(q.s, q.skk, q.c, A, j, j + 1, j + 1, end[0], q.t,rnd,rnd_re)
        if j > start[0]:
            RotateRight_(q.s_t, q.skk_t, q.c_t, A, j - 1, j, start[0], j + 1, q.t,rnd,rnd_re)	
        _mpc_set(&q.s_t, q.s,rnd_re)
        _mpc_set(&q.skk_t, q.skk,rnd_re)
        mpfr_set(q.c_t, q.c, rnd_re)
    for j from start[0] <= j < end[0]:
        mpc_add(A[j][j], A[j][j], mu,rnd) #_re)

#/* l1 = a, l2 = d*/
#
cdef int Eigenvalue2(mpc_t * l1, mpc_t * l2, mpc_t a, mpc_t b, mpc_t c, mpc_t d, mpc_t t[3],mpc_rnd_t rnd,mpfr_rnd_t rnd_re):
    
    if (mpc_zero_p(b) or mpc_zero_p(c)):
        _mpc_set(l1, a,rnd_re)
        _mpc_set(l2, d,rnd_re)
        return 0
    if (mpc_zero_p(a) and mpc_zero_p(d)) :
        _mpc_mul(l2, b, c,t,rnd, rnd_re)
        mpc_sqrt(l1[0], l2[0], rnd)
        mpc_neg(l2[0], l1[0],rnd)
        return 2
    _mpc_sub(l2, a, d,rnd_re)
    #_mpc_div_ui(l2, l2[0], 2,rnd)
    mpc_div_ui(l2[0], l2[0], 2,rnd)
    # pc_mul(l1, l2[0], l2[0],t,rnd, rnd_re)
    mpc_mul(l1[0],l2[0], l2[0],rnd)
    #_mpc_mul(l2, b, c,t,rnd, rnd_re)
    mpc_mul(l2[0], b, c,rnd)
    _mpc_add(l1, l1[0], l2[0],rnd_re)
    mpc_sqrt(l2[0], l1[0], rnd)
    _mpc_add(l1, a, d,rnd_re)
    mpc_div_ui(l1[0], l1[0], 2,rnd)
    _mpc_add(&t[0], l1[0], l2[0],rnd_re)
    _mpc_sub(l2, l1[0], l2[0],rnd_re)
    _mpc_set(l1, t[0],rnd_re)
    return 3

cdef void _qr_double_step(mpc_t ** A, int *nrows, int *start, int *end, QR_set * q, mpc_rnd_t rnd, mpfr_rnd_t rnd_re):
    r"""
    Apply the QR-step wth shifts given by both eigenvalues of
    the lower right 2x2 submatrix.
    """
    cdef int n
    n = end[0] - 1
    cdef int r
    cdef mpc_t mu
    cdef int prec = mpc_get_prec(A[0][0])
    if n<=0:
        return
    r = Eigenvalue2(&q.mu1, &q.mu2, A[n - 1][n - 1], A[n - 1][n], A[n][n - 1], A[n][n], q.t,rnd,rnd_re)
    #print "mu1 =",print_mpc(q.mu1)
    #print "mu2 =",print_mpc(q.mu2)
    _qr_step(A, nrows, start, end, q, q.mu2, rnd, rnd_re)
    if (r == 3):
        _qr_step(A, nrows, start, end, q, q.mu1, rnd, rnd_re)


cdef void _qr_single_step(mpc_t ** A, int *nrows, int *start, int *end, QR_set * q, mpc_rnd_t rnd, mpfr_rnd_t rnd_re):
    r"""
    A is an nrows x nrows complex matrix.
    """
    cdef int n
    n = end[0] - 1
    cdef int r
    cdef int prec = mpc_get_prec(A[0][0])
    if n<=0:
        return
    _wilkinson_shift(&q.mu2, &A[n - 1][n - 1], &A[n - 1][n], &A[n][n - 1], &A[n][n],prec,rnd,rnd_re, q.t,q.x)
    #print "mu=",print_mpc(q.mu2)
    _qr_step(A, nrows, start, end, q, q.mu2, rnd, rnd_re)





cdef void _wilkinson_shift(mpc_t *mu,mpc_t* a,mpc_t* b,mpc_t* c,mpc_t* d,int prec,mpc_rnd_t rnd,mpfr_rnd_t rnd_re, mpc_t *zz, mpfr_t *xx):
    r"""
    Computes the eigenvalue of [[a,b],[c,d]] closest to d.
    Algorithm 2.7 in Matrix Algorithms by G. W. Stewart

    zz should  be an allocated and initialized array of 4 mpc_t's
    """
    cdef mpfr_t *s,*x,*y,*t
    cdef mpc_t *p,*q,*r,*z  #, *w
    #cdef mpc_t w[2]
    p = zz
    q=zz+1;    z=zz+2;
    x = xx; y=xx + 1; t=xx + 2; s=xx + 3
    _mpc_set(mu,d[0],rnd_re)
    mpc_abs(s[0],a[0],rnd_re)
    #mpc_init2(res.t[3],prec)
    mpc_abs(x[0],b[0],rnd_re)
    mpfr_add(s[0],s[0],x[0],rnd_re)
    #print "|a|+|b|=",nu
    mpc_abs(x[0],c[0],rnd_re)
    mpfr_add(s[0],s[0],x[0],rnd_re)
    mpc_abs(x[0],d[0],rnd_re)
    mpfr_add(s[0],s[0],x[0],rnd_re)
    if mpfr_zero_p(s[0])<>0:
        return 
    _mpc_div_fr(p,b[0],s[0],rnd_re)
    _mpc_div_fr(q,c[0],s[0],rnd_re)
    _mpc_mul(q,q[0],p[0],zz+3,rnd,rnd_re)
    #q=(b/s)*(c/s)
    mpc_abs(x[0],q[0],rnd_re)
    #mpc_imag(y[0],q[0],rnd_re)
    #print "s=",print_mpfr(s[0])
    #print "x=",print_mpfr(x[0]),mpfr_zero_p(x[0])
    if mpfr_zero_p(x[0])==0:
        #print "here 4!"
        _mpc_div_fr(z,a[0],s[0],rnd_re)
        #print "a/s=",print_mpc(z[0])
        _mpc_div_fr(z+1,d[0],s[0],rnd_re)
        #print "d/s=",print_mpc(z[1])
        _mpc_sub(p,z[0],z[1],rnd_re)
        _mpc_div_ui(p,p[0],2,rnd_re)
        #p=((a/s)-(d/s))/two
        _mpc_mul(z,p[0],p[0],zz+3,rnd,rnd_re)
        _mpc_add(zz+3,z[0],q[0],rnd_re)
        mpc_sqrt(zz[3],zz[3],rnd)  
        #zz[3]=r=(p*p+q).sqrt()
        mpfr_mul(x[0],p[0].re,zz[3].re,rnd_re)
        mpfr_mul(y[0],p[0].im,zz[3].im,rnd_re)
        mpfr_add(t[0],x[0],y[0],rnd_re)
        #if r.real()*p.real()+p.imag()*r.imag()<0:
        if mpfr_sgn(t[0])<0:
            mpc_neg(zz[3],zz[3],rnd)
        _mpc_add(z,p[0],zz[3],rnd_re)
        _mpc_div(p,q[0],z[0],zz+3,rnd_re)
        _mpc_mul_fr(z,p[0],s[0],rnd,rnd_re)
        _mpc_sub(mu,mu[0],z[0],rnd_re)
        #mu = mu - s*(q/(p+r))

cdef void set_zero(mpc_t** A, int *nrows, mpfr_t delta, mpc_t * t,mpc_rnd_t rnd, mpfr_rnd_t rnd_re):
    cdef int j
    cdef mpfr_t x
    mpfr_init(x)
    for j from  1 <= j < nrows[0]:
        if (mpc_zero_p(A[j][j - 1])):
            continue
        #/*if (mpfr_cmpabs(A[j][j - 1][0].re,delta) < 0):
        #   if(mpfr_cmpabs(A[j][j - 1][0].im,delta) < 0):
        #     mpc_set_si(&A[j][j - 1], 0, 0)
        #     continue                            
        mpc_abs(t[0].re, A[j][j],rnd_re)
        mpc_abs(t[0].im, A[j - 1][j - 1],rnd_re)
        mpfr_add(t[0].re, t[0].re, t[0].im, rnd_re)
        #		/*if(mpfr_cmp_si(t[0].re,1) < 0)
        #		   mpfr_set_si(t[0].re,1,rnd_re) */
        ### One can do a more accurate estimate
        ### but it doesn't seem to help much... 
        #mpfr_mul(t[0].re, t[0].re, delta, rnd_re)
        #if mpfr_cmp_ui(t[0].re,1)<=0:
        #    mpfr_mul(t[0].re, t[0].re, delta, rnd_re)
        #else:
        #    mpfr_set(t[0].re,delta,rnd_re)
        #mpc_abs(t[0].im,A[j-1][j],rnd_re)
        #mpc_abs(x,A[j][j-1],rnd_re)
        #mpfr_sub(x,x,t[0].im,rnd_re)
        #mpfr_abs(x,x,rnd_re)
        #if mpfr_cmp_ui(x,1)>0:
        #    mpfr_div(t[0].re,t[0].re,x,rnd_re)
        mpfr_mul(t[0].re, t[0].re, delta, rnd_re)
        mpc_abs(t[0].im, A[j][j - 1],rnd_re)
        if (mpfr_cmp(t[0].im, t[0].re) <= 0):
            mpc_set_si(A[j][j - 1], 0, rnd)


cdef void get_submatrix(int *s, int *e, mpc_t** A, int *nrows):
    cdef int j
    e[0] = 0
    s[0] = 0
    #for j from  nrows[0] - 2 >=  j > 0 by -1:
    for j in range(nrows[0] - 2,0,-1):
        if not mpc_zero_p(A[j + 1][j]):
            if not mpc_zero_p(A[j][j - 1]):
                e[0] = j + 2
                break
            j-=1
    #print "e[0]=",e[0]
    
    for j in range(e[0]-1,0,-1):
    #for j from e[0] - 1 >= j > 0 by -1:
        if mpc_zero_p(A[j][j - 1]):
            s[0] = j
            break
    



cdef void get_eigenvalues(mpc_t* res, mpc_t** A,int n,int prec, mpc_t t[2], mpc_rnd_t rnd,mpfr_rnd_t rnd_re):
    cdef int i
    #cdef mpc_t t[2]
    #mpc_init2(t[0],prec)
    #mpc_init2(t[1],prec)
    #t = mpc_init()
    #print "A[",n-1,n-1,"]=",print_mpc(A[n-1][n-1])
    #print "n=",n
    for i from 0 <= i < n:
        if (i + 1 < n):
            #print "A[",i+1,i,"]=",print_mpc(A[i+1][i])
            if (not mpc_zero_p(A[i + 1][i])):
                Eigenvalue2(&res[i],&res[i + 1],A[i][i], A[i][i + 1], A[i + 1][i],A[i + 1][i + 1], t,rnd,rnd_re)
                #print "res2[",i,"]=",print_mpc(res[i])
                #print "res2[",i+1,"]=",print_mpc(res[i+1])
                i+=1
                continue
        _mpc_set(&res[i], A[i][i],rnd_re)
        #print "A[",i,i,"]=",print_mpc(A[i][i])
        #print "res[",i,"]=",print_mpc(res[i])


cdef void _eigenvalues(mpc_t* res, mpc_t** A,int nrows,int prec,mpc_rnd_t rnd,mpfr_rnd_t rnd_re):
    r"""
    A is a square nrows x nrows Matrix 
    
    """
    #/*int i = 0 */
    cdef mpfr_t delta, delta2
    cdef int ii,j
    cdef QR_set q 
    cdef int start, end    
    cdef double pl
    # see if we have to add to the precision
    cdef int prec_add =_get_additional_precision(A, nrows, nrows, delta, prec,rnd, rnd_re)
    #prec = prec + prec_add 
    #print "adding precision",prec_add
    q  = init_QR(prec)
    mpfr_init2(delta,prec)
    mm_get_delta(&delta,rnd_re)
    #print "delta=",print_mpfr(delta)
    _hessenberg_reduction(A,nrows, q, prec, rnd, rnd_re)
    #print "Hessenberg(A)="
    #for i from 0 <= i < nrows: # we never need more than n steps
    #    for ii from 0 <= ii < nrows: # we never need more than n steps
    #        print "A[{0}][{1}]=".format(i,ii),print_mpc(A[i][ii])
    #begin_pline("computing eigenvalues")
    end = 1
    for i from 0 <= i < 1000*nrows*nrows:
        # How many steps do we really need?
        set_zero(A, &nrows, delta, &q.r,rnd,rnd_re)
        get_submatrix(&start, &end, A, &nrows)
        #_qr_double_step(A, &nrows, &start, &end, &q, rnd, rnd_re)
        _qr_single_step(A, &nrows, &start, &end, &q, rnd, rnd_re)
        #print "After qr_single step: ",i
        #for j from 0 <= j < nrows: # we never need more than n steps
        #    for ii from 0 <= ii < nrows: # we never need more than n steps
        #print "A[{0}][{1}]=".format(j,ii),print_mpc(A[j][ii])
        if end == 0:
            #print "Stopping!"
            break
        #for ii from 0 <= ii <10:
        #    QR_double_step(A, &nrows, &start, &end, &q, rnd, rnd_re)
    #print "A=",print_mat(A,nrows,nrows)
    if end<>0:
        clear_QR(&q)
        mpfr_clear(delta)
        raise ArithmeticError,"QR-algorithm did not cinverge in {0} steps!".format(i)
    get_eigenvalues(res, A,nrows,prec,q.t,rnd,rnd_re)
    clear_QR(&q)
    mpfr_clear(delta)



cdef void _hessenberg_reduction(mpc_t** A, int nrows, QR_set  q,int prec, mpc_rnd_t rnd,mpfr_rnd_t rnd_re, int rt = 0):
    r"""
    Input a square nrows x nrows complex matrix.
    At return the matrix A is in upper Hessenberg form.
    If rt=1 then the matrix transforming A to Hessenberg form can be
    reconstructed as a sequence of Givens rotations given by parameters
    't' in the part of A below the second subdiagonal.
    """
    cdef int i,j
    cdef mpc_t t
    mpc_init2(t,prec)
    #print "FIX!!!!!"
    for j from 0<= j < nrows - 2:
        for i from j + 2 <=  i < nrows:
            if (mpc_zero_p(A[i][j])):
                continue
            givens(&q.s, &q.skk, &q.c, &q.r,&A[j + 1][j], &A[i][j],q.t, rnd,rnd_re)
            #print "A[j+1][j]=",print_mpc(A[j+1][j])
            #print "A[j][j]=",print_mpc(A[j][j])
            #print "q.s=",print_mpc(q.s)
            #print "q.skk=",print_mpc(q.skk)
            #print "q.c=",print_mpfr(q.c)
            #print "q.r=",print_mpc(q.r)
            if rt:
                _mpc_div(&t,A[i][j],A[j+1][j],q.t,rnd_re)
                mpc_conj(t,t,rnd)
            mpc_set_si(A[i][j], 0, rnd)
            mpc_set(A[j + 1][j], q.r,rnd)
            RotateLeft_( q.s, q.skk, q.c, A, j + 1, i, j + 1, nrows, q.t,rnd,rnd_re)
            # we have to rotate from right too so that we keep the eigenvalues
            RotateRight_(q.s, q.skk, q.c, A, j + 1, i, 0, nrows, q.t,rnd,rnd_re)
            if rt:
                #print "Setting A[",i,"][",j,"]=",print_mpc(t)
                mpc_set(A[i][j], t,rnd)

cdef _reconstruct_matrix(mpc_t** Q, mpc_t** A, int m, int n, int k, int prec, mpc_rnd_t rnd, mpfr_rnd_t rnd_re):
    r"""
    Use parameters in A[i,j] to reconstruct the matrix Q
    as  a product of Givens rotations.
    Here n>=j>=0 and m>=i>=j+k
    """
    cdef int i,j
    cdef mpc_t t,s,sc,tt[4]
    cdef mpfr_t x,c
    mpc_init2(s,prec); mpc_init2(sc,prec)
    mpfr_init2(x,prec); mpfr_init2(c,prec)
    mpc_init2(tt[0],prec); mpc_init2(tt[1],prec)
    mpc_init2(tt[2],prec); mpc_init2(tt[3],prec)
    mpc_init2(t,prec)
    for i from 0<= i < m:
        for j from 0 <= j < n:
            if i<>j:
                mpc_set_ui(Q[i][j],0,rnd)
        mpc_set_ui(Q[i][i],1,rnd)

    for j in xrange(n-k-1,-1,-1):
        for i in xrange(m-1,j+k-1,-1):
            mpc_set(t,A[i][j],rnd)
            mpc_set_ui(A[i][j],0,rnd)
            if mpfr_inf_p(t.re) or mpfr_inf_p(t.im):
                mpfr_set_ui(c,0,rnd_re)
                mpc_set_ui(s,1,rnd)
                mpc_set_ui(sc,1,rnd)
            else:
                mpc_abs(x,t,rnd_re)
                mpfr_mul(x,x,x,rnd_re)
                mpfr_add_ui(x,x,1,rnd_re)
                mpfr_sqrt(x,x,rnd_re)
                mpfr_ui_div(c,1,x,rnd_re)
                mpc_mul_fr(s,t,c,rnd)
            mpc_conj(sc,s,rnd)
            # multiuply from left with inverse rotation...
            mpc_neg(s,s,rnd)
            mpc_neg(sc,sc,rnd)
            RotateLeft_(s, sc, c, Q, j+k-1 ,  i, j+k-1, n, tt,rnd,rnd_re)

cdef void solve_upper_triangular(mpc_t** R,mpc_t* b,int n, int prec, mpc_rnd_t rnd, mpfr_rnd_t rnd_re):
    r"""
    Solve Rx=b where R is n x n upper-triangular matrix
    and b is a vector of length n.
    At exit b is overwritten with the solution x.
    """
    cdef mpc_t s[2]
    cdef mpc_t t[2]
    cdef mpc_t *v
    mpc_init2(s[0],prec)
    mpc_init2(s[1],prec)
    mpc_init2(t[0],prec)
    mpc_init2(t[1],prec)
    #print "b[",n-1,"]=",print_mpc(b[n-1])
    #print "R[",n-1,n-1,"]=",print_mpc(R[n-1][n-1])
    _mpc_div(&s[0],b[n-1],R[n-1][n-1], t, rnd_re)
    _mpc_set(&b[n-1],s[0],rnd_re)
    #mpc_div(s[0],b[n-1],R[n-1][n-1], rnd)
    #mpc_set(b[n-1],s[0],rnd)
    #print "b/R=",print_mpc(b[n-1])
    for i in xrange(n-2,-1,-1):
        mpc_set_si(s[0],0,rnd)
        v =<mpc_t *> sage_malloc(sizeof(mpc_t) * n)
        for j from i <= j < n:
            mpc_init2(v[j],prec)
            mpc_set(v[j],R[i][j],rnd)
        for j from i+1<= j < n:
            #_mpc_mul(&R[i][j],R[i][j],b[j],t,rnd,rnd_re)
            #_mpc_add(&s[0],s[0],R[i][j],rnd_re)
            _mpc_mul(&v[j],v[j],b[j],t,rnd,rnd_re)
            _mpc_add(&s[0],s[0],v[j],rnd_re)

        _mpc_sub(&b[i],b[i],s[0],rnd_re)
        #_mpc_div(&b[i],b[i],R[i][i],s,rnd_re)
        _mpc_div(&b[i],b[i],v[i],s,rnd_re)
        for j from i <= j < n:
            mpc_clear(v[j])
        sage_free(v)
    #for j from 0 <= j <n:
    #        print "b[",j,"]",print_mpc(b[j])
    mpc_clear(s[0])
    mpc_clear(s[1])
    mpc_clear(t[0])
    mpc_clear(t[1])
    


cpdef my_div(x,y):
    cdef mpc_t xt,yt
    cdef MPComplexNumber z
    mpc_init2(xt,x.parent().prec())
    mpc_init2(yt,x.parent().prec())
    mpc_div(xt,xt,yt,MPC_RNDNN)
    z = x.parent().one()
    mpc_set(z.value,xt,MPC_RNDNN)
    return z


cdef void _norm(mpfr_t norm,mpc_t** A, int nrows,int ncols, int prec,mpfr_rnd_t rnd_re,int ntype):
    r"""
    """
    cdef mpfr_t x
    cdef int i,j
    mpfr_init2(x,prec)
    mpfr_set_ui(norm,0,rnd_re)
    if ntype==2:
        for i from 0 <= i < nrows:
            for j from 0 <= j < ncols:
                mpc_abs(x,A[i][j],rnd_re)
                mpfr_sqr(x,x,rnd_re)
                mpfr_add(norm,norm,x,rnd_re)
        mpfr_sqrt(norm,norm,rnd_re)
    elif ntype==0: ## infinity or max norm
        for i from 0 <= i < nrows:
            for j from 0 <= j < ncols:
                #print "|A[",i*ncols+j,"]"
                mpc_abs(x,A[i][j],rnd_re)
                if mpfr_cmp(x,norm)>0:
                    mpfr_set(norm,x,rnd_re)
    else:
        mpfr_set_si(norm,-1,rnd_re)
        # raise NotImplementedError,"Only infinty and 2-norm is currently implemented"
    mpfr_clear(x)







cdef _get_additional_precision(mpc_t** A,int m, int n, mpfr_t mpeps, int prec, mpc_rnd_t rnd, mpfr_rnd_t rnd_re):
    r"""
    Computes the number of extra bits of precision needed
    to compute eigenvalues accuraetly.
    """
    cdef mpfr_t amax,x,v,y
    cdef int i,j
    mpfr_init2(x,prec); mpfr_init2(y,prec); mpfr_init2(v,prec)
    mpfr_init2(amax,prec)
    _norm(amax,A,m,n,prec,rnd_re,0)
    #for i from 0 <= i < n:
    #    mpc_abs(x,A[i][i],rnd_re)
    #    for j from i+1 <= j < n:
    #        mpc_abs(y,A[i][j],rnd_re)
    ##        if mpfr_cmp(v,mpeps)<0:
    #            continue
    #        mpfr_div(v,y,x,rnd_re)
    #        if mpfr_cmp(v,amax)>0:
    #            mpfr_set(amax,v,rnd_re)
    mpfr_mul(amax,amax,amax,rnd_re)
    mpfr_mul_ui(amax,amax,n,rnd_re)
    mpfr_set_d(x,3.01,rnd_re)
    mpfr_mul(amax,amax,x,rnd_re)
    #if deb:
    #    mpfr_set(tmpr.value,amax,rnd_re)
    #    print "amax=",tmpr
    mpfr_log2(amax,amax,rnd_re)
    mpfr_abs(amax,amax,rnd_re)
    mpfr_add_ui(amax,amax,1,rnd_re)
    return  mpfr_get_ui(amax,GMP_RNDD)+30 # truncating


###
### Helper functions for complex arithmetic
###


cdef void mm_get_delta(mpfr_t *delta,mpfr_rnd_t rnd_re):
    mpfr_set_si(delta[0], 1, rnd_re)
    mpfr_nextabove(delta[0])
    mpfr_sub_si(delta[0], delta[0], 1, rnd_re)


## cdef inline void mpc_z_to_pol(mpc_t * res, mpc_t z,mpfr_rnd_t rnd_re):
##     #assert_nsame(res[0].re, z.re)
##     mpfr_hypot(res[0].re, z.re, z.im, rnd_re)
##     mpfr_atan2(res[0].im, z.im, z.re, rnd_re)



## ## This is about 20% faster than standard mpc_mul
## cdef inline void _mpc_mul(mpc_t* z, mpc_t a, mpc_t b, mpc_t *t, mpc_rnd_t rnd, mpfr_rnd_t rnd_re):
##     # __inline__ void compl_mul(complex_t * z, complex_t a, complex_t b)
##     #assert_nsame(z->re, a.re);
##     #assert_nsame(z->re, b.re);
##     if mpfr_zero_p(a.re):
##         if mpfr_zero_p(a.im):
##             mpc_set_si(z[0], 0, rnd)
##             return
##         if mpfr_zero_p(b[0].im):
##             mpfr_set_si(z[0].re, 0, rnd_re)
##         else:
##             mpfr_mul(z[0].re, a.im, b[0].im, rnd_re)
##             mpfr_neg(z[0].re, z[0].re, rnd_re)            
##         if mpfr_zero_p(b[0].re):
##             mpfr_set_si(z[0].im, 0, rnd_re)
##         else:
##             mpfr_mul(z[0].im, a.im, b[0].re, rnd_re)
##         return
##     if mpfr_zero_p(a.im):
##         if mpfr_zero_p(b[0].re):
##             mpfr_set_si(t[0].re, 0, rnd_re)
##         else:
##             mpfr_mul(t[0].re, a.re, b[0].re, rnd_re)
##         if mpfr_zero_p(b[0].im):
##             mpfr_set_si(t[0].im, 0, rnd_re)
##         else:
##             mpfr_mul(t[0].im, a.re, b[0].im, rnd_re)
##         mpfr_set(z[0].re,t[0].re,rnd_re)
##         mpfr_set(z[0].im,t[0].im,rnd_re)
##         return
##     if mpfr_zero_p(b[0].re):
##         if mpfr_zero_p(b[0].im):
##             mpc_set_si(z[0], 0, rnd)
##             return
##         mpfr_mul(t[0].re, a.im, b[0].im, rnd_re)
##         mpfr_neg(t[0].re, t[0].re, rnd_re)
##         mpfr_mul(z[0].im, a.re, b[0].im, rnd_re)
##         mpfr_set(z[0].re,t[0].re,rnd_re)
##         return
##     if mpfr_zero_p(b[0].im):
##         mpfr_mul(z[0].re, a.re, b[0].re, rnd_re)
##         mpfr_mul(z[0].im, a.im, b[0].re, rnd_re)
##         return
##     mpfr_mul(t[0].re, a.re, b[0].re, rnd_re)
## #/* z[0].re = a.re*b[0].re */
##     mpfr_mul(t[0].im, a.im, b[0].im, rnd_re)
## #/* z[0].im = a.im*b[0].im */
##     mpfr_sub(t[0].re, t[0].re, t[0].im, rnd_re)
## #/* z->re = a.re*b[0].re - a.im*b[0].im */

##     mpfr_mul(t[0].im, a.im, b[0].re, rnd_re)
## #/* z->im = a.im*b[0].re */
##     mpfr_fma(z[0].im, a.re, b[0].im, t[0].im, rnd_re)	#/* z->im = a.re*b[0].im + a.im*b[0].re */
##     mpfr_set(z[0].re,t[0].re,rnd_re)




## cdef inline void _mpc_mul_fr(mpc_t* z, mpc_t a, mpfr_t b, mpc_rnd_t rnd, mpfr_rnd_t rnd_re):
##     cdef int kk
##     mpfr_mul(z[0].re, a.re, b, rnd_re)
##     mpfr_mul(z[0].im, a.im, b, rnd_re)
##     return


## cdef inline void _mpc_div(mpc_t * z, mpc_t a, mpc_t b, mpc_t t[2], mpfr_rnd_t rnd_re):
##     r"""

    
    
##     z = a/b
##     This works if z is a pointer to a but not b.
##     """
##     # note that this is  if we apply divsion by self
##     # t is a temporary variable
##     #	assert_nsame(z->re, a.re);
##     #	assert_nsame(z->re, b.re);
##     # if 1
##     #	m_assert(!mpfr_zero_p(b.re) || !mpfr_zero_p(b.im))
##     if mpfr_zero_p(a.re):
##         if mpfr_zero_p(a[0].im): 
##             mpc_set_si(z[0], 0, rnd)
##             return
##         if mpfr_zero_p(b[0].re):
##             mpfr_set_si(z[0].im, 0, rnd_re)
##             mpfr_div(z[0].re, a[0].im, b[0].im, rnd_re)
##             return
            
##         if mpfr_zero_p(b[0].im):
##             mpfr_set_si(z[0].re, 0, rnd_re)
##             mpfr_div(z[0].im, a[0].im, b[0].re, rnd_re)
##             return
##         mpfr_mul(t[0].re, a[0].im, b[0].im, rnd_re)
##         mpfr_mul(t[0].im, a[0].im, b[0].re, rnd_re)
        
##         mpfr_mul(t[1].re, b[0].re, b[0].re, rnd_re)
##         mpfr_fma(t[1].re, b[0].im, b[0].im, t[0].re, rnd_re)
##         mpfr_div(z[0].re, t[0].re, t[1].re, rnd_re)
##         mpfr_div(z[0].im, t[0].im, t[1].re, rnd_re)
##         return
    

##     if mpfr_zero_p(a[0].im):
##         if mpfr_zero_p(b[0].re):
##             mpfr_set_si(z[0].re, 0, rnd_re)
##             mpfr_div(z[0].im, a[0].re, b[0].im, rnd_re)
##             mpfr_neg(z[0].im, z[0].im, rnd_re)
##             return
##         if mpfr_zero_p(b[0].im):
##             mpfr_set_si(z[0].im, 0, rnd_re)
##             mpfr_div(z[0].re, a[0].re, b[0].re, rnd_re)
##             return
##         mpfr_mul(t[0].re, a[0].re, b[0].re, rnd_re)
##         mpfr_mul(t[0].im, a[0].re, b[0].im, rnd_re)
##         mpfr_neg(t[0].im, t[0].im, rnd_re)
        
##         mpfr_mul(t[1].re, b[0].re, b[0].re, rnd_re)
##         mpfr_fma(t[1].re, b[0].im, b[0].im, t[0].re, rnd_re)
##         mpfr_div(z[0].re, t[0].re, t[1].re, rnd_re)
##         mpfr_div(z[0].im, t[0].im, t[1].re, rnd_re)
##         return
##     if mpfr_zero_p(b[0].re):
##         mpfr_div(t[0].re, a[0].im, b[0].im, rnd_re)
        
##         mpfr_div(t[0].im, a[0].re, b[0].im, rnd_re)
##         mpfr_neg(t[0].im, t[0].im, rnd_re)
##         mpfr_set(z[0].re,t[0].re,rnd_re)
##         mpfr_set(z[0].im,t[0].im,rnd_re)
##         return
##     if mpfr_zero_p(b[0].im):
##         _mpc_div_fr(&z[0], a, b[0].re,rnd_re)
##         return
##     # endif
##     mpfr_mul(t[0].im, a[0].im, b[0].re, rnd_re)
## #/* t[0].im = a.im*b.re */
##     mpfr_mul(t[0].re, a[0].re, b[0].im, rnd_re)
## #/* t[0].re = a.re*b.im */
##     mpfr_sub(t[0].im, t[0].im, t[0].re, rnd_re)
## #/* z.im = a.im*b.re - a.re*b.im */
    
##     mpfr_mul(t[0].re, a[0].re, b[0].re, rnd_re)
## #/* t[0].re = a.re*b.re */
##     mpfr_fma(t[0].re, a[0].im, b[0].im, t[0].re, rnd_re)
## #/* z.im = a.re*b.re + a.im*b.im */
##     mpfr_mul(t[1].re, b[0].re, b[0].re, rnd_re)
##     mpfr_fma(t[1].re, b[0].im, b[0].im, t[1].re, rnd_re)
    
##     assert(not mpfr_zero_p(t[0].re))
##     mpfr_div(z[0].re, t[0].re, t[1].re, rnd_re)
##     mpfr_div(z[0].im, t[0].im, t[1].re, rnd_re)




## cdef inline void _mpc_add(mpc_t *res, mpc_t a, mpc_t b,mpfr_rnd_t rnd_re):
##     mpfr_add(res[0].re, a[0].re, b[0].re, rnd_re)
##     mpfr_add(res[0].im, a[0].im, b[0].im, rnd_re)


## cdef inline void _mpc_sub(mpc_t *res, mpc_t a, mpc_t b,mpfr_rnd_t rnd_re):
##     mpfr_sub(res[0].re, a[0].re, b[0].re, rnd_re)
##     mpfr_sub(res[0].im, a[0].im, b[0].im, rnd_re)


## cdef inline void _mpc_div_fr(mpc_t * r, mpc_t a, mpfr_t d, mpfr_rnd_t rnd):
##     mpfr_div(r[0].re, a[0].re, d, rnd)
##     mpfr_div(r[0].im, a[0].im, d, rnd)


## cdef inline void _mpc_div_ui(mpc_t *res,mpc_t z, unsigned int i, mpfr_rnd_t rnd_re):
##     mpfr_div_ui(res[0].re,z.re,i,rnd_re)
##     mpfr_div_ui(res[0].im,z.im,i,rnd_re)



## cdef inline void _mpc_set(mpc_t *res, mpc_t z, mpfr_rnd_t rnd_re):
##     mpfr_set(res[0].re, z[0].re,rnd_re)
##     mpfr_set(res[0].im, z[0].im,rnd_re)

## cdef inline void _mpc_conj(mpc_t *res, mpc_t z, mpfr_rnd_t rnd_re):
##     mpfr_set(res[0].re, z[0].re,rnd_re)
##     mpfr_set(res[0].im, z[0].im,rnd_re)
##     mpfr_neg(res[0].im, res[0].im,rnd_re)

## # Note: mpc_abs and mpc_sqrt are more or less optimal.
## #       we don't need to make any new versions of these


## cdef inline void _pochammer(mpc_t *res, mpc_t z, int n,mpc_rnd_t rnd,mpfr_rnd_t rnd_re):
##     r"""
##     Pochammer symbol
##     """
##     cdef int j
##     cdef int prec= mpc_get_prec(z)
##     cdef mpc_t x
##     cdef mpc_t t[2]
##     mpc_init2(t[0],prec+10)
##     mpc_init2(t[1],prec+10)
##     mpc_init2(x,prec+10)
##     mpc_set_ui(res[0],1,rnd)
##     for j from 0 <= j < n:
##         #mpc_add_ui(z.re,z.re,1,rnd_re)
##         mpc_add_ui(x,z,j,rnd)
##         _mpc_mul(res,res[0],x,t,rnd,rnd_re)
##     mpc_clear(t[0])
##     mpc_clear(t[1])
##     mpc_clear(x)

## def mpc_pochammer(MPComplexNumber z,n):
##     assert isinstance(z,MPComplexNumber)
##     cdef mpc_t x
##     cdef MPComplexNumber res
##     #assert ZZ(n)==n
##     F = z.parent()
##     mpc_init2(x,F.prec())
##     res = F.one()
##     rnd_re = _mpfr_rounding_modes.index(F.rounding_mode_real())
##     rnd = _mpc_rounding_modes.index(F.rounding_mode())
##     _pochammer(&x,z.value,n,rnd,rnd_re)
##     mpc_set(res.value,x,rnd)
##     mpc_clear(x)
##     return res

## cdef print_mat(mpc_t **A,int m,int n): #, mpc_rnd_t rnd):
##     r""" """

##     from sage.rings.complex_mpc import MPComplexField
##     cdef MPComplexNumber z
##     cdef int i,j,prec
##     prec = mpc_get_prec(A[0][0])
##     z = MPComplexField(prec)(1)
##     s = "["
##     for i from 0 <= i < m:
##         s+="["
##         for j from 0 <= j < m:        
##             mpc_set(z.value,A[i][j],MPC_RNDNN)
##             s+=str(z)
##             if j<m-1:
##                 s+=", "
##         s+="]"
##         if i<m-1:
##             s+=",\n"
##     s+="]"
##     print s

## cdef print_mpc(mpc_t z):
##     from sage.rings.complex_mpc import MPComplexField
##     cdef int prec
##     cdef MPComplexNumber zz
##     prec = mpc_get_prec(z)
##     zz = MPComplexField(prec)(1)
##     mpc_set(zz.value,z,MPC_RNDNN)
##     return str(zz)

## def print_mpfr(mpfr_t x):
##     cdef int prec
##     cdef RealNumber xx
##     prec = mpfr_get_prec(x)
##     xx = RealField(prec)(1)
##     mpfr_set(xx.value,x,GMP_RNDN)
##     return str(xx)


## cdef void mpc_sign(mpc_t *res, mpc_t *z,mpfr_rnd_t rnd_re):
##     r"""
##     The complex sign function. Returns z/abs(z)=exp(i*Arg(z))
##     """
##     mpc_abs(res[0].re, z[0],rnd_re)
##     if mpfr_zero_p(res[0].re):
##         mpfr_set_ui(res[0].re,1,rnd_re)
##         mpfr_set_ui(res[0].im,0,rnd_re)
##         #assert not 
##     else:
##         mpfr_div(res[0].im, z[0].im, res[0].re, rnd_re)
##         mpfr_div(res[0].re, z[0].re, res[0].re, rnd_re)


## ##
## ## Test the standard mpc_mul and mpc_add against the
## ## above versions.
## ## 
## cpdef mul1(F):
##     cdef MPComplexNumber a,b
##     cdef RealNumber d
##     a = F.random_element()
##     b = F.random_element()
##     d = F.base().random_element()
##     cdef mpc_t x,y,s
##     cdef mpfr_t xx
##     mpc_init2(x,F.prec())
##     mpc_init2(y,F.prec())
##     mpc_init2(s,F.prec())
##     mpfr_init2(xx,F.prec())
##     mpc_set(x,a.value,MPC_RNDNN)
##     mpfr_set(xx,d.value,GMP_RNDN)
##     #mpc_div(s,x,y,MPC_RNDNN)
##     cdef int i
##     for i from 0 <= i <= 10000:
##         mpc_div_ui(s,x,2,MPC_RNDNN)
##         #mpc_conj(x,y,MPC_RNDNN)
##         #mpc_set(x,y,MPC_RNDNN)


## cpdef mul3(F):
##     cdef MPComplexNumber a,b,c
##     cdef mpc_t x, y,z
##     cdef mpc_t t[3]
    
##     a = F.random_element()
##     b = F.random_element()
##     mpc_init2(x,F.prec())
##     mpc_init2(t[0],F.prec())
##     mpc_init2(t[1],F.prec())
##     mpc_init2(t[2],F.prec())
##     mpc_init2(y,F.prec())
##     mpc_init2(z,F.prec())

##     mpc_set(y,a.value,MPC_RNDNN)
##     _mpc_mul(&z,x,y,t,MPC_RNDNN,GMP_RNDN)
##     print "x*y=",print_mpc(z)
##     mpc_mul(z,x,y,<mpc_rnd_t>MPC_RNDNN)
##     print "x*y=",print_mpc(z)

## cpdef mul2(F):
##     cdef MPComplexNumber a,b,c
##     cdef RealNumber d
##     a = F.random_element()
##     b = F.random_element()
##     c = F.random_element()
##     d = F.base().random_element()
##     cdef mpc_t x, y, s,z
##     cdef mpfr_t xx
##     cdef mpc_rnd_t rnd
##     rnd = MPC_RNDNN
##     cdef mpfr_rnd_t rnd_re
##     rnd_re = GMP_RNDN
##     cdef int i
##     mpc_init2(x,F.prec())
##     mpc_init2(y,F.prec())
##     mpc_init2(s,F.prec())
##     mpc_init2(z,F.prec())
##     mpfr_init2(xx,F.prec())
##     mpc_set(x,a.value,MPC_RNDNN)
##     mpfr_set(xx,d.value,GMP_RNDN)
##     #print "a=",a
##     #print "b=",b
##     #_mpc_div(&s,x,y,&z,rnd_re)
##     #mpc_set(c.value,s,MPC_RNDNN)
##     #print "ab=",c
##     #mpc_div(s,x,y,MPC_RNDNN)    
##     #mpc_set(c.value,s,MPC_RNDNN)    
##     #print "ab=",c
##     #return 

##     for i from 0 <= i <= 10000:
##         _mpc_div_ui(&s,x,2,rnd_re)
##         #_mpc_abs(&xx,&x,rnd_re)
##         #mpc_conj(x,y,rnd_re)
##         #mpfr_set(x.re, y.re,rnd_re)
##         #mpfr_set(x.im, y.im,rnd_re)
##         #_mpc_set(&x, y,rnd_re)
##         #mpfr_neg(x.im, x.im,rnd_re)
    


## cpdef add1(F):
##     cdef MPComplexNumber a,b
##     a = F.random_element()
##     b = F.random_element()
##     cdef mpc_t x,y,s
##     cdef mpfr_rnd_t rnd_re
##     rnd_re = GMP_RNDN
##     mpc_init2(x,F.prec())
##     mpc_init2(y,F.prec())
##     mpc_init2(s,F.prec())
##     mpc_set(x,a.value,MPC_RNDNN)
##     mpc_set(y,b.value,MPC_RNDNN)
##     cdef int i
##     for i from 0 <= i <= 10000:
##         mpc_sub(s,x,y,MPC_RNDNN)

## cpdef add2(F):
##     cdef MPComplexNumber a,b
##     a = F.random_element()
##     b = F.random_element()
##     cdef mpc_t x, y, s
##     mpc_init2(x,F.prec())
##     mpc_init2(y,F.prec())
##     mpc_init2(s,F.prec())
##     mpc_set(x,a.value,MPC_RNDNN)
##     mpc_set(y,b.value,MPC_RNDNN)
##     cdef mpfr_rnd_t rnd
##     rnd = GMP_RNDN
##     cdef int i 
##     for i from 0 <= i <= 10000:
##         _mpc_sub(&s,x,y,rnd)


