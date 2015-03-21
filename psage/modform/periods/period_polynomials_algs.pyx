"""
Cython routines needed for computation of period integrals associated to spaces of modular symbols.
"""

#############################################################################
#       Copyright (C) 2012 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

cdef extern:
    complex cexp(complex)
    double exp(double)
    
import math
cdef complex c0 = complex(0,2*math.pi)  # 2*pi*i

from sage.rings.complex_number cimport ComplexNumber
from sage.rings.real_mpfr cimport RealNumber
from sage.rings.real_mpfr import RealField
from sage.rings.complex_field import ComplexField



cpdef complex exp_z_integral(complex alpha, unsigned long n, unsigned int m):
    """
    Numerically compute the integral

       integral_{alpha}^{infty} exp(2*pi*i*n*z) z^m dz

    where alpha is complex, n>=1 an integer, and m>=0 an integer.
    This uses Lemma 10.4 from Stein's 2007 AMS modular forms book.

    Note that m is typically small but n gets large.
    """
    cdef complex t = n*c0                  # = 2*pi*i*n
    cdef complex one_over_t = 1/t
    cdef double j_prod = 1                 # = prod_{(m+1)-s}^m j
    cdef complex denom = 1/t               # = 1/(2*pi*i*n)^(s+1)
    cdef complex alpha_pow = alpha**m      # = alpha^(m-s)
    cdef complex alpha_inv = 1/alpha       
    cdef int sgn = -1                       # = (-1)^s

    cdef unsigned int s

    cdef complex summation = sgn*alpha_pow*denom*j_prod
    
    for s in range(1,m+1):
        j_prod *= m+1-s
        denom *= one_over_t
        sgn *= -1
        alpha_pow *= alpha_inv
        summation += sgn*alpha_pow*denom*j_prod

    return cexp(t*alpha) * summation
        
cpdef complex extended_period_integral(unsigned int m, complex alpha, list v):
    """
    Entries of v = [a0,a1,a2,...] are assumed to be complex.
    """
    # There are many nicer ways that this code could be written, e.g., using
    # sum, range, enumerate, etc. -- don't bother, as they are all way slower,
    # at least with Cython 0.15.1.
    cdef complex summation = 0
    cdef unsigned long n = 1
    cdef complex z
    for z in v[1:]:
        summation += z * exp_z_integral(alpha, n, m)
        n += 1
    return summation * c0

    
                              
#### Extensions by me

cpdef complex exp_z_integral_w(complex alpha, unsigned long n, unsigned long width,unsigned int m):
    """
    Numerically compute the integral

       integral_{alpha}^{infty} exp(2*pi*i*n*z) z^m dz

    where alpha is complex, n>=1 an integer, and m>=0 an integer.
    This uses Lemma 10.4 from Stein's 2007 AMS modular forms book.

    Note that m is typically small but n gets large.
    """
    cdef complex t = n*c0/<double>width                  # = 2*pi*i*n
    cdef complex one_over_t = 1/t
    cdef double j_prod = 1                 # = prod_{(m+1)-s}^m j
    cdef complex denom = 1/t               # = 1/(2*pi*i*n)^(s+1)
    cdef complex alpha_pow = alpha**m      # = alpha^(m-s)
    cdef complex alpha_inv = 1/alpha       
    cdef int sgn = -1                       # = (-1)^s

    cdef unsigned int s

    cdef complex summation = sgn*alpha_pow*denom*j_prod
    
    for s in range(1,m+1):
        j_prod *= m+1-s
        denom *= one_over_t
        sgn *= -1
        alpha_pow *= alpha_inv
        summation += sgn*alpha_pow*denom*j_prod

    return cexp(t*alpha) * summation


cpdef ComplexNumber exp_z_integral_wmp(ComplexNumber alpha, unsigned long n, unsigned long width,unsigned int m):
    """
    Numerically compute the integral

       integral_{alpha}^{infty} exp(2*pi*i*n*z) z^m dz

    where alpha is complex, n>=1 an integer, and m>=0 an integer.
    This uses Lemma 10.4 from Stein's 2007 AMS modular forms book.

    NOTE: There is a sign wrong in the formula of the Lemma
    Here the sign is corrected!
    
    Note that m is typically small but n gets large.
    """
    CF = alpha.parent()
    cdef int prec = CF.prec()
    RF = RealField(prec)
    cdef ComplexNumber t
    t = n*CF(0,2*RF.pi())/RF(width)                  # = 2*pi*i*n
    cdef ComplexNumber one_over_t
    one_over_t = CF(t)**-1
    cdef RealNumber j_prod
    j_prod = RF(1)                 # = prod_{(m+1)-s}^m j
    cdef ComplexNumber denom
    denom = CF(t)**-1              # = 1/(2*pi*i*n)^(s+1)
    cdef ComplexNumber alpha_pow
    alpha_pow = CF(alpha)**CF(m)      # = alpha^(m-s)
    cdef ComplexNumber alpha_inv
    alpha_inv = CF(alpha)**-1
    cdef int sgn = -1                       # = (-1)^s

    cdef unsigned int s

    cdef ComplexNumber summation
    summation = CF(sgn)*alpha_pow*denom*j_prod
    
    for s in range(1,m+1):
        j_prod *= m+1-s
        denom *= one_over_t
        sgn *= -1
        alpha_pow *= alpha_inv
        summation += sgn*alpha_pow*denom*j_prod

    return (t*alpha).exp() * summation

cpdef int get_truncation(int k,int L,int prec):
    r"""
    Compute a lower bound on M s.t. if we truncate the q-expansion
    of a holomorphic modular form of weight k:
    f(q)=\Sum_{n=1}^{oo} a(n) q^n, q=e(z)  at n=M then the error is less than 2**-prec for all z=x+iy wuth y>=1/N
    
    """
    cdef int maxit,M,M0
    cdef double eps,argk,f1,f2,Lf
    cdef double twopi = 6.28318530717959
    maxit = 1000*max(k,2)  ## To avoid infinite loops..
    eps = 2.0**-<double>prec
    argk = <double>(k-1)/<double>(2)
    Lf = <double>L
    f1 = (Lf/twopi)**argk
    argk = argk-2.0 #RF(3)
    if argk>10.0:
        M0 = math.ceil(argk)
    else:
        M0 = 10
    #M0 = max(10,argk)    
    argk = argk-1.0 #RF(3)
    for M in range(M0,maxit):
        argM = twopi*<double>(M)/Lf
        ## We are to the right of the critical point so we can use the asymptotic estimate.
        f2 = argM**argk*exp(-argM)  
        #f2 = incomplete_gamma(argk,argM)
        if f1*f2 < eps:
            return M
    raise ArithmeticError,"Could not find truncation for prec={0}".format(prec)


    
