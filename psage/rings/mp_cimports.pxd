r"""
Standard cimports for working with multiprecision numbers.


"""
from sage.libs.mpfr cimport *
from sage.libs.mpc cimport *

from sage.rings.complex_mpc cimport * 
from sage.rings.complex_mpc cimport MPComplexNumber
from sage.rings.real_mpfr cimport RealNumber,RealField_class
from sage.rings.complex_number cimport ComplexNumber
from sage.rings.integer cimport Integer
from sage.rings.complex_mpc cimport MPComplexNumber
from sage.rings.complex_mpc cimport MPComplexField_class
from sage.rings.real_mpfr cimport RealNumber

## Define the rounding modes
cdef mpc_rnd_t rnd
cdef mpfr_rnd_t rnd_re
rnd = MPC_RNDNN
rnd_re = MPFR_RNDN
