#include "sage/rings/mpc.pxi"
from psage.modules.vector_complex_dense cimport Vector_complex_dense
from psage.modules.vector_real_mpfr_dense cimport Vector_real_mpfr_dense
from psage.matrix.matrix_complex_dense cimport Matrix_complex_dense
#from psage.rings.mpc_extras cimport *
from psage.rings.mpc_extras cimport _mpc_mul,print_mpc,print_mpfr,_mpc_sub,_mpc_set,_mpc_add
