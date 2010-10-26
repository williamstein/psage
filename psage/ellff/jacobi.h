#ifndef JACOBI_H
#define JACOBI_H

#include <NTL/lzz_p.h>
#include <NTL/lzz_pE.h>
#include <NTL/ZZX.h>

NTL_CLIENT

void jacobi_sum(ZZ ***sum, long d, long *chi_of_16);

#endif	// JACOBI_H
