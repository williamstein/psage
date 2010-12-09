#ifndef _JACORDER_INCLUDE_
#define _JACORDER_INCLUDE_

#include <limits.h>

/*
    Copyright 2007 Andrew V. Sutherland

    This file is part of smalljac.

    smalljac is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    smalljac is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with smalljac.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
	This module contains generic order computations for Jacobians using BSGS.
	The detail are described in [SutherlandThesis] and [KedlayaSutherland2007].
	
	Most of smalljac's life is spent in the function jac_parallel_search.
*/

#define JAC_INVALID_A1			LONG_MAX

struct a2tab_entry {
	double a1;
	double a2median;
	double a2mae;
};

int jac_order (unsigned long *pP1, unsigned long Min, unsigned long Max,  long a1, int *constraints, curve_t c[1], int fExponentOnly);
unsigned long jac_parallel_search (jac_t a[1],  unsigned long M, unsigned long W, unsigned long Min, unsigned long Max, int tiny, int parity, curve_t c[1]);
int jac_search (mpz_t e[2], jac_t a[1],  unsigned m, mpz_t Min, mpz_t Max, curve_t c[1]);

#endif
