#ifndef _JAC_INCLUDE_
#define _JAC_INCLUDE_

#include <stdint.h>
#include "hecurve.h"
#include "ff.h"
#include "ffpoly.h"

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
	Basic group operations for Jacobians are defined here, plus some generic
	fastorder computation functions (computations of element orders from
	a known exponent).
*/	

#define JAC_MAX_CURVES			100
#define JAC_MAX_GAPS				64
#define JAC_BUFSIZE				4096
#define JAC_MAX_FACTORS			100
#define JAC_U_DEGREE				HECURVE_GENUS
#define JAC_V_DEGREE				HECURVE_GENUS-1
#define JAC_CONFIDENCE_LEVEL		20			// larger value used internally in jac_sylow - note that smalljac results are (except where noted) provably correct regardless of this setting
#define JAC_MAX_GENERATORS		20			// this should be 2g+1 (not 2g, we need room for one extra) but for debugging purposes we set it higher then needed
#define JAC_CYCLIC_TESTS			4			// number of times to attempt to prove p-Sylow subgroup is cyclic before computing it

unsigned long jac_gops;

typedef struct {
	ff_t u[JAC_U_DEGREE+1];
	ff_t v[JAC_V_DEGREE+1];
} jac_t;

struct _jac_vector_struct {
	jac_t a[JAC_MAX_CURVES];
};
typedef struct _jac_vector_struct jac_vec_t;

struct _curve_struct {
	ff_t f[HECURVE_DEGREE+1];
	int d;
};
typedef struct _curve_struct curve_t;


#define _curve_copy(c1,c2)		((c1)=(c2))
#define _curve_print(c)			ff_poly_print((c).f, (c).d)
#define _curve_sprint(s,c)		ff_poly_sprint(s,(c).f, (c).d)

#define _jac_set(o,a)			_hecurve_set((o).u,(o).v,(a).u,(a).v)
#define __jac_mult(o,a,b,c,ctx)	hecurve_compose ((o).u, (o).v, (a).u, (a).v, (b).u, (b).v, (c).f, ctx)
#define __jac_square(o,a,c,ctx)	hecurve_square ((o).u, (o).v, (a).u, (a).v, (c).f, ctx)
#define _jac_mult(o,a,b,c)		{ __jac_mult (o,a,b,c,0);  jac_gops++; }
#define _jac_square(o,a,c)		{ __jac_square (o,a,c,0);  jac_gops++; }
#define _jac_invert(o)			hecurve_invert((o).u,(o).v)
#define _jac_random(o,c)		hecurve_random((o).u,(o).v,(c).f)

#define _jac_cmp(a,b)			hecurve_cmp ((a).u,(a).v,(b).u,(b).v)  	// 1 if equal, -1 if inverses, 0 o.w.
#define _jac_set_identity(a)		_hecurve_set_identity ((a).u,(a).v)
#define _jac_is_identity(a)		_hecurve_is_identity ((a).u,(a).v)
#define _jac_2tor(a)			_hecurve_2tor((a).u,(a).v)
#define _jac_bits(a,c)			hecurve_bits ((a).u,(a).v,(c).f)

#define _jac_print(a)			hecurve_print((a).u, (a).v)
#define _jac_sprint(s,a)			hecurve_sprint(s,(a).u,(a).v)

static inline void jac_parallel_set (jac_t a[], jac_t b[], int n) { register int i;  for ( i = 0 ; i < n ; i++ ) { _jac_set(a[i],b[i]); }}

void curve_random (curve_t c[1]);

void jac_parallel_mult (jac_t o[], jac_t a[], jac_t b[], curve_t c[], int n);
void jac_parallel_mult_c (jac_t o[], jac_t a[], jac_t b[], curve_t c[1], int n);		// multiple a's and b's but same c
void jac_parallel_mult_1 (jac_t o[], jac_t a[], jac_t b[1], curve_t c[1], int n);
void jac_parallel_square (jac_t o[], jac_t a[], curve_t c[], int n);

void jac_square_mult (jac_t a[1], jac_t b[1], curve_t c[1]);
void jac_mult2 (jac_t o1[1], jac_t a1[1], jac_t b1[1], jac_t o2[1], jac_t a2[1], jac_t b2[1], curve_t c[1]);	// o1 can't overlap a2 or b2

void jac_exp_ui (jac_t o[1], jac_t a[1], unsigned long e, curve_t c[1]);
void jac_exp2_ui (jac_t o[2], jac_t a[1], unsigned long e1, unsigned long e2, curve_t c[1]);
void jac_exp_ui32_powers (jac_t o[1], jac_t a[], uint32_t e, curve_t c[1]);
void jac_exp_mpz (jac_t o[1], jac_t a[1], mpz_t e, curve_t c[1]);
int jac_verify_group_exponent (mpz_t e, curve_t c[1]);
unsigned long jac_pp_order_ui (jac_t a[1], unsigned long p, unsigned long h, curve_t c[1]);
int jac_fastorder_ui (unsigned long *po, jac_t a[1], unsigned long e, curve_t c[1]);
int jac_factored_order_ui (unsigned long *po, jac_t a[1], unsigned long p[], unsigned long h[], int w, curve_t c[1]);
unsigned long jac_fastorder_powersmooth (mpz_t o, jac_t a[1], unsigned long L, curve_t c[1]);
int jac_sylow (jac_t a[JAC_MAX_GENERATORS], unsigned long ords[JAC_MAX_GENERATORS], unsigned p, unsigned long E, unsigned long M, unsigned long maxgops, curve_t c[1]);
int jac_structure (long m[], curve_t c[1], unsigned long order, int fExponent);

#endif
