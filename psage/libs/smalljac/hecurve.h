#ifndef _HECURVE_INCLUDE_
#define _HECURVE_INCLUDE_

#include "gmp.h"
#include "smalljac.h"
#include "ff.h"

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

// some counters used for performance diagnostics
extern unsigned long hecurve_fastorders;
extern unsigned long hecurve_steps;
extern unsigned long hecurve_fastorder_expbits;
extern unsigned long hecurve_expbits;
extern unsigned long hecurve_retries;
extern unsigned long hecurve_g1_tor3_counter[10];
extern unsigned long hecurve_g1_tor4_counter[17];
extern unsigned long hecurve_g1_a_counter;
extern unsigned long hecurve_g1_s2known_counter;


/*
	Implementation of Jacobian group operation for genus 1, 2 and 3 hyperelliptic curves over prime fields
	Currently the genus is fixed at compile time.  This is planned to change.  This module is the common
	entry point for all genera and calls the appropriate genus specific code in hecurve1.c, hecurve2.c, hecurve3.c
*/

#define HECURVE_GENUS		SMALLJAC_GENUS
#define HECURVE_DEGREE		(2*HECURVE_GENUS+1)

#define HECURVE_RANDOM		0		// generate truly random divisors
#define HECURVE_VERIFY		0		// verifies that all inputs and outputs are valid elements of the Jacobian - ONLY USE FOR DEBUGGING THIS IS VERY SLOW
#define HECURVE_FAST			1		// define to turn off certain error checking.  see also FF_FAST (only marginally effects performance)
#define HECURVE_SPARSE		0		// assumes f2 through f_2g are zero (makes only about 1-2 % difference)

#define HECURVE_RANDOM_POINT_RETRIES	100			// this fails with probability 1/2, and we really need it to succeed
#define HECURVE_G1_ORDER_RETRIES		20
#define HECURVE_G1_SHORT_INTERVAL		5			// don't bother with full BSGS search for very short intervals
#define HECURVE_G1_4TOR_MINP			(1<<23)		// don't compute 4-torsion for p smaller than this, its not worth the effort
#define HECURVE_G1_8TOR_MINP			(1<<26)		// don't check any 8-torsion for p smaller than this, its not worth the effort

#if HECURVE_GENUS == 1
// note that rep in genus 1 accomodates projective coords, so we use u[1]=z and have (x,y,z) => u[0]=-x, u[1]=z, v[0]=y (note the sign difference of u[0] and x)
#define _hecurve_init_uv(u,v)			{ _ff_init(u[0]);  _ff_init(v[0]);  }
#define _hecurve_set(u1,v1,u2,v2)		{_ff_set(u1[1],u2[1]);_ff_set(u1[0],u2[0]);_ff_set(v1[0],v2[0]);}
#define _hecurve_is_identity(u,v)		_ff_zero((u)[1])
#define _hecurve_set_identity(u,v) 		{_ff_set_zero((u)[1]);  _ff_set_zero((u)[0]);  _ff_set_zero((v)[0]); }
#define _hecurve_2tor(u,v)			(_ff_zero((v)[0]))
#endif
#if HECURVE_GENUS == 2
#define _hecurve_init_uv(u,v)			{ _ff_init(u[0]);  _ff_init(u[1]);  _ff_init(u[2]);  _ff_init(v[0]);  _ff_init(v[1]);  }
#define _hecurve_set(u1,v1,u2,v2)		{_ff_set(u1[2],u2[2]);_ff_set(u1[1],u2[1]);_ff_set(u1[0],u2[0]);_ff_set(v1[1],v2[1]);_ff_set(v1[0],v2[0]);}
#define _hecurve_is_identity(u,v)		(_ff_zero((u)[1]) && _ff_zero((u)[2]))			// ignores v
#define _hecurve_set_identity(u,v) 		{_ff_set_one((u)[0]); _ff_set_zero((u)[1]);  _ff_set_zero((u)[2]);  \
								  _ff_set_zero((v)[0]);  _ff_set_zero((v)[1]);}
#define _hecurve_2tor(u,v)			(_ff_zero((v)[1]) && _ff_zero((v)[0]))
#endif
#if HECURVE_GENUS == 3
#define _hecurve_init_uv(u,v)			{ _ff_init(u[0]);  _ff_init(u[1]);  _ff_init(u[2]);  _ff_init(u[3]);  _ff_init(v[0]);  _ff_init(v[1]);  _ff_init(v[2]);  }
#define _hecurve_set(u1,v1,u2,v2)		{_ff_set(u1[3],u2[3]);_ff_set(u1[2],u2[2]);_ff_set(u1[1],u2[1]);_ff_set(u1[0],u2[0]);_ff_set(v1[2],v2[2]);_ff_set(v1[1],v2[1]);_ff_set(v1[0],v2[0]);}
#define _hecurve_is_identity(u,v)		(_ff_zero((u)[1]) && _ff_zero((u)[2]) && _ff_zero((u)[3]))			// ignores v
#define _hecurve_set_identity(u,v) 		{_ff_set_one((u)[0]); _ff_set_zero((u)[1]);  _ff_set_zero((u)[2]); _ff_set_zero((u)[3]);  \
								  _ff_set_zero((v)[0]);  _ff_set_zero((v)[1]); _ff_set_zero((v)[2]);}
#define _hecurve_2tor(u,v)			(_ff_zero((v)[2]) && _ff_zero((v)[1]) && _ff_zero((v)[0]))
#endif

struct _hecurve_ctx_struct {
	ff_t invert;
	ff_t s0;
	ff_t s1;
	ff_t s2;
	ff_t r;
	ff_t inv1;
	int state;					// must be 0 on first call, 1 on second call
};
typedef struct _hecurve_ctx_struct hecurve_ctx_t;
	
// Jacobian coordinate rep for genus 1, represents the affine point s(x/z,y/z), in Mumford rep: u(t)=t-x, v(t)=y.
struct _hec_j_struct {
	ff_t x, y, z;		
};
typedef struct _hec_j_struct hec_j_t;
	
// reduced Chudnovsky Jacobian coordinate representation for genus 1
// Represents the affine point (x/z2,y/z3), or in Mumford rep: u(t)=t-x, v(t)=y.
struct _hec_rc_struct {
	ff_t x, y, z2, z3;				// note that we don't maintain z, this saves a field multiplication when composing elements
};
typedef struct _hec_rc_struct hec_jc_t;

#define _hecurve_init_ctx(ctx)		{_ff_init(ctx.invert);_ff_init(ctx.s0);_ff_init(ctx.s1);_ff_init(ctx.r);_ff_init(ctx.inv1);ctx.state=0;}
#define _hecurve_clear_ctx(ctx)	{_ff_clear(ctx.invert);_ff_clear(ctx.s0);_ff_clear(ctx.s1);_ff_clear(ctx.r);_ff_clear(ctx.inv1);ctx.state=0;}

// return 1 if equal, -1 if unequal inverses, 0 otherwise
#if HECURVE_GENUS == 1
// IMPORTANT - point comparison is only supported for points in affice coordinates!
// we could compare projective points by cross multiplying, but currently comparisons are only used during table lookup, in which case
// we expect to use a unique (affine) representation anyway.
static inline int hecurve_cmp (ff_t u1[3], ff_t v1[2], ff_t u2[3], ff_t v2[2])
{
	if ( _ff_zero(u1[1]) ) return ( _ff_zero(u2[1]) ? 1 : 0 );
	if ( _ff_zero(u2[1]) ) return 0;
#if ! HECURVE_FAST
	if ( ! _ff_one(u1[1]) || ! _ff_one(u2[1]) ) { puts ("Attempt to compare genus 1 points in non-affine coordinates!");  exit (0); }
#endif
	if ( _ff_equal(u1[0],u2[0]) ) {
		if ( _ff_equal (v1[0], v2[0]) ) return 1;
		if ( ff_is_negation(v1[0],v2[0]) ) return -1;
	}
	return 0;
}
#endif
#if HECURVE_GENUS == 2
static inline int hecurve_cmp (ff_t u1[3], ff_t v1[2], ff_t u2[3], ff_t v2[2])
{
	if ( _ff_equal(u1[0],u2[0]) && _ff_equal(u1[1],u2[1]) && _ff_equal(u1[2],u2[2]) ) {
		if ( _ff_equal (v1[0], v2[0]) && _ff_equal(v1[1], v2[1]) ) return 1;
		if ( ff_is_negation(v1[0],v2[0]) && ff_is_negation (v1[1],v2[1]) ) return -1;
	}
	return 0;
}
#endif
#if HECURVE_GENUS == 3
static inline int hecurve_cmp (ff_t u1[4], ff_t v1[3], ff_t u2[4], ff_t v2[3])
{
	if ( _ff_equal(u1[0],u2[0]) && _ff_equal(u1[1],u2[1]) && _ff_equal(u1[2],u2[2]) && _ff_equal(u1[3],u2[3]) ) {
		if ( _ff_equal (v1[0], v2[0]) && _ff_equal(v1[1], v2[1]) && _ff_equal(v1[2], v2[2]) ) return 1;
		if ( ff_is_negation(v1[0],v2[0]) && ff_is_negation (v1[1],v2[1]) && ff_is_negation (v1[2],v2[2]) ) return -1;
	}
	return 0;
}
#endif
static inline void hecurve_invert (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS])
{
#if HECURVE_GENUS == 3
	ff_negate(v[2]);
#endif
#if HECURVE_GENUS >= 2
	ff_negate(v[1]);
#endif
	ff_negate(v[0]);
}


#if HECURVE_GENUS == 1
#define hecurve_compose		hecurve_g1_compose
#define hecurve_square			hecurve_g1_square
#endif
#if HECURVE_GENUS == 2
#define hecurve_compose		hecurve_g2_compose
#define hecurve_square			hecurve_g2_square
#endif
#if HECURVE_GENUS == 3
#define hecurve_compose		hecurve_g3_compose
#define hecurve_square			hecurve_g3_square
#endif

void hecurve_setup (mpz_t p);
int hecurve_init_ctx (hecurve_ctx_t *ctx, ff_t f[HECURVE_DEGREE+1]);
void hecurve_clear_ctx (hecurve_ctx_t *ctx);

void hecurve_compose_cantor (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS],
						  ff_t u1[HECURVE_GENUS+1], ff_t v1[HECURVE_GENUS],
						  ff_t u2[HECURVE_GENUS+1], ff_t v2[HECURVE_GENUS], ff_t f[HECURVE_DEGREE+1]);

// ctx->state should be set by the caller to compose and square
//	0=inital call, return to invert		1=2nd call with inverted element		-1=initial call, don't return to invert
// note that inputs and outputs may overlap

int hecurve_g1_compose (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS], ff_t u1[HECURVE_GENUS+1],
				         ff_t v1[HECURVE_GENUS], ff_t u2[HECURVE_GENUS+1],
				         ff_t v2[HECURVE_GENUS+1], ff_t f[HECURVE_DEGREE+1], hecurve_ctx_t *ctx);
int hecurve_g1_square (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS],
				      ff_t u1[HECURVE_GENUS+1], ff_t v1[HECURVE_GENUS], ff_t f[HECURVE_DEGREE+1], hecurve_ctx_t *ctx);

int hecurve_g2_compose (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS], ff_t u1[HECURVE_GENUS+1],
				         ff_t v1[HECURVE_GENUS], ff_t u2[HECURVE_GENUS+1],
				         ff_t v2[HECURVE_GENUS+1], ff_t f[HECURVE_DEGREE+1], hecurve_ctx_t *ctx);
int hecurve_g2_square (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS],
				      ff_t u1[HECURVE_GENUS+1], ff_t v1[HECURVE_GENUS], ff_t f[HECURVE_DEGREE+1], hecurve_ctx_t *ctx);

int hecurve_g3_compose (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS], ff_t u1[HECURVE_GENUS+1],
				         ff_t v1[HECURVE_GENUS], ff_t u2[HECURVE_GENUS+1],
				         ff_t v2[HECURVE_GENUS+1], ff_t f[HECURVE_DEGREE+1], hecurve_ctx_t *ctx);
int hecurve_g3_square(ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS],
				     ff_t u1[HECURVE_GENUS+1], ff_t v1[HECURVE_GENUS], ff_t f[HECURVE_DEGREE+1], hecurve_ctx_t *ctx);

void hecurve_g1_exp_ui (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS], ff_t u1[HECURVE_GENUS+1], ff_t v1[HECURVE_GENUS], unsigned long e, ff_t f[HECURVE_DEGREE+1]);
void hecurve_g1_AJ_powers (hec_j_t p[], ff_t x0, ff_t y0, int k, ff_t f1);
int hecurve_g1_AJC_steps_1 (hec_jc_t s[], ff_t x0, ff_t y0, int n, ff_t f1);
int hecurve_g1_AJC_steps_2 (hec_jc_t s[], ff_t x0, ff_t y0, ff_t x1, ff_t y1, int n, ff_t f1);
long hecurve_g1_order (long *pd, ff_t f[4]);
int hecurve_g1_group_structure (long n[2], long N, long d, ff_t f[4]);
int hecurve_g1_test_exponent (long e, ff_t f[4]);


// The function below should are used in both genus 2 and 3 although the compression support functions (bits/unbits)
// are currently only supported in genus 2

int hecurve_verify (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS], ff_t f[HECURVE_DEGREE+1]);
void hecurve_random (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS], ff_t f[HECURVE_DEGREE+1]);
int  hecurve_random_point (ff_t *px, ff_t *py, ff_t f[HECURVE_DEGREE+1]);
void hecurve_invert (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS]);
unsigned hecurve_bits (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS], ff_t f[HECURVE_DEGREE+1]);
int hecurve_unbits (ff_t v[HECURVE_GENUS], ff_t u[HECURVE_GENUS+1], unsigned bits, ff_t f[HECURVE_DEGREE+1]);

void hecurve_copy (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS], ff_t u1[HECURVE_GENUS+1], ff_t v1[HECURVE_GENUS]);
void hecurve_print (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS]);
int hecurve_sprint (char *s, ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS]);

#endif
