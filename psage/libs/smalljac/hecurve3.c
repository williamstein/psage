#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "gmp.h"
#include "mpzutil.h"
#include "ffwrapper.h"
#include "ffpoly.h"
#include "hecurve.h"
#include "cstd.h"


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
    along with GenericGroupDemo.  If not, see <http://www.gnu.org/licenses/>.
*/

#if HECURVE_GENUS == 3

#define _dbg_print_uv(s,u,v)	if ( dbg_level >= DEBUG_LEVEL ) { printf (s);  hecurve_print(u,v); }

#define sub		_ff_sub
#define subf		_ff_subfrom
#define add 		_ff_add
#define addt 		_ff_addto
#define mult		_ff_mult
#define multx(a,b)	ff_mult(a,a,b)
#define sqr		_ff_square
#define neg		_ff_neg
#define dbl		_ff_x2

// For simplicity, we assume no initilization of field elements is required.  Debug code assume ff_t fits in an unsigned long.

int hecurve_g3_compose (ff_t u[4], ff_t v[3], ff_t u1[4], ff_t v1[3], ff_t u2[4], ff_t v2[3], ff_t f[8], hecurve_ctx_t *ctx)
{
	_ff_t_declare_reg d0, d1, d2, t0, w0, w1, w2, w3, w4, w5, w6, w7, r, x0, x1, x2;

	if ( ctx && ctx->state == 1 ) {
		// get invert result and restore saved variables
		_ff_set (w1, ctx->invert);		// inverted value of r*s1
		_ff_set (r, ctx->r);
		_ff_set (d0, ctx->s0);
		_ff_set (d1, ctx->s1);
		_ff_set (d2, ctx->s2);
///		dbg_printf ("Restored r = %lu, s2 = %lu, s1 = %lu, s0 = %lu, w1 = %lu\n", r, d2, d1, d0, w1);
		goto hecurve_g3_compose_inverted;
	}

#if ! HECURVE_FAST
	if ( dbg_level >= 2 ) {
		printf ("Compose_3_3 inputs\n");  hecurve_print(u1,v1);  hecurve_print(u2,v2);
	}
#endif
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u1, v1, f) )  { err_printf ("hecurve_g3_compose: invalid input\n");  exit (0); }
	if ( ! hecurve_verify (u2, v2, f) ) { err_printf ("hecurve_g3_compose: invalid input\n");  exit (0); }
#endif

	if ( ! _ff_one(u1[3]) || ! _ff_one(u2[3]) ) { hecurve_compose_cantor(u,v,u1,v1,u2,v2,f);  return 1; }

	// 1. Compute r = Resultant(u1,u2)
	sub(d0,u2[0],u1[0]);					// d0 = u20-u10
	sub(d1,u2[1],u1[1]);					// d1 = u21-u11
	sub(d2,u2[2],u1[2]);					// d2 = u22-u12
	mult(w1,u1[1],u2[0]);					// w1 = t3 = u11u20
	mult(t0,u1[0],u2[1]);					// t0 = t4 = u10u21
	subf(w1,t0);							// w1 = t3-t4
	mult(w2,u1[2],u2[0]);					// w2 = t5 = u12u20
	mult(t0,u1[0],u2[2]);					// t0 = t6 = u10u22
	subf(w2,t0);							// w2 = t5-t6
	mult(w3,u1[2],u2[1]);					// w3 = t1 = u12u21
	mult(t0,u1[1],u2[2]);					// t0 = t2 = u11u22
	subf(w3,t0);							// w3 = t1-t2
	mult(w4,d1,d0);						// w4 = t11
	sqr(w5,d1);							// w5 = t8
	mult(t0,d2,w1);						// t0 = t9
	sqr(w6,d0);							// w6 = t7
	subf(w6,t0);							// w6 = t7-t9
	mult(w7,d2,w2);						// w7 = t10 = (u22-u12)(t5-t6)
	subf(w7,w4);							// w7 = t10-t11
	add(t0,d0,w3);
	mult(r,t0,w6);
	sub(t0,w7,w4);
	multx(t0,w2);
	addt(r,t0);
	mult(t0,w5,w1);
	addt(r,t0);							// r = (d0+w3)w6+w2(w7-w4)+w5w1
///	dbg_printf ("r = %lu\n", r);
//	if ( _ff_zero(r) ) { hecurve_compose_cantor(u,v,u1,v1,u2,v2,f);  return 1; }		// handled below
	
	// 2. Compute inv = r/u1 mod u2
	add(t0,w3,d0);						// t0 = t1-t2+u20-u10
	mult(w1,t0,d2);
	subf(w1,w5);							// w1 = inv2 = (t1-t2+u20-u10)(u22-u12)-t8
	mult(w2,u2[2],w1);
	subf(w2,w7);							// w2 = inv1 = inv2u22 -t10+t11
	mult(t0,u2[2],w7);
	addt(t0,w6);							// t0 = u22(t10-t11)+t7-t9
	mult(w3,w1,u2[1]);
	subf(w3,t0);							// w3 = inv0 = inv2u21 - (u22(t10-t11)+t7-t9)
///	dbg_printf ("inv = %lux^2 + %lux + %lu\n", w1, w2, w3);
	
	// 3. Compute s' = rs = (v2-v1)inv (mod u2) - note that we don't use Karatusba here since field mults aren't much slower than add/sub
	sub(d0,v2[0],v1[0]);					// do = v20-v10
	sub(d1,v2[1],v1[1]);					// d1 = v21-v11
	sub(d2,v2[2],v1[2]);					// d2 = v22-v12
	mult(w5,d0,w3);						// w5 = t15 = r'0
	mult(w4,d2,w1);						// w4 = t17 = r'4
	mult(w6,d1,w3);
	mult(t0,d0,w2);
	addt(w6,t0);							// w6 = r'1 = d1w3+d0w2
	mult(w7,d1,w2);
	mult(t0,d2,w3);
	addt(w7,t0);
	mult(t0,d0,w1);
	addt(w7,t0);							// w7 = r'2 = d1w2+d2w3+d0w1
	mult(w3,d2,w2);
	mult(t0,d1,w1);
	addt(t0,w3);
	mult(w3,u2[2],w4);
	subf(w3,t0);							// w3 = t18 = u22w4 - (d2w2+d1w2)
	mult(w2,u2[0],w3);						// w2 = t15 = u20w3
	mult(w1,u2[1],w4);						// w1 = t16 = u21w4
	add(d0,w2,w5);						// d0 = s'0 = w2+w5
	add(d1,u2[0],u2[1]);
	sub(t0,w4,w3);
	multx(t0,d1);
	sub(d1,w6,t0);
	addt(d1,w1);
	subf(d1,w2);							// d1 = s'1 = w6 - (u20+u21)(w4-w3)+w1-w2
	mult(d2,u2[2],w3);
	subf(d2,w1);
	addt(d2,w7);							// d2 = s'2 = w7-w1+u22w3
///	dbg_printf ("s' = %lux^2 + %lux + %lu\n", d2, d1, d0);
	
	mult(t0,r,d2);
	if ( _ff_zero(t0) ) {
		int sts;
		// check for squaring and inverse cases before resorting to Cantor - caller really should use hecurve_square
		sts =  hecurve_cmp (u1,v1,u2,v2);
		if ( sts == 1 ) return hecurve_g3_square(u,v,u1,v1,f,ctx);
		if ( sts == -1 ) { _hecurve_set_identity(u,v);  return 1; }
		hecurve_compose_cantor(u,v,u1,v1,u2,v2,f);
		return 1;
	}
	if ( ctx && ! ctx->state ) {
		_ff_set (ctx->s0, d0);
		_ff_set (ctx->s1, d1);
		_ff_set (ctx->s2, d2);
		_ff_set (ctx->r, r);
		_ff_set (ctx->invert, t0);	// return to let caller invert
		ctx->state = 1;
///		dbg_printf("Returning to let caller invert %lu\n", t0);
		return 0;
	}
	_ff_invert(w1,t0);
	
hecurve_g3_compose_inverted:		// caller returns here after inverting
								// we assume w1, r, d0, d1, and d2 are set, all others free
	// 4. Compute s = s'/r and make monic
	mult(w2,r,w1);							// w2 = rw1 = 1/s2
	sqr(w7,d2);
	multx(w7,w1);
	ff_negate(w7);						// w7 = -s2/r
	mult(w6,r,w2);							// w6 = r/s2 (w4 in HECHECC)
	sqr(w5,w6);							// w5  = (r/s2)^2
	multx(d0,w2);							// d0 = s0'/s2 = s0
	multx(d1,w2);							// d1 = s1'/s2 = s1
/// 	dbg_printf ("s = x^2 + %lux + %lu\n", d1, d0);
	
	// 5. Compute z = su1
	mult (w0,d0,u1[0]);						// w0 = z0 = s0u10
	mult (w1,d1,u1[0]);
	mult (t0,d0,u1[1]);
	addt(w1,t0);							// w1 = z1 = s1u10+s0u11
	mult (w2,d0,u1[2]);
	mult(t0,d1,u1[1]);
	addt(w2,t0);
	addt(w2,u1[0]);						// w2 = z2 = s0u12+s1u11+u10
	mult(w3,d1,u1[2]);
	addt(w3,d0);
	addt(w3,u1[1]);						// w3 = z3 = s1u12+s0+u11
	add(w4,d1,u1[2]);						// w4 = z4 = u1[2],s1
///	dbg_printf ("z = x^5 + %lux^4 + %lux^3 + %lux^2 + %lux + %lu\n", w4, w3, w2, w1, w0);

	// 6. Compute u' = [s(z+2w6v1) - w5((u20-v1^2)]/u2		note h = 0, w6 is w4 in handbook
	add(r,w4,d1);
	subf(r,u2[2]);							// r = u'3 = z4+s1-u22
	mult(x2,d1,w4);
	mult(t0,r,u2[2]);
	subf(x2,t0);
	subf(x2,u2[1]);
	addt(x2,w3);
	addt(x2,d0);							// x2 = u'2 = -u22u'3 - u21 + z3 + s0 + s1z4
	mult(x1,w6,v1[2]);
	dbl(x1);
	mult(t0,d1,w3);
	addt(x1,t0);
	mult(t0,d0,w4);
	addt(x1,t0);
	addt(x1,w2);
	subf(x1,w5);
	mult(t0,x2,u2[2]);
	subf(x1,t0);
	mult(t0,r,u2[1]);
	subf(x1,t0);
	subf(x1,u2[0]);						// x1 = u'1 = 2w6v12 + s1z3 + s0z4 + z2 - w5 - u22u'2 - u21u'3 - u2[0]
	mult(x0,d1,v1[2]);
	addt(x0,v1[1]);
	multx(x0,w6);
	dbl(x0);
	mult(t0,d1,w2);
	addt(x0,t0);
	addt(x0,w1);
	mult(t0,d0,w3);
	addt(x0,t0);
	mult(t0,w5,u1[2]);
	addt(x0,t0);
	mult(t0,u2[2],x1);
	subf(x0,t0);
	mult(t0,u2[1],x2);
	subf(x0,t0);
	mult(t0,u2[0],r);
	subf(x0,t0);							// x0 = u'0 = 2w6(v11+s1v12) + s1z2 + z1 + s0z3 + w5u12 - u22u'1 - u21u'2 - u20u'3
///	dbg_printf ("u' = x^4 + %lux^3 + %lux^2 + %lux + %lu\n", r, x2, x1, x0);

	// 7. Compute v' = w7z - v1 mod u'			note h = 0 and w7 is -w3 in handbook
	sub(w5,r,w4);							// w5 = t1 = u'3 - z4
	mult(d0,x0,w5);
	addt(d0,w0);
	multx(d0,w7);
	subf(d0,v1[0]);						// d0 = v'0 = w7(u'0t1+z0) - v10
	mult(d1,x1,w5);
	subf(d1,x0);
	addt(d1,w1);
	multx(d1,w7);
	subf(d1,v1[1]);						// d1 =-v'1 = w7(u'1t1-u'0+z1) - v11
	mult(d2,x2,w5);
	subf(d2,x1);
	addt(d2,w2);
	multx(d2,w7);
	subf(d2,v1[2]);						// d2 = v'2 = w7(u'2t1-u'1+z2) - v12
	mult(w0,r,w5);
	subf(w0,x2);
	addt(w0,w3);
	multx(w0,w7);							// w0 = v'3 = w7(u'3t1 - u'2 + z3)
///	dbg_printf ("v' = %lux^3 + %lux^2 + %lux + %lu\n", w0, d2, d1, d0);

	// 8. reduce u' to u'' = (f-v'^2)/u'
	_ff_set_one(u[3]);
	sqr(t0,w0);
	addt(t0,r);
	neg(u[2],t0);							// u[2] = u''2 = -(u'3+v'3^2)		note f6 = 0
	mult(t0,r,u[2]);						// t0 = u'3u''2
	addt(t0,x2);							// t0 = u'3u''2+u'2
	mult(w7,d2,w0);
	dbl(w7);								// w7 = 2v'2v'3
	addt(t0,w7);
	neg(u[1],t0);
#if ! HECURVE_SPARSE
	addt(u[1],f[5]);						// u[1] = u''1 = -(u'2 + u''2u'3 + 2v'2v'3) + f5
#endif
	mult(t0,u[2],x2);						// t0 = u''2u'2
	addt(t0,x1);							// t0 = u''2u'2+u'1
	mult(w7,u[1],r);						// w7 = u''1u'3
	addt(t0,w7);
	mult(w7,d1,w0);
	dbl(w7);								// w7 = 2v'1v'3
	addt(t0,w7)
	sqr(w7,d2);							// w7 = v'2^2
	addt(t0,w7);
	neg(u[0],t0);
#if ! HECURVE_SPARSE
	addt(u[0],f[4]);						// u[0] = u''0 = -(u'1 + u''2u'2 + u''1u'3 + 2v'1v'3 + v'2^2) + f4
#endif
///	dbg_printf ("u = x^3 + %lux^2 + %lux + %lu\n", u[2], u[1], u[0]);

	// 9. compute v'' = -v' mod u3				note h = 0
	mult(t0,w0,u[2]);
	sub(v[2],t0,d2);						// v''2 = v'3u''2 - v'2
	mult(t0,w0,u[1]);
	sub(v[1],t0,d1);						// v'1 = v'3u''1 - v'1
	mult(t0,w0,u[0]);
	sub(v[0],t0,d0);						// v''0 = v'3u''0 - v'0
///	dbg_printf ("v = %lux^2 + %lux + %lu\n", v[2], v[1], v[0]);

#if HECURVE_VERIFY
	if ( ! hecurve_verify (u, v, f) ) {
		err_printf ("hecurve_compose_genus3 output verification failed for inputs:\n");
		err_printf ("      ");  hecurve_print(u1,v1);
		err_printf ("      ");  hecurve_print(u2,v2);
		err_printf ("note that inputs may have been modified if output overlapped.\n");
		exit (0);
	}
#endif

	return 1;	
}


int hecurve_g3_square (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS],
				 ff_t u1[HECURVE_GENUS+1], ff_t v1[HECURVE_GENUS], ff_t f[HECURVE_DEGREE+1], hecurve_ctx_t *ctx)
{
	_ff_t_declare_reg t0, w0, w1, w2, w3, w4, w5, w6, w7, r, z0, z1, z2, g0, g1, g2, g3,g4;

	if ( ctx && ctx->state == 1 ) {
		// get invert result and restore saved variables
		_ff_set (w1, ctx->invert);		// inverted value of r*s1
		_ff_set (r, ctx->r);
		_ff_set (z0, ctx->s0);
		_ff_set (z1, ctx->s1);
		_ff_set (z2, ctx->s2);
///		dbg_printf ("Restored r = %lu, s2 = %lu, s1 = %lu, s0 = %lu, w1 = %lu\n", r, z2, z1, z0, w1);
		goto hecurve_g3_square_inverted;
	}

#if ! HECURVE_FAST
	if ( dbg_level >= 2 ) {
		printf ("hecurve_g3_square input\n");  hecurve_print(u1,v1);
	}
#endif
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u1, v1, f) )  { err_printf ("hecurve_g3_square: invalid input\n");  exit (0); }
#endif

	if ( ! _ff_one(u1[3]) ) { hecurve_compose_cantor(u,v,u1,v1,u1,v1,f);  return 1; }

	// 1. Compute r = Resultant(u,2v)
	add(z0,v1[0],v1[0]);					// z0 = h'0 = 2v0
	add(z1,v1[1],v1[1]);					// z1 = h'1 = 2v1
	add(z2,v1[2],v1[2]);					// z2 = h'2 = 2v2
	mult(t0,u1[1],z2);
	mult(w1,u1[2],z1);
	subf(w1,t0);							// w1 = t1-t2 = u2z1-u1z2
	mult(t0,u1[0],z1);
	mult(w2,u1[1],z0);
	subf(w2,t0);							// w2 = t3-t4 = u1z0-u0z1
	mult(t0,u1[0],z2);
	mult(w3,u1[2],z0);
	subf(w3,t0);							// w3 = t5-t6 = u2z0-u0z2
	mult(t0,z2,w2);
	sqr(w4,z0);
	subf(w4,t0);							// w4 = t7-t9 = z0^2 = z2w2
	mult(w6,z1,z0);						// w6 = t11 = z1z0
	mult(w5,z2,w3);
	subf(w5,w6);							// w5 = t10-t11 = z2w3 - z1z0
	sqr(w0,z1);							// w0 = t8 = z1^2
	add(r,z0,w1);
	multx(r,w4);
	sub(t0,w5,w6);
	multx(t0,w3);
	addt(r,t0);
	mult(t0,w0,w2);
	addt(r,t0);							// r = (z0+w1)w4 + w3(w5+z1z0) + w0w2
///	dbg_printf ("r = %lu\n", r);

	// 2. Compute inv = r/2v mod u
	add(t0,w1,z0);
	multx(t0,z2);
	sub(w2,w0,t0);						// w2 = inv2 = w0-z2(w1+z0)
	mult(w1,w2,u1[2]);
	addt(w1,w5);							// w1 = inv1 = inv2u2 + w5
	mult(w0,w2,u1[1]);
	mult(t0,w5,u1[2]);
	addt(w0,t0);
	addt(w0,w4);							// w0 = inv0 = inv2u1+u2w5+w4
///	dbg_printf ("inv = %lux^2 + %lux + %lu\n", w2, w1, w0);

	// 3. Compute z = (f-v^2)/u mod u
	mult(w5,u1[1],u1[2]);					// w5 = -t13 = u1u2				note that f6 = 0 so z'3 = -u2
	sqr(w4,u1[2]);
	subf(w4,u1[1]);
#if ! HECURVE_SPARSE
	addt(w4,f[5]);							// w4 = z'2 = u2^2-u1+f5		note h and f6 = 0
#endif
	sqr(w3,v1[2]);
	subf(w3,w5);
	mult(t0,w4,u1[2]);
	addt(w3,t0);
	addt(w3,u1[0]);
	ff_negate(w3);
#if ! HECURVE_SPARSE
	addt(w3,f[4]);							// w3 = z'1 = -(v2^2+w5+w4u2+u0) + f4
#endif
///	dbg_printf ("z' = x^4 - %lux^3 + %lux^2 + %lux + ?, w5 = -t13 = u1u2 = %lu\n", u1[2], w4, w3, w5);
	add(z2,u1[2],u1[2]);
	addt(z2,u1[2]);
	multx(z2,u1[2]);
	add(t0,u1[1],u1[1]);
	subf(z2,t0);
#if ! HECURVE_SPARSE
	addt(z2,f[5]);							// z2 = 3u2^2-2u1+f5
#endif
	add(z1,w5,w5);
	addt(z1,w3);
	subf(z1,u1[0]);						// z1 = 2u1u2+w3-u0
	mult(t0,v1[1],v1[2]);
	dbl(t0);
	mult(z0,w4,u1[1]);
	addt(t0,z0);
	mult(z0,w3,u1[2]);
	addt(t0,z0);
	add(z0,u1[2],u1[2]);
	addt(z0,u1[2]);
	multx(z0,u1[0]);
	subf(z0,t0);
#if ! HECURVE_SPARSE
	addt(z0,f[3]);							// z0 = 3u0u2 - (2v1v2+z'2u1+z'1u2) + f3
#endif
///	dbg_printf ("z = %lux^2 + %lux + %lu\n", z2, z1, z0);

	// 4. Compute s' = (z inv) mod u
	mult(t0,w1,z0);
	mult(w3,w0,z1);
	addt(w3,t0);							// w3 = r'1 = w0z1+w1z0
	mult(w4,w1,z1);
	mult(t0,w0,z2);
	addt(w4,t0);
	mult(t0,w2,z0);
	addt(w4,t0);							// w4 = r'2 = w1z1+w0z2+w2z0
	mult(t0,w2,z1);
	mult(w5,w1,z2);
	addt(w5,t0);							// w5 = r'3 = w1z2+w2z1
	multx(z2,w2);							// z2 = r'4
	mult(w6,z2,u1[2]);
	subf(w6,w5);							// w6 = t18 = u2z2-w5
	mult(w1,w6,u1[0]);						// w1 = t15
	mult(w2,z2,u1[1]);						// w2 = t16
	multx(z0,w0);
	addt(z0,w1);							// z0 = s'0 = z0w0+w1
	sub(t0,w6,z2);
	add(z1,u1[0],u1[1]);
	multx(z1,t0);
	addt(z1,w3);
	addt(z1,w2);
	subf(z1,w1);							// z1 = s'1 = w3+(u1+u0)(w6-z2)+w2-w1
	mult(z2,w6,u1[2]);
	addt(z2,w4);
	subf(z2,w2);							// z2 = s'2 = w4-w2+u2w6
///	dbg_printf ("s' = %lux^2 + %lux + %lu\n", z2, z1, z0);

	mult(t0,r,z2);
	if ( _ff_zero(t0) ) { hecurve_compose_cantor(u,v,u1,v1,u1,v1,f);  return 1; }
	if ( ctx && ! ctx->state ) {
		_ff_set (ctx->s0, z0);
		_ff_set (ctx->s1, z1);
		_ff_set (ctx->s2, z2);
		_ff_set (ctx->r, r);
		_ff_set (ctx->invert, t0);	// return to let caller invert
		ctx->state = 1;
///		dbg_printf("Returning to let caller invert %lu\n", t0);
		return 0;
	}
	_ff_invert(w1,t0);
	
hecurve_g3_square_inverted:		// caller returns here after inverting
							// we assume w1, r, z0, z1, and z2 are set, all others free
	// 5. Compute s = s'/r and make s monic
	mult(w2,w1,r);							// w2 = w1r
	sqr(w3,z2);
	multx(w3,w1);
	ff_negate(w3);							// w3 = -w1s'2^2 (note sign change)
	mult(w4,w2,r);							// w4 = r/s'2
	sqr(w5,w4);
	multx(z0,w2);							// z0 = s0 = w2s'0
	multx(z1,w2);							// z1 = s1 = w2s'1
///	dbg_printf ("s = x^2 + %lux + %lu, w3 = %lu (-), w4 = %lu, w5 = %lu\n", z1, z0, w3, w4, w5);
	
	// 6. Compute G=su
	mult(g0,z0,u1[0]);						// g0 = z0u1
	mult(t0,z1,u1[0]);
	mult(g1,z0,u1[1]);
	addt(g1,t0);							// g1 = z1u0+z0u1
	mult(g2,z0,u1[2]);
	mult(t0,z1,u1[1]);
	addt(g2,t0);
	addt(g2,u1[0]);						// g2 = z0u2+z1u1+u0
	mult(g3,z1,u1[2]);
	addt(g3,z0);
	addt(g3,u1[1]);						// g3 = z1u2+z0+u1
	add(g4,z1,u1[2]);						// g4 = z1+u2
///	dbg_printf ("G = su = x^5 + %lux^4 + %lux^3 + %lux^2 + %lux + %lu\n", g4, g3, g2, g1, g0);

	// 7. Compute u' = [(G+w4v)^2 + w5f]/u^2
	add(w6,z1,z1);							// w6 = u'3  = 2z1
	sqr(w2,z1);
	add(t0,z0,z0);
	addt(w2,t0);							// w2 = u'2 = z1^2+2z0
	mult(t0,w4,v1[2]);
	mult(w1,z0,z1);
	addt(w1,t0);
	dbl(w1);
	subf(w1,w5);							// w1 = u'1 = 2(z0z1+w4v2)-w5
	sub(t0,z1,u1[2]);
	multx(t0,v1[2]);
	addt(t0,v1[1]);
	multx(t0,w4);
	multx(w5,u1[2]);
	addt(t0,w5);
	dbl(t0);
	sqr(w0,z0);
	addt(w0,t0);							// w0 = u'0 = z0^2 + 2(w4(v1+v2(z1-u2))+w5u2)
///	dbg_printf ("u' = x^4 + %lux^3 + %lux^2 + %lux + %lu\n", w6, w2, w1, w0);

	// 8. Compute v' = -(Gw3 +v) mod u'
	sub(w4,w6,g4);						// w4 = t1 = w6-g4
	subf(g3,w2);
	mult(t0,w4,w6);
	addt(g3,t0);
	multx(g3,w3);							// g3 = v'3 = w3(g3-w2+w4w6)	(note w3 sign changed from HECHECC)
	subf(g2,w1);
	mult(t0,w4,w2);
	addt(g2,t0);
	multx(g2,w3);
	subf(g2,v1[2]);						// g2 = v'2 = w3(g2-w1+w4w2) - v2
	subf(g1,w0);
	mult(t0,w4,w1);
	addt(g1,t0);
	multx(g1,w3);
	subf(g1,v1[1]);						// g1 = v'1 = w3(g1-w0+w4w1) - v1
	mult(t0,w4,w0);
	addt(g0,t0);
	multx(g0,w3);
	subf(g0,v1[0]);						// g0 = v'0 = w3(g0+w4w0) - v0
///	dbg_printf ("v' = %lux^3 + %lux^2 + %lux + %lu\n", g3, g2, g1, g0);

	// 9. Compute u'' = (f-v'^2)/u'
	_ff_set_one(u[3]);
	sqr(t0,g3);
	addt(t0,w6);
	neg(u[2],t0);							// u[2] = u''2 = -(w6+g3^2)
	mult(t0,g2,g3);
	dbl(t0);
	mult(z0,w6,u[2]);
	addt(z0,t0);
	addt(z0,w2);
	neg(u[1],z0);
#if ! HECURVE_SPARSE
	addt(u[1],f[5]);						// u[1] = u''1 = -(w2+w6u[2]+2g1g2) + f5
#endif
	sqr(t0,g2);
	mult(z0,g1,g3);
	dbl(z0);
	addt(z0,t0);
	mult(t0,w6,u[1]);
	addt(z0,t0);
	mult(t0,w2,u[2]);
	addt(z0,t0);
	addt(z0,w1);
	neg(u[0],z0);
#if ! HECURVE_SPARSE
	addt(u[0],f[4]);						// u[0] = u''0 = -(w1+w2u[2]+w6u[1] + 2g1g3+g2^2) + f4
#endif
///	dbg_printf ("u = x^3 + %lux^2 + %lux + %lu\n", u[2], u[1], u[0]);

	// 10. Compute v'' = -v' mod u''
	mult(t0,g3,u[2]);
	sub(v[2],t0,g2);						// v[2] = v''2 = g2+g3u[2]
	mult(t0,g3,u[1]);
	sub(v[1],t0,g1);						// v[1] = v''1 = g1+g3u[1]
	mult(t0,g3,u[0]);
	sub(v[0],t0,g0);						// v[0] = v'0 = g0+g3u[0]
///	dbg_printf ("v = %lux^2 + %lux + %lu\n", v[2], v[1], v[0]);
	
	
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u, v, f) ) {
		err_printf ("hecurve_g3_square output verification failed for input:\n");
		err_printf ("      ");  hecurve_print(u1,v1);
		err_printf ("note that input may have been modified if output overlapped.\n");
		exit (0);
	}
#endif

	return 1;	
}

#endif
