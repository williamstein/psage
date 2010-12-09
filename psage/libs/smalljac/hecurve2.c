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
    along with smalljac.  If not, see <http://www.gnu.org/licenses/>.
*/

#if HECURVE_GENUS == 2

static char buf[4096];

#define _deg_u(u)			(_ff_nonzero(u[2])?2:(_ff_nonzero(u[1])?1:0))
#define _deg_v(v)			(_ff_nonzero(v[1])?1:(_ff_nonzero(v[0])?0:-1))
#define _dbg_print_uv(s,u,v)	if ( dbg_level >= DEBUG_LEVEL ) { printf (s);  hecurve_print(u,v); }

void hecurve_g2_square_1 (ff_t u[3], ff_t v[2], ff_t u0, ff_t v0, ff_t f[6]);
void hecurve_g2_square_2_r0 (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t f[6]);
void hecurve_g2_compose_1_2 (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t u2[3], ff_t v2[2], ff_t f[6]);
int hecurve_g2_compose_special (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t u2[3], ff_t v2[3], ff_t f[6]);

/*
	TODO: Clean up code to not use quite so many variables.  Useful for debugging at the moment.
*/
// note that (u,v) may be equal to (u1,v1) or (u2,v2) but we assume overlap is aligned.
int hecurve_g2_compose (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t u2[3], ff_t v2[2], ff_t f[6], hecurve_ctx_t *ctx)
{
	_ff_t_declare_reg r, inv0, inv1, w0, w1, w2, w3, w4, w5, s0, s1, t1, L0, L1, L2;

	if ( ctx && ctx->state == 1 ) {
		// get invert result and restore saved variables
		_ff_set (w1, ctx->invert);		// inverted value of r*s1
		_ff_set (r, ctx->r);
		_ff_set (s0, ctx->s0);
		_ff_set (s1, ctx->s1);
		_ff_set (inv1, ctx->inv1);
///		dbg_printf ("Restored r = %Zd, s1 = %Zd, s0 = %Zd, inv1 = %Zd, w1 = %Zd\n",_ff_wrap_mpz(r,0), _ff_wrap_mpz(s1,0), _ff_wrap_mpz(s0,0), _ff_wrap_mpz(inv1,0), _ff_wrap_mpz(w1,0));
		goto hecurve_g2_compose_2_2_inverted;
	}
	if ( _ff_nonzero(f[4]) ) { hecurve_compose_cantor (u, v, u1, v1, u2, v2, f);  return 1; }		// revert to Cantor if f4 is non-zero - should only happen over F_5
	if ( _ff_zero(u1[2]) || _ff_zero(u2[2]) ) { hecurve_g2_compose_special (u, v, u1, v1, u2, v2, f);  return 1; }
#if ! HECURVE_FAST
	dbg_printf ("Compose_2_2: (%Zdx^2+%Zdx+%Zd, %Zdx+%Zd) * (%Zdx^2+%Zdx+%Zd, %Zdx+%Zd)\n",
			    _ff_wrap_mpz(u1[2],0), _ff_wrap_mpz(u1[1],1), _ff_wrap_mpz(u1[0],2), _ff_wrap_mpz(v1[1],3), _ff_wrap_mpz(v1[0],4),
			    _ff_wrap_mpz(u2[2],5), _ff_wrap_mpz(u2[1],6), _ff_wrap_mpz(u2[0],7), _ff_wrap_mpz(v2[1],8), _ff_wrap_mpz(v2[0],9));
	if ( !_ff_one(u1[2]) || !_ff_one(u2[2]) ) { err_printf ("u1 or u2 not monic deg 2 in hecurve_compose_2_2!\n");  exit (0); }
#endif
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u1, v1, f) )  { err_printf ("hecurve_compose: invalid input\n");  exit (0); }
	if ( ! hecurve_verify (u2, v2, f) ) { err_printf ("hecurve_compose: invalid input\n");  exit (0); }
#endif
		
	// 1. compute r = Resultant(u1,u2) and save z1=inv1 and z3=inv0
	_ff_sub(inv1,u1[1],u2[1]);			// reduce inv1
	_ff_qsub(w0,u2[0],u1[0]);			// w0 = z2 (unreduced)
	_ff_multadd(inv0,u1[1],inv1,w0);		// z3 = inv0 = u1[1]*inv1+w0 = u1[1]*z1+z2
	_ff_mult(w1,inv1,inv1); 
	ff_mult(w1,w1,u1[0]);
	_ff_multadd(r,w0,inv0,w1);
///	dbg_printf ("r = %Zd\n",_ff_wrap_mpz(r,0));
//	if ( _ff_zero(r) ) { hecurve_compose_special (u, v, u1, v1, u2, v2, f);  return 1; }		this is checked below
	
	// 2. compute almost inverse of u2 mod u1, i.e. inv = (r/u2) mod u1
	// nothing to do, simply note that inv = inv1x + inv0, i.e. inv1 = z1 and inv0 = z3
///	dbg_printf ("inv = %Zdx + %Zd\n",_ff_wrap_mpz(inv1,0), _ff_wrap_mpz(inv0,1));
	
	// 3. compute s' = rs = ((v1-v2)inv) mod u1
	_ff_qsub(w0,v1[0],v2[0]);
	_ff_mult(w2,inv0,w0);
	_ff_qsub(w1,v1[1],v2[1]);
	_ff_mult(w3,inv1,w1); 
	_ff_incmult(t1,u1[1],w3);
	_ff_qadd (w4,inv0,inv1);
	_ff_qadd (w5,w0,w1);
	_ff_qnegsum(w0,w2,t1);
	_ff_multadd(s1,w4,w5,w0);
	_ff_mult(t1,u1[0],w3);
	_ff_qsub(s0,w2,t1);
///	dbg_printf ("w0 = %Zd, w1 = %Zd, w2 = %Zd, w3 = %Zd, w4 = %Zd, w5 = %Zd, s\' = rs = %Zdx + %Zd\n",
///			  _ff_wrap_mpz(w0,0), _ff_wrap_mpz(w1,1), _ff_wrap_mpz(w2,2), _ff_wrap_mpz(w3,3), _ff_wrap_mpz(w4,4), _ff_wrap_mpz(w5,5), _ff_wrap_mpz(s1,6), _ff_wrap_mpz(s0,7));
	
	// Note that we need to save w3, w4, and w5, but w0, w1, and w2 are free
//	if ( _ff_zero(s1) ) goto hecurve_compose_s1_zero;								this is checked below

	// 4. compute s'' = x+s0/s1 = x + s'0/s'1 and s1
	_ff_mult(w2,r,s1);			// assume r and s1 both non-zero - this is the usual case
	if ( _ff_zero(w2) ) {
		if ( _ff_nonzero(r) ) goto hecurve_g2_compose_s1_zero;
		hecurve_g2_compose_special (u, v, u1, v1, u2, v2, f);  return 1;
	}
	if ( ctx && ! ctx->state ) {
		_ff_set (ctx->s0, s0);
		_ff_set (ctx->s1, s1);
		_ff_set (ctx->r, r);
		_ff_set(ctx->inv1, inv1);
		_ff_set (ctx->invert, w2);	// return to let caller invert
		ctx->state = 1;
		return 0;
	}
	_ff_invert(w1,w2);
	
hecurve_g2_compose_2_2_inverted:	// caller returns here after inverting
								// we assume w1, r, s0, s1, and inv1 are set, all others free, s0 is unreduced, the rest are reduced
	_ff_mult(w2,w1,r);
	_ff_square(w3,s1); 
	ff_mult(w3,w3,w1);
	_ff_mult(w4,r,w2); 
	_ff_square(w5,w4);
	_ff_mult(s0,s0,w2);
///	dbg_printf ("s\'\' = x + %Zd, s1 = %Zd\n",_ff_wrap_mpz(s0,0), _ff_wrap_mpz(w3,1));
	
	// 5. compute L = s''0*u2 = x^3 + L2*x^2 + L1*x + L0
	_ff_qadd(L2,u2[1],s0);		// ok to leave unreduced
	_ff_multadd(L1,u2[1],s0,u2[0]);
	_ff_mult(L0,u2[0],s0);
///	dbg_printf ("L = x^3 + %Zdx^2 + %Zdx + %Zd\n",_ff_wrap_mpz(L2,0), _ff_wrap_mpz(L1,1), _ff_wrap_mpz(L0,2));
	
	// 6. compute u = (s(L+2v2)-t)/u1 = x^2 + u[1]*x + u[0]		note h = 0
	_ff_qsub(w0,s0,u1[1]);
	_ff_qsub(w1,s0,inv1);
	_ff_mult(w2,w0,w1);			// w2 = (s0-u1[1])(s0-inv1)
	_ff_qadd(w1,u2[1],u2[1]);
	_ff_qaddto(w1,inv1);
	ff_mult(w1,w1,w5);
	_ff_qaddto (w2,w1);
	_ff_qneg(w1,u1[0]);
	_ff_qaddto(w2,w1);
	_ff_addto(w2,L1);
	_ff_qadd(w0,v2[1],v2[1]);
	_ff_multadd(u[0],w0,w4,w2);				// potentially destroys u1[0] and u2[0]
	// no mults to work with so we are better off reducing as we go to compute u[1]
	_ff_add(w0,s0,s0);
	_ff_subfrom(w0,inv1);
	_ff_sub(u[1],w0,w5);
	_ff_set_one(u[2]);
///	dbg_printf ("u = %Zdx^2 + %Zdx + %Zd\n",_ff_wrap_mpz(u[2],0), _ff_wrap_mpz(u[1],1), _ff_wrap_mpz(u[0],2));
	
	// 7. compute v = (-(L+v2)mod u = v[1]*x + v[0]					note h = 0
	_ff_qsub(w1,L2,u[1]);
	_ff_qsub(w0,u[0],L1);
	_ff_multadd(w2,w1,u[1],w0);
	_ff_qneg(w0,v2[1]);
	_ff_multadd (v[1],w2,w3,w0);
	_ff_qneg(w0,L0);
	_ff_multadd(w2,w1,u[0],w0);
	_ff_qneg(w0,v2[0]);
	_ff_multadd(v[0],w2,w3,w0);
///	dbg_printf ("v = %Zdx + %Zd\n",_ff_wrap_mpz(v[1],0), _ff_wrap_mpz(v[0],1));
#if ! HECURVE_FAST
	_dbg_print_uv ("Compose_2 Result: ", u, v);
#endif
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u, v, f) ) {
		err_printf ("hecurve_compose output verification failed for inputs:\n");
		err_printf ("      ");  hecurve_print(u1,v1);
		err_printf ("      ");  hecurve_print(u2,v2);
		err_printf ("note that inputs have been modified if output overlapped.\n");
		exit (0);
	}
#endif
	return 1;

	// Special case is not optimized
hecurve_g2_compose_s1_zero:
///		dbg_printf ("Special case, s\'\'1 == 0\n");

	// 4'. compute s
	_ff_invert(inv0,r);
	ff_mult(s0,s0,inv0); 
///	dbg_printf ("s = %Zd\n",_ff_wrap_mpz(s0,0));
	
	// Part of 6' moved here incase u2[0] gets overwritten below
	_ff_mult(w2,u2[0],s0);
	_ff_addto(w2,v2[0]);

	// 5'. compute u = (t-s(L+2v2))/u1 = x + u[0]		note h = 0
	_ff_square(w1,s0);
	_ff_addto(w1,u2[1]);
	_ff_addto(w1,u1[1]);
	_ff_neg(u[0],w1);				// this potentially destroys u1[0] and u2[0], we assume u2[1] is ok
///	dbg_printf ("u = x + %Zd\n",_ff_wrap_mpz(u[0],0));
	
	// 6. compute v = -(L+v2) mod u = v[0]				note h = 0
	_ff_sub (w1,u2[1],u[0]);				// BUG FIX - handbook has w1 = s0(u21+u[0]) (+ SHOULD BE -)
	ff_mult(w1,w1,s0);
	_ff_addto(w1,v2[1]);
	ff_mult(w1,w1,u[0]);
	_ff_sub(v[0],w1,w2);
	_ff_set_zero(v[1]);

	_ff_set_zero(u[2]);  _ff_set_one(u[1]);		// set these last to avoid overlap issues

///	dbg_printf ("v = %Zd\n",_ff_wrap_mpz(v[0],0));
#if ! HECURVE_FAST
	_dbg_print_uv ("Compose_2 Result: ", u, v);
#endif
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u, v, f) ) {
		err_printf ("hecurve_compose output verification failed for inputs:\n");
		err_printf ("      ");  hecurve_print(u1,v1);
		err_printf ("      ");  hecurve_print(u2,v2);
		err_printf ("note that inputs have been modified if output overlapped.\n");
		exit (0);
	}
#endif
	return 1;
}


int hecurve_g2_square (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t f[6], hecurve_ctx_t *ctx)
{
	_ff_t_declare_reg w0, w1, w2, w3, w4, w5, inv0, inv1, L0, L1, L2, r, s0, s1, t0, t1;
		
#if ! HECURVE_FAST
	_dbg_print_uv ("Square_2: ", u, v);
#endif	
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u1, v1, f) )  { err_printf ("hecurve_square: invalid input\n");  exit (0); }
#endif
	
	if ( ctx && ctx->state == 1 ) {
		// get invert result and restore saved variables
		_ff_set (w1, ctx->invert);		// inverted value of r*s1
		_ff_set (r, ctx->r);
		_ff_set (s0, ctx->s0);
		_ff_set (s1, ctx->s1);
///		dbg_printf ("Restored r = %Zd, s1 = %Zd, s0 = %Zd, w1 = %Zd\n",_ff_wrap_mpz(r,0), _ff_wrap_mpz(s1,1), _ff_wrap_mpz(s0,2), _ff_wrap_mpz(w1,3));
		goto hecurve_g2_square_2_2_inverted;
	}
	if ( _ff_nonzero(f[4]) ) { hecurve_compose_cantor (u, v, u1, v1, u1, v1, f);  return 1; }		// revert to Cantor if f4 is non-zero - should only happen over F_5
	if ( _ff_zero(u1[2]) ) { hecurve_g2_square_1 (u, v, u1[0], v1[0], f);  return 1; }

	// 1. compute v' = 2v mod u = v1*x + v0			note h = 0, we use w5=~v1, inv0=~v0
	_ff_add(w5,v1[1],v1[1]);				// w5 gets negated below, must be reduced
	_ff_qadd(inv0,v1[0],v1[0]);
///	dbg_printf ("2v = %Zdx + %Zd\n",_ff_wrap_mpz(w5,0), _ff_wrap_mpz(inv0,1));
	
	// 2. compute resultant r = Resultant(2v,u)
	_ff_square(w0,v1[1]);
	_ff_square(w1,u1[1]);
	_ff_qadd(w2,w0,w0); 
	_ff_qadd(w2,w2,w2);
	_ff_mult(w3,u1[1],w5); 
	_ff_qsub(w4,inv0,w3);
	ff_mult(w4,w4,inv0);
	_ff_multadd(r,u1[0],w2,w4);
//	_ff_addto(r,w4);
///	dbg_printf ("w0 = %Zd, w1 = %Zd, w2 = %Zd, w3 = %Zd, w4 = %Zd r = res(2v,u) = %Zd\n",
///			  _ff_wrap_mpz(w0,0), _ff_wrap_mpz(w1,1), _ff_wrap_mpz(w2,2), _ff_wrap_mpz(w3,3), _ff_wrap_mpz(w4,4), _ff_wrap_mpz(r,5));
//	if ( _ff_zero(r) ) { hecurve_square_2_r0 (u, v, u1, v1, f);  return 1; }		handled below
	
	// 3. compute almost inverse invv = r*inv
	_ff_neg(inv1,w5);
	_ff_subfrom(inv0,w3);
///	dbg_printf ("inv = %Zdx + %Zd\n", _ff_wrap_mpz(inv1,0), _ff_wrap_mpz(inv0,1));

	// 4. compute t = ((f-v^2)/u) mod u = tt1*x+tt0	note h = 0
	// this code could be improved further when f[3] == 0
#if HECURVE_SPARSE
	_ff_set (w3,w1);	// this could be further optimized
#else
	_ff_add(w3,f[3],w1);
#endif
	_ff_add(w4,u1[0],u1[0]);
	_ff_qadd(t1,w1,w1);
	_ff_qaddto(t1,w3);
	_ff_qsubfrom(t1,w4);
	_ff_qadd(t0,w4,w4);
	_ff_qsubfrom(t0,w3);
	_ff_qneg(w5,w0);
	_ff_multadd(t0,t0,u1[1],w5);
#if ! HECURVE_SPARSE
	_ff_addto(t0,f[2]);
#endif
///	dbg_printf ("t = %Zdx + %Zd\n", _ff_wrap_mpz(t1,0), _ff_wrap_mpz(t0,1));  

	// 5. compute ss = (tt*invv) mod u
	_ff_mult(w0,t0,inv0);
	_ff_mult(w1,t1,inv1);
	_ff_add(s1,inv0,inv1);		// need s1 reduced to keep mult from getting too big
	_ff_add(w2,t0,t1);
	_ff_qneg(w5,w0);
	_ff_multadd(s1,s1,w2,w5);
	_ff_incmult(w2,u1[1],w1);
	_ff_qsubfrom(s1,w2);
	_ff_mult(w5,w1,u1[0]);
	_ff_qsub(s0,w0,w5);
///	dbg_printf ("s' = %Zdx + %Zd\n",_ff_wrap_mpz(s1,0), _ff_wrap_mpz(s0,1));
//	if ( _ff_zero(s1) ) goto hecurve_square_s1_zero;		handled below

	// 6. compute s'' = x + s0/s1 and s1
	_ff_mult(w0,r,s1);		// assume non-zero r and s1 is usual case
	if ( _ff_zero(w0) ) {
		if ( _ff_nonzero(r) ) goto hecurve_g2_square_s1_zero;
		hecurve_g2_square_2_r0 (u, v, u1, v1, f);
		return 1;
	}
	if ( ctx && ! ctx->state  ) {
		_ff_set (ctx->s0, s0);
		_ff_set (ctx->s1, s1);
		_ff_set (ctx->r, r);
		_ff_set (ctx->invert, w0);
		ctx->state = 1; 								// return to let caller invert
		return 0;
	}
	_ff_invert(w1,w0);

hecurve_g2_square_2_2_inverted:	// caller returns here after inverting
							// we assume w0, s0, s1, and r are set and everything else is free
	_ff_mult(w2,w1,r); 
	_ff_square(w3,s1);
	ff_mult(w3,w3,w1);
	_ff_mult(w4,w2,r);
	_ff_square(w5,w4);
	ff_mult(s0,s0,w2);
///	dbg_printf ("s\'\'0 = %Zd, w0 = %Zd, w1 = %Zd, w2 = %Zd, w3 = %Zd, w4 = %Zd, w5 = %Zd\n",
///			  _ff_wrap_mpz(s0,0), _ff_wrap_mpz(w0,1), _ff_wrap_mpz(w1,2), _ff_wrap_mpz(w2,3), _ff_wrap_mpz(w3,4), _ff_wrap_mpz(w4,5), _ff_wrap_mpz(w5,6));
	// note that w3, w4, and w5, need to be saved, but w0, w1 and w2 are free
		
	// 7. compute LL = sss*u = x^3 + LL2*x^2 + LL1*x + LL0
	_ff_qadd(L2,s0,u1[1]);
	_ff_multadd(L1,s0,u1[1],u1[0]);
	_ff_mult(L0,s0,u1[0]);
///	dbg_printf ("L' = s''u = x^3 + %Zdx^2 + %Zdx + %Zd\n",_ff_wrap_mpz(L2,0), _ff_wrap_mpz(L1,1), _ff_wrap_mpz(L0,2));
	
	// 8. compute u = s^2 + 2vs/u + (v^2-f)/u1^2		note h = 0
	_ff_square(w0,s0);
	_ff_qadd(w1,v1[1],v1[1]);
	ff_mult(w1,w1,w4);
	_ff_qaddto(w0,w1);
	_ff_qadd(w1,u1[1],u1[1]);
	_ff_multadd(u[0],w1,w5,w0);						// potentially destroys u1[0] and u2[0]
	_ff_add(w0,s0,s0);								// no mult to absorb so stay reduced
	_ff_sub(u[1],w0,w5);							// potentially destroys u1[1] and u2[1]
	_ff_set_one(u[2]);
///	dbg_printf ("u = x^2 + %Zdx + %Zd\n",_ff_wrap_mpz(u[1],0), _ff_wrap_mpz(u[0],1));
	
	// 9. compute v = -(L+v1) mod u = v[1]x + v[0]			note h = 0
	_ff_qsub(w1,L2,u[1]);
	_ff_qsub (w0,u[0],L1);
	_ff_multadd(w2,w1,u[1], w0);
	_ff_qneg(w0,v1[1]);		// v1[1] could overlap v[1]
	_ff_multadd(v[1],w2,w3,w0);
	_ff_qneg(w0,L0);
	_ff_multadd(w2,w1,u[0],w0);
	_ff_qneg(w0,v1[0]);		// v1[0] could overlap v[0]
	_ff_multadd(v[0],w2,w3,w0);
///	dbg_printf ("v = %Zdx + %Zd\n",_ff_wrap_mpz(v[1],0), _ff_wrap_mpz(v[0],1));
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u, v, f) ) {
		err_printf ("hecurve_square: output verification failed for inputs:\n");
		err_printf ("      ");  hecurve_print(u1,v1);
		err_printf ("note that input has been modified if output overlapped.\n");
		exit (0);
	}
#endif
	return 1;

hecurve_g2_square_s1_zero:
///	dbg_printf ("Special case, s'1 == 0\n");
	// 6'. compute s and precomputations
	_ff_invert(w1,r);
	ff_mult(s0,s0,w1);
	_ff_mult(w2,s0,u1[0]);
	_ff_addto(w2,v1[0]);
///	dbg_printf ("w1 = %Zd, w2 = %Zd, s0 = %Zd\n",_ff_wrap_mpz(w1,0), _ff_wrap_mpz(w2,1), _ff_wrap_mpz(s0,2));

	// 7'. compute u' = (f-v^2)/u^2 - 2vs/u - s^2			note h =0
	_ff_square(w3,s0);
	_ff_add(w4,u1[1],u1[1]);
	_ff_addto(w3,w4);					
	_ff_neg(u[0],w3);							// potentially destroys u1[0] and u2[0]
///	dbg_printf ("u = x + %Zd\n",_ff_wrap_mpz(u[0],0));

	// 8'. compute v' = -(su+v) mod u
	_ff_sub(w1,u1[1],u[0]);
	ff_mult(w1,w1,s0);
	_ff_addto(w1,v1[1]);
	ff_mult(w1,w1,u[0]);
	_ff_sub(v[0],w1,w2);
	_ff_set_zero(v[1]);
///	dbg_printf ("v = %Zd\n",_ff_wrap_mpz(v[0],0));

	_ff_set_one(u[1]);							// set these last in case of overlap
	_ff_set_zero(u[2]);
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u, v,f) ) {
		err_printf ("hecurve_square output: verification failed for inputs:\n");
		err_printf ("      ");  hecurve_print(u1,v1);
		err_printf ("note that input has been modified if output overlapped.\n");
		exit (0);
	}
#endif
	return 1;
}

// Note that (u,v) == (u1,v1) or (u2,v2) (or both) are possible.  We do assume that any overlap is aligned.
// Thus modifying u[1] shouldn't effect u1[0] or u2[0] but may destory u1[1] and/or u2[1]
//
// This algorithm handles all the unusual cases - it assumes that either u1 or u2 has degree < 2
// or that Resultant(u1,u2) = 0, otherwise hecurve_compose would have handled it
int hecurve_g2_compose_special (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t u2[3], ff_t v2[3], ff_t f[6])
{
	_ff_t_declare_reg t1, t2, t3, r;
	_ff_t_declare *swap;
	int d1, d2, swapd, sts;

#if ! HECURVE_FAST
	dbg_printf ("Compose Special: (%Zdx^2+%Zdx+%Zd, %Zdx+%Zd) * (%Zdx^2+%Zdx+%Zd, %Zdx+%Zd)\n",
			    _ff_wrap_mpz(u1[2],0), _ff_wrap_mpz(u1[1],1), _ff_wrap_mpz(u1[0],2), _ff_wrap_mpz(v1[1],3), _ff_wrap_mpz(v1[0],4),
			    _ff_wrap_mpz(u2[2],5), _ff_wrap_mpz(u2[1],6), _ff_wrap_mpz(u2[0],7), _ff_wrap_mpz(v2[1],8), _ff_wrap_mpz(v2[0],9));
#endif
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u1, v1, f) )  { err_printf ("hecurve_compose_special: invalid input\n");  exit (0); }
	if ( ! hecurve_verify (u2, v2, f) ) { err_printf ("hecurve_compose_special: invalid input\n");  exit (0); }
#endif
	d1 = _deg_u(u1);
	d2 = _deg_u(u2);
	if ( d1 > d2 ) {
		swap = u2;  u2 = u1;  u1 = swap;
		swap = v2;  v2 = v1;  v1 = swap;
		swapd = d2;  d2 = d1;  d1 = swapd;
///		dbg_printf ("swapped u1 = x + %Zd, v1 = %Zd.\n",_ff_wrap_mpz(u1[0],0), _ff_wrap_mpz(v1[0],1));
///		dbg_printf ("swapped u2 = x^2 + %Zdx + %Zd, v2 = %Zdx + %Zd.\n",
///				   _ff_wrap_mpz(u2[1],0), _ff_wrap_mpz(u2[0],1), _ff_wrap_mpz(v2[1],2), _ff_wrap_mpz(v2[0],3));
	}
	switch(d1){
	case 0:		// deg u1 == 0, i.e. u1 is the identity
///		dbg_printf ("Case 1, deg(u1) == 0\n");
		_hecurve_set (u, v, u2, v2);
		break;
	case 1: 	// deg u1 == 1
		// case 2
///		dbg_printf ("Case 2, deg(u1) == 1\n");
		if ( d2 == 1 ) { // deg u1 == dege u2 == 1
///			dbg_printf ("Case 2.A, deg(u1) == deg(u2) == 1\n");
			// case 2.A
			if ( _ff_equal(u1[0],u2[0]) && _ff_equal(u1[1],u2[1]) ) {
				if ( _ff_zero(v1[0]) && _ff_zero(v2[0]) ) {
///					dbg_printf ("Case 2.A.i.a D1 == -D2 == D1\n");
					_hecurve_set_identity(u,v);
					break;
				}
				_ff_neg(t1,v2[0]);
				if ( _ff_equal(v1[0],t1) ) {
///					dbg_printf ("Case 2.A.i.b D1 == -D2\n");
					_hecurve_set_identity(u,v);
					break;
				}
				if ( _ff_equal(v1[0],v2[0]) ) {
///					dbg_printf ("Case 2.A.ii D1 == D2\n");
					hecurve_g2_square_1 (u, v, u1[0], v1[0], f);
					break;
				}
			}
///			dbg_printf ("Case 2.A.iii D1 != +- D2\n");
			if ( _ff_equal(u1[0],u2[0]) ) { err_printf ("u1[0] == u2[0] in case 2.A.iii of hecurve_compose!\n");  exit(0); }
			
			// u = u1*u2 and v = ((v2-v1)x+v2u1[0]-v1u2[0])/(u1[0]-u2[0])		note that v1 and v2 are constants (deg 0)
			_ff_sub(t1,u1[0],u2[0]);
			ff_invert(t1,t1);
			_ff_sub(t2,v2[0],v1[0]);
			_ff_mult(v[1],t2,t1);				// we assume writing [1] doesn't hurt [0]
			_ff_mult(t2,v2[0],u1[0]);
			_ff_mult(t3,v1[0],u2[0]);
			_ff_subfrom(t2,t3);
			_ff_mult(v[0],t2,t1);
			_ff_set_one(u[2]);					// update u last to handle overlap
			_ff_add(u[1],u1[0],u2[0]); 		// we assume writing [1] doesn't hurt [0]
			ff_mult(u[0],u1[0],u2[0]);		// mult should be safe
			break;
		} else { // deg u1 == 1, deg u2 == 2
///			dbg_printf ("Case 2.B, deg(u1) == 1, deg(u2) == 2\n");
			// compute u2(-u1[0]) = (-u1[0])(-u1[0] + u2[1]) + u2[0]
			_ff_neg(t1,u1[0]);				// save t1 = -u1[0] for later
			_ff_add(t2,t1,u2[1]);
			ff_mult(t2,t2,t1);
			_ff_addto(t2,u2[0]);
			if ( _ff_nonzero(t2) ) {
///			dbg_printf ("Case 2.B.i, u2(-u1[0]) != 0\n");
				hecurve_g2_compose_1_2 (u, v, u1, v1, u2, v2, f);
				break;
			}
///			dbg_printf ("Case 2.B.ii\n");
			_ff_add(t2,u1[0],u1[0]);
			if ( _ff_equal(u2[1],t2) ) {
				_ff_square(t2,u1[0]);
				if ( _ff_equal(u2[0],t2) ) {
///					dbg_printf ("Case 2.B.ii.a, u2[1] == 2u1[0] and u2[0] == u1[0]^2\n");
					// compute v2(-u1[0])
					_ff_mult(t2,v2[1],t1);
					_ff_addto(t2,v2[0]);
					if ( _ff_equal(t2,v1[0]) ) {
						_ff_t_declare_dynamic u3[3], v3[2], v4[2];
						
						// this is expensive, but this is a special case - may optimize later to handle recursion more elegantly
						_ff_init(u3[0]);  _ff_init(u3[1]);  _ff_init(u3[2]);
						_ff_init(v3[0]);  _ff_init(v3[1]);  _ff_init(v4[0]);  _ff_init(v4[1]);
						
///						dbg_printf ("Case 2.B.ii.a.I D2 = 2D1\n");
						// Double D2 then subtract D1
						// this is slightly inefficient but easier and is a rare case anyway
						hecurve_square (u3, v3, u2, v2, f, 0);  // double D2
						// Bug fix - need to avoid infinite recursion on order 3 elements, this occurs when 2D2 == -D2
						if ( _ff_equal(u3[0],u2[0]) && _ff_equal(u3[1],u2[1]) && _ff_equal(u3[2],u2[2]) ) {
							_ff_neg(t2,v2[0]);
							_ff_neg(t3,v2[1]);
							if ( _ff_equal(v3[0],t2) && _ff_equal(v3[1],t3) ) {
								hecurve_compose_cantor (u, v, u1, v1, u2, v2, f);
								break;
//								err_printf ("Attempt to compute 2D+D for D with order 6 - unable to handle this case\n");
//								exit (0);
							}
						}
///						dbg_printf ("2D2 != -2D1, computing 2D2-D1\n");
						_ff_neg(v4[0],v1[0]);
						_ff_neg(v4[1],v1[1]);						// D4 = invert D1, note u4 = u1
						hecurve_compose (u, v, u1, v4, u3, v3, f, 0);		// compute D2+D2-D1
						_ff_clear(u3[0]);  _ff_clear(u3[1]);  _ff_clear(u3[2]);
						_ff_clear(v3[0]);  _ff_clear(v3[1]);  _ff_clear(v4[0]);  _ff_clear(v4[1]);
					} else {
///						dbg_printf ("Case 2.B.ii.a.II, D2 = -2P1\n");
						_hecurve_set (u, v, u1, v1);
						hecurve_invert (u, v);						// compute inverse of D1
					}
					break;
				}
			} else {
///				dbg_printf ("Case 2.B.ii.b, u2[1] != 2u1[0] or u2[0] != u1[0]^2\n");
				// compute -v2(-u1[0]) - this is a bug fix - p. 314/315 of handbook uses v2(-u1[0]) (or not)
				_ff_mult(t3,v2[1],t1);		// t1 = -u1[0]
				_ff_addto(t3,v2[0]);
				_ff_neg(t2,t3); 
				// we need to compute v2(u1[0]-u2[1]) either way
				_ff_sub(t3,u1[0],u2[1]);
				ff_mult(t3,t3,v2[1]);
				_ff_addto(t3,v2[0]);
///				dbg_printf ("t3 = v2[1]*(u1[0]-u2[1]) + v2[0] = %Zd\n",_ff_wrap_mpz(t3,0));
				if ( _ff_equal(t2,v1[0]) ) {
///					dbg_printf ("Case 2.B.ii.b.I, -P1 occurs in D2\n");
					_ff_add(u[0],u2[1],t1);						// t1 = -u1[0]
					_ff_set_zero(u[2]);  _ff_set_one(u[1]);				// u[1] must get set here, could destory u2[1]
					_ff_set_zero(v[1]);
					_ff_set(v[0],t3);
				} else {
					_ff_t_declare_dynamic u3[3], v3[2], u4[3], v4[2];
					
					// this is expensive, but this is a special case - may optimize later to handle recursion more elegantly
					_ff_init(u3[0]);  _ff_init(u3[1]);  _ff_init(u3[2]);  _ff_init(v3[0]);  _ff_init(v3[1]);
					_ff_init(u4[0]);  _ff_init(u4[1]);  _ff_init(u4[2]);  _ff_init(v4[0]);  _ff_init(v4[1]);

///					dbg_printf ("Case 2.B.ii.b.II, Doubling D1 and composing with (x+u21-u10, v2(u10-u21)),\n");
					hecurve_g2_square_1 (u3, v3, u1[0], v1[0], f);
					_ff_set_zero(u4[2]);  _ff_set_one(u4[1]);
					_ff_sub(u4[0],u2[1],u1[0]);
					_ff_set_zero(v4[1]);
					_ff_set(v4[0],t3);
					hecurve_compose (u, v, u4, v4, u3, v3, f, 0);
					_ff_clear(u3[0]);  _ff_clear(u3[1]);  _ff_clear(u3[2]);  _ff_clear(v3[0]);  _ff_clear(v3[1]);
					_ff_clear(u4[0]);  _ff_clear(u4[1]);  _ff_clear(u4[2]);  _ff_clear(v4[0]);  _ff_clear(v4[1]);
				}
			}
		}
		break;
	case 2: // deg u1 == deg u2 == 2 (case 3)
///		dbg_printf ("Case 3, deg(u1) == deg(u2) == 2\n");
		if ( _ff_equal(u1[0],u2[0]) && _ff_equal(u1[1],u2[1]) ) {
///			dbg_printf ("Case 3.A, u1 == u2\n");
			if ( _ff_equal(v1[1],v2[1]) && _ff_equal(v1[0],v2[0]) ) {
				// catch inverse case for order 2 elements
				if ( _ff_zero(v1[1]) && _ff_zero(v1[0]) ) {
///					dbg_printf ("Case 3.A.i*, v1 == -v2 == 0\n");
					_hecurve_set_identity (u, v);
					break;
				}
///				dbg_printf ("Case 3.A.ii, v1 == v2\n");
				hecurve_square (u, v, u1, v1, f, 0);				// handles case 3.A.ii.b (zero resultant) also
				break;
			}
			_ff_neg(t1,v2[0]);
			_ff_neg(t2,v2[1]);
			if ( _ff_equal(v1[0],t1) && _ff_equal(v1[1],t2) ) {
///				dbg_printf ("Case 3.A.i, v1 == -v2\n");
				_hecurve_set_identity (u, v);
				break;
			}
///			dbg_printf ("Case 3.A.iii v1 != +- v2\n");
			// case 3.A.iii
			_ff_sub(t1,v2[1],v1[1]);
#if ! HECURVE_FAST
			if ( _ff_zero(t1) ) { err_printf ("v2[1] == v1[1] in case 3.A.iii of hecurve_compose_special\n");  exit (0); }
#endif
			_ff_invert(t2,t1);
			_ff_sub(t1,v1[0],v2[0]);
			ff_mult(t1,t1,t2);
			_ff_neg(t2,t1);
			_ff_mult(t3,v1[1],t1);
			_ff_addto(t3,v1[0]);
			hecurve_g2_square_1 (u, v, t2, t3, f);
		} else { // u1 != u2
///			dbg_printf ("Case 3.B, u1 != u2\n");
			// case 3.B
			// We can assume Resultant(u1,u2) == 0, since otherwise we wouldn't have been called
			_ff_t_declare_dynamic x2, y2, x3, y3;
			
			// this is expensive, but this is a special case - may optimize later to handle recursion more elegantly
			_ff_init(x2);  _ff_init(x3);  
			_ff_init(y2);  _ff_init(y3);  
		
///			dbg_printf ("Case 3.B.ii, res(u1,u2) == 0\n");
			// need to compute gcd(u1,u2) = x-x1 where x1 = (u1[1]-u2[1])/(u2[0]-u1[0]))
			_ff_sub(t1,u1[1],u2[1]);
#if ! HECURVE_FAST
			if ( _ff_zero(t1) ) { err_printf ("u1[1] == u2[1] in case 3.B.ii of hecurve_compose_special\n");  exit (0); }
#endif
			_ff_invert(t2,t1);
			_ff_sub(t1,u2[0],u1[0]);
			_ff_mult(t3,t1,t2);
///			dbg_printf ("gcd(u1,u2) = x-%Zd\n", _ff_wrap_mpz(t3,0));
			// compute v1(x1) and v2(x1) and compare
			_ff_mult(t1,v1[1],t3);
			_ff_addto(t1,v1[0]); 
			_ff_mult(t2,v2[1],t3);
			_ff_addto(t2,v2[0]); 
			// Need to extract coords of P2 and P3 in either case, so compute them here
			// This code could be cleaned up a bit - there is some double negation, but it is a special case...
			_ff_add(y2,u1[1],t3); 
			_ff_neg(x2,y2);  
			_ff_mult(y2,v1[1],x2);
			_ff_addto(y2,v1[0]);
			_ff_add(y3,u2[1],t3);
			_ff_neg(x3,y3); 
			_ff_mult(y3,v2[1],x3);
			_ff_addto(y3,v2[0]); 
///			dbg_printf ("P2 = (%Zd,%Zd), P3 = (%Zd,%Zd)\n", _ff_wrap_mpz(x2,0), _ff_wrap_mpz(y2,1), _ff_wrap_mpz(x3,2), _ff_wrap_mpz(y3,3));
			if ( _ff_equal(t1,t2) ) {
				_ff_t_declare_dynamic u3[3], v3[2], u4[3], v4[2], u5[3], v5[2];
				
				// this is expensive, but this is a special case - may optimize later to handle recursion more elegantly
				_ff_init(u3[0]);  _ff_init(u3[1]);  _ff_init(u3[2]);  _ff_init(v3[0]);  _ff_init(v3[1]);
				_ff_init(u4[0]);  _ff_init(u4[1]);  _ff_init(u4[2]);  _ff_init(v4[0]);  _ff_init(v4[1]);
				_ff_init(u5[0]);  _ff_init(u5[1]);  _ff_init(u5[2]);  _ff_init(v5[0]);  _ff_init(v5[1]);
			
///					dbg_printf ("Case 3.B.ii.a, v1(x1) == v2(x1)\n");
				_ff_neg(t1,t3);
				_ff_mult(t2,v1[1],t3);
				_ff_add(t2,t2,v1[0]);
				hecurve_g2_square_1 (u4, v4, t1, t2, f);		// D' = 2(P1)
///				dbg_printf ("D' = (%Zdx^2 + %Zdx + %Zd, %Zdx + %Zd)\n",_ff_wrap_mpz(u4[2],0), _ff_wrap_mpz(u4[1],1), _ff_wrap_mpz(u4[0],2), _ff_wrap_mpz(v4[1],3), _ff_wrap_mpz(v4[0],4));
				// could we use hecurve_compose_1_2 here?
				_ff_set_zero(u3[2]);  _ff_set_one(u3[1]);
				_ff_neg(u3[0],x2);
				_ff_set_zero(v3[1]);
				_ff_set(v3[0],y2);
				hecurve_compose (u5, v5, u3, v3, u4, v4, f, 0);		// D'' = D' + (P2)
///				dbg_printf ("D'' = (%Zdx^2 + %Zdx + %Zd, %Zdx + %Zd)\n",_ff_wrap_mpz(u5[2],0), _ff_wrap_mpz(u5[1],1), _ff_wrap_mpz(u5[0],2), _ff_wrap_mpz(v5[1],3), _ff_wrap_mpz(v5[0],4));
				_ff_neg(u3[0],x3);
				_ff_set(v3[0],y3);
				hecurve_compose (u, v, u3, v3, u5, v5, f, 0);		// D = D'' + (P3)
				_ff_clear(u3[0]);  _ff_clear(u3[1]);  _ff_clear(u3[2]);  _ff_clear(v3[0]);  _ff_clear(v3[1]);
				_ff_clear(u4[0]);  _ff_clear(u4[1]);  _ff_clear(u4[2]);  _ff_clear(v4[0]);  _ff_clear(v4[1]);
				_ff_clear(u5[0]);  _ff_clear(u5[1]);  _ff_clear(u5[2]);  _ff_clear(v5[0]);  _ff_clear(v5[1]);
			} else {
///				dbg_printf ("Case 3.B.ii.a, v1(x1) != v2(x1)\n");
				hecurve_make_2 (u, v, x2, y2, x3, y3);
			}
			_ff_clear(x2);  _ff_clear(x3);  
			_ff_clear(y2);  _ff_clear(y3);  				
		}
		break;
	default:
		err_printf ("Invalid degree in hecurve_compose!\n");  exit (0);
	}
#if ! HECURVE_FAST
	_dbg_print_uv ("Compose Special Result: ", u, v);
#endif
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u, v, f) ) {
		err_printf ("hecurve_compose_special output verification failed for inputs:\n");
		err_printf ("      ");  hecurve_print(u1,v1);
		err_printf ("      ");  hecurve_print(u2,v2);
		err_printf ("note that inputs have been modified if output overlapped.\n");
		exit (0);
	}
#endif
	return 1;
}


// Algorithm 14.20 p.318
// As above, handle aligned overlap of (u,v) with (u1,v1) and/or (u2,v2)
void hecurve_g2_compose_1_2 (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t u2[3], ff_t v2[2], ff_t f[6])
{
	_ff_t_declare_reg r, inv, t1, t2, w1, L0, L1, s0;

///	dbg_printf ("Compose_1_2\n");	

	// 1. compute r = u2 mod u1
	_ff_sub(w1,u2[1],u1[0]);
	ff_mult(w1,w1,u1[0]);
	_ff_sub(r,u2[0],w1);
///	dbg_printf ("r = %Zd\n",_ff_wrap_mpz(r,0));
	
	// 2. compute inverse of u2 mod u1
	_ff_invert(inv, r);
///	dbg_printf ("inv = %Zd\n",_ff_wrap_mpz(inv,0));
	
	// 3. compute s = ((v1-v2)inv) mod u1
	_ff_mult(w1,v2[1],u1[0]);			// BUG FIX - handbook incorrectly uses -v2[1]u1[0] - should be positive
	_ff_addto(w1,v1[0]);		
	_ff_subfrom(w1,v2[0]); 
	_ff_mult(s0,w1,inv);
///	dbg_printf ("s = %Zd\n", _ff_wrap_mpz(s0,0));
	
	// 4. compute L = su2 = s0*x^2 + L1*x + L0
	_ff_mult(L1,s0,u2[1]);
	_ff_mult(L0,s0,u2[0]);
	// no reduction required for L coeffs
///	dbg_printf ("L = su2 = %Zdx^2 + %Zdx + %Zd\n",_ff_wrap_mpz(s0,0), _ff_wrap_mpz(L1,1), _ff_wrap_mpz(L0,2));
	
	// 5. compute t = (f - v2^2)/u2 = x^3 + t2*x2 + t1*x + t0   note h = 0
	_ff_neg(t2,u2[1]);
	_ff_square(t1,u2[1]);
	_ff_subfrom(t1,u2[0]);
	if ( _ff_nonzero(f[3]) ) _ff_addto(t1,f[3]);
///	dbg_printf ("t = x^3 + %Zdx^2 + %Zdx + t0\n",_ff_wrap_mpz(t2,0), _ff_wrap_mpz(t1,1));
	
	// 6. compute u = (t-s(L+2v2))/u1 = x^2 + u[1]*x + u[0]	 	note h = 0
	_ff_set_one(u[2]);
	_ff_square(r,s0);
	_ff_sub(w1,t2,r);
	_ff_sub(u[1],w1,u1[0]);			// ther should be no overlap - potentially destroys u1[1] and u2[1]
	_ff_add(w1,v2[1],v2[1]);
	_ff_addto(w1,L1);
	ff_mult(w1,w1,s0);
	_ff_sub(r,t1,w1);
	_ff_mult(t2,u1[0],u[1]);
	_ff_sub(u[0],r,t2);								// potentially destroys u1[0] and u2[0]
///	dbg_printf ("u = x^2 + %Zdx + %Zd\n",_ff_wrap_mpz(u[1],0), _ff_wrap_mpz(u[0],1));
	
	// 7. compute v = (-(L+v2)) mod u = v[1]x + v[0]			note h = 0
	_ff_mult(w1,s0,u[1]);
	_ff_subfrom(w1,L1);
	_ff_subfrom(w1,v2[1]);			// v2[1] could overlap v[1]
	_ff_set (v[1],w1);
	_ff_mult(w1,s0,u[0]);
	_ff_subfrom(w1,L0);
	_ff_subfrom(w1,v2[0]);			// v2[0] could overlap v[0]
	_ff_set(v[0],w1);
///	dbg_printf ("v = %Zdx + %Zd\n",_ff_wrap_mpz(v[1],0), _ff_wrap_mpz(v[0],1));
}


// This code would benefit from precomputing f' - todo later
// Note that overlap of u[0] and u0 (and v[0] and v0) is possible!
void hecurve_g2_square_1 (ff_t u[3], ff_t v[2], ff_t u0, ff_t v0, ff_t f[6])
{
	_ff_t_declare_reg x, y, t1, t2, t3;
	
	if ( _ff_zero(v0) ) { _hecurve_set_identity (u, v);  return; }
	// The code below corresponds to (14.9) on p. 314 of Handbook of E+HE Crypto
///	dbg_printf ("Square_1 (x + %Zd,%Zd)\n",_ff_wrap_mpz(u0,0),_ff_wrap_mpz(v0,0));
	// u = u1^2
	_ff_set (t3, u0);				// save t3=u0 since we need it later and may overwrite it here.
	_ff_set_one(u[2]);
	_ff_add(u[1],t3,t3);
	_ff_square(u[0],t3);
///	dbg_printf ("u = x^2 + %Zdx + %Zd\n",_ff_wrap_mpz(u[1],0), _ff_wrap_mpz(u[0],1));
	// v = (f'(-u1[0])x + f'(-u1[0])u1[0])/(2v1) + v1
	// Compute y = f'(-u1[0])  - note that this code assumes that 5* max coeff of f fits in an ff_t
	_ff_neg(x,t3);
	_ff_set(y,f[1]);
	if ( _ff_nonzero(f[2]) ) {
		_ff_add(t2,f[2],f[2]);
		ff_mult(t2,t2,x);
		_ff_addto(y,t2);
	}
	_ff_square(t1,x);
	if ( _ff_nonzero(f[3]) ) {
		_ff_add(t2,f[3],f[3]); 
		_ff_addto(t2,f[3]);
		ff_mult(t2,t2,t1);
		_ff_addto(y,t2);
	}
#if FF_WORDS == 1
	if ( _ff_p != 5 ) {
#endif
		ff_square(t1,t1);
		_ff_add(t2,t1,t1); 
		_ff_x2(t2); 
		_ff_addto(t2,t1);
		_ff_addto(y,t2); 
#if FF_WORDS == 1
	}
#endif
///	dbg_printf ("f'(-u0) = %Zd\n",_ff_wrap_mpz(y,0));
	_ff_add(t2,v0,v0);
	_ff_invert(t1,t2);
	_ff_mult(v[1],t1,y);
	_ff_mult(t1,v[1],t3);
	_ff_add(v[0],t1,v0);
///	dbg_printf ("v = %Zdx + %Zd\n",_ff_wrap_mpz(v[1],0), _ff_wrap_mpz(v[0],0));
}

/*
	This algorithm handles squaring a degree 2 element where Resultant(2v,u) = 0.
	In this case the divisor contains a point with order 2 and we simply want to square the other point.
	This corresponds to case 3.A.ii.b in the general algorithm, but may also be detected in hecurve_square
*/
void hecurve_g2_square_2_r0 (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t f[6])
{
	_ff_t_declare_reg t1, t2;

	// Since h = 0, we simply need gcd(2v1, u1) which must be (x+v1[0]/v1[1]) hence x1 = -v1[0]/v1[1] is the x-coord of P1
	// We want to then double P2 = (-u1[1]-x1,v1(-u1[1]-x1))  see p. 315 of the handbook
	if ( _ff_zero(v1[1]) ) {
		if ( _ff_zero(v1[0]) ) { _hecurve_set_identity (u, v);  return; }
		// this should be impossible if Results(u,2v) = 0
		err_printf ("v1[1] = 0 and v1[0] != 0 in hecurve_square_2_r0!");
		exit (0);		
	}
	_ff_invert(t1,v1[1]);
	_ff_neg(t2,v1[0]);
	ff_mult(t1,t1,t2);
	_ff_addto(t1,u1[1]);
	_ff_neg(t2,t1);
	ff_mult(t2,t2,v1[1]);
	_ff_addto(t2,v1[0]);
	hecurve_g2_square_1 (u, v, t1, t2, f);
}


#endif
