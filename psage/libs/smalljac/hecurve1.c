#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "gmp.h"
#include "mpzutil.h"
#include "hecurve.h"
#include "ffwrapper.h"
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

/*
	This module implements generic Jacobian group operation in genus 1, i.e.
	point addition (and doubling) on elliptic curves.

	For consistency with higher genera we use a Mumford representation,
	the point (x,y) is represented by u(z) = z-x and v(z) = y.
	Thus x = -u[0] and y = v[0].  The identity has u(z) = 1 and v(z) = 0.
	(*** note the sign difference on x and u[0] ***)
	
	There is minimal overhead in doing this .  u[1] is, in fact, never accessed
	and need not be allocated.
	
	We use an affine representation, as elsewhere, so that we can benefit
	from parallel operation, which is achieved via the ctx parameter.  If
	ctx is non-null and a field inversion is required, the functions will fill
	ctx with the entry to be inverted along with other state information and
	return 0 to the caller.  The caller should (eventually) perform the required
	inversion (along with a bunch of others) and call back to finish the operation.	
*/

#define _dbg_print_uv(s,u,v)	if ( dbg_level >= 2 ) { printf (s);  hecurve_print(u,v); }

int hecurve_g1_compose (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS], ff_t u1[HECURVE_GENUS+1],
				         ff_t v1[HECURVE_GENUS], ff_t u2[HECURVE_GENUS+1],
				         ff_t v2[HECURVE_GENUS+1], ff_t f[HECURVE_DEGREE+1], hecurve_ctx_t *ctx)
{
	_ff_t_declare_reg m, t1, t2;
#ifdef FF_INIT_REQUIRED
	static int init;
	
	if ( ! init ) { _ff_init (m);  _ff_init (t1);  _ff_init (t2);  init = 1; }
#endif
	
#if ! HECURVE_FAST
	_dbg_print_uv ("Genus 1 Compose A: ", u1, v1);
	_dbg_print_uv ("Genus 1 Compose B: ", u2, v2);
#endif	
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u1, v1, f) ) { err_printf ("hecurve_compose_g1 input verification failed for input: "); hecurve_print(u1,v1);  exit (0); }
	if ( ! hecurve_verify (u2, v2, f) ) { err_printf ("hecurve_compose_g1 input verification failed for input: "); hecurve_print(u2,v2);  exit (0); }
#endif

	if ( _hecurve_is_identity (u1,v2) ) { _hecurve_set(u,v,u2,v2);  return 1; }
	if ( _hecurve_is_identity (u2,v2) ) { _hecurve_set(u,v,u1,v1);  return 1; }
	
	if ( ctx && ctx->state ) {
		if ( ctx->state == 1 ) {
			// get invert result and restore saved variables
			_ff_set (t2, ctx->invert);		// inverted value of t1 = x1-x2 = u2[0]-u1[0]
			goto hecurve_g1_compose_inverted;
		} else if ( ctx->state == 2 ) {
			// state 2 indicates square completing
			return hecurve_g1_square (u, v, u1, v1, f, ctx);
		}
	}
	
	_ff_sub (t1, u2[0], u1[0]);
	if ( _ff_zero(t1) ) {
		if ( _ff_equal(v1[0],v2[0]) ) {
			return hecurve_g1_square (u, v, u1, v1, f, ctx);
		} else {
			_hecurve_set_identity (u,v);
			return 1;
		}
	}

	if ( ctx && ! ctx->state ) {
		_ff_set (ctx->invert, t1);	// return to let caller invert
		ctx->state = 1;
		return 0;
	}
	_ff_invert(t2,t1);
	
hecurve_g1_compose_inverted:	// caller returns here after inverting - we have t2 = 1/(x1-x2) = 1/(u2[0]-u1[0])

	_ff_sub (t1,v1[0],v2[0]);
	_ff_mult (m,t1,t2);
	_ff_square (t1,m);
	_ff_add (t2,u1[0],u2[0]);
	_ff_addto (t2,t1);
	_ff_neg (t1,t2);			// t1 = -x3 = -(m^2+u1[0]+u2[0]) = -(m^2-x1-x2)	
	_ff_sub (t2, t1, u1[0]);
	_ff_set_one (u[1]);
	_ff_set (u[0], t1);			// must wait to set u[0] until here due to possible overlap with u1 or u2
	_ff_mult (t1, t2, m);
	_ff_subfrom (t1, v1[0]);
	_ff_set (v[0], t1);
#if ! HECURVE_FAST
	_dbg_print_uv ("Genus 1 Compose Result: ", u, v);
#endif
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u, v, f) ) {
		err_printf ("hecurve_compose output verification failed for inputs:\n");
		err_printf ("      ");  hecurve_print(u1,v1);
		err_printf ("      ");  hecurve_print(u2,v2);
		err_printf ("note that inputs have been modified if output overlapped input.\n");
		exit (0);
	}
#endif
	return 1;
}


int hecurve_g1_square (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS],
				      ff_t u1[HECURVE_GENUS+1], ff_t v1[HECURVE_GENUS], ff_t f[HECURVE_DEGREE+1], hecurve_ctx_t *ctx)
{
	_ff_t_declare_reg m, t1, t2;
#ifdef FF_INIT_REQUIRED
	static int init;
	
	if ( ! init ) { _ff_init (m);  _ff_init (t1);  _ff_init (t2);  init = 1; }
#endif

#if HECURVE_VERIFY
	if ( ! hecurve_verify (u1, v1, f) ) { err_printf ("hecurve_square_g1 input verification failed for input: "); hecurve_print(u1,v1);  exit (0); }
#endif

	if ( ctx && ctx->state == 2 ) {
		// get invert result and restore saved variables
		_ff_set (t2, ctx->invert);		// inverted value of t1 = x1-x2 = u2[0]-u1[0]
		goto hecurve_g1_square_inverted;
	}

#if ! HECURVE_FAST
	_dbg_print_uv ("Genus 1 Square: ", u1, v1);
#endif	
	
	if ( _hecurve_2tor (u1,v1) ) { _hecurve_set_identity (u,v);  return 1; }

	_ff_add (t1, v1[0], v1[0]);
	if ( ctx && ! ctx->state ) {
		_ff_set (ctx->invert, t1);
		ctx->state = 2;
		return 0;
	}
	_ff_invert (t2,t1);

hecurve_g1_square_inverted:			// caller returns here after inverting - t2 = 1/(2y1) = 1/(2v1[0])
	_ff_square (t1, u1[0]);
	_ff_add (m, t1, t1);
	_ff_addto (t1, m);
	_ff_addto (t1, f[1]);
	_ff_mult (m, t1, t2);				// m = (3x1^2 + f1)/2y1 = 3(u1[0]^2 + f1)/(2v1[0])
	_ff_square (t1,m);
	_ff_add (t2,u1[0],u1[0]);
	_ff_addto (t2,t1);
	_ff_neg (t1, t2);
	_ff_sub (t2,t1,u1[0]);
	ff_mult (t2, t2, m);
	_ff_subfrom (t2, v1[0]);
	_ff_set_one (u[1]);
	_ff_set (u[0], t1);				// note that we must wait until here to set u[0] due to possible overlap with u1 or u2
	_ff_set (v[0], t2);
#if ! HECURVE_FAST
	_dbg_print_uv ("Genus 1 square Result: ", u, v);
#endif
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u, v, f) ) {
		err_printf ("hecurve_square output verification failed for input:\n");
		err_printf ("      ");  hecurve_print(u1,v1);
		err_printf ("note that inputs have been modified if output overlapped input.\n");
		exit (0);
	}
#endif
	return 1;
}
