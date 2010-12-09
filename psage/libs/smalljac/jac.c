#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include "gmp.h"
#include "mpzutil.h"
#include "ff.h"
#include "ffpoly.h"
#include "hecurve.h"
#include "jac.h"
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
	Basic group operations for Jacobians are implemented here, plus some generic
	fastorder computation functions (computations of element orders from
	a known exponent).
*/	


#define JAC_MAX_FASTORDER_UI_W		20


// o1 should not overlap a2 or b2
void jac_mult2 (jac_t o1[1], jac_t a1[1], jac_t b1[1], jac_t o2[1], jac_t a2[1], jac_t b2[1], curve_t c[1])
{
	hecurve_ctx_t ctx[2];
	ff_t inverts[2];
	
	ctx[0].state = ctx[1].state = 0;
	if ( ! __jac_mult (o1[0], a1[0], b1[0], c[0], ctx) ) {
		if ( ! __jac_mult (o2[0], a2[0], b2[0], c[0], ctx+1) ) {
			_ff_set (inverts[0], ctx[0].invert);
			_ff_set (inverts[1], ctx[1].invert);
			ff_parallel_invert (inverts, inverts, 2);
			_ff_set (ctx[0].invert, inverts[0]);
			_ff_set (ctx[1].invert, inverts[1]);
			__jac_mult (o1[0], a1[0], b1[0], c[0], ctx) ;
			__jac_mult (o2[0], a2[0], b2[0], c[0], ctx+1) ;
		} else {
			_ff_invert (ctx[0].invert, ctx[0].invert);
			__jac_mult (o1[0], a1[0], b1[0], c[0], ctx);
		}
	} else {
		__jac_mult (o2[0], a2[0], b2[0], c[0], 0);
	}
	jac_gops += 2;
}


// o1 should not overlap a2 or b2
void jac_square2 (jac_t o1[1], jac_t a1[1], jac_t o2[1], jac_t a2[1], curve_t c[1])
{
	hecurve_ctx_t ctx[2];
	ff_t inverts[2];
	
	ctx[0].state = ctx[1].state = 0;
	if ( ! __jac_square (o1[0], a1[0], c[0], ctx) ) {
		if ( ! __jac_square (o2[0], a2[0], c[0], ctx+1) ) {
			_ff_set (inverts[0], ctx[0].invert);
			_ff_set (inverts[1], ctx[1].invert);
			ff_parallel_invert (inverts, inverts, 2);
			_ff_set (ctx[0].invert, inverts[0]);
			_ff_set (ctx[1].invert, inverts[1]);
			__jac_square (o1[0], a1[0], c[0], ctx) ;
			__jac_square (o2[0], a2[0], c[0], ctx+1) ;
		} else {
			_ff_invert (ctx[0].invert, ctx[0].invert);
			__jac_square (o1[0], a1[0], c[0], ctx);
		}
	} else {
		__jac_square (o2[0], a2[0], c[0], 0);
	}
	jac_gops += 2;
}


// replaces a by a^2 and b by ab - a and b should not overlap 
void jac_square_mult (jac_t a[1], jac_t b[1], curve_t c[1])
{
	hecurve_ctx_t ctx[2];
	jac_t t;
	ff_t inverts[2];
	
	ctx[0].state = ctx[1].state = 0;
	if ( ! __jac_mult (b[0], a[0], b[0], c[0], ctx) ) {
		if ( ! __jac_square (t, a[0], c[0], ctx+1) ) {			// can't overwrite a - needed to finish ab and square might not return for inversion
			_ff_set (inverts[0], ctx[0].invert);
			_ff_set (inverts[1], ctx[1].invert);
			ff_parallel_invert (inverts, inverts, 2);
			_ff_set (ctx[0].invert, inverts[0]);
			_ff_set (ctx[1].invert, inverts[1]);
			__jac_mult (b[0], a[0], b[0], c[0], ctx) ;
			__jac_square (a[0], a[0], c[0], ctx+1) ;			// but we can overwrite a here, in which case we didn't actually use t - nice
		} else {
			_ff_invert (ctx[0].invert, ctx[0].invert);
			__jac_mult (b[0], a[0], b[0], c[0], ctx);
			_jac_set (a[0], t);							// here's where we need to use t
		}
	} else {
		__jac_square (a[0], a[0], c[0], 0);
	}
	jac_gops += 2;
}


void jac_parallel_mult (jac_t o[], jac_t a[], jac_t b[], curve_t c[], int n)
{
	int inv[JAC_MAX_CURVES];
	ff_t inverts[JAC_MAX_CURVES];
	hecurve_ctx_t ctx[JAC_MAX_CURVES];
	register int i, j;
	
	j = 0;
	for ( i = 0 ; i < n ; i++ ) {
		ctx[i].state = 0;
		if ( ! __jac_mult (o[i], a[i], b[i], c[i], ctx+i) ) {
			_ff_set (inverts[j], ctx[i].invert);
			inv[j++] =i;
		}
	}
	ff_parallel_invert (inverts, inverts, j);	// we could modify this to avoid overlap if needed
	for ( i = 0 ; i < j ; i++ ) {
		_ff_set (ctx[inv[i]].invert, inverts[i]);
		__jac_mult (o[inv[i]], a[inv[i]], b[inv[i]], c[inv[i]], ctx+inv[i]);
	}
	jac_gops += n;
}


void jac_parallel_mult_c (jac_t o[], jac_t a[], jac_t b[], curve_t c[1], int n)
{
	int inv[JAC_MAX_CURVES];
	ff_t inverts[JAC_MAX_CURVES];
	hecurve_ctx_t ctx[JAC_MAX_CURVES];
	register int i, j;
	
	j = 0;
	for ( i = 0 ; i < n ; i++ ) {
		ctx[i].state = 0;
		if ( ! __jac_mult (o[i], a[i], b[i], c[0], ctx+i) ) {
			_ff_set (inverts[j], ctx[i].invert);
			inv[j++] =i;
		}
	}
	ff_parallel_invert (inverts, inverts, j);	// we could modify this to avoid overlap if needed
	for ( i = 0 ; i < j ; i++ ) {
		_ff_set (ctx[inv[i]].invert, inverts[i]);
		__jac_mult (o[inv[i]], a[inv[i]], b[inv[i]], c[0], ctx+inv[i]);
	}
	jac_gops += n;
}


void jac_parallel_mult_1 (jac_t o[], jac_t a[], jac_t b[1], curve_t c[1], int n)
{
	int inv[JAC_MAX_CURVES];
	ff_t inverts[JAC_MAX_CURVES];
	hecurve_ctx_t ctx[JAC_MAX_CURVES];
	register int i, j;
	
	j = 0;
	for ( i = 0 ; i < n ; i++ ) {
		ctx[i].state = 0;
		if ( ! __jac_mult (o[i], a[i], b[0], c[0], ctx+i) ) {
			_ff_set (inverts[j], ctx[i].invert);
			inv[j++] =i;
		}
	}
	ff_parallel_invert (inverts, inverts, j);	// we could modify this to avoid overlap if needed
	for ( i = 0 ; i < j ; i++ ) {
		_ff_set (ctx[inv[i]].invert, inverts[i]);
		__jac_mult (o[inv[i]], a[inv[i]], b[0], c[0], ctx+inv[i]);
	}
	jac_gops += n;
}


void jac_parallel_square (jac_t o[], jac_t a[], curve_t c[], int n)
{
	int inv[JAC_MAX_CURVES];
	ff_t inverts[JAC_MAX_CURVES];
	hecurve_ctx_t ctx[JAC_MAX_CURVES];
	register int i, j;

	j = 0;
	for ( i = 0 ; i < n ; i++ ) {
		ctx[i].state = 0;
		if ( ! __jac_square (o[i], a[i], c[i], ctx+i) ) {
			_ff_set (inverts[j], ctx[i].invert);
			inv[j++] =i;
		}
	}
	ff_parallel_invert (inverts, inverts, j);	// we could modify this to avoid overlap if needed
	for ( i = 0 ; i < j ; i++ ) {
		_ff_set ( ctx[inv[i]].invert, inverts[i]);
		__jac_square (o[inv[i]], a[inv[i]], c[inv[i]], ctx+inv[i]);
	}
	jac_gops += n;
}

// o cannot overlap a[]
void jac_exp_ui32_powers (jac_t o[1], jac_t a[], uint32_t e, curve_t c[1])
{
	register int i;
	
	if ( ! e ) { _jac_set_identity (o[0]);  return; }
	// could use parallel mult here but don't bother for now
	for ( i = 0 ; ! (e & 0x1) ; e >>= 1, i++ );
	_jac_set (o[0], a[i]);
	for ( e >>= 1, i++ ; e ; e >>= 1, i++ ) {
		if ( (e & 0x1) ) _jac_mult (o[0], o[0], a[i], c[0]);
	}
}

/*
	Generic exponentiation optimized to take advantage of a fast parallel group operation.
	Splits exponent and does two smaller exponentiations in parallel, then shifts (squares) and combines.
*/
void jac_exp_ui (jac_t o[1], jac_t a[1], unsigned long e, curve_t c[1])
{
	register int i, j;
	unsigned long e1, e2;
	jac_t b[2];

#if SMALLJAC_GENUS == 1
	hecurve_g1_exp_ui(o[0].u,o[0].v,a[0].u,a[0].v,e,c[0].f);
	jac_gops += (3*ui_len(e))/2;
	return;
#endif
		
	// hardwire common cases
	switch (e) {
	case 0: _jac_set_identity (o[0]);  return;
	case 1: _jac_set (o[0], a[0]);  return;
	case 2: _jac_square (o[0], a[0], c[0]);  return;
	case 3: _jac_square (b[0], a[0], c[0]); _jac_mult(o[0], a[0], b[0], c[0]);  return;
	case 4: _jac_square (b[0], a[0], c[0]); _jac_square (o[0], b[0], c[0]);  return;
	case 5: _jac_square (b[0], a[0], c[0]); _jac_square (b[0], b[0], c[0]);  _jac_mult(o[0], a[0], b[0], c[0]);  return;
	case 6: _jac_square (b[0], a[0], c[0]); _jac_square (b[1], b[0], c[0]);  _jac_mult(o[0], b[0], b[1], c[0]);  return;
	case 7: _jac_square (b[0], a[0], c[0]); _jac_square (b[1], b[0], c[0]);  _jac_mult(b[1], b[1], b[0], c[0]);  _jac_mult(o[0], b[1], a[0], c[0]);  return;
	case 8: _jac_square (b[0], a[0], c[0]); _jac_square (b[0], b[0], c[0]); _jac_square (o[0], b[0], c[0]);  return;
	case 9: _jac_square (b[0], a[0], c[0]); _jac_square (b[0], b[0], c[0]); _jac_square (b[0], b[0], c[0]);  _jac_mult(o[0], b[0], a[0], c[0]);  return;
	}

	// In theory this uses about the same number of gops as binary exponentiation (about 3n/2), but performs 2/3 of them in parallel at a 3/2 speedup, giving an effective 7n/9 serial gops
	j = ui_len(e)/2;
	e2 = e >> j;
	e1 = e & ((1UL<<j)-1);
	jac_exp2_ui (b, a, e1, e2, c);
	for ( i = 0 ; i < j ; i++ ) _jac_square (b[1], b[1], c[0]);
	_jac_mult (o[0], b[0], b[1], c[0]);
}

// simple simultaneous exponentiation of a common base to two exponents using a 1-bit window size
void jac_exp2_ui (jac_t o[2], jac_t a[1], unsigned long e1, unsigned long e2, curve_t c[1])
{
	register unsigned long m, max;
	jac_t b[2], b01, b10, b11;

	if ( e1 < 2 ) { if ( e1 ) { _jac_set(o[0],a[0]); } else { _jac_set_identity(o[0]); } jac_exp_ui (o+1, a, e2, c);  return; }
	if ( e2 < 2 ) { if ( e2 ) { _jac_set(o[1],a[0]); } else { _jac_set_identity(o[1]); }  jac_exp_ui (o, a, e1, c);  return; }
	max = _ui_max (e1,e2);
	m = 1;
	_jac_set (b[0], a[0]);
	if ( m&e2 ) {
		if ( m&e1) { _jac_set (b11, b[0]);  _jac_set_identity (b10);  _jac_set_identity (b01);  } else {  _jac_set (b10, b[0]);  _jac_set_identity (b11);  _jac_set_identity (b01);  }
	} else {
		if ( m&e1) { _jac_set (b01, b[0]);  _jac_set_identity (b10);  _jac_set_identity (b11);  } else {  _jac_set_identity (b10);  _jac_set_identity (b11);  _jac_set_identity (b01); }
	}
	_jac_square (b[0], b[0], c[0]);
	m <<= 1;
	max >>= 1;
	while ( m <= max ) {
		if ( m&e2 ) {
			if ( m&e1 ) { jac_square_mult (b, &b11, c);  } else {  jac_square_mult (b, &b10, c); }
		} else {
			if ( m&e1 ) { jac_square_mult (b, &b01, c); } else { _jac_square (b[0], b[0], c[0]); }
		}
		m <<= 1;
	}
	if ( m&e2 ) {
		if ( m&e1 ) { _jac_mult (b11, b11, b[0], c[0]); } else { _jac_mult (b10, b10, b[0], c[0]); }
	} else {
		if ( m&e1 ) { _jac_mult (b01, b01, b[0], c[0]); } // else should be impossible
	}
	jac_mult2 (o, &b01, &b11, o+1, &b10, &b11, c);
}

/*
	This function isn't currently used but it could be.  It assumes that b[] contains at least 3 powers of a, indexed by exponent.
	This is useful if you have already computed the powers for some other reason, e.g. baby-steps.
*/
void jac_exp2_powers_ui (jac_t o[2], jac_t a[1], unsigned long e1, unsigned long e2, jac_t b[], unsigned n, int parity, curve_t c[1])
{
	register unsigned long mask, j1, j2;
	register int i, j, k;
	
	// we assume b[0] holds the identity and that powers > 2.
	// we make no attempt to combine operations other than squarings, which should be most of them
	
	k = ui_lg_floor(n);
	mask = (1UL<<k)-1;	
	i = _ui_max (ui_len(e1), ui_len(e2));
	if ( i > k ) {
		j1 = (e1&(mask<<(i-k)))>>(i-k);  j2 = (e2&(mask<<(i-k)))>>(i-k);
	} else {
		j1 = e1;  j2 = e2;
	}
	if ( parity && ! (j1&0x1) ) {											// only use odd values if parity is set
		_jac_set (o[0], b[j1-1]);
		_jac_mult (o[0], o[0], a[0], c[0]);
	} else {
		_jac_set (o[0], b[j1]);
	}
	if ( parity && ! (j2&0x1) ) {											// only use odd values if parity is set
		_jac_set (o[1], b[j2-1]);
		_jac_mult (o[1], o[1], a[0], c[0]);
	} else {
		_jac_set (o[1], b[j2]);
	}
	i -= k;
	while ( i > k ) {
		for ( j = 0 ; j < k ; j++ ) jac_square2 (o, o, o+1, o+1, c);				// this is where the time is spent
		j1 = (e1&(mask<<(i-k)))>>(i-k);
		if ( j1 ) {
			if ( parity && ! (j1&0x1) ) {									// only use odd values if parity is set
				_jac_mult (o[0], o[0], b[j1-1], c[0]);
				_jac_mult (o[0], o[0], a[0], c[0]);
			} else {
				_jac_mult (o[0], o[0], b[j1], c[0]);
			}
		}
		j2 = (e2&(mask<<(i-k)))>>(i-k);
		if ( j2 ) {
			if ( parity && ! (j2&0x1) ) {									// only use odd values if parity is set
				_jac_mult (o[1], o[1], b[j2-1], c[0]);
				_jac_mult (o[1], o[1], a[0], c[0]);
			} else {
				_jac_mult (o[1], o[1], b[j2], c[0]);
			}
		}
		i -= k;
	}
	for ( j = 0 ; j < i ; j++ ) jac_square2 (o, o, o+1, o+1, c);					// and here
	j1 = e1&(mask>>(k-i));
	if ( j1 ) {
		if ( parity && ! (j1&0x1) ) {										// only use odd values if parity is set
			_jac_mult (o[0], o[0], b[j1-1], c[0]);
			_jac_mult (o[0], o[0], a[0], c[0]);
		} else {
			_jac_mult (o[0], o[0], b[j1], c[0]);
		}
	}
	j2 = e2&(mask>>(k-i));
	if ( j2 ) {
		if ( parity && ! (j2&0x1) ) {										// only use odd values if parity is set
			_jac_mult (o[1], o[1], b[j1-1], c[0]);
			_jac_mult (o[1], o[1], a[0], c[0]);
		} else {
			_jac_mult (o[1], o[1], b[j2], c[0]);
		}
	}
}


// quick and dirty binary exponentiation - don't worry about efficiency, this is only used for validation
void jac_exp_mpz (jac_t o[1], jac_t a[1], mpz_t e, curve_t c[1])
{
	register int i, t;
	jac_t b;
	
	if ( ! mpz_sgn(e) ) { _jac_set_identity (o[0]);  return; }

	_jac_set (b, a[0]);			// copy a to handle overlap
	_jac_set (o[0], a[0]) ;
	i = mpz_sizeinbase(e,2) - 1;
	while ( i > 0 ) {
		_jac_square (o[0], o[0], c[0]);
		i--;
		if ( mpz_tstbit(e,i) ) _jac_mult (o[0], o[0], b, c[0]);
	}
}


int jac_verify_group_exponent (mpz_t e, curve_t c[1])
{
	jac_t a, b;
	int i;

	for ( i = 0 ; i < JAC_CONFIDENCE_LEVEL ; i++ ) {
		_jac_random (a, c[0]);
		jac_exp_mpz (&b, &a, e, c);
		if ( ! _jac_is_identity (b) ) return 0;
	}
	return 1;
}


/*
	Computess the order of a given that a^{p^h} = 1
*/
unsigned long jac_pp_order_ui (jac_t a[1], unsigned long p, unsigned long h, curve_t c[1])
{
	jac_t x;
	unsigned long o;
	int i;
	
	o = 1;
	if ( _jac_is_identity (a[0]) ) return o;
	_jac_set (x, a[0]);
	for ( i = 1 ; i < h ; i++ ) {
		o *= p;
		jac_exp_ui (&x, &x, p, c);
		if ( _jac_is_identity (x) ) break;
	}
	if ( ! _jac_is_identity (x) ) o *= p;
	return o;
}


/*
	Recursive fastorder algorithm optimized for small exponents (takes advantage of fast dual exponentiation).
	This algorithm is, to my knowledge, new and not written up anywhere - maybe it should be, it's quite fast (AVS).
*/
int jac_fastorder_ui (unsigned long *po, jac_t a[1], unsigned long e, curve_t c[1])
{
	jac_t t[JAC_MAX_FASTORDER_UI_W];
	unsigned long q[JAC_MAX_FASTORDER_UI_W];
	unsigned long h[JAC_MAX_FASTORDER_UI_W];
	unsigned long o, qp, m, n;
	register int i, j, k, w;

	if ( _jac_is_identity (a[0]) ) { *po = 1;  return 1; }
	if ( ! e ) { err_printf ("%7u: zero exponent for non-identity element in fastorder!\n", _ff_p);  exit (0); }
	w = ui_factor (q, h, e);
	if ( w > JAC_MAX_FASTORDER_UI_W ) { err_printf ("w=%d too large for jac_fastorder_ui, exceeded %d\n", w, JAC_MAX_FASTORDER_UI_W);  exit (0); }
	return jac_factored_order_ui (po, a, q, h, w, c);

}

int jac_factored_order_ui (unsigned long *po, jac_t a[1], unsigned long p[], unsigned long h[], int w, curve_t c[1])
{
	jac_t t[JAC_MAX_FASTORDER_UI_W];
	register unsigned long m, n;
	register int i, j, k;
	
	if ( w == 1 ) {
		n = p[0];
		_jac_set (t[0], a[0]);
		for ( i = 1 ; i < h[0] ; i++ ) {
			jac_exp_ui (t, t, p[0], c);
			if ( _jac_is_identity (t[0]) ) break;
			n *= p[0];
		}
		*po = n;
		return 1;
	}
	k = w-1;	// index of biggest prime factor
	for ( i = 1, n = p[k] ; i < h[k] ; i++ ) n *= p[k];
	m = 1;
	for ( i = 0 ; i < k ; i++ )
		for ( j = 1, m*=p[i] ; j < h[i] ; j++ ) m *= p[i];
		
	jac_exp2_ui (t, a, m, n, c); 
	
	if ( _jac_is_identity(t[0]) ) return jac_factored_order_ui (po, a, p, h, k, c);	// don't need biggest prime power at all
	for ( i = 1 ; i < h[k] ; i++ ) {
		jac_exp_ui (t, t, p[k], c);
		if ( _jac_is_identity(t[0]) ) break;
	}
	if ( i < h[k] ) {												// don't need all the powers of biggest prime factor (rare)
		_jac_set (t[1], a[0]);									// adjust t1 to include only the powers required
		for ( j = 0, n=1 ; j < i ; j++ ) { jac_exp_ui (t+1, t+1, p[k], c);  n*= p[k]; }
	}
	if ( _jac_is_identity (t[1]) ) { *po = n;  return 1; }
	if ( ! jac_factored_order_ui (po, t+1, p, h, k, c) ) return 0;
	*po *= n;
	return 1;
}


/*
    This function is not used by smalljac.

    An implementation of Algorithm 9.1 in [SutherlandThesis].  It is slightly out-of-date as it doesn't include
    the refinements of Algorithm 3 in [Sutherland2007].  These don't impact performance significantly,
    but make the algorithm simpler - should update eventually.

    For now don't try to save space - the space required is roughly L/log(L) which isn't bad compared to the baby/giant tables
*/
unsigned long jac_fastorder_powersmooth (mpz_t o, jac_t a[1], unsigned long L, curve_t c[1])
{
	static mpz_t e;
	static int init;
	prime_enum_ctx_t *ctx;
	jac_t *b, temp;
	unsigned *pp;
	unsigned q[JAC_MAX_FACTORS];
	unsigned p[JAC_MAX_FACTORS];
	unsigned h[JAC_MAX_FACTORS];
	unsigned i, j, k, t, w, low, mid, high;
	unsigned long size;
	int sts;
	
	if ( ! init ) { mpz_init(e);  init = 1; }
	
	mpz_set_ui(o,1);
	if ( _jac_is_identity (a[0]) ) return 1;
		
	ctx = fast_prime_enum_start (0, L, 1);
	w = (unsigned) ui_pi_bound(L) + 1;
	b = mem_alloc (w*sizeof(*b));
	pp = mem_alloc (w*sizeof(*pp));
	size = (unsigned long)w*(sizeof(*b)+sizeof(*pp));
	_jac_set (b[1], a[0]);
	pp[0] = 1;
	for ( i = 1 ; i < w-1 ; i++ ) {			// w is larger than needed, but we should never reach it (assuming |a| is really L-powersmooth)
		if ( _jac_is_identity (b[i]) ) break;
		pp[i] = fast_prime_enum (ctx);
		if ( pp[i] > L ) break;
		jac_exp_ui (b+i+1, b+i, pp[i], c);
	}
	fast_prime_enum_end (ctx);
	if ( ! _jac_is_identity (b[i]) ) { sts = 0; goto cleanup; }
	sts = 1;
	t = i-1;
	w = 1;
	q[0] = pp[t];
	mpz_mul_ui (o, o, pp[t]);
	for (;;) {
		if ( t == 1 ) break;
		low = 1;  high = t-1;
		while ( high > low ) {
			mid = (high + low + 1)/2;
			jac_exp_mpz (&temp, b+mid, o, c);
			if ( _jac_is_identity (temp) ) {
				high = mid-1;
			} else {
				low = mid;
			}
		}
		t = low;
		jac_exp_mpz (&temp, b+t, o, c);
		if ( _jac_is_identity(temp) ) { if ( t != 1 ) printf ("breaking with t = %d\n", t);  break; }
		if ( w >= JAC_MAX_FACTORS ) { err_printf ("Excceded JAC_MAX_FACTORS in jac_fastorder_powersmooth!\n");  exit (0); }
		mpz_mul_ui (o, o, pp[t]);
		q[w++] = pp[t];
	}
	
	// We now have a reasonably small factored exponent in q[0]...q[w] - do a standard classical order algorithm to finish
	for ( i = 0 ; i < w ; i++ ) {
		p[i] = ui_pp_base (q[i]);
		for ( h[i] = 1, t = p[i] ; t < q[i] ; h[i]++, t *= p[i] );
	}
	jac_exp_mpz (b, a, o, c);
	if ( ! _jac_is_identity (b[0]) ) { err_printf ("Exponent validation failed in jac_fastorder_powersmooth\n");  exit (0); }
	
	for ( i = 0 ; i < w ; i++ ) {
		mpz_divexact_ui (e, o, q[i]);
		jac_exp_mpz (b+i, a, e, c);
	}
	mpz_set_ui (o, 1);
	for ( i = 0 ; i < w ; i++ ) {
		for ( j = 0 ; j < h[i]-1 ; j++ ) {
			if ( _jac_is_identity (b[i]) ) break;
			jac_exp_ui (b+i, b+i, p[i], c);
		}
		if ( ! _jac_is_identity (b[i]) ) j++;
		for ( k = 0 ; k < j ; k++ ) mpz_mul_ui (o, o, p[i]);
	}
cleanup:
	mem_free (b);
	mem_free (pp);
	return ( sts ? size : 0 );
}
