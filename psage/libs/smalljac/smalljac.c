#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "gmp.h"
#include "mpzutil.h"
#include "ffwrapper.h"
#include "ffpoly.h"
#include "jac.h"
#include "jacorder.h"
#include "hecurve.h"
#include "smalljac.h"
#include "smalljactab.h"
#if SMALLJAC_GENUS > 1
#include "smalljac_g23.h"
#endif
#include "smalljac_internal.h"
#include "pointcount.h"
#include "cstd.h"

/*
    Copyright 2007=2008 Andrew V. Sutherland

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
	Main smalljac module. 
	Still some work to be done in genus > 1 but the interface should remain unchanged.
	Genus 2 support for curves with a single pt at infinity (y^2=f(x) with f deg 5) is essentially complete,
	but still need to add support for curves with zero or two pts at infty (f deg 6).
	
	Lot's of work still to do in genus 3, currently only handles generic curves (minimal endomorphism ring)
	with a single pt at infinity reliably (results are always provably correct, but it may fail frequently for
	exceptional curves).

	Not currently recommend for general use in genus 3 to compute Lpolys, but will compute a_p reliably.
*/

#define SMALLJAC_MAX_CONSTRAINTS		64

int cornacchia (long *x, long *y, long p, long d);				// used to handle curves of the form y^2=x^6+a

static int _smalljac_initted;
void smalljac_init (void)
{
	if ( _smalljac_initted ) return;
	smalljac_table_alloc (SMALLJAC_TABBITS);
	pointcount_init (SMALLJAC_INIT_COUNT_P);
	_ff_wrapper_init();
	_smalljac_initted = 1;
}


int smalljac_Lpoly (long a[], char *curve, unsigned long q, unsigned long flags)
{
	static mpz_t P;
	static int init;
	smalljac_Qcurve *qc;
	unsigned long p;
	long t, t1, tnm1, tn;
	int i, n, h, error;

	if ( ! init ) { smalljac_init();  mpz_init (P);  init = 1; }
	if ( (flags&SMALLJAC_A1_ONLY) && (flags&SMALLJAC_GROUP) ) return SMALLJAC_INVALID_FLAGS;

	qc = smalljac_Qcurve_init (curve, &error);
	if (  ! qc ) return error;
	if ( qc->genus != SMALLJAC_GENUS ) { error = SMALLJAC_WRONG_GENUS;  goto done; }
	error = 0;
	
	// check for prime powers
	mpz_set_ui (P, q);
	h = mpz_pp_base (P, P);
	if ( ! h ) { error = SMALLJAC_INVALID_PP;  goto done; }
	if ( h > 1 ) {
		if ( (flags&SMALLJAC_GROUP) ) { error = SMALLJAC_INVALID_PP;  goto done; }			// group computation only supported for prime fields at the moment
		if ( (flags&SMALLJAC_A1_ONLY) ) { flags &= ~SMALLJAC_A1_ONLY; n = qc->genus; }		// turn off A1_ONLY flag if q is a prime power
	}
	
	p = mpz_get_ui (P);
	if ( p <= SMALLJAC_TINY_P ) {
		n = smalljac_tiny_Lpoly (a, qc, p, flags);
		goto done;
	}
	if ( mpz_divisible_ui_p (qc->D, p) ) { n = 0;  goto done; }
	
	// Precompute delta values for pointcounting if needed
	if ( smalljac_Qcurve_degree(qc) == 6 || p <= SMALLJAC_COUNT_P ) {
		smalljac_Qcurve_init_Deltas(qc);
		pointcount_precompute (qc->Deltas, qc->f, qc->degree);
		qc->flags |= SMALLJAC_QCURVE_DELTA;
	}
	n = smalljac_internal_Lpoly (a, qc, p, flags);
	if ( n <= 0 ) { error = SMALLJAC_INTERNAL_ERROR;  goto done; }
done:
	smalljac_Qcurve_clear (qc);
	if ( h > 1 ) if ( ! smalljac_Lpoly_extend (a, n, p, h) ) error = SMALLJAC_INTERNAL_ERROR;
	return ( error ? error : n);
}


long smalljac_Lpolys (smalljac_Qcurve_t curve, unsigned long start, unsigned long end, unsigned long flags,
				   int (*callback)(smalljac_Qcurve_t curve, unsigned long p, int good, long a[], int n, void *arg), void *arg)
{
	static mpz_t P;
	static int init;
	smalljac_Qcurve *qc;
	prime_enum_ctx_t *ctx;
	unsigned long h[SMALLJAC_MAX_BAD_PRIMES];
	register unsigned long p, pbitmask, pbits;
	int i, filter, good, good_only,  error;
	
	if ( ! init ) { smalljac_init();  mpz_init (P);  init = 1; }
	if ( (flags&SMALLJAC_A1_ONLY) && (flags&SMALLJAC_GROUP) ) return SMALLJAC_INVALID_FLAGS;
	if ( end < start || end > SMALLJAC_MAX_END ) return SMALLJAC_INVALID_INTERVAL;
	qc = (smalljac_Qcurve *)curve;
	
	// Check that curve genus matches compiled genus - this restriction should go away eventually
	if ( qc->genus < 4 ) {
		if ( qc->genus != SMALLJAC_GENUS ) return SMALLJAC_WRONG_GENUS;
	} else {
		if ( ! (flags&SMALLJAC_A1_ONLY) ) return SMALLJAC_WRONG_GENUS;
	}
	
	good_only = (flags&SMALLJAC_GOOD_ONLY);
	filter = (flags&SMALLJAC_FILTER);
	
	// Precompute delta values for pointcounting, if needed
	if ( (smalljac_Qcurve_degree(qc)==6||start <= SMALLJAC_COUNT_P) && ! (qc->flags&SMALLJAC_QCURVE_DELTA) ) {
		smalljac_Qcurve_init_Deltas(qc);
		pointcount_precompute (qc->Deltas, qc->f, qc->degree);
		qc->flags |= SMALLJAC_QCURVE_DELTA;
	}
	
	pbitmask = (flags&SMALLJAC_SPLIT)>>SMALLJAC_SPLIT_SHIFT;
	pbits = (flags&SMALLJAC_HIGH)>>SMALLJAC_HIGH_SHIFT;
	// First deal with tiny primes
	if ( start <= SMALLJAC_TINY_P ) {
		if ( ! start ) p = 2; else p = ui_small_prime(ui_small_prime_index(start-1)+1);			// set start to first prime >= start
		for ( ; p <= SMALLJAC_TINY_P ; p = ui_small_prime(ui_small_prime_index(p)+1) ) {
			if ( p > end ) return end;
			if ( (p&pbitmask) != pbits ) continue;
			if ( filter && ! (*callback) (curve, p, -1, 0, 0, arg) ) continue;
			qc->n = smalljac_tiny_Lpoly (qc->a, qc, qc->p=p, flags);						// note tiny_Lpoly handles bad reduction check
			if ( qc->n < 0 ) return SMALLJAC_INTERNAL_ERROR;
			if ( qc->n > 0 || ! good_only ) {
				if ( ! (*callback) (curve, p, (qc->n > 0 ? 1 : 0), qc->a, qc->n, arg) ) return p;
			}
		}
		if ( p > end ) return end;
		start = p;
	}
	error = 0;
	// use fast prime enumeration (based on a wheeled sieve) if possible - we should extend this to cover all cases
	// the code duplication below is intentional
	if ( end <= MPZ_MAX_ENUM_PRIME ) {
		ctx = fast_prime_enum_start (start, end, 0);
		while ( (p = fast_prime_enum(ctx)) ) {
			if ( (p&pbitmask) != pbits ) continue;
			if ( filter && ! (*callback) (curve, p, -1, 0, 0, arg) ) continue;
			if ( mpz_divisible_ui_p (qc->D,p) ) {
				if ( ! good_only ) if ( ! (*callback) (curve, p, 0, 0, 0, arg) ) break;
				continue;
			}
			qc->n = smalljac_internal_Lpoly (qc->a, qc, qc->p=p, flags);
			if ( qc->n <= 0 ) { error = 1;  break; }
			if ( ! (*callback) (curve, p, 1, qc->a, qc->n, arg) ) break;
		}
		fast_prime_enum_end (ctx);
		if ( ! p ) p = end;
	} else {
		// using nextprime here is a very slow -- not an issue in genus >1, but in genus 1 we could do better by extending fast prime enum code
		mpz_set_ui (P, start-1);
		for ( mpz_nextprime (P,P) ; (p = mpz_get_ui(P)) <= end ; mpz_nextprime (P,P) ) {
			if ( (p&pbitmask) != pbits ) continue;
			if ( filter && ! (*callback) (curve, p, -1, 0, 0, arg) ) continue;
			if ( mpz_divisible_ui_p (qc->D,p) ) {
				if ( ! good_only ) if ( ! (*callback) (curve, p, 0, 0, 0, arg) ) break;
				continue;
			}
			qc->n = smalljac_internal_Lpoly (qc->a, qc, qc->p=p, flags);
			if ( qc->n <= 0 ) { error = 1;  break; }
			if ( ! (*callback) (curve, p, 1, qc->a, qc->n, arg) ) break;
		}
		if ( p > end ) p = end;
	}
	if ( error ) { printf ("smalljac internal error at p=%lu\n", p);  return SMALLJAC_INTERNAL_ERROR; }
	return (long)p;
}

/*
	smalljac_internal_Lpoly assumes p > qc->degree and good reduction.
        It computes the coefficients of L_p(T) for the curve specified by F/denom and degree.
	The array D0 holds precomputed delta values used for pointcounting which will be
	used to determine a1 if p < SMALLJAC_COUNT_P.
	
	If p >= SMALLJAC_PADIC_P, a p-adic computation will be used to compute the 
	coefficients of L_p(T) modulo p (or p^N in genus > 3), using Harvey's frobenius
	code, followed by a generic group computation to get the exact function.
	
	Otherwise, a generic computation is performed, possibly using the a1 value
	obtained from pointcounting.
	
	The return value is either 0 (error) or the number of coefficients computed.
	If a1_only is set AND if pointcounting is performed (not necessarily the case) then
	this will be 1, otherwise it will be g.   Note that for large p it may be more efficient
	to compute a1 via a method that computes all the coefficients rather than
	pointcounting, in this case all the coefficients are returned (may as well...).
*/

int smalljac_internal_Lpoly (long a[], smalljac_Qcurve *qc, unsigned long p, unsigned long flags)
{
	curve_t c;

	qc->pts = 0;
	if ( qc->genus >2 && (flags&SMALLJAC_A1_ONLY) ) {
		qc->pts = smalljac_pointcount (qc, p);
		if ( ! qc->pts ) return 0;
		a[0] = (long)qc->pts - (long)(p+1); 
		return 1;
	}
// avoid requiring p-adic code to be compiled in genus < 3
#if SMALLJAC_GENUS > 2
	if ( p >= SMALLJAC_PADIC_P ) return smalljac_padic_Lpoly (a, qc, p, flags);
#endif
	
	// invoke special code for y^2=x^6+a curves
	if ( !(flags&SMALLJAC_GROUP) && smalljac_Qcurve_special(qc) == SMALLJAC_SPECIAL_X6PA ) return smalljac_x6pa_Lpoly (a, qc, p, flags);
	
	if ( smalljac_Qcurve_degree(qc)==6 || p <= SMALLJAC_COUNT_P ) {
		qc->pts = smalljac_pointcount (qc, p);
		if ( ! qc->pts ) { printf ("smalljac_pointcount failed at p=%lu\n", p); return 0; }
		if ( (flags&SMALLJAC_A1_ONLY) ) { a[0] = (long)qc->pts - (long)(p+1);  return 1; }
		if ( qc->genus == 1 ) {										// in genus 1, we're basically done, but we might need to do a group structure computation
			// this code is used in genus 1 only for very small p,  in fact if SMALLJAC_COUNT_P <= SMALLJAC_TINY_P we never get here
			if ( ! (flags&SMALLJAC_GROUP) ) { a[0] = (long)qc->pts - (long)(p+1);  return 1; }
			ff_setup_ui (p);
			if ( ! ff_poly_set_rational_mpz (c.f, qc->f, (c.d=qc->degree), qc->denom)  ) return 0;
			return hecurve_g1_group_structure (a,qc->pts,0,c.f);			// we don't have any torsion info so specify d=0
		}
	}
	return smalljac_generic_Lpoly (a, qc, p, qc->pts, flags);
}


unsigned long smalljac_pointcount (smalljac_Qcurve *qc,  unsigned long p)
{
	unsigned long Deltaf0[SMALLJAC_MAX_DEGREE+1];

	if ( ! ui_poly_set_rational_mpz_mod_p (Deltaf0, qc->Deltas, qc->degree, qc->denom, p)  ) return 0;
	if ( p < SMALLJAC_BIG_COUNT_P ) {
		switch (qc->degree) {
		case 3: return pointcount_g1 (Deltaf0,p);
		case 5: return pointcount_g2 (Deltaf0,p);
		case 6:
			if ( mpz_cmp_ui(qc->denom,1) != 0 ) { err_printf ("Non-integral coefficients in degree 6 poly unexpected, denom=%Zd\n",qc->denom);  return 0; }
			return pointcount_g2d6 (Deltaf0,p,mpz_fdiv_ui(qc->f[6],p));
		case 7: return pointcount_g3 (Deltaf0,p);
		case 9: return pointcount_g4 (Deltaf0,p);
		}
	} else {
		switch (qc->degree) {
		case 3: return pointcount_big_g1 (Deltaf0,p);
		case 5: return pointcount_big_g2 (Deltaf0,p);
		case 6:
			if ( mpz_cmp_ui(qc->denom,1) != 0 ) { printf ("Non-integral coefficients in degree 6 poly unexpected\n");  return 0; }
			return pointcount_big_g2d6 (Deltaf0,p,mpz_fdiv_ui(qc->f[6],p));
		case 7: return pointcount_big_g3 (Deltaf0,p);
		case 9: return pointcount_big_g4 (Deltaf0,p);
		}
	}
	return 0;
}

// Computes L-poly coefficients or group structure (or both) using generic group algorithms
// Returns 0 for errors, otherwise the number of entries in a which will be g for L-poly coefficients, or the rank of the group.
int smalljac_generic_Lpoly (long a[], smalljac_Qcurve *qc, unsigned long p, unsigned long pts, unsigned long flags)
{
	curve_t c, twist;
	unsigned long Min, Max;
	unsigned long e, te, o, P1, PN1;
	unsigned long m, gops;
	long a1, d, tcon;
	int constraints[SMALLJAC_MAX_CONSTRAINTS+1];
	int i, j, k, tk, n, sts;
	double x, t;

	ff_setup_ui (p);
	if ( ! ff_poly_set_rational_mpz (c.f, qc->f, (c.d=qc->degree), qc->denom)  ) return 0;
#if SMALLJAC_GENUS == 1
	// Use faster genus 1 code in hecurve1x instead of general jac functions
	if ( pts ) { a[0] = (long)pts - (long)(p+1);  return qc->genus; }
	P1 = hecurve_g1_order (&d,c.f);
	if ( (flags&SMALLJAC_GROUP) ) {
		return hecurve_g1_group_structure (a,P1,d,c.f);
	} else {
		a[0] = (long)(P1) - (long)(p+1);
		return 1;
	}
#endif
	x = sqrt((double)p);
	Min = (unsigned long) (floor(pow(x-1.0, 2.0*SMALLJAC_GENUS)));
	Max = (unsigned long) (ceil(pow(x+1.0, 2.0*SMALLJAC_GENUS)));
	if ( pts ) {
		a1 = (long)pts - (long)(p+1);
	} else {
		a1 = JAC_INVALID_A1;
	}
	k = jac_order (&P1, Min, Max, a1, 0, &c, flags&SMALLJAC_SGROUP);
	if ( k <= 0 ) return 0;
	if ( k > SMALLJAC_MAX_CONSTRAINTS ) { printf ("%7lu: Exceeded SMALLJAC_MAX_CONSTRAINTS", p);  return 0; }
#if SMALLJAC_GENUS == 2
	if ( k > 1 ) { printf ("%lu: Ambiguous result in genus 2 not handled\n", p);  return 0; }		// should be impossible provided pointcounting is used for p < SMALLJAC_TINY_P
	if ( (flags&SMALLJAC_GROUP) ) return jac_structure (a, &c, P1, flags&SMALLJAC_SGROUP);
	if ( pts ) {
		a[0] = (long)pts - (long)(p+1);
		a[1] = (long)(P1) - ((long)p*p+1)-(long)(p+1)*a[0];
		return qc->genus;
	} else {
		if ( ! smalljac_genus2_charpoly_from_P1 (a, P1, Min, Max, &c) ) {
			// If we can't deduce the charpoly from P1 (very rare) go ahead and compute PN1 from scratch
			// This is a bit silly, since we know the value of PN1 mod 2(p+1) and could use this knowledge to speed things up,
			// but it happens so rarely that we don't bother
			ff_poly_twist (twist.f, c.f, c.d);
			twist.d = c.d;
			k = jac_order (&PN1, Min, Max, JAC_INVALID_A1, 0, &twist, 0);
			if ( k <= 0 ) { printf ("%lu: Attempted twist order computation failed\n", p);  return 0; }
			if ( k > 1 ) { printf ("%lu: Ambiguous result in genus 2 not handled\n", p);  return 0; }		// should be impossible
			if ( ! smalljac_genus2_charpoly (a, p, x, P1, PN1) ) { err_printf ("%lu: smalljac_genus2_charpoly failed\n", p);  return 0; }
			return qc->genus;
		}
		return qc->genus;
	}
#endif
#if SMALLJAC_GENUS == 3
	if ( k == 1 ) {
		if ( (flags&SMALLJAC_GROUP) ) return jac_structure (a, &c, P1, flags&SMALLJAC_SGROUP);
		return ( smalljac_genus3_charpoly_from_P1 (a, P1, pts, Min, Max, &c) ? qc->genus : 0 );
	} else {
		unsigned long unique_P1, unique_PN1;
		
		ff_poly_twist (twist.f, c.f, c.d);
		twist.d = c.d;
		m = 2*(p+1);
		e = P1;
		P1 = _ui_ceil_ratio(Min,e)*e;
		for ( i = 0 ; i < k ; i++ ) {
			tcon = 2*(p*p*p+1) - P1;
			constraints[i] = (tcon < 0 ? (int) (m - ((-tcon)%m)) : (int) (tcon%m) );
			P1 += e;
		}
		constraints[k] = -1;
		tk = jac_order (&PN1, Min, Max, (pts ? -a1 : JAC_INVALID_A1), constraints, &twist, 0);
		if ( tk <= 0 ) { printf ("%lu: Attempted twist order computation failed\n", p);  return 0; }
		P1 = _ui_ceil_ratio(Min,e)*e;
		te = PN1;
		n = 0;
		for ( i = 0 ; i < k ; i++ ) {
			PN1 = _ui_ceil_ratio(Min,te)*te;
			for ( j = 0 ; j < tk ; j++ ) {
				sts = smalljac_genus3_charpoly (a, p, x, P1, PN1, pts);
				if ( sts ) {
					if ( ! n ) { unique_P1 = P1;  unique_PN1 = PN1; }
					n++;
				}
				PN1 += te;
			}
			P1 += e;
		}
		if ( n != 1 ) { printf ("%lu: Unable to resolve ambiguity, primary subgroup %lu, twist subgroup %lu, pts = %lu, n = %d, giving up\n", p, e, te, pts, n);  return 0; }
		if ( (flags&SMALLJAC_GROUP) ) return jac_structure (a, &c, unique_P1, 0);
		return ( smalljac_genus3_charpoly (a, p, x, unique_P1, unique_PN1, pts) ? qc->genus : 0 );
	}
#endif
	err_printf ("Fell through to unreachable code - compiled genus %d\n", SMALLJAC_GENUS);  exit (0);
}

#if SMALLJAC_GENUS > 2
// Computes Lpoly coefficients or group structure using p-adic computations to determine L_p(T) mod p
// by calling David Harvey's frobenius() function.  If group structure is requested, uses #J(C/F) = L_p(1)
// to compute the group order.  Returns 0 if an error occurs, the number of entries in a otherwise (either g or the rank of the group)
int smalljac_padic_Lpoly (long a[], smalljac_Qcurve *qc, unsigned long p, unsigned long flags)
{
	curve_t c;
	long f[SMALLJAC_MAX_DEGREE+1];
	int i;
	
	ff_setup_ui (p);
	if ( ! ff_poly_set_rational_mpz (c.f, qc->f, qc->degree, qc->denom)  ) return 0;
	c.d = qc->degree;
	// this is slightly annoying, but we need to convert to standard integers (mod p) and ff uses Montogemery
	for ( i = 0 ; i <= qc->degree ; i++ ) f[i] = _ff_get_ui (c.f[i]);
	if ( ! padic_charpoly (a, f, qc->degree, p) ) { err_printf ("%lu: padic_charpoly failed\n", p);  return 0; }
	switch ( qc->genus) {
	case 2:
		if ( ! smalljac_genus2_charpoly_from_Pmodp (a, &c) ) return 0;
		break;
	case 3:
		if ( ! smalljac_genus3_charpoly_from_Pmodp (a, &c) ) return 0;
		break;
	default:
		err_printf ("Invalid genus %d in smalljac_padic_lpoly\n", qc->genus);
		exit (0);
	}
	if ( (flags&SMALLJAC_GROUP) ) {
		return jac_structure (a, &c, smalljac_Lp1_ui (a, qc->genus, p), flags&SMALLJAC_SGROUP);
	} else {
		return qc->genus;
	}
}
#endif

// table of Lpoly a1 values for all possible genus 1 Q curve reductions over F_2.  -9 denotes singular cases
// table is indexed by n = a1+2*a2+4*a3+8*a4+16*a6 with a_i reduced mod 2.
char ws2a1tab[32] = { -9, -9, -9, -9, 0, -9, 2, 1,
	                          -9, 1, -9, -1, 2, -9, 0, 1,
	                          -9, 1, -9, -1, 0, -1, -2, -9,
	                          -9, -9, -9, -9, -2, -1, 0, -9 };

// Handle general Weierstrass case for p=3
int smalljac_Lpoly3_ws (long a[], smalljac_Qcurve *qc)
{
	register int a1, a2, a3, a4, a6, x, y, z, f0, f1, f2, f3;
	register unsigned ws;
	register int i, pts;

	pts = 1;
	ws = qc->ws3;
	a1 = ws&0x3;  ws >>= 2;  a2 = ws&0x3;  ws >>= 2;  a3 = ws&0x3; ws >>= 2;  a4 = ws&0x3; ws >>= 2;  a6 = ws&0x3;
	f0 = a2+a1*a1;   f1 = a4+2*a1*a3;  f2 = a6+a3*a3;			// f0 = c2', f1 = c4', f2 = c6' after completing the square in char 3
	z = f0*f0*f1*f1 - f0*f0*f0*f2 - f1*f1*f1;						// discriminant of transformed poly (in char 3)
	if ( z%3 == 0 ) return 0;
	for ( x = 0 ; x <= 2 ; x++ ) {
		for ( y = 0 ; y <= 2 ; y++ ) {
			z = x*x*x + a2*x*x + a4*x +a6 - y*y - a1*x*y - a3*y;
			if ( z%3 == 0 ) pts++;
		}
	}
	a[0] = pts - 4;
	return 1;
}

// Computes either Lpoly coefficients or group structure for tiny values of p.  Returns 0 for error, number of entries in a otherwise
int smalljac_tiny_Lpoly (long a[], smalljac_Qcurve *qc, int p, unsigned long flags)
{
	unsigned long f[SMALLJAC_MAX_DEGREE+1];
	curve_t c;
	long P1,pts;
	int i;

	// handle special cases for 2 and 3 first
	if ( p == 2 ) {
		// note that except for genus 1 curves in general Weierstrass form, we only support curves y^2 = f(x) which all have bad reduction at 2
		if ( ! (qc->flags&SMALLJAC_QCURVE_WS) ) return 0;
		if ( ws2a1tab[qc->ws2] < -2 ) return 0;
		a[0] = ws2a1tab[qc->ws2];
		if ( (flags&SMALLJAC_GROUP) ) a[0] += 3;		// Happily, the Jacobian over F_2 is always cyclic
		return 1;
	} else if ( p == 3 && (qc->flags&SMALLJAC_QCURVE_WS) ) {
		// there are 8 reduced general Weierstrass forms with non-cyclic groups over F_3 (Z_2 x Z_2 is the only possibility)
		// we handle these cases explicitly, and otherwise assume the group is cyclic.
		if ( (flags&SMALLJAC_GROUP) ) {
			switch (qc->ws3) {
			case 0x22a:	// [2,2,2,0,2]		note the encoding is a1 in bits 0-1, a2 in bits 2-3, ..., a6 in bits 8-9.
			case 0x219:	// [1,2,1,0,2]		
			case 0x269:	// [1,2,2,1,2]		
			case 0x25a:	// [2,2,1,1,2]		
			case 0x8a: 	// [2,2,0,2,0]		
			case 0x89:	// [1,2,0,2,0]		
			case 0x290:	// [0,0,1,2,2]		
			case 0x2a0:	// [0,0,2,2,2]		
			case 0x80:	// [0,0,0,2,0]		
			a[0] = a[1] = 2;  return 2;
			}
		}
		if ( ! smalljac_Lpoly3_ws (a, qc) ) return 0;
		if ( (flags&SMALLJAC_GROUP) ) a[0] += 4;
		return 1;
	}
	// Check reduction at p - we need to be careful here, since we may have modified the curve to clear the 2g coefficient
	// leaving the curve undefined over F_p if p divides (2g+1), even though the original curve may have had good reduction at p
	if ( /*(qc->flags & SMALLJAC_QCURVE_2G) &&*/ !(qc->degree%p) ) {					// actually, we may as well always do this and not worry about setting SMALLJAC_QCURVE_2G
		if ( mpz_divisible_ui_p(qc->disc, p) ) return 0;
		if ( ! ui_poly_expr_mod_p (f, qc->degree, qc->str, p) ) return 0;					// use original curve string reduced mod p
	} else {
		if ( mpz_divisible_ui_p(qc->D, p) ) return 0;
		if ( ! ui_poly_set_rational_mpz_mod_p (f, qc->f, qc->degree, qc->denom, p) ) return 0;
	}
	qc->pts = pointcount_tiny (f, qc->degree, p);
	if ( qc->genus == 1 ) {
		if ( (flags&SMALLJAC_GROUP) ) {
			ff_setup_ui (p);
			for ( i = 0, c.d = qc->degree ; i <= c.d ; i++ ) _ff_set_ui (c.f[i], f[i]);
			return hecurve_g1_group_structure (a,qc->pts,0,c.f);						// we don't have any torsion info so specify d=0
		}
		a[0] = qc->pts-p-1;
		return 1;
	}
	a[0] = qc->pts-p-1;
	if ( (flags&SMALLJAC_A1_ONLY) ) return 1;
	pts = pointcount_tiny_d2 (f, qc->degree, p);
	a[1] = (pts - p*p - 1+a[0]*a[0]) / 2;
	if ( qc->genus == 2 ) {
		if ( (flags&SMALLJAC_GROUP) ) {
			P1 = smalljac_Lp1_ui (a, qc->genus, p);
			ff_setup_ui (p);
			for ( i = 0, c.d = qc->degree ; i <= c.d ; i++ ) _ff_set_ui (c.f[i], f[i]);
			return jac_structure (a, &c, P1, flags&SMALLJAC_SGROUP);
		} else {
			return 2;
		}
	}

	pts = pointcount_tiny_d3 (f, qc->degree, p);
	a[2] = (pts - (long)p*p*p - 1 - a[0]*a[0]*a[0] + 3*a[0]*a[1]) / 3;
	if ( qc->genus == 3 ) {
		if ( (flags&SMALLJAC_GROUP) ) {
			P1 = smalljac_Lp1_ui (a, qc->genus, p);
			ff_setup_ui (p);
			for ( i = 0, c.d = qc->degree ; i <= c.d ; i++ ) _ff_set_ui (c.f[i], f[i]);
			return jac_structure (a, &c, P1, flags&SMALLJAC_SGROUP);
		} else {
			return 3;
		}
	}
	err_printf ("Unhandled genus in smalljac_Lpoly_tiny for p=%d\n", p);
	return 0;
}


// Computes either Lpoly coefficients for curves y^2=x^6+a, doesn't compute group structure
int smalljac_x6pa_Lpoly (long a[], smalljac_Qcurve *qc, long p, unsigned long flags)
{
	long x, y, k, m, n, ap, a2;
	int i, split;
	ff_t s, t, u;
	ff_t f[7];

	// f(x) must be of the form x^6+a
	if ( qc->degree != 6 ) { err_printf ("Invalid curve degree in smalljac_x6pa_Lpoly\n");  return 0; }
	if ( mpz_cmp_ui(qc->f[6],1) != 0 )  { err_printf ("Invalid curve in smalljac_x6pa_Lpoly\n");  return 0; }
	for ( i = 1 ; i < 6 ; i++ ) if ( mpz_sgn(qc->f[i]) )  { err_printf ("Invalid curve in smalljac_x6pa_Lpoly\n");  return 0; }
	
	// If p is not 1 mod 3, the Hasse-witt matrix is 0, so a_1=a_2=0 mod p, and the curve is supersingular
	// by Thm 6.1 of Galbraith "Supersingular curves in cryptography", and the L-poly must be (pT^2+1)^2.
	if ( p%3 != 1 ) { a[0] = 0;  a[1] = 2*p;  return 2; }
	
	// p = 6m+1
	// the Hasse Witt matrix is diagonal, with entries binom(3m,m)a^{2m} and binom(3m,m)a^m
	// we compute (3m,m) by applying formula 9.1 (p. 453) of "Binomial coefficients and Jacobi Sums"
	// by Hudson and Williams.
	ff_setup_ui (p);
	if ( ! cornacchia (&x, &y, p, 3) ) { printf ("%d is not a prime of the form x^2 + ny^2\n", p);  return 0; }
	if ( x*x+3*y*y != p ) { printf ("cornacchia failed: %d != %d^2 +3*d^2\n", p, x, y);  return 0; }
	if ( (x%3) != 1 ) x = -x;	// make sure x is 1 mod 3
	// binom(3m,m) = 2x
	_ff_set_i (u, mpz_get_i(qc->f[0]));
	k = (p-1)/6;
	ff_exp_ui (&s, &u, k);
	_ff_square (t, s);
	_ff_set_i (u, x);
	_ff_x2 (u);
	_ff_mult(s,s,u);
	_ff_mult(t,t,u);
	_ff_add(u,s,t);
	ap = _ff_get_ui(u);
	_ff_mult(u,s,t);
	a2 = _ff_get_ui(u);
	if ( ap > p/2 ) ap -= p;
	// We now know the L-poly mod p, we just need to nail a_2

	m = 2*p + (ap*ap)/4;
	if ( m > a2 ) m -= (m-a2)%p; else m -= p - ((a2-m)%p);
	n = p*p+1-(p+1)*ap;
	if ( mpz_cmp_ui(qc->f[0],1)==0 ) {
		k=0;
	} else {
		_ff_set_one (f[6]);  _ff_set_mpz(f[0],qc->f[0]);
		for ( i = 1 ; i < 6 ; i++ ) _ff_set_zero (f[i]);
		k = ( ff_poly_factors (f, 6) > 2 ? 0 : 3);		// NOTE: 3-3 split does NOT yield 2-torsion.  1-5 and 2-4 splits can't happen
	}
	
	// We now know Lp(1)=k mod  6 which uniquely determines a_2
	
	while ( (m+n-k)%6 ) m -= p;
	a[0] = -ap;
	a[1] = m;
	return 2;
}



smalljac_Qcurve *smalljac_Qcurve_alloc (void)
{
	smalljac_Qcurve *qc;
	
	qc = mem_alloc (sizeof(*qc));
	mpz_init (qc->disc);
	mpz_init (qc->denom);
	mpz_init (qc->D);
	return qc;
}


smalljac_Qcurve_t smalljac_Qcurve_init (char *str, int *err)
{
	smalljac_Qcurve *qc;

	qc = smalljac_Qcurve_alloc();
	if ( ! smalljac_Qcurve_set_str (qc, str, err) ) { smalljac_Qcurve_clear ((smalljac_Qcurve_t)qc);  return (smalljac_Qcurve_t)0; }
	return (smalljac_Qcurve_t) qc;
}


int smalljac_Qcurve_set_str (smalljac_Qcurve *qc, char *str, int *err)
{
	static mpz_t W[5];
	static mpq_t F[SMALLJAC_MAX_DEGREE+1], DQ;
	static int init;
	int i, degree, ws;
	
	if ( ! init ) { smalljac_init();  for ( i = 0 ; i < 5 ; i++ ) mpz_init (W[i]);  for ( i = 0 ; i <= SMALLJAC_MAX_DEGREE ; i++ ) mpq_init (F[i]);  mpq_init(DQ);  init = 1; }
	if ( err ) *err = 0;
	if ( strlen(str)+1 > sizeof(qc->str) ) { if ( err ) *err = SMALLJAC_PARSE_ERROR;  return 0; }
	strcpy (qc->str, str);

	degree = poly_expr_degree (str, &ws);
	if ( degree > SMALLJAC_MAX_DEGREE ) { if ( err ) *err = SMALLJAC_UNSUPPORTED_CURVE;  return 0; }

	qc->flags = qc->ws2 = qc->ws3 = 0;

	if ( ws ) {
		if ( gmp_sscanf (str, "[%Zd,%Zd,%Zd,%Zd,%Zd]", W[0], W[1], W[2], W[3], W[4]) != 5 ) { if ( err ) *err = SMALLJAC_PARSE_ERROR;  return 0; }
		mpq_poly_weierstrass (F, W);
		for ( i = 4 ; i >= 0 ; i-- ) {
			if ( i < 4 ) { qc->ws2 <<= 1;  qc->ws3 <<= 2; }
			qc->ws2 |= mpz_tstbit(W[i],0);
			qc->ws3 |= mpz_fdiv_ui (W[i],3);
		}
		qc->flags |= SMALLJAC_QCURVE_WS;
	} else {
		if ( mpq_poly_expr (F, degree, str) != degree ) { if ( err ) *err = SMALLJAC_PARSE_ERROR;  return 0; }
		if ( (degree&1) && mpq_poly_standardize (F, degree) ) qc->flags |= SMALLJAC_QCURVE_2G;			// note if curve was modified to make 2g coefficient zero - only used for odd degree
	}

	// Filter out unsupported curves here rather than above so we can be a little more informative about parsing errors
	if ( degree < 3 || (degree!=6 && ! (degree&1)) ) { if ( err ) *err = SMALLJAC_UNSUPPORTED_CURVE;  return 0; }		// degree 6 is currently only supported even degree
	smalljac_Qcurve_init_f (qc, degree);
	mpq_poly_discriminant (DQ, F, degree);
	if ( ! mpq_sgn (DQ) )  {if ( err ) *err = SMALLJAC_SINGULAR_CURVE;  return 0; }
	mpz_poly_set_mpq (qc->denom, qc->f, F, degree);
	mpz_set (qc->disc, mpq_numref(DQ));
	mpz_mul (qc->D, qc->denom, qc->disc);		// could use lcm here
	if ( degree == 6 ) {
		if ( mpz_cmp_ui(qc->f[6],1)==0 ) {
			for ( i = 1 ; i < 6 ; i++ ) if ( mpz_sgn(qc->f[i]) ) break;
			if ( i == 6 ) qc->special = SMALLJAC_SPECIAL_X6PA;
		}
	}
	return 1;
}

char *smalljac_Qcurve_str (smalljac_Qcurve_t curve) { return ((smalljac_Qcurve*)curve)->str; }
int smalljac_Qcurve_genus (smalljac_Qcurve_t curve) { return ((smalljac_Qcurve*)curve)->genus; }

int smalljac_Qcurve_set_mpz (smalljac_Qcurve *qc, mpz_t f[], int degree, mpz_t denom, mpz_t disc, char *str)
{
	char *t;
	int i;
	
	if ( degree < 3 || degree > SMALLJAC_MAX_DEGREE || (! (degree&1)&&degree!=6) || ! mpz_sgn(denom) || ! mpz_sgn(disc) ) return 0;
	qc->flags = qc->ws2 = qc->ws3 = 0;
	smalljac_Qcurve_init_f (qc, degree);
	for ( i = 0 ; i <= degree ; i++ ) mpz_set (qc->f[i], f[i]);
	mpz_set (qc->disc, disc);
	mpz_set (qc->denom, denom);
	mpz_mul (qc->D, disc, denom);
	if ( degree == 6 ) {
		if ( mpz_cmp_ui(f[6],1)==0 ) {
			for ( i = 1 ; i < 6 ; i++ ) if ( mpz_sgn(f[i]) ) break;
			if ( i == 6 ) qc->special = SMALLJAC_SPECIAL_X6PA;
		}
	}
	// trust that str actually matches the curve - the idea is to avoid recomputing the discriminant - may not be worth the trouble...
	strcpy (qc->str, str);
	return 1;
}

int smalljac_Qcurve_set_i (smalljac_Qcurve *qc, long f[], int degree, char *str)
{
	char *t;
	int i;
	
	if ( degree < 3 || degree > SMALLJAC_MAX_DEGREE || (! (degree&1)&&degree!=6)  ) return 0;
	qc->flags = qc->ws2 = qc->ws3 = 0;
	smalljac_Qcurve_init_f (qc, degree);
	for ( i = 0 ; i <= degree ; i++ ) mpz_set_i (qc->f[i], f[i]);
	mpz_set_ui (qc->denom, 1);
	mpz_poly_discriminant (qc->disc, qc->f, degree);
	mpz_mul (qc->D, qc->disc, qc->denom);
	if ( degree == 6 ) {
		if ( mpz_cmp_ui(qc->f[6],1)==0 ) {
			for ( i = 1 ; i < 6 ; i++ ) if ( mpz_sgn(qc->f[i]) ) break;
			if ( i == 6 ) qc->special = SMALLJAC_SPECIAL_X6PA;
		}
	}
	strcpy (qc->str, str);
	return 1;
}

void smalljac_Qcurve_clear (smalljac_Qcurve_t curve)
{
	smalljac_Qcurve *qc;
	int i;
	
	qc = (smalljac_Qcurve *)curve;
	mpz_clear (qc->D);
	mpz_clear (qc->denom);
	for ( i = 0 ; i < qc->f_inits ; i++ ) mpz_clear (qc->f[i]);
	for ( i = 0 ; i < qc->Delta_inits ; i++ ) mpz_clear (qc->Deltas[i]);
	mem_free (qc);
}


int smalljac_Lpoly_extend (long a[], int n, unsigned long p, int h)
{
	register long t, t1, tn, tnm1;
	register int i;

	if ( h < 1 ) return 0;
	if ( h == 1 ) return 1;
	switch (n) {
	case 1:
		tnm1 = t1 = -a[0];  tn = t1*t1-2*(long)p; 
		for ( i = 2 ; i < h ; i++ ) { t = t1*tn - p*tnm1;  tnm1 = tn;  tn = t; }
		a[0] = -tn;
		break;
	default:
		return 0;		// add support for genus > 1 later
	}
	return 1;
}

// implementation of Cornacchia's algorithm to find (x,y) s.t. x^2+dy^2 = p  see Cohen, Alg. 1.5.2
// returns 0 if no solution exists.
int cornacchia (long *x, long *y, long p, long d)
{
	ff_t k, t;
	long a, b, c, r, L;
	
	_ff_set_i (t, -d);
	if ( ! ff_sqrt (&k, &t) ) return 0;
	b = _ff_get_ui(k);
	if ( b < (p>>1) ) b = p-b;
	a = p;
	L = sqrt(p);
	while ( b > L ) { r = a%b;  a=b;  b=r; }
	r = p-b*b;
	c = r/d;
	if ( c*d != r ) return 0;
	if ( (*y=ui_sqrt(c)) < 0 ) return 0;
	*x=b;
	return 1;
}
