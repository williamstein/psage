#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <stdint.h>
#include <memory.h>
#include <math.h>
#include "gmp.h"
#include "mpzutil.h"
#include "ffwrapper.h"
#include "ffpoly.h"
#include "hecurve.h"
#include "jac.h"
#include "smalljac.h"
#include "smalljactab.h"
#include "jacorder.h"
#include "cstd.h"

/*
    Copyright 2007-2008 Andrew V. Sutherland

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
	
	This code is still subject to occasional tweaking, so some debugging code has been
	left in the source (but commented out).
	
	6/1/2008: As of version 3.0, genus 1 order and structure computations are now handled in hecurve1x.c,
	so none of the code here gets used in genus 1 anymore, but for the moment we won't change things.
*/

// SMALLJAC_M determines the number of group operations to perform simultaneously
// when doing a parallel search
#if HECURVE_GENUS == 1
unsigned long SMALLJAC_M = 4;				// increase to 8 for large p
#else
unsigned long SMALLJAC_M = 32;
#endif

#define SMALLJAC_BABY_STASHSIZE		4096
jac_t baby_stash[SMALLJAC_BABY_STASHSIZE];

/*
	The tables below are used to select parameters for BSGS searches in genus 2 and 3
	when a1 is known.  They contain the median and expected mean absolute error
	for a2 conditioned on a1, precomputed using the Haar distribution on USp(2g).
	See [KedlayaSutherland2007] for details.
*/

struct a2tab_entry g2a2tab[] = {
{ -4.00, 5.36, 0.09 },
{ -3.60, 4.72, 0.15 },
{ -3.20, 4.08, 0.17 },
{ -2.80, 3.44, 0.19 },
{ -2.40, 2.80, 0.22 },
{ -2.00, 2.24, 0.29 },
{ -1.60, 1.68, 0.38 },
{ -1.20, 1.20, 0.49 },
{ -0.80, 0.80, 0.62 },
{ -0.40, 0.48, 0.74 },
{ 0.00, 0.48, 0.74 },
{ 0.40, 0.80, 0.62 },
{ 0.80, 1.20, 0.49 },
{ 1.20, 1.68, 0.38 },
{ 1.60, 2.24, 0.29 },
{ 2.00, 2.80, 0.22 },
{ 2.40, 3.44, 0.19 },
{ 2.80, 4.08, 0.17 },
{ 3.20, 4.72, 0.15 },
{ 3.60, 5.36, 0.09 },
{ 4.00, 0.00 , 0.00 }		// can't happen, used to terminate
};

struct a2tab_entry g3a2tab[] = {
{ -6.00, 12.92, 0.15 },
{ -5.40, 10.84, 0.26 },
{ -4.80, 8.92, 0.34 },
{ -4.20, 7.16, 0.40 },
{ -3.60, 5.56, 0.45 },
{ -3.00, 4.12, 0.52 },
{ -2.40, 2.84, 0.61 },
{ -1.80, 1.72, 0.66 },
{ -1.20, 0.92, 0.61 },
{ -0.60, 0.60, 0.55 },
{ 0.00, 0.60, 0.55 },
{ 0.60, 0.92, 0.61 },
{ 1.20, 1.72, 0.66 },
{ 1.80, 2.84, 0.61 },
{ 2.40, 4.12, 0.52 },
{ 3.00, 5.56, 0.45 },
{ 3.60, 7.16, 0.40 },
{ 4.20, 8.92, 0.34 },
{ 4.80, 10.84, 0.26 },
{ 5.40, 12.92, 0.15 },
{ 6.00, 0.00 , 0.00 }		// can't happen, used to terminate
};

// various counters used for tuning and testing - none of these are functionally necessary
unsigned long jac_curve_count, jac_prebaby_gops, jac_baby_gops, jac_pregiant_gops, jac_giant_gops, jac_fastorder_gops, jac_ambexp_gops, jac_exp_gops, jac_charpoly_gops, jac_order_count;

/*
    The function jac_order computes #J(C/F_p) = P(1) given that Min <= #J(C/F_p) <= Max and (optionally in genus < 3)
    given the coefficient a1 of P(T) = L_p(T) (the numerator of Z(C/F_p;T) and (optionally) a list of constraints on the
    possible value of P(1) modulo 2(p+1). 

    The return value is 0 for an error, otherwise the number of multiples of the order of largest subgroup found 
    that lie in  the interval [Min,Max].  If there is only one multiple in the interval, *pP1 is set to the unique
    value which is then (provably) equal to P(1) = #J(C/F_p).  Otherwise *pP1 is set to the size of the largest,
    which will typically be less than Min.  In this case it is up to the caller to resolve the ambiguity.

    See [KedlayaSutherland2007] for more details.
*/

int jac_order (unsigned long *pP1, unsigned long Min, unsigned long Max, long a1, int *constraints, curve_t c[1], int fExponentOnly)
{
	static mpz_t Z;
	static int init;
	jac_t g, gen[JAC_MAX_GENERATORS];
	unsigned long ords[JAC_MAX_GENERATORS];
	unsigned long q[MPZ_MAX_UI_PP_FACTORS], h[MPZ_MAX_UI_PP_FACTORS];
	int i, j, k, two_rank, q_rank, cnt;
	unsigned long  m, n, e, e2, o, p, order, M, W, M2, W2;
	unsigned long SylowLimit, abandoned_q;
	double x, z;
	long t;
	int sts, full, ambiguous_exponent,tor3;
	
	if ( ! init ) { mpz_init (Z);  init = 1; }
	full = 0;															// default is to do a fastorder computation, not a full search of the interval
	p = (unsigned long)_ff_p;
	x = sqrt((double)p);
	if ( a1 != JAC_INVALID_A1 ) {
#if HECURVE_GENUS == 2
		t = (long)(p*p+1) + (long)(p+1)*a1 + (a1*a1)/2 - (long)2*p;			// lower bound on a2 given a1 is (a1^2)/2 - g*p
		if ( t > (long)Min ) Min = (unsigned long) t;
		t = (long)(p*p+1) + (long)(p+1)*a1 + (a1*a1)/4 + (long)2*p + 1;		// upper bound on a2 given a1 is (a1^2)/2 + g*p - (a1^2)/(2g)
		if ( t < (long)Max ) Max = (unsigned long)t;
		z = (double)a1/x;
		for ( i = 0 ; g2a2tab[i].a1 < z ; i++ );
		i--;
		M = (p*p+1) + (p+1)*a1 + g2a2tab[i].a2median*p ;					// Median value of a2 conditioned on a1 from table
		W = (unsigned long) ceil(g2a2tab[i].a2mae*p);						// Conditional expected distance of a2 from conditional median from table
#endif
#if HECURVE_GENUS == 3
		t = (long)(p*p*p +1) + (long)(p*p+1)*a1 + ((a1*a1)/2 - (long)3*p)*(p+1) - 20*p*(long)ceil(x);		// lower bound on a2 given a1 is (a1^2)/2 - g*p
		if ( t > (long)Min ) Min = (unsigned long) t;
		t = (long)(p*p*p+1) + (long)(p*p+1)*a1 + ((a1*a1)/3 + (long)3*p  + 1)*(p+1) + 20*p*(long)ceil(x);  // upper bound on a2 given a1 is (a1^2)/2 + g*p - (a1^2)/(2g)
		if ( t < (long)Max ) Max = (unsigned long)t;
		z = (double)a1/x;
		for ( i = 0 ; g3a2tab[i].a1 < z ; i++ );
		i--;
		M =(p*p*p +1) +(p*p+1)*a1 + g3a2tab[i].a2median*p*(p+1);			// Median value of a2 conditioned on a1 from table
		W = (unsigned long) ceil(g3a2tab[i].a2mae*p*p);					// Conditional expected distance of a2 from conditional median from table
#endif
	} else {
		M = Min + (Max-Min)/2;											// median value of a1- symmetrically distributed
#if HECURVE_GENUS == 1
		// In genus 1 we may want to do a "full" search of the interval to avoid a fastorder computation - see [KedlayaSutherland2007] section 4.2
		if ( p < SMALLJAC_FULL_P ) {
			W = (unsigned long) ceil(x);									//  (1+4/(3pi))sqrt(p) ~ 1.4244 but with divisor optimization, 1 is about right
			full = 1;
		} else {
			W = (unsigned long) ceil(0.8488*x);							// expected distance of a1 from the median is 8pi/3~0.8488...
		}
		if ( p > (1UL<<28) ) SMALLJAC_M = 8;								// tweak parallelism
#endif
#if HECURVE_GENUS == 2
		W = (unsigned long) ceil(0.7905*p*x);								// expected distance of a1 from the median is (exactly) 4096/525pi^2 = 0.790498...
#endif
#if HECURVE_GENUS == 3
		err_printf ("%lu: Invalid a1 value in jac_order in genus 3!\n", _ff_p);  exit (0);
#endif
	}
//printf ("%u: M = %ld, W = %ld, a1 = %ld, Min = %ld, Max = %ld\n", p, M, W, a1, Min, Max);
	ambiguous_exponent = 0;

	tor3 = 0;
#if HECURVE_GENUS == 1
	two_rank = (ff_poly_roots_d3(0,c[0].f)+1)/2;			// this takes < 1 microsecond (AMD Athlon 2.5GHz), equivalent to less than 10 gops
	if ( ! _ff_p1mod3 ) {								// for now only use 3tor when p=2mod3 and it is fast and easy
		tor3 = ff_poly_g1_3tor(c[0].f);
		if ( tor3 !=3 && tor3 != 9 ) tor3=1;			// don't deal with special cases
	}
#else
	two_rank = ff_poly_factors (c[0].f, c[0].d)-1;			// this takes 10 to 20 microseconds (AMD Athlon 2.5GHz)
	if ( a1 == JAC_INVALID_A1 ) 
		tor3 = ff_poly_g2_3tor(c[0].f);					// this takes 300-500 microseconds, don't use unless we are past pointcounting range
#endif
	e = ( two_rank >0 ? 2 : 1);
	if ( tor3 ) e*=tor3;
	for(;;) {
		M2 = M/e;
		W2 = _ui_ceil_ratio (W,e);
		for ( i = 0 ; i < SMALLJAC_RETRIES ; i++ ) {
			_jac_random (g, c[0]);
			if ( e > 1 ) jac_exp_ui (&g, &g, e, c);
			if ( ! _jac_is_identity (g) ) break;
		}
		if ( i == SMALLJAC_RETRIES ) { ambiguous_exponent = 1;  break; }
		o = jac_parallel_search (&g, M2, W2, Min/e, _ui_ceil_ratio(Max,e), full,  (two_rank?0:1), c);
		if ( ! o ) { err_printf ("%lu: search failed with e=%lu, M=%lu, W=%lu, m=%u for element: ", p, e, M, W, SMALLJAC_M);  _jac_print(g);  return 0; }
		e *= o;
		if ( fExponentOnly ) continue;
		e2 = e;
		if ( two_rank > 1&& ! full ) e2 <<= (two_rank-1);		// use knowledge of 2-rank but only when doing an order computation rather than a full search
		order = _ui_ceil_ratio(Min,e2) * e2;
		if ( ! two_rank && ! (order&0x1) ) order += e2;		// if 2-rank is 0, order must be odd
		if ( order > Max ) { err_printf ("%7u: No multiple of exponent %lu in interval (%lu,%lu) with 2-rank %d after computing order %lu with full = %d for element:   ", p, e, Min, Max, two_rank, o, full);  _jac_print(g);  return 0; }
		if ( order+e2 > Max ) break;
		if ( ! two_rank && !((order+e2)&0x1) && order + 2*e2 > Max ) break;
	}
	if ( ! two_rank && !(e&0x1) ) { err_printf ("%lu: Even exponent %ld found for group with trivial 2-rank\n", p, e);  return 0; }
	if ( fExponentOnly ) {
		*pP1 = e;
		return 1;
	}
	if ( ambiguous_exponent ) {
		n = jac_gops;
		e2 = e;
		if ( two_rank > 1 && ! full ) e2 <<= (two_rank-1);		// use knowledge of 2-rank, but only when not doing a full search
		dbg_printf ("%9lu: Ambiguous exponent %d, %d possible orders in [%ld,%ld]\n", p, e2, ui_multiples_in_range (e2, Min, Max), Min, Max);
		SylowLimit = ceil(sqrt(Max));
		k = ui_factor (q, h, e2);
		abandoned_q = 0;
		for ( i = 0 ; i < k ; i++ ) {
			if ( constraints ) {
				order = 0;
				cnt = 0;
				for ( o = _ui_ceil_ratio(Min,e2)*e2 ; o <= Max ; o += e2 ) {
					register int *pcon;
					for ( pcon = constraints ; *pcon >= 0 ; pcon++ ) {
						if ( (o % (2*(p+1))) == *pcon ) {
							cnt++;
							if ( cnt == 1 ) { order = o; }
						}
					}
				}
				if ( cnt == 1 ) { printf ("%lu: Ambiguity resolved via modular constraint.\n", p);  break; }
				if ( ! cnt ) { err_printf ("%lu: No orders satisfy modularity constraints!\n", p);  return 0; }
			}
			if ( q[i]*e2 <= Max ) {
				for ( j = 0, o = e2 ; j < h[i] ; j++, o/= q[i] );
				m = Max/o;
				// set M to the largest power of q[i] <= m, but use mpz to avoid overflow
				mpz_set_ui (Z, q[i]);
				while ( mpz_cmp_ui(Z,m) <= 0  ) { M = mpz_get_ui(Z);  mpz_mul_ui (Z, Z, q[i]); }
				for ( m = M ; m > SylowLimit ; m /= q[i] );
				q_rank = jac_sylow (gen, ords, q[i], e2, M, m, c);
				if ( q_rank < 0 ) {
					info_printf ("%lu: Abandoned large %d-Sylow subgroup computation, m= %lu,  SylowLimit = %lu.\n", p, q[i], m, SylowLimit);
					if ( abandoned_q ) { err_printf ("%lu: Abandoned second %d-Sylow subgroup computation - this should be impossible Max = %lu, m = %lu, SylowLimit = %lu\n", p, q[i], Max, m, SylowLimit);  return 0; }
					abandoned_q = q[i];
				}
				for ( j = 0 ; j < q_rank ; j++ ) o *= ords[j];
				if ( o > e2 ) {
					dbg_printf ("%7u: %d-Sylow computation enlarged subgroup to %lu\n", p, q[i], o);
					e2 = o;
				} else {
					dbg_printf ("%7u: %d-Sylow computation failed to enlarge subgroup\n", p, q[i]);
				}
				order = _ui_ceil_ratio(Min,e2)*e2;
				if ( ! two_rank && ! (order&0x1) ) order += e2;		// if 2-rank is 0, order must be odd
				if ( order > Max ) { err_printf ("%lu: No multiple of subgroup order %lu in interval (%lu,%lu) after %d-Sylow computation, 2-rank = %d\n", p, e2, Min, Max, q[i], two_rank); return 0; }
				if ( order + e2 > Max ) break;
				if ( ! two_rank && !((order+e2)&0x1) && order + 2*e2 > Max ) break;
				order = 0;
			}
		}
		jac_ambexp_gops += jac_gops-n;
		if ( i == k ) { 
			// in genus 3 don't rely on Sylow subgroup computations - we can't be sure we generated the whole group due to poor random element generation
#if HECURVE_GENUS == 3
			info_printf ("%lu: Maximal subgroup order %lu is ambiguous in (%lu,%lu)\n", p, e2, Min, Max); *pP1 = e2;  return ui_multiples_in_range (e2, Min, Max);
#else
			if ( ! abandoned_q ) { err_printf ("%lu: Scanned all Sylow subgroups and maximal order %lu is still ambiguous in (%lu,%lu)\n", p, e2, Min, Max);  *pP1 = e2;  return ui_multiples_in_range (e2, Min, Max); }
			for ( order = e2 ; order < Min ; order *= abandoned_q );
			if ( order > Max ) { err_printf ("%lu: No %lu-multitple of subgroup order %lu in interval (%lu,%lu) after all but one Sylow subgroups computed, 2-rank = %d\n",
									p, abandoned_q, e2, Min, Max, two_rank);  return 0; }
			info_printf ("%lu: Ambiguity resolved by computing all but 1 Sylow subgroup using %d group operations\n", p, jac_gops-n);
#endif
		}
		info_printf ("%lu: Ambiguity resolved using %d group operations\n", p, jac_gops-n);
	}
	*pP1 = order;
	return 1;
}


/*
	Parallel baby/giant search on [Min,Max] from starting point M with expected width W to find an exponent of the element a.
	The parity flag, if set, indicates that an odd exponent of a is known to lie in the interval.

	The returned value is either the order of a (the default), or if "full" is non-zero it may be either the order of a or the unique
	exponent of a in the interval [Min,Max].  Note that in either case, if [Min,Max] is known to contain a divisor of the
	group order, then the returned value divides the group order (and, in fact, the group exponent).

	The parameter full is usually 0 or 1, but can be set to 2 to force a true full search of the interval.  By default, the algorithm
	attempts to use divisors of the first exponent found to terminate the search early, and discontinues the search in the
	direction of the first divisor found - see section 4.2 of [KedlayaSutherland2007] for a discussion of these optimizations.

	6/1/2008: Many of the optimizations below are primarily relevant in genus 1, but these are now handled by the new code in hecurve1x.c.
			We may eventually want to take these out to simplify things, but for the moment we won't change anything.
*/
unsigned long jac_parallel_search (jac_t a[1],  unsigned long M, unsigned long W, unsigned long Min, unsigned long Max, int full, int parity, curve_t c[1])
{
	jac_t babys[SMALLJAC_M], giants[SMALLJAC_M], gsteps[SMALLJAC_M];
	jac_t b[64];		// needs to be at least as big as 2*SMALLJAC_M and SMALLJAC_VBITS
	jac_t s1, s2;
	register unsigned long GS, dvalue, uvalue;
	long e, n, o, o1;
	register unsigned i, j, k, cnt;
	unsigned baby_steps, baby_span;
	register uint32_t value, vstep;
	uint32_t matches[SMALLJAC_MAX_MATCHES];
	register int sign, stash, fast;
	unsigned long tabbits;
	int w;
	unsigned long p[MPZ_MAX_FACTORS], h[MPZ_MAX_FACTORS];

	
	n = jac_gops;
//printf ("Begin search p = %lu, M = %lu, W = %lu, Min = %lu, Max = %lu, parity = %d, full = %d\n", _ff_p, M, W, Min, Max, parity, full);  _jac_print(a[0]);
	if ( W == 0 ) { err_printf ("%9u: W == 0 in smalljac_parallel_search, increasing to 1\n", _ff_p);  W = 1; }
	if ( _jac_is_identity (a[0]) ) { return 1; }
	
	// We can check inverses for free, so we effectively double the giant step size and in s+s steps  we cover 2*s^2 = 2W, hence we want s=sqrt(W)
	// If we know the exponent is odd, we double the size of the baby steps and in s+s steps  we cover 4*s^2 = 2W, hence we want s=sqrt(W/2)
	if ( parity ) {	// odd exponent
		parity = 1;
		baby_steps = (unsigned)ceil(sqrt((double)W/2.0));
		if ( M&0x1 ) M--;			// M must be even
	} else {
		if ( full < 2 ) {	
			baby_steps = (unsigned)ceil(sqrt((double)W));
		} else {				
			baby_steps = (unsigned)ceil(sqrt((double)2.0*W));		// don't use inverse optimization - only for comparison
		}
	}
	baby_span = (parity ? 2*SMALLJAC_M*_ui_ceil_ratio(baby_steps,SMALLJAC_M) : SMALLJAC_M*_ui_ceil_ratio(baby_steps,SMALLJAC_M));
	if ( baby_span >= (1<<SMALLJAC_VBITS) ) { err_printf ("SMALLJAC_VBITS = %d is too small for baby span %d, increase SMALLJAC_VBITS\n", SMALLJAC_VBITS, baby_span);  exit (0); }
	stash = baby_span;
	if ( stash >= SMALLJAC_BABY_STASHSIZE ) stash = 0;

	_jac_set (babys[0], a[0]);
	for ( i = 1 ; i < SMALLJAC_M ; i += i ) jac_parallel_mult_1 (babys+i, babys, babys+i-1, c, i);
	if ( ! stash ) {									// if stash isn't big enough, precompute baby powers to use for reconstruction
		for ( i = 0 ; (1<<i) <= SMALLJAC_M ; i++ ) _jac_set (b[i], babys[(1<<i)-1]);
		for ( i-- ; i < SMALLJAC_VBITS ; i++ ) _jac_square (b[i+1], b[i], c[0]);
	}
	if ( parity ) {
		for ( i = 2 ; i < SMALLJAC_M ; i += 2 ) _jac_set (babys[(i+1)/2], babys[i]);
		_jac_set (s1, babys[SMALLJAC_M-1]);			// overlap is probably ok, but no need to push it
		jac_parallel_mult_1(babys+SMALLJAC_M/2, babys, &s1, c, SMALLJAC_M/2);
		_jac_mult (s1, babys[0], babys[SMALLJAC_M-1], c[0]);	// m-th baby has odd index 2m-1, we add one to get s1 = 2m: (1,3,...,2m-1) -> (2m+1,2m+3,...,4m-1)
		vstep = 2;
	} else {
		_jac_set (s1, babys[SMALLJAC_M-1]);			// s1 = m-th baby:  (1,2,...,m) -> (m+1,m+2,...,2m)
		vstep = 1;
	}

	// put identity in table
	_jac_set_identity (baby_stash[0]);
	for ( tabbits = 8 ; (1<<tabbits) < 2*baby_steps ; tabbits++ );
	smalljac_table_init (tabbits);
	smalljac_table_insert (baby_stash, 0);
	// take baby steps
	o = 0;
	value = 1;
	n = jac_gops;
	if ( stash ) {
		for (;; ) {
			for ( i = 0 ; i < SMALLJAC_M ; i++ ) {
				if ( _jac_is_identity (babys[i]) ) { o = value;  goto found; }
//printf ("%lu: inserting baby value %d: ", _ff_p, value);  _jac_print(babys[i]);
				smalljac_table_insert (babys+i, value);
				_jac_set (baby_stash[value], babys[i]);					// we could avoid copying here, but parity makes it awkward
				value += vstep;
			}
			if ( value > baby_span ) break;									// value is now the index of the next baby we will compute - one (or two) greater than the last
			jac_parallel_mult_1 (babys, babys, &s1, c, SMALLJAC_M);
		}
	} else {
		for (;; ) {
			for ( i = 0 ; i < SMALLJAC_M ; i++ ) {
				if ( _jac_is_identity (babys[i]) ) { o = value;  goto found; }
				smalljac_table_insert (babys+i, value);
				value += vstep;
			}
			if ( value > baby_span ) break;									// value is now the index of the next baby we will compute - one (or two) greater than the last
			jac_parallel_mult_1 (babys, babys, &s1, c, SMALLJAC_M);
		}
	}
	value -= vstep;													// set value to index of last baby computed
	jac_baby_gops += jac_gops-n;
//printf ("babys done, used %lu gops\n", jac_gops-n);	

	e = M/value;														// the cost of the exponentiation can be quite significant in genus 1 - > 1/4 of the total time for p~2^24
	if ( parity && (e&1) ) e++;											// we use the largest baby to get the first giant by tweaking M to be a multiple of value (with the right parity)
	M = e * value;														// this adds at most one giant step.
	n = jac_gops;
	jac_exp_ui (giants, babys+SMALLJAC_M-1, e, c);							// exponentiate biggest baby to get first giant
	jac_exp_gops += jac_gops-n;
	
	_jac_square (s1, babys[SMALLJAC_M-1], c[0]);								// s1 is double last baby
	if ( parity ) {
		if ( full < 2 ) GS = 2*value;
		else GS = ( (value&1) ? value+1 : value );							// giant spacing must be even if babys are odd
	} else {
		_jac_mult (s1, s1, a[0], c[0]);										// we can add one otherwise
		if ( full < 2 ) GS = 2*value+1; else GS = value;
	}
	// compute m/2 giants - currently done serially
	for ( i = 1 ; i < SMALLJAC_M/2 ; i++ ) _jac_mult (giants[i], giants[i-1], s1, c[0]);
	for ( i = 1 ; i < SMALLJAC_M/2 ; i += i ) _jac_square (s1, s1, c[0]);				// power up s1 so that each giant will step m/2 giant spacings
	_jac_set (s2, s1);
	if ( _jac_is_identity (s2) ) { o = GS*(SMALLJAC_M/2);  goto found; }			// unlikely but possible and the giants won't get very far if this happens
	_jac_invert (s2);													// s2 = -s1 is used for downward steps
	jac_parallel_mult_1 (giants+SMALLJAC_M/2, giants, &s2, c, SMALLJAC_M/2);		// have m/2 up giants step down once to form m/2 down giants - this means down giants are in reverse order but we can cope
	for ( i = 0 ; i < SMALLJAC_M/2 ; i++ ) {									// gsteps array is just m/2 copies of s1 and s2, used to perform up and down steps in parallel
		_jac_set (gsteps[i],s1);
		_jac_set (gsteps[i+SMALLJAC_M/2],s2);
	}

	uvalue = M;														// index of first up giant
	dvalue = M-GS;													// index of "last" down giant giants[SMALLJAC_M-1] - this is the one with the largest index

//printf ("%lu: Beginning giant steps GS = %ld, uvalue = %ld, dvalue = %ld, used %lu gops\n", _ff_p, GS, uvalue, dvalue, jac_gops-n);	

	//giant steps 
	// the code duplication below is intentional.  it could be avoided, but at the cost of more conditional code inside the loops
	o = o1 = 0;
	for ( j = 0 ; ; j++ ) {												// search until we find it - could add Min/Max check
		if ( uvalue ) {
			for ( i = 0 ; i < SMALLJAC_M/2 ; i++ ) {						// up giants
//printf ("%lu: checking up giant = %lu ", _ff_p, uvalue);  _jac_print(giants[i]);
				if ( (cnt = smalljac_table_lookup (matches, giants+i)) > 0 ) {
//printf ("%9lu: matched %d up\n", _ff_p, cnt);
					for ( k = 0 ; k < cnt ; k++ ) {
//printf ("%d\n", matches[k]);
						if ( stash ) {
							sign = _jac_cmp (baby_stash[matches[k]], giants[i]);
						} else {
							jac_exp_ui32_powers (&s1, b, matches[k], c); 			// reconstruct baby for comparison
							sign = _jac_cmp (s1, giants[i]);
						}
						if ( sign ) {
//printf ("%9lu: Found a match at uvalue = %lu, k = %d, matches[k] = %u, sign = %d, %d gops\n", _ff_p, uvalue, k, matches[k], sign, jac_gops-n);
							if ( o ) o1 = o;
							o = ( sign > 0 ? uvalue - matches[k] : uvalue + matches[k] );
							if ( ! o ) { o = o1;  continue; }						// can happen for small groups
							if ( o < 0 ) o = -o;									// can happen for small groups
							if ( ! full ) goto found;
							if ( o == o1 ) continue;								// yes this can happen - up and down steps can overlap in odd parity case
							if ( o1  ) { o = ( o > o1 ? o-o1 : o1-o);  goto found; }
							if ( full < 2 && j && 2*o > Max-Min ) {					// make sure we are actually in the top half of the interval - we might not have started at the center
								if ( ! dvalue ) return o;							// if we already hit the other end, we're done
								w = ui_factor (p, h, o);
								if ( w == 1 && h[0] == 1 ) { jac_order_count++; jac_giant_gops += jac_gops-n; return o; }
								if ( p[w-1] > o-Min && (vstep*o/p[w-1])+dvalue+baby_span < o ) {
//printf ("%lu: [%lu,%lu], exp %lu, q=%lu, o-dvalue-baby_span=%lu, dvalue=%lu, baby_span = %d, parity =%d\n", _ff_p, Min, Max, o, p[w-1], o-dvalue-baby_span, dvalue, baby_span, parity); 
									jac_order_count++; jac_giant_gops += jac_gops-n; return o;
								}
								uvalue = 0;
								goto udone;
							}
						}
					}
				}
				if ( uvalue ) uvalue += GS;
			}
		}
udone:
		if ( dvalue ) {
			for ( i = SMALLJAC_M-1 ; dvalue && i >= SMALLJAC_M/2 ; i-- ) {			// down giants - process largest dvalue first
//printf ("%lu: checking down giant = %lu ", _ff_p, dvalue);  _jac_print(giants[i]);
				if ( (cnt = smalljac_table_lookup (matches, giants+i)) > 0 ) {
//printf ("%9lu: matched %d down\n", _ff_p, cnt);
					for ( k = 0 ; k < cnt ; k++ ) {
						if ( stash ) {
							sign = _jac_cmp (baby_stash[matches[k]], giants[i]);
						} else {
							jac_exp_ui32_powers (&s1, b, matches[k], c); 			// reconstruct baby for comparison
							sign = _jac_cmp (s1, giants[i]);
						}
						if ( sign ) {
//printf ("%9lu: Found a match at dvalue = %lu, k = %d, matches[k] = %u, sign = %d, %d gops\n", _ff_p, dvalue, k, matches[k], sign, jac_gops-n);
							if ( o ) o1 = o;
							o = ( sign > 0 ? dvalue - matches[k] : dvalue + matches[k] );
							if ( ! o ) { o = o1;  continue;	}						// can happen for small groups
							if ( o < 0 ) o = -o;									// can happen for small groups
							if ( ! full ) goto found;
							if ( o == o1 ) continue;								// yes this can happen - up and down steps can overlap in odd parity case
							if ( o1 ) { o = ( o > o1 ? o-o1 : o1-o);  goto found; }
							if ( full < 2 && 2*o < Max-Min ) {
								if ( ! uvalue ) return o;							// if we already hit the other end, we're done
								w = ui_factor (p, h, o);
								if ( w == 1 && h[0] == 1 ) { jac_order_count++; jac_giant_gops += jac_gops-n; return o; }
								if ( p[w-1] > Max-o && (vstep*o/p[w-1])+o+baby_span < uvalue ) {
//printf ("%lu: [%lu,%lu], exp %lu, q=%lu, uvalue-o-baby_span=%lu, parity =%d\n", _ff_p, Min, Max, o, p[w-1], uvalue-o-baby_span, parity); 		
									jac_order_count++; jac_giant_gops += jac_gops-n;  return o;
								}
								dvalue = 0;
								goto ddone;
							}
						}
					}
				}
				if ( dvalue <= GS ) dvalue = 0; else dvalue -= GS;
			}
		}
ddone:
		if ( uvalue > Max+baby_span ) uvalue = 0;
		if ( dvalue+baby_span < Min ) dvalue = 0;
		if ( ! dvalue && ! uvalue ) break;
		if ( ! dvalue ) {
			jac_parallel_mult_c (giants, giants, gsteps, c, SMALLJAC_M/2);
		} else if ( ! uvalue ) {
			jac_parallel_mult_c (giants+SMALLJAC_M/2, giants+SMALLJAC_M/2, gsteps+SMALLJAC_M/2, c, SMALLJAC_M/2);
		} else {
			jac_parallel_mult_c (giants, giants, gsteps, c, SMALLJAC_M);
		}
	}
	jac_giant_gops += jac_gops-n;
//printf ("giant steps in full computation, used %lu gops\n", jac_gops-n);	
	if ( ! o ) {
		err_printf ("%lu: No exponent found in smalljac_parallel_search, full=%d, parity=%d, baby_span = %d, GS = %ld, dvalue = %d, uvalue = %ld, Min = %lu, Max = %lu!\n",
				 _ff_p, full, parity, baby_span, GS, dvalue, uvalue, Min, Max);  _curve_print(c[0]);
		_jac_print(a[0]);
		return 0;
	}
//if ( full ) jac_full_not_found_count++;
//printf ("%lu: found one exponent %lu\n", _ff_p, o);
	if ( o < 0 ) { err_printf("%lu: Negative exponent %ld found with full=%d, o1=%ld for element ", _ff_p, o, full, o1);  _jac_print(a[0]);  return 0; }
	return o;
found:
	if ( ! o ) { err_printf("%lu: Zero exponent found with full=%d, o1=%ld for element ", _ff_p, full, o1);  _jac_print(a[0]);  exit (0); }
	if ( o < 0 ) { err_printf("%lu: Negative exponent %ld found with full=%d, o1=%ld for element ", _ff_p, o, full, o1);  _jac_print(a[0]);  return 0; }

//printf ("found, used %lu gops\n", jac_gops-n);	
//if ( j > 3*baby_s ) printf ("%9u: Long search - %d baby steps, %d giant steps, W = %ld\n", _ff_p, baby_s*SMALLJAC_M, j*SMALLJAC_M, W); 
//printf ("giant steps in order computation, used %lu gops\n", jac_gops-n);	
	jac_giant_gops += jac_gops-n;
	n = jac_gops;
	w = ui_factor (p, h, o);
/*
	There is no reason not to do this, but we need to be sure the caller handles it - currently it is assumed that the order is returned unless the full flag is set.
	
	if ( j && p[w-1] > o-Min && (vstep*o/p[w-1])+dvalue+baby_span < o ) { jac_order_count++;  return o; }
	if ( p[w-1] > Max-o && (vstep*o/p[w-1])+o+baby_span < uvalue ) { jac_order_count++;  return o; }
*/
	jac_factored_order_ui (&o, a, p, h, w, c);
	jac_fastorder_gops += jac_gops-n;
//jac_exp_ui (&s1, a, o, c);
//if ( ! _jac_is_identity(s1) ) { err_printf ("%d: %d is not an exponent of ", _ff_p, o);  _jac_print(s1);  exit (0); }
//printf ("%u: %lu is the order of ", _ff_p, o); _jac_print(a[0]);
	if ( o < 0 ) { err_printf("%lu: Negative exponent %ld after fastorder with full=%d, o1=%ld for element ", _ff_p, o, full, o1);  _jac_print(a[0]);  return 0; }
	return o;
}


// Straight baby/giant search over values Min + km up to Max as k ranges from 0 to ceil((Max-Min)/m)
// Finds least two matches in interval, returns number of matches found
int jac_search (mpz_t e[2], jac_t a[1],  unsigned m, mpz_t Min, mpz_t Max, curve_t c[1])
{
	static int init;
	static mpz_t o, o1, o2, G;
	jac_t b[32], baby, giant;
	jac_t s1, step, babystep;
	unsigned long S, tabbits;
	unsigned i, j, k, s, cnt;
	uint32_t matches[SMALLJAC_MAX_MATCHES];
	int sign;

	if ( ! init ) { mpz_init (o);  mpz_init(o1);  mpz_init (o2);  mpz_init(G); init = 1; }
	mpz_sub (G, Max, Min);
	S = mpz_get_ui (G);
	s = (unsigned) ceil(sqrt(_ui_ceil_ratio(S, m)));
//out_printf ("Small search, search p = %d, m = %d, Min = %Zd, Max = %Zd, s = %d\n", _ff_p, m, Min, Max, s);  _jac_print(a[0]);	
	for ( tabbits = 8 ; (1<<tabbits) < 2*s ; tabbits++ );
	smalljac_table_init (tabbits);
	// put identity in table
	_jac_set_identity (step);
	smalljac_table_insert (&step, 0);
	// Compute babystep = first baby
	jac_exp_ui (&babystep, a, m, c);
	_jac_set (baby, babystep);
	// Compute baby powers used to recompute baby matches
	_jac_set (b[0], babystep);
	for ( i = 1 ; i < 32 ; i++ ) _jac_square(b[i], b[i-1], c[0]);
	// Compute first giant
//	G = Min + (unsigned long)s*(unsigned long)m;							// First giant starts one babyspan into the interval
	mpz_set(G, Min);
	mpz_add_ui (G, G, (unsigned long)s*(unsigned long)m);
	jac_exp_mpz (&giant, a, G, c);
	// Compute giant step - twice the span of the babys
	S = 2*(unsigned long)s*(unsigned long)m;							// giant step size - should fit in 64 bits
	jac_exp_ui (&step, a, S, c);											// could be clever and use last baby, but only saves a few gops - not needed here

	// baby steps
	for ( i = 1 ; ; i++ ) {
		smalljac_table_insert (&baby, i);
		if ( i > s ) break;
		_jac_mult (baby, baby, babystep, c[0]);
	}
	
	// trim s to avoid taking extra giant steps - it may be larger than necessary
	for (;;) {
		mpz_set (o, G);
		mpz_add_ui (o, o, (unsigned long)s*(S-(unsigned long)m));
		if ( mpz_cmp (o, Max) <= 0 ) break;
		s--;
	}
	
	//giant steps 
	mpz_set_ui (o, 0);
	mpz_set_ui (o1, 0);  mpz_set_ui (o2, 0);
	for ( j = 0 ; j < s ; j++ ) {
		if ( (cnt = smalljac_table_lookup (matches, &giant)) > 0 ) {
			for ( k = 0 ; k < cnt ; k++ ) {
				if ( matches[k] ) {
					// reconstruct baby for comparison
					jac_exp_ui32_powers (&s1, b, matches[k], c);
					if ( (sign = _jac_cmp(s1,giant)) ) {
						mpz_set(o, G);
						mpz_add_ui (o, o, j*S);
						if ( sign > 0 ) { mpz_sub_ui (o, o, (unsigned long)matches[k]*(unsigned long)m); } else { mpz_add_ui (o, o, (unsigned long)matches[k]*(unsigned long)m); }
					} else {
						mpz_set_ui (o, 0);
					}
				} else {
					mpz_set (o, G);
					mpz_add_ui (o, o, j*S);
				}
//jac_exp_mpz (&s1, b, o, c);
//if ( ! _jac_is_identity (s1) ) { err_printf ("%Zd is not an exponent of ", o); _jac_print(b[0]);  exit (0); }
				// track the smallest two matches - note that we don't necessarily find them in order.
				if ( mpz_sgn(o) && mpz_cmp(o,Min) >= 0 && mpz_cmp (o,Max) <= 0 ) {		// ignore matches outside the interval - these can happen and can cause confusion when they do
					if ( ! mpz_sgn(o1) ) {
						mpz_set (o1, o);;
					} else {
						if ( mpz_cmp(o,o1) < 0 ) {
							mpz_set (o2, o1);  mpz_set (o1, o);
						} else if ( mpz_cmp(o,o1) > 0 && (! mpz_sgn(o2) || mpz_cmp(o,o2) < 0) ) {
							mpz_set (o2, o);
						}
					}
				}
			}
		}
		if ( j == s ) break;
		_jac_mult (giant, giant, step, c[0]);
	}
//out_printf ("o1=%Zd, o2=%Zd\n", o1, o2);
	mpz_set (e[0], o1);
	mpz_set (e[1], o2);
	return ( mpz_sgn(o2) ? 2 : ( mpz_sgn(o1) ? 1 : 0 ));
}
