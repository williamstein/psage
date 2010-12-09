#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gmp.h"
#include "mpzutil.h"
#include "ffwrapper.h"
#include "ffpoly.h"
#include "hecurve.h"
#include "jac.h"
#include "jacorder.h"
#include "smalljac.h"
#include "smalljac_g23.h"
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
    This module contains genus 2 and genus 3 specific code for deriving the
    numerator of the zeta function, P(T) = L_p(T), when given certain
    partial information such as the coefficients mod p or the values P(1) and P(-1).
*/

/*
	Derive P(T) from P(1) and P(-1) in genus 2.  No group operations required.  Assumes p < 2^31
	Sanity check coefficients - caller may not be sure about P(1) and P(-1) and may use this to narrow down the possibilities
*/
int smalljac_genus2_charpoly (long a[3], long p, double sqrtp, long P1, long PN1)
{
	if ( p > (1L<<31) ) { err_printf ("Overflow in smalljac_genus3_charpoly, p = %ld\n", p);  return 0; }
	a[0] = (P1-PN1)/(2*(p+1));
	a[1] = ((P1+PN1) - 2*(p*p+1)) / 2;
	if ( (p*p+1) + (p+1)*a[0] + a[1] !=  P1 ) return 0;
	if ( (p*p+1) - (p+1)*a[0] + a[1] != PN1 ) return 0;
	if ( abs(a[0]) > 4*sqrtp ) return 0;										// standard Weil bound
	if ( a[1] > 2*p + a[0]*a[0]/4 || a[1] < -2*p+a[0]*a[0]/2 ) return 0;			// these bounds follow from Prop. 4 of [KedlayaSutherland2007]
	if ( abs(a[2]) > 20*p*sqrtp ) return 0;
	return 1;
}


/*
	Derive P(T) from P(1), P(-1) and #C/F_p in genus 3.  No group operations required.  Assumes p<2^20
	Sanity check coefficients - caller may not be sure about P(1) and P(-1) and may use this to narrow down the possibilities
*/
int smalljac_genus3_charpoly (long a[3], long p, double sqrtp, long P1, long PN1, unsigned long pts)
{
	if ( p > (1L<<20) ) { err_printf ("Overflow in smalljac_genus3_charpoly, p = %ld\n", p);  return 0; }
	a[0] = (long)pts - (long)(p+1);
	a[2] = ((P1-PN1) - 2*(p*p+1)*a[0]) / 2;
	a[1] = (P1 - (p*p*p+1) - (p*p+1)*a[0] - a[2])/(p+1);
	if ( (p*p*p+1) + (p*p+1)*a[0] + (p+1)*a[1] + a[2] !=  P1 ) return 0;
	if ( (p*p*p+1) - (p*p+1)*a[0] + (p+1)*a[1] - a[2] != PN1 ) return 0;
	if ( abs(a[0]) > 6*sqrtp ) return 0;												// standard Weil bound
	if ( a[1] > 3*p + a[0]*a[0]/3 || a[1] < -3*p +a[0]*a[0]/2 ) return 0;					// these bounds on a2 follow from Prop. 4 of [KedlayaSutherland2007]
	if ( abs(a[2]) > 20*p*sqrtp ) return 0;
	return 1;
}


/*
	Derive P(T) from P(1) in genus 2, using O(\lg(p)) gops in the twist.
*/

int smalljac_genus2_charpoly_from_P1 (long a[2], long P1, long Min, long Max, curve_t c[1])
{
	curve_t twist;
	jac_t h, g;
	unsigned long n;
	long e, o, to, to0, a0, p;
	int i,k;
	
	n = jac_gops;
	ff_poly_twist (twist.f, c[0].f, c[0].d);
	twist.d = 5;

	// We want to compute a and b such that P(T) = 1 + a1T + a2T^2 + pa1T^3 + p^2T^4
	// We are given P(1) and we know that P(-1) is the order of the twist
	// We also know that a1 >= (P(1) - (p^2+1) - 6p)/(p+1) and a <= (P(1) - (p^2+1) + 6p)/(p+1)
	// which gives a range of at most 12 values for a1, and thus at most 12 possibilities for P(-1) = P(1) - 2(p+1)a1
	// all of which differ by a multiple of 2(p+1).  We just need to distinguish P(-1) among these 11 possibilities.
	p = _ff_p;
	a0 = (P1 - p*p + 6*p - 1)/(p+1);						// max value for a1 = (P(1) - (p^2+6p+1)) / (p+1) gives min value for twist order to
	to = P1 - 2*(p+1)*a0;
	to0 = to;
	i = 0;
	while ( to < Min ) { i++; to += 2*(p+1); }				// Skip twist orders not in Weil interval
	for ( k = 0 ; k < SMALLJAC_RETRIES ; k++ ) {
		_jac_random(g, twist);
		jac_exp_ui (&h, &g, 2*(p+1), &twist);
		if ( ! _jac_is_identity(h) ) break;
	}
	if ( k == SMALLJAC_RETRIES ) {
		// 2*(p+1) is an exponent for the group - hopefully only one candidate twist order is compatible, since they are all exponents
		o = 0;
		for ( ; i <= 12 ; i++, to += 2*(p+1) ) {
			if ( to > Max ) break;
			if ( ui_compatible (to, 2*(p+1)) ) {
				if ( o ) { /*err_printf ("%7u: More than one twist order implied by P(1) = %ld is a 2(p+1) compatible exponent of the twist (%ld,%ld).\n", _ff_p, P1, o, to);*/  return 0; }
				o = to;
			}
		}
	} else {
		// Ok we have an element whose order contains a prime not in 2*(p+1).  This should almost always uniquely determine the twist
		// since two candidate twist orders differ by i*2*(p+1) where i <= 10
		o = 0;
		jac_exp_ui (&g, &g, to, &twist);
		for ( ; i <= 12 ; i++, to += 2*(p+1) ) {
			if ( to > Max ) break;
			if ( _jac_is_identity (g) ) {
				if ( o ) {	// we know that e=(to-o) is an exponent of g - try to find an element for which e is not an exponent
					e = to-o;
					for ( k = 0 ; k < SMALLJAC_RETRIES ; k++ ) {
						_jac_random(g, twist);
						jac_exp_ui (&h, &g, e, &twist);
						if ( ! _jac_is_identity(h) ) break;
					}
					if ( k == SMALLJAC_RETRIES ) {/*err_printf ("%7u: More than one twist order implied by P(1) = %ld is an exponent of the twist (%ld,%ld).\n", _ff_p, P1, o, to);*/  return 0; }
					jac_exp_ui (&h, &g, o, &twist);
					if ( _jac_is_identity(h) ) {
						// to can't be an exponent, and o is probably the twist order, reset g and h to continue the search
						jac_exp_ui (&h, &g, 2*(p+1), &twist);
						jac_exp_ui (&g, &g, to, &twist);
					} else {
						o = 0;
						jac_exp_ui (&h, &g, 2*(p+1), &twist);
						jac_exp_ui (&g, &g, to, &twist);
						if ( _jac_is_identity (g) ) o = to;
					}
				} else {
					o = to;
				}
			}
			_jac_mult (g, g, h, twist);
		}
	}
	if ( ! o ) { err_printf ("%7u: None of the implied twist orders for P(1) = %ld are exponents with to0 = %ld in Weil interval [%ld,%ld]\n", _ff_p, P1, to0, Min, Max);  return 0; }
	a[0] = (P1 - o) / (2*(p+1));		// a1 = P(1)-P(-1)/(2(p+1))
	a[1] = (P1+o-2*(p*p+1))/2;		// a2 = (P(1)+P(-1) - 2(p^2+1))/2
	smalljac_charpoly_gops += jac_gops-n;
	return 1;
}


/*
	Derive P(T) from P(1) in genus 3, using O(p^{1/4}) gops in the twist.
*/

int smalljac_genus3_charpoly_from_P1 (long o[3], long P1, long pts, long Min, long Max, curve_t c[1])
{
	static int init;
	static mpz_t Tmin, Tmax, e[2];
	curve_t twist;
	jac_t g;
	long p, tmin, tmax, to, spc, bmin, bmax, k;
	unsigned long n, e0;
	double x;
	int constraints[2];
	int i, cnt;
	
	n = jac_gops;
	if ( ! init ) { mpz_init (Tmin);  mpz_init (Tmax);  mpz_init (e[0]);  mpz_init (e[1]);  init = 1; }
	
	ff_poly_twist (twist.f, c[0].f, c[0].d);
	twist.d = c[0].d;

	p = (long) _ff_p;
	o[0] = pts - 1 - p;
	x = sqrt((double)p);
	spc = (long)ceil(x);
	bmin = (P1 - (p*p*p+1) - (p*p+1)*o[0] - 20*p*spc) / (p+1);
	tmin = 2*(p*p*p+1) - P1+ 2*(p+1)*bmin;
	bmax = (P1 - (p*p*p+1) - (p*p+1)*o[0] + 20*p*spc) / (p+1) + 1;	
	tmax = 2*(p*p*p+1) - P1+ 2*(p+1)*bmax;
	if ( tmin < Min ) tmin += _ui_ceil_ratio((Min-tmin),2*(p+1)) * 2*(p+1);
	tmax = 2*(p*p*p+1) - P1+ 2*(p+1)*15*p;
	if ( tmax > Max ) tmax -= (tmax-Max)/(2*(p+1)) * 2*(p+1);
	_jac_random (g,twist);
	mpz_set_ui (Tmin, tmin);  mpz_set_ui (Tmax, tmax);
	cnt = jac_search (e, &g, 2*(p+1), Tmin, Tmax, &twist);
	if ( ! cnt ) { err_printf ("%7u: Search failed in genus3_charpoly with P(1) = %ld in [%ld,%ld]\n", _ff_p, P1, Min, Max);  return 0; }
	if ( cnt > 1 ) {
		k =  2*(p*p*p+1) - P1;		// note k may be big
		constraints[0] = ( k < 0 ? (int) (2*(p+1) - ((-k)%(2*(p+1)))) : (int) (k%(2*(p+1))) );
		constraints[1] = -1;
		if ( jac_order (&to, Min, Max, -o[0], constraints, &twist, 0) != 1 ) {
			err_printf ("%7u: Ambiguous result in genus3_charpoly with P(1) = %ld in [%ld,%ld] (%ld,%ld), order computation for twist failed\n", _ff_p, P1, Min, Max, e[0], e[1]);
			return 0;
		}
	} else {
		e0 = mpz_get_ui (e[0]);
		to = _ui_ceil_ratio(tmin,e0)*e0;
	}
	
	if ( ! smalljac_genus3_charpoly (o, p, x, P1, to, pts) ) { err_printf ("%7u: Unable to derive consistent Lpolynomial from P(1) = %lu and P(-1) = %lu given a1 = %ld\n", p, P1, to, o[0]);  return 0; }
	smalljac_charpoly_gops += jac_gops-n;
	return 1;
}

/*
	Derive P(T) given P(T) mod p in genus 2, using O(lg(p)) group operations.
*/
int smalljac_genus2_charpoly_from_Pmodp (long a[2], curve_t c[1])
{
	jac_t g, h;
	unsigned long n;
	long p, e0, e1, N, Min, Max;
	
	p = (long) _ff_p;

	n = jac_gops;
	if ( a[1] < 0 ) a[1] += p;
	Min = (p*p+1) + (p+1)*a[0] - 2*p;
	Max = (p*p+1) + (p+1)*a[0] + 6*p;
	N = ((p*p+1)+(p+1)*a[0] + a[1]) % p;						// N = P(1) mod p
	Min = p*(Min/p) + N;									// congruent to P(1) mod p
	e0 = Min;
	_jac_random(g, c[0]);
	jac_exp_ui (&h, &g, e0, c);
	jac_exp_ui (&g, &g, p, c);
	e1 = 0;
	while ( e0 <= Max ) {
		if ( _jac_is_identity (h) ) {
			if ( e1 ) { err_printf ("%lu: Ambiguous result in genus2_charpoly_from_Pmodp %lu and %lu possible exponents\n", p, e1, e0);  return 0; }
			e1 = e0;
		}
		_jac_mult (h, h, g, c[0]);
		e0 += p;
	}
	if ( ! e1 ) { err_printf ("%9u: Search failed in genus2_charpoly_from_Pmodp, a[1] = %ld, Min = %ld, Max = %ld\n", _ff_p, a[1], Min, Max);  return 0; }
	a[1] = e1 - (p*p+1) - (p+1)*a[0];
	smalljac_charpoly_gops += jac_gops-n;
	return 1;
}

/*
	Derive P(T) given P(T) mod p in genus 3, using O(p^{1/4}) group operations.

	We assume p is big enough to avoid any unpleasant special cases but less than 2^31
*/
int smalljac_genus3_charpoly_from_Pmodp (long a[3], curve_t c[1])
{
	static int init;
	static mpz_t Min, Max, e[2], z, O;
	jac_t g;
	long p, M, N, T, a2, a2Min, a2Max, good_a2, spc;
	unsigned long n;
	double x;
	int i, cnt;

	if ( ! init ) { mpz_init (Min);  mpz_init (Max);  mpz_init (e[0]);  mpz_init (e[1]); mpz_init (z); mpz_init (O);  init = 1; }

	p = (long) _ff_p;
	x = sqrt((double)p);
	spc = (long)ceil(x);

	cnt = 0;
	n = jac_gops;
	if ( a[1] < 0 ) a[1] += p;
	// Set a2 bounds based on Proposition 4 of KedlayaSutherland 2007
	a2Min = -p+a[1];
	a2Max = 3*p + _ui_ceil_ratio(a[0]*a[0], 3);
	for ( a2 = a2Min ; a2 <= a2Max ; a2 += p ) {
		mpz_set_ui (Min, p*p);  mpz_mul_ui (Min, Min, p);  mpz_add_ui (Min, Min, 1);
		mpz_set_i (z, a[0]);  mpz_mul_ui (z, z, p*p+1);  mpz_add (Min, Min, z);
		mpz_set_i (z, a2);  mpz_mul_ui (z, z, p+1); mpz_add (Min, Min, z);
		mpz_set_i (e[0], a[2]);
		mpz_add (z, Min, e[0]);
		mpz_mod_ui (z, z, p);
		N = mpz_get_ui (z);								// N is P(1)=p^3+1 + (p^2+1)a1 + (p+1)a2 +a3 mod p
		mpz_add_ui (Max, Min, 20*p*spc);
		mpz_sub_ui (Min, Min, 20*p*spc);
		mpz_tdiv_q_ui (Min, Min, p);
		mpz_mul_ui (Min, Min, p);
		mpz_add_ui (Min, Min, N);
		_jac_random (g, c[0]);
		i = jac_search (e, &g, p, Min, Max, c);
		if ( ! cnt && i ) { mpz_set(O,e[0]);  good_a2 = a2; }
		cnt += i;
		if ( cnt > 1 ) {
			err_printf ("%lu: Ambiguous result in genus3_charpoly_from_Pmodp %Zd and %Zd possible exponents (i=%d)\n", _ff_p, O, (i>1?e[1]:e[0]),i);
			return 0;
		}
	}
	if ( cnt != 1 ) {
		err_printf ("%9u: Search failed in genus3_charpoly_from_Pmodp a2Min = %ld, a2Max = %ld, Min = %Zd, Max = %Zd\n", _ff_p, a2Min, a2Max, Min, Max);
		err_printf ("%9u: padic_lpoly coefficients: a1=%ld, a2=%ld, a3=%ld\n", _ff_p, a[0], a[1], a[2]);
		return 0;
	}
	a[1] = good_a2;
	mpz_set_ui (z, p*p);  mpz_mul_ui (z, z, p);  mpz_add_ui (z, z, 1);
	mpz_sub (O, O, z);
	mpz_set_i (z, a[0]);  mpz_mul_ui (z, z, p*p+1);
	mpz_sub (O, O, z);
	mpz_set_i (z, a[1]);  mpz_mul_ui (z, z, p+1);
	mpz_sub (O, O, z);
	a[2] = mpz_get_i (O);
	smalljac_charpoly_gops += jac_gops-n;
	return 1;
}
