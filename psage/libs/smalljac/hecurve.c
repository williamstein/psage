#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "gmp.h"
#include "mpzutil.h"
#include "ffwrapper.h"
#include "ffext.h"
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

/*
	General module for hyperelliptic curve group operations.
	Fast genus specific group operations for the common cases
	are in hecurve1.c, hecurve2.c, and hecurve3.c.  This module
	handles everything else.  It also includes a non-genus specific
	implementation of Cantor's algorithm.
*/

static char buf[4096];

#if HECURVE_GENUS == 1
#define _deg_u(u)			(_ff_nonzero(u[1])?0:1)
#define _deg_v(v)			(_ff_nonzero(v[0])?0:-1)
#define _hecurve_set_zero(u,v)		{_ff_set_zero(u[0]);_ff_set_zero(v[0]);}
#endif
#if HECURVE_GENUS == 2
#define _deg_u(u)			(_ff_nonzero(u[2])?2:(_ff_nonzero(u[1])?1:0))
#define _deg_v(v)			(_ff_nonzero(v[1])?1:(_ff_nonzero(v[0])?0:-1))
#define _hecurve_set_zero(u,v)		{_ff_set_zero(u[2]);_ff_set_zero(u[1]);_ff_set_zero(u[0]);_ff_set_zero(v[1]);_ff_set_zero(v[0]);}
#endif
#if HECURVE_GENUS == 3
#define _deg_u(u)			(_ff_nonzero(u[3])?3:(_ff_nonzero(u[2])?2:(_ff_nonzero(u[1])?1:0)))
#define _deg_v(v)			(_ff_nonzero(v[2])?2:(_ff_nonzero(v[1])?1:(_ff_nonzero(v[0])?0:-1)))
#define _hecurve_set_zero(u,v)		{_ff_set_zero(u[3]);_ff_set_zero(u[2]);_ff_set_zero(u[1]);_ff_set_zero(u[0]);_ff_set_zero(v[2]);_ff_set_zero(v[1]);_ff_set_zero(v[0]);}
#endif


#define _dbg_print_uv(s,u,v)	if ( dbg_level >= DEBUG_LEVEL ) { printf (s);  hecurve_print(u,v); }

// The make functions are independent of genus
void hecurve_make_2 (ff_t u[3], ff_t v[2], ff_t x1, ff_t y1, ff_t x2, ff_t y2);
void hecurve_make_3 (ff_t u[4], ff_t v[3], ff_t x1, ff_t y1, ff_t x2, ff_t y2, ff_t x3, ff_t y3);

// needed to double pts in random 
void hecurve_g2_square_1 (ff_t u[3], ff_t v[2], ff_t u0, ff_t v0, ff_t f[6]);

// note that caller must clear higher coefficients of u and v
static inline void hecurve_make_1 (ff_t u[2], ff_t v[1], ff_t x1, ff_t y1)
{
	_ff_set_one(u[1]);		// this also works in genus 1 Jacobian coords (z=u[1]=1)
	_ff_neg(u[0],x1);
	_ff_set(v[0],y1);
}

void hecurve_setup (mpz_t p)
{
#if HECURVE_VERIFY
	err_printf ("HECURVE VERIFY IS ON\n");
#endif
	ff_setup (p);
}


void hecurve_print (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS])
	{  hecurve_sprint(buf, u, v);  out_printf ("%s\n", buf); }

	
int hecurve_sprint (char *s, ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS])
{
#if HECURVE_GENUS==3
	if ( _ff_nonzero(u[3]) ) {
		if ( _ff_one(u[3]) ) {
			return gmp_sprintf (s, "(x^3 + %Zd*x^2 + %Zd*x + %Zd, %Zd*x^2 + %Zd*x + %Zd)",
							_ff_wrap_mpz(u[2],0),_ff_wrap_mpz(u[1],1), _ff_wrap_mpz(u[0],2), _ff_wrap_mpz(v[2],3), _ff_wrap_mpz(v[1],4), _ff_wrap_mpz(v[0],5));
		} else {
			err_printf ("Warning, non-monic u in hecurve_sprint\n");
			return gmp_sprintf (s, "(%Zd*x^2 + %Zd*x + %Zd, %Zd*x + %Zd)",
							_ff_wrap_mpz(u[3],0),_ff_wrap_mpz(u[2],1),_ff_wrap_mpz(u[1],2), _ff_wrap_mpz(u[0],3), _ff_wrap_mpz(v[2],4),_ff_wrap_mpz(v[1],5), _ff_wrap_mpz(v[0],6));
		}
	}
#endif
#if HECURVE_GENUS >= 2
	if ( _ff_nonzero(u[2]) ) {
#if HECURVE_GENUS == 3
		if ( _ff_nonzero(v[2]) ) err_printf ("Warning, deg v == deg u detected in hecurve_sprint, v[2] = %Zd\n", _ff_wrap_mpz(v[2],0));
#endif
		if ( _ff_one(u[2]) ) {
			return gmp_sprintf (s, "(x^2 + %Zd*x + %Zd, %Zd*x + %Zd)",
							_ff_wrap_mpz(u[1],1), _ff_wrap_mpz(u[0],2), _ff_wrap_mpz(v[1],3), _ff_wrap_mpz(v[0],4));
		} else {
			err_printf ("Warning, non-monic u in hecurve_sprint\n");
			return gmp_sprintf (s, "(%Zd*x^2 + %Zd*x + %Zd, %Zd*x + %Zd)",
							_ff_wrap_mpz(u[2],0),_ff_wrap_mpz(u[1],1), _ff_wrap_mpz(u[0],2), _ff_wrap_mpz(v[1],3), _ff_wrap_mpz(v[0],4));
		}
	} else if ( _ff_nonzero(u[1]) ) {
		if ( _ff_nonzero(v[1]) ) err_printf ("Warning, deg v == deg u = 1 detected in hecurve_sprint, v = %Zdx+%Zd\n", _ff_wrap_mpz(v[1],0),_ff_wrap_mpz(v[0],1));
		if ( _ff_one(u[1]) ) {
			return gmp_sprintf (s, "(x + %Zd, %Zd)", _ff_wrap_mpz(u[0],0), _ff_wrap_mpz(v[0],1));
		} else {
			err_printf ("Warning, non-monic u in hecurve_sprint\n");
			return gmp_sprintf (s, "(%Zdx + %Zd, %Zd)", _ff_wrap_mpz(u[1],0), _ff_wrap_mpz(u[0],1), _ff_wrap_mpz(v[0],2));
		}
	} else if ( _ff_nonzero(u[0]) ) {
		if ( _ff_nonzero(v[1]) || _ff_nonzero(v[0]) ) err_printf ("Warning, deg v >= deg u = 0 detected in hecurve_sprint, v = %Zdx+%Zd\n", _ff_wrap_mpz(v[1],0),_ff_wrap_mpz(v[0],1));
		if ( _ff_one(u[0]) ) {
			return gmp_sprintf (s, "(1, 0)");
		} else {
			err_printf ("Warning, non-monic u in hecurve_sprint\n");
			return gmp_sprintf (s, "(%Zd, 0)", _ff_wrap_mpz(u[0],0));
		}			
	} else {
		err_printf ("u is zero in hecurve_sprint!\n");  exit (0);
	}
#else
	if ( ! _ff_zero(u[1]) && ! _ff_one(u[1]) ) { printf ("non-monic u in hecurve print, lc=%d(raw=%u,one=%u(\n", _ff_get_ui(u[1]),u[1],_ff_mont_R); exit(0); }
	if ( _hecurve_is_identity(u,v) ) {
		return gmp_sprintf (s, "(1, 0)");
	} else {
		return gmp_sprintf(s, "(x + %Zd, %Zd)", _ff_wrap_mpz(u[0],1), _ff_wrap_mpz(v[0],2));
	}
#endif
}


int hecurve_verify (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS], ff_t f[HECURVE_DEGREE+1])
{
	_ff_t_declare t[2*HECURVE_GENUS-1], g[HECURVE_DEGREE+1], x, y, z;
	int i, d_u, d_v, d_t, d_g;

	if ( ! _ff_one (f[HECURVE_DEGREE]) ) { err_printf ("Invalid f polynomial in hecurve_verify, not monic degree %d\n", HECURVE_DEGREE);  exit (0); }
	if ( ! u ) { err_printf ("Null pointer in hecurve_verification!\n");  return 0; }
	
#if HECURVE_GENUS == 1
	if ( _ff_zero(u[1]) ) {
		if ( _ff_zero(v[0]) && _ff_zero(u[0]) ) return 1;
		out_printf ("hecurve verification failed - invalid identity element in genus 1\n");
		return 0;
	}
if ( ! _ff_one(u[1]) ) { out_printf ("hecurve verification failed - non-affine point in genus 1, u[1]=%d\n", _ff_get_ui(u[1])); return 0; }
	_ff_neg (x, u[0]);
	ff_poly_eval (&y, f, HECURVE_DEGREE, &x);
	+ff_square (z, v[0]);
	if ( ! _ff_equal (y,z) ) {
		out_printf ("hecurve verification failed for element:\n ");  hecurve_print(u,v);
		out_printf ("point not on curve\n");
		return 0;
	}
	return 1;
#endif
	
	d_u = _deg_u(u);
	if ( d_u == 0 ) {
		if ( ! _hecurve_is_identity(u,v) ) {
			out_printf ("hecurve verification failed for element:\n  ");  hecurve_print(u,v);
			out_printf ("degree of u is zero, but (u,v) is not the identity.\n");
			return 0;
		}
		return 1;
	}
	if ( ! _ff_one(u[d_u]) ) {
		out_printf ("hecurve verification failed for element:\n  ");  hecurve_print(u,v);
		out_printf ("u is not monic.\n");
		printf ("d_u = %d, u[d_u] = %d (raw) %Zd\n", d_u, u[d_u], _ff_wrap_mpz(u[d_u], 0));
		return 0;
	}
	d_v = _deg_v(v);
	ff_poly_mult (t, &d_t, v, d_v, v, d_v);						// t = v^2
	ff_poly_sub (g, &d_g, f, HECURVE_DEGREE, t, d_t);				// g = f-v^2
	ff_poly_mod (t, &d_t, g, d_g, u, d_u);
	if ( d_t >= 0 ) {
		out_printf ("hecurve verification failed for element: \n    ");  hecurve_print (u,v);
		out_printf ("f-v^2 mod u = ");  ff_poly_print (t, d_t);
		out_printf("degree u = %d, degree v = %d\n", d_u, d_v);
		out_printf("f-v^2 = "); ff_poly_print(g,d_g);
		return 0;
	}
	return 1;
}


int hecurve_random_point (ff_t px[1], ff_t py[1], ff_t f[HECURVE_DEGREE+1])
{
	_ff_t_declare x, y, z, t;
	int i;
#if INIT_REQUIRED
	static int init;
	if ( ! init ) { _ff_init(x);  _ff_init(y);  _ff_init(z);  _ff_init(t);  init = 1; }
#endif
	i = 0;
	do {
		if ( i++ > HECURVE_RANDOM_POINT_RETRIES ) return 0;				// it is possible for very small p that there are no rational points on the curve other than the point at infinity
		_ff_random(x);
///		dbg_printf ("Random x = %Zd in F_%Zd\n",_ff_wrap_mpz(x,0), _ff_mpz_p);
		ff_poly_eval (&y, f, HECURVE_DEGREE, &x);
///		dbg_printf ("y^2 = f(%Zd) = %Zd mod %Zd, computing square root\n",_ff_wrap_mpz(x,0),_ff_wrap_mpz(y,1), _ff_mpz_p);
	} while ( ! ff_sqrt (&z, &y) );
	if ( ui_randomb(1) ) ff_negate (z);	// flip square root 50% of the time
///	dbg_printf ("Found square root %Zd\n", _ff_wrap_mpz(z,0));
	_ff_set(px[0],x);
	_ff_set(py[0],z);
#if ! HECURVE_FAST
	if ( dbg_level >= 2 ) {
		ff_poly_sprint (buf, f, HECURVE_DEGREE);
		dbg_printf ("Random Point (%Zd,%Zd) on y^2 = %s over F_%Zd\n",_ff_wrap_mpz(px[0],0), _ff_wrap_mpz(py[0],1), buf, _ff_mpz_p);
	}
#endif
	return 1;
}

void hecurve_random (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS], ff_t f[HECURVE_DEGREE+1])
{
	_ff_t_declare x1, y1, x2, y2, x3, y3;
	_ff_t_declare alpha[3], beta[3];
#if HECURVE_GENUS == 3
	_ff_t_declare alpha_p[3], beta_p[3], temp1[3], temp2[3], temp3[3];
	_ff_t_declare_reg gamma, delta, b, c, d, D;
#endif
	register unsigned long t;
	unsigned bits;
	int r;

#if HECURVE_GENUS == 1
	if ( ! hecurve_random_point(&x1,&y1, f) ) { _hecurve_set_identity(u,v);  return; }
	hecurve_make_1 (u,v,x1,y1);
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u, v, f) ) {
		err_printf ("Verification failed in hecurve_random!\n");
		err_printf ("P1 = (%Zd,%Zd)\n", _ff_wrap_mpz(x1,0), _ff_wrap_mpz(y1,1));
	}
#endif
#if ! HECURVE_FAST
	if ( dbg_level >= 2 ) { printf ("Random element: ");  hecurve_print (u, v); }
#endif
	return;
#endif
	
#if HECURVE_GENUS == 2
#if HECURVE_RANDOM == 1
	do {
		_ff_random(u[2]);
		if ( _ff_zero(u[2]) ) {						// use weight one divisor with probability 1/p
			_hecurve_set_zero(u,v);
			if ( !  hecurve_random_point(&x1,&y1, f) ) { _hecurve_set_identity(u,v);  return; }
			hecurve_make_1 (u, v, x1, y1);
			break;
		}
		_ff_set_one(u[2]);
		_ff_random(u[1]);
		_ff_random(u[0]);
		bits = (unsigned) ui_randomm(4);
	} while ( ! hecurve_unbits (v,u,bits,f) );
#else
#if HECURVE_RANDOM == 2
	_hecurve_set_zero (u, v);
	// Flip a coin to decide how to split u(x)
	if ( ui_randomb(1) ) {
		// u will be reducible.  pick two points on the curve
		if ( ! hecurve_random_point (&x1, &y1, f) ) { _hecurve_set_identity(u,v);  return; }		// identity return only happens for curves with no rational pts
		if ( ! hecurve_random_point (&x2, &y2, f) ) { _hecurve_set_identity(u,v);  return; }
		if ( _ff_equal(x1,x2) ) {
			if ( _ff_equal(y1,y2) ) {
				// if the points are identical, double the pt, unless it has 2-torsion
				if ( _ff_zero(y1) ) {
					hecurve_make_1 (u, v, x1, y1);		// for 2-torsion points don't double (means we never return identity here)
//					puts ("singleton 2-torsion");
				} else {
					ff_negate(x1);						// for single pt divisor, u(x)=x-x1 so u0=-x1
					hecurve_g2_square_1 (u, v, x1, y1, f);
//					puts ("double pt");
				}
			} else {
				// Treat the case of a pt and its opposite by just making a divisor with 1 pt.
				// This will occur with probability \approx 1/p, which is right ratio
				hecurve_make_1 (u, v, x1, y1);
//				puts ("singleton non-torsion");
			}
		} else {
			hecurve_make_2 (u, v, x1, y1, x2, y2);
//			puts ("pair of pts");
		}
	} else {
		// make u irreducible by picking a root alpha of degree 2 over F_p and setting u(x)=(x-alpha)(x-phi(alpha))
		// where phi is the frobenius map.  If alpha = a_0+a_1z  \in F_p[z]/(z^2-nr), then phi(alpha) = a_0-a_1z
		// note that if beta^2 = f(alpha) then phi(beta)^2 = f(phi(alpha)), so the divisor necessarily contains the
		// two points (alpha,+/-beta) and (phi(alpha),+/-phi(beta)).  The signs must match, otherwise v(x) won't be rational.
		for(;;) {
			do { _ff_random(alpha[1]); } while ( _ff_zero(alpha[1]) );
			_ff_random(alpha[0]);
//printf ("alpha=%ldz+%ld\n", _ff_get_ui(alpha[1]), _ff_get_ui(alpha[0]));
			ff2_poly_eval(beta, f, HECURVE_DEGREE, alpha);
			if ( ! ff2_sqrt(beta, beta) ) continue;
			if ( ui_randomb(1) ) ff2_negate (beta);	// flip square root 50% of the time
//printf ("beta=%ldz+%ld\n", _ff_get_ui(beta[1]), _ff_get_ui(beta[0]));
			_ff_set_one(u[2]);
			_ff_add(x1,alpha[0],alpha[0]);
			_ff_neg(u[1],x1);
			_ff_square(x1,alpha[0]);
			_ff_square(x2,alpha[1]);
			_ff_mult(x3,x2,_ff_nr);
			_ff_sub(u[0],x1,x3);
			_ff_invert(x1,alpha[1]);
			_ff_mult(v[1],x1,beta[1]);
			_ff_mult(x2,alpha[0],v[1]);
			_ff_sub(v[0],beta[0],x2);
			// cost to get (u,v) from (alpha,beta) is I+5M+4A, counting squares as mults and subs/negs as additions.
//			puts ("2 non-rational conjugates");
			break;
		}
	}
	
#else
	_hecurve_set_zero (u, v);
	if ( ! hecurve_random_point (&x1, &y1, f) ) { _hecurve_set_identity(u,v);  return; }
	if ( !  hecurve_random_point (&x2, &y2, f) ) { _hecurve_set_identity(u,v);  return; }
	if ( _ff_equal(x1,x2) ) {
		hecurve_make_1 (u, v, x1, y1);
	} else {
		hecurve_make_2 (u, v, x1, y1, x2, y2);
	}
#endif
#endif
#endif
#if HECURVE_GENUS == 3
	_hecurve_set_zero (u, v);
#if HECURVE_RANDOM == 2
	
	// Flip a 6-sided coin to decide how to split u(x)
retry3:
	r = ui_randomm(6);
	if ( !r ) {
		// with probabilty 1/6, split u completely
		for(;;) {
			_hecurve_set_zero (u, v);
			// u will be split completely.  pick three points on the curve
			if ( ! hecurve_random_point (&x1, &y1, f) ) goto retry3;								// only happens for curves with no rational pts
			if ( ! hecurve_random_point (&x2, &y2, f) ) goto retry3;
			if ( ! hecurve_random_point (&x3, &y3, f) ) goto retry3;
			if ( _ff_equal(x1,x2) ) {
				if ( _ff_equal(x1,x3) || ! _ff_equal(y1,y2) ) {
					// If all three x values are equal, let P=(x1,y1) and we create P, 2P, or 3P with equal probabilities
					// Each case should have probability 1/p^3, however the probability we have reached this point with a particular P is (1/6)*(4/p^3) = 2/(3p^3)
					// and we really want it to be 3/p^3.  To fix things up, we include some of the cases where (x1,y1) and (x2,y2) are opposite points
					if ( ! _ff_equal(y1,y2) ) {
						r = ui_randomm(3*_ff_p);
						if ( r >= 7 ) continue;
						// we get here with probability 7/(3p^3), which added to 2/(3p^3) (all three x values equal) gives the desired 3/p^3.
					}
					r = ui_randomm(3);
					switch (r) {
					case 0:
						hecurve_make_1 (u, v, x1, y1);
//						puts ("single point");
						break;
					case 1:
						hecurve_make_1 (u, v, x1, y1);
						hecurve_g3_square (u, v, u, v, f, 0);
//						puts ("double point");
						break;
					case 2:
						hecurve_make_1 (u, v, x1, y1);
						hecurve_g3_square (alpha, beta, u, v, f, 0);
						hecurve_g3_compose (u, v, alpha, beta, u, v, f, 0);
//						puts ("triple point");
					}
				} else {
					// we have two copies of a pt and a third distinct pt
					_hecurve_set_zero(alpha,beta);
					hecurve_make_2 (alpha, beta, x1, y1, x3, y3);	// compose distinct pts
					hecurve_make_1 (u, v, x2, y2);				// make one of the dups singleton
					hecurve_compose (u, v, alpha, beta, u, v, f, 0);	// put it all together
//					puts ("2+1 point");
				}
			} else if ( _ff_equal (x1, x3) || _ff_equal(x2,x3) ) {
				// points 1 and 2 are distinct, but 3 is a duplicate or an opposite
				// We have handled all cases except 2 distinct rational points, which we handle here if (x1,y1)==(x3,y3)
				// otherwise we retry in order to avoid overweighting
				if (!  (_ff_equal(x1,x3) && _ff_equal(y1,y3)) ) continue;
				hecurve_make_2 (u, v, x1, y1, x2, y2);
//				puts ("2 distinct rational points");
			} else {
				hecurve_make_3 (u, v, x1, y1, x2, y2, x3, y3);		// usual case
//				puts ("3 distinct rational points");
			}
			break;
		}
	} else if ( r < 4 ) {
		// with probability 1/2, split u(x) into a linear and a quadratic, in general,
		// but with probaility 1/p (within 1/2) don't include a linear factor (or if there are none)
		if ( ! ui_randomm(_ff_p+1) || ! hecurve_random_point (&x1, &y1, f) ) {
			// this is an exact duplicate of the genus 2 code above
			for(;;) {
				do { _ff_random(alpha[1]); } while ( _ff_zero(alpha[1]) );
				_ff_random(alpha[0]);
//printf ("alpha=%ldz+%ld\n", _ff_get_ui(alpha[1]), _ff_get_ui(alpha[0]));
				ff2_poly_eval(beta, f, HECURVE_DEGREE, alpha);
				if ( ! ff2_sqrt(beta, beta) ) continue;
				if ( ui_randomb(1) ) ff2_negate (beta); 		// flip square root 50% of the time
//printf ("beta=%ldz+%ld\n", _ff_get_ui(beta[1]), _ff_get_ui(beta[0]));
				_ff_set_one(u[2]);
				_ff_add(x1,alpha[0],alpha[0]);
				_ff_neg(u[1],x1);
				_ff_square(x1,alpha[0]);
				_ff_square(x2,alpha[1]);
				_ff_mult(x3,x2,_ff_nr);
				_ff_sub(u[0],x1,x3);
				_ff_invert(x1,alpha[1]);
				_ff_mult(v[1],x1,beta[1]);
				_ff_mult(x2,alpha[0],v[1]);
				_ff_sub(v[0],beta[0],x2);
				// cost to get (u,v) from (alpha,beta) is I+5M+4A, counting squares as mults and subs/negs as additions.
//				puts ("2 non-rational conjugates");
				break;
			}
		} else {
		// split u(x) into a linear and a quadratic
//printf("x1=%ld, y1=%ld, s=%ld, p=%ld\n", _ff_get_ui(x1), _ff_get_ui(y1), _ff_get_ui(_ff_nr), _ff_p);
			for(;;) {
				do { _ff_random(alpha[1]); } while ( _ff_zero(alpha[1]) );
				_ff_random(alpha[0]);
	//printf ("alpha=%ldz+%ld\n", _ff_get_ui(alpha[1]), _ff_get_ui(alpha[0]));
				ff2_poly_eval(beta, f, HECURVE_DEGREE, alpha);
				if ( ! ff2_sqrt(beta, beta) ) continue;
				if ( ui_randomb(1) ) ff2_negate (beta);	// flip square root 50% of the time
	//printf ("beta=%ldz+%ld\n", _ff_get_ui(beta[1]), _ff_get_ui(beta[0]));
				_ff_add(b,alpha[0],alpha[0]);
				ff_negate(b);						// b = -tr(alpha) is coeff of x in quadratic factor
	//printf("b=%ld\n", _ff_get_ui(b));
				_ff_square(c,alpha[0]);
				_ff_square(y2,alpha[1]);
				_ff_mult(x2,y2,_ff_nr);
				_ff_add(d,c,x2);					// d = a0^2+s*a1^2 used below to compute delta
				_ff_subfrom(c,x2);					// c = N(alpha) is const. coeff in quadratic factor
	//printf("c=%ld\n", _ff_get_ui(c));
				_ff_mult(D,b,x1);
				_ff_set_one(u[3]);					// u[3] = 1
				_ff_sub(u[2],b,x1);					// u[2] = b-x1
				_ff_sub(u[1],c,D);					// u[1] = c-bx1
				_ff_mult(y3,c,x1);	
				_ff_neg(u[0],y3);					// u[0] = -cx1
				_ff_square(x2,x1);
				_ff_addto(D,x2);
				_ff_addto(D,c);
				ff_mult(y3,D,alpha[1]);
	//printf("D=%ld\n", _ff_get_ui(D));
				ff_invert(D,y3);						// D = 2z/ [(x1^2+bx1+c)(alpha-phi(alpha))] 		note the 2's all cancel
				_ff_mult(gamma,alpha[1],beta[0]);
				_ff_mult(delta,alpha[0],gamma);
				_ff_mult(y3,beta[1],alpha[0]);
				_ff_subfrom(gamma,y3);				// gamma = alpha*phi(beta),  only need gamma[1]
	//printf("gamma=%ld\n", _ff_get_ui(gamma));
				_ff_x2(delta);
				_ff_mult(y3,beta[1],d);
				_ff_subfrom(delta,y3);				// delta = alpha^2*phi(beta)	only need delta[1]
	//printf("delta=%ld\n", _ff_get_ui(delta));
				ff_mult(y1,alpha[1],y1);				// all use of y1 from here on out will include the factor alpha[1]=(alpha-phi(alpha))/2
				_ff_sub(y2,y1,gamma);
				_ff_mult(y3,beta[1],x1);
				_ff_subfrom(y2,y3);
				_ff_mult(v[2],D,y2);					// v[2] = [alpha1*y1 - gamma1 - beta1*x1] / D
				_ff_mult(y2,x2,beta[1]);
				_ff_addto(y2,delta);
				_ff_mult(y3,b,y1);
				_ff_addto(y2,y3);
				_ff_mult(v[1],D,y2);					// v[1] = [alpha1*y1*b + delta1+ beta1*x1^2] / D
				_ff_mult(y2,c,y1);
				_ff_mult(y3,x1,delta);
				_ff_subfrom(y2,y3);
				_ff_mult(y3,gamma,x2);
				_ff_addto(y2,y3);
				_ff_mult(v[0],D,y2);					// v[0] = [alpha1*y1*c - delta1*x1 + gamm1*x1^2] / D
				// Cost to get (u,v) from (alpha,beta) is I+21M+18A
//				puts ("2,1 non-rational conjugates and rational pt");
				break;
			}
		}
		// construct quadratic factor as in genus 2 case above.
	} else {
		// with probability 1/3 make u(x) irreducible
		for(;;) {
			do { ff3_random(alpha); } while (  _ff_zero(alpha[1]) && _ff_zero(alpha[2]) );
//printf ("alpha=%ldz^2+%ldz+%ld\n", _ff_get_ui(alpha[2]), _ff_get_ui(alpha[1]), _ff_get_ui(alpha[0]));
			ff3_poly_eval(beta, f, HECURVE_DEGREE, alpha);
//printf ("f(alpha)=%ldz^2+%ldz+%ld\n", _ff_get_ui( beta[2]), _ff_get_ui( beta[1]), _ff_get_ui( beta[0]));
			if ( ! ff3_sqrt(beta, beta) ) continue;
			if ( ui_randomb(1) ) ff3_negate (beta);	// flip square root 50% of the time
//printf ("beta=%ldz^2+%ldz+%ld\n", _ff_get_ui( beta[2]), _ff_get_ui( beta[1]), _ff_get_ui( beta[0]));
			ff3_exp_p(alpha_p,alpha);
			ff3_mult(temp1,alpha,alpha_p);
			ff3_trace(&x1,alpha);
			ff3_norm(&x2,alpha);
			_ff_set_one(u[3]);				// u[3] = 1
			_ff_neg(u[2],x1);				// u[2]=-tr(alpha)
			ff3_trace(u+1,temp1);			// u[1]=tr(alpha^(p+1))
			_ff_neg(u[0],x2);				// u[0]=-norm(alpha)
			ff3_mult(temp2,temp1,alpha_p);	// temp2=alpha^(2p+1)
			ff3_mult(temp1,temp1,alpha);		// temp1=alpha^(p+2)
			ff3_trace(&x1,temp2);
			ff3_trace(&x2,temp1);
			_ff_sub(x3,x1,x2);				
			ff_invert(D,x3);					// D = 1 / [(alpha-alpha^p)(alpha^p-alpha^(p^2))(alpha^(p^2)-alpha)]
			ff3_exp_p(beta_p,beta);
			ff3_exp_p(beta,beta_p);			// beta now set to beta^2 (we don't need beta anymore)
			ff3_mult(temp3,temp2,beta);
			ff3_trace(&x1,temp3);
			ff3_mult(temp3,temp1,beta);
			ff3_trace(&x2,temp3);
			_ff_sub(x3,x1,x2);
			_ff_mult(v[0],D,x3);				// v[0] = D * [tr(alpha^(2p+1)*beta^(p^2)) - tr(alpha^(p+2)*beta^(p^2))]
			ff3_square(temp1,alpha);
			ff3_mult(temp2,temp1,beta);
			ff3_trace(&x1,temp2);
			ff3_mult(temp3,temp1,beta_p);
			ff3_trace(&x2,temp3);
			_ff_sub(x3,x1,x2);
			_ff_mult(v[1],D,x3);				// v[1] = D * [tr(alpha^2*(beta^(p^2)) - tr(alpha^2*beta^p)]
			ff3_mult(temp2,alpha,beta_p);
			ff3_trace(&x1,temp2);
			ff3_mult(temp3,alpha,beta);
			ff3_trace(&x2,temp3);
			_ff_sub(x3,x1,x2);
			_ff_mult(v[2],D,x3);				// v[2] = D * [tr(alpha*beta^p) = tr(alpha*beta^(p^2))]
//			puts ("3 non-rational conjugates");
			break;
		}
	}
#else
	if ( ! hecurve_random_point (&x1, &y1, f) ) { _hecurve_set_identity(u,v);  return; }
	if ( !  hecurve_random_point (&x2, &y2, f) ) { _hecurve_set_identity(u,v);  return; }
	if ( ! hecurve_random_point (&x3, &y3, f) ) { _hecurve_set_identity(u,v);  return; }
	if ( _ff_equal(x1,x2) ) {
		if ( _ff_equal (x1,x3) ) {
			hecurve_make_1 (u, v, x1, y1);
		} else {
			hecurve_make_2 (u, v, x1, y1, x3, y3);
		}
	} else {
		if ( _ff_equal(x1,x3) || _ff_equal(x2,x3) ) {
			hecurve_make_2 (u, v, x1, y1, x2, y2);
		} else {
			hecurve_make_3 (u, v, x1, y1, x2, y2, x3, y3);			// usual case
		}
	}
#endif
#endif
#if ! HECURVE_FAST
	if ( dbg_level >= 2 ) { printf ("Random element: ");  hecurve_print (u, v); }
#endif
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u, v, f) ) {
		err_printf ("Verification failed in hecurve_random!\n");
#if HECURVE_GENUS == 3
		err_printf ("P1 = (%Zd,%Zd), P2 = (%Zd,%Zd), P3 = (%Zd,%Zd)\n",
				_ff_wrap_mpz(x1,0), _ff_wrap_mpz(y1,1), _ff_wrap_mpz(x2,2), _ff_wrap_mpz(y2,3), _ff_wrap_mpz(x3,4), _ff_wrap_mpz(y3,5));
#else
		err_printf ("P1 = (%Zd,%Zd), P2 = (%Zd,%Zd)\n", _ff_wrap_mpz(x1,0), _ff_wrap_mpz(y1,1), _ff_wrap_mpz(x2,2), _ff_wrap_mpz(y2,3));
#endif
		exit (0);
	}
#endif
}

//  See formulas on p.307 of [HECHECC]
void hecurve_make_3 (ff_t u[4], ff_t v[3], ff_t x1, ff_t y1, ff_t x2, ff_t y2, ff_t x3, ff_t y3)
{
	_ff_t_declare_reg s1, s2, s3, t1, t2, t3, t, w, z;
#if ! HECURVE_FAST
	if ( _ff_equal(x1,x2) ) { err_printf ("Call to hecurve_make_3 with x1 = x2!\n");  exit (0); }
	if ( _ff_equal(x2,x3) ) { err_printf ("Call to hecurve_make_3 with x2 = x3!\n");  exit (0); }
	if ( _ff_equal(x3,x1) ) { err_printf ("Call to hecurve_make_3 with x3 = x1!\n");  exit (0); }
#endif
	// Construct u(x) = (x-x1)(x-x2)(x-x3) with roots x1, x2, and x3 then construct v(x) so that v(x_i) = y_i ensuring u|v^2-f (note h =0)
	_ff_set_one(u[3]);						// u[3] = 1
	_ff_sub(t1,x1,x2);						// t1=x1-x2
	_ff_sub(t2,x2,x3);						// t2=x2-x3
	_ff_sub(t3,x3,x1);						// t3=x3-x1
	_ff_mult(t,t1,t2);
	_ff_mult(z,t,t3);
	_ff_invert(z,z);
	ff_negate(z);							// z = -1/[(x1-x2)(x2-x3)(x3-x1)]
	ff_mult(t1,t1,y3);						// t1=(x1-x2)y3
	ff_mult(t2,t2,y1);						// t2=(x2-x3)y1
	ff_mult(t3,t3,y2);						// t3=(x3-x1)y2
	_ff_add(t,t1,t2);
	_ff_addto(t,t3);
	ff_mult(v[2],z,t);						// v[2] = -[(x1-x2)y3+(x2-x3)y1+(x3-x1)y2]/[(x1-x2)(x2-x3)(x3-x1)]
	_ff_add(s1,x1,x2);						// s1=x1+x2
	_ff_add(s2,x2,x3);						// s2=x2+x3
	_ff_add(s3,x3,x1);						// s3=x3+x1
	_ff_add(t,s1,x3);
	_ff_neg(u[2],t);						// u[2] = -(x1+x2+x3)
	_ff_mult(t,s1,t1);
	_ff_mult(w,s2,t2);
	_ff_addto(t,w);
	_ff_mult(w,s3,t3);
	_ff_addto(t,w);
	ff_negate(t);
	ff_mult(v[1],z,t);						// v[1] = [(x1-x2)(x1+x2)y3+(x2-x3)(x2+x3)y1+(x3-x1)(x3+x1)y2]/[(x1-x2)(x2-x3)(x3-x1)]
	_ff_mult(s1,x1,x2);						// s1=x1x2
	_ff_mult(s2,x2,x3);						// s2=x2x3
	_ff_mult(s3,x3,x1);						// s3=x3x1
	_ff_mult(t,s1,x3);
	_ff_neg(u[0],t);						// u[0] = -x1x2x3
	_ff_add(t,s1,s2);
	_ff_add(u[1],t,s3);						// u[1] = x1x2+x2x3+x3x1
	_ff_mult(t,s1,t1);
	_ff_mult(w,s2,t2);
	_ff_addto(t,w);
	_ff_mult(w,s3,t3);
	_ff_addto(t,w);
	_ff_mult(v[0],z,t);						// v[0] = -[x1x2(x1-x2)y3+x2x3(x2-x3)y1)+x3x1(x3-x1)y2]/[(x1-x2)(x2-x3)(x3-x1)]
	// Total cost is I+18M+16A (counting negations and subtractions as additions)
}


// note this function is not genus specific, however caller is responsible for zeroing higher degree coefficients
void hecurve_make_2 (ff_t u[3], ff_t v[2], ff_t x1, ff_t y1, ff_t x2, ff_t y2)
{
	_ff_t_declare_reg t1, t2, t3;
#if ! HECURVE_FAST
	if ( _ff_equal(x1,x2) ) { err_printf ("Call to hecurve_make_2 with x1 = x2!\n");  exit (0); }
#endif
	// Construct u(x) = (x-x_1)(x-x_2) has roots x1 and x2 then construct v(x) so that v(x_i) = y_i ensuring u|v^2-f (note h =0)
	_ff_set_one(u[2]);						// u[2] = 1
	_ff_add(t1,x1,x2);
	_ff_neg(u[1],t1);						// u[1] = -(x1+x2)
	_ff_mult(u[0],x1,x2);						
	_ff_sub(t1,x1,x2);
	_ff_invert(t1,t1);						// t1 = 1/(x1-x2) 
	_ff_sub(t2,y1,y2);
	_ff_mult(v[1],t2,t1);						// v1 = (y1-y1)/(x1-x2)
	_ff_mult(t3,x1,y2);
	_ff_mult(t2,x2,y1);
	_ff_subfrom(t3,t2);
	_ff_mult(v[0],t3,t1);						// v0 = (x1y2-x2y1)/(x1-x2)
	// total cost is I+5M+5A (counting subtractions and negations as additions)
}

/*
	Basic implementation of Cantor's algorithm.  Not optimized in the least - only intended for use
	in exceptional cases.  These generally are quite infrequent, however for very small p this
        algorithm sees heavy use (in fact, smalljac typically runs faster with p > 1000 than p < 1000 for this reason).

	Based on Algorithm 14.7 of [HECHECC].
*/
void hecurve_compose_cantor (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS],
						  ff_t u1[HECURVE_GENUS+1], ff_t v1[HECURVE_GENUS],
						  ff_t u2[HECURVE_GENUS+1], ff_t v2[HECURVE_GENUS], ff_t f[HECURVE_DEGREE+1])
{
	_ff_t_declare d1[HECURVE_GENUS+1], e1[HECURVE_GENUS+1], e2[HECURVE_GENUS+1];
	_ff_t_declare d[HECURVE_GENUS], c1[HECURVE_GENUS], c2[HECURVE_GENUS];
	_ff_t_declare s1[POLY_MAX_DEGREE+1], s2[POLY_MAX_DEGREE+1], s3[POLY_MAX_DEGREE+1];
	_ff_t_declare q[POLY_MAX_DEGREE+1], r[POLY_MAX_DEGREE+1], s[POLY_MAX_DEGREE+1], t[POLY_MAX_DEGREE+1];
	int i, d_d1, d_e1, d_e2, d_d, d_c1, d_c2, d_s, d_t, d_s1, d_s2, d_s3, d_u1, d_u2, d_v1, d_v2, d_q, d_r, d_f, d_u, d_v;

	// assume no initialization of ff_t's is required - too much of a hassle
#if FF_MPZ
	err_printf ("Cantor's algorithm not supported for MPZ field implementation\n");
	exit (0);
#endif

//	out_printf ("Cantor composition:\n");
//	out_printf ("    (u1,v1) = ");  hecurve_print (u1,v1);
//	out_printf ("    (u2,v2) = ");  hecurve_print (u2,v2);
	
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u1, v1, f) )  { err_printf ("hecurve_compose_cantor: invalid input\n");  exit (0); }
	if ( ! hecurve_verify (u2, v2, f) ) { err_printf ("hecurve_compose_cantor: invalid input\n");  exit (0); }
#endif
	
	d_u1 = _deg_u(u1);    d_u2 = _deg_u(u2);    d_v1 = _deg_v(v1);   d_v2 = _deg_v(v2);  d_f = HECURVE_DEGREE;
	
	_poly_gcdext (d1, e1, e2, u1, u2);
	_poly_add (t, v1, v2);
	_poly_gcdext (d, c1, c2, d1, t);
	_poly_mult(s1,c1,e1);  _poly_mult(s2,c1,e2);  _poly_set(s3,c2);
	_poly_mult (q, u1, u2);
	_poly_mult (t, d, d);
	_poly_div (s, r, q, t);			// s = u1u2/d^2		(s = u)
	if ( d_r != -1 ) { err_printf ("Inexact division in Cantor's algorithm, r = ");  _poly_print(r);  exit (0); }
	_poly_mult (r, u1, v2);
	_poly_mult (t, r, s1);
	_poly_mult (r, u2, v1);
	_poly_mult (q, r, s2);
	_poly_addto (t, q);
	_poly_mult (r, v1, v2);
	_poly_addto (r, f);
	_poly_mult (q, r, s3);
	_poly_addto (t, q);
	_poly_div (q, r, t, d);
	if ( d_r != -1 ) { err_printf ("Inexact division in Cantor's algorithm, r = ");  _poly_print(r);  exit (0); }
	_poly_mod (t, q, s);			// t = (s1u1v1+s2u2v1+s3(v1v2+f))/d mod u
	do {
		_poly_mult (s3, t, t);
		_poly_sub (s2, f, s3);
		_poly_div (s1, r, s2, s);	// s1 = (f-v^2)/u	(s1 = u')
		_poly_set (s, s1);
		_poly_mod (s2, t, s1);
		_poly_neg (t, s2);
	} while ( d_s > HECURVE_GENUS );
	_poly_monic (u, s);
	_poly_set (v, t);
	for ( i = d_u+1 ; i <= HECURVE_GENUS ; i++ ) _ff_set_zero(u[i]);
	for ( i = d_v+1 ; i < HECURVE_GENUS ; i++ ) _ff_set_zero(v[i]);
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u, v, f) ) {
		err_printf ("hecurve_compose_cantor output verification failed for inputs:\n");
		err_printf ("      ");  hecurve_print(u1,v1);
		err_printf ("      ");  hecurve_print(u2,v2);
		err_printf ("note that inputs have been modified if output overlapped.\n");
		exit (0);
	}
#endif	
}

/*
	compression/decompression currently only supported for genus 1 and 2.  Note that
	decompression is used for random element generation.

	See p. 312 of [HECHECC] for a description of compression in genus 2.
*/


#if HECURVE_GENUS <= 2

unsigned hecurve_bits (ff_t u[3], ff_t v[2], ff_t f[6])
{
	_ff_t_declare_reg s0, t1, t2;
	register unsigned long bits;
#if HECURVE_GENUS == 1
	return ( _ff_parity(v[0]) ? 1 : 0 );
#endif

	// In genus 2, compute low order bit of s0 and low order bit of v0 (or v1 when v0 = 0)
	
	_ff_square (t1, u[1]);
	_ff_add(t2, u[0], u[0]);
	_ff_subfrom (t2, t1);
	_ff_mult (t1, t2, u[1]);
	_ff_square (s0, v[1]);
	_ff_subfrom (s0, t1);
#if ! HECURVE_SPARSE
	if ( _ff_nonzero(f[3]) ) {
		_ff_mult (t1, f[3], u[1]);
		_ff_addto (s0, t1);
	}
	if ( _ff_nonzero(f[2]) ) _ff_subfrom (s0, f[2]);
#endif
	bits = _ff_parity(s0);
	if ( v[0] ) {
		if ( _ff_parity(v[0]) ) bits |= 0x2;
	} else {
		if ( _ff_parity(v[1]) ) bits |= 0x2;
	}
	return bits;
}


// Used for genus 2 decompression and random element generation - attempts to reconstruct v from u and 2 bits
// No comparable support exists in genus 3 - decompression is much more difficult for g > 2
// This code needs to handle f4 != 0, since it may get used over F_5
int hecurve_unbits (ff_t v[2], ff_t u[3], unsigned bits, ff_t f[6])
{
	_ff_t_declare_reg s, A, B, a, b, c, t0, t1, t2,w1,w2,w3,w4;
	_ff_t_declare T0, T1;
		
#if HECURVE_GENUS == 1
	_ff_square (t1, u[0]);
	_ff_addto (t1, f[1]);
	_ff_sub (T0, f[0], t1);
	if ( ! ff_sqrt(v,&T0) ) { /*dbg_printf ("unbits: sqrt(f0-f1+u0^2) = sqrt(%Zd) failed\n", _ff_wrap_mpz(t0,0)); */ return 0; }
	if ( (bits&0x1) != _ff_parity(v[0]) ) ff_negate(v[0]);
	return 1;
#endif
		
	_ff_square(w1,u[1]);			// save w1 = u[1]^2
	_ff_add(t1,u[0],u[0]);
	_ff_add(w4,t1,t1);				// save w4 = 4u[0]
	_ff_subfrom(t1,w1);
	_ff_mult(w2,u[1],t1);			// save w2 = u[1](2u[0]-u[1]^2)
	_ff_sub(t1,u[0],w1);
	_ff_mult(w3,u[0],t1);			// save w3 = u[0](u[0]-u[1]^2)
#if ! HECURVE_SPARSE
	if ( _ff_nonzero(f[3]) ) {
		_ff_mult(t1,f[3],u[1]);
		_ff_sub(A,f[2],t1);
		_ff_addto(A,w2);
		_ff_mult(t1,f[3],u[0]);
		_ff_sub(B,f[1],t1);
		_ff_addto(B,w3);
	} else {
		_ff_add(A,w2,f[2]);
		_ff_add(B,w3,f[1]);
	}
	if ( _ff_nonzero(f[4]) ) {
		_ff_sub(t1,u[0],w1);
		ff_mult(t1,t1,f[4]);
		_ff_subfrom(A,t1);
		_ff_mult(t1,u[0],u[1]);
		ff_mult(t1,t1,f[4]);
		_ff_addto(B,t1);
	}
#else
	_ff_set(A,w2);
	_ff_add(B,f[1],w3);
#endif
	// (14.5) on p. 312 is: (s+A)x^2 + (u[1]s+B)x + (u[0]s+f[0]) which must have discriminant 0 (note s = s0)
	// this implies (u[1]s+B)^2 - 4(s+A)(u[0]s+f[0]) = 0
	// thus (u[1]^2-4u[0])s^2 + (2Bu[1]-4(Au[0]+f[0]))s + (B^2-4Af[0]) = as^2 + bs + c = 0
	// we now compute a, b, and c
	_ff_sub(a,w1,w4);				// a = u[1]^2-4u[0]
	_ff_multadd (t0, A, u[0], f[0]);
	_ff_x2(t0);  _ff_x2(t0);
	_ff_mult(t1,B,u[1]);
	_ff_add(b,t1,t1);
	_ff_subfrom(b,t0);				// b = 2Bu[1] - 4(Au[0]+f[0])
	_ff_square(c,B);
	_ff_mult(t0,A,f[0]);
	_ff_x2(t0);  _ff_x2(t0);
	_ff_subfrom(c,t0);				// c = B^2 - 4Af[0]
	if ( _ff_zero(a) ) {
		if ( _ff_zero(b) ) { /*dbg_printf ("unbits: a=b=0\n"); */ return 0; }
		_ff_invert(t0,b);
		_ff_neg(t1,c);
		_ff_mult(s,t0,t1);
	} else {
		_ff_mult(t0,a,c);
		_ff_add(t1,t0,t0);  _ff_x2(t1);
		_ff_square(T0,b);
		_ff_subfrom(T0,t1);
		if ( ! ff_sqrt(&T1,&T0) ) { /*dbg_printf ("unbits: sqrt(b^2-4ac) = sqrt(%Zd) failed\n", _ff_wrap_mpz(t1,0));*/  return 0; }
		_ff_add(t0,a,a);
		ff_invert(t0,t0);
		_ff_sub(s,T1,b);
		ff_mult(s,s,t0);					// s = (-b + sqrt(b^2-4ac))/2a
		if ( bits&0x1 != _ff_parity(s) ) {	// 50/50 chance here
			_ff_neg(s,T1);
			_ff_subfrom(s,b);
			ff_mult(s,s,t0);
		}
	}
	// we have s, now we use (14.6), (14.7), and (14.8) on p.312 of HHECC to compute v[0] and v[1]
	// First use (14.6) to compute v0 = sqrt(u0s0+f0)
	_ff_multadd(T0,s,u[0],f[0]);
	if ( ! ff_sqrt(v,&T0) ) { /*dbg_printf ("unbits: sqrt(u[1]s0+f[0]) = sqrt(%Zd) failed\n", _ff_wrap_mpz(t0, 0));*/   return 0; }
	if ( _ff_zero(v[0]) ) {
		// use (14.8) to calculate v[1] when v[0] is zero - note that 14.8 is just v[1]^2 = s+A where A was computed above
		_ff_add(T0,s,A);
		if ( ! ff_sqrt(v+1,&T0) ) { /*dbg_printf ("unbits: sqrt(s0+f[2]-f[3]u[1]+u[1](2u[0]-u[1]^2) = sqrt(%Zd) failed\n", _ff_wrap_mpz(t0,0));*/  return 0; }
		if ( ((bits&0x2)>>1) != _ff_parity(v[1]) ) ff_negate(v[1]);
	} else {
		// got v[0], set sign appropriately and use (14.7) to calculuate v[1] = (u[1]s+B)/(2v[0]) where B was computed above
		if ( ((bits&0x2)>>1) != _ff_parity(v[0]) ) ff_negate(v[0]);
		_ff_add(t0,v[0],v[0]);
		ff_invert(t0,t0);
		_ff_multadd(t1,u[1],s,B);
		_ff_mult(v[1],t0,t1);
	}
	return 1;
}

#endif

