#include <stdlib.h>
#include <stdio.h>
#include "ffext.h"
#include "ffpoly.h"
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
	This module contains performance-oriented code for operations on polynomials of low degree,
	and support functions for computations on genus 1 and genus 2 curves.
	
	It is assumed here that the ff_t type does not require initialization (no support for GMP's mpz type),
	and in a number of cases it is further assumed that the ff_t type fits in an unsigned long (this assumption is verified when made).
*/


void ff_poly_xpan_mod_d2 (ff_t g[2], ff_t a, unsigned long n, ff_t f[1]);
void ff_poly_xn_mod_d3 (ff_t g[3], unsigned long n, ff_t f[4]);
void ff_poly_xpan_mod_d3 (ff_t g[2], ff_t a, unsigned long n, ff_t f[4]);
void ff_poly_mult_mod_d3 (ff_t w[3], ff_t u[3], ff_t v[3], ff_t f[4]);
void ff_poly_square_mod_d3 (ff_t w[3], ff_t u[3], ff_t f[4]);
void ff_poly_exp4_mod_d3 (ff_t g[3], ff_t p[12], unsigned long n, ff_t f[4]);
void ff_poly_exp8_mod_d3 (ff_t g[3], ff_t p[24], unsigned long n, ff_t f[4]);

void ff_poly_xn_mod_d4 (ff_t g[4], unsigned long n, ff_t f[5]);
void ff_poly_xpan_mod_d4 (ff_t g[3], ff_t a, unsigned long n, ff_t f[5]);
void ff_poly_mult_mod_d4 (ff_t w[4], ff_t u[4], ff_t v[4], ff_t f[5]);
void ff_poly_square_mod_d4 (ff_t w[4], ff_t u[4], ff_t f[5]);
void ff_poly_exp4_mod_d4 (ff_t g[4], ff_t p[16], unsigned long n, ff_t f[5]);
void ff_poly_exp8_mod_d4 (ff_t g[4], ff_t p[32], unsigned long n, ff_t f[5]);


// assumes f monic, char > 3, sets t to translation value: depressed cubic is f(x-t)
void ff_depress_cubic (ff_t t[1], ff_t f[4])
{
	_ff_t_declare_reg t1, t2, f1;
	
	_ff_mult(t1,_ff_third,f[2]);														// t1 = f2/3
	_ff_set(t[0],t1);
	_ff_mult(f1,t1,f[2]);															// f1 = f2^2/3
	_ff_square(t2,t1); ff_mult(t2,t2,t1); _ff_x2(t2); ff_mult(t1,t1,f[1]);						// t2 = f2^3/27, t1 = f1f2/3
	_ff_addto(f[0],t2);  _ff_subfrom(f[0],t1);
	_ff_subfrom(f[1],f1);
	_ff_set_zero(f[2]);
	// 5M+3A
}

// computes cubic resolvent of f=x^4+f2x^2+f1x+f0 and depresses it
void ff_depressed_cubic_resolvent (ff_t t[1], ff_t g[4], ff_t f[5])
{
	_ff_t_declare_reg t0,t1,t2;

	_ff_set_one(g[3]);
	_ff_add(t1,f[2],f[2]);
	_ff_neg(g[2],t1);															// g2 = -2f2
	_ff_add(t1,f[0],f[0]); _ff_x2(t1);
	_ff_square(t2,f[2]);
	_ff_sub(g[1],t2,t1);															// g1 = f2^2-4f0
	_ff_square(g[0],f[1]);														// g0 = f1^2
	ff_depress_cubic (t,g);
	// total 7M+5A
}


// assumes f monic, char > 3, sets t to translation value: depressed cubic is f(x-t)
void ff_depress_quartic (ff_t t[1], ff_t f[5])
{
	_ff_t_declare_reg t0,t1, t2, t3, t4, w1, w2;
	
	_ff_mult(t1,f[3],_ff_half);														// t1 = f3/2
	ff_mult(t0,t1,_ff_half);														// t0 = f3/4
	_ff_set(t[0],t0);			
	_ff_square(t2,t0);  _ff_square(t4,t3);											// t2 = f3^2/16, t4=f3^4/256
	_ff_add(t3,t4,t4);  _ff_addto(t3,t4);												// t3 = 3f3^4/256
	_ff_mult(w1,f[2],t2);															// w1 = f2f3^2/16
	_ff_mult(w2,f[1],t0);															// w2 = f1f3/4
	_ff_subfrom(w1,w2);
	_ff_subfrom(w1,t3);
	_ff_addto(f[0],w1);															// new f0 = f0 + f2f3^2/16 - f1f3/4 - 3f^4/256
	_ff_mult(t3,t0,t2);															// t3 = f3^2/8
	_ff_mult(w1,t3,f[3]);
	_ff_mult(w2,f[2],t1);
	_ff_subfrom(w1,w2);
	_ff_addto(f[1],w1);															// new f1 = f1 + f3^3/8 - f2f3/2
	_ff_add(w1,t3,t3);  _ff_addto(w1,t3);											// w1 = 3f3^2/8
	_ff_subfrom(f[2],w1);														// new f2 = f2-3f3^2/8
	_ff_set_zero(f[3]);
	// 9M+10A
}

// computes the gcd of f and g where f=x^3+ax+b 
void ff_poly_gcd_g1 (ff_t h[4], int *d_h, ff_t f[4], ff_t g[4], int d_g)
{
	_ff_t_declare_reg t1,t2;
	ff_t r[4], s[4];
	int d_r, d_s;
	
	if ( d_g == 3 ) {
		_ff_neg(r[2],g[2]);
		_ff_mult(t1,f[1],g[3]);
		_ff_sub(r[1],t1,g[1]);
		_ff_mult(t1,f[0],g[3]);
		_ff_sub(f[0],t1,g[0]);
		for ( d_r = 2 ; d_r >= 0 && _ff_zero(r[d_r]) ; d_r-- );
	} else {
		ff_poly_copy (r,&d_r,g,d_g);
	}
	if ( d_r < 0 ) { ff_poly_copy (h, d_h, f, 3);  return; }
	if ( ! d_r ) { _ff_set_one(h[0]); *d_h = 0;  return; }
	_ff_neg(s[2],r[d_r-1]);
	_ff_mult(t1,f[1],r[d_r]);
	if ( d_r > 1 ) { _ff_sub(s[1],t1,r[d_r-2]); } else { _ff_set(s[1],t1); }
	_ff_mult(s[0],f[0],r[d_r]);
	for ( d_s = 2 ; d_s >= 0 && _ff_zero(s[d_s]) ; d_s-- );
	if ( d_s == 2 && d_r == 2) {
		_ff_mult(t1,s[1],r[2]);
		_ff_mult(t2,r[1],s[2]);
		_ff_sub(s[1],t1,t2);
		_ff_mult(t1,s[0],r[2]);
		_ff_mult(t2,r[0],s[2]);
		_ff_sub(s[0],t1,t2);
		for ( d_s = 1 ; d_s >= 0 && _ff_zero(s[d_s]) ; d_s-- );
	}
	if ( d_s < 0 ) { ff_poly_copy(h,d_h,r,d_r);  return; }
	if ( ! d_s ) { _ff_set_one(h[0]); *d_h = 0;  return; }
	if ( d_r == 2 && d_s == 1 ) {
		_ff_mult(t1,r[1],s[1]);
		_ff_mult(t2,s[0],r[2]);
		_ff_sub(r[1],t1,t2);
		ff_mult(r[0],r[0],s[1]);
		for ( d_r = 1 ; d_r >= 0 && _ff_zero(r[d_r]) ; d_r-- );
		if ( d_r < 0 ) { ff_poly_copy (h, d_h, s, d_s);  return; }
		if ( ! d_r ) { _ff_set_one(h[0]); *d_h = 0;  return; }
	}
	if ( d_s == 2 && d_r == 1 ) {
		_ff_mult(t1,s[1],r[1]);
		_ff_mult(t2,r[0],s[2]);
		_ff_sub(s[1],t1,t2);
		ff_mult(s[0],s[0],r[1]);
		for ( d_s = 1 ; d_s >= 0 && _ff_zero(s[d_s]) ; d_s-- );
		if ( d_s < 0 ) { ff_poly_copy (h, d_h, r, d_r);  return; }
		if ( ! d_s ) { _ff_set_one(h[0]); *d_h = 0;  return; }
	}
	// we must have d_r==d_s==1 at this point
	_ff_mult(t1,r[0],s[1]);
	_ff_mult(t2,s[0],r[1]);
	_ff_sub(s[0],t2,t1);
	if ( _ff_zero(s[0]) ) { ff_poly_copy(h,d_h,r,d_r);  return; }
	_ff_set_one(h[0]); *d_h = 0;
	// worst case 13M+7A and some copies
	return;	
}

/*
	Standard root finding algorithm specialized to degree 3, note that roots may be null if only the count is required
	Returns the number of roots found.  This function is superseded by ff_poly_roots_d3 below which solves the cubic with radicals.
*/
int ff_old_poly_roots_d3 (ff_t roots[3], ff_t f[4])
{
	_ff_t_declare g[3], h[3], hh[3], a;
	_ff_t_declare_reg t0, t1, t2;
	int i, j, d_g, d_h, d_hh;

#if FF_WORDS > 1
	err_printf ("ff_poly_factors_g1 only implemented for single word fields\n");  exit (0);
#endif

	// this is wasteful if the caller already knows the discriminant!
	_ff_square(t0,f[0]);  _ff_set_ui(t1,27);  _ff_mult(t2,t0,t1);			// t2 = 27f1^2
	_ff_square(t0,f[1]);  _ff_mult(t1,t0,f[1]);  _ff_x2(t1); _ff_x2(t1);		// t1 = 4f0^3
	_ff_addto(t1,t2);
	if ( _ff_zero(t1) ) {											// D=0 => curve has a repeated root
		if ( roots ) {											
			if ( _ff_zero(f[1]) ) {									// if f[1]=0 then f[0]=0 and we have a triple root at 0
				_ff_set_zero(roots[0]);  _ff_set_zero(roots[1]); _ff_set_zero(roots[2]);
			} else {											// otherwise return the distinct root is 3b/a and the double root is (-3b/(2a)
				_ff_invert(t0,f[1]);								// 1/a
				_ff_add(t1,f[0],f[0]);  _ff_addto(t1,f[0]);				// 3b
				_ff_mult(roots[0],t0,t1);							// 3b/a
				_ff_neg(t0,roots[0]);
				_ff_mult(roots[1],t0,_ff_half);						// -3b/(2a)
				_ff_set(roots[2],roots[1]);
			}
		}
		return 3;											
	}
	
	// compute x^p mod f
	ff3_zn_mod (g, _ff_p, f);	// g = x^p mod f
	_ff_dec(g[1]);			// g = x^p-x mod f
	for ( d_g = 2 ; d_g >= 0 && _ff_zero(g[d_g]) ; d_g-- );
	ff_poly_gcd_red (h, &d_h, f, 3, g, d_g);
	switch (d_h) {
	case 0:
		return 0;
	case 1:
		if ( roots ) {
			_ff_invert(t1,h[1]);
			_ff_mult(t2,t1,h[0]);
			_ff_neg(roots[0], t2);
		}
		return 1;
	case 2:
		err_printf ("Impossible, non-singular cubic with 2 roots in ff_old_poly_roots_d3: ");  ff_poly_print (f,3);  ff_poly_print(h,d_h);  exit (0);
	case 3:
		if ( ! roots ) return 3;
		// We know there are three roots, and separate the roots using the standard probabilistic algorithm
		// Pick a random a and compute g=(x+a)^((p-1)/2).  Half of the elements in F_p* will be roots of g, 
		// and with probability 3/4 gcd(f,g) will have degree 1 or 2 allowing us to pick out a root
		// As a minor optimization, we pick a=0 as our first "random" choice, since we can compute g faster in this case.
		for ( i = 0 ; ; i++ ) {
			if ( ! i ) {
				ff3_zn_mod (g, (_ff_p-1)/2, f);
			} else {
				ff_random(&a);
				ff_poly_xpan_mod_d3 (g, a, (_ff_p-1)/2, f);
			}
			_ff_dec (g[0]);		// g = (x+a)^((p-1)/2) mod f
			for ( d_g = 2 ; d_g >= 0 && _ff_zero(g[d_g]) ; d_g-- );
			ff_poly_gcd_red (h, &d_h, f, 3, g, d_g);
			if ( d_h > 0 && d_h < 3 ) break;
		}
		i = ff_poly_roots_d2(roots,h,d_h);
		if ( i != d_h ) { err_printf ("poly_roots_d2 failed in ff_old_poly_roots_d3!");  exit (0); }
		ff_poly_div (h, &d_h, g, &d_g, f, 3, h, d_h);
		if ( d_g != -1 ) { err_printf ("inexact poly division in ff_old_poly_factors_g1!\n");  exit(0); }
		j = ff_poly_roots_d2(roots+i,h,d_h);
		if ( j != d_h ) { err_printf ("poly_roots_d2 failed in ff_old_poly_roots_d3!");  exit (0); }
		if ( i+j != 3 ) { err_printf ("missing root in ff_old_poly_roots_d3!");  exit (0); }
		return 3;
	}
	err_printf ("Unreachable code in ff_old_poly_factors_g1!\n");
	exit (0);
}




/*
	Finds the roots of a monic depreseed cubic x^3+ax+b over F_p.
	If roots is NULL only the # of roots is determined, which can be significantly quicker (e.g. when D is not a QR)
	The value d is optional.  If non-null it must contain sqrt(-D/3) (this is used when computing 3-torsion)

	Assumes _ff_p fits in an unsigned long
*/
int _ff_poly_roots_d3 (ff_t roots[3], ff_t f[4], ff_t *pd)
{
	_ff_t_declare g[3], h[2], w[2], v[2], vk[2], d, t, r, s;
	_ff_t_declare_reg a, t0, t1, t2;
	int i, k, sts;

#if FF_WORDS > 1
	err_printf ("ff_poly_factors_g1 only implemented for single word fields\n");  exit (0);
#endif
	
	// handle p=3 separately
	if ( _ff_p == 3 ) {
		if ( _ff_zero(f[0]) ) { if ( roots ) _ff_set_zero(roots[0]); if ( _ff_one(f[1]) ) return 1; if ( roots ) {  _ff_set_one(roots[1]); _ff_set(roots[2],_ff_negone); } return 3; }
		if ( _ff_zero(f[1]) ) { if ( roots ) { _ff_neg(roots[0],f[0]); _ff_neg(roots[0],f[0]); _ff_neg(roots[0],f[0]); } return 3; }
		if ( _ff_one(f[1]) ) { if ( roots ) _ff_set(roots[0],f[0]);  return 1; } else return 0;
	}
	
	// f=x^3+ax
	if ( _ff_zero(f[0]) ) {
		_ff_neg(t,f[1]);  
		if ( ! roots ) { return (ff_residue(t)?3:1); }
		_ff_set_zero(roots[0]);
		if ( ! ff_sqrt(roots+1,&t) ) return 1;
		_ff_neg(roots[2],roots[1]);
		return 3;
	}
	
	// f = x^3+b
	if ( _ff_zero(f[1]) ) {
		_ff_neg(t,f[0]);
		if ( _ff_p1mod3 ) {
			if ( ! ff_cbrt(&d,&t) ) return 0;
			if ( roots ) { _ff_set(roots[0],d);  _ff_mult(roots[1], roots[0], _ff_cbrt_unity);  _ff_mult(roots[2], roots[1], _ff_cbrt_unity); }
			return 3;
		} else {
			if ( roots ) if ( ! ff_cbrt(roots,&t) ) { printf ("Impossible, ff_cbrt failed on input %ld when p=%ld is 2 mod 3\n", _ff_get_ui(t), _ff_p); exit(0); }
			return 1;
		}
	}

	if ( ! pd ) {
		_ff_square(t0,f[0]);  _ff_set_ui(t1,27);  _ff_mult(t2,t0,t1);			// t2 = 27f1^2
		_ff_square(t0,f[1]);  _ff_mult(t1,t0,f[1]);  _ff_x2(t1); _ff_x2(t1);		// t1 = 4f0^3
		_ff_add(t0,t1,t2);  _ff_mult(t,t0,_ff_third);						// t = -D/3 = (27f1^2+4f0^3)/3
		if ( _ff_zero(t) ) {											// t=0 => D=0 => curve is singular
			if ( roots ) {											// distinct root is 3b/a (we know a!=0), and (-3b/(2a) is a double root).
				_ff_invert(t0,f[1]);									// 1/a
				_ff_add(t1,f[0],f[0]);  _ff_addto(t1,f[0]);					// 3b
				_ff_mult(roots[0],t0,t1);								// 3b/a
				_ff_neg(t0,roots[0]);
				_ff_mult(roots[1],t0,_ff_half);							// -3b/(2a)
				_ff_set(roots[2],roots[1]);
			}
			return 3;											
		}
		sts = ff_sqrt_ext(&d,&t);										// use extended sqrt so we get a solution in F_p^2 if necessary
	} else {
		_ff_set(d,*pd);
	}
	_ff_mult(t2,_ff_half,_ff_third);									// t2=1/6 (used later)
	_ff_mult(t0, _ff_half, f[0]);
	_ff_neg(a,t0);												// a = -f[0]/2 (used later)
	if ( _ff_p1mod3 ) {
		if ( ! sts ) {											// p=1mod3 => -1/3 is a QR => (t not a QR => D not a QR => f has an even # of factors (by Stickleberger), so 1 root
			if ( ! roots ) return 1;								// if all the caller wants is the number of roots, we are done
			_ff_square(t0,t2);									// 1/36
			_ff_mult(s,t,t0);									// s = -D/108 is a not a QR, we will now work in the field F_p^2 = F_p[z]/(z^2-s) to compute (z+a)^1/3
			// compute (z+t1)^((p+1)*m) in F_p^2, this will be in the 3-Sylow subgroup of F_p^2, necessarily in F_p (here p=3^e*m+1, m not divisible by 3)
			// in the process, we compute g=(z+1)^((m-i)/3) where i is congruent to m mod 3 and h=(z+a)^m
			if ( _ff_p3_m1mod3 ) {
				ff_poly_xpan_mod_d2(g,a,(_ff_p3_m-1)/3,&s);			// g = (z+a)^((m-1)/3)
				ff2_square_s(h,g,s);  ff2_mult_s (h,h,g,s);			// h = g^3 = (z+a)^(m-1)
				ff2_mult_zpa_s(h,h,a,s);							// h = g^3*(z+a) = (z+a)^m
			} else {
				ff_poly_xpan_mod_d2(g,a,(_ff_p3_m-2)/3,&s);			// g = (z+a)^((m-2)/3)
				ff2_square_s (h, g, s);  ff2_mult_s (h, h, g, s);		// h = g^3 = (z+a)^(m-2)
				_ff_square(t0,a); _ff_add(w[1],a,a);  _ff_add(w[0],t0,s);	// w = (z+a)^2
				ff2_mult_s (h,h,w,s);							// h = g^3*(z+a)^2 = (z+a)^m
			}
			ff2_norm_s (&t,h,s);								// norm(h)=((z+a)^(m*(p+1)) is in the 3-Sylow of F_p
			if ( ! ff_3Sylow_invcbrt(&r,&t) ) {						// (z+a) is not a cubic residue - this should be impossible
				printf ("(z+%ld)^%ld = ", _ff_get_ui(a), _ff_p3_m); ff_poly_print(h,1);
				printf ("norm(h)=%ld is not a cube in F_%ld\n", _ff_get_ui(t), _ff_p);
				printf ("z+%ld is not a CR in F_p[z]/(z^2-%ld)\n", _ff_get_ui(a), _ff_get_ui(s));
				ff_poly_print(f,3);
				exit (0);
			}
			// we now know (z+a)^-((p+1)m)/3
			ff2_norm_s(&t,g,s);									// t = norm(g) = g^(p+1)
			ff2_exp_ui_s(h,h,(_ff_p+1)/(3*_ff_p3_m),s);				// compute h^(3^(e-1))
			if ( !_ff_p3_m1mod3 ) {
				// We need to construct (z+a)^n where n = (2(p+1)m+1)/3 = 2(p+1)(m-2)/3 + 2(2p+1)/3 + 1 = 2(p+1)(m-2)/3 + 2(2*3^(e-1)*m+1) + 1
				// We then have (z+a)^n = norm(g)^2*((h^(3^(e-1))^2*(z+a))^2*(z+a)   (note that we also need to square r to get (z+a)^(-2(p+1)m/3) here)
				ff_square(t,t);  ff_square(r,r);
				ff2_square_s(h,h,s);  ff2_mult_zpa_s (h,h,a,s);  ff2_square_s(h,h,s);
			}// when m is 1 mod 3:
				// We need to construct (z+a)^n where n = ((p+1)m+1)/3 = (p+1)(m-1)/3 + (p-1)/3 + 1 = (p+1)(m-1)/3 + 3^(e-1)*m + 1
				// We then have (z+a)^n = norm(g)*h^(3^(e-1))*(z+a) 
			ff2_scalar_mult(g,t,h);								
			ff2_mult_zpa_s(g,g,a,s);								// g = (z+a)^n
			ff2_scalar_mult(g,r,g);								// g = (z+a)^(1/3)
			// We know that (g-f1/(3g)) is a root of f.  To get others we multiply by cube roots of unity (which are conveniently in F_p).
			// We use this to find the root in F_p (there must be exactly one, since we know f has 2 factors by Stickleberger)
			// The multiple h of g yielding a root in F_p must have the property that -f1/3=norm(h), which we can test without inverting, and we then have tr(h) as our root
			_ff_mult(t0,f[1],_ff_third);
			ff_negate(t0);
			for ( i = 0 ; i < 3 ; i++ ) {  ff2_norm_s(&t,g,s);  if ( _ff_equal(t,t0) ) break;  ff2_scalar_mult(g,_ff_cbrt_unity,g); }
			if ( i == 3 ) {
				printf ("g=%ldz+%ld is a cbrt of (z+%ld) in F_p[z]/(z^2-%ld), but couldn't find an F_%ld root of 2-factor f = ", _ff_get_ui(g[1]), _ff_get_ui(g[0]),_ff_get_ui(a), _ff_get_ui(s), _ff_p);
				ff_poly_print(f,3);
				printf ("t0=%ld, cube root of unity is %ld\n", _ff_get_ui(t0), _ff_get_ui(_ff_cbrt_unity));
				exit (0);
			}
			ff2_trace(roots,g);
			return 1;
		} else {												// in this case have a sqrt of t and we know that f has 1 or 3 factors
			ff_mult(d,d,t2);									// d = sqrt(t)/6 = sqrt(f1^3/27+f0^2/4) (this is the square root in Cardano's formula for depressed cubics)
			_ff_add(t,a,d);										// t = -f0/2 + d (first cube to check)
			if ( ! ff_cbrt(&r,&t) ) return 0;
			_ff_sub(t,a,d);										// t = -f0/2 - d (second cube to check - we really ought to be able to compute this from the first one)
			if ( ! ff_cbrt(&s,&t) ) { printf("Impossible cube root failure in ff_poly_roots_d3 for p=%ld, f = ", _ff_p); ff_poly_print(f,3);  exit(0); }
			if ( ! roots ) return 3;								// we know the roots exist but we don't need to find them
			// We have to cycle through the choices of one of our cube roots to verify that we have the correct s and t, unforutnately this means we may need to evaluate f twice
			_ff_add(t,r,s);
			_ff_square(t0,t);  _ff_addto(t0,f[1]);  _ff_multadd(t1,t0,t,f[0]);
			if ( ! _ff_zero(t1) ) {
				ff_mult(s,s,_ff_cbrt_unity);						// try second s
				_ff_add(t,r,s);
				_ff_square(t0,t);  _ff_addto(t0,f[1]);  _ff_multadd(t1,t0,t,f[0]);
				if ( ! _ff_zero(t1) ) ff_mult(s,s,_ff_cbrt_unity);			// must be the third s
			}
			_ff_add(roots[0],r,s);								// third choice must work since first two didn't
			_ff_square(t0,_ff_cbrt_unity);							// could cache this in ff.c
			_ff_mult(t1,r,_ff_cbrt_unity);
			_ff_mult(t2,s,t0);
			_ff_add(roots[1],t1,t2);								// second root is wr+w^2s where w is a primitive cbrt of unity
			_ff_mult(t1,r,t0);
			_ff_mult(t2,s,_ff_cbrt_unity);
			_ff_add(roots[2],t1,t2);								// third root is wr^2+ws
			return 3;
		}
	} else {													// p=2mod3
		if ( sts ) {											// p=2mod3 => -1/3 not a QR => (-D/3 a QR => D not a QR => f has an even # of factors (by Stickleberger))
			if ( ! roots ) return 1;								// if we only need the # of factors, we are done.
			ff_mult(d,d,t2);									// d = sqrt(t)/6 = sqrt(f1^3/27+f0^2/4) (this is the square root in Cardano's formula for depressed cubics)
			_ff_add(t,a,d);										// t = -f0/2 + d (first cube)
			ff_cbrt(&r,&t);										// p is 2 mod 3 so we know there is exactly 1 cube root
			_ff_sub(t,a,d);										// t = -f0/2 - d (second cube)
			ff_cbrt(&s,&t);
			_ff_add(t,r,s);										// there was no choice of cuberoots (and changing sign of d doesn't change r+s=t), so this has got to be the root
			_ff_set(roots[0],t);
			return 1;
		} else {												// in this case t = -D/3 is not a QR, D=d^2 is a QR, and we know that f has 1 or 3 factors
			// for p=2mod3, handle roots=NULL separately, since we can be much more efficient
			if ( ! roots ) {
				_ff_square(t0,t2);								// 1/36
				_ff_mult(s,t,t0);								// s = -D/108 is a not a QR, we will now work in the field F_p^2 = F_p[z]/(z^2-s)
				ff_poly_xpan_mod_d2(g,a,(_ff_p+1)/3,&s);			// (z+a)^((p+1)/3) is in F_p iff (z+a)^((p+1)(p-1)/3)=1 iff (z+a) is a cubic residue
				return ( _ff_zero(g[1]) ? 3 : 0 );
			}

			// in this case we will be slightly less efficient and work in the standard basis for F_p^2=F_p[z](z^2-s) because we
			// want don't want to have to compute a new generator of the Sylow 3-subgroup in F_p^2 every time 
			_ff_set(v[0],a);  _ff_mult(v[1],d,t2);						// v is the element of F_p^2 (in the standard basis) whose cuberoot we need
			// this code is slightly inefficient, but we are trying to compute h=v^((p+1)/3) as quickly as possible (if h is not in F_p, f is irreducible and we are done)
			// while also saving values we will need to compute a cube root of v when f splits
			if ( _ff_p3_m1mod3 ){ k=1; ff2_set(vk,v); }
			else { k=2; ff2_square(vk,v); }						// vk = v^k
			ff2_exp_ui(g,v,(_ff_p3_m-k)/3);							// g = v^((m-k)/3)
			ff2_square(h,g);  ff2_mult(h,h,g);						// h = v^(m-k)
			ff2_mult(h,h,vk);									// h = v^m
			if ( _ff_p3_e > 1 ) {
				ff2_exp_ui(w,h,(_ff_p+1)/(3*_ff_p3_m)-1);				// w = v^((3^(e-1)-1)m) = v^(((p+1)-3m)/3)
				_ff_mult(t1,w[1],h[0]);  _ff_mult(t2,w[0],h[1]);			// compute z coefficient of w*h
				_ff_add(t0,t1,t2);
				if ( ! _ff_zero(t0) ) return 0;						// if w*h=v^((p+1)/3) is not in F_p, then v^((p^2-1)/3)=h^(p-1) is not 1 and v is not a cube, and o.w. it is
			} else {
				if ( ! _ff_zero(h[1]) ) return 0;						// here h=v^m=v^((p+1)/3) since e=1
			}
			// We now know v is a cube and there are three roots
			if ( ! roots ) return 3;
			ff2_norm(&t,g);									// t = N(g)=g^(p+1)=v^((p+1)(m-k)/3) where k=m mod 3
			ff2_scalar_mult(g,t,g);								// g = v^((p+2)(m-k)/3)
			// When e==1, life is pretty simple: we know v is a cube, and there is only one cube in the Sylow 3-subgroup, name 1, and its inverse cube root is 1
			// It follows that we can assume v^((p-1)m/3)=v^((p-2)m/3)=1.
			if ( _ff_p3_e == 1 ) {
				if ( k==2 ) {
					ff2_mult(g,g,h);							// g = v^((p+2)(m-2)/3+m)=v^((pm+2m-2p-4+3m)/3)=v^((pm+2m-2(3m-1)-4+3m)/3)=v^(((p-1)m-2)/3)
				} else {										// note that e=1 and k=1 implies g=v^((p+1)(m-1)/3) = v^(((p-1)m-1)/3)
					ff2_square(g,g);							// g=v^((2(p-1)m-2)/3)
				}
				ff2_mult(g,g,v);								// g = v^(((3-k)(p-1)m+1)/3) = v^(1/3}
				ff2_trace(roots,g);								// Every choice of cube root of v corresponds to a root of f, so we can just take the trace
				ff2_setup_cbrt();								// we still need to get the cube root of unity (this is annoying)
				ff2_mult(g,g,_ff2_cbrt_unity);
				ff2_trace(roots+1,g);
				ff2_mult(g,g,_ff2_cbrt_unity);
				ff2_trace(roots+2,g);
				return 3;
			}
			// now we have e>1
			if ( k==2 ) { ff2_square(w,w);  ff2_mult(w,w,h); }			// w = v^((k*3^(e-1)-1)m)
			ff2_mult(g,g,w);									// g = v^(((p-1)m-k)/3)  (note that (p+2)(m-k)/3 + (k*3^(e-1)-1)m = (pm+2m-kp-2k+k3^em-3m)/3
															//                                                                                                     = (pm-m-kp-2k+k(p+1))/3=((p-1)m-k)/3
			ff2_square(h,g);  ff2_mult(h,h,g);						// h = v^((p-1)m-k)
			ff2_mult(h,h,vk);									// h = v^((p-1)m) is in the 3-Sylow subgroup, and it must be a cubic residue
			if ( ! ff2_3Sylow_invcbrt(h,h) ) {						// h = v^(-(p-1)m)/3)  if e=1 
				printf ("Impossible situation in ff_poly_factors_g1 with p=%ld,m=%ld:  v^((p-1)m)=(%ldz+%ld) is not a cube in 3-Sylow of F_p^2=F_p[z]/(z^2-%ld) ",
					   _ff_p, _ff_p3_m, _ff_get_ui(a), _ff_get_ui(h[1]), _ff_get_ui(h[0]), _ff_get_ui(s));
				ff_poly_print(f,3);
				exit (0);
			}
			if ( k==1 ) {
				ff2_square(h,h);								// h = v^(-2(p-1)m/3)
				ff2_square(g,g);								// g = v^(2((p-1)m-1)/3)
				ff2_mult(g,g,v);								// g = v^(2((p-1)m-1)/3+1) = v^((2(p-1)m+1)/3
			} else {
				ff2_mult(g,g,v);								// g = v^((p-1)m-2)/3+1 = v^((p-1)m+1/3
			}
			ff2_mult(g,g,h);									// g = v^(1/3)
			ff2_trace(roots,g);									// Every choice of cube root of v corresponds to a root of f, so we can just take the trace of each
			ff2_mult(g,g,_ff2_cbrt_unity);
			ff2_trace(roots+1,g);
			ff2_mult(g,g,_ff2_cbrt_unity);
			ff2_trace(roots+2,g);
			return 3;
		}
	}
	puts ("Unreachable code in ff_poly_factors_g1!");
	exit (0);
}

// Compute roots of poly of degree 2, or less
int ff_poly_roots_d2 (ff_t r[2], ff_t f[3], int d_f)
{
	ff_t t1,t2,D;
	
	switch (d_f) {
	case 0:  return 0;
	case 1: 
		if ( _ff_one(f[1]) ) { _ff_set_one(t1); } else { _ff_invert(t1,f[1]); }
		ff_negate(t1); _ff_mult(r[0],t1,f[0]);  return 1;
	case 2:
		_ff_square(t1,f[1]);
		_ff_mult(t2,f[2],f[0]);
		_ff_x2(t2);  _ff_x2(t2);
		_ff_subfrom(t1,t2);
		if ( ! ff_sqrt(&D,&t1) ) return 0;
		if ( _ff_one(f[2]) ) {					
			_ff_set(t1,_ff_half);
		} else {
			_ff_add(t2,f[2],f[2]);
			_ff_invert(t1,t2);
		}
		_ff_sub(t2,D,f[1]);
		_ff_mult(r[0],t1,t2);
		ff_negate(D);
		_ff_sub(t2,D,f[1]);
		_ff_mult(r[1],t1,t2);
		return 2;
	}
}


/*
	ff_poly_g1_3tor computes the size of the 3-torsion subgroup of y^2=f(x) = x^3+ax+b.
	If the 3-torsion is trivial, it will instead attempt to determine the 3-torsion subgroup of the twist, and if it is non-trivial
	will return a negative value whose absolute value is the size of the 3-torsion subgroup in the twist.

        This means that a return value of 1 indicates that neither the curve nor its twist have non-trivial 3-torsion.
	In this situation, if p=1mod3, then a_p=0mod3 and the group order must be congruent to 2 mod 3.

	If the curve is singular (discriminant zero), 0 is returned.
*/
int ff_poly_g1_3tor (ff_t f[4])
{
	ff_t div[5],r[4],g[4],t0,t1, t2, d;
	int j, k, tor;
	
#if FF_WORDS > 1
	err_printf ("ff_poly_g1_3tor only supports single word finite fields\n"); exit (0);
#endif

	if ( _ff_p == 3 ) { if ( _ff_zero(f[1]) ) return 0; else return 1; }					// no curve y^2=x^3+ax+b over F_3 has 3 torsion
	
	_ff_square(t1,f[0]); _ff_set_ui(t2,27); _ff_mult(t0,t1,t2);							// t0 = 27f0^2
	_ff_square(t2,f[1]); _ff_mult(t1,t2,f[1]); _ff_x2(t1); _ff_x2(t1);						// t1 = 4f1^3, t2 = f1^2
	_ff_add(d, t0, t1);														// d = 4f1^3+27f0^2 = -discriminant of f
	if ( _ff_zero(d) ) return 0;
	_ff_add(t0,_ff_third, _ff_third); _ff_x2(t0); _ff_square(t1,t0);						// t0 = (4/3)^2 = 16/9
	ff_mult(d,d,t0);														// d^2 = 256/81*D_f^2 = -D/3 where D is discriminant of the 3-div poly

	ff_negate(t2);
	ff_mult(div[0],_ff_third,t2);	
	_ff_add(div[1],f[0],f[0]); _ff_x2(div[1]);
	_ff_add(div[2],f[1],f[1]);
	_ff_set_zero(div[3]);
	_ff_set_one(div[4]);
	// 3-division poly is x^4 + 2f1x^2 +4f0x - f1^2/3

	// handle p=2mod3 ourselves since we can easily save some time.  note that in this case, either both the curve and its twist have non-trivial 3-torsion, or neither do.
	// we also know that the 3-rank is 1 and there can be at most one root of the 3-division poly (2 pts of order 3).
	if ( ! _ff_p1mod3 ) {
		if ( _ff_zero(div[1]) ) {												// biquadratic case is handled quickly by ff_poly_roots
			k = _ff_poly_roots_d4(r,div,&d);
			return ( k ? 3 : 1);
		} else {
			ff_depressed_cubic_resolvent(&t0,g,div);
			k = _ff_poly_roots_d3(r,g,&d);
			if ( k != 1 ) { printf ("p=%ld is 2mod3 but cubic resolvent of 3div poly has k=%d!=1 roots!\n  f=", _ff_p, k); ff_poly_print(f,3); exit(0); }
		}
		_ff_sub(t1,t0,r[0]);
		return ( ff_residue(t1) ? 3 : 1);
	}
	
	// Now handle p=1mod3 by computing the roots of the 3-division polynomial (we actually only need to know the number of roots and the value of one of them)
	k = _ff_poly_roots_d4(r,div,&d);
	if ( ! k ) return 1;
	_ff_square(t1,r[0]);  _ff_addto(t1,f[1]); _ff_mult(t2,t1,r[0]); _ff_addto(t2,f[0]);		// t2 = f(r[0])
	if ( _ff_zero(t2) ) { printf ("impossible - non-id point (%lu,0) with 2-tor and 3-tor!\n", _ff_get_ui(r[j]));  ff_poly_print(f,3);  ff_poly_print(div,4); exit(0); }
	
	// We now rely on the fact that the 3-division poly of the curve and its twist are intimately related.  They have the same factorization pattern,
	// and either all the roots of the 3-division poly correspond to points on the curve, or none of them do, so we only need to check one.
	// In the latter case, it must be that all the roots of the 3-division poly of the twist correspond to points on the twist.
	if ( ff_residue(t2) ) return ( k==1?3:9); else return (k==1?-3:-9);
}

#if SMALLJAC_GENUS==2
int ff_poly_g2_3tor(ff_t f[6])
{
	ff_t g[POLY_G2TOR3_DEGREE+1], h[POLY_G2TOR3_DEGREE], r[POLY_G2TOR3_DEGREE];
	int k, d_h, d_r;
	
	ff_poly_g2tor3_modpoly(g,f);
	
	// first make sure modular poly is square free
	ff_poly_derivative (h,&d_h,g,POLY_G2TOR3_DEGREE);
	ff_poly_gcd (r,&d_r,g,POLY_G2TOR3_DEGREE,h,d_h);
	if ( d_r > 0 ) return 0;
	
	k = ff_poly_roots(g,40);
	if ( ! k ) return 1;
	if ( k ==2 ||k==5 ||k==8 ) return 3;
	return 0;
}
#endif

/*
	Given a point (x,y) on y^2=f(x)=x^3+f1x+f0, replaces (x,y) with a point (u,v) s.t. (u,v) composed with itself yields (x,y) or its inverse (we don't distinguish these cases!)
	Returns 1 for success, 0 if no such point exists.
*/
int ff_poly_g1_halve (ff_t *x, ff_t *y, ff_t f[4])
{
	ff_t g[5], r[4];
	register ff_t t0, t1, t2;
	int i, k;
	
	// construct g(z) = z^4 - 2(f1+3x^2)z^2 - 8y^2z + f1^2 - 3x(y^2+f1x + 3f0).  If (u,v)*(u,v)=(x,y) then g(u-x)=0 (note translation by x to kill degree 3 coeff of g).
	_ff_set_one(g[4]);  _ff_set_zero(g[3]);
	_ff_square(t0,*x);  _ff_add(t1,t0,t0); _ff_addto(t1,t0);  _ff_addto(t1,f[1]); _ff_x2(t1); _ff_neg(g[2],t1);			// g[2] = -2(f1+3x^2)
	_ff_add(t0,*y,*y); _ff_square(t1,t0); _ff_x2(t1); _ff_neg(g[1],t1);										// g[1] = -8y^2
	_ff_square(t0,*y); _ff_mult(t1,f[1],*x); _ff_addto(t0,t1); _ff_add(t1,f[0],f[0]); _ff_addto(t1,f[0]); _ff_addto(t0,t1);	// t0 = y^2 + f1x + 3f0
	_ff_add(t1,*x,*x); _ff_addto(t1,*x); _ff_mult(t2,t0,t1); _ff_square(t1,f[1]);  _ff_sub(g[0],t1,t2);				// g[0] = f1^2 - 3x(y^2+f1x+3f0)
//printf ("(%d,%d) havling poly: ", _ff_get_ui(*x), _ff_get_ui(*y)); ff_poly_print(g,4);
	
	k  = ff_poly_roots_d4 (r, g);
//printf ("%d roots\n", k);
	for ( i = 0 ; i < k ; i++ ) {
		_ff_add(t2,r[i],*x);
		_ff_square(t0,t2);  _ff_addto(t0,f[1]);  _ff_mult(t1,t0,t2); _ff_add(g[1],t1,f[0]);						// g[1] = f(r[i]+x)
		if ( ff_sqrt(g,g+1) ) break;
	}
	if ( i == k ) return 0;
//printf ("found (%d,%d)^2 = (%d,%d)\n", _ff_get_ui(t2), _ff_get_ui(g[0]), _ff_get_ui(*x), _ff_get_ui(*y));
	_ff_set(*x,t2); _ff_set(*y,*g);
	return 1;
}


/*
	Returns the size of the 2-Sylow subgroup of the elliptic curve y^2=f(x)=x^3+f1x+f0, provided it is cyclic.
	Otherwise the return value is -1, which indicates that Z/2Z x Z/2Z is a subgroup (and the group order is divisible by 4)
        (we could compute the entire 2-Sylow in this case, but it would be much more time consuming)
*/
int ff_poly_g1_2Sylow (ff_t f[4])
{
	ff_t r[3];
	ff_t x, y;
	int n;
	
	n = ff_poly_roots_d3(r,f);
	if ( ! n ) return 1;
	if ( n > 1) return 0;
	
	_ff_set(x,r[0]);  _ff_set_zero(y);
	for ( n = 2 ; ff_poly_g1_halve (&x, &y, f) ; n<<= 1 );
	return n;
}


/*
	Computes the order and rank of the 4-torsion subgroup of the elliptic curve y^2=f(x)=x^3+f1x+f0. 
	Set o to the order and returns d to indicate the group Z/dZ x Z/eZ where d*e = o, d divides e (and may be 1)

	If flag8 is set, also check for Z/8Z (for a random curve with full Galois image in GL(2,Z/8Z), these occurs with probability 1/8).
	We could also check for Z/2ZxZ/8Z, or even Z/4ZxZ/8Z or Z/8ZxZ/8Z (these occur with prob. 3/64, 3/512 and 1/1536 resp);
	
	The latter two are rare enough not to be worth the trouble, but Z/2xZ/8Z might be worth doing.  Note, however, that Z/2xZ/4Z
	contains 4 points of order 4 (two pairs of inverses) and we would need to compute two of these.
*/
int ff_poly_g1_4tor (int *o, ff_t f[4], int flag8)
{
	ff_t r[3];
	ff_t u, v, v2;
	register ff_t t0,t1,t2;
	register int i, n;
	
	n = ff_poly_roots_d3(r,f);
	for ( i = 0 ; i < n ; i++ ) {
		_ff_square(t0,r[i]); 	_ff_add(t1,t0,t0); _ff_addto(t1,t0); _ff_add(u,t1,f[1]);			// u = 3x^2+a where (x,0) is a 2-torsion point (x=r[i] was a root of f)
		if ( ff_sqrt(&v,&u) ) {													// v is a root of the translated halving poly
			_ff_add(t2,r[i],v);													// t2=x+v is a root of the halving poly
			_ff_square(t0,t2);  _ff_addto(t0,f[1]);  _ff_mult(t1,t0,t2); _ff_add(u,t1,f[0]);		// u = f(t2)
			if ( ff_sqrt(&v,&u) ) break;											// if successful, (t2,v) is a point of order 4
			_ff_sub(t2,r[i],v);													// t2=x-u is also a root of the halving poly
			_ff_square(t0,t2);  _ff_addto(t0,f[1]);  _ff_mult(t1,t0,t2); _ff_add(u,t1,f[0]);		// u = f(t2)
			if ( ff_sqrt(&v,&u) ) break;											// if successful, (t2,v) is a point of order 4
		}
	}
	if ( i == n ) { *o = n+1;  return (n==3?2:1); }									// no pts of order 4, so 4-torsion subgroup = 2-torsion subgroup
//printf ("%ld: found point(%ld,%ld) with order 4 (i=%d,n=%d) on curve y^2 = ", _ff_p, _ff_get_ui(t2), _ff_get_ui(v), i, n); ff_poly_print(f,3);
	if ( n==1 ) {																// rank 1 case, with a point of order 4
		if ( ! flag8 ) { *o = 4;  return 1; }
//printf("%ld: flag8 set (Z/4Z), attempting to halve order 4 point (%ld,%ld) on curve y^2 = ", _ff_p, _ff_get_ui(t2), _ff_get_ui(v)); ff_poly_print(f,3);
		_ff_set(u,t2);
		*o = ( ff_poly_g1_halve(&u,&v,f) ? 8 : 4 );									// check for a point of order 8 (this takes us outside the 4-torsion subgroup)
//if ( *o == 8 ) puts ("Succeeded"); else puts ("Failed");
		return 1;
	}
	if ( i > 0 ) { *o = 8; return 2; }	
		
	// We are here if the 2-rank is 2 and we succesfully halved the first point of order 2 that we tried.  We need to check one more to distinguish Z/2Z x Z/4Z from Z/4Z x Z/4Z.
	// Don't change v or t2 from above, we may need them below
	i = 1;
	_ff_square(t0,r[i]); 	_ff_add(t1,t0,t0); _ff_addto(t1,t0); _ff_add(u,t1,f[1]);			// u = 3x^2+a where (x,0) is a 2-torsion point (x=r[0] was a root of f)
	if ( ff_sqrt(&v2,&u) ) {													// v2 is a root of the translated halving poly
		_ff_add(t1,r[i],v2);													// t1is a root of the halving poly
		_ff_square(t0,t1);  _ff_addto(t0,f[1]);  ff_mult(t1,t1,t0); _ff_addto(t1,f[0]);		// t1 = f(t1)
		if ( ff_residue(t1) ) { *o = 16; return 4; }
	}
	*o = 8;
	return 2;
}

/*
	Computes the roots of f(x)=x^4+ax^2+bx+c, solving by radicals.
	Handles singular f, will return repeated roots correctly.

	Currently r is required.  We could easily extend to make it optional and only return a count.
	The parameter pd is optional.  If specified it points to sqrt(-D/3) (used for 3-torsion).
*/
int _ff_poly_roots_d4 (ff_t r[4], ff_t f[5], ff_t *pd)
{
	ff_t g[4], h[3], s[3], x, y, t;
	register ff_t t0,t1,t2,w1,w2;
	register int i,k,res;

	if ( _ff_zero(f[1]) ) {
		// f = x^4+f2x^2+f0 = g(x^2) where g is monic quadratic
		_ff_set_one(g[2]);  _ff_set(g[1],f[2]);  _ff_set(g[0],f[0]);
		k = ff_poly_roots_d2(s,g,2);
		if ( ! k ) return 0;
		if ( ff_sqrt(r,s) ) {
			_ff_neg(r[1],r[0]);
			// could optimize for double root here, but why bother
			if ( ff_sqrt(r+2,s+1) ) { _ff_neg(r[3],r[2]); return 4; }
			return 2;
		}
		if ( ff_sqrt(r,s+1) ) { _ff_neg(r[1],r[0]); return 2; }
		return 0;
	}
	ff_depressed_cubic_resolvent(&t,g,f);
	k = _ff_poly_roots_d3(r,g,pd);				// put the roots of g into r for now, they will get replaced by roots of f below.
	if ( k ) {
		res = 0;
		for ( i = 0 ; i < k ; i++ ) {
//			ff_poly_eval (&y, g, 3, r+i);
//			if ( ! _ff_zero(y) ) { printf ("(%d roots) %ld is not a root of ", k, _ff_get_ui(r[i]));  ff_poly_print(g,3); exit(0); }
			_ff_sub(x,r[i],t);
			ff_negate(x);
			if ( ! ff_sqrt(s+i,&x) ) break;
			if ( ++res == 2 ) break;
		}
		if ( res==2 ) {
			_ff_mult(t1,s[0],s[1]);
			ff_negate(t1);
			_ff_invert(t0,t1);
			_ff_mult(s[2],t0,f[1]);
			_ff_square(t1,s[2]);
			_ff_sub(x,r[2],t);
			ff_negate(x);
			if ( ! _ff_equal(t1,x) ) { printf("%ld: division failed s[0]=%d, s[1]=%d, d=%d\n", _ff_p, _ff_get_ui(s[0]), _ff_get_ui(s[1]), _ff_get_ui(t0));  ff_poly_print(f,4); ff_poly_print(g,3); exit(0); }
			_ff_add(t1,s[0],s[1]);  _ff_addto(t1,s[2]); ff_mult(t1,t1,_ff_half);
			_ff_set(r[0],t1);
			ff_negate(s[0]); ff_negate(s[2]);
			_ff_add(t1,s[0],s[1]);  _ff_addto(t1,s[2]); ff_mult(t1,t1,_ff_half);
			_ff_set(r[1],t1);
			ff_negate(s[0]); ff_negate(s[1]);
			_ff_add(t1,s[0],s[1]);  _ff_addto(t1,s[2]); ff_mult(t1,t1,_ff_half);
			_ff_set(r[2],t1);
			ff_negate(s[0]); ff_negate(s[2]);
			_ff_add(t1,s[0],s[1]);  _ff_addto(t1,s[2]); ff_mult(t1,t1,_ff_half);
			_ff_set(r[3],t1);
			return 4;
		} else if ( k==1 && res==1 ) {
			_ff_set_one(h[2]);
			_ff_set(h[1],s[0]);
			_ff_square(t1,s[0]);						// t1 = s^2
			_ff_add(w1,t1,f[2]);
			ff_mult(w1,w1,_ff_half);					// w1 = (s^2+f2)/2
			_ff_square(t0,f[2]);
			_ff_add(t2,f[0], f[0]); _ff_x2(t2);			// t2 = 4f0
			_ff_subfrom(t0,t2);						// t0 = f2^2-4f0
			_ff_add(t2,f[2],f[2]);					// t2 = 2f2
			_ff_add(w2,t1,t2);						// w2 = s^2+2f2
			ff_mult(w2,w2,t1);						// w2 = s^4+2f2s^2
			_ff_addto(w2,t0);
			ff_mult(w2,w2,s[0]);						// w2 = s(s^4+2f2s^2+f2^2-4f0)				
			_ff_invert(t0,f[1]);
			ff_mult(t1,_ff_half,t0);
			ff_mult(w2,w2,t1);
			_ff_sub(h[0],w1,w2);
			if ( ff_poly_roots_d2(r,h,2) != 2 ) {
				ff_negate(h[1]);
				_ff_add(h[0],w1,w2);
				if ( ff_poly_roots_d2(r,h,2) != 2 ) { printf("%ld: fail 2-1-1 split with s=%ld\n", _ff_p, _ff_get_ui(s[0]));  ff_poly_print(f,4);  ff_poly_print(g,3);  ff_poly_print(h,2); exit(0); }
			}
			return 2;
		} else {
			return 0;
		}
	} else {
		// this is the annoying case, we know we have a 1,3 split, but getting the value of the root is expensive
		// we need sqrt(-z) in F_p[z]/h(z) where h(z) is the undepressed resolvent
		// so we want sqrt(t-z) in F_p[z]/g(z) which is sqrt(z+t) in F_p[z]/k(z) where k(z)=-g(-z)=z^3+g1z-g0z.  In terms of z^3-rz-s, r=-g1 and s=g0
		_ff_neg(t1,g[1]);
		if ( ff3_trsqrt_zpa_mod_rs(&x,t,t1,g[0]) ) res++; else { printf ("%ld: impossible failure of ff3_trsqrt_zpa_mod_rs in ff_poly_roots_d4!", _ff_p); ff_poly_print(f,4); ff_poly_print(g,3); exit(0); }
		ff_mult(x,x,_ff_half);
		ff_poly_eval(&y,f,4,&x);
		if ( ! _ff_zero(y) ) ff_negate(x);
		_ff_set(r[0],x);
		return 1;
	}
	puts ("Unreachable code in ff_poly_factors_g1!");
	exit (0);
}

/*
	Hard-wired code for computing the roots of f(x)=x^4+ax^2+bx+c.  This is used
	primarily for testing for 3-torsion by computing roots of the 3-division poly (made monic).
	If f is singular, it only returns distinct roots.

	It is generally very fast, although when f splits completely more work is required.
*/
int ff_old_poly_roots_d4 (ff_t r[4], ff_t f[5])
{
	_ff_t_declare g[5], h[5], A[5], t, a;
	int i,j,k,d_g, d_h, d_A;

#if FF_WORDS > 1
	err_printf ("ff_poly_roots_d4 only implemented for single word fields\n");  exit (0);
#endif
	ff_poly_xn_mod_d4 (g, _ff_p, f);		// compute x^p mod f
	_ff_set_one(t);
	_ff_subfrom (g[1],t);				// x^p-x mod f
	for ( d_g = 3 ; d_g >= 0 && _ff_zero(g[d_g]) ; d_g-- );
	ff_poly_gcd (h, &d_h, g, d_g, f, 4);
	switch (d_h) {
	case 0: return 0;
	case 1:
	case 2: return ff_poly_roots_d2 (r, h, d_h);
	case 3:
	case 4:
		ff_poly_monic (A, &d_A, h, d_h);
		i = j = k = 0;
		_ff_set_zero(t);
		do {
			// if A was just made cubic, it may need to be translated to make x^2 coeff zero
			// if this is done, we remember the translation so we can shift roots back at the end
			if ( d_A == 3 && ! _ff_zero(A[2]) ) { j = i; ff_depress_cubic(&t, A); }
			if ( k ) ff_random(&a);		// use a=0 first time through
			if ( d_A == 4 ) {
				if ( !k ) {
					ff_poly_xn_mod_d4 (g, (_ff_p-1)/2, A);
				} else {
					ff_poly_xpan_mod_d4 (g, a, (_ff_p-1)/2, A);
				}
				_ff_dec (g[0]);
				for ( d_g = 3 ; d_g >= 0 && _ff_zero(g[d_g]) ; d_g-- );
			} else {
				if ( !k ) {
					ff3_zn_mod (g, (_ff_p-1)/2, A);
				} else {
					ff_poly_xpan_mod_d3 (g, a, (_ff_p-1)/2, A);
				}
				_ff_dec (g[0]);
				for ( d_g = 2 ; d_g >= 0 && _ff_zero(g[d_g]) ; d_g-- );
			}
			ff_poly_gcd (h, &d_h, g, d_g, A, d_A);
			if ( d_h > 0 && d_h < d_A ) {
				ff_poly_div (A, &d_A, g, &d_g, A, d_A, h, d_h);
				if ( d_g != -1 ) { err_printf ("inexact poly division in ff_poly_roots_d4!\n");  exit(0); }
				if ( d_h < 3 ) {
					i += ff_poly_roots_d2 (r+i,h,d_h);
					ff_poly_copy (h,&d_h,A,d_A);
					ff_poly_monic (A,&d_A,h,d_h);
				} else {
					i += ff_poly_roots_d2 (r+i,A,d_A);
					ff_poly_monic (A,&d_A,h,d_h);
				}
			}
			k++;
		} while ( d_A > 2 );
		i += ff_poly_roots_d2 (r+i,A,d_A);
		// shift any roots found after translation occurred (if it did)
		if ( ! _ff_zero(t) ) while ( j < i ) { _ff_subfrom(r[j],t); j++; }
		return i;
	}
	return d_h;
}


void ff_poly_xn_mod_d4 (ff_t g[4], unsigned long n, ff_t f[5])
{
	_ff_t_declare_reg t1,t2,t3,nf0,f02,f12,f22,w1,w2;
	register unsigned long j, m;
	register int i;

	for ( i = 0 ; i < 4 ; i++ ) _ff_set_zero(g[i]);
	if ( !n ) { _ff_set_one(g[0]);  return; }
	i = _asm_highbit(n);
	if ( i&1 ) i--;
	m = 3UL<<i;
	j = (n&m)>>i;
	_ff_set_one(g[j]);
	m >>= 2;  i -= 2;
	if ( m > 1 ) {
		_ff_neg(nf0,f[0]);
		// these won't get used if no digits of n are 3, we could check for this...
		_ff_mult(f02,f[0],f[2]);
		_ff_mult(f12,f[1],f[2]);
		_ff_square(f22,f[2]);
	}
	while (m) {
		ff_poly_square_mod_d4 (g,g,f);
		ff_poly_square_mod_d4 (g,g,f);
		j = ((n&m)>>i);
		switch (j) {
		case 1:	// 3M+2A
			_ff_set(t3,g[3]);
			_ff_set(g[3],g[2]);
			_ff_mult(w2,f[2],t3);
			_ff_sub(g[2],g[1],w2);
			_ff_mult(w2,f[1],t3);
			_ff_sub(g[1],g[0],w2);
			_ff_mult(g[0],t3,nf0);
			break;
		case 2: // 6M+3A
			_ff_set(t3,g[3]);
			_ff_set(t2,g[2]);
			_ff_mult(w1,f[2],t3);
			_ff_sub(g[3],g[1],w1);
			_ff_mult(w1,f[2],t2);
			_ff_mult(w2,f[1],t3);
			_ff_addto(w2,w1);
			_ff_sub(g[2],g[0],w2);
			_ff_mult(w1,t3,nf0);
			_ff_mult(w2,f[1],t2);
			_ff_sub(g[1],w1,w2);
			_ff_mult(g[0],t2,nf0);
			break;
		case 3:	// 9M+8A
			_ff_set(t3,g[3]);
			_ff_set(t2,g[2]);
			_ff_set(t1,g[1]);
			_ff_mult(w1,t2,f[2]);
			_ff_mult(w2,t3,f[1]);
			_ff_addto(w1,w2);
			_ff_sub(g[3],g[0],w1);
			_ff_sub(w2,f22,f[0]);
			_ff_mult(w1,t3,w2);
			_ff_mult(w2,t1,f[2]);
			_ff_subfrom(w1,w2);
			_ff_mult(w2,t2,f[1]);
			_ff_sub(g[2],w1,w2);
			_ff_mult(w1,t1,f[1]);
			_ff_mult(w2,t2,f[0]);
			_ff_addto(w1,w2);
			_ff_mult(w2,t3,f12);
			_ff_sub(g[1],w2,w1);
			_ff_mult(w1,t3,f02);
			_ff_mult(w2,t1,f[0]);
			_ff_sub(g[0],w1,w2);			
		}
		m >>= 2;  i -= 2;
	}
}

// Computes (x+a)^n modulo x^2-f0 (note the sign).  Assumes n < 2^63
void ff_poly_xpan_mod_d2 (ff_t g[2], ff_t a, unsigned long n, ff_t f[1])
{
	_ff_t_declare_reg t0,t1,t2,t3,s0,s1;
	register unsigned long m;
	int i, j;

	if ( ! n ) { _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	_ff_set_one(g[1]);  _ff_set(g[0],a);
	if ( n==1 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	_ff_add(s0,a,a);
	_ff_add(s1,f[0],f[0]);
	while (m) {
		j = ((n&m)>>i);
		if ( j ) {
			// square g and multiply by (x+a) to get (g0^2+g1^2f0+a*2g0g1)*x + a(g0^2+g1^2f0) + 2g0g1f0
			_ff_square(t0,g[0]);				// t0=g0^2
			_ff_square(t1,g[1]);				// t1=g1^2
			_ff_mult(t2,f[0],t1);				// t2 = g1^2f0
			_ff_add(t1,t0,t2);				// t1 = g0^2+g1^2f0
			_ff_mult(t0,g[0],g[1]);			// t0 = g0g1
			_ff_mult(t2,a,t1);				// t2 = a(g0^2+g1^2f0)
			_ff_mult(t3,s1,t0);				// t3 = 2f0g0g1
			_ff_add(g[0],t2,t3);				// g0 = a(g0^2+g1^2f0) + 2g0g1f0
			_ff_mult(t2,t0,s0);				// t2 = 2ag0g1
			_ff_add(g[1],t1,t2);				// g1 = g0^2+g1^2f0+a*2g0g1			
			// 7M + 3A
		} else {
			// square g mod f to get  2g0g1*x + (g0^2+g1^2f0)
			_ff_square(t0, g[0]);  _ff_square(t1, g[1]);
			_ff_mult(t2,g[0],g[1]);
			_ff_add(g[1],t2,t2);
			_ff_mult(t2,t1,f[0]);
			_ff_add(g[0],t0,t2);
			// 4M + 2A
		}
		m >>= 1;  i -= 1;
	}
}

// Computes (x+a)^n modulo a cubic of the form x^3+f1x+f0.  Assumes n < 2^63
void ff_poly_xpan_mod_d3 (ff_t g[3], ff_t a, unsigned long n, ff_t f[4])
{
	_ff_t_declare_reg w1,w2,t2;
	register unsigned long m;
	register int i, j;
	
	if ( ! n ) { _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	_ff_set_zero(g[2]); _ff_set_one(g[1]);  _ff_set(g[0],a);
	if ( n==1 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	while (m) {
		ff_poly_square_mod_d3 (g,g,f);
		j = ((n&m)>>i);
		if ( j ) {	// 5M+4A
			_ff_set(t2,g[2]);
			_ff_mult(w1,t2,a);
			_ff_add(g[2],g[1],w1);
			_ff_mult(w1,g[1],a);
			_ff_mult(w2,t2,f[1]);
			_ff_subfrom(w1,w2);
			_ff_add(g[1],g[0],w1);
			_ff_mult(w1,g[0],a);
			_ff_mult(w2,t2,f[0]);
			_ff_sub(g[0],w1,w2);
		}
		m >>= 1;  i -= 1;
	}
}

// squares a quadratic (arbitrary coeff) modulo a cubic of the form x^3+f1x+f0
// no initialization or validation is done, overlap is ok. 
void ff_poly_square_mod_d3 (ff_t w[3], ff_t u[3], ff_t f[4])
{
	_ff_t_declare_reg t0,t1,s1,s2, w1, w2;

	_ff_mult(t0,u[2],f[0]);
	_ff_mult(t1,u[2],f[1]);
	_ff_mult(s1,t0,u[2]);
	_ff_add(s2,u[1],u[1]);
	
	_ff_add(w1,u[0],u[0]);
	_ff_subfrom(w1,t1);
	_ff_mult(w2,w1,u[2]);
	_ff_square(w1,u[1]);
	_ff_add(w[2],w1,w2);			// w2 = u1^2+(2u0-f1u2)u2 = 2u0u2+u1^2 - f1u2^2
	
	_ff_sub(w1,u[0],t1);
	_ff_mult(w2,w1,s2);
	_ff_sub(u[1],w2,s1);			// w1 = 2u1(u0-f1u2)-(f0u2)u2) = 2u0u2+u1^2 - f1u2^2
	
	_ff_square(w1,u[0]);
	_ff_mult(w2,s2,t0);
	_ff_sub(u[0],w1,w2);			// w0 = u0^2-(2u1)(f0u2)
	 // 8M+7A
}

// Computes (x+a)^n modulo a quartic of the form x^4+f2x^2+f1x+f0.  Assumes n < 2^63
void ff_poly_xpan_mod_d4 (ff_t g[4], ff_t a, unsigned long n, ff_t f[5])
{
	_ff_t_declare_reg w1,w2,t3;
	ff_t h[4];
	register unsigned long m;
	int i, j;
	
	if ( ! n ) { _ff_set_zero(g[1]); _ff_set_one(g[0]); return; }
	_ff_set_zero(g[3]);  _ff_set_zero(g[2]); _ff_set_one(g[1]);  _ff_set(g[0],a);
	if ( n==1 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	while (m) {
		ff_poly_square_mod_d4 (g,g,f);
		j = ((n&m)>>i);
		if ( j ) {
			_ff_set(t3,g[3]);
			_ff_mult(w1,t3,a);
			_ff_add(g[3],g[2],w1);
			_ff_mult(w1,g[2],a);
			_ff_mult(w2,t3,f[2]);
			_ff_subfrom(w1,w2);
			_ff_add(g[2],g[1],w1);
			_ff_mult(w1,g[1],a);
			_ff_mult(w2,t3,f[1]);
			_ff_subfrom(w1,w2);
			_ff_add(g[1],g[0],w1);
			_ff_mult(w1,g[0],a);
			_ff_mult(w2,t3,f[0]);
			_ff_sub(g[0],w1,w2);
			// 7M+6A
		}
		m >>= 1;  i -= 1;
	}
}

// squares a cubic (arbitrary coeff) modulo a quartic of the form x^4+ax^2+bx+c
// no initialization or validation is done, overlap of u,w is OK.  Uses 18M+21A
void ff_poly_square_mod_d4 (ff_t w[4], ff_t u[4], ff_t f[5])
{
	_ff_t_declare_reg t0,t1,t2,s1,s2,s3,s4,s5,w1,w2,w3;
	
	// compute t_i = u_3*f_i
	_ff_mult(t0,u[3],f[0]);  _ff_mult(t1,u[3],f[1]);  _ff_mult(t2,u[3],f[2]);
	_ff_mult(s1,u[2],f[2]);			// s1=u2f2
	_ff_square(s2,u[2]);				// s2=u2^2
	_ff_mult(s3,u[3],t0);			// s3=u3^2f0	(not reused, but need to save since u[3] may get clobbered)
	_ff_mult(s4,u[2],t0);
	_ff_x2(s4);					// s4 = 2u2u3f0 (not reused, but need to save since u[2] may get clobbered)
	_ff_mult(s5,u[1],u[3]);
	_ff_x2(s5);					// s5 = 2u1u3 (not reused, but need to save since u[1],u[3] may get clobbered)

	_ff_sub(w1,u[0],s1);
	_ff_x2(w1);
	_ff_subfrom(w1,t1);
	_ff_mult(w2,w1,u[3]);
	_ff_mult(w1,u[1],u[2]);
	_ff_x2(w1);
	_ff_add(w[3],w1,w2);			// w3 = (2(u0-u2f2)-u3f2)u3+2u1u2 = 2*u0*u3 + 2*u1*u2 - 2*u2*u3*f2 - u3^2*f1
	// note u[3] is now (potentially) clobbered

	_ff_sub(w2,u[1],t2);
	_ff_square(w1,w2);
	_ff_sub(w2,u[0],t1);
	_ff_x2(w2);
	_ff_subfrom(w2,s1);
	_ff_mult(w3,u[2],w2);
	_ff_addto(w1,w3);
	_ff_sub(w[2],w1,s3);			// w2 = (u1-u3f2)^2+u2(2(u0-u3f1)-u2f2)-u3(u3f0) = 2*u0*u2 + u1^2 - 2*u1*u3*f2 - u2^2*f2 - 2*u2*u3*f1 - u3^2*f0 + u3^2*f2^2
	// note u[2] is now (potentially) clobbered
	
	_ff_sub(w2,u[0],t1);
	_ff_mult(w1,u[1],w2);
	_ff_x2(w1);
	_ff_mult(w2,s2,f[1]);
	_ff_subfrom(w1,w2);
	_ff_subfrom(w1,s4);
	_ff_mult(w2,t1,t2);
	_ff_add(w[1],w1,w2);			// w1 = 2u1(u0-u3f1)-u2^2f1-2u2(u3f0)+(u3f1)(u3f2) = 2*u0*u1 - 2*u1*u3*f1 - u2^2*f1 - 2*u2*u3*f0 + u3^2*f1*f2
	// note u[1] is now (potentially) clobbered

	_ff_square(w1,u[0]);
	_ff_add(w2,s5,s2);
	_ff_mult(w3,w2,f[0]);
	_ff_subfrom(w1,w3);
	_ff_mult(w2,t0,t2);
	_ff_add(w[0],w1,w2);			// w0 = u0^2-f0(2u1u3+u2^2)+(u3f0)(u3f2) = u0^2 - 2*u1*u3*f0 - u2^2*f0 + u3^2*f0*f2
	 // 18M+21A
}


/*

The code below has all been superseded by more efficient routines above or in ffext.c
Keep it around in case we need it for debugging.

// squares a quadratic (arbitrary coeff) modulo a cubic of the form x^3+f1x+f0
// no initialization or validation is done, overlap is ok.  Uses 10M+9A
void old_ff_poly_square_mod_d3 (ff_t w[3], ff_t u[3], ff_t f[4])
{
	_ff_t_declare_reg t1,t3, t4, w0, w1, w2;

	_ff_mult(t1,u[0],u[2]);
	_ff_add(w2,t1,t1);
	_ff_square(t1,u[1]);
	_ff_addto(w2,t1);
	_ff_square(t4,u[2]);				// t4 = u2^2
	_ff_mult(t1,f[1],t4);
	_ff_subfrom(w2,t1);				// w2 = 2u0u2+u1^2 - f1u2^2
	_ff_mult(t1,u[0],u[1]);
	_ff_add(w1,t1,t1);				
	_ff_mult(t1,f[0],t4);
	_ff_subfrom(w1,t1);
	_ff_mult(t1,u[1],u[2]);
	_ff_add(t3,t1,t1);				// t3 = 2u1u2
	_ff_mult(t1,f[1],t3);
	_ff_subfrom(w1,t1);				// w1 = 2u0u1 - f0u2^2 - 2f1u1u2
	_ff_square(w0,u[0]);
	_ff_mult(t1,f[0],t3);
	_ff_sub(w[0],w0,t1);				// w0 = u0^2 - 2u1u2f0
	_ff_set(w[1],w1);
	_ff_set(w[2],w2);
}

// Computes (x+a)^n modulo a cubic of the form x^3+f1x+f0.  Assumes n < 2^63
void ff_poly_xpan_mod_d3 (ff_t g[3], ff_t a, unsigned long n, ff_t f[5])
{
	_ff_t_declare p[12], *q;
	_ff_t_declare_reg t1,a2, ca;
	register unsigned long m;
	int i, d;
	
	memset (p, 0, sizeof(p));
	q = p;
	_ff_set_one(q[0]);		// (x+a)^0 = 1
	q += 3;
	_ff_set_one(q[1]);
	_ff_set(q[0],a);		// (x+a)^1 = x + a
	q += 3;
	_ff_set_one(q[2]);
	_ff_add(ca,a,a);
	_ff_set(q[1],ca);
	_ff_square(a2,a);
	_ff_set(q[0],a2);		// (x+a)^2 = x^2 + 2ax + a^2
	q += 3;
	_ff_addto(ca,a);
	_ff_set(q[2],ca);
	_ff_mult(t1,ca,a);
	_ff_sub(q[1],t1,f[1]);
	_ff_mult(t1,a2,a);
	_ff_sub(q[0],t1,f[0]);	// (x+a)^3 = 3ax^2+(3a^2-f1)x + (a^3-f0)

	// could move up to exp8 but we'll settle for a window size of 2 bits for now
	ff_poly_exp4_mod_d3 (g, p, n, f);
}


// Computes x^n modulo a cubic of the form x^3+ax+b.  Assumes n < 2^63
void ff_poly_xn_mod_d3 (ff_t g[3], unsigned long n, ff_t f[4])
{
	_ff_t_declare p[24], *q;
	_ff_t_declare_reg t0, t1, t3;
	register unsigned long m;
	int i, d;
	
	memset (p, 0, sizeof(p));
	q = p;
	_ff_set_one(q[0]);		// x^0 = 1
	q += 3;
	_ff_set_one(q[1]);		// x^1 = x
	q += 3;
	_ff_set_one(q[2]);		// x^2 = x^2
	q += 3;
	_ff_neg(t1,f[1]);
	_ff_neg(t0,f[0]);
	_ff_set(q[1],t1);
	_ff_set(q[0],t0);		// x^3 = -f1x-f0
	q += 3;
	_ff_set(q[2],t1);
	_ff_set(q[1],t0);		// x^4 = -f1x^2-f0x
	q += 3;
	_ff_set(q[2],t0);
	_ff_square(t3,f[1]);
	_ff_set(q[1],t3);		// t3 = f1^2
	_ff_mult(t1,f[0],f[1]);
	_ff_set(q[0],t1);		// x^5 = -f0x^2+f1x^2+f0f1
	q += 3;
	_ff_set(q[2],t3);
	_ff_add(t0,t1,t1);		// t0 = 2f0f1
	_ff_set(q[1],t0);
	_ff_square(t1,f[0]);		// t1 = f0^2
	_ff_set(q[0],t1);		// x^6 = f1^2x^2+2f0f1x+f0^2
	q += 3;
	_ff_set(q[2],t0);
	_ff_mult(t0,f[1],t3);
	_ff_sub(q[1],t1,t0);
	_ff_mult(t0,f[0],t3);
	_ff_neg(q[0],t0);		// x^7 = 2f0f1x^2 + (f0^2-f1^3)x - f0f1^2
	
	ff_poly_exp8_mod_d3 (g, p, n, f);
}


void ff_poly_exp4_mod_d3 (ff_t g[3], ff_t p[12], unsigned long n, ff_t f[4])
{
	register ff_t *q;
	register unsigned long m;
	int i, d;
	
	for ( m = (0x3UL<<60), i = 60 ; ! (n&m) ; m >>= 2, i-=2 );
	q = p + 3*((n&m)>>i);
	_ff_set (g[2],q[2]);  _ff_set(g[1],q[1]);  _ff_set(g[0],q[0]);
	m >>= 2;  i -= 2;
	while (m) {
		ff_poly_square_mod_d3 (g,g,f);
		ff_poly_square_mod_d3 (g,g,f);
		q = p + 3*((n&m)>>i);
		ff_poly_mult_mod_d3 (g, g, q, f);
		m >>= 2;  i -= 2;
	}
}

void ff_poly_exp8_mod_d3 (ff_t g[3], ff_t p[24], unsigned long n, ff_t f[4])
{
	register ff_t *q;
	register unsigned long m;
	int i, d;
	
	for ( m = (0x7UL<<60), i = 60 ; ! (n&m) ; m >>= 3, i-=3 );
	q = p + 3*((n&m)>>i);
	_ff_set (g[2],q[2]);  _ff_set(g[1],q[1]);  _ff_set(g[0],q[0]);

	m >>= 3;  i -= 3;
	while (m) {
		ff_poly_square_mod_d3 (g,g,f);
		ff_poly_square_mod_d3 (g,g,f);
		ff_poly_square_mod_d3 (g,g,f);
		q = p + 3*((n&m)>>i);
		ff_poly_mult_mod_d3 (g, g, q, f);
		m >>= 3;  i -= 3;
	}
}


// multiplies two quadratics (arbitrary coeff) modulo a cubic of the form x^3+ax+b
// no initialization or validation is done, overlap of u,v,w is ok.  Uses 10M+11A
void ff_poly_mult_mod_d3 (ff_t w[3], ff_t u[3], ff_t v[3], ff_t f[4])
{
//	_ff_t_declare_reg t1, t3, t4, w0, w1, w2;
	_ff_t_declare_reg t0, t1, t2, s01, s02, s12, w0, w1;
	
	_ff_mult(t0,u[0],v[0]);  _ff_mult(t1,u[1],v[1]);  _ff_mult(t2,u[2],v[2]);  
	_ff_add(w1,u[0],u[1]); _ff_add(w0,v[0],v[1]); _ff_mult(s01,w1,w0); _ff_subfrom(s01,t0);  _ff_subfrom(s01,t1);
	_ff_add(w1,u[0],u[2]); _ff_add(w0,v[0],v[2]); _ff_mult(s02,w1,w0); _ff_subfrom(s02,t0);  _ff_subfrom(s02,t2);
	_ff_add(w1,u[1],u[2]); _ff_add(w0,v[1],v[2]); _ff_mult(s12,w1,w0); _ff_subfrom(s12,t1);  _ff_subfrom(s12,t2);
	
	_ff_add(w1,s02,t1);
	_ff_mult(w0,f[1],t2);
	_ff_sub(w[2],w1,w0);			// w2 = s02 + t1 - t2f1 = u2v0+u1v1+ u0v2 - f1u2v2
	
	_ff_mult(w1,f[0],t2);
	_ff_mult(w0,f[1],s12);
	_ff_addto(w0,w1);
	_ff_sub(w[1],s01,w0);			// w1= s01 - t2f0 - s12f1= u1v0+u0v1-f0u2v2-f1(u1v2+u2v1)
	
	_ff_mult (w0, f[0],s12);
	_ff_sub(w[0],t0,w0);				// w0 = t0 - s12f0 = u0v0 - (u1v2+u2v1)f0
	return;
}


*/

/*
// multiplies two cubics (arbitrary coeff) modulo a quartic of the form x^4+ax^2+bx+c
// no initialization or validation is done, overlap of u,v,w is ok.  Uses 23M+38A
void ff_poly_mult_mod_d4 (ff_t w[4], ff_t u[4], ff_t v[4], ff_t f[5])
{
	_ff_t_declare_reg t0, t1, t2, t3, s01,s02,s03,s12,s13,s23,w0, w1;		// don't need so many variables

	// compute t_i = u_i*v_i
	_ff_mult(t0,u[0],v[0]); _ff_mult(t1,u[1],v[1]); _ff_mult(t2,u[2],v[2]); _ff_mult(t3,u[3],v[3]);
	
	// compute s_ij = u_i*v_j + u_j*v_i
	_ff_add(w0,u[0],u[3]); _ff_add(w1,v[0],v[3]); _ff_mult(s03,w0,w1); _ff_subfrom(s03,t0);  _ff_subfrom(s03,t3);
	_ff_add(w0,u[0],u[2]); _ff_add(w1,v[0],v[2]); _ff_mult(s02,w0,w1); _ff_subfrom(s02,t0);  _ff_subfrom(s02,t2);
	_ff_add(w0,u[0],u[1]); _ff_add(w1,v[0],v[1]); _ff_mult(s01,w0,w1); _ff_subfrom(s01,t0);  _ff_subfrom(s01,t1);
	_ff_add(w0,u[1],u[3]); _ff_add(w1,v[1],v[3]); _ff_mult(s13,w0,w1); _ff_subfrom(s13,t1);  _ff_subfrom(s13,t3);
	_ff_add(w0,u[1],u[2]); _ff_add(w1,v[1],v[2]); _ff_mult(s12,w0,w1); _ff_subfrom(s12,t1);  _ff_subfrom(s12,t2);
	_ff_add(w0,u[2],u[3]); _ff_add(w1,v[2],v[3]); _ff_mult(s23,w0,w1); _ff_subfrom(s23,t2);  _ff_subfrom(s23,t3);
	
	// s13 is special, add t2 to it
	_ff_addto(s13,t2);
	
	_ff_mult(w0,s23,f[2]);
	_ff_mult(w1,t3,f[1]);
	_ff_addto(w0,w1);
	_ff_add(w1,s03,s12);
	_ff_sub(w[3],w1,w0);			// w3 = s03 + s12 - s23f2 - t3f1 = u0v3 + u1v2 + u2v1 - u2v3f2 + u3v0 - u3v2f2 - u3v3f1

	_ff_mult(w1,s13,f[2]);
	_ff_mult(w0,s23,f[1]);
	_ff_addto(w0,w1);
	_ff_square(w1,f[2]);
	_ff_subfrom(w1,f[0]);
	ff_mult(w1,w1,t3);
	_ff_addto(w1,s02);
	_ff_addto(w1,t1);
	_ff_sub(w[2],w1,w0);			// w2 = s02 + t1 - s13f2 - s23f1 + t3(f2^2-f0) = u0v2 + u1v1 - u1v3f2 + u2v0 - u2v2f2 - u2v3f1 - u3v1f2 - u3v2f1 - u3v3f0 + u3v3f2^2
	
	_ff_mult(w1,s13,f[1]);
	_ff_mult(w0,s23,f[0]);
	_ff_addto(w0,w1);
	_ff_mult(w1,f[1],f[2]);
	ff_mult(w1,w1,t3);
	_ff_addto(w1,s01);
	_ff_sub(w[1],w1,w0);			// w1 = s01 - s13f1 - s23f0 + t3f1f2  = u0v1 + u1v0 - u1v3f1 - u2v2f1 - u2v3f0 - u3v1f1 - u3v2f0 + u3v3f1f2

	_ff_mult(w0,s13,f[0]);
	_ff_mult(w1,f[0],f[2]);
	ff_mult(w1,w1,t3);
	_ff_addto(w1,t0);
	_ff_sub(w[0],w1,w0);			// w0 = t0 - s13f0 + t3f0f2 = u0v0 - u1v3f0 - u2v2f0 - u3v1f0 + u3v3f0f2
}


// squares a cubic (arbitrary coeff) modulo a quartic of the form x^4+ax^2+bx+c
// no initialization or validation is done, overlap of u,w is OK.  Uses 23M+21A
void old_ff_poly_square_mod_d4 (ff_t w[4], ff_t u[4], ff_t f[5])
{
	_ff_t_declare_reg t0,t1,t2,t3,s03,s02,s01,s13,s12,s23,w0,w1;

	// compute t_i = u_i^2
	_ff_square(t0,u[0]);	_ff_square(t1,u[1]);	_ff_square(t2,u[2]);	_ff_square(t3,u[3]);
	
	// compute s_ij = u_i*u_j
	_ff_mult(s13,u[1],u[3]);  _ff_mult(s23,u[2],u[3]);
	
	_ff_mult(w1,u[0],u[3]);
	_ff_mult(w0,u[1],u[2]);
	_ff_addto(w1,w0);
	_ff_mult(w0,s23,f[2]);
	_ff_subfrom(w1,w0);
	_ff_x2(w1);
	_ff_mult(w0,t3,f[1]);
	_ff_sub(w[3],w1,w0);			// w3 = 2(s03 + s12 - s23f2) - t3f1 = 2*u0*u3 + 2*u1*u2 - 2*u2*u3*f2 - u3^2*f1
	
	_ff_mult(w0,s13,f[2]);
	_ff_mult(w1,s23,f[1]);
	_ff_addto(w0,w1);
	_ff_mult(w1,u[0],u[2]);
	_ff_subfrom(w1,w0);
	_ff_x2(w1);
	_ff_addto(w1,t1);
	_ff_square(w0,f[2]);
	_ff_subfrom(w0,f[0]);
	ff_mult(w0,w0,t3);
	_ff_addto(w1,w0);
	_ff_mult(w0,t2,f[2]);
	_ff_sub(w[2],w1,w0);			// w2 = 2(s02 - s12f2 - s23f1) + t1 - t2f2 + t3(f2^2-f0) = 2*u0*u2 + u1^2 - 2*u1*u3*f2 - u2^2*f2 - 2*u2*u3*f1 - u3^2*f0 + u3^2*f2^2
	
	_ff_mult(w0,s23,f[0]);
	_ff_mult(w1,u[0],u[1]);
	_ff_subfrom(w1,w0);
	_ff_x2(w1);
	_ff_mult(w0,f[1],f[2]);
	ff_mult(w0,w0,t3);
	_ff_addto(w1,w0);
	_ff_add(w0,s13,s13);
	_ff_addto(w0,t2);
	ff_mult(w0,w0,f[1]);
	_ff_sub(w[1],w1,w0);			// w1 = 2(s01 - s23f0) - (t2 + 2s13)f1 + t3f1f2 = 2*u0*u1 - 2*u1*u3*f1 - u2^2*f1 - 2*u2*u3*f0 + u3^2*f1*f2
	
	_ff_mult(w0,t3,f[2]);
	_ff_add(w1,s13,s13);
	_ff_addto(w1,t2);
	_ff_subfrom(w0,w1);
	_ff_mult(w1,w0,f[0]);
	_ff_add(w[0],w1,t0);			// w0 = t0 + (t3f2 - 2s13 - t2)f0 = u0^2 - 2*u1*u3*f0 - u2^2*f0 + u3^2*f0*f2
}
*/
/*
void old_ff_poly_xpan_mod_d4 (ff_t g[4], ff_t a, unsigned long n, ff_t f[5])
{
	_ff_t_declare p[16], *q;
	_ff_t_declare_reg a2, ca;
	register unsigned long m;
	int i, d;
	
	memset (p, 0, sizeof(p));
	q = p;
	_ff_set_one(q[0]);		// (x+a)^0 = 1
	q += 4;
	_ff_set_one(q[1]);
	_ff_set(q[0],a);		// (x+a)^1 = x + a
	q += 4;
	_ff_set_one(q[2]);
	_ff_add(ca,a,a);
	_ff_set(q[1],ca);
	_ff_square(a2,a);
	_ff_set(q[0],a2);		// (x+a)^2 = x^2 + 2ax + a^2
	q += 4;
	_ff_set_one(q[3]);		
	_ff_addto(ca,a);
	_ff_set(q[2],ca);
	_ff_mult(q[1],ca,a);
	_ff_mult(q[0],a2,a);		// (x+a)^3 = x^3 + 3ax^2+3a^2x+a^3

	// could move up to exp8 but we'll settle for a window size of 2 bits for now
	ff_poly_exp4_mod_d4 (g, p, n, f);
}


void ff_poly_exp4_mod_d4 (ff_t g[4], ff_t p[16], unsigned long n, ff_t f[5])
{
	register ff_t *q;
	register unsigned long m;
	int i, d;
	
	for ( m = (0x3UL<<60), i = 60 ; ! (n&m) ; m >>= 2, i-=2 );
	q = p + 4*((n&m)>>i);
	_ff_set(g[3],q[3]);  _ff_set (g[2],q[2]);  _ff_set(g[1],q[1]);  _ff_set(g[0],q[0]);
	m >>= 2;  i -= 2;
	while (m) {
		ff_poly_square_mod_d4 (g,g,f);
		ff_poly_square_mod_d4 (g,g,f);
		q = p + 4*((n&m)>>i);
		ff_poly_mult_mod_d4 (g, g, q, f);
		m >>= 2;  i -= 2;
	}
}

void ff_poly_exp8_mod_d4 (ff_t g[4], ff_t p[32], unsigned long n, ff_t f[5])
{
	register ff_t *q;
	register unsigned long m;
	int i, d;
	
	for ( m = (0x7UL<<60), i = 60 ; ! (n&m) ; m >>= 3, i-=3 );
	q = p + 4*((n&m)>>i);
	_ff_set(g[3],q[3]);  _ff_set (g[2],q[2]);  _ff_set(g[1],q[1]);  _ff_set(g[0],q[0]);
	m >>= 3;  i -= 3;
	while (m) {
		ff_poly_square_mod_d4 (g,g,f);
		ff_poly_square_mod_d4 (g,g,f);
		ff_poly_square_mod_d4 (g,g,f);
		q = p + 4*((n&m)>>i);
		ff_poly_mult_mod_d4 (g, g, q, f);
		m >>= 3;  i -= 3;
	}
}
*/

/*
// Computes x^n modulo a quartic of the form x^4+f2x^2+f1x+f0.  Assumes n < 2^63
void old_ff_poly_xn_mod_d4 (ff_t g[4], unsigned long n, ff_t f[5])
{
	_ff_t_declare p[32], *q;
	_ff_t_declare_reg t0, t1, t2, t3;
	register unsigned long m;
	int i, d;
	
	memset (p, 0, sizeof(p));
	q = p;
	_ff_set_one(q[0]);		// x^0 = 1
	q += 4;
	_ff_set_one(q[1]);		// x^1 = x
	q += 4;
	_ff_set_one(q[2]);		// x^2 = x^2
	q += 4;
	_ff_set_one(q[3]);		// x^3 = x^3
	q += 4;
	_ff_neg(t2,f[2]);
	_ff_neg(t1,f[1]);
	_ff_neg(t0,f[0]);
	_ff_set(q[2],t2);
	_ff_set(q[1],t1);
	_ff_set(q[0],t0);		// x^4 = -f2x^2 - f1x-f0
	q += 4;
	_ff_set(q[3],t2);
	_ff_set(q[2],t1);
	_ff_set(q[1],t0);		// x^5 = -f2x^3 - f1x^2 - f0x
	q += 4;
	_ff_set(q[3],t1);
	_ff_square(t3,f[2]);
	_ff_subfrom(t3,f[0]);
	_ff_set(q[2],t3);
	_ff_mult(t1,f[1],f[2]);
	_ff_set(q[1],t1);
	_ff_mult(t0,f[0],f[2]);
	_ff_set(q[0],t0);		// x^6 = -f1x^3 + (f2^2-f0)x^2 + f1f2x^2 + f0f2
	q += 4;
	_ff_set(q[3],t3);
	_ff_add(q[2],t1,t1);
	_ff_square(t2,f[1]);
	_ff_add(q[1],t0,t2);
	_ff_mult(q[0],f[0],f[1]);	// x^7 = (f2^2-f0)x^3 + 2f1f2x^2 + (f0f2+f1^2)x + f0f1
	
	ff_poly_exp8_mod_d4 (g, p, n, f);
}
*/
