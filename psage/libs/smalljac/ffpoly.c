#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include "gmp.h"
#include "mpzutil.h"
#include "ffwrapper.h"
#include "polydisc.h"
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

// This module consists mostly of standard polynomial arithmetic functions that are better implemented in NTL
// but are handy for a standalone C implementation.

// The discriminant code is woefully incomplete and resultants aren't implemented.

int mpq_poly_jinv_i (mpq_t j, long a[5])
{
	static mpz_t t1, t2, a0, a1, a2, a3, a4, a6, b2, b4, b6, b8, c4, D;  // more variables than is really necessary
	static mpq_t d;
	static int init;

	if ( ! init ) {
		mpz_init (a0); mpz_init (a1);  mpz_init (a2);  mpz_init (a3);  mpz_init (a4);  mpz_init (a6);
		mpz_init (b2);  mpz_init (b4);  mpz_init (b6);  mpz_init (b8);  mpz_init (c4);  mpz_init (D);
		mpz_init (t1); mpz_init (t2);  mpq_init(d);  init = 1;
	}
	mpz_set_i (a1, a[0]);  mpz_set_i (a2, a[1]);  mpz_set_i (a3, a[2]);  mpz_set_i (a4, a[3]);  mpz_set_i (a6, a[4]);
	
	// compute b2 = a1^2+4a2
	mpz_mul(t2,a1,a1);			// save t2 = a1^2
	mpz_mul_2exp (t1,a2,2);
	mpz_add(b2,t2,t1);
	// compute b4 = 2a4+a1a3
	mpz_mul_2exp (b4,a4,1);
	mpz_mul (t1,a1,a3);
	mpz_add (b4,b4,t1);
	// compute b6 = a3^2+4a6
	mpz_mul (b6,a3,a3);
	mpz_mul_2exp (t1,a6,2);
	mpz_add (b6, b6, t1);
	// compute b8 = (a1^2+4a2)a6 + (a2a3-a1a4)a3 - a4^2
	mpz_mul_2exp (b8,a2,2);
	mpz_add (b8,b8,t2);		// t2 = a1^2 from above
	mpz_mul (b8,b8,a6);
	mpz_mul (t1,a2,a3);
	mpz_mul (t2,a1,a4);
	mpz_sub (t1,t1,t2);
	mpz_mul (t1,t1,a3);
	mpz_add (b8,b8,t1);
	mpz_mul (t1,a4,a4);
	mpz_sub (b8,b8,t1);
	// compute D = (9b4b6-b2b8)b2 - 8b4^3 - 27b6^2
	mpz_mul (D,b4,b6);
	mpz_mul_ui (D,D,9);
	mpz_mul (t1,b2,b8);
	mpz_sub (D,D,t1);
	mpz_mul (D,D,b2);
	mpz_mul (t1,b4,b4);
	mpz_mul (t1,t1,b4);
	mpz_mul_2exp (t1,t1,3);
	mpz_sub (D,D,t1);
	mpz_mul (t1,b6,b6);
	mpz_mul_ui (t1,t1,27);
	mpz_sub (D,D,t1);
	if ( ! mpz_sgn(D) ) return 0;
	// compute t1 = numerator of j = c4^3 = (b2^2 - 24b4)^3
	mpz_mul (t1,b2,b2);
	mpz_mul_ui (t2,b4,24);
	mpz_sub (t1,t1,t2);
	mpz_mul (t2, t1, t1);
	mpz_mul (t1,t1,t2);
	mpq_set_z (j, t1);
	mpq_set_z (d,D);
	mpq_div (j,j,d);
	return 1;
}


/*
	For monic f, substitute x <- (x - f_{d-1}/[d*f_d]) to make f_{d-1} term zero if necessary
	Returns 1 if poly modified, 0 if not
*/
int mpq_poly_standardize (mpq_t f[], int d)
{
	static mpq_t e[POLY_MAX_DEGREE+1], g[POLY_MAX_DEGREE+1], h[POLY_MAX_DEGREE+1], c;
	static int init;
	int i, j, d_e;
	
	if ( ! init ) { for ( i = 0 ; i <= POLY_MAX_DEGREE ; i++ ) { mpq_init (e[i]);  mpq_init (g[i]);  mpq_init (h[i]); } mpq_init (c);  init = 1; }
	if ( d < 1 ) return 0;
	if ( mpz_cmp_ui (mpq_numref(f[d]),1) != 0 || mpz_cmp_ui(mpq_denref(f[d]),1) != 0 ) return 0;
	if ( ! mpq_sgn(f[d-1]) ) return 0;
	mpq_set_ui (c, 1, d);
	mpq_div (c, c, f[d]);
	mpq_mul (c, c, f[d-1]);
	mpq_neg (c, c);				// c = -f_{d-1}/(d*f_d)
	for ( i = 0 ; i <= d ; i++ ) { mpq_set_ui (e[i], 0, 1);  mpq_set_ui (g[i], 0, 1);  mpq_set_ui (h[i], 0, 1); }
	mpq_set_ui (e[0], 1, 1);  d_e = 0;	// e = 1
	mpq_set (g[0], f[0]);			// g = f0
	for ( j = 1 ; j <= d ; j++ ) {
		for ( i = 0 ; i <= d_e ; i++ ) mpq_mul (h[i], c, e[i]);		// h = c*e
		for ( i = d_e ; i >= 0 ; i-- ) mpq_set (e[i+1], e[i]);			// e *= x
		mpq_set (e[0], h[0]);
		for ( i = 1 ; i <= d_e ; i++ ) mpq_add (e[i], e[i], h[i]);		// e += h
		d_e++;											// e is now e*(x+c)
		for ( i = 0 ; i <= d_e ; i++ ) mpq_mul (h[i], f[j], e[i]);		// h = f[j]*e
		for ( i = 0 ; i <= d_e ; i++ ) mpq_add (g[i], g[i], h[i]);		// g += h
	}
	for ( i = 0 ; i <= d ; i++ ) mpq_set (f[i], g[i]);					// copy g into f
	if ( mpq_sgn(f[d-1]) ) { err_printf ("mpq_poly_standardize failed!  f[d-1] = %Qd \n", f[d-1]);  exit (0); }		// should be impossible
	return 1;
}


/*
	Convert long weierstrass form [a1,a2,a3,a4,a6] to short form x^3+f1*x + f0
	Weierstrass coefficients are specified by W[0] = a1, ..., W[3] = a4, w[4] = a6

	This is done by first applying the substitution y -> y - a1/2*x - a3/2 which kills a1 and a3 
	and sends a2 -> a2+ a1^2/4, a4 -> a4 + a1a3/2 and a6 -> a6 + a3^2/4.
	
	Then apply the substitution x -> x - a2/3 to the resulting curve, which kills a2, keeps
	a1 and a3 zero if they have already been killed, and sends a4 -> a4 - a2^2/3 and
	a6 -> a6 + 2/27a2^3 - a2a4/3
*/
void mpq_poly_weierstrass (mpq_t f[4], mpz_t W[5])
{
	static mpq_t t1, t2, c2, c4, c6, three;
	static int init;
	
	if ( ! init ) { mpq_init(t1);  mpq_init(t2);  mpq_init(c2);  mpq_init(c4);  mpq_init (c6);  mpq_init(three);  init = 1; }

	mpq_set_ui (f[3], 1, 1);
	mpq_set_ui (f[2], 0, 1);
	
	// set c2 = a2 + a1^2/4
	mpq_set_z (t2, W[0]);
	mpq_div_2exp (t1, t2, 1);						// t1 = a1/2
	mpq_mul (c2, t1, t1);						// c2 = a1^2/4
	mpq_set_z (t2, W[1]);
	mpq_add (c2, c2, t2);						// c2 = a2+a1^2/4 is the new a2

	// set c4 = a4 + a1a3/2
	mpq_set_z (t2, W[2]);
	mpq_mul (c4, t1, t2);						// c4 = a1a3/2
	mpq_set_z (t2, W[3]);
	mpq_add (c4, c4, t2);						// c4 = a4 + a1a3/2 is the new a4
	
	// set c6 = a6 + a3^2/4 this is rolled in below
	mpq_set_z (t2, W[2]);
	mpq_div_2exp (c6, t2, 1);
	mpq_mul (c6, c6, c6);						// c6 = a3^2/4
	mpq_set_z (t2, W[4]);
	mpq_add (c6, c6, t2);						// c6 = a6 + a3^2/4 is the new a6 
	
	// set f[1] = a4' = c4 - c2^2/3
	mpq_set_ui (three, 3, 1);
	mpq_div (c2, c2, three);
	mpq_mul (t1, c2, c2);						// t1 = c2^2/9
	mpq_mul (t2, t1, three);						// t2 = c2^2/3
	mpq_sub (f[1], c4, t2);						// f[1] = c4 - c2^2/3
	
	// set f[0] = a6' = a6 + a3^2/4 + 2(c2/3)^3 - c2c4/3
	mpq_mul (t2, t1, c2);						// t2 = c2^3/27
	mpq_mul_2exp (t2, t2, 1);					// t2 = 2/27*c2^3
	mpq_mul (t1, c4, c2);						// t1 = c2c4/3
	mpq_sub (t2, t2, t1);						// t2 = 2/27*c2^3 - c2c4/3
	mpq_add (f[0], c6, t2);						// f[0] = c6 + 2/27*c2^3 - c2c4/3

	return;
}	

void mpq_poly_weierstrass_mpq (mpq_t f[4], mpq_t a1, mpq_t a2, mpq_t a3, mpq_t a4, mpq_t a6)
{
	static mpq_t t1, t2, c2, c4, c6, three;
	static int init;
	
	if ( ! init ) { mpq_init(t1);  mpq_init(t2);  mpq_init(c2);  mpq_init(c4);  mpq_init (c6);  mpq_init(three);  init = 1; }
	
	mpq_set_ui (f[3], 1, 1);
	mpq_set_ui (f[2], 0, 1);
	
	// set c2 = a2 + a1^2/4
	mpq_div_2exp (t1, a1, 1);					// t1 = a1/2
	mpq_mul (c2, t1, t1);						// c2 = a1^2/4
	mpq_add (c2, c2, a2);						// c2 = a2+a1^2/4 is the new a2

	// set c4 = a4 + a1a3/2
	mpq_mul (c4, t1, a3);						// c4 = a1a3/2
	mpq_add (c4, c4, a4);						// c4 = a4 + a1a3/2 is the new a4
	
	// set c6 = a6 + a3^2/4 this is rolled in below
	mpq_div_2exp (c6, a3, 1);
	mpq_mul (c6, c6, c6);						// c6 = a3^2/4
	mpq_add (c6, c6, a6);						// c6 = a6 + a3^2/4 is the new a6 
	
	// set f[1] = a4' = c4 - c2^2/3
	mpq_set_ui (three, 3, 1);
	mpq_div (c2, c2, three);
	mpq_mul (t1, c2, c2);						// t1 = c2^2/9
	mpq_mul (t2, t1, three);						// t2 = c2^2/3
	mpq_sub (f[1], c4, t2);						// f[1] = c4 - c2^2/3
	
	// set f[0] = a6' = a6 + a3^2/4 + 2(c2/3)^3 - c2c4/3
	mpq_mul (t2, t1, c2);						// t2 = c2^3/27
	mpq_mul_2exp (t2, t2, 1);					// t2 = 2/27*c2^3
	mpq_mul (t1, c4, c2);						// t1 = c2c4/3
	mpq_sub (t2, t2, t1);						// t2 = 2/27*c2^3 - c2c4/3
	mpq_add (f[0], c6, t2);						// f[0] = c6 + 2/27*c2^3 - c2c4/3
	return;
}	


int poly_expr_degree (char *str, int *ws)
{
	register char *s;
	int d, degree;
	
	if ( ws ) *ws = 0;
	for ( s = str ; *s && *s != 'x' && *s != ',' ; s++ );
	if ( *s == ',' && str[0] == '[' ) { if ( ws ) *ws = 1;  return 3; }
	if ( *s != 'x' ) return 0;	// don't try to distinguish bad or zero polys from constants
	if ( *(s+1) == '^' ) degree = atoi (s+2); else degree = 1;
	s++;
	for (;;) {
		while ( *s && *s != 'x' ) s++;
		if ( ! *s ) break;
		if ( *(s+1) == '^' ) { d = atoi (s+2);  if ( d > degree ) degree = d; }
		s++;
	}
	return degree;
}


static inline int isqdigit(char c) { return (isdigit(c) || c == '/' ); }

int mpq_poly_expr (mpq_t f[], int maxd, char *expr)
{
	static int init;
	static mpq_t c;
	static mpz_t W[5];
	char cbuf[4096];
	register char *s, *t;
	int i;
	int sign;
	
	if ( ! init ) { mpq_init (c);  for ( i = 0 ; i < 5 ; i++ ) mpz_init (W[i]); init = 1; }
	
	s = expr;
	// Check if expr looks like it is in weierstrass form "[a1,a2,a3,a4,a6]" and handle this seperately
	while ( isspace(*s) ) s++;
	if ( *s == '[' && (isdigit (*(s+1)) || *(s+1) == '-') ) {
		for ( t = s ; *t && *t != ',' && *t != ']' ; t++ );
		if ( *t == ',' ) {
			if ( gmp_sscanf (s, "[%Zd,%Zd,%Zd,%Zd,%Zd]", W[0], W[1], W[2], W[3], W[4]) != 5 ) return -2;
			mpq_poly_weierstrass (f, W);
			return 3;
		}
	}
	
	for ( i = 0 ; i <= maxd ; i++ ) mpq_set_ui(f[i], 0, 1);
	if ( *s == '[' || *s == '(' ) s++;
	for(;;) {
		while ( isspace(*s) ) s++;
		// handle signs explicitly to deal with things like '+ -3' which gmp doesn't like and we need to handle
		sign = 1;
		if ( *s == '-' ) {
			sign = -1;
			for ( s++ ; isspace(*s) ; s++);
		} else if ( *s == '+' ) {
			for ( s++ ; isspace(*s) ; s++);
			if ( *s == '-' ) {
				sign = -1;
				for ( s++ ; isspace(*s) ; s++);
			}
		}
		if ( isdigit(*s) ) {
			for ( t = cbuf ; isqdigit(*s) ; *t++ = *s++ );
			*t = '\0';
			if ( mpq_set_str (c, cbuf, 0) < 0 ) return -2;
			mpq_canonicalize (c);
			while ( isqdigit(*s) ) s++;
		} else if ( *s == 'x' || *s == 'X' ) {
			mpq_set_ui (c, 1, 1);
		} else {
			break;
		}
		if ( sign < 0 ) mpq_neg(c,c);
		if ( *s == '*' ) s++;
		if ( *s == 'x' || *s == 'X' ) {
			s++;
			if ( *s == '^' ) s++;
			if ( isdigit(*s) ) {
				i = atoi (s);
				while ( isdigit(*s) ) s++;
				if ( i > maxd ) return -2;
			} else {
				i = 1;
			}
		} else {
			i = 0;
		}
		mpq_add (f[i], f[i], c);
	}
	for ( i = maxd ; i && ! mpq_sgn(f[i]) ; i-- );
	if ( ! i && ! mpq_sgn(f[i]) ) return -1;
	maxd = i;
	if  ( *s == ']' && *(s+1) == '/' ) {
		for ( t = cbuf, s+=2 ; isdigit(*s) ; *t++ = *s++ );
		*t = '\0';
		if ( mpq_set_str(c, cbuf, 0) < 0 ) return -2;
		mpq_canonicalize (c);
		for ( i = 0 ; i <= maxd ; i++ ) mpq_div (f[i], f[i], c);
	}
	return maxd;
}


int mpz_poly_expr (mpz_t f[], int maxd, char *expr)
{
	static int init;
	static mpz_t c;
	char cbuf[4096];
	char *s, *t;
	int i;
	int sign;
	
	if ( ! init ) { mpz_init (c);  init = 1; }
	for ( i = 0 ; i <= maxd ; i++ ) mpz_set_ui(f[i], 0);
	s = expr;
	if ( *s == '[' || *s == '(' ) s++;
	for(;;) {
		while ( isspace(*s) ) s++;
		// handle signs explicitly to deal with things like '+ -3' which gmp doesn't like and we need to handle
		sign = 1;
		if ( *s == '-' ) {
			sign = -1;
			for ( s++ ; isspace(*s) ; s++);
		} else if ( *s == '+' ) {
			for ( s++ ; isspace(*s) ; s++);
			if ( *s == '-' ) {
				sign = -1;
				for ( s++ ; isspace(*s) ; s++);
			}
		}
		if ( isdigit(*s) ) {
			for ( t = cbuf ; isdigit(*s) ; *t++ = *s++ );
			*t = '\0';
			if ( mpz_set_str (c, cbuf, 0) < 0 ) return -2;
			while ( isdigit(*s) ) s++;
		} else if ( *s == 'x' || *s == 'X' ) {
			mpz_set_ui (c, 1);
		} else {
			break;
		}
		if ( sign < 0 ) mpz_neg(c,c);
		if ( *s == '*' ) s++;
		if ( *s == 'x' || *s == 'X' ) {
			s++;
			if ( *s == '^' ) s++;
			if ( isdigit(*s) ) {
				i = atoi (s);
				while ( isdigit(*s) ) s++;
				if ( i > maxd ) return -2;
			} else {
				i = 1;
			}
		} else {
			i = 0;
		}
		mpz_add (f[i], f[i], c);
	}
	for ( i = maxd ; i && ! mpz_sgn(f[i]) ; i-- );
	if ( ! i && ! mpz_sgn(f[i]) ) return -1;
	return i;
}


int ui_poly_expr (unsigned long f[], int maxd, char *expr)
{
	static int init;
	static mpz_t F[POLY_MAX_DEGREE+1];
	int i, d;
	
	if ( ! init ) { for ( d = 0 ; d <= POLY_MAX_DEGREE ; d++ ) mpz_init (F[d]); init = 1; }
	d = mpz_poly_expr (F, maxd, expr);
	if ( d < -1 ) return d;
	for ( i = 0 ; i <= d ; i++ ) {
		if ( mpz_sgn(F[i]) < 0 ) { err_printf("Attempt to parse poly with negative coefficients in ui_poly_expr\n");  exit (0); }
		f[i] = mpz_get_ui(F[i]);
	}
	for ( ; i <= maxd ; i++ ) f[i] = 0;
	return d;
}


int ui_poly_expr_mod_p (unsigned long f[], int maxd, char *expr, unsigned long p)
{
	static int init;
	static mpq_t F[POLY_MAX_DEGREE+1];
	static mpz_t t;
	unsigned long x;
	int i, d;
	
	if ( ! init ) { for ( d = 0 ; d <= POLY_MAX_DEGREE ; d++ ) mpq_init (F[d]); mpz_init (t); init = 1; }
	d = mpq_poly_expr (F, maxd, expr);
	if ( d < -1 ) return d;
	for ( i = 0 ; i <= d ; i++ ) {
		x = ui_inverse (mpz_fdiv_ui(mpq_denref(F[i]), p), p);
		mpz_mul_ui (t, mpq_numref(F[i]), x);
		f[i] = mpz_fdiv_ui (t, p);
	}
	for ( ; i <= maxd ; i++ ) f[i] = 0;
	return d;
}


int i_poly_expr (long f[], int maxd, char *expr)
{
	static int init;
	static mpz_t F[POLY_MAX_DEGREE+1];
	int i, d;
	
	if ( ! init ) { for ( d = 0 ; d <= POLY_MAX_DEGREE ; d++ ) mpz_init (F[d]); init = 1; }
	d = mpz_poly_expr (F, maxd, expr);
	if ( d < -1 ) return d;
	for ( i = 0 ; i <= d ; i++ ) f[i] = mpz_get_i(F[i]);
	for ( ; i <= maxd ; i++ ) f[i] = 0;
	return d;
}


int ff_poly_expr (ff_t f[], int maxd, char *expr)
{
	static int init;
	static mpz_t F[POLY_MAX_DEGREE+1];
	int i, d;
	
	if ( ! init ) { for ( d = 0 ; d <= POLY_MAX_DEGREE ; d++ ) mpz_init (F[d]); init = 1; }
	d = mpz_poly_expr (F, maxd, expr);
	if ( d < -1 ) return d;
	for ( i = 0 ; i <= d ; i++ ) _ff_set_mpz(f[i], F[i]);
	for ( ; i <= maxd ; i++ ) _ff_set_zero (f[i]);
	return d;
}


void ff_poly_print (ff_t f[], int d_f)
{
	char buf[4096];

	ff_poly_sprint (buf, f, d_f);
	out_printf ("%s\n", buf);
}


void ui_poly_print (unsigned long f[], int d_f)
{
	char buf[4096];

	ui_poly_sprint (buf, f, d_f);
	out_printf ("%s\n", buf);
}


int ui_poly_sprint (char *s, unsigned long f[], int d_f)
{
	char *t;
	int i, n;
	
	if ( d_f < 0 ) { strcpy (s, "[zero polynomial]");  return strlen(s); }
	t = s;
	if ( d_f >= 2 ) {
		if ( f[d_f] != 1 ) t += sprintf (t, "[%lux^%d", f[d_f], d_f); else  t += sprintf (t, "[x^%d", d_f);
	} else if ( d_f == 1 ) {
		if ( f[d_f] != 1 ) t += sprintf (t, "[%lux", f[d_f]); else  t += sprintf (t, "[x");
	} else {
		t += sprintf (t, "[%lu", f[d_f]);
	}
	for ( i = d_f-1 ; i >= 0 ; i-- ) {
		if ( f[i] ) {
			if ( i >= 2 ) {
				t += sprintf (t, " + %lux^%d", f[i], i);
			} else if ( i == 1 ) {
				t += sprintf (t, " + %lux", f[i]);
			} else {
				t += sprintf (t, " + %lu", f[i]);
			}
		}
	}
	*t++ = ']';
	*t= '\0';
	return t-s;
}


void i_poly_print (long f[], int d_f)
{
	char buf[4096];

	i_poly_sprint (buf, f, d_f);
	out_printf ("%s\n", buf);
}

int i_poly_sprint (char *s, long f[], int d_f)
{
	char *t;
	int i, n;
	
	if ( d_f < 0 ) { strcpy (s, "[zero polynomial]");  return strlen(s); }
	t = s;
	if ( d_f >= 2 ) {
		if ( f[d_f] != 1 ) t += sprintf (t, "[%ldx^%d", f[d_f], d_f); else  t += sprintf (t, "[x^%d", d_f);
	} else if ( d_f == 1 ) {
		if ( f[d_f] != 1 ) t += sprintf (t, "[%ldx", f[d_f]); else  t += sprintf (t, "[x");
	} else {
		t += sprintf (t, "[%ld", f[d_f]);
	}
	for ( i = d_f-1 ; i >= 0 ; i-- ) {
		if ( f[i] > 1) {
			t += sprintf (t, " + %ld", f[i]);
		} else if ( f[i] == 1 ) {
			t += sprintf (t, " + ");
		} else if ( f[i] == -1 ) {
			t += sprintf (t, " - ");
		} else if ( f[i] < -1 ) {
			t += sprintf (t, " - %ld", -f[i]);
		}
		if ( f[i] ) {
			if ( i >= 2 ) {
				t += sprintf (t, "x^%d", i);
			} else if ( i == 1 ) {
				t += sprintf (t, "x", f[i]);
			} else {
				if ( f[i] == 1 || f[i] == -1 ) t += sprintf (t, "1");
			}
		}
	}
	*t++ = ']';
	*t= '\0';
	return t-s;
}


void mpz_poly_print (mpz_t f[], int d_f)
{
	char buf[65536];

	mpz_poly_sprint (buf, f, d_f);
	out_printf ("%s\n", buf);
}

int mpz_poly_sprint (char *s, mpz_t f[], int d_f)
{
	char *t;
	int i, n;
	
	if ( d_f < 0 ) { strcpy (s, "[zero polynomial]");  return strlen(s); }
	t = s;
	if ( d_f >= 2 ) {
		if ( mpz_cmp_ui(f[d_f],1) != 0 ) {
			t += gmp_sprintf (t, "[%Zdx^%d", f[d_f], d_f);
		} else {
			t += gmp_sprintf (t, "[x^%d", d_f);
		}
	} else if ( d_f == 1 ) {
		if ( mpz_cmp_ui(f[d_f],1) != 0 ) {
			t += gmp_sprintf (t, "[%Zdx", f[d_f]);
		} else {
			t += gmp_sprintf (t, "[x");
		}
	} else {
		t += gmp_sprintf (t, "[%Zd", f[d_f]);
	}
	for ( i = d_f-1 ; i >= 0 ; i-- ) {
		if ( mpz_sgn(f[i]) ) {
			if ( i >= 2 ) {
				if ( mpz_cmp_ui(f[i],1) != 0 ) {
					t += gmp_sprintf (t, " + %Zdx^%d", f[i], i);
				} else {
					t += gmp_sprintf (t, " + x^%d", i);
				}
			} else if ( i == 1 ) {
				if ( mpz_cmp_ui(f[i],1) != 0 ) {
					t += gmp_sprintf (t, " + %Zdx", f[i]);
				} else {
					t += gmp_sprintf (t, " + x");
				}
			} else {
				t += gmp_sprintf (t, " + %Zd", f[i]);
			}
		}
	}
	*t++ = ']';
	*t= '\0';
	return t-s;
}


int ff_poly_sprint (char *s, ff_t f[], int d_f)
{
	char *t;
	int i, n;
	
	if ( d_f < 0 ) { strcpy (s, "[zero polynomial]");  return strlen(s); }
	t = s;
	if ( d_f >= 2 ) {
		t += gmp_sprintf (t, "[%Zdx^%d", _ff_wrap_mpz(f[d_f],0), d_f);
	} else if ( d_f == 1 ) {
		t += gmp_sprintf (t, "[%Zdx", _ff_wrap_mpz(f[d_f],0));
	} else {
		t += gmp_sprintf (t, "[%Zd", _ff_wrap_mpz(f[d_f],0));
	}
	for ( i = d_f-1 ; i >= 0 ; i-- ) {
		if ( _ff_nonzero(f[i]) ) {
			if ( i >= 2 ) {
				t += gmp_sprintf (t, " + %Zdx^%d", _ff_wrap_mpz(f[i],0), i);
			} else if ( i == 1 ) {
				t += gmp_sprintf (t, " + %Zdx", _ff_wrap_mpz(f[i],0));
			} else {
				t += gmp_sprintf (t, " + %Zd", _ff_wrap_mpz(f[i],0));
			}
		}
	}
	*t++ = ']';
	*t= '\0';
	return t-s;
}

// Converts rational coefficients to common denominator, sets f[i] to numerator of F[i] and returns common denominator
unsigned long i_poly_set_mpq (long f[], mpq_t F[], int d)
{
	static mpz_t Q, x;
	static int init;
	unsigned long q;
	int i;
	
	if ( ! init ) { mpz_init (Q);  mpz_init (x);  init = 1; }
	mpq_get_den (Q, F[0]);
	for ( i = 1 ; i <= d ; i++ ) if ( mpq_sgn(F[i]) && mpz_cmp_ui(mpq_denref(F[i]),1) != 0 ) mpz_lcm (Q, Q, mpq_denref(F[i]));
	q = mpz_get_ui (Q);
	// handle integer case quickly
	if ( q == 1 ) {
		for ( i = 0 ; i <= d ; i++ ) {
			if ( mpz_cmpabs_ui (mpq_numref(F[i]), LONG_MAX) > 0 ) return 0;
			f[i] = mpz_get_i (mpq_numref(F[i]));
		}
		return 1;
	}
	for ( i = 0 ; i <= d ; i++ ) {
		if ( mpq_sgn (F[i]) ) {
			mpz_divexact (x, Q, mpq_denref(F[i]));
			mpz_mul (x, x, mpq_numref(F[i]));
			if ( mpz_cmpabs_ui (x, LONG_MAX) > 0 ) return 0;
			f[i] = mpz_get_i (x);
		} else {
			f[i] = 0;
		}
	}
	return q;
}

// Given rational poly g, computes integer poly f s.t. f/denom = F
void mpz_poly_set_mpq (mpz_t denom, mpz_t f[], mpq_t F[], int d)
{
	int i;
	
	mpq_get_den (denom, F[0]);
	for ( i = 1 ; i <= d ; i++ ) if ( mpq_sgn(F[i]) && mpz_cmp_ui(mpq_denref(F[i]),1) != 0 ) mpz_lcm (denom, denom, mpq_denref(F[i]));

	// handle integer case quickly
	if ( mpz_cmp_ui(denom,1) == 0 ) {
		for ( i = 0 ; i <= d ; i++ ) mpz_set (f[i], mpq_numref(F[i]));
	} else {
		for ( i = 0 ; i <= d ; i++ ) {
			if ( mpq_sgn (F[i]) ) {
				mpz_divexact (f[i], denom, mpq_denref(F[i]));
				mpz_mul (f[i], f[i], mpq_numref(F[i]));
			} else {
				mpz_set_ui (f[i], 0);
			}
		}
	}
	return;
}


int ff_poly_set_rational_mpz (ff_t f[], mpz_t F[], int d, mpz_t denom)
{
	_ff_t_declare_reg z;
	int i;
#if FF_INIT_REQUIRED
	static int init;
	if ( ! init ) { _ff_init (z);  init = 1; }
#endif
	
	_ff_set_mpz (z, denom);
	if  ( _ff_zero(z) ) return 0;
	ff_invert (z, z);
	for ( i = 0 ; i <= d ; i++ ) {
		_ff_set_mpz (f[i], F[i]);
		if ( _ff_nonzero(f[i]) ) ff_mult (f[i], f[i], z);
	}
	return 1;
}


int ui_poly_set_rational_mpz_mod_p (unsigned long f[], mpz_t F[], int d, mpz_t denom, unsigned long p)
{
	static mpz_t t;
	static int init;
	unsigned long z;
	int i;
	
	if ( ! init ) { mpz_init (t);  init = 1; }
	z = mpz_fdiv_ui (denom, p);
	z = ui_inverse (z, p);
	if ( ! z ) return 0;
	for ( i = 0 ; i <= d ; i++ ) {
		if ( mpz_sgn(F[i]) ) {
			mpz_mul_ui (t, F[i], z);
			f[i] = mpz_fdiv_ui (t, p);
		} else {
			f[i]  = 0;
		}
	}
	return 1;
}


void ff_poly_eval (ff_t o[1], ff_t f[], int d, ff_t x[1])
{
	_ff_t_declare_reg t,y;
	register int i;
#if FF_INIT_REQUIRED
	static int init;
	if ( ! init ) { _ff_init (t);  _ff_init (y); init = 1; }
#endif
	if ( d < 0 ) { _ff_set_zero (o[0]);  return; }
	
	_ff_set (y,f[d]);
	_ff_set (t,x[0]);
	for ( i = d-1 ; i >= 0 ; i-- ) ff_multadd(y,y,t,f[i]);
	_ff_set(o[0],y);
	return;
}


void ff_poly_twist_mpz (ff_t g[], mpz_t F[], int d)
{
	_ff_t_declare f[POLY_MAX_DEGREE+1];
	int i;
#if FF_INIT_REQUIRED
	static int init;
	if ( ! init ) { for ( i = 0 ; i <= POLY_MAX_DEGREE ; i++ ) _ff_init(f[i]);  init = 1; }
#endif
	if ( d > POLY_MAX_DEGREE ) { err_printf ("ff_poly_twist_mpz: exceeded max degree %d > %d\n", d, POLY_MAX_DEGREE);  exit (0); }
	for ( i = 0 ; i <= d ; i++ ) _ff_set_mpz (f[i], F[i]);
	ff_poly_twist (g, f, d);
	return;
}

/*
	Computes the quadratic twist of y^2=f(x) by a non-residue a in F_p, giving the curve a*y^2=f(x), assuming f(x) has odd degree.

	We put this in a more convenient isomorphic form by replacing y with y/a^{g+1} and x with x/a and multiplying
	both sides by a^{2g+1}.  This will make the result monic if f is, and keep any zero coefficients zero.

	Overlap of f and g is ok.
*/
void ff_poly_twist (ff_t g[], ff_t f[], int d)
{
	_ff_t_declare a, b;
	int i;
#if FF_INIT_REQUIRED
	static int init;
	if ( ! init ) {  _ff_init (a); _ff_init(b);  init = 1; }
#endif
	if ( d > POLY_MAX_DEGREE ) { err_printf ("ff_poly_twist: exceeded max degree %d > %d\n", d, POLY_MAX_DEGREE);  exit (0); }
	if ( ! (d&0x1) ) { err_printf ("ff_poly_twist: degree must be odd.\n");  exit (0); }
	_ff_set (g[d], f[d]);					// leading coefficient unchanged - means twist is monic if f is
	_ff_set (b, _ff_2g);					// use the generator of the 2-Sylow as our non-residue
	_ff_set (a,b);
	for ( i = d-1 ; i >= 0 ; i-- ) {
		ff_mult (g[i], b, f[i]);				// g[i] = a^{d-i}f[i] for i = d-1 down to 0
		ff_mult (b, b, a);
	}
	return;
}

// overlap not ok
void ff_poly_derivative (ff_t g[], int *pd_g, ff_t f[], int d_f)
{
	_ff_t_declare a;
	int i;
#if FF_INIT_REQUIRED
	static int init;
	if ( ! init ) {  _ff_init (a); init = 1; }
#endif
	
	_ff_set_one (a);
	for ( i = 0 ; i < d_f ; i++ ) {
		_ff_mult (g[i], a, f[i+1]);
		_ff_inc(a);
	}
	for ( i = d_f-1 ; i >= 0 && _ff_zero(g[i]) ; i-- );
	*pd_g = i;
}


// overlap not ok - use poly_addto instead
void ff_poly_add (ff_t c[], int *pd_c, ff_t a[], int d_a, ff_t b[], int d_b)
{
	int d_c;
	int i;

	if ( d_a < 0 ) { ff_poly_copy (c, pd_c, b, d_b);  return; }
	if ( d_b < 0 ) { ff_poly_copy (c, pd_c, a, d_a);  return; }
	
	d_c = d_a;
	if ( d_b > d_a ) d_c = d_b;
	for ( i = 0 ; i <= d_c ; i++ ) {
		_ff_set_zero (c[i]);
		if ( i <= d_a ) _ff_addto (c[i], a[i]);
		if ( i <= d_b ) _ff_addto (c[i], b[i]);
	}
	for ( i = d_c ; i >= 0 && _ff_zero(c[i]) ; i-- );
	*pd_c = i;
}

// overlap ok
void ff_poly_addto (ff_t c[], int *pd_c, ff_t b[], int d_b)
{
	int i;

	if ( d_b < 0 ) return;
	for ( i = 0 ; i <= d_b ; i++ ) {
		if ( i > *pd_c ) {
			_ff_set (c[i], b[i]);
			*pd_c = i;
		} else {
			_ff_addto (c[i], b[i]);
		}
	}
	for ( i = *pd_c ; i >= 0 && _ff_zero(c[i]) ; i-- );
	*pd_c = i;
}


// overlap not ok
void ff_poly_neg (ff_t c[], int *pd_c, ff_t b[], int d_b)
{
	int i;

	*pd_c = d_b;
	for ( i = 0 ; i <= d_b ; i++ ) _ff_neg(c[i], b[i]);
}


// overlap not ok
void ff_poly_monic (ff_t c[], int *pd_c, ff_t b[], int d_b)
{
	_ff_t_declare z;
	int i;
#if FF_INIT_REQUIRED
	static int init;
	if ( ! init ) { _ff_init (z);  init = 1; }
#endif
	
	if ( d_b < 0 || _ff_one (b[d_b]) ) { ff_poly_copy (c, pd_c, b, d_b);  return; }
	_ff_invert (z, b[d_b]);
	*pd_c = d_b;
	for ( i = 0 ; i <= d_b ; i++ ) _ff_mult(c[i],z,b[i]);
}


void ff_poly_sub (ff_t c[], int *pd_c, ff_t a[], int d_a, ff_t b[], int d_b)
{
	int d_c;
	int i;

	if ( d_a < 0 ) { ff_poly_neg (c, pd_c, b, d_b);  return; }
	if ( d_b < 0 ) { ff_poly_copy (c, pd_c, a, d_a);  return; }
	
	d_c = d_a;
	if ( d_b > d_a ) d_c = d_b;
	for ( i = 0 ; i <= d_c ; i++ ) {
		_ff_set_zero (c[i]);
		if ( i <= d_a ) _ff_addto (c[i], a[i]);
		if ( i <= d_b ) _ff_subfrom (c[i], b[i]);
	}
	for ( i = d_c ; i >= 0 && _ff_zero(c[i]) ; i-- );
	*pd_c = i;
}


void ff_poly_mult (ff_t c[], int *pd_c, ff_t a[], int d_a, ff_t b[], int d_b)
{
	_ff_t_declare s, t[POLY_MAX_DEGREE+1];
	int d_t;
	int i, j;
#if FF_INIT_REQUIRED
	static int init;
	if ( ! init ) { _ff_init (s);  for ( i = 0 ; i <= POLY_MAX_DEGREE ; i++ ) _ff_init (t[i]); init = 1; }
#endif
	
	if ( d_a < 0 || d_b < 0 ) { *pd_c = -1;  return; }
	if ( _ff_zero(a[d_a]) || _ff_zero(b[d_b]) ) { err_printf ("zero leading coefficient in ff_poly_mult!\n");  ff_poly_print (a, d_a); puts("");  ff_poly_print (b, d_b); exit (0); }
	
	d_t = d_a + d_b;
	for ( i = 0 ; i <= d_t ; i++ ) _ff_set_zero (t[i]);
	
	// multiply into temporary storage to ensure we handle duplicated inputs
	for ( i = 0 ; i <= d_a ; i++ ) {
		for ( j = 0 ; j <= d_b ; j++ ) {
			_ff_mult (s, a[i], b[j]);
			_ff_addto (t[i+j], s);
		}
	}
	ff_poly_copy (c, pd_c, t, d_t);
}


// overlap is ok.  Does not require any inversions if b is monic
void ff_poly_div (ff_t q[], int *pd_q, ff_t r[], int *pd_r, ff_t a[], int d_a, ff_t b[], int d_b)
{
	_ff_t_declare_reg s, t, x;
	ff_t c[POLY_MAX_DEGREE+1];
	register int i, j;
	int d_c, d_q, d_r, d_s;
#if FF_INIT_REQUIRED
	static int init;
	if ( ! init ) { for ( i = 0 ; i <= POLY_MAX_DEGREE ; i++ ) { _ff_init (c[i]); } _ff_init (s);  _ff_init (t);  _ff_init (x); init = 1; }
#endif
	if ( d_b < 0 ) { err_printf ("poly division by zero in ff_poly_div!\n");  exit (0); }
	if ( _ff_zero(b[d_b]) ) { err_printf ("zero leading coefficient in ff_poly_div!\n");  exit (0); }
	if ( d_a >= 0 && _ff_zero(a[d_a]) ) { err_printf ("zero leading coefficient in ff_poly_div!\n");  exit (0); }
	
	ff_poly_copy (c, &d_c, b, d_b);	// copy b to ensure we handle overlapping input/output case
	ff_poly_copy (r, &d_r, a, d_a);
	d_q = -1;
	// the code duplication below saves a mult and/or a test inside the loop
	if ( ! _ff_one(c[d_c]) ) {
		_ff_invert (x, c[d_c]);
		while ( d_r >= d_c ) {
			_ff_mult (s, x, r[d_r]);
			d_s = d_r - d_c;
			if ( d_s > d_q ) {
				for ( j = d_q+1 ; j < d_s ; j++ ) _ff_set_zero (q[j]);
				d_q = d_s;
				_ff_set (q[d_q], s);
			} else {
				_ff_addto (q[d_s], s);
			}
			for ( i = 0 ; i < d_c ; i++ ) {
				_ff_mult (t, s, c[i]);
				_ff_subfrom (r[i+d_s], t);
			}
			_ff_set_zero(r[d_r]);
			while ( d_r && _ff_zero (r[d_r]) ) d_r--;
			if ( _ff_zero (r[d_r]) ) d_r = -1;
		}
	} else {
		while ( d_r >= d_c ) {
			_ff_set (s, r[d_r]);
			d_s = d_r - d_c;
			if ( d_s > d_q ) {
				for ( j = d_q+1 ; j < d_s ; j++ ) _ff_set_zero (q[j]);
				d_q = d_s;
				_ff_set (q[d_q], s);
			} else {
				_ff_addto (q[d_s], s);
			}
			for ( i = 0 ; i < d_c ; i++ ) {
				_ff_mult (t, s, c[i]);
				_ff_subfrom (r[i+d_s], t);
			}
			_ff_set_zero(r[d_r]);
			while ( d_r && _ff_zero (r[d_r]) ) d_r--;
			if ( _ff_zero (r[d_r]) ) d_r = -1;
		}
	}
	*pd_q = d_q;
	*pd_r = d_r;
}


void ff_poly_mod (ff_t g[], int *pd_g, ff_t a[], int d_a, ff_t f[], int d_f)
{
	_ff_t_declare q[POLY_MAX_DEGREE+1];
	int i;
	int d_q;
#if FF_INIT_REQUIRED
	static int init;
	if ( ! init ) { for ( i = 0 ; i <= POLY_MAX_DEGREE ; i++ ) { _ff_init (q[i]); }  init = 1; }
#endif
	ff_poly_div (q, &d_q, g, pd_g, a, d_a, f, d_f);
}

// computes h=af-bg of degree less than max(d_f,d_g) (f and g must both be non-zero).  One of a or b (or both) is a unit
// h is divisible by gcd(f,g) but could be (a lot) bigger than f mod g or g mod f
// this function is most useful when d_f~d_g or both are small
void ff_poly_reduce (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_t g[], int d_g)
{
	_ff_t_declare_reg t0,t1,t2,t3;
	register int i,j;
#if FF_INIT_REQUIRED	
	static int init;
	if ( ! init ) { _ff_init (t0);  _ff_init (t1);  _ff_init(t2);  _ff_init (t3);   init = 1; }
#endif
	
	if ( d_f < 0 || d_g < 0 ) { puts ("Zero polynomial in ff_poly_reduce!");  exit (0); }
	if ( ! d_f || ! d_g ) { _ff_set_one(h[0]);  *pd_h=0; return; }

	// avoid unnecessary copying, but allow overlap of h with f or g
	_ff_set(t0,f[d_f]);
	_ff_set(t1,g[d_g]);
	if ( d_f > d_g ) {
		j = d_f-d_g;
		for ( i = d_f-1 ; i >= j ; i-- ) {
			_ff_mult(t2,t1,f[i]);
			_ff_mult(t3,t0,g[i-j]);
			_ff_sub(h[i],t2,t3);
		}
		for ( ; i >= 0 ; i-- ) ff_mult(h[i],f[i],t1);
		for ( i = d_f-1 ; i >= 0 ; i-- ) if ( ! _ff_zero(h[i]) ) break;
		*pd_h = i;
	} else {
		j = d_g-d_f;
		for ( i = d_g-1 ; i >= j ; i-- ) {
			_ff_mult(t2,t0,g[i]);
			_ff_mult(t3,t1,f[i-j]);
			_ff_sub(h[i],t2,t3);
		}
		for ( ; i >= 0 ; i-- ) ff_mult(h[i],g[i],t0);
		for ( i = d_g-1 ; i >= 0 ; i-- ) if ( ! _ff_zero(h[i]) ) break;
		*pd_h = i;
	}
}

// note that g is not made monic (and need not be, even if a and b both are)
void ff_poly_gcd (ff_t g[], int *pd_g, ff_t a[], int d_a, ff_t b[], int d_b)
{
	_ff_t_declare q[POLY_MAX_DEGREE+1], r[POLY_MAX_DEGREE+1], s[POLY_MAX_DEGREE+1], t[POLY_MAX_DEGREE+1];
	int d_q, d_r, d_s, d_t;
	int i, j;
#if FF_INIT_REQUIRED	
	static int init;
	if ( ! init ) { for ( i = 0 ; i <= POLY_MAX_DEGREE ; i++ ) { _ff_init (q[i]);  _ff_init (r[i]);  _ff_init(s[i]);  _ff_init (t[i]); }  init = 1; }
#endif
	if ( d_a > POLY_MAX_DEGREE || d_b > POLY_MAX_DEGREE ) { err_printf ("Exceeded POLY_MAX_DEGREE in ff_poly_gcd!\n");  exit (0); }
	ff_poly_copy (s, &d_s, a, d_a);
	ff_poly_copy (t, &d_t, b, d_b);
	while ( d_t >= 0 ) {
		ff_poly_div (q, &d_q, r, &d_r, s, d_s, t, d_t);
		ff_poly_copy (s, &d_s, t, d_t);
		ff_poly_copy (t, &d_t, r, d_r);
	}
	ff_poly_copy (g, pd_g, s, d_s);
	if ( *pd_g==0 ) _ff_set_one(g[0]);
	return;
}

// compute gcd using poly_reduce - faster than ff_poly_gcd for small or roughly equal size polys but slower if one poly is much bigger than the other
// uses 2n^2+O(n) M+A and *no* inversions
void ff_poly_gcd_red (ff_t g[], int *pd_g, ff_t a[], int d_a, ff_t b[], int d_b)
{
	_ff_t_declare s[POLY_MAX_DEGREE+1], t[POLY_MAX_DEGREE+1];
	int d_s, d_t;
	int i, j;
#if FF_INIT_REQUIRED	
	static int init;
	if ( ! init ) { for ( i = 0 ; i <= POLY_MAX_DEGREE ; i++ ) { _ff_init (q[i]);  _ff_init (r[i]);  _ff_init(s[i]);  _ff_init (t[i]); }  init = 1; }
#endif
	if ( d_a > POLY_MAX_DEGREE || d_b > POLY_MAX_DEGREE ) { err_printf ("Exceeded POLY_MAX_DEGREE in ff_poly_gcd!\n");  exit (0); }
	if ( d_a < 0 ) { ff_poly_copy(g,pd_g, b,d_b); return; }
	if ( d_b < 0 ) { ff_poly_copy(g,pd_g, a,d_a); return; }

	// reduce into s first, avoid copying
	ff_poly_reduce (s, &d_s, a, d_a, b, d_b);
	if ( d_s < 0 ) { if ( d_a < d_b ) ff_poly_copy (g,pd_g,a,d_a); else ff_poly_copy(g,pd_g,b,d_b); return; }
	if ( d_s == 0 ) { *pd_g = 0;  _ff_set_one(g[0]); return; }
	if ( d_a < d_b ) {
		ff_poly_reduce (t, &d_t, s, d_s, a, d_a);
		if ( d_t < 0 ) { if ( d_a < d_s ) ff_poly_copy (g,pd_g,a,d_a); else ff_poly_copy(g,pd_g,s,d_s); return; }
		if ( d_t == 0 ) { *pd_g = 0;  _ff_set_one(g[0]); return; }
		if ( d_s > d_a ) ff_poly_copy(s,&d_s,a,d_a);		// can't easily avoid a copy in this case
	} else {
		ff_poly_reduce (t, &d_t, s, d_s, b, d_b);
		if ( d_t < 0 ) { if ( d_b < d_s ) ff_poly_copy (g,pd_g,b,d_b); else ff_poly_copy(g,pd_g,s,d_s); return; }
		if ( d_t == 0 ) { *pd_g = 0;  _ff_set_one(g[0]); return; }
		if ( d_s > d_b ) ff_poly_copy(s,&d_s,b,d_b);		// can't easily avoid a copy in this case
	}

	// we now have polys in s and t of degree > 0
	while ( d_s > 0 && d_t > 0 ) {
		if ( d_s < d_t ) {
			ff_poly_reduce (t, &d_t, s, d_s, t, d_t);
		} else {
			ff_poly_reduce (s, &d_s, s, d_s, t, d_t);
		}
	}
	if ( d_s < d_t ) {
		if ( d_s < 0 ) { ff_poly_copy (g,pd_g,t,d_t); return; }
	} else {
		if ( d_t < 0 ) { ff_poly_copy (g,pd_g,s,d_s); return; }
	}
	*pd_g = 0;  _ff_set_one(g[0]);
	return;
}

// [CCANT] algorithm 3.2.2.  None of the polynomials can overlap
void ff_poly_gcdext (ff_t D[], int *pd_D, ff_t u[], int *pd_u, ff_t v[], int *pd_v, ff_t a[], int d_a, ff_t b[], int d_b)
{
	_ff_t_declare q[POLY_MAX_DEGREE+1], r[POLY_MAX_DEGREE+1], t[POLY_MAX_DEGREE+1];
	_ff_t_declare t2[POLY_MAX_DEGREE+1], v1[POLY_MAX_DEGREE+1], v3[POLY_MAX_DEGREE+1];
	int d_q, d_r, d_t, d_t2, d_v1, d_v3;
	int i, j;
#if FF_INIT_REQUIRED	
	static int init;
	if ( ! init ) {
		for ( i = 0 ; i <= POLY_MAX_DEGREE ; i++ )
			{ _ff_init (q[i]);  _ff_init (r[i]);  _ff_init(s[i]);  _ff_init (t[i]); _ff_init (t2[i]); _ff_init (v1[i]); _ff_init (v3[i]); }
		init = 1;
	}
#endif
	// GCd(0,b) = b = 0*a+1*b
	if ( d_a < 0 ) {
		ff_poly_copy (D, pd_D, b, d_b);
		*pd_v = 0;  _ff_set_one(v[0]);
		*pd_u = -1;
		return;
	}
	// GCD(a,0) = a = 1*a+0*b
	if ( d_b < 0 ) {
		ff_poly_copy (D, pd_D, a, d_a);
		*pd_u= 0;  _ff_set_one(u[0]);
		*pd_v = -1;
		return;
	}
	
	if ( d_a > POLY_MAX_DEGREE || d_b > POLY_MAX_DEGREE ) { err_printf ("Exceeded POLY_MAX_DEGREE in ff_poly_gcd!\n");  exit (0); }
	*pd_u = 0;  	_ff_set_one(u[0]);			// u = 1
	ff_poly_copy (D, pd_D, a, d_a);			// D = a
	d_v1 = -1;							// v1 = 0
	ff_poly_copy (v3, &d_v3, b, d_b);			// v3 = b	

	while ( d_v3 >= 0 ) {
		ff_poly_div (q, &d_q, r, &d_r, D, *pd_D, v3, d_v3);
		ff_poly_mult (t2, &d_t2, v1, d_v1, q, d_q);
		ff_poly_sub (t, &d_t, u, *pd_u, t2, d_t2);	// t = u-v1q
		ff_poly_copy (u, pd_u,  v1, d_v1);		// u = v1
		ff_poly_copy (D, pd_D, v3, d_v3);		// D = v3
		ff_poly_copy (v1, &d_v1, t, d_t);		// v1 = t
		ff_poly_copy (v3, &d_v3, r, d_r);		// v3 = r
	}
	ff_poly_mult (v3, &d_v3, a, d_a, u, *pd_u);
	ff_poly_sub (v1, &d_v1, D, *pd_D, v3, d_v3);
	ff_poly_div (v, pd_v, r, &d_r, v1, d_v1, b, d_b);	// v = (D-au)/b 		(exact division)
	if ( d_r != -1 ) { err_printf ("Division not exact in ff_poly_gcdext!  r = ");  ff_poly_print(r,d_r);  exit (0); }
	return;
}


void ff_poly_exp_mod (ff_t g[], int *pd_g, ff_t a[], int d_a, mpz_t e, ff_t f[], int d_f)
{
	_ff_t_declare h[POLY_MAX_DEGREE+1];
	int d_h;
	int i, j, t;
#if FF_INIT_REQUIRED	
	static int init;
	if ( ! init ) { for ( i = 0 ; i <= POLY_MAX_DEGREE ; i++ ) { _ff_init (h[i]); } init = 1; }
#endif
	if ( d_f <= 1 ) { err_printf ("Invalid modulus in ff_poly_exp_mod, degree must be > 1\n");  exit (0); }
	if ( _ff_zero(f[d_f]) ) { err_printf ("zero leading coefficient in ff_poly_exp_mod!\n");  exit (0); }
	if ( ! mpz_sgn (e) ) { _ff_set_one(g[0]);  *pd_g = 0;  return; }
	if ( d_a < 0 ) { *pd_g = 0;  return; }
	if ( _ff_zero(a[d_a]) ) { err_printf ("zero leading coefficient in ff_poly_mult!\n");  exit (0); }
	
	_ff_set_one (h[0]);
	d_h = 0;
	t = mpz_sizeinbase (e, 2) - 1;
	while ( t >= 0 ) {
		if ( mpz_tstbit (e, t) ) {
			ff_poly_mult (h, &d_h, h, d_h, a, d_a);
			ff_poly_mod (h, &d_h, h, d_h, f, d_f);
		}
		if ( t > 0 ) {
			ff_poly_mult (h, &d_h, h, d_h, h, d_h);
			ff_poly_mod (h, &d_h, h, d_h, f, d_f);
		}
		t--;
	}
	ff_poly_copy (g, pd_g, h, d_h);
}

// computes x^n mod f for monic f
void ff_poly_xn_mod (ff_t g[], int *pd_g, mpz_t e, ff_t f[], int d_f)
{
	_ff_t_declare h[POLY_MAX_DEGREE+1], w[2*POLY_MAX_DEGREE];
	_ff_t_declare_reg t0;
	register int i, j, m, n, t;
#if FF_INIT_REQUIRED	
	static int init;
	if ( ! init ) { for ( i = 0 ; i <= POLY_MAX_DEGREE ; i++ ) { _ff_init (h[i]); } init = 1; }
#endif
	if ( d_f <= 1 ) { err_printf ("Invalid modulus in ff_poly_exp_mod, degree must be > 1\n");  exit (0); }
	if ( ! _ff_one(f[d_f]) ) { err_printf ("nonmonic in ff_poly_exp_mod!\n");  exit (0); }
	if ( ! mpz_sgn (e) ) { _ff_set_one(g[0]);  *pd_g = 0;  return; }
	
	n = d_f;
	for ( m = d_f-1 ; _ff_zero(f[m]) ; m-- );						// m is the index of first nonzero nonleading coefficient of f
	_ff_set_one(h[1]);  _ff_set_zero(h[0]);
	for ( i = 2 ; i < n ; i++ ) _ff_set_zero(h[i]);
	t = mpz_sizeinbase (e, 2) - 2;
	while ( t >= 0 ) {
		if ( mpz_tstbit (e, t) ) {
			// square h and multiply it by x, into w
			for ( i = 0 ; i < n ; i++ ) { _ff_set_zero(w[2*i]); _ff_square(w[2*i+1],h[i]); }
			for ( i = 0 ; i < n ; i++ ) for ( j = 0 ; j < i ; j++ ) { _ff_mult(t0,h[i],h[j]); _ff_x2(t0);  _ff_addto(w[i+j+1],t0); }
		} else {
			// square h into w
			for ( i = 0 ; i < n ; i++ ) { _ff_square(w[2*i],h[i]); _ff_set_zero(w[2*i+1]); }
			for ( i = 0 ; i < n ; i++ ) for ( j = 0 ; j < i ; j++ ) { _ff_mult(t0,h[i],h[j]); _ff_x2(t0);  _ff_addto(w[i+j],t0); }
		}
		// reduce w mod f into h
		for ( i = 2*n-1 ; i >= n ; i-- ) if ( ! _ff_zero(w[i]) ) for ( j = 0 ; j <= m ; j++ ) { _ff_mult(t0,w[i],f[j]); _ff_subfrom(w[i-n+j],t0); }
		for ( i = 0 ; i < n ; i++ ) _ff_set(h[i],w[i]);
		t--;
	}
	for ( i = n-1 ; i >= 0 && _ff_zero(h[i]) ; i-- );
	for ( *pd_g = i ; i >= 0 ; i-- ) _ff_set(g[i],h[i]);
}


// Algorithm IPT from Shoup "A Computational Introduction to Number Theory and Algebra" p. 463
int ff_poly_irreducible (ff_t f[], int d_f, int *pnroots)
{
	static mpz_t q;
	static int init;
	_ff_t_declare h[POLY_MAX_DEGREE+1], X[2], t[POLY_MAX_DEGREE+1];
	int d, d_h, d_t, d_X;
	int i;

	if ( ! init ) { mpz_init(q);  init = 1; }

	if ( d_f < 1 ) { err_printf ("Invalid input to ff_poly_irreducible, must have degree > 0\n");  exit (0); }
	if ( ! _ff_one(f[d_f]) ) { err_printf ("Input to ff_poly_irreducible must be monic.\n");  exit (0); }
	if ( d_f == 1 ) {
		if ( pnroots ) *pnroots = 1;
		return 1;
	}
	
	_ff_set_one (h[1]);  _ff_set_zero (h[0]);
	d_h = 1;
	_ff_set_zero (X[0]);  _ff_set_one (X[1]);  ff_negate (X[1]);
	d_X = 1;
	mpz_set(q,_ff_mpz_p);
	for ( i = 1 ; i <= d_f/2 ; i++ ) {
		if ( i > 1 ) mpz_mul(q,q,_ff_mpz_p);
		ff_poly_xn_mod (h, &d_h, q, f, d_f);
		ff_poly_add (t, &d_t, h, d_h, X, d_X);
		ff_poly_gcd (t, &d_t, t, d_t, f, d_f);
		if ( d_t > 0 ) break;
	}
	if ( i <= d_f/2 ) {
		if ( pnroots ) {
			if ( i == 1 ) {
				*pnroots = d_t;
			} else {
				*pnroots = 0;
			}
		}
		return 0;
	}
	return 1;

}

/*
	Counts the distinct factors of f using the standard distinct degree factorization algorithm (just a straight-forward generalization of the irreducibiilty test)
*/
int ff_poly_factors (ff_t f[], int d_f)
{
	static mpz_t q;
	static int init;
	_ff_t_declare g[POLY_MAX_DEGREE+1], h[POLY_MAX_DEGREE+1], X[2], t[POLY_MAX_DEGREE+1], r[POLY_MAX_DEGREE+1], s[POLY_MAX_DEGREE+1], v[POLY_MAX_DEGREE];
	int d, d_g, d_r, d_s, d_h, d_t, d_X, d_v;
	int i, n;

	if ( ! init ) { mpz_init(q);  init = 1; }
	if ( d_f < 1 ) { err_printf ("Invalid input to ff_poly_factors, must have degree > 0\n");  exit (0); }
	if ( ! _ff_one(f[d_f]) ) { err_printf ("Input to ff_poly_factors must be monic.\n");  exit (0); }
	if ( d_f == 1 ) return 1;

	_ff_set_one (h[1]);  _ff_set_zero (h[0]);
	d_h = 1;
	_ff_set_zero (X[0]);  _ff_set_one (X[1]);  ff_negate (X[1]);
	d_X = 1;
	ff_poly_monic (g, &d_g, f, d_f);
	n = 0;
	mpz_set(q,_ff_mpz_p);
	for ( i = 1 ; i <= d_g/2 ; i++ ) {
		if ( i > 1 ) mpz_mul(q,q,_ff_mpz_p);
		ff_poly_xn_mod (h, &d_h, q, g, d_g);					// h = x^{p^i}
		ff_poly_add (t, &d_t, h, d_h, X, d_X);					// t = x^{p^i}-x
		do {
			ff_poly_gcd_red (s, &d_s, t, d_t, g, d_g);
			if ( i > 1 && d_s % i ) { err_printf ("Unexpected gcd in ff_poly_factors, degree not a multiple of i\n");  exit (0); }
			n += d_s / i;
			if ( d_s > 0 ) {
				ff_poly_div (v, &d_v, r, &d_r, g, d_g, s, d_s);
				if ( d_r > -1 ) { err_printf ("Error in ff_poly_factors - division by gcd not exact, remainder:");  ff_poly_print(r, d_r); exit(0); }
				ff_poly_monic(g,&d_g,v,d_v);
			}
		} while ( d_s > 0 );
	}
	if ( d_g ) n++;
	return n;

}

// returns the degree of gcd(f,x^p-x) which will be the number of roots if f is square-free and will always be non-zero if f has a root over F_p
int ff_poly_roots (ff_t f[], int d_f)
{
	_ff_t_declare g[POLY_MAX_DEGREE+1], h[POLY_MAX_DEGREE+1], X[2], t[POLY_MAX_DEGREE+1], r[POLY_MAX_DEGREE+1], s[POLY_MAX_DEGREE+1];
	int d, d_g, d_r, d_s, d_h, d_t, d_X;
	int i, n;

	if ( d_f < 1 ) { err_printf ("Invalid input to ff_poly_has_root, must have degree > 0\n");  exit (0); }
	if ( ! _ff_one(f[d_f]) ) { err_printf ("Input to ff_poly_roots must be monic.\n");  exit (0); }
	if ( d_f == 1 ) return 1;

	_ff_set_one (h[1]);  _ff_set_zero (h[0]);
	d_h = 1;
	_ff_set_zero (X[0]);  _ff_set_one (X[1]);  ff_negate (X[1]);
	d_X = 1;
	ff_poly_monic (g, &d_g, f, d_f);
	ff_poly_xn_mod (h, &d_h, _ff_mpz_p, g, d_g);				// h = x^p
	ff_poly_add (t, &d_t, h, d_h, X, d_X);					// t = x^p-x
	ff_poly_gcd_red (s, &d_s, t, d_t, g, d_g);
	return d_s;
}


// This function only returns bad primes that fit in an unsigned long - larger primes are ignored.
int mpz_poly_bad_primes (unsigned long bad[], int maxbad, mpz_t f[], int d)
{
	static int init;
	static mpz_t D, P;
	unsigned long h[MPZ_MAX_FACTORS];
	int w;
	
	if ( ! init ) { mpz_init (D);  mpz_init (P);  init = 1; }

	mpz_poly_discriminant(D, f, d);
	if ( ! mpz_sgn(D) ) return -1;
	w = mpz_factor_small (bad, h, P, D, _ui_min(maxbad,MPZ_MAX_FACTORS), 128);
	if ( ! w ) { err_printf ("unable to factor_small %Zd\n", D);  return -2; }
	return w;
}


// This function only returns bad primes that fit in an unsigned long - larger primes are ignored.
int mpq_poly_bad_primes (unsigned long bad[], int maxbad, mpq_t f[], int d)
{
	static int init;
	static mpq_t D;
	static mpz_t P, Q;
	unsigned long p1[MPZ_MAX_FACTORS], p2[MPZ_MAX_FACTORS];
	unsigned long h[MPZ_MAX_FACTORS];
	int i, j, k, w1, w2;
	
	if ( ! init ) { mpq_init (D);  mpz_init (P);  mpz_init (Q);  init = 1; }
	mpq_poly_discriminant (D, f, d);
	if ( ! mpq_sgn(D) ) return -1;
	w1 = mpz_factor_small (p1, h, P, mpq_numref(D), MPZ_MAX_FACTORS, 128);
	if ( ! w1 ) { err_printf ("unable to small_factor %Zd\n", mpq_numref(D));  return -2; }
	
	mpq_get_den(Q, f[0]);
	for ( i = 1 ; i <= d ; i++ ) mpz_lcm (Q, Q, mpq_denref(f[i]));
	if ( mpz_cmp_ui (Q,1) == 0 ) {
		for ( i = 0 ; i < w1 ; i++ ) bad[i] = p1[i];
		return w1;
	}
	w2 = mpz_factor_small (p2, h, P, Q, MPZ_MAX_FACTORS, 128);
	if ( ! w2 ) { err_printf ("unable to small_factor %Zd\n", Q);  return -2; }
	
	for ( i = j = k = 0 ; i < w1 || j < w2 ; k++ ) {
		if ( k >= maxbad ) { err_printf ("too many bad primes, increase maxbad value\n");  return -2; }
		if ( i >= w1 ) { bad[k] = p2[j++];  continue; }
		if ( j >= w2 ) { bad[k] = p1[i++];  continue; }
		if ( p1[i] < p2[j] ) {
			bad[k] = p1[i++];
		} else {
			bad[k] = p2[j++];
		}
	}
	return k;
}

/*
	The discriminant functions below are optimized for computing the discriminant of monic polys of degree 3,  5, or 7
	whose d-1 term is 0.  They use an explicit formula generated via Magma with coefficients in polydisc.h.
*/

void mpz_poly_discriminant (mpz_t D, mpz_t f[], int degree)
{
	static int init, warn;
	static mpz_t temp1, temp2, fpow[6][8];
	struct disc_terms *dt;
	int i, j, k;
	
	if ( ! init ) { mpz_init (temp1);  mpz_init (temp2);  for ( i = 0 ; i < 6 ; i++ ) for ( j = 0 ; j <= 7 ; j++ ) mpz_init(fpow[i][j]);  init = 1; }
	if ( ! (degree&0x1) ) { err_printf ("Don't know how to compute discriminant for even degree\n");  exit (0); }
	if ( mpz_cmp_ui (f[degree], 1) != 0 ) { err_printf ("Don't know how to compute discriminant for non-monic poly, degree %d, leading coeff %Zd\n", degree, f[degree]);  exit (0); }
	if ( mpz_sgn(f[degree-1]) ) { err_printf ("Don't know how to compute discriminant for non-zero d-1 term, please standardize poly\n");  exit (0); }
	for ( i = degree-1 ; i > 1 ; i-- ) if ( mpz_sgn(f[i]) ) break;
	if ( i == 1 ) {
		mpz_ui_pow_ui (temp1, (degree-1), (degree-1));
		mpz_pow_ui (D, f[1], degree);
		mpz_mul (D, D, temp1);
		mpz_ui_pow_ui (temp1, degree, degree);
		mpz_pow_ui (temp2, f[0], degree-1);
		mpz_mul (temp2, temp2, temp1);
		mpz_add (D, D, temp2);
		mpz_neg (D, D);
	} else {
		if ( degree != 5 && degree != 7) { err_printf ("Don't know how to compute discriminant for degree != 5 or 7 and not x^d+ax+b\n");  exit (0); }
		
		// not the most efficient way to do this...
		for ( i = 0 ; i < degree-1 ; i++ ) {
			mpz_set_ui (fpow[i][0], 1);
			for ( j = 1 ; j <= degree ; j++ ) mpz_mul (fpow[i][j], fpow[i][j-1], f[i]);
		}
		
		mpz_set_ui (temp1, 0);
		if ( degree == 5 ) {
			for ( i = 0 ; i < PD_D5_TERMS ; i++ ) {
				if ( PD_D5[i].c < 0 ) {
					mpz_mul_ui (temp2, fpow[0][PD_D5[i].f[0]], -PD_D5[i].c);
					mpz_neg(temp2, temp2);
				} else {
					mpz_mul_ui (temp2, fpow[0][PD_D5[i].f[0]], PD_D5[i].c);
				}
				for ( j = 1 ; j < degree-1 ; j++ ) mpz_mul (temp2, temp2, fpow[j][PD_D5[i].f[j]]);
				mpz_add (temp1, temp1, temp2);
			}
		} else if ( degree == 7 ) {
			for ( i = 0 ; i < PD_D7_TERMS ; i++ ) {
				if ( PD_D7[i].c < 0 ) {
					mpz_mul_ui (temp2, fpow[0][PD_D7[i].f[0]], -PD_D7[i].c);
					mpz_neg(temp2, temp2);
				} else {
					mpz_mul_ui (temp2, fpow[0][PD_D7[i].f[0]], PD_D7[i].c);
				}
				for ( j = 1 ; j < degree-1 ; j++ ) mpz_mul (temp2, temp2, fpow[j][PD_D7[i].f[j]]);
				mpz_add (temp1, temp1, temp2);
			}
		}
		mpz_set (D, temp1);
	}
}


// For degree 6 (only) we handle the fully general case.
void mpq_poly_discriminant (mpq_t D, mpq_t f[], int degree)
{
	static int init, warn;
	static mpq_t temp1, temp2, temp3, fpow[7][8];
	static mpz_t z;
	struct disc_terms *dt;
	int i, j, k;
	
	if ( ! init ) { mpz_init (z); mpq_init (temp1);  mpq_init (temp2);  mpq_init (temp3);  for ( i = 0 ; i < 7 ; i++ ) for ( j = 0 ; j <= 7 ; j++ ) mpq_init(fpow[i][j]);  init = 1; }
	if ( degree != 6 ) {
		if ( degree < 2 || ! (degree&0x1) ) { err_printf ("Degree must be odd > 1 to compute discriminant\n");  exit (0); }
		if ( mpq_cmp_ui (f[degree], 1, 1) != 0 ) { err_printf ("Don't know how to compute discriminant for non-monic poly, degree %d, leading coeff %Zd\n", degree, f[degree]);  exit (0); }
		if ( mpq_sgn(f[degree-1]) ) { err_printf ("Don't know how to compute discriminant for non-zero d-1 term, please standardize poly\n");  exit (0); }
		for ( i = degree-1 ; i > 1 ; i-- ) if ( mpq_sgn(f[i]) ) break;
	} else {
		i = 5;
	}
	if ( i == 1 ) {
		mpq_set (D, f[1]);
		for ( i = 1 ; i < degree ; i++ ) mpq_mul (D, D, f[1]);
		mpz_ui_pow_ui (z, (degree-1), (degree-1));
		mpq_set_z (temp1, z);
		mpq_mul (D, D, temp1);
		mpq_set (temp2, f[0]);
		for ( i = 1 ; i < degree-1 ; i++ ) mpq_mul (temp2, temp2, f[0]);
		mpz_ui_pow_ui (z, degree, degree);
		mpq_set_z (temp1, z);
		mpq_mul (temp2, temp2, temp1);
		mpq_add (D, D, temp2);
		mpq_neg (D, D);
	} else {
		if ( degree < 5 || degree > 7) { err_printf ("Don't know how to compute discriminant for degree != 5, 6, or 7 and not x^d+ax+b\n");  mpq_set_ui(D,1,1); return; }
		
		// not the most efficient way to do this...
		k = (degree==6?7:degree-1);
		for ( i = 0 ; i < k ; i++ ) {
			mpq_set_ui (fpow[i][0], 1, 1);
			for ( j = 1 ; j <= degree ; j++ ) mpq_mul (fpow[i][j], fpow[i][j-1], f[i]);
		}
		
		mpq_set_ui (temp1, 0, 1);
		if ( degree == 5 ) {
			for ( i = 0 ; i < PD_D5_TERMS ; i++ ) {
				mpq_set_si (temp3, PD_D5[i].c, 1);
				mpq_mul (temp2, fpow[0][PD_D5[i].f[0]], temp3);
				for ( j = 1 ; j < k ; j++ ) mpq_mul (temp2, temp2, fpow[j][PD_D5[i].f[j]]);
				mpq_add (temp1, temp1, temp2);
			}
		} else if ( degree == 6 ) {
			for ( i = 0 ; i < PD_D6_TERMS ; i++ ) {
				mpq_set_si (temp3, PD_D6[i].c, 1);
				mpq_mul (temp2, fpow[0][PD_D6[i].f[0]], temp3);
				for ( j = 1 ; j < k ; j++ ) mpq_mul (temp2, temp2, fpow[j][PD_D6[i].f[j]]);
				mpq_add (temp1, temp1, temp2);
			}			
		} else if ( degree == 7 ) {
			for ( i = 0 ; i < PD_D7_TERMS ; i++ ) {
				mpq_set_si (temp3, PD_D7[i].c, 1);
				mpq_mul (temp2, fpow[0][PD_D7[i].f[0]], temp3);
				for ( j = 1 ; j < k ; j++ ) mpq_mul (temp2, temp2, fpow[j][PD_D7[i].f[j]]);
				mpq_add (temp1, temp1, temp2);
			}
		}
		mpq_set (D, temp1);
	}
}


/*
	Assumes f is monic, degree 3, 5, or, 7, and has d-1 coefficient zero.
*/
void _ff_poly_discriminant (ff_t disc[1], ff_t f[], int degree)
{
	 _ff_t_declare temp1, temp2, temp3, fpow[6][8];
	struct disc_terms *dt;
	int i, j, k;
#if FF_INIT_REQUIRED
	static int init;	
	if ( ! init ) { _ff_init (temp1);  _ff_init (temp2);  _ff_init(temp3); for ( i = 0 ; i < 6 ; i++ ) for ( j = 0 ; j <= 7 ; j++ ) _ff_init(fpow[i][j]);  init = 1; }
#endif
	if ( degree != 5 && degree != 7) { err_printf ("Don't know how to compute discriminant for degree != 5 or 7 and not x^d+ax+b\n");  exit (0); }
	if ( ! _ff_one(f[degree]) ) { err_printf ("Don't know how to compute discriminant when poly is not monic\n");  exit (0); }
	if ( _ff_nonzero(f[degree-1]) ) { err_printf ("Don't know how to compute discriminant when degree-1 coefficient is non-zero, please standardize poly\n");  exit (0); }
	
	// not the most efficient way to do this...
	for ( i = 0 ; i < degree-1 ; i++ ) {
		_ff_set_one (fpow[i][0]);
		for ( j = 1 ; j <= degree ; j++ ) _ff_mult (fpow[i][j], fpow[i][j-1], f[i]);
	}
	
	_ff_set_zero(temp1);
	if ( degree == 5 ) {
		for ( i = 0 ; i < PD_D5_TERMS ; i++ ) {
			_ff_set_i (temp3, PD_D5[i].c);
			_ff_mult (temp2, temp3, fpow[0][PD_D5[i].f[0]]);
			for ( j = 1 ; j < degree-1 ; j++ ) _ff_mult (temp2, temp2, fpow[j][PD_D5[i].f[j]]);
			_ff_addto (temp1, temp2);
		}
	} else if ( degree == 7 ) {
		for ( i = 0 ; i < PD_D7_TERMS ; i++ ) {
			_ff_set_i (temp3, PD_D7[i].c);
			_ff_mult (temp2, temp3, fpow[0][PD_D7[i].f[0]]);
			for ( j = 1 ; j < degree-1 ; j++ ) _ff_mult (temp2, temp2, fpow[j][PD_D7[i].f[j]]);
			_ff_addto (temp1, temp2);
		}
	}
	_ff_set (disc[0], temp1);
}

// Requires f monic, hardwired for degrees 2, 3, and 4, and also x^d+ax+b cases for d=5 or 7
void ff_poly_discriminant (ff_t disc[1], ff_t f[], int degree)
{
	_ff_t_declare_reg A, B, temp1, temp2, temp3;
	int i;
#if FF_INIT_REQUIRED
	static int init;
	if ( ! init ) { _ff_init(A);  _ff_init(B);  _ff_init (temp1);  _ff_init (temp2);  _ff_init(temp3);  init = 1; }
#endif
#ifndef FF_FAST
	if ( ! _ff_one (f[degree]) ) { err_printf ("Don't know how to compute discriminant for non-monic poly, leading coeff = %Zd\n", _ff_wrap_mpz(f[degree],0));  exit (0); }
#endif
	switch (degree) {
	case 1: _ff_set_one (disc[0]);  return;
	case 2: _ff_square(temp1,f[1]);  _ff_add(temp2,f[0],f[0]);  _ff_x2(temp2); _ff_sub(disc[0],temp1,temp2);  return;
	case 3:
		_ff_square(temp1,f[1]);  ff_mult(temp1,temp1,f[1]);  _ff_x2(temp1);  _ff_x2(temp1);								// temp1 = 4f1^3
		_ff_square(temp2,f[0]); _ff_set_i(temp3,27); ff_mult(temp2,temp2,temp3);										// temp2 = 27f2^2
		_ff_add(disc[0],temp1,temp2);										
		ff_negate(disc[0]);																					// disc = -4f1^3-27f0^2
		if ( ! _ff_zero(f[2]) ) {
			_ff_square(temp2,f[2]);  _ff_mult(temp3,f[2],f[0]);  ff_mult(temp1,temp2,temp3); _ff_x2(temp1);  _ff_x2(temp1);		// temp1 = 4f2^3f0
			_ff_square(temp3,f[1]);  _ff_mult(temp3,temp2,temp3);  _ff_subfrom(temp3,temp1);							// temp3 = f2^2f1^2 - 4f2^3f0
			_ff_set_ui(temp2,18);  _ff_mult (temp1,f[2],temp2);  _ff_mult(temp2,f[1],f[0]); ff_mult(temp1,temp1,temp2);		// temp1 = 18f0f1f2
			_ff_addto(disc[0],temp1);
			_ff_addto(disc[0],temp3);																		// disc = -4f2^3f0+f2^2f1^2+18f2f1f0-4f1^3-27f0^2
		}
		return;
	case 4:
		if ( ! _ff_zero(f[3]) ) {
			_ff_square(A,f[3]);  ff_mult(A,A,f[0]);  _ff_square(temp1,f[1]); _ff_addto(A,temp1);								// A = f3^2f0 + f1^2
			_ff_mult(temp1,f[2],f[0]);  _ff_x2(temp1);  _ff_x2(temp1); _ff_subfrom(A,temp1);								// A =  f3^2f0 + f1^2 - 4f2f0
			_ff_mult(B,f[3],f[1]);  _ff_add(temp1,f[0],f[0]);  _ff_x2(temp1);  _ff_subfrom(B,temp1);							// B = f3f1 - 4f0
		} else {
			 _ff_square(A,f[1]);  _ff_mult(temp1,f[2],f[0]);  _ff_x2(temp1);  _ff_x2(temp1); _ff_subfrom(A,temp1);				// A =  f1^2 - 4f2f0
			 _ff_add(B,f[0],f[0]);  _ff_x2(B);  ff_negate(B);															// B = -4f0
		}
		_ff_square(temp3,f[2]);  _ff_mult(temp1,temp3,f[2]);  ff_mult(temp1,temp1,A); _ff_x2(temp1); _ff_x2(temp1);			// temp1 = 4f2^3A
		_ff_square(temp2,B);  ff_mult(disc[0],temp2,temp3); _ff_subfrom (disc[0],temp1);								// disc = -4f2^3A + f2^2B^2
		ff_mult(temp2,temp2,B);  _ff_x2(temp2);  _ff_x2(temp2); _ff_subfrom(disc[0],temp2);								// disc = -4f2^3A + f2^2B^2 - 4*B^3
		_ff_set_ui(temp1,18); _ff_mult(temp3,temp1,f[0]); _ff_mult(temp1,A,B); ff_mult(temp1,temp1,temp3);					// temp1 = 18f0BA
		_ff_set_ui(temp2,27);  _ff_square(temp3,A);  ff_mult(temp2,temp2,temp3); _ff_subfrom(temp1,temp2);				// temp1 = 18f0BA - 27A^2
		_ff_addto(disc[0],temp1);																			 // disc = -4f2^3A + f2^2B^2 + 18f0BA - B^3 - 27A^2
		return;
	case 5:
		if ( ! _ff_zero(f[4]) || _ff_zero(f[3]) || ! _ff_zero(f[2]) ) {  _ff_poly_discriminant (disc, f, degree);  return; }
		_ff_square(temp1,f[1]);
		ff_square(temp1,temp1);
		ff_mult(temp1,temp1,f[1]);
		_ff_set_i(temp2,-256);
		_ff_mult(disc[0],temp1,temp2);
		_ff_square(temp1,f[0]);
		ff_square(temp1,temp1);
		_ff_set_i(temp2,-3125);
		ff_mult(temp1,temp1,temp2);
		_ff_addto(disc[0],temp1);
		return;
	case 7:
		if ( ! _ff_zero(f[6]) || ! _ff_zero(f[5]) || ! _ff_zero(f[4]) || _ff_zero(f[3]) || ! _ff_zero(f[2]) ) {  _ff_poly_discriminant (disc, f, degree);  return; }
		_ff_square(temp1,f[1]);
		ff_square(temp2,temp1);
		ff_mult(temp1,temp1,temp2);
		ff_mult(temp1,temp1,f[1]);
		_ff_set_i(temp2,-46656);
		_ff_mult(disc[0],temp1,temp2);
		_ff_square(temp1,f[0]);
		ff_square(temp2,temp1);
		ff_mult(temp1,temp1,temp2);
		_ff_set_i(temp2,-823543);
		ff_mult(temp1,temp1,temp2);
		_ff_addto(disc[0],temp1);
		return;
	default:
		err_printf ("Unhandled poly in ff_poly_discriminant\n");  exit (0);
	}
}

